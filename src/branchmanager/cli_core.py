"""BranchManager CLI — lightweight entrypoint for the branchmanager package.

Implements the BranchManager office workflow and technical reporting utilities.
"""
import argparse
import csv
import hashlib
import json
import logging
import os
import re
import shutil
import sys
from pathlib import Path

# When this file is invoked directly, the package root (src)
# may not be on sys.path. Ensure the parent of this `branchmanager` package is
# available so absolute imports like `branchmanager.db.interface` work.
try:
    here = Path(__file__).resolve().parent
    src_root = str(here.parent)
    if src_root not in sys.path:
        sys.path.insert(0, src_root)
except Exception:
    pass  # could not determine package path; continuing without modification

from branchmanager.db.interface import Database
from branchmanager.pipeline import (
    classify,
    itol,
    neighbourhood,
    novelty,
    paper_trail as _paper_trail_module,
    qc,
    quarterly_review,
    tree,
)
from branchmanager.pipeline import cluster_report as _cluster_report
from branchmanager.pipeline import mwl as _mwl
from branchmanager.pipeline import selection_sets as _selection_sets
from branchmanager.pipeline.classify import _derive_db_name as _classify_derive_db_name
from branchmanager.pipeline.collapse import collapse_fasta_within_taxa
from branchmanager.pipeline.workflow_helpers import (
    _assignment_source_is_fasta,
    build_orig_to_short_map as _build_orig_to_short_map_helper,
    build_placement_warning_rows,
    build_selection_decision,
    build_sequence_assessment_rows,
    classification_ids_matching_kingdom as _classification_ids_matching_kingdom_helper,
    collect_db_taxonomy_rows,
    iter_assignment_rows,
    iter_classification_rows,
    load_classification_results_for_dataset,
    load_taxonomy_entries_from_assignments,
    merge_combined_taxonomy_rows,
    prune_dataset_by_kingdom as _prune_dataset_by_kingdom_helper,
    read_combined_taxonomy_ids,
    write_combined_taxonomy_tsv,
    write_baseline_hits_tsv,
    write_placement_warning_tsv,
    write_selection_summary_tsv,
    write_sequence_assessment_tsv,
)
from branchmanager.taxonomy import (
    canonicalise_sequence_id, normalise_domain_query, taxonomy_matches_kingdom,
    compare_expected_taxon, parse_taxon_string, domain_designation,
)
from branchmanager.partner_metadata import load_partner_sequencing_metadata
from branchmanager.run_manifest import RunManifest, utc_now, sha256_file
from branchmanager.project_config import (
    apply_project_config_defaults, config_path_from_argv, load_project_config,
)

def _find_tree_file_in_dir(d: str):
    """Return path to a tree file in directory d if present, preferring current_tree.nwk."""
    p = Path(d)
    for cand in (p / 'current_tree.nwk', p / 'tree' / 'current_tree.nwk'):
        if cand.exists():
            return str(cand)
    # otherwise search for any .nwk or .tree file in the root or organised tree dir
    for base in (p, p / 'tree'):
        for ext in ('*.nwk', '*.tree', '*.tre'):
            found = next(base.glob(ext), None)
            if found:
                return str(found)
    return None


def _configure_logging(outdir: str):
    """Configure root logger: console INFO + optional file DEBUG in outdir.

    The root logger level must be set to DEBUG so that INFO/DEBUG records
    actually reach the handlers — handler-level filtering alone is not enough
    because Python's logging framework gates records at the logger level first.
    """
    logger = logging.getLogger()
    # remove existing handlers to avoid duplicate messages on repeated calls
    for h in list(logger.handlers):
        logger.removeHandler(h)
    # Set root logger to DEBUG so all records flow through to handlers;
    # each handler then applies its own level filter.
    logger.setLevel(logging.DEBUG)
    fmt = logging.Formatter("%(asctime)s %(levelname)s: %(message)s")
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    logger.addHandler(ch)
    try:
        fh = logging.FileHandler(os.path.join(outdir, 'branchmanager.log'))
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    except Exception:
        logger.warning("Could not write log file to %s", outdir)


def _build_orig_to_short(alias_entries):
    return _build_orig_to_short_map_helper(alias_entries, Database(':memory:'))


def _write_id_map_tsv(path: str | Path, entries, *, short_header: str = 'short_id', original_header: str = 'original_header'):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, 'w') as handle:
        handle.write(f'{short_header}\t{original_header}\n')
        for short, original in entries:
            handle.write(f'{short}\t{original}\n')
    return str(p)


def _dataset_sequence_ids(db: Database, dataset: str) -> set[str]:
    with db.connect() as conn:
        cur = conn.cursor()
        cur.execute("SELECT id FROM sequences WHERE dataset = ?", (dataset,))
        return {str(iid) for (iid,) in cur.fetchall()}


def _filter_fasta_to_ids(src: str | Path, dst: str | Path, allowed_ids: set[str]) -> int:
    from branchmanager.utils.fasta import read_fasta, write_fasta
    records = [(h, s) for h, s in read_fasta(str(src)) if str(h) in allowed_ids]
    write_fasta(records, str(dst))
    return len(records)


def _load_qc_rejection_reasons(outdir: str | Path, extra_manifest: str | Path | None = None) -> dict:
    """Look up why a sequence ID never made it into the current run's accepted
    ID set, using the paper-trail/merge-meeting QC manifest one stage up from
    `outdir` (typically `<root>/02_paper_trail_merge_meeting/failed_qc_sequences/
    failed_qc_manifest.tsv`). Returns {sequence_id: 'reason(;reason...)'} for
    every sequence that was rejected before assembly/insertion, so downstream
    warnings can report the real upstream cause instead of a generic
    "not found" message.
    """
    reasons: dict[str, str] = {}
    root = Path(outdir).parent
    candidates = [
        root / '02_paper_trail_merge_meeting' / 'failed_qc_sequences' / 'failed_qc_manifest.tsv',
        Path(outdir) / 'failed_qc_sequences' / 'failed_qc_manifest.tsv',
    ]
    manifest_path = next((c for c in candidates if c.exists()), None)
    if manifest_path is not None:
        try:
            with open(manifest_path, newline='') as handle:
                for row in csv.DictReader(handle, delimiter='\t'):
                    seq_id = str(row.get('SequenceID') or '').strip()
                    if not seq_id:
                        continue
                    status = str(row.get('Status') or '').strip()
                    rec = str(row.get('Recommendation') or '').strip()
                    seq_reasons = str(row.get('Reasons') or '').strip()
                    detail = status or 'failed_qc'
                    if rec:
                        detail = f'{detail}_recommend_{rec.lower()}'
                    if seq_reasons:
                        detail = f'{detail}:{seq_reasons}'
                    reasons[seq_id] = detail
        except Exception:
            logging.getLogger(__name__).warning(
                "[PERFORMANCE REVIEW] Could not read QC rejection manifest %s to annotate partner metadata warnings",
                manifest_path,
            )
    direct_qc_path = Path(outdir) / 'qc_rejections.tsv'
    if direct_qc_path.is_file():
        try:
            with open(direct_qc_path, newline='') as handle:
                for row in csv.DictReader(handle, delimiter='\t'):
                    seq_id = str(row.get('ID') or '').strip().split()[0]
                    if not seq_id:
                        continue
                    detail = str(row.get('Reasons') or 'failed_qc').strip()
                    reasons[seq_id] = f'{reasons[seq_id]}; {detail}' if seq_id in reasons else detail
        except Exception:
            logging.getLogger(__name__).warning(
                "[PERFORMANCE REVIEW] Could not read QC rejection table %s to register failed isolates",
                direct_qc_path,
            )
    if extra_manifest and Path(extra_manifest).is_file():
        try:
            with open(extra_manifest, newline='') as handle:
                for row in csv.DictReader(handle, delimiter='\t'):
                    seq_id = str(row.get('ID') or row.get('SequenceID') or '').strip()
                    if not seq_id:
                        continue
                    detail = str(row.get('Reasons') or row.get('Reason') or 'failed_qc').strip()
                    reasons[seq_id] = f'{reasons[seq_id]}; {detail}' if seq_id in reasons else detail
        except Exception:
            logging.getLogger(__name__).warning(
                '[PERFORMANCE REVIEW] Could not read cohort QC failures from %s', extra_manifest,
            )
    return reasons


def _write_partner_metadata_warnings(path: str | Path, warning_rows):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, 'w') as handle:
        handle.write('SourceID\tReason\n')
        for source_id, reason in warning_rows:
            handle.write(f'{source_id}\t{reason}\n')
    return str(p)


def _write_rejected_before_assembly(path: str | Path, rejected_rows):
    """Write a sidecar of partner-metadata IDs that never produced a marker
    sequence (e.g. every read failed QC before assembly), so downstream
    reports can still list them with an explicit resequence-required status
    instead of silently omitting them.
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, 'w') as handle:
        handle.write('SourceID\tPartnerID\tAlreadySequenced\tSelectedForSequencing\tReason\n')
        for source_id, partner_id, already_sequenced, selected_for_sequencing, reason in rejected_rows:
            handle.write(f'{source_id}\t{partner_id}\t{already_sequenced}\t{selected_for_sequencing}\t{reason}\n')
    return str(p)


_MARKER_REPLACEMENT_REPORT_FIELDS = [
    'sequence_id', 'outcome', 'retained_marker', 'comparison', 'qc_comparison',
    'previous_length', 'incoming_length', 'previous_qc_class', 'incoming_qc_class',
    'previous_qc_reasons', 'incoming_qc_reasons',
    'previous_read_count', 'incoming_read_count', 'primer_read_change',
    'previous_read_ids', 'incoming_read_ids', 'previous_used_read_ids',
    'incoming_used_read_ids', 'previous_read_files',
    'incoming_read_files', 'previous_sequence_sha256', 'incoming_sequence_sha256',
    'details',
]

_MARKER_QC_RANK = {
    'FAIL_QC': 0,
    'QUALITY_UNVERIFIED': 1,
    'PASS_WITH_WARNINGS': 2,
    'PASS_HIGH_CONFIDENCE': 3,
}


def _used_marker_read_count(used_read_ids, read_ids, qc_reasons):
    used = [value.strip() for value in str(used_read_ids or '').split(';') if value.strip()]
    if used:
        return len(used)
    if 'only_one_primer_read_passed_qc' in str(qc_reasons or ''):
        return 1
    return len([value for value in str(read_ids or '').split(';') if value.strip()])


def _merge_identical_marker_provenance(previous, incoming):
    """Record new reads for identical bases without downgrading accepted QC."""
    merged = dict(incoming)
    for key in ('source_read_ids', 'source_used_read_ids', 'source_read_files'):
        values = []
        for item in (previous.get(key), incoming.get(key)):
            values.extend(value.strip() for value in str(item or '').split(';') if value.strip())
        merged[key] = ';'.join(dict.fromkeys(values))

    try:
        checksums = json.loads(previous.get('source_read_checksums_json') or '{}')
        checksums.update(json.loads(incoming.get('source_read_checksums_json') or '{}'))
        merged['source_read_checksums_json'] = json.dumps(checksums, sort_keys=True)
    except (TypeError, ValueError):
        merged['source_read_checksums_json'] = incoming.get('source_read_checksums_json') or (
            previous.get('source_read_checksums_json') or '{}'
        )

    previous_qc = str(previous.get('marker_qc_class') or '').upper()
    incoming_qc = str(incoming.get('marker_qc_class') or '').upper()
    if (
        previous_qc in _MARKER_QC_RANK
        and incoming_qc in _MARKER_QC_RANK
        and _MARKER_QC_RANK[incoming_qc] < _MARKER_QC_RANK[previous_qc]
    ):
        for key in (
            'source_manifest', 'source_sequence_file', 'source_sequence_sha256',
            'marker_qc_class', 'marker_qc_recommendation', 'marker_qc_reasons',
            'manual_review_status', 'chimera_call', 'chimera_score',
        ):
            merged[key] = previous.get(key)
    return merged


def _marker_replacement_report_row(sequence_id, result=None, **overrides):
    result = dict(result or {})
    result.update(overrides)
    previous_qc = str(result.get('previous_qc_class') or '').upper()
    incoming_qc = str(result.get('incoming_qc_class') or '').upper()
    if previous_qc and incoming_qc:
        old_rank = _MARKER_QC_RANK.get(previous_qc, 0)
        new_rank = _MARKER_QC_RANK.get(incoming_qc, 0)
        qc_comparison = 'INCOMING_QC_HIGHER' if new_rank > old_rank else (
            'EXISTING_QC_HIGHER' if old_rank > new_rank else 'SAME_QC_CLASS'
        )
    else:
        qc_comparison = 'QC_COMPARISON_UNAVAILABLE'
    old_reads = int(result.get('previous_read_count') or 0)
    new_reads = int(result.get('incoming_read_count') or 0)
    if old_reads == 1 and new_reads >= 2:
        primer_change = 'ADDITIONAL_PRIMER_READS_ASSEMBLED'
    elif new_reads >= 2:
        primer_change = 'TWO_PRIMER_MARKER'
    elif new_reads == 1:
        primer_change = 'SINGLE_PRIMER_MARKER'
    else:
        primer_change = 'NO_INCOMING_MARKER'
    status = str(result.get('status') or result.get('outcome') or 'NOT_PROCESSED')
    retained = (
        'INCOMING' if status == 'REPLACED'
        else 'EXISTING' if status in {
            'TOO_SHORT', 'FAILED_QC', 'NO_INCOMING_CONSENSUS', 'LOWER_QC_DECLINED',
        } or status == 'NO'
        else 'NONE'
    )
    return {
        'sequence_id': str(sequence_id),
        'outcome': status,
        'retained_marker': result.get('retained_marker') or retained,
        'comparison': result.get('comparison') or '',
        'qc_comparison': result.get('qc_comparison') or qc_comparison,
        'previous_length': result.get('previous_length', ''),
        'incoming_length': result.get('incoming_length', ''),
        'previous_qc_class': previous_qc,
        'incoming_qc_class': incoming_qc,
        'previous_qc_reasons': result.get('previous_qc_reasons') or '',
        'incoming_qc_reasons': result.get('incoming_qc_reasons') or '',
        'previous_read_count': old_reads,
        'incoming_read_count': new_reads,
        'primer_read_change': result.get('primer_read_change') or primer_change,
        'previous_read_ids': result.get('previous_read_ids') or '',
        'incoming_read_ids': result.get('incoming_read_ids') or '',
        'previous_used_read_ids': result.get('previous_used_read_ids') or '',
        'incoming_used_read_ids': result.get('incoming_used_read_ids') or '',
        'previous_read_files': result.get('previous_read_files') or '',
        'incoming_read_files': result.get('incoming_read_files') or '',
        'previous_sequence_sha256': result.get('previous_sequence_sha256') or '',
        'incoming_sequence_sha256': result.get('incoming_sequence_sha256') or '',
        'details': result.get('detail') or result.get('details') or '',
    }


def _write_marker_replacement_report(path: str | Path, rows):
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, 'w', newline='') as handle:
        writer = csv.DictWriter(
            handle, fieldnames=_MARKER_REPLACEMENT_REPORT_FIELDS,
            delimiter='\t', extrasaction='ignore',
        )
        writer.writeheader()
        writer.writerows(rows)
    return str(output)


def _load_rejected_before_assembly(path: str | Path) -> list[dict]:
    p = Path(path)
    if not p.exists():
        return []
    rows = []
    try:
        with open(p, newline='') as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                source_id = str(row.get('SourceID') or '').strip()
                if not source_id:
                    continue
                rows.append({
                    'id': source_id,
                    'partner_id': str(row.get('PartnerID') or source_id).strip() or source_id,
                    'already_sequenced': str(row.get('AlreadySequenced') or '').strip().lower() == 'true',
                    'selected_for_sequencing': str(row.get('SelectedForSequencing') or '').strip().lower() == 'true',
                    'reason': str(row.get('Reason') or '').strip(),
                })
    except Exception:
        logging.getLogger(__name__).warning(
            "[PERFORMANCE REVIEW] Could not read rejected-before-assembly sidecar %s", p
        )
    return rows


def _load_partner_metadata_for_run(args, db: Database, outdir: str, orig_to_short: dict, run_ids):
    metadata_path = getattr(args, 'partner_metadata', None)
    command = getattr(args, 'command', None)
    if not metadata_path:
        if command == 'performance-review':
            raise SystemExit(
                '[PERFORMANCE REVIEW] --partner-metadata / --sequencing-metadata is required: '
                'a cumulative CSV/TSV ledger with sequence IDs, partner acronyms, optional selection commitments, and already-sequenced status.'
            )
        return {}

    log = logging.getLogger(__name__)
    try:
        metadata_rows = load_partner_sequencing_metadata(metadata_path)
    except Exception as e:
        raise SystemExit(f'[PERFORMANCE REVIEW] Failed to read partner metadata {metadata_path}: {e}')

    run_id_set = {str(x) for x in run_ids}
    qc_rejection_reasons = _load_qc_rejection_reasons(
        outdir, getattr(args, 'cohort_qc_failures', None),
    )
    with db.connect() as conn:
        sequence_datasets = {
            str(sequence_id): str(dataset or '')
            for sequence_id, dataset in conn.execute('SELECT id, dataset FROM sequences')
        }
    dataset_roles = db.get_dataset_roles()
    isolate_statuses = db.get_isolate_statuses()
    maldi_sequence_ids = {
        str(item.get('sequence_id') or '')
        for item in db.get_maldi_identifications()
        if item.get('sequence_id')
    }
    resolved_rows = []
    warnings = []
    rejected_before_assembly = []
    markerless_qc_failures = {}
    matched_sources = set()
    for row in metadata_rows:
        source_id = str(row.get('source_id') or '').strip()
        if not source_id:
            continue
        qc_reason = qc_rejection_reasons.get(source_id)
        if qc_reason and source_id not in run_id_set:
            existing = db.resolve_sequence_id(source_id)
            if existing:
                with db.connect() as conn:
                    sequence_row = conn.execute(
                        'SELECT sequence FROM sequences WHERE id=?', (existing,),
                    ).fetchone()
                has_marker = bool(sequence_row and str(sequence_row[0] or '').strip())
                if has_marker:
                    db.record_marker_replacement_attempt(
                        existing, 'FAILED_QC', qc_reason,
                        source_manifest=str(Path(outdir).parent / '02_paper_trail_merge_meeting' / 'run_manifest.json'),
                        source_fasta=str(Path(outdir).parent / '02_paper_trail_merge_meeting' / 'assembled.fasta'),
                        incoming_qc_class='FAIL_QC',
                        incoming_qc_reasons=qc_reason,
                    )
                    warnings.append((source_id, f'incoming_marker_failed_qc_existing_marker_retained:{qc_reason}'))
                else:
                    markerless_qc_failures[existing] = qc_reason
                rejected_before_assembly.append((
                    source_id,
                    row.get('partner_id') or source_id,
                    bool(row.get('selected_for_wgs')),
                    bool(row.get('selected_for_sequencing')),
                    qc_reason,
                ))
                if has_marker:
                    # Keep the accepted marker and its taxonomy; the failed
                    # incoming attempt is recorded separately below.
                    pass
            else:
                run_dataset = getattr(args, 'dataset', 'user')
                db.upsert_dataset_role(run_dataset, 'candidate', genomes_available=False)
                db.insert_sequences([(source_id, '')], dataset=run_dataset)
                sequence_datasets[source_id] = run_dataset
                dataset_roles[run_dataset] = {'role': 'candidate'}
                markerless_qc_failures[source_id] = qc_reason
                rejected_before_assembly.append((
                    source_id,
                    row.get('partner_id') or source_id,
                    bool(row.get('selected_for_wgs')),
                    bool(row.get('selected_for_sequencing')),
                    qc_reason,
                ))
        mapped = orig_to_short.get(source_id)
        if not mapped:
            try:
                cid = canonicalise_sequence_id(source_id)
            except Exception:
                cid = None
            if cid:
                mapped = orig_to_short.get(cid)
        if not mapped and source_id in run_id_set:
            mapped = source_id
        if not mapped:
            existing = db.resolve_sequence_id(source_id)
            if existing:
                mapped = existing
        if not mapped:
            if qc_reason:
                # The failed isolate is still a database entity. Store an
                # empty marker row so a later submission can attach its first
                # passing marker using the same isolate ID.
                run_dataset = getattr(args, 'dataset', 'user')
                db.upsert_dataset_role(run_dataset, 'candidate', genomes_available=False)
                db.insert_sequences([(source_id, '')], dataset=run_dataset)
                sequence_datasets[source_id] = run_dataset
                dataset_roles[run_dataset] = {'role': 'candidate'}
                mapped = source_id
                markerless_qc_failures[mapped] = qc_reason
                rejected_before_assembly.append((
                    source_id, row.get('partner_id') or source_id,
                    bool(row.get('selected_for_wgs')),
                    bool(row.get('selected_for_sequencing')), qc_reason,
                ))
            else:
                warnings.append((source_id, 'metadata_id_not_found_in_project_database_or_current_run'))
                continue

        mapped_dataset = sequence_datasets.get(str(mapped), '')
        if dataset_roles.get(mapped_dataset, {}).get('role') == 'baseline':
            warnings.append((source_id, 'metadata_id_belongs_to_baseline_dataset'))
            continue

        matched_sources.add(str(mapped))
        existing_metadata = db.get_sequencing_metadata_for_ids([mapped]).get(str(mapped), {})
        from branchmanager.partner_metadata import merge_partner_project_metadata
        project_metadata, catalogue_is_gold = merge_partner_project_metadata(
            existing_metadata.get('project_metadata') or {},
            row.get('project_metadata') or {},
        )
        expected_taxon = str(
            (project_metadata.get('partner_taxon_raw') if catalogue_is_gold else '')
            or project_metadata.get('expected_taxon')
            or row.get('expected_taxon') or ''
        ).strip()
        if expected_taxon:
            project_metadata['expected_taxon'] = expected_taxon
        elif not catalogue_is_gold and domain_designation((row.get('project_metadata') or {}).get('partner_taxon_raw')):
            # A newly supplied generic domain label supersedes any stale
            # specific-looking value previously stored for this isolate.
            project_metadata.pop('expected_taxon', None)
        project_metadata['expected_domain'] = str(row.get('expected_domain') or 'Bacteria').strip() or 'Bacteria'
        addition_method = row.get('addition_method')
        if str(mapped) in run_id_set and str(mapped) in maldi_sequence_ids:
            addition_method = 'both'
        if addition_method:
            project_metadata['addition_method'] = addition_method
        lifecycle_status = str(isolate_statuses.get(str(mapped), {}).get('status') or '').upper()
        if lifecycle_status == 'SAB_APPROVED':
            selected_commitment = True
        elif lifecycle_status in {'SAB_DECLINED', 'SAB_DEFERRED'}:
            selected_commitment = False
        else:
            # A catalogue's stale 'no' must not erase a later SAB approval.
            selected_commitment = bool(row.get('selected_for_sequencing')) or bool(existing_metadata.get('selected_for_sequencing'))
        resolved_rows.append({
            'id': mapped,
            'partner_id': row.get('partner_id') or source_id,
            # Preserve the dataset in which an earlier isolate entered the
            # project instead of relabelling it as part of the current batch.
            'dataset': mapped_dataset or getattr(args, 'dataset', ''),
            'selected_for_sequencing': selected_commitment,
            # A registered genome is factual evidence and must not be erased
            # by a stale catalogue row. New genome availability is monotonic.
            'selected_for_wgs': bool(row.get('selected_for_wgs')) or bool(existing_metadata.get('selected_for_wgs')),
            'source_id': source_id,
            'source_file': str(metadata_path),
            'raw_selected_value': row.get('raw_selected_value', ''),
            'raw_commitment_value': (
                existing_metadata.get('raw_commitment_value', '')
                if existing_metadata.get('selected_for_sequencing') and not row.get('selected_for_sequencing')
                else row.get('raw_commitment_value', '')
            ),
            'project_metadata': project_metadata,
        })

    for run_id in sorted(run_id_set - matched_sources):
        warnings.append((run_id, 'run_sequence_missing_from_partner_metadata'))

    inserted = db.upsert_sequencing_metadata(resolved_rows)
    if markerless_qc_failures:
        provenance_rows = []
        for sequence_id, reason in markerless_qc_failures.items():
            provenance_rows.append({
                'sequence_id': sequence_id,
                'source_manifest': str(Path(outdir).parent / '02_paper_trail_merge_meeting' / 'run_manifest.json'),
                'source_sequence_file': '',
                'marker_qc_class': 'FAIL_QC',
                'marker_qc_recommendation': 'RESEQUENCE',
                'marker_qc_reasons': reason,
                'manual_review_status': 'NOT_REVIEWED',
                'chimera_call': 'NOT_RUN',
            })
        db.upsert_sequence_provenance(provenance_rows)
        for sequence_id, reason in markerless_qc_failures.items():
            # A markerless isolate must not carry a classifier assignment from
            # an earlier incomplete/legacy import. MALDI evidence is separate.
            with db.connect() as conn:
                sequence_row = conn.execute(
                    'SELECT sequence FROM sequences WHERE id=?', (sequence_id,),
                ).fetchone()
                if not sequence_row or str(sequence_row[0] or '').strip():
                    continue
                conn.execute('DELETE FROM taxonomy WHERE id=?', (sequence_id,))
                conn.execute('DELETE FROM taxonomy_alt WHERE id=?', (sequence_id,))
                conn.execute('DELETE FROM distances WHERE id=?', (sequence_id,))
                conn.execute('DELETE FROM classification_evidence WHERE sequence_id=?', (sequence_id,))
            db.update_isolate_status(
                sequence_id, 'MARKER_QC_FAILED', detail=reason,
                source_file=str(Path(outdir) / 'qc_rejections.tsv'),
            )
    log.info(
        '[PERFORMANCE REVIEW] Loaded partner sequencing metadata from %s: %d matched rows, %d warning(s)',
        metadata_path,
        inserted,
        len(warnings),
    )
    if warnings:
        warn_path = _write_partner_metadata_warnings(
            Path(outdir) / 'partner_metadata_warnings.tsv',
            warnings,
        )
        log.warning('[PERFORMANCE REVIEW] Partner metadata warnings written to %s', warn_path)
    if rejected_before_assembly:
        rejected_path = _write_rejected_before_assembly(
            Path(outdir) / 'rejected_before_assembly.tsv',
            rejected_before_assembly,
        )
        log.warning(
            '[PERFORMANCE REVIEW] %d partner isolate(s) produced no marker sequence and will be carried into '
            'the final reports as resequence-required placeholders: %s',
            len(rejected_before_assembly), rejected_path,
        )

    return db.get_sequencing_metadata_for_ids(run_ids)


def _load_id_map_from_tsv(path: str | Path, db: Database | None = None):
    mapping = {}
    with open(path) as handle:
        next(handle, None)
        for line in handle:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            short, orig = parts[0], parts[1]
            if not short or not orig:
                continue
            mapping[orig] = short
            mapping[short] = short
            if db is not None:
                try:
                    cid = db._canonical_from_header(orig)
                    if cid:
                        mapping[cid] = short
                except Exception:
                    pass  # could not resolve canonical ID for alias entry
    return mapping


def _find_preferred_id_map(directory: str | Path):
    p = Path(directory)
    preferred = [
        p / 'filing_cabinet_id_map.tsv',
        p / 'user_id_map.tsv',
        p / 'user_id_map.csv',
        p / 'ids' / 'filing_cabinet_id_map.tsv',
        p / 'ids' / 'user_id_map.tsv',
        p / 'ids' / 'user_id_map.csv',
        p / 'ids' / 'baseline_id_map.tsv',
    ]
    for cand in preferred:
        if cand.exists():
            return cand
    try:
        for base in (p, p / 'ids'):
            for cand in base.glob('*_id_map.tsv'):
                if cand.name != 'id_map.tsv':
                    return cand
    except Exception:
        pass  # could not glob ID map files in directory
    return None


def _combined_taxonomy_candidates(directory: str | Path):
    p = Path(directory)
    return [
        p / 'filing_cabinet_combined_taxonomy.tsv',
        p / 'combined_taxonomy.tsv',
        p / 'taxonomy' / 'tree_taxonomy.tsv',
        p / 'baseline' / 'loaded_baseline_taxonomy.tsv',
    ]


def _read_qc_stats(outdir: Path) -> dict:
    stats = {}
    for cand in (outdir / 'assessment' / 'qc.stats', outdir / 'qc.stats'):
        if not cand.exists():
            continue
        try:
            with open(cand) as handle:
                for line in handle:
                    line = line.rstrip('\n')
                    if not line or '\t' not in line:
                        continue
                    key, value = line.split('\t', 1)
                    stats[key] = value
            stats['_path'] = str(cand)
            return stats
        except Exception:
            return {}
    return {}


def _rel_output_path(root: Path, path: str | Path) -> str:
    try:
        return str(Path(path).relative_to(root))
    except Exception:
        return str(path)


def _write_detailed_output_guide(outdir: Path, rows):
    """Write a detailed Markdown guide for files and core metrics."""
    guide = outdir / 'OUTPUT_GUIDE.md'
    qc_stats = _read_qc_stats(outdir)
    qc_rejections = outdir / 'assessment' / 'qc_rejections.tsv'
    if not qc_rejections.exists():
        qc_rejections = outdir / 'qc_rejections.tsv'

    def stat(name: str, default: str = 'NA') -> str:
        return qc_stats.get(name, default)

    file_rows = []
    for name, path, expl in rows:
        rel = _rel_output_path(outdir, path)
        if rel in ('OUTPUT_GUIDE.md',):
            continue
        safe_expl = str(expl).replace('|', '/').replace('\t', ' ')
        file_rows.append(f"| `{rel}` | {safe_expl} |")

    qc_section = [
        "## QC Filtering",
        "",
        "Generic QC happens before classification, novelty scoring, and tree building. "
        "By default, sequences below min_len or above max_n_percent are excluded. "
        "However, an explicit APPROVE_FAIL_QC decision from the validated cohort marker-review file "
        "can retain a failed marker for downstream assessment. Such a sequence remains technically FAIL_QC; "
        "the reviewer override and original QC reasons are recorded separately.",
        "",
        "Filtering rules:",
        "",
        f"- `too_short`: sequence length is less than `min_len` (`{stat('min_len')}` bp for this run).",
        f"- `too_many_n`: ambiguous `N` bases exceed `max_n_percent` (`{stat('max_n_percent')}`% for this run).",
        "- A sequence can have more than one reason in `qc_rejections.tsv`.",
        "",
    ]
    if qc_stats:
        qc_section.extend([
            "This run:",
            "",
            f"- Input sequences: `{stat('total_input')}`",
            f"- Kept after QC: `{stat('kept')}`",
            f"- Passed the generic length/N filter: `{stat('passed_generic_qc')}`",
            f"- Retained because of manual failed-QC approval: `{stat('retained_by_manual_override')}`",
            f"- Rejected total: `{stat('rejected_total', str(max(0, int(stat('total_input', '0')) - int(stat('kept', '0'))) if stat('total_input', '0').isdigit() and stat('kept', '0').isdigit() else 'NA'))}`",
            f"- Rejected as too short: `{stat('rejected_too_short')}`",
            f"- Rejected for too many `N` bases: `{stat('rejected_too_many_n')}`",
            "",
        ])
        if qc_rejections.exists():
            qc_section.append(
                f"See `{_rel_output_path(outdir, qc_rejections)}` for every sequence failing the generic filter, "
                "including its disposition (`EXCLUDED` or `INCLUDED_BY_MANUAL_REVIEW`), length, N count, and exact reason."
            )
        else:
            qc_section.append("Older runs may only include the summary file `qc.stats`; rerun to generate per-sequence `qc_rejections.tsv`.")
        override_audit = outdir / 'assessment' / 'qc_manual_overrides.tsv'
        if override_audit.is_file():
            qc_section.append(
                f"Explicit cohort approvals are listed separately in `{_rel_output_path(outdir, override_audit)}`."
            )
        qc_section.append("")
    else:
        qc_section.extend([
            "No `qc.stats` file was found in this output directory. If this was Filing Cabinet only or a partially completed Performance Review, QC may not have been executed.",
            "",
        ])

    lines = [
        "# BranchManager Output Guide",
        "",
        "This guide explains the output files and metrics produced by the BranchManager Performance Review.",
        "",
        "## Run Scope",
        "",
        "Assessment, selection, per-database taxonomy, novelty, baseline-hit, and neighbourhood target rows are limited to the dataset supplied to this run. Prior partner datasets and cultured baselines remain available as comparison context. Neighbourhood context leaves, the combined tree taxonomy, MSA, tree, and iTOL membership files are cumulative project-context outputs by design.",
        "",
        "## Recommended Reading Order",
        "",
        "1. `assessment/sequencing_sets.tsv` - proposed primary, backup, and alternate isolates grouped by GTDB taxon rank and local tree clade.",
        "1a. `assessment/pangenome_hierarchy.tsv` - independent species-to-domain pangenome summaries wherever catalogue lineage ranks are available.",
        "2. `assessment/selection_summary.tsv` - decision-facing per-isolate recommendations and evidence.",
        "3. `assessment/neighbourhoods/` - labelled local-clade figures linked from both reports.",
        "4. `assessment/sequence_assessment.tsv` - full per-sequence audit table.",
        "5. `baseline/baseline_hits.tsv`, `taxonomy/`, `tree/`, and `itol/` - supporting assignments and phylogeny.",
        "6. `assessment/novelty_metrics.tsv` - detailed component calculations.",
        "7. `performance_review_dashboard.html` - interactive, offline-capable dashboard for collection-level interpretation and per-isolate evidence exploration.",
        "",
        "## Directory Overview",
        "",
        "- `assessment/`: selection/audit reports, local-clade figures, novelty metrics, cluster reports, and warnings.",
        "- `baseline/`: nearest-hit reports and loaded baseline taxonomy/sequences.",
        "- `taxonomy/`: per-database classification outputs and combined taxonomy used for tree metadata.",
        "- `quality/`: marker provenance and reference-UCHIME results carried into selection evidence.",
        "- `tree/`: MSA, Newick tree, and tree/alignment warning files.",
        "- `itol/`: one iTOL metadata dataset per metadata type, usually colour-strip files.",
        "- `ids/`: short ID to original FASTA header maps.",
        "- `intermediate/`: QC, optional collapse, debug FASTA files, and scratch pools.",
        "- `logs/`: pipeline log.",
        "",
        *qc_section,
        "## Main Assessment Metrics",
        "",
        "`assessment/selection_summary.tsv` is the decision-facing table. It separates scientific value from redundancy and marker-evidence quality:",
        "",
        "- `Recommendation`: primary, backup, secondary, review, target-met, or already-sequenced action; the reason is explicit in `RecommendationReason`.",
        "- `SequencingSetID`, `SelectionGroupType`, `SequencingSetRole`, `SequencingSetRank`: connect each isolate to its nine-member diversity panel. `BMEXT_*` identifies a baseline-pangenome extension; `BMSET_*` identifies a candidate-only group. `PRIMARY` rows fill the target, `BACKUP` rows protect against DNA-extraction failure, and redundancy/boundary rows are retained without a rank.",
        "- `BaselineExtensionStatus`: records whether a species-level candidate also meets the cultured-baseline identity and query-coverage criteria for a baseline extension; this secondary evidence does not override the GTDB-defined pangenome taxon.",
        "- `PartnerTaxon`, `PartnerTaxonRank`, `PangenomeRank`, `PangenomeTaxon`: the partner-catalogue declaration and the GTDB-derived rank/taxon used for pangenome grouping. `PartnerTaxon` controls pathway routing only; `GTDBGenomeTaxonomy` is the pangenome taxonomy source. `ClassifierGTDBDisagreement` flags marker-classifier versus GTDB discordance.",
        "- `PangenomeCommittedGenomeCount`, `PangenomeBaselineGenomeCount`, `PangenomeTarget`, and `PangenomeGap`: coverage at the stated GTDB rank, subject to tree/clade concordance. See `pangenome_hierarchy.tsv` for parallel rank-level summaries.",
        "- `EvidenceQuality`: whether the marker classification and placement warnings are strong enough to support selection.",
        "- `MarkerQC`, `MarkerReview`: Paper Trail/Merge Meeting evidence and any explicit manual decision. Unreviewed warning/unverified markers cannot become PRIMARY/BACKUP.",
        "- `CulturedGap`: distance from the cultured baseline: large below 97 percent, moderate from 97 to below 98.65 percent, small at or above 98.65 percent.",
        "- `ProjectCoverage`: number of other partner candidates at or above 97 percent identity and 95 percent query coverage, including selected-for-sequencing partner isolates, the current rolling collection, and excluding the query itself.",
        "- `ReferenceContext`: divergence from the external reference using 94.5 and 98.65 percent full-length 16S heuristic boundaries.",
        "- `LocalTreeFigure`: relative path to the grouped local-clade PNG for visual inspection.",
        "",
        "These boundaries are decision-support heuristics, not declarations of a new species or genus. Genome analysis remains necessary for formal novelty claims.",
        "",
        "`assessment/sequence_assessment.tsv` is the main table. Important column groups:",
        "",
        "- `GTDBTaxonomy`, `GTDBClassificationHit`, `GTDBClassificationIdentity`, `GTDBQueryCoverage`, `GTDBTargetCoverage`, `GTDBAlignmentLength`: GTDB's sequence-derived classifier assignment and the extent of its supporting alignment. These fields are classifier evidence, not the partner-isolate identity. `GTDBAssignmentEvidenceLevel`, `GTDBAlignmentEvidenceStatus`, `GTDBEvidenceTier`, and `GTDBEvidenceInterpretation` qualify that evidence.",
        "- `Taxonomy_<DB>`, `ClassificationHit_<DB>`, `Identity_<DB>`, and `SourceClassifierConfidence_<DB>`: assignment from additional databases such as GG2, SILVA, or NCBI. Source classifier confidence is retained only when supplied by that classifier; it is never substituted with percentage identity.",
        "- `CulturedProjectNearestHit`, `CulturedProjectNearestHitDataset`, `CulturedProjectNearestHitTaxonomy`, `CulturedProjectNearestIdentity`, `CulturedProjectNearestQueryCoverage`: closest sequence in the combined cultured-project pool; inspect identity together with query coverage.",
        "- `ProjectNearestHit`, `ProjectNearestIdentity`, `ProjectNoveltyScore`: comparison against all partner candidates in the rolling project collection, excluding self.",
        "- `GTDBReferenceNearestHit`, `GTDBReferenceNearestIdentity`, `GTDBReferenceNoveltyScore`: nearest hit and context against the GTDB reference FASTA.",
        "- `PartnerTaxon`, `PartnerTaxonRank`, `PangenomeRank`, `PangenomeTaxon`: partner-catalogue declaration and GTDB-derived pangenome rank/taxon. `PartnerTaxon` controls domain/pathway routing only; `GTDBGenomeTaxonomy` supplies pangenome identity. `ClassifierGTDBDisagreement` and `ClassifierDatabaseConcordanceStatus` retain independent marker-classifier diagnostics. The legacy `*SameAssessmentSpecies` columns remain classifier-derived context; use `PangenomeCommittedGenomeCount`, `PangenomeTarget`, and `PangenomeGap` for GTDB-rank coverage.",
        "- `SequencingSet*`: proposed primary/backup membership based on the GTDB rank-specific pangenome gap, marker evidence quality, and phylogenetic spread. `pangenome_hierarchy.tsv` contains parallel species-to-domain summaries.",
        "- `InTree`, `ClusterRepresentative`, `ClusterSize`, `ClusteredMembers`: whether the sequence itself entered the tree or was represented by another clustered sequence.",
        "- `PlacementFlags`: warnings such as low classification identity, low nearest identity, or novelty/classification disagreement.",
        "",
        "## Novelty Metrics",
        "",
        "`assessment/novelty_metrics.tsv` contains cultured-baseline novelty, rolling project novelty, GTDB-reference context, and genome-collection coverage.",
        "",
        "- Leading `Nearest*`, `NoveltyScore`, `Crowding`, and `SequencingPriority` columns mirror the cultured-project comparison when a baseline pool exists; otherwise they mirror project novelty.",
        "- `Hungate*`, `SecondaryBaseline*`, and `CulturedProject*` columns keep priority, secondary, and combined cultured-project evidence separate.",
        "- `Project*` columns compare against all partner candidate datasets, including partner isolates selected for sequencing and the current rolling collection, while excluding each query's self-hit. Baseline datasets are not mixed into this pool.",
        "- `Reference*` columns compare against the chosen external reference FASTA supplied with `--ref`, usually GTDB. These are separate from the baseline/project novelty scores.",
        "- `GenomeCollection*` columns compare against every baseline genome and every partner isolate with a genome already available. Exact same-GTDB-species counts drive the nine-genome target, while baseline identity and query coverage prevent near-identical cultured isolates from consuming diversity-panel ranks.",
        "- `NearestIdentity`: vsearch global-alignment percent identity to the highest-identity sequence with at least 95% query coverage; if no sequence meets that coverage, the nearest hit is `None` and unsupported partial matches do not drive novelty.",
        "- `NearestQueryCoverage`, `NearestAlignmentLength`: how much marker evidence supports that identity. Nearest hits are selected by highest identity among alignments covering at least 95% of the query, then by coverage and alignment length for ties.",
        "- `Novel`: `True` when no hit reaches 95 percent query coverage or the best qualifying hit is below 97 percent identity.",
        "- `MatchesGE99`, `MatchesGE97`, `MatchesGE95`: number of pool sequences at or above 99, 97, and 95 percent identity with at least 95 percent query coverage. These describe how busy the supported local neighbourhood is.",
        "- `Crowding`: `isolated` when there is at most one hit at both 99 and 97 percent; `sparse` when <=3 hits at 97 percent; `moderate` when <=10 hits at 97 percent; otherwise `crowded`. Only alignments covering at least 95 percent of the query count.",
        "- `NoveltyScore`: 0-100 score where higher means more novel and less crowded. It combines distance from the nearest hit with density bonuses for sparse neighbourhoods.",
        "- `SequencingPriority`: `HIGH` for <97 percent identity among hits with at least 95 percent query coverage and few close neighbours, `MEDIUM` for moderately novel/sparse cases, otherwise `LOW`.",
        "- `DensitySource`: names the pool used, for example `baseline:Hungate`, `project_collection`, `target_fasta`, or `reference_fasta` fallback.",
        "",
        "Interpretation: a sequence distant from one named baseline but close to another cultured-project sequence may be novel relative to the first baseline, but not relative to all project evidence.",
        "",
        "## Taxonomy Metrics",
        "",
        "- `ClassificationIdentity` is the vsearch percent identity to the best reference hit used for taxonomy assignment.",
        "- `SourceClassifierConfidence` is retained only when the source classifier explicitly provides it. It is not inferred from percentage identity and should not be interpreted as a calibrated probability unless the source classifier documents that meaning.",
        "- `QueryCoverage`, `TargetCoverage`, `AlignmentLength`, `QueryLength`, `TargetLength`, `Mismatches`, and `Gaps` expose the alignment behind each best hit. Identity is calculated only across the aligned region; coverage fields show how much of the query and reference that region spans. `AssignmentEvidenceLevel` and `EvidenceInterpretation` summarise the limitations without converting identity into a probability.",
        "- `AssignmentEvidenceLevel` is an evidence category, not a probability: `HIGH` requires high identity (>=98%), >=95% query coverage, >=90% reference coverage, and a confident taxonomic consensus; `MODERATE` indicates intermediate alignment or provisional consensus; `LOW` indicates weak coverage/identity or contested taxonomy; `TENTATIVE` is diagnostic-only evidence; `UNRESOLVED` means alignment evidence is insufficient. These are reporting bands, not species-identification cut-offs.",
        "- `TaxonomySupportPct` must be read with `TaxonomyCandidateCount`: 100% support from one reference is not independent confirmation. `EvidenceInterpretation` states the actual identity, query/reference coverage, search tier, and candidate count.",
        "- `taxonomy/all_databases.tsv` shows assignments across all configured databases in one place.",
        "- `taxonomy/input_warnings.tsv` flags mismatches between reference FASTA IDs and supplied taxonomy tables.",
        "",
        "## MWL Metrics",
        "",
        "When `--mwl` is supplied, MWL columns are added to `sequence_assessment.tsv` and `assessment/mwl_matches.tsv` is written.",
        "",
        "- `MWLMatch`: whether the GTDB-derived taxonomy matched a Most Wanted List taxon.",
        "- `MWLMatchedRank` and `MWLMatchedTaxon`: deepest rank/taxon that matched the MWL entry.",
        "- `MWLTaxonomicScore`, `MWLIdentity`, `MWLScore`: taxonomic and identity contributions to MWL priority.",
        "- `EvaluationScore`: combined candidate score after MWL contribution is considered.",
        "",
        "## Tree And iTOL Files",
        "",
        "- `tree/current_alignment.fasta`: MSA used to build the tree.",
        "- `tree/current_tree.nwk`: tree to upload to iTOL.",
        "- `itol/phylum.itol`, `itol/family.itol`, `itol/genus.itol`: taxonomy colour strips.",
        "- `itol/dataset_membership.itol`: dataset-of-origin colour strip.",
        "- `itol/marker_replacement_status.itol`: red tree-tip markers for isolates whose new marker failed QC and whose previous marker remains in use.",
        "- `itol/novelty.itol`: novelty/nearest-identity colour strip.",
        "- `assessment/neighbourhoods/clade_*.png`: compact labelled subtrees with selection-set role/rank, MSA pident to the P1 isolate, and red triangles for retained old markers after failed resequencing QC. Nearest baseline hits are forced into the displayed context.",
        "- `assessment/neighbourhoods/clade_*_pairwise_pident.tsv`: all displayed leaf-to-leaf MSA percent identities and compared-column counts.",
        "- `assessment/neighbourhoods/neighbourhood_manifest.tsv`: maps each assessed ID to its figure/pident table, P1 identity anchor, forced baseline hits, and displayed context.",
        "- `assessment/marker_replacement_report.tsv`: per-isolate replacement outcome, previous/incoming QC, marker lengths, read support, and which marker was retained. `assessment/replaced_markers.fasta` contains the exact incoming sequences applied. Replacement consensus minimum is the normal 800 bp QC threshold.",
        "- `assessment/marker_replacement_history.tsv`: append-only replacement attempt history with outcome, sequence hashes, QC details, and source paths; successful replacements are also retained in the project database.",
        "",
        "BranchManager keeps iTOL `DATASET_COLORSTRIP` files only. Older branch `TREE_COLORS` and symbol-strip variants duplicated the same metadata and are removed to keep outputs readable.",
        "",
        "## File Catalogue",
        "",
        "| File | Explanation |",
        "|---|---|",
        *file_rows,
        "",
    ]
    try:
        guide.write_text('\n'.join(lines))
    except Exception:
        logging.getLogger(__name__).debug('Could not write OUTPUT_GUIDE.md')


def _write_output_explanations(outdir: str):
    """Write short human-readable explanation files for common outputs under outdir.

    For each file discovered under outdir we create a companion file named
    <filename>.explain.txt containing a brief description of the file and how
    it is produced. This helps downstream users and CI to understand outputs.
    """
    p = Path(outdir)
    if not p.exists():
        return
    # mapping of filename substrings (lowercase) to explanation text
    patterns = [
        ('tree_taxonomy.tsv', 'Combined taxonomy TSV used for tree/iTOL metadata. Columns: ID\tTaxon\tConfidence.'),
        ('combined_taxonomy.tsv', 'Combined taxonomy TSV. Columns: ID\tTaxon\tConfidence. Used to generate iTOL colour/legend files.'),
        ('loaded_baseline_taxonomy.tsv', 'Combined taxonomy TSV for baseline/provided datasets loaded before Performance Review. Columns: ID\tTaxon\tConfidence.'),
        ('filing_cabinet_combined_taxonomy.tsv', 'Combined taxonomy TSV for the Filing Cabinet baseline. Columns: ID\tTaxon\tConfidence.'),
        ('all_databases.tsv', 'Per-sequence taxonomic assignments from every configured reference database.'),
        ('baseline_hits.tsv', 'Nearest-hit report against baseline/provided datasets such as Hungate.'),
        ('nearest_project_hits_raw.tsv', 'Raw nearest-hit table produced for the newest submission; rolling project and baseline comparisons are summarised in novelty_metrics.tsv.'),
        ('qc_rejections.tsv', 'Per-sequence generic QC decision table. Columns include ID, Length, NCount, NPercent, Reasons, MinLength, MaxNPercent, ManualReviewDecision, and Disposition. Manually approved failed-QC markers are retained and marked INCLUDED_BY_MANUAL_REVIEW; excluded records are marked EXCLUDED.'),
        ('qc_manual_overrides.tsv', 'Audit table of explicit APPROVE_FAIL_QC decisions passed through the generic QC filter, including sequence length, N content, automatic QC reasons, and whether a manual override was needed or the sequence already passed the generic filter.'),
        ('partner_metadata_warnings.tsv', 'Warnings from --partner-metadata mapping. Rows indicate metadata IDs not found in the project database/current run, baseline IDs supplied as partner metadata, or current-run sequences missing from the metadata table.'),
        ('expected_taxonomy_validation.tsv', 'Per-sequence comparison of the immutable partner-catalogue taxon with the VSEARCH-derived lineage. Classifier disagreements are reported but do not alter the catalogue taxon or exclude the isolate.'),
        ('marker_qc_provenance.tsv', 'Per-isolate bridge from Paper Trail/Merge Meeting raw reads and marker QC into Performance Review, including source checksums, manual review, and chimera call.'),
        ('manual_qc_override_warnings.tsv', 'High-severity audit of manually approved isolates that still carry failed automated marker QC or chimera evidence. Confirmed chimera overrides are critical and remain excluded from automated selection.'),
        ('marker_replacement_report.tsv', 'Per-isolate replacement summary comparing prior and incoming QC class, length, read sources, final outcome, and which marker was retained.'),
        ('replaced_markers.fasta', 'Exact incoming marker sequences applied as replacements in this Performance Review run.'),
        ('marker_replacement_history.tsv', 'Append-only record of marker replacement attempts and successful replacements, including outcomes, old/new sequence hashes, QC reasons, and source paths. The project database retains this history across runs.'),
        ('marker_replacement_status.itol', 'iTOL colour strip: red marks a tree leaf whose attempted replacement failed QC and whose previous marker is retained.'),
        ('chimera_screen.tsv', 'Reference-UCHIME marker screen. CHIMERA, INDETERMINATE, SKIPPED, or NOT_RUN calls force selection evidence to review.'),
        ('decision_changes.tsv', 'Difference between the two latest stored assessment snapshots, including recommendation, set role, pangenome gap, and tracked evidence changes.'),
        ('run_manifest.json', 'Machine-readable workflow manifest with input/output SHA256 checksums, software/tool versions, stage status, timestamps, warnings, and failure details.'),
        ('performance_review_dashboard.html', 'Interactive, self-contained dashboard for collection-level interpretation, taxonomy evidence, pangenome coverage, and per-isolate exploration.'),
        ('qc.stats', 'QC summary with total input, kept count, generic-QC passes, manually retained failed-QC markers, actual exclusions, min_len, max_n_percent, and pointers to the per-sequence QC audit tables.'),
        ('user_id_map.tsv', 'Mapping of runtime sequence IDs to original headers produced when inserting user sequences into the DB. When shortening is disabled these usually match.'),
        ('filing_cabinet_id_map.tsv', 'Mapping of Filing Cabinet runtime IDs back to original FASTA headers. Use this to trace tree labels back to source records.'),
        ('*_id_map.tsv', 'ID map mapping original headers to runtime DB ids. Useful for iTOL and metadata tracing.'),
        ('collapsed_map.tsv', 'Cluster map for collapsed sequences: rep_id\ttaxonomy\tcount.'),
        ('filing_cabinet_collapsed_map.tsv', 'Cluster map for collapsed Filing Cabinet sequences: rep_id\ttaxonomy\tcount.'),
        ('collapsed_members.tsv', 'Member->representative mapping (member\trep) for collapsed clusters.'),
        ('filing_cabinet_collapsed_members.tsv', 'Member-to-representative mapping for collapsed Filing Cabinet clusters.'),
        ('current_dataset_sequences.fasta', 'QC-passing current-dataset FASTA with the runtime IDs used by the database. Exact marker duplicates retain separate isolate IDs.'),
        ('current_dataset_collapsed.fasta', 'Optional collapsed tree-view FASTA; assessment and database rows still retain every current-dataset isolate ID.'),
        ('filing_cabinet_collapsed.fasta', 'Collapsed Filing Cabinet FASTA; representatives retained for tree building.'),
        ('novelty_matches.tsv', 'vsearch BLAST-like output used to compute nearest-neighbour novelty identities.'),
        ('novelty_supplementary_reference_summary.tsv', 'Nearest-hit checks against supplementary 16S references such as GG2, SILVA, and NCBI. These confirm marker-sequence recognition only; they do not contribute to genome coverage or the primary novelty score.'),
        ('novelty_marker_density_coverage_comparison.tsv', 'GTDB-reference and project density counts compared at 95% and 80% query coverage.'),
        ('novelty_baseline_density_coverage_comparison.tsv', 'Hungate and secondary baseline density counts compared at 95% and 80% query coverage.'),
        ('vsearch_search_audit.tsv', 'Per-invocation VSEARCH parameters, reference checksum, sequence counts, runtime, version, and raw hit count.'),
        ('novelty_metrics.tsv', (
            'Per-sequence novelty metrics. The leading Nearest*/NoveltyScore columns mirror the '
            'cultured-project comparison when available; explicit Hungate*, SecondaryBaseline*, and '
            'CulturedProject* columns separate baseline-tier evidence, and Project* columns compare against all partner candidates '
            'stored in the DB. Reference* columns compare against the external reference FASTA, usually '
            'GTDB. GenomeCollection* and Pangenome* columns add rolling genome-coverage context. PartnerTaxon is the immutable partner-catalogue gold standard; classifier disagreement is reported separately. '
            'DensitySource columns name the comparison pool.'
        )),
        ('sequence_assessment.tsv', (
            'Unified per-sequence assessment. '
            'COLUMN GROUPS: '
            '(1) TAXONOMY/CLASSIFICATION — GTDBTaxonomy, GTDBClassificationHit, '
            'GTDBClassificationIdentity, GTDBQueryCoverage, GTDBTargetCoverage, GTDBAlignmentLength, '
            'GTDBSourceClassifierConfidence, GTDBAssignmentEvidenceLevel, '
            'GTDBAlignmentEvidenceStatus, GTDBEvidenceTier, GTDBEvidenceInterpretation: authoritative GTDB assignment and qualified alignment evidence. '
            'Repeated as Taxonomy_<DB>, ClassificationHit_<DB>, Identity_<DB>, SourceClassifierConfidence_<DB> '
            'for each additional --alt-ref database. '
            '(2) NOVELTY — NearestHit, NearestIdentity, MatchesGE*, NoveltyScore, Crowding, '
            'SequencingPriority: baseline/cultured novelty when a baseline pool exists. '
            'Project* columns repeat the same metrics against the rolling partner collection. '
            'Reference* columns repeat the same metrics against the external taxonomy reference '
            'FASTA, usually GTDB. GenomeCollection* and Pangenome* columns report available '
            'baseline genomes, already-sequenced partner genomes, and remaining same-species coverage. '
            '(3) TREE/CLUSTER — InTree, ClusterRepresentative, ClusterSize, ClusteredMembers: '
            'records whether the sequence entered the phylogenetic tree directly or was '
            'represented by a cluster representative after --collapse. '
            '(4) MWL — when --mwl is supplied, EvaluationScore and MWL* columns describe '
            'GTDB-based Most Wanted List matches and MWL priority contribution.'
        )),
        ('selection_summary.tsv', (
            'Concise scientific-advisory-board selection table. Contains one row per assessed '
            'sequence with a transparent decision, marker evidence quality, cultured gap, project '
            'coverage, external-reference context, taxonomy/MWL evidence, available-genome coverage, '
            'local-tree figure, and a short rationale. '
            'Use sequence_assessment.tsv for the full audit trail.'
        )),
        ('pangenome_hierarchy.tsv', 'Rank-stratified pangenome assessment from species through domain using GTDB taxonomy and rank-matched tree/clade concordance. Marker-classifier disagreements are reported but never replace GTDB assignments.'),
        ('sequencing_sets.tsv', (
            'Rolling GTDB-rank genome-sequencing diversity panel. Every candidate is '
            'placed in a BMEXT baseline-pangenome extension or a BMSET candidate-only group. A baseline '
            'where classifier and catalogue species agree, baseline-extension evidence also considers marker identity and query coverage. '
            'Classifier discordance is reported but does not reassign or remove the catalogue-defined group. PRIMARY rows '
            'fill the pangenome target; BACKUP and DIVERSITY_CANDIDATE rows preserve extraction-failure '
            'resilience and additional within-group spread. BASELINE_REDUNDANT rows remain auditable but '
            'are excluded using the reported identity and query-coverage thresholds; '
            'PANGENOME_BOUNDARY_REVIEW rows require lineage review before grouping.'
        )),
        ('neighbourhood_manifest.tsv', 'Maps each assessed sequence to its grouped local-clade image and pairwise-pident table, including the P1 identity anchor, forced nearest-baseline hits, displayed leaf count, assessed peers, baseline leaves, and already-sequenced leaves.'),
        ('visual_report_manifest.tsv', 'Page index for Paper Trail PNG reports. Records each page file, read/isolate range, dimensions, configured height ceiling, and any split-isolate continuation.'),
        ('read_error_profiles_page_', 'Height-bounded page of per-read Phred quality, trim windows, and internally masked low-quality regions.'),
        ('trace_chromatograms_page_', 'Height-bounded page of AB1 dye-channel traces, retained windows, and mixed-peak markers.'),
        ('assembly_overview_page_', 'Height-bounded page aligning every read to consensus coordinates, with read contribution, assembly status, and QC decision.'),
        ('assembly_read_placements.tsv', 'Auditable per-read consensus coordinates and contribution status used by the assembly overview figures.'),
        ('pairwise_pident.tsv', 'Long-form pairwise MSA percent identity for every displayed leaf pair. Identity is identical A/C/G/T bases divided by jointly unambiguous A/C/G/T MSA columns; ComparableACGTColumns reports the overlap denominator. Terminal gaps and ambiguous bases are excluded.'),
        ('clade_', 'Local phylogenetic-neighbourhood image with full sequence IDs, recommendation role/rank labels, MSA pident to P1, dataset/taxonomy context, forced nearest-baseline markers, already-sequenced markers, and primary/backup recommendation stars.'),
        ('clusters.csv', 'Consolidated cluster-membership table containing every cluster and member isolate, representative status, backup rank, taxonomy, novelty, crowding, phylogenetic isolation, and placement fields.'),
        ('cluster_summary.tsv', 'One row per cluster with consensus taxonomy, novelty/crowding summaries, phylogenetic isolation, investigation score, and backup availability.'),
        ('mwl_matches.tsv', 'Most Wanted List match report. Contains sequences whose GTDB taxonomy matched an MWL taxon, with matched rank, MWL score, evaluation score, and functional role.'),
        ('taxonomy_input_warnings.tsv', 'Warnings about inconsistencies between the classifier reference FASTA and supplied taxonomy table.'),
        ('tree_build_warnings.tsv', 'Warnings about weak phylogenetic signal, missing anchors, or poor alignment quality.'),
        ('tree_orientation_summary.tsv', 'Sequence-level audit of tree-input orientation checks. Reports which sequences were kept forward, reverse-complemented, or lacked orientation evidence before alignment.'),
        ('placement_warnings.tsv', 'Warnings about low-support placements, low identity matches, or potentially artefactual novelty assignments.'),
        ('taxa_assignments_classout.tsv', 'Synthetic classification-like TSV created when --taxa-assignments provided. Columns: id\tbest\tidentity\ttaxon\tconfidence.'),
        ('dataset_membership.itol', 'iTOL DATASET_COLORSTRIP mapping sequence IDs to dataset colours (membership).'),
        ('novelty.itol', 'iTOL colour strip showing novelty (nearest identity) for run sequences.'),
        ('filing_cabinet_dataset.itol', 'iTOL colour strip for the Filing Cabinet dataset.'),
        ('itol_dataset_membership.itol', 'iTOL DATASET_COLORSTRIP mapping sequence IDs to dataset colours (membership).'),
        ('itol_novelty.itol', 'iTOL colour strip showing novelty (nearest identity) for run sequences.'),
        ('.nwk', 'Newick tree file (phylogenetic tree). Commonly named current_tree.nwk.'),
        ('.itol', 'iTOL dataset file (text format) describing colours/strips/legends for visualisation in iTOL.'),
        ('.log', 'Log file produced by the pipeline (branchmanager.log) containing debug/info messages')
    ]

    rows = []
    # iterate files in outdir (non-recursive and recursive) and collect explanations
    for fp in sorted(p.rglob('*')):
        if fp.is_dir():
            continue
        if fp.name.endswith('.explain.txt'):
            try:
                fp.unlink()
            except Exception:
                pass  # could not remove old .explain.txt file
            continue
        name = fp.name
        lname = name.lower()
        expl = None
        for patt, text in patterns:
            # treat wildcard at start/end
            if patt.startswith('*') and patt.endswith('*'):
                key = patt.strip('*').lower()
                if key in lname:
                    expl = text
                    break
            elif patt.startswith('*'):
                key = patt.lstrip('*').lower()
                if lname.endswith(key):
                    expl = text
                    break
            elif patt.endswith('*'):
                key = patt.rstrip('*').lower()
                if lname.startswith(key):
                    expl = text
                    break
            elif patt.startswith('.'):
                if lname.endswith(patt):
                    expl = text
                    break
            else:
                if patt in lname:
                    expl = text
                    break
        if expl is None:
            # generic descriptions based on extension
            if lname.endswith('.fasta') or lname.endswith('.fa'):
                expl = 'FASTA file containing sequences. May be dereplicated, collapsed, or exported from the DB.'
            elif lname.endswith('.tsv') or lname.endswith('.csv'):
                expl = 'Tab/CSV-delimited table used by the pipeline.'
            elif lname.endswith('.uc'):
                expl = 'vsearch UC membership file describing clusters.'
            else:
                expl = 'Pipeline output file.'

        rows.append((fp.name, str(fp), expl))

    try:
        manifest = p / 'OUTPUT_EXPLANATIONS.tsv'
        with open(manifest, 'w') as ef:
            ef.write('File\tPath\tDescription\n')
            for name, path, expl in rows:
                ef.write(f"{name}\t{path}\t{expl}\n")
    except Exception:
        logging.getLogger(__name__).debug('Could not write OUTPUT_EXPLANATIONS.tsv')
    _write_detailed_output_guide(p, rows)


def _replace_path(src: Path, dst: Path):
    if not src.exists():
        return None
    try:
        if src.resolve() == dst.resolve():
            return str(dst)
    except Exception:
        pass  # could not resolve paths for same-file check
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        if dst.is_dir():
            shutil.rmtree(dst)
        else:
            dst.unlink()
    shutil.move(str(src), str(dst))
    return str(dst)


def _move_if_exists(outdir: Path, name: str, target_dir: str, target_name: str | None = None):
    return _replace_path(outdir / name, outdir / target_dir / (target_name or name))


def _move_glob(outdir: Path, pattern: str, target_dir: str, rename=None):
    moved = []
    for src in sorted(outdir.glob(pattern)):
        if not src.exists() or src.is_dir():
            continue
        target_name = rename(src) if rename else src.name
        dst = outdir / target_dir / target_name
        got = _replace_path(src, dst)
        if got:
            moved.append(got)
    return moved


def _unlink_glob(outdir: Path, pattern: str):
    for src in sorted(outdir.glob(pattern)):
        try:
            if src.is_dir():
                shutil.rmtree(src)
            else:
                src.unlink()
        except Exception:
            pass  # could not remove file during glob cleanup


def _organise_run_outputs(outdir: str, *, primary_db_name: str | None = None):
    """Clean and organise Performance Review outputs into report folders."""
    out = Path(outdir)
    if not out.exists():
        return

    # Remove transient classifier outputs that are represented in final reports.
    for patt in (
        'novelty_matches.tsv',
        'novelty_density_matches.tsv',
        'novelty_*_matches.tsv',
        'ref_uncompressed*.fasta',
    ):
        _unlink_glob(out, patt)

    # Overall assessment and prioritisation reports.
    for name in (
        'sequence_assessment.tsv',
        'marker_replacement_report.tsv',
        'marker_replacement_history.tsv',
        'replaced_markers.fasta',
        'cluster_summary.tsv',
        'clusters.csv',
        'backup_candidates.tsv',
        'placement_warnings.tsv',
        'novelty_metrics.tsv',
        'selection_summary.tsv',
        'sequencing_sets.tsv',
        'pangenome_hierarchy.tsv',
        'rumen_functions_draft.tsv',
        'qc.stats',
        'qc_rejections.tsv',
        'qc_manual_overrides.tsv',
        'manual_qc_override_warnings.tsv',
        'partner_metadata_warnings.tsv',
        'expected_taxonomy_validation.tsv',
        'rejected_before_assembly.tsv',
    ):
        _move_if_exists(out, name, 'assessment')
    for name in (
        'marker_qc_provenance.tsv', 'chimera_screen.tsv', 'chimera_uchime.tsv',
        'chimera_flagged.fasta', 'chimera_passed.fasta',
    ):
        _move_if_exists(out, name, 'quality')
    _replace_path(out / 'neighbourhoods', out / 'assessment' / 'neighbourhoods')
    _move_if_exists(out, 'mwl_matches.tsv', 'assessment')

    # Direct baseline/provided-dataset nearest-hit reports.
    _move_if_exists(out, 'baseline_hits.tsv', 'baseline')
    _move_if_exists(out, 'novelty.tsv', 'assessment', 'nearest_project_hits_raw.tsv')

    # Taxonomy assignment reports.
    primary_name = primary_db_name or 'primary'
    safe_primary = ''.join(c if c.isalnum() or c in ('_', '-') else '_' for c in str(primary_name)).strip('_') or 'primary'
    _move_if_exists(out, 'taxonomy.tsv', 'taxonomy', f'{safe_primary}.tsv')
    _move_if_exists(out, 'taxonomy_all_dbs.tsv', 'taxonomy', 'all_databases.tsv')
    _move_if_exists(out, 'combined_taxonomy.tsv', 'taxonomy', 'tree_taxonomy.tsv')
    _move_if_exists(out, 'taxonomy_input_warnings.tsv', 'taxonomy', 'input_warnings.tsv')
    _move_if_exists(out, 'taxa_assignments_classout.tsv', 'taxonomy', 'input_assignments.tsv')
    # Preserve the raw VSEARCH evidence, labelled by database, alongside the
    # corresponding taxonomy assignment tables instead of deleting it as a temp file.
    _move_if_exists(out, 'matches.tsv', 'taxonomy')
    _move_glob(out, 'matches_*.tsv', 'taxonomy')

    # Keep the database-qualified filenames intact in the final report folder.
    _move_glob(out, 'taxonomy_*.tsv', 'taxonomy')
    _move_glob(out, 'taxonomy_conflicts_*.tsv', 'taxonomy')
    _move_glob(out, 'classification_supporting_hits_*.tsv', 'taxonomy')

    # Tree/MSA deliverables for iTOL upload.
    for name in (
        'current_tree.nwk',
        'current_tree_labelled.nwk',
        'current_alignment.fasta',
        'tree_build_warnings.tsv',
        'tree_orientation_summary.tsv',
    ):
        _move_if_exists(out, name, 'tree')
    _move_glob(out, 'tree_sequences_phylum_*.fasta', 'tree')

    # One iTOL metadata dataset per metadata type.
    itol_renames = {
        'itol_phylum_colours.itol': 'phylum.itol',
        'itol_family_colours.itol': 'family.itol',
        'itol_genus_colours.itol': 'genus.itol',
        'itol_dataset_membership.itol': 'dataset_membership.itol',
        'itol_baseline_tier.itol': 'baseline_tier.itol',
        'itol_marker_replacement_status.itol': 'marker_replacement_status.itol',
        'itol_novelty.itol': 'novelty.itol',
        'itol_user_colours.itol': 'user_colours.itol',
        'filing_cabinet_dataset.itol': 'filing_cabinet_dataset.itol',
    }
    for src_name, dst_name in itol_renames.items():
        _move_if_exists(out, src_name, 'itol', dst_name)

    def _rename_function_itol(src: Path):
        stem = src.stem
        if stem.startswith('itol_func_'):
            stem = 'function_' + stem[len('itol_func_'):]
        return stem + src.suffix

    _move_glob(out, 'itol_func_*.itol', 'itol', rename=_rename_function_itol)

    # Stable ID maps and run-sequence processing artefacts.
    for name in ('user_id_map.tsv', 'filing_cabinet_id_map.tsv'):
        _move_if_exists(out, name, 'ids')

    # Inline Filing Cabinet internals: keep useful reports, drop raw matches/ref cache.
    baseline_filing_cabinet = out / 'filing_cabinet_baseline'
    if baseline_filing_cabinet.exists():
        _replace_path(baseline_filing_cabinet / 'baseline_id_map.tsv', out / 'ids' / 'baseline_id_map.tsv')
        _replace_path(
            baseline_filing_cabinet / 'baseline_combined_taxonomy.tsv',
            out / 'baseline' / 'loaded_baseline_taxonomy.tsv',
        )
        _replace_path(
            baseline_filing_cabinet / 'taxonomy.tsv',
            out / 'baseline' / 'loaded_baseline_classification.tsv',
        )
        for src in sorted(baseline_filing_cabinet.glob('filing_cabinet_*_sequences.fasta')):
            _replace_path(src, out / 'baseline' / src.name)
        try:
            shutil.rmtree(baseline_filing_cabinet)
        except Exception:
            pass  # could not remove old baseline_filing_cabinet directory

    for name in (
        'qc.fasta',
        'current_dataset_sequences.fasta',
        'current_dataset_collapsed.fasta',
        'collapsed_map.tsv',
        'collapsed_members.tsv',
        'submitted_sequences.fasta',
        'project_collection_reference.fasta',
        'density_query_db.fasta',
        'db_sequences.fasta',
        'combined_input.fasta',
        'new_sequences.fasta',
        'combined_aln.fasta',
        'tree_orientation_ref.fasta',
        'id_map.tsv',
    ):
        _move_if_exists(out, name, 'intermediate')
    _move_glob(out, '*_oriented.fasta', 'intermediate')
    _move_glob(out, '*_unique.fasta', 'intermediate')
    _move_glob(out, 'novelty_*_pool.fasta', 'intermediate')
    _move_glob(out, '*.uc', 'intermediate')

    _move_if_exists(out, 'branchmanager.log', 'logs')


def _classification_ids_matching_kingdom(classification_tsv: str, kingdom: str):
    return _classification_ids_matching_kingdom_helper(classification_tsv, kingdom)


def _validate_partner_taxonomy(classification_tsv: str, partner_metadata_path: str,
                               kingdom: str, outdir: str, input_fasta: str | None = None):
    """Retain unclassified partner isolates, but reject explicit taxonomic conflict.

    The partner catalogue's ``taxon`` column is an expected identification, not
    a replacement for sequence-derived taxonomy. A missing VSEARCH assignment
    therefore remains unverified rather than being treated as non-bacterial.
    """
    metadata_rows = load_partner_sequencing_metadata(partner_metadata_path)
    metadata_by_id = {}
    for row in metadata_rows:
        sid = str(row.get('source_id') or '').strip()
        if not sid:
            continue
        metadata_by_id[sid] = row
        canonical = canonicalise_sequence_id(sid)
        if canonical:
            metadata_by_id.setdefault(canonical, row)

    retained = set()
    report_rows = []
    requested_domain = normalise_domain_query(kingdom)
    results_by_id = {}
    for result in iter_classification_rows(classification_tsv):
        result_id = str(result.get('qid') or '')
        if result_id:
            results_by_id[result_id] = result
            canonical_result_id = canonicalise_sequence_id(result_id)
            if canonical_result_id:
                results_by_id.setdefault(canonical_result_id, result)

    # VSEARCH/classifier outputs may omit queries with no accepted hit. Those
    # sequences must still be evaluated: for partner isolate collections, a
    # missing hit is unclassified evidence, not evidence of a non-bacterial
    # sequence. Iterate the input FASTA so no-hit queries are retained under
    # the declared bacterial expectation and explicitly reported.
    if input_fasta:
        from branchmanager.utils.fasta import read_fasta
        query_ids = [str(qid) for qid, _sequence in read_fasta(input_fasta)]
    else:
        query_ids = list(dict.fromkeys(
            str(row.get('qid') or '') for row in results_by_id.values() if row.get('qid')
        ))

    for qid in query_ids:
        result = results_by_id.get(qid) or results_by_id.get(canonicalise_sequence_id(qid) or '')
        result = result or {'qid': qid, 'tax': None}
        qid = str(qid)
        metadata = metadata_by_id.get(qid) or metadata_by_id.get(canonicalise_sequence_id(qid) or '')
        expected_taxon = str((metadata or {}).get('expected_taxon') or '').strip()
        expected_domain = str((metadata or {}).get('expected_domain') or 'Bacteria').strip() or 'Bacteria'
        # Defensive normalisation for metadata already persisted in the project
        # database before generic domain labels were normalised on import.
        taxon_domain = domain_designation(expected_taxon)
        if taxon_domain:
            expected_domain = taxon_domain
            expected_taxon = ''
        observed_taxonomy = str(result.get('tax') or '').strip()
        if observed_taxonomy.lower() in {'', 'na', 'none', 'unresolved'}:
            status = 'UNCLASSIFIED_EXPECTED_BACTERIAL'
            detail = 'No usable VSEARCH taxonomy; retained under the partner-isolate bacterial expectation'
            retain = True
        else:
            observed_ranks = parse_taxon_string(observed_taxonomy)
            explicit_domain = observed_ranks.get('d') or observed_ranks.get('k')
            explicit_domain_matches = (
                taxonomy_matches_kingdom(observed_taxonomy, requested_domain)
                if explicit_domain else None
            )
            if explicit_domain and not explicit_domain_matches:
                status = 'DOMAIN_DISAGREEMENT'
                detail = (
                    f'Classifier lineage conflicts with the requested {requested_domain} domain; '
                    'partner-catalogue identity remains authoritative and the sequence is retained'
                )
                retain = True
            elif expected_taxon:
                comparison = compare_expected_taxon(expected_taxon, observed_taxonomy)
                if comparison == 'DISAGREEMENT':
                    status = 'TAXON_DISAGREEMENT'
                    detail = (
                        'Sequence-derived lineage conflicts with the partner-catalogue taxon; '
                        'the catalogue value is retained as the gold-standard identity'
                    )
                    retain = True
                elif comparison == 'LESS_SPECIFIC':
                    status = 'COMPATIBLE_LESS_SPECIFIC'
                    detail = 'Observed lineage is compatible but less specific than the partner-catalogue taxon'
                    retain = True
                else:
                    status = 'CONFIRMED_COMPATIBLE' if explicit_domain_matches else 'TAXON_COMPATIBLE_DOMAIN_UNRESOLVED'
                    detail = 'Observed lineage is compatible with the partner-catalogue taxon'
                    retain = True
            elif explicit_domain_matches:
                status = 'ASSIGNED_REQUESTED_DOMAIN'
                detail = 'Observed lineage explicitly belongs to the requested domain; no partner taxon was supplied'
                retain = True
            else:
                status = 'UNRESOLVED_DOMAIN_EXPECTED_BACTERIAL'
                detail = 'Taxonomy is present but does not provide an explicit conflicting domain; retained for review'
                retain = True
        if retain:
            retained.add(qid)
        report_rows.append({
            'sequence_id': qid,
            'expected_domain': expected_domain,
            'expected_taxon': expected_taxon or 'NA',
            'observed_taxonomy': observed_taxonomy or 'NA',
            'status': status,
            'retained': 'YES' if retain else 'NO',
            'detail': detail,
        })

    report_path = Path(outdir) / 'expected_taxonomy_validation.tsv'
    with open(report_path, 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            'sequence_id', 'expected_domain', 'expected_taxon',
            'observed_taxonomy', 'status', 'retained', 'detail',
        ], delimiter='\t')
        writer.writeheader()
        writer.writerows(report_rows)
    logging.getLogger(__name__).info(
        '[PERFORMANCE REVIEW] Expected-versus-observed taxonomy: %d retained, %d rejected; report=%s',
        len(retained), len(report_rows) - len(retained), report_path,
    )
    return retained


def _normalise_sequence_domain(value: str | None) -> str | None:
    if value is None:
        return None
    got = normalise_domain_query(str(value))
    if got in ('', 'none', 'all', 'mixed'):
        return None
    return got


def _sequence_domain_filter(args, *, default: str | None = 'bacteria') -> str | None:
    """Return the domain filter implied by CLI flags."""
    profile = getattr(args, 'sequence_domain', None)
    if profile:
        return _normalise_sequence_domain(profile)
    return _normalise_sequence_domain(default)


def _sequence_domain_label(domain_filter: str | None) -> str:
    return domain_filter if domain_filter else 'mixed'


def _prune_dataset_by_kingdom(db: Database, dataset: str, kingdom: str | None, log_prefix: str):
    deleted = _prune_dataset_by_kingdom_helper(db, dataset, kingdom)
    if deleted:
        logging.getLogger(__name__).info(
            "%s Removed %d sequences from dataset %s not matching kingdom %s",
            log_prefix,
            deleted,
            dataset,
            kingdom,
        )
    return deleted


def _load_performance_review_baseline(
    args,
    db: Database,
    outdir: str,
    effective_ref: str | None,
    effective_taxa_tsv: str | None,
    threads: int,
):
    """Load an optional baseline FASTA before evaluating the current run."""
    baseline_fasta = getattr(args, 'baseline_fasta', None)
    if not baseline_fasta:
        return None

    baseline_dataset = getattr(args, 'baseline_dataset', None) or 'Baseline'
    run_dataset = getattr(args, 'dataset', None)
    if run_dataset and str(baseline_dataset) == str(run_dataset):
        raise SystemExit(
            "[BASELINE] --baseline-dataset must differ from --dataset so novelty can exclude the current run correctly."
        )

    log = logging.getLogger(__name__)
    baseline_tier = _normalise_cli_baseline_tier(
        baseline_dataset, getattr(args, 'baseline_tier', None),
    )
    db.upsert_dataset_role(
        baseline_dataset, 'baseline', genomes_available=True,
        baseline_tier=baseline_tier,
    )
    baseline_out = Path(outdir) / 'filing_cabinet_baseline'
    baseline_out.mkdir(parents=True, exist_ok=True)

    log.info(
        "[BASELINE] Loading baseline FASTA %s into dataset=%s before evaluating %s",
        baseline_fasta,
        baseline_dataset,
        run_dataset or '(current run)',
    )
    alias_entries, mapped_fasta = db.register_filing_cabinet(
        baseline_fasta,
        taxa_tsv=None,
        colour_csv=getattr(args, 'baseline_colours', None),
        source='baseline',
        dataset=baseline_dataset,
        outdir=str(baseline_out),
        shorten_ids=bool(getattr(args, 'baseline_shorten_ids', False)),
        baseline_tier=baseline_tier,
    )

    try:
        if alias_entries:
            map_path = _write_id_map_tsv(baseline_out / 'baseline_id_map.tsv', alias_entries)
            log.info("[BASELINE] Wrote baseline id mapping to %s", map_path)
    except Exception as e:
        log.warning("[BASELINE] Could not write baseline id mapping file: %s", e)

    orig_to_short = _build_orig_to_short(alias_entries)
    for short, _orig in alias_entries or []:
        orig_to_short[short] = short

    baseline_ref, baseline_taxa_tsv, baseline_assignment_tsv = _resolve_reference_inputs(
        effective_ref,
        effective_taxa_tsv,
        getattr(args, 'baseline_taxa_assignments', None),
        source_fasta_path=baseline_fasta,
        log_prefix='[BASELINE]',
    )

    if baseline_assignment_tsv:
        log.info("[BASELINE] Using baseline taxonomy assignments from %s", baseline_assignment_tsv)
        try:
            tax_entries = load_taxonomy_entries_from_assignments(
                baseline_assignment_tsv,
                orig_to_short,
                db,
                baseline_dataset,
                source_fasta_path=baseline_fasta,
            )
        except Exception as e:
            raise RuntimeError(f'Failed to read baseline taxonomy assignments {baseline_assignment_tsv}: {e}') from e
        if not tax_entries:
            raise RuntimeError(
                f'No baseline taxonomy rows mapped to {baseline_dataset}; check baseline FASTA/table IDs.'
            )
        if tax_entries:
            db.insert_taxonomy(tax_entries)
            log.info("[BASELINE] Inserted/updated taxonomy for %d baseline ids", len(tax_entries))
    elif not bool(getattr(args, 'baseline_skip_classify', False)):
        if not baseline_ref:
            log.warning("[BASELINE] No reference FASTA available; baseline loaded without taxonomy classification")
        else:
            input_for_classify = str(mapped_fasta) if mapped_fasta else baseline_fasta
            log.info("[BASELINE] Classifying baseline %s against %s", input_for_classify, baseline_ref)
            try:
                class_out = classify.run_classification(
                    input_for_classify,
                    str(baseline_out),
                    ref_fasta=baseline_ref,
                    taxa_tsv=baseline_taxa_tsv,
                    threads=threads,
                    db_name=getattr(args, 'ref_name', None) or _classify_derive_db_name(baseline_ref),
                )
                tax_entries, dist_entries = load_classification_results_for_dataset(
                    class_out,
                    orig_to_short,
                    db,
                    baseline_dataset,
                )
                if tax_entries:
                    db.insert_taxonomy(tax_entries)
                    log.info("[BASELINE] Inserted/updated taxonomy for %d baseline ids", len(tax_entries))
                if dist_entries:
                    db.insert_distances(dist_entries)
                    log.info("[BASELINE] Inserted/updated nearest-reference distances for %d baseline ids", len(dist_entries))
            except Exception as e:
                raise RuntimeError(f'Baseline classification failed: {e}') from e
    else:
        log.info("[BASELINE] --baseline-skip-classify set; baseline loaded without taxonomy classification")

    baseline_domain = _sequence_domain_filter(args)
    if baseline_domain:
        try:
            _prune_dataset_by_kingdom(db, baseline_dataset, baseline_domain, '[BASELINE]')
        except Exception as e:
            raise RuntimeError(f'Baseline domain filtering failed: {e}') from e

    try:
        combined_tax = baseline_out / 'baseline_combined_taxonomy.tsv'
        with db.connect() as conn:
            cur = conn.cursor()
            cur.execute(
                "SELECT s.id, t.taxonomy, t.confidence "
                "FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id "
                "WHERE s.dataset = ?",
                (baseline_dataset,),
            )
            rows = cur.fetchall()
        write_combined_taxonomy_tsv(combined_tax, rows)
        log.info("[BASELINE] Wrote baseline combined taxonomy for %d ids to %s", len(rows), combined_tax)
    except Exception as e:
        log.warning("[BASELINE] Failed to write baseline combined taxonomy: %s", e)

    return str(baseline_out)


def _same_path(path_a: str | None, path_b: str | None) -> bool:
    if not path_a or not path_b:
        return False
    try:
        return os.path.samefile(path_a, path_b)
    except Exception:
        try:
            return Path(path_a).resolve() == Path(path_b).resolve()
        except Exception:
            return str(path_a) == str(path_b)


def _build_alt_databases(args) -> list:
    """Return a list of (ref_fasta, taxa_tsv_or_None, db_name) for alt references."""
    alt_refs = getattr(args, 'alt_ref', None) or []
    alt_taxa = getattr(args, 'alt_taxa', None) or []
    alt_names = getattr(args, 'alt_ref_name', None) or []
    result = []
    for i, ref in enumerate(alt_refs):
        taxa = alt_taxa[i] if i < len(alt_taxa) else None
        raw_name = alt_names[i] if i < len(alt_names) else None
        name = raw_name if raw_name else _classify_derive_db_name(ref)
        result.append((ref, taxa, name))
    return result


def _store_alt_taxonomy_in_db(db: Database, all_results: dict, main_db_name: str):
    """Persist alt-db classification results to taxonomy_alt table.

    Iterates over *all_results* ``{db_name: {qid: (hit, pct, tax, conf)}}``,
    skipping the primary database (its results go into the standard taxonomy table).
    """
    log = logging.getLogger(__name__)
    for db_name, results in all_results.items():
        if db_name == main_db_name:
            continue
        alt_entries = []
        for qid, row in results.items():
            if not isinstance(row, (tuple, list)) or len(row) < 4:
                continue
            hit, pct, tax, conf = row[:4]
            if tax == 'NA' and hit == 'NA':
                continue  # skip unclassified rows to keep table lean
            alt_entries.append((qid, db_name, tax, conf, hit, pct))
        if alt_entries:
            try:
                db.insert_taxonomy_alt(alt_entries)
                log.info("[DB] Stored %d alt-db taxonomy entries for ref_db=%s", len(alt_entries), db_name)
            except Exception as e:
                log.warning("[DB] Failed to store alt-db taxonomy for %s: %s", db_name, e)
        evidence = []
        for qid, row in results.items():
            if not isinstance(row, (tuple, list)) or len(row) < 11:
                continue
            evidence.append({
                'sequence_id': qid, 'ref_db': db_name, 'best_hit': row[0], 'identity': row[1],
                'query_coverage': row[4], 'target_coverage': row[5], 'alignment_length': row[6],
                'query_length': row[7], 'target_length': row[8], 'mismatches': row[9], 'gaps': row[10],
                'source_classifier_confidence': row[3],
                'assignment_status': row[11] if len(row) > 11 else None,
                'taxonomy_support_pct': row[12] if len(row) > 12 else None,
                'taxonomy_candidate_count': row[16] if len(row) > 16 else None,
                'evidence_tier': row[15] if len(row) > 15 else None,
                'alignment_evidence_status': row[17] if len(row) > 17 else None,
                'assignment_evidence_level': row[18] if len(row) > 18 else None,
                'evidence_interpretation': row[19] if len(row) > 19 else None,
            })
        db.upsert_classification_evidence(evidence)


def _resolve_reference_inputs(
    ref_fasta: str | None,
    taxa_tsv: str | None,
    taxa_assignments: str | None,
    *,
    source_fasta_path: str | None,
    log_prefix: str,
):
    """Resolve whether --taxa-assignments is a TSV assignment file or a reference FASTA.

    External FASTA/FASTA.gz files provided via --taxa-assignments are treated as
    the effective classifier reference database. TSV inputs (and the special case
    where the assignments file is the same file as the source FASTA) are treated
    as direct input-sequence taxonomy assignments.
    """
    effective_ref = ref_fasta
    effective_taxa = taxa_tsv
    assignments_tsv = None

    if not taxa_assignments:
        return effective_ref, effective_taxa, assignments_tsv

    use_as_reference = False
    try:
        use_as_reference = (
            not _same_path(taxa_assignments, source_fasta_path)
            and _assignment_source_is_fasta(taxa_assignments, source_fasta_path=source_fasta_path)
        )
    except Exception:
        use_as_reference = False

    if use_as_reference:
        log = logging.getLogger(__name__)
        if effective_ref and not _same_path(effective_ref, taxa_assignments):
            log.info(
                "%s Treating --taxa-assignments=%s as the GTDB/reference FASTA and using it instead of --ref=%s (--ref is the preferred flag for reference databases)",
                log_prefix,
                taxa_assignments,
                effective_ref,
            )
        else:
            log.info(
                "%s Treating --taxa-assignments=%s as the GTDB/reference FASTA for classification (--ref is the preferred flag for reference databases)",
                log_prefix,
                taxa_assignments,
            )
        if not effective_taxa:
            log.info(
                "%s No --taxa TSV/CSV provided; taxonomy will be parsed directly from reference FASTA headers",
                log_prefix,
            )
        effective_ref = taxa_assignments
        return effective_ref, effective_taxa, None

    assignments_tsv = taxa_assignments
    return effective_ref, effective_taxa, assignments_tsv


def _normalise_cli_baseline_tier(dataset: str, explicit: str | None = None) -> str:
    value = str(explicit or '').strip().lower()
    if value in {'primary', 'hungate'}:
        value = 'priority'
    if value in {'priority', 'secondary'}:
        return value
    dataset_name = str(dataset or '').strip().lower()
    return 'priority' if 'hungate' in dataset_name else 'secondary'


def _resolve_novelty_baseline_datasets(db, args) -> list[str]:
    """Return every registered cultured baseline plus explicit current-run additions."""
    roles = db.get_dataset_roles()
    names = [
        name for name, details in roles.items()
        if details.get('role') == 'baseline'
    ]
    default_name = getattr(args, 'baseline_dataset', None)
    if default_name and (
        getattr(args, 'baseline_fasta', None)
        or roles.get(default_name, {}).get('role') == 'baseline'
    ):
        names.append(default_name)
    names.extend(getattr(args, 'novelty_baseline_datasets', None) or [])
    return list(dict.fromkeys(str(name) for name in names if str(name).strip()))


def _cmd_filing_cabinet_impl(args):
    db = Database(args.db)
    db.initialise()
    outdir = args.out or '.'
    threads = int(getattr(args, 'threads', 4) or 4)
    requested_kingdom = _sequence_domain_filter(args)
    os.makedirs(outdir, exist_ok=True)
    _configure_logging(outdir)
    logging.getLogger(__name__).info("[FILING CABINET] Starting baseline registration in %s", args.db)
    logging.getLogger(__name__).info(
        "[FILING CABINET] Sequence-domain profile: %s",
        _sequence_domain_label(requested_kingdom),
    )
    baseline_tier = _normalise_cli_baseline_tier(
        getattr(args, 'dataset', 'Baseline'), getattr(args, 'baseline_tier', None),
    )
    db.upsert_dataset_role(
        getattr(args, 'dataset', 'Baseline'),
        'baseline',
        genomes_available=True,
        baseline_tier=baseline_tier,
    )

    alias_entries, mapped_fasta = db.register_filing_cabinet(
        args.fasta,
        taxa_tsv=getattr(args, 'taxa', None),
        colour_csv=getattr(args, 'colours', None),
        source='filing_cabinet',
        dataset=getattr(args, 'dataset', 'Baseline'),
        outdir=outdir,
        shorten_ids=bool(getattr(args, 'shorten_ids', False)),
        baseline_tier=baseline_tier,
    )
    try:
        from branchmanager.marker_evidence import ensure_marker_evidence_schema, register_baseline_marker_records
        ensure_marker_evidence_schema(db)
        # The mapped FASTA contains the stable IDs actually inserted into the DB.
        if mapped_fasta and Path(mapped_fasta).is_file():
            evidence_count = register_baseline_marker_records(
                db, mapped_fasta, dataset=str(getattr(args, 'dataset', 'Baseline')),
            )
            logging.getLogger(__name__).info(
                '[FILING CABINET] Registered %d curated baseline marker-evidence record(s)', evidence_count,
            )
    except Exception as exc:
        raise RuntimeError(f'Could not register baseline marker evidence: {exc}') from exc
    try:
        if alias_entries:
            id_map_path = _write_id_map_tsv(Path(outdir) / 'filing_cabinet_id_map.tsv', alias_entries)
            logging.getLogger(__name__).info('[FILING CABINET] Wrote baseline ID mapping to %s', id_map_path)
    except Exception as e:
        logging.getLogger(__name__).warning('[FILING CABINET] Could not write baseline ID mapping file: %s', e)
    effective_ref, effective_taxa_tsv, assignment_tsv = _resolve_reference_inputs(
        getattr(args, 'ref', None),
        getattr(args, 'taxa', None),
        getattr(args, 'taxa_assignments', None),
        source_fasta_path=args.fasta,
        log_prefix='[FILING CABINET]',
    )
    classification_requested = bool(getattr(args, 'classify', False) or (getattr(args, 'taxa_assignments', None) and not assignment_tsv))

    # If the user provided a table of predetermined taxa assignments, use it
    # instead of running the classifier. The table should have at least ID and
    # taxonomy columns; TSV, CSV, and .gz variants are supported.
    if assignment_tsv:
        taxa_file = assignment_tsv
        logging.getLogger(__name__).info("[FILING CABINET] Using taxa assignments from %s (skipping classifier)", taxa_file)
        # build mapping orig->short from alias_entries
        orig_to_short = _build_orig_to_short(alias_entries)
        try:
            tax_entries = load_taxonomy_entries_from_assignments(
                taxa_file,
                orig_to_short,
                db,
                getattr(args, 'dataset', 'Baseline'),
                source_fasta_path=args.fasta,
            )
        except Exception as e:
            raise RuntimeError(f'Could not read Filing Cabinet taxonomy assignments {taxa_file}: {e}') from e

        if not tax_entries:
            raise RuntimeError(
                f'Filing Cabinet taxonomy assignments {taxa_file} did not map to any baseline sequence IDs'
            )
        db.insert_taxonomy(tax_entries)
        logging.getLogger(__name__).info("[FILING CABINET] Inserted/updated taxonomy for %d baseline IDs from taxa_assignments", len(tax_entries))

    # If classification requested, run classifier on the mapped fasta (short ids)
    # unless the user supplied a taxa assignments table, in which case use that
    # instead and skip running the external classifier.
    if classification_requested and not assignment_tsv:
        if not effective_ref:
            raise RuntimeError(
                'Filing Cabinet classification was requested but no reference FASTA was provided '
                'via --ref or --taxa-assignments'
            )
        else:
            input_for_classify = str(mapped_fasta) if mapped_fasta else args.fasta
            logging.getLogger(__name__).info("[FILING CABINET] Classifying baseline FASTA %s against %s", input_for_classify, effective_ref)

            alt_databases = _build_alt_databases(args)
            ref_name = getattr(args, 'ref_name', None) or _classify_derive_db_name(effective_ref)
            main_db = getattr(args, 'main_ref', None) or ref_name
            all_results: dict = {}

            if alt_databases:
                logging.getLogger(__name__).info(
                    "[FILING CABINET] Multi-database classification: primary=%s, alt=%s, main=%s",
                    ref_name, [n for _, _, n in alt_databases], main_db,
                )
                class_out, all_results = classify.run_all_classifications(
                    input_for_classify, outdir,
                    primary_ref=effective_ref,
                    primary_taxa=effective_taxa_tsv,
                    primary_name=ref_name,
                    alt_refs=alt_databases,
                    threads=threads,
                    main_db=main_db,
                )
                _store_alt_taxonomy_in_db(db, all_results, main_db)
            else:
                class_out = classify.run_classification(input_for_classify, outdir, ref_fasta=effective_ref, taxa_tsv=effective_taxa_tsv, threads=threads, db_name=ref_name)

            # parse classification output and persist taxonomy/distances only for
            # the registered baseline sequence IDs (short ids present in the DB)
            orig_to_short = _build_orig_to_short(alias_entries)
            try:
                tax_entries, dist_entries = load_classification_results_for_dataset(
                    class_out,
                    orig_to_short,
                    db,
                    getattr(args, 'dataset', 'Baseline'),
                )

                if tax_entries:
                    db.insert_taxonomy(tax_entries)
                    logging.getLogger(__name__).info("[FILING CABINET] Inserted/updated taxonomy for %d baseline IDs", len(tax_entries))
                if dist_entries:
                    db.insert_distances(dist_entries)
                    logging.getLogger(__name__).info("[FILING CABINET] Inserted/updated distances for %d baseline IDs", len(dist_entries))
            except Exception as e:
                raise RuntimeError(f'Could not persist Filing Cabinet classification output: {e}') from e

    # Apply the active sequence-domain profile before tree building so the
    # backbone tree is domain-specific.
    try:
        _prune_dataset_by_kingdom(
            db,
            getattr(args, 'dataset', 'Baseline'),
            requested_kingdom,
            '[FILING CABINET]',
        )
    except Exception as e:
        raise RuntimeError(f'Filing Cabinet domain filtering failed: {e}') from e

    # Optionally build a baseline tree/alignment from the registered baseline sequences
    if getattr(args, 'build_tree', False):
        try:
            # prefer mapped fasta (short ids) if present
            user_fasta = str(mapped_fasta) if mapped_fasta else args.fasta
            if requested_kingdom:
                try:
                    allowed_ids = _dataset_sequence_ids(db, getattr(args, 'dataset', 'Baseline'))
                    domain_fasta = Path(outdir) / f'filing_cabinet_{requested_kingdom}_sequences.fasta'
                    kept = _filter_fasta_to_ids(user_fasta, domain_fasta, allowed_ids)
                    if kept:
                        logging.getLogger(__name__).info(
                            "[FILING CABINET] Filtered baseline tree FASTA to %d %s sequence(s): %s",
                            kept,
                            requested_kingdom,
                            domain_fasta,
                        )
                        user_fasta = str(domain_fasta)
                    else:
                        raise RuntimeError(
                            f'Filing Cabinet domain filter {requested_kingdom} left no sequences for tree building'
                        )
                except Exception as e:
                    raise RuntimeError(f'Could not prepare domain-filtered Filing Cabinet tree FASTA: {e}') from e

            # optionally collapse registered baseline sequences before building tree
            if getattr(args, 'collapse', False):
                # require classification to be run for safe taxon-based collapsing
                if not (classification_requested or assignment_tsv):
                    raise RuntimeError('--collapse requires Filing Cabinet taxonomy assignments or --classify')
                else:
                    try:
                        threshold = float(getattr(args, 'collapse_threshold', 99.9))
                    except Exception:
                        threshold = 99.9

                    # build qid -> tax mapping from DB taxonomy for this dataset
                    qid_to_tax = {}
                    try:
                        with db.connect() as conn:
                            cur = conn.cursor()
                            cur.execute("SELECT s.id, t.taxonomy FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id WHERE s.dataset = ?", (getattr(args, 'dataset', 'Baseline'),))
                            for rid, tax in cur.fetchall():
                                qid_to_tax[rid] = tax
                    except Exception:
                        qid_to_tax = {}

                    # group records by tax and cluster within groups
                    from branchmanager.utils.fasta import read_fasta
                    taxa_groups = {}
                    try:
                        for h, s in read_fasta(user_fasta):
                            tax = qid_to_tax.get(h)
                            taxa_groups.setdefault(tax, []).append((h, s))
                    except Exception:
                        taxa_groups = {None: []}
                    try:
                        artefacts = collapse_fasta_within_taxa(
                            taxa_groups,
                            outdir,
                            'filing_cabinet_collapsed.fasta',
                            'filing_cabinet_collapsed_map.tsv',
                            'filing_cabinet_collapsed_members.tsv',
                            threshold=threshold,
                            threads=threads,
                            log_prefix='[FILING CABINET COLLAPSE]',
                            strict=True,
                        )
                        user_fasta = artefacts.collapsed_path
                        logging.getLogger(__name__).info("[FILING CABINET COLLAPSE] Wrote collapsed baseline FASTA %s (reps=%d)", user_fasta, len(artefacts.collapsed_records))
                    except Exception as e:
                        raise RuntimeError(f'Filing Cabinet collapse failed: {e}') from e

            logging.getLogger(__name__).info("[FILING CABINET] Building baseline tree/alignment in %s from registered baseline sequences", outdir)
            tree.initialise_or_update_tree(
                ref_fasta=effective_ref,
                user_fasta=user_fasta,
                outdir=outdir,
                db=None,
                threads=threads,
                anchor_file=getattr(args, 'anchors', None),
                tree_method=getattr(args, 'tree_method', 'fasttree'),
            )
            logging.getLogger(__name__).info("[FILING CABINET] Baseline tree/alignment written to %s", outdir)
            try:
                warning_rows = tree.collect_tree_build_warnings(user_fasta=str(user_fasta), anchor_file=getattr(args, 'anchors', None), db=None)
                warning_rows.extend(tree.summarise_alignment_quality(str(Path(outdir) / 'current_alignment.fasta')))
                if warning_rows:
                    warn_path = tree.write_tree_warning_tsv(outdir, warning_rows)
                    logging.getLogger(__name__).warning("[FILING CABINET] Tree/alignment warnings written to %s", warn_path)
            except Exception as e:
                logging.getLogger(__name__).warning("[FILING CABINET] Failed to summarise tree/alignment quality: %s", e)
        except Exception as e:
            raise RuntimeError(f'Filing Cabinet tree construction failed: {e}') from e

    # If a sequence-domain filter is requested for the Filing Cabinet, remove any registered baseline sequences
    # whose assigned taxonomy indicates they are not the requested kingdom.
    try:
        _prune_dataset_by_kingdom(
            db,
            getattr(args, 'dataset', 'Baseline'),
            requested_kingdom,
            '[FILING CABINET]',
        )
    except Exception as e:
        raise RuntimeError(f'Final Filing Cabinet domain validation failed: {e}') from e

    # Build combined taxonomy and generate iTOL colour files for the Filing Cabinet dataset.
    try:
        out_p = Path(outdir)
        combined_tax = out_p / 'filing_cabinet_combined_taxonomy.tsv'
        if not combined_tax.exists():
            try:
                with db.connect() as conn:
                    cur = conn.cursor()
                    cur.execute("SELECT s.id, t.taxonomy, t.confidence FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id WHERE s.dataset = ?", (getattr(args, 'dataset', 'Baseline'),))
                    rows = cur.fetchall()
                write_combined_taxonomy_tsv(combined_tax, rows)
                logging.getLogger(__name__).info("[FILING CABINET] Wrote combined taxonomy for %d ids to %s", len(rows), combined_tax)
            except Exception as e:
                logging.getLogger(__name__).warning("[FILING CABINET] Failed to build combined taxonomy: %s", e)

        # Build id_map from alias_entries so the generator can emit short ids
        id_map = {}
        try:
            if alias_entries:
                for short, orig in alias_entries:
                    id_map[orig] = short
                    id_map[short] = short
                    try:
                        cid = db._canonical_from_header(orig)
                        if cid:
                            id_map[cid] = short
                    except Exception:
                        pass  # could not resolve canonical ID for alias entry
        except Exception:
            id_map = {}

        # Call generator (fallback if id_map not accepted)
        try:
            tree_path = Path(outdir) / 'current_tree.nwk'
            tfile = str(tree_path) if tree_path.exists() else _find_tree_file_in_dir(outdir)
            itol.generate_itol_colours(str(combined_tax), outdir, user_colour_csv=getattr(args, 'colours', None), id_map=id_map, tree_file=tfile, phylum_groups=getattr(args, 'group_phyla', None))
            logging.getLogger(__name__).info("[FILING CABINET] Generated iTOL colour files in %s", outdir)
        except TypeError:
            try:
                itol.generate_itol_colours(str(combined_tax), outdir, user_colour_csv=getattr(args, 'colours', None), phylum_groups=getattr(args, 'group_phyla', None))
                logging.getLogger(__name__).info("[FILING CABINET] Generated iTOL colour files in %s (no id_map)", outdir)
            except Exception as e:
                logging.getLogger(__name__).warning("[FILING CABINET] Failed to generate iTOL colours: %s", e)
        except Exception as e:
            logging.getLogger(__name__).warning("[FILING CABINET] Failed to generate iTOL colours: %s", e)

        try:
            func_tsv = getattr(args, 'functional', None)
            if func_tsv:
                try:
                    written = itol.write_functional_annotations(str(func_tsv), outdir, id_map=id_map)
                    logging.getLogger(__name__).info("[FILING CABINET] Wrote functional annotation iTOL files: %s", ','.join(written) if written else '(none)')
                except Exception as e:
                    logging.getLogger(__name__).warning("[FILING CABINET] Functional annotations generation failed: %s", e)
        except Exception:
            pass  # functional annotation argument parsing failed

        if getattr(args, 'draft_rumen_functions', False):
            try:
                combined_tax_path = str(out_p / 'combined_taxonomy.tsv')
                tsv_out, itol_out = itol.generate_rumen_function_draft(
                    combined_tax_path, outdir, id_map=id_map
                )
                if tsv_out:
                    logging.getLogger(__name__).info(
                        "[FILING CABINET] Draft rumen functional annotation: %s", tsv_out
                    )
                if itol_out:
                    logging.getLogger(__name__).info(
                        "[FILING CABINET] Rumen functional iTOL file: %s", itol_out
                    )
            except Exception as e:
                logging.getLogger(__name__).warning(
                    "[FILING CABINET] Draft rumen functions generation failed: %s", e
                )

        try:
            ds_colour = itol._name_to_dataset_colour(getattr(args, 'dataset', 'Baseline'))
            itol_path = out_p / 'filing_cabinet_dataset.itol'
            dataset_label = getattr(args, 'dataset', 'Baseline')
            with db.connect() as conn:
                cur = conn.cursor()
                cur.execute("SELECT id FROM sequences WHERE dataset = ?", (getattr(args, 'dataset', 'Baseline'),))
                id_to_colour = {iid: ds_colour for (iid,) in cur.fetchall()}
            itol.write_dataset_colourstrip(str(itol_path), dataset_label, id_to_colour, legend_title=f"{dataset_label} legend")
            logging.getLogger(__name__).info("[FILING CABINET] Wrote dataset iTOL colour strip to %s", itol_path)
            tier_path = out_p / 'baseline_tier.itol'
            itol.write_baseline_tier_strip(str(tier_path), list(id_to_colour), {iid: baseline_tier for iid in id_to_colour})
            logging.getLogger(__name__).info("[FILING CABINET] Wrote baseline-tier iTOL colour strip to %s", tier_path)
        except Exception as e:
            logging.getLogger(__name__).warning("[FILING CABINET] Failed to write dataset iTOL colour strip: %s", e)
    except Exception:
        pass  # Filing Cabinet iTOL generation failed

    # write brief explanations for files produced by the Filing Cabinet
    try:
        _write_output_explanations(outdir)
    except Exception:
        pass  # could not write output explanations


def cmd_filing_cabinet(args):
    """Register the baseline Filing Cabinet without publishing a partial database."""
    import fcntl
    import sqlite3

    outdir = Path(args.out or '.').expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    original = Path(args.db).expanduser().resolve()
    original.parent.mkdir(parents=True, exist_ok=True)
    if getattr(args, 'build_tree', False):
        tree.preflight_tree_tools(getattr(args, 'tree_method', 'fasttree'))
    manifest = RunManifest(outdir, 'filing_cabinet')
    manifest.add_input(args.fasta, role='baseline_fasta')
    for role, source in (
        ('taxonomy', getattr(args, 'taxa', None)),
        ('taxonomy_assignments', getattr(args, 'taxa_assignments', None)),
        ('reference', getattr(args, 'ref', None)),
        ('colours', getattr(args, 'colours', None)),
        ('anchors', getattr(args, 'anchors', None)),
        ('functional_annotations', getattr(args, 'functional', None)),
    ):
        if source:
            manifest.add_input(source, role=role)
    for source in getattr(args, 'alt_ref', None) or []:
        manifest.add_input(source, role='alternate_reference')
    for source in getattr(args, 'alt_taxa', None) or []:
        if source:
            manifest.add_input(source, role='alternate_reference_taxonomy')

    staged_path = outdir / '.filing_cabinet_project.sqlite'
    lock_path = original.with_name(original.name + '.lock')
    with open(lock_path, 'a+') as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        try:
            staged_path.unlink(missing_ok=True)
            if original.is_file():
                manifest.add_input(original, role='project_database_before')
                with sqlite3.connect(str(original)) as source, sqlite3.connect(str(staged_path)) as target:
                    source.backup(target)
            else:
                Database(str(staged_path)).initialise()

            staged_args = argparse.Namespace(**vars(args))
            staged_args.db = str(staged_path)
            manifest.add_stage('filing_cabinet', 'RUNNING')
            _cmd_filing_cabinet_impl(staged_args)

            staged_db = Database(str(staged_path))
            staged_db.initialise()
            with staged_db.connect() as conn:
                integrity = conn.execute('PRAGMA integrity_check').fetchone()[0]
                baseline_count = conn.execute(
                    'SELECT COUNT(*) FROM sequences WHERE dataset = ?',
                    (getattr(args, 'dataset', 'Baseline'),),
                ).fetchone()[0]
            if integrity != 'ok':
                raise RuntimeError(f'staged Filing Cabinet database failed integrity_check: {integrity}')
            if baseline_count < 1:
                raise RuntimeError('Filing Cabinet registration produced no baseline sequences')

            for role, path, required in (
                ('id_map', outdir / 'filing_cabinet_id_map.tsv', False),
                ('combined_taxonomy', outdir / 'filing_cabinet_combined_taxonomy.tsv', True),
                ('alignment', outdir / 'current_alignment.fasta', bool(getattr(args, 'build_tree', False))),
                ('tree', outdir / 'current_tree.nwk', bool(getattr(args, 'build_tree', False))),
            ):
                manifest.add_output(path, role=role, required=required)
            missing = manifest.verify_required_outputs()
            if missing:
                raise RuntimeError('required Filing Cabinet outputs are missing: ' + ', '.join(missing))

            manifest.add_stage('filing_cabinet', 'COMPLETE', detail=f'{baseline_count} baseline sequences')
            manifest.finish('COMPLETE')
            staged_db.record_project_run(
                f'filing-cabinet:{manifest.data["started_at"]}', 'filing_cabinet', 'COMPLETE',
                dataset=getattr(args, 'dataset', 'Baseline'), manifest_path=str(manifest.json_path),
                started_at=manifest.data['started_at'], completed_at=manifest.data['completed_at'],
            )

            publish_tmp = original.with_name(original.name + '.publishing')
            publish_tmp.unlink(missing_ok=True)
            with sqlite3.connect(str(staged_path)) as source, sqlite3.connect(str(publish_tmp)) as target:
                source.backup(target)
                if target.execute('PRAGMA integrity_check').fetchone()[0] != 'ok':
                    raise RuntimeError('published Filing Cabinet database failed integrity_check')
            os.replace(publish_tmp, original)
            staged_path.unlink(missing_ok=True)
            manifest.add_output(original, role='project_database')
        except BaseException as exc:
            manifest.finish('FAILED', error=exc)
            raise
        finally:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
    print(
        f'[filing-cabinet] Registered {baseline_count} baseline sequence(s).\n'
        f'  Project database: {original}\n'
        f'  Backbone output : {outdir}'
    )


def cmd_filing_cabinet_configured_baselines(args):
    """Register every baseline in a profile, then make one combined backbone."""
    import fcntl
    import sqlite3
    import re
    from branchmanager.utils.fasta import write_fasta

    baselines = list(getattr(args, 'configured_baselines', []) or [])
    if not baselines:
        return cmd_filing_cabinet(args)
    outdir = Path(args.out or '.').expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    original = Path(args.db).expanduser().resolve()
    original.parent.mkdir(parents=True, exist_ok=True)
    if getattr(args, 'build_tree', False):
        tree.preflight_tree_tools(getattr(args, 'tree_method', 'fasttree'))
    manifest = RunManifest(outdir, 'filing_cabinet')
    for baseline in baselines:
        manifest.add_input(baseline['fasta'], role=f"baseline:{baseline.get('dataset', 'unnamed')}")
        for field, role in (('taxonomy_assignments', 'taxonomy_assignments'), ('taxa', 'taxonomy')):
            if baseline.get(field):
                manifest.add_input(baseline[field], role=f"{role}:{baseline.get('dataset', 'unnamed')}")
    if args.ref:
        manifest.add_input(args.ref, role='taxonomy_reference')
    if args.anchors:
        manifest.add_input(args.anchors, role='tree_anchors')
    staged = outdir / '.configured_baselines.sqlite'
    lock_path = original.with_name(original.name + '.lock')
    with open(lock_path, 'a+') as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        try:
            staged.unlink(missing_ok=True)
            if original.exists():
                with sqlite3.connect(str(original)) as source, sqlite3.connect(str(staged)) as target:
                    source.backup(target)
            else:
                Database(str(staged)).initialise()

            datasets = []
            for index, baseline in enumerate(baselines, start=1):
                required = {'dataset', 'fasta', 'tier'}
                missing = required - set(baseline)
                if missing:
                    raise ValueError(f"Baseline profile entry {index} is missing: {', '.join(sorted(missing))}")
                dataset = str(baseline['dataset'])
                datasets.append(dataset)
                slug = re.sub(r'[^A-Za-z0-9_.-]+', '_', dataset).strip('._') or f'baseline_{index}'
                child = argparse.Namespace(**vars(args))
                child.db = str(staged)
                child.fasta = baseline['fasta']
                child.dataset = dataset
                child.baseline_tier = baseline['tier']
                child.taxa = baseline.get('taxa')
                child.taxa_assignments = baseline.get('taxonomy_assignments')
                child.classify = bool(baseline.get('classify', False))
                child.build_tree = False
                child.out = str(outdir / 'baselines' / f'{index:02d}_{slug}')
                cmd_filing_cabinet(child)

            project_db = Database(str(staged))
            if getattr(args, 'build_tree', False):
                placeholders = ','.join('?' for _ in datasets)
                with project_db.connect() as conn:
                    rows = conn.execute(
                        f"SELECT id, sequence FROM sequences WHERE dataset IN ({placeholders}) "
                        "AND sequence IS NOT NULL AND sequence != '' ORDER BY dataset, id",
                        tuple(datasets),
                    ).fetchall()
                    taxonomy_rows = conn.execute(
                        f"SELECT s.id, t.taxonomy, t.confidence FROM sequences s "
                        f"LEFT JOIN taxonomy t ON s.id=t.id WHERE s.dataset IN ({placeholders}) "
                        "ORDER BY s.dataset, s.id",
                        tuple(datasets),
                    ).fetchall()
                if not rows:
                    raise RuntimeError('No marker sequences were available across configured baselines to build a tree')
                combined_fasta = outdir / 'configured_baselines.fasta'
                write_fasta(rows, str(combined_fasta))
                write_combined_taxonomy_tsv(outdir / 'filing_cabinet_combined_taxonomy.tsv', taxonomy_rows)
                tree.initialise_or_update_tree(
                    ref_fasta=getattr(args, 'ref', None), user_fasta=str(combined_fasta),
                    outdir=str(outdir), db=None, threads=int(getattr(args, 'threads', 4) or 4),
                    anchor_file=getattr(args, 'anchors', None),
                    tree_method=getattr(args, 'tree_method', 'fasttree'),
                )
                manifest.add_output(combined_fasta, role='combined_baseline_markers')
                manifest.add_output(outdir / 'current_alignment.fasta', role='alignment')
                manifest.add_output(outdir / 'current_tree.nwk', role='tree')

            with sqlite3.connect(str(staged)) as conn:
                integrity = conn.execute('PRAGMA integrity_check').fetchone()[0]
                total = conn.execute(
                    f"SELECT COUNT(*) FROM sequences WHERE dataset IN ({','.join('?' for _ in datasets)})",
                    tuple(datasets),
                ).fetchone()[0]
            if integrity != 'ok':
                raise RuntimeError(f'configured baseline database failed integrity_check: {integrity}')
            publish = original.with_name(original.name + '.publishing')
            publish.unlink(missing_ok=True)
            with sqlite3.connect(str(staged)) as source, sqlite3.connect(str(publish)) as target:
                source.backup(target)
            os.replace(publish, original)
            staged.unlink(missing_ok=True)
            manifest.add_stage('filing_cabinet', 'COMPLETE', detail=f'{total} records across {len(datasets)} baselines')
            manifest.finish('COMPLETE')
            manifest.add_output(original, role='project_database')
            print(f'[filing-cabinet] Registered {total} sequence record(s) across {len(datasets)} configured baselines.\n  Project database: {original}\n  Backbone output : {outdir}')
        except BaseException as exc:
            manifest.finish('FAILED', error=exc)
            raise
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


def _cmd_performance_review_impl(args):
    db = Database(args.db)
    db.initialise()
    # Marker evidence registry is additive to the legacy sequence/provenance tables.
    # Backfill once before the first v1 assessment so existing curated baselines and
    # accepted partner markers are classified under the explicit evidence policy.
    try:
        from branchmanager.marker_evidence import ensure_marker_evidence_schema, backfill_legacy_marker_evidence
        ensure_marker_evidence_schema(db)
        backfill_legacy_marker_evidence(db)
    except Exception as exc:
        raise RuntimeError(f'Could not initialise marker evidence registry: {exc}') from exc
    outdir = args.out
    threads = int(getattr(args, 'threads', 4) or 4)
    os.makedirs(outdir, exist_ok=True)
    _configure_logging(outdir)
    logging.getLogger(__name__).info("[PERFORMANCE REVIEW] Starting run pipeline (input=%s)", args.input)
    domain_filter = _sequence_domain_filter(args)
    logging.getLogger(__name__).info(
        "[PERFORMANCE REVIEW] Sequence-domain profile: %s",
        _sequence_domain_label(domain_filter),
    )
    effective_ref, effective_taxa_tsv, assignment_tsv = _resolve_reference_inputs(
        getattr(args, 'ref', None),
        getattr(args, 'taxa', None),
        getattr(args, 'taxa_assignments', None),
        source_fasta_path=args.input,
        log_prefix='[PERFORMANCE REVIEW]',
    )
    if not effective_ref:
        raise SystemExit("[PERFORMANCE REVIEW] A reference FASTA is required via --ref, or `--taxa-assignments` must point to a GTDB/reference FASTA rather than a taxonomy assignment table.")
    if getattr(args, 'command', None) == 'performance-review' and not getattr(args, 'partner_metadata', None):
        raise SystemExit(
            '[PERFORMANCE REVIEW] --partner-metadata / --sequencing-metadata is required: '
            'a cumulative CSV/TSV ledger with sequence IDs, partner acronyms, optional selection commitments, and already-sequenced status.'
        )

    _load_performance_review_baseline(args, db, outdir, effective_ref, effective_taxa_tsv, threads)
    if (
        hasattr(args, 'skip_chimera_check')
        and getattr(args, 'command', None) == 'performance-review'
    ):
        baseline_datasets = db.get_dataset_names_by_role('baseline')
        if not baseline_datasets:
            raise SystemExit(
                '[PERFORMANCE REVIEW] No cultured baseline is registered. Supply --baseline-fasta for this run '
                'or establish it first with `branchmanager filing-cabinet`.'
            )

    # QC. An explicit APPROVE_FAIL_QC decision must survive this generic
    # length/N prefilter; otherwise a reviewer-approved short/high-N sequence
    # would be silently removed before the marker-provenance workflow can use it.
    manual_approved_ids = set()
    review_path = getattr(args, 'marker_review', None)
    if review_path:
        from branchmanager.marker_provenance import _review_decisions
        reviewed = _review_decisions(review_path)
        manual_approved_ids = {
            sequence_id for sequence_id, review in reviewed.items()
            if str(review.get('decision') or '').upper() == 'APPROVE_FAIL_QC'
        }
    qc_out = qc.run_qc(
        args.input,
        outdir,
        min_len=getattr(args, 'min_len', 800),
        max_n_percent=getattr(args, 'max_n_percent', 5.0),
        manual_approved_ids=manual_approved_ids,
    )

    # Preserve every isolate ID, including exact marker-sequence duplicates.
    # Distinct isolates can still contain different genomes or provide useful
    # DNA-extraction backups; optional tree collapsing handles visual density.
    candidate_input = qc_out

    # Map user-provided IDs to runtime IDs and insert them into the DB.
    from branchmanager.utils.fasta import read_fasta, write_fasta
    current_dataset_fasta = Path(outdir) / 'current_dataset_sequences.fasta'
    used_ids = set(db.get_all_ids())
    orig_to_short = {}
    mapped_records = []
    initial_marker_records = []
    pending_initial_markers = {}
    replaced_marker_records = []
    skipped_existing = 0
    resequencing_decisions = []
    replacement_report_rows = []
    provenance_for_replacements = {}
    declined_resequencing_ids = set()
    failed_resequencing_ids = set()

    def _incoming_ab1_files(sequence_id):
        """Find raw-read paths recorded alongside the incoming marker FASTA."""
        from branchmanager.marker_provenance import _read_table
        files = set()
        input_parent = Path(args.input).expanduser().resolve().parent
        for read_qc in (input_parent / 'read_qc.tsv', Path(outdir) / 'read_qc.tsv'):
            if not read_qc.is_file():
                continue
            for row in _read_table(read_qc):
                row_id = row.get('sequenceid') or row.get('sequence_id') or ''
                source = row.get('sourcefile') or row.get('source_file') or ''
                if row_id == sequence_id and source and Path(source).suffix.lower() in {'.ab1', '.abi'}:
                    files.add(str(Path(source).expanduser()))
        return sorted(files)

    def _write_resequencing_decisions():
        if not resequencing_decisions:
            return
        decision_path = Path(outdir) / 'resequencing_conflicts.tsv'
        with open(decision_path, 'w', newline='') as handle:
            writer = csv.DictWriter(
                handle, fieldnames=[
                    'sequence_id', 'incoming_ab1_files', 'decision', 'qc_class', 'qc_reasons',
                    'qc_review_decision', 'detail',
                ], delimiter='\t',
            )
            writer.writeheader(); writer.writerows(resequencing_decisions)

    # If the user requested a kingdom filter, classify the QC-passing input
    # before insertion so only sequences assigned to that domain are retained.
    early_class_out = None
    allowed_qids = None
    kingdom = domain_filter
    if kingdom:
        kingdom_text = str(kingdom)
        if not effective_ref:
            logging.getLogger(__name__).warning("[PERFORMANCE REVIEW] --sequence-domain specified but no reference FASTA was available; cannot classify to filter; proceeding without domain filtering")
            allowed_qids = None
        else:
            try:
                logging.getLogger(__name__).info("[PERFORMANCE REVIEW] Running pre-insert classification to filter by domain=%s", kingdom)
                early_class_out = classify.run_classification(
                    str(candidate_input), outdir, ref_fasta=effective_ref,
                    taxa_tsv=effective_taxa_tsv, threads=threads,
                    db_name=getattr(args, 'ref_name', None) or _classify_derive_db_name(effective_ref),
                )
                if getattr(args, 'partner_metadata', None):
                    allowed_qids = _validate_partner_taxonomy(
                        early_class_out, args.partner_metadata, kingdom_text, outdir,
                        input_fasta=str(candidate_input),
                    )
                else:
                    allowed_qids = _classification_ids_matching_kingdom(early_class_out, kingdom_text)
                logging.getLogger(__name__).info("[PERFORMANCE REVIEW] Domain/taxon validation: %d QC-passing sequences retained for %s", len(allowed_qids), kingdom)
            except Exception as e:
                raise SystemExit(f'[PERFORMANCE REVIEW] Domain classification/filtering failed: {e}')

    for h, s in read_fasta(candidate_input):
        # if kingdom filtering is active, skip sequences that did not match
        if allowed_qids is not None:
            # try a few candidate header forms for matching
            candidates = [h, h.split('|')[-1], h.split()[-1] if h.split() else h]
            hit = False
            for c in candidates:
                if c in allowed_qids:
                    hit = True
                    break
            if not hit:
                continue

        # If DB already contains this sequence ID (or another alias), we still
        # want to include it in the tree for visualisation, but we use the
        # existing DB ID so downstream taxonomy/metadata resolve consistently.
        # The "INSERT OR IGNORE" during insert_sequences will prevent DB duplicates.
        existing_id = db.resolve_sequence_id(h)
        if existing_id and existing_id in used_ids:
            with db.connect() as conn:
                stored_row = conn.execute(
                    'SELECT sequence FROM sequences WHERE id=?', (existing_id,),
                ).fetchone()
            if stored_row and not str(stored_row[0] or '').strip():
                previous_pending = pending_initial_markers.get(existing_id)
                if previous_pending is not None and previous_pending != s:
                    raise SystemExit(
                        f'[PERFORMANCE REVIEW] Multiple different marker sequences were submitted for '
                        f'markerless isolate {existing_id!r} in one run.'
                    )
                if previous_pending is None:
                    pending_initial_markers[existing_id] = s
                    initial_marker_records.append((existing_id, s))
                    mapped_records.append((existing_id, s))
                orig_to_short[h] = existing_id
                try:
                    canonical_id = canonicalise_sequence_id(h)
                    if canonical_id:
                        orig_to_short[canonical_id] = existing_id
                except Exception:
                    pass  # could not canonicalise sequence ID for orig_to_short
                continue
            try:
                db.assert_sequence_compatible(existing_id, s)
            except ValueError as e:
                incoming_files = _incoming_ab1_files(existing_id)
                print(f'[DB] {e}')
                if incoming_files:
                    print(f'[DB] Incoming AB1 files for {existing_id}:')
                    for incoming_file in incoming_files:
                        print(f'  {incoming_file}')
                else:
                    print(f'[DB] No AB1 paths were found in read_qc.tsv for {existing_id}.')
                try:
                    response = input(f'Treat this as a resequenced marker and re-examine isolate {existing_id}? [y/N] ').strip().lower()
                except EOFError:
                    response = ''
                approved = response in {'y', 'yes'}
                decision_row = {
                    'sequence_id': existing_id,
                    'incoming_ab1_files': ';'.join(incoming_files),
                    'decision': 'YES' if approved else 'NO',
                    'qc_class': '', 'qc_reasons': '', 'qc_review_decision': '',
                    'detail': str(e),
                }
                resequencing_decisions.append(decision_row)
                _write_resequencing_decisions()
                provenance_row = {}
                if approved:
                    from branchmanager.marker_provenance import load_marker_provenance
                    if existing_id not in provenance_for_replacements:
                        replacement_provenance, _qc_source = load_marker_provenance(
                            args.input,
                            qc_path=getattr(args, 'marker_qc', None),
                            review_path=getattr(args, 'marker_review', None),
                            id_map={h: existing_id},
                            accept_unverified=bool(getattr(args, 'accept_unverified_marker_qc', False)),
                            unverified_override_authorizer=getattr(args, 'unverified_override_authorizer', ''),
                            unverified_override_reason=getattr(args, 'unverified_override_reason', ''),
                        )
                        provenance_for_replacements[existing_id] = next(
                            (row for row in replacement_provenance if row.get('sequence_id') == existing_id),
                            {'sequence_id': existing_id},
                        )
                    provenance_row = provenance_for_replacements[existing_id]
                    qc_class = str(provenance_row.get('marker_qc_class') or 'QUALITY_UNVERIFIED').upper()
                    qc_reasons = str(provenance_row.get('marker_qc_reasons') or 'no_marker_qc_sidecar')
                    decision_row['qc_class'] = qc_class
                    decision_row['qc_reasons'] = qc_reasons
                    previous_provenance = db.get_sequence_provenance([existing_id]).get(existing_id, {})
                    previous_qc_class = str(
                        previous_provenance.get('marker_qc_class') or '',
                    ).upper()
                    if (
                        len(s) >= 800
                        and previous_qc_class in _MARKER_QC_RANK
                        and qc_class in _MARKER_QC_RANK
                        and _MARKER_QC_RANK[qc_class] < _MARKER_QC_RANK[previous_qc_class]
                    ):
                        print(
                            f'[DB] Incoming marker for {existing_id} has lower QC '
                            f'({qc_class}) than the stored marker ({previous_qc_class}).'
                        )
                        try:
                            qc_response = input(
                                f'Use the lower-QC incoming marker for {existing_id}? [y/N] ',
                            ).strip().lower()
                        except EOFError:
                            qc_response = ''
                        approved = qc_response in {'y', 'yes'}
                        decision_row['qc_review_decision'] = (
                            'USER_APPROVED_LOWER_QC_REPLACEMENT'
                            if approved else 'USER_DECLINED_LOWER_QC_REPLACEMENT'
                        )
                        if not approved:
                            decision_row['detail'] = (
                                f'{decision_row["detail"]}; lower-QC replacement declined '
                                f'({qc_class} versus {previous_qc_class})'
                            )
                            db.record_marker_replacement_attempt(
                                existing_id, 'LOWER_QC_DECLINED',
                                f'user retained existing {previous_qc_class} marker over incoming {qc_class} marker',
                                source_manifest=str(Path(outdir) / 'run_manifest.json'),
                                source_fasta=str(args.input),
                                incoming_qc_class=qc_class,
                                incoming_qc_reasons=qc_reasons,
                                incoming_sequence=s,
                            )
                    else:
                        decision_row['qc_review_decision'] = 'USER_APPROVED_RESEQUENCING'
                    _write_resequencing_decisions()
                if approved:
                    replacement = {**provenance_row, 'sequence_id': existing_id, 'sequence': s}
                    result = db.replace_selected_markers(
                        [replacement], source_manifest=str(Path(outdir) / 'run_manifest.json'),
                        source_fasta=str(args.input), allow_unselected=True,
                    )[0]
                    replacement_report_rows.append(
                        _marker_replacement_report_row(existing_id, result)
                    )
                    if result.get('status') != 'REPLACED':
                        detail = result.get('detail') or result.get('status')
                        decision_row['detail'] = f'{decision_row["detail"]}; replacement not applied: {detail}'
                        _write_resequencing_decisions()
                        if result.get('status') != 'TOO_SHORT':
                            raise SystemExit(
                                f'[PERFORMANCE REVIEW] Approved resequencing for {existing_id!r} could not be applied '
                                f'({result.get("status")}: {detail}); inspect resequencing_conflicts.tsv and resolve '
                                'the replacement eligibility issue before rerunning the review.'
                            )
                        with db.connect() as conn:
                            existing_row = conn.execute(
                                'SELECT sequence FROM sequences WHERE id=?', (existing_id,),
                            ).fetchone()
                        if not existing_row or not existing_row[0]:
                            raise SystemExit(
                                f'[PERFORMANCE REVIEW] Replacement for {existing_id!r} failed QC and no stored marker '
                                'is available to retain.'
                            )
                        mapped_records.append((existing_id, existing_row[0]))
                        if result.get('status') == 'TOO_SHORT':
                            failed_resequencing_ids.add(existing_id)
                            print(
                                f'[DB] Replacement for {existing_id} did not meet the 800 bp minimum ({detail}); '
                                'retaining and flagging the old marker.'
                            )
                    else:
                        provenance_for_replacements[existing_id] = replacement
                        db.upsert_sequence_provenance([replacement])
                        db.update_isolate_status(
                            existing_id, 'MARKER_REPLACED', detail=result.get('detail', ''),
                            source_file=str(args.input),
                        )
                        mapped_records.append((existing_id, s))
                        replaced_marker_records.append((existing_id, s))
                        short = existing_id
                    _write_resequencing_decisions()
                    short = existing_id
                else:
                    short = existing_id
                    declined_resequencing_ids.add(existing_id)
                    with db.connect() as conn:
                        existing_row = conn.execute(
                            'SELECT sequence FROM sequences WHERE id=?', (existing_id,),
                        ).fetchone()
                    previous_provenance = db.get_sequence_provenance([existing_id]).get(existing_id, {})
                    existing = existing_row[0] if existing_row else None
                    if existing:
                        mapped_records.append((short, existing))
                    else:
                        raise SystemExit(f'[DB] Existing sequence {existing_id!r} could not be loaded after declining resequencing.')
                    lower_qc_declined = (
                        decision_row.get('qc_review_decision') == 'USER_DECLINED_LOWER_QC_REPLACEMENT'
                    )
                    replacement_report_rows.append(_marker_replacement_report_row(
                        existing_id,
                        {
                            'status': 'LOWER_QC_DECLINED' if lower_qc_declined else 'USER_DECLINED',
                            'detail': (
                                f'user retained existing {previous_provenance.get("marker_qc_class") or "QC"} '
                                f'marker over incoming {provenance_row.get("marker_qc_class") or "lower-QC"} marker'
                                if lower_qc_declined
                                else 'user declined replacement; existing marker retained'
                            ),
                            'retained_marker': 'EXISTING',
                            'previous_length': len(existing),
                            'previous_qc_class': previous_provenance.get('marker_qc_class') or '',
                            'previous_qc_reasons': previous_provenance.get('marker_qc_reasons') or '',
                            'previous_read_ids': previous_provenance.get('source_read_ids') or '',
                            'previous_read_files': previous_provenance.get('source_read_files') or '',
                            'incoming_length': len(s) if lower_qc_declined else '',
                            'incoming_qc_class': provenance_row.get('marker_qc_class') or '',
                            'incoming_qc_reasons': provenance_row.get('marker_qc_reasons') or '',
                            'incoming_read_ids': provenance_row.get('source_read_ids') or '',
                            'incoming_read_files': ';'.join(incoming_files),
                        },
                    ))
                orig_to_short[h] = short
                try:
                    cid = canonicalise_sequence_id(h)
                    if cid:
                        orig_to_short[cid] = short
                except Exception:
                    pass  # could not canonicalise sequence ID for orig_to_short
                skipped_existing += 1
                continue
            incoming_files = _incoming_ab1_files(existing_id)
            if incoming_files:
                from branchmanager.marker_provenance import load_marker_provenance
                incoming_rows, _ = load_marker_provenance(
                    args.input,
                    qc_path=getattr(args, 'marker_qc', None),
                    review_path=getattr(args, 'marker_review', None),
                    id_map={h: existing_id},
                    accept_unverified=bool(getattr(args, 'accept_unverified_marker_qc', False)),
                    unverified_override_authorizer=getattr(args, 'unverified_override_authorizer', ''),
                    unverified_override_reason=getattr(args, 'unverified_override_reason', ''),
                )
                incoming = next(
                    (row for row in incoming_rows if row.get('sequence_id') == existing_id),
                    {'sequence_id': existing_id},
                )
                previous = db.get_sequence_provenance([existing_id]).get(existing_id, {})
                merged_provenance = _merge_identical_marker_provenance(previous, incoming)
                db.upsert_sequence_provenance([merged_provenance])
                db.record_marker_replacement_attempt(
                    existing_id, 'SEQUENCE_UNCHANGED',
                    'incoming sequence is identical to the stored marker; existing sequence retained',
                    source_manifest=str(incoming.get('source_manifest') or ''),
                    source_fasta=str(incoming.get('source_sequence_file') or args.input),
                    incoming_qc_class=incoming.get('marker_qc_class') or '',
                    incoming_qc_reasons=incoming.get('marker_qc_reasons') or '',
                    incoming_sequence=s,
                )
                replacement_report_rows.append(_marker_replacement_report_row(
                    existing_id,
                    {
                        'status': 'SEQUENCE_UNCHANGED',
                        'detail': 'incoming sequence is identical to the stored marker; incoming read provenance recorded',
                        'retained_marker': 'EXISTING',
                        'comparison': 'SEQUENCE_IDENTICAL',
                        'previous_length': len(s),
                        'incoming_length': len(s),
                        'previous_qc_class': previous.get('marker_qc_class') or '',
                        'incoming_qc_class': incoming.get('marker_qc_class') or '',
                        'previous_qc_reasons': previous.get('marker_qc_reasons') or '',
                        'incoming_qc_reasons': incoming.get('marker_qc_reasons') or '',
                        'previous_read_ids': previous.get('source_read_ids') or '',
                        'incoming_read_ids': incoming.get('source_read_ids') or '',
                        'previous_used_read_ids': previous.get('source_used_read_ids') or '',
                        'incoming_used_read_ids': incoming.get('source_used_read_ids') or '',
                        'previous_read_count': _used_marker_read_count(
                            previous.get('source_used_read_ids'),
                            previous.get('source_read_ids'),
                            previous.get('marker_qc_reasons'),
                        ),
                        'incoming_read_count': _used_marker_read_count(
                            incoming.get('source_used_read_ids'),
                            incoming.get('source_read_ids'),
                            incoming.get('marker_qc_reasons'),
                        ),
                        'previous_read_files': previous.get('source_read_files') or '',
                        'incoming_read_files': ';'.join(incoming_files),
                        'previous_sequence_sha256': previous.get('source_sequence_sha256') or '',
                        'incoming_sequence_sha256': incoming.get('source_sequence_sha256') or '',
                    },
                ))
            skipped_existing += 1
            short = existing_id
            orig_to_short[h] = short
            try:
                cid = canonicalise_sequence_id(h)
                if cid:
                    orig_to_short[cid] = short
            except Exception:
                pass  # could not canonicalise sequence ID for orig_to_short
            # Still add to mapped_records so it appears in the tree
            mapped_records.append((short, s))
        else:
            try:
                short = db.choose_effective_sequence_id(
                    h,
                    used_ids,
                    shorten_ids=bool(getattr(args, 'shorten_ids', False)),
                )
            except ValueError as e:
                raise SystemExit(f"[PERFORMANCE REVIEW] {e}")
            mapped_records.append((short, s))
            orig_to_short[h] = short
            try:
                cid = canonicalise_sequence_id(h)
                if cid:
                    orig_to_short[cid] = short
            except Exception:
                pass  # could not canonicalise sequence ID for orig_to_short
    if skipped_existing:
        logging.getLogger(__name__).info("[DB] Found %d sequences already in DB; keeping them for tree inclusion", skipped_existing)
    _write_marker_replacement_report(
        Path(outdir) / 'marker_replacement_report.tsv',
        replacement_report_rows,
    )
    write_fasta(
        replaced_marker_records,
        str(Path(outdir) / 'replaced_markers.fasta'),
    )
    write_fasta(mapped_records, str(current_dataset_fasta))
    if resequencing_decisions:
        decision_path = Path(outdir) / 'resequencing_conflicts.tsv'
        _write_resequencing_decisions()
        logging.getLogger(__name__).info('[DB] Resequencing decisions written to %s', decision_path)
    logging.getLogger(__name__).info("[DB] Mapped %d user sequence IDs to runtime IDs and wrote %s", len(mapped_records), current_dataset_fasta)
    if not mapped_records:
        logging.getLogger(__name__).warning(
            '[PERFORMANCE REVIEW] No marker sequences passed QC/domain filtering; failed isolates will be '
            'registered as markerless QC-failed records with no classifier taxonomy.'
        )
    # insert mapped records into DB (dataset provided by user)
    run_dataset = getattr(args, 'dataset', 'user')
    db.upsert_dataset_role(run_dataset, 'candidate', genomes_available=False)
    db.insert_sequences(mapped_records, dataset=run_dataset)
    for sequence_id, sequence in initial_marker_records:
        db.attach_initial_marker_sequence(sequence_id, sequence)
        replacement_report_rows.append(_marker_replacement_report_row(
            sequence_id,
            {
                'status': 'INITIAL_MARKER_ADDED',
                'retained_marker': 'INCOMING',
                'previous_length': 0,
                'incoming_length': len(sequence),
                'detail': 'first passing marker sequence attached to the existing markerless isolate ID',
            },
        ))
    # write mapping file
    try:
        map_path = Path(outdir) / 'user_id_map.tsv'
        _write_id_map_tsv(map_path, ((short, orig) for orig, short in orig_to_short.items()))
        logging.getLogger(__name__).info("[DB] Wrote user id mapping to %s", map_path)
    except Exception as e:
        logging.getLogger(__name__).warning("[DB] Could not write user id mapping file: %s", e)

    for short, _seq in mapped_records:
        orig_to_short[short] = short

    run_metadata = _load_partner_metadata_for_run(
        args,
        db,
        outdir,
        orig_to_short,
        [short for short, _seq in mapped_records],
    )

    # Carry Paper Trail/Merge Meeting quality and raw-read provenance into the persistent
    # project record. Direct FASTA inputs are intentionally QUALITY_UNVERIFIED
    # unless the user supplies an explicit acceptance or review decision.
    try:
        from branchmanager.marker_provenance import load_marker_provenance, write_marker_qc_bridge
        marker_provenance, marker_qc_source = load_marker_provenance(
            args.input,
            qc_path=getattr(args, 'marker_qc', None),
            review_path=getattr(args, 'marker_review', None),
            id_map=orig_to_short,
            accept_unverified=bool(getattr(args, 'accept_unverified_marker_qc', False)),
            unverified_override_authorizer=getattr(args, 'unverified_override_authorizer', ''),
            unverified_override_reason=getattr(args, 'unverified_override_reason', ''),
        )
        mapped_id_set = {short for short, _sequence in mapped_records}
        marker_provenance = [
            row for row in marker_provenance if row.get('sequence_id') in mapped_id_set
        ]
        if declined_resequencing_ids:
            prior_provenance = db.get_sequence_provenance(declined_resequencing_ids)
            marker_provenance = [
                prior_provenance[row['sequence_id']]
                if row['sequence_id'] in declined_resequencing_ids and row['sequence_id'] in prior_provenance
                else row
                for row in marker_provenance
                if row['sequence_id'] not in declined_resequencing_ids or row['sequence_id'] in prior_provenance
            ]
        if provenance_for_replacements:
            marker_provenance = [
                provenance_for_replacements.get(row['sequence_id'], row)
                for row in marker_provenance
            ]
        if failed_resequencing_ids:
            prior_provenance = db.get_sequence_provenance(failed_resequencing_ids)
            marker_provenance = [
                prior_provenance[row['sequence_id']]
                if row['sequence_id'] in prior_provenance
                else row
                for row in marker_provenance
                if row['sequence_id'] not in failed_resequencing_ids or row['sequence_id'] in prior_provenance
            ]
        db.upsert_sequence_provenance(marker_provenance)
        marker_bridge_path = write_marker_qc_bridge(
            Path(outdir) / 'marker_qc_provenance.tsv', marker_provenance,
        )
        existing_statuses = db.get_isolate_statuses(mapped_id_set)
        terminal_states = {
            'SAB_APPROVED', 'DNA_EXTRACTION_PENDING', 'DNA_EXTRACTION_FAILED',
            'LIBRARY_PENDING', 'SEQUENCED', 'GENOME_QC_FAILED', 'GENOME_QC_PASSED', 'WITHDRAWN',
        }
        for row in marker_provenance:
            sequence_id = row['sequence_id']
            if existing_statuses.get(sequence_id, {}).get('status') in terminal_states:
                continue
            qc_class = str(row.get('marker_qc_class') or '').upper()
            review = str(row.get('manual_review_status') or '').upper()
            status = 'MARKER_QC_PASSED' if qc_class == 'PASS_HIGH_CONFIDENCE' or review == 'APPROVED' else 'TRACE_REVIEW'
            db.update_isolate_status(
                sequence_id, status,
                detail=str(row.get('marker_qc_reasons') or ''),
                source_file=marker_qc_source or str(Path(args.input).resolve()),
            )
        logging.getLogger(__name__).info(
            '[PERFORMANCE REVIEW] Marker provenance for %d sequence(s) written to %s',
            len(marker_provenance), marker_bridge_path,
        )
    except Exception as exc:
        raise SystemExit(f'[PERFORMANCE REVIEW] Marker-QC provenance could not be established: {exc}')

    try:
        if (
            not mapped_records
            or not hasattr(args, 'skip_chimera_check')
            or getattr(args, 'skip_chimera_check', False)
        ):
            chimera_results = {
                sequence_id: {'call': 'SKIPPED', 'score': None}
                for sequence_id, _sequence in mapped_records
            }
            logging.getLogger(__name__).warning(
                '[PERFORMANCE REVIEW] Chimera screening explicitly skipped; affected marker evidence will require review.'
            )
        else:
            from branchmanager.pipeline import chimera as _chimera
            chimera_report, chimera_results = _chimera.run_reference_screen(
                str(current_dataset_fasta),
                getattr(args, 'chimera_ref', None) or effective_ref,
                outdir,
                threads=threads,
            )
            logging.getLogger(__name__).info('[PERFORMANCE REVIEW] Chimera report written to %s', chimera_report)
        for row in marker_provenance:
            result = chimera_results.get(row['sequence_id'], {'call': 'INDETERMINATE', 'score': None})
            row['chimera_call'] = result.get('call', 'INDETERMINATE')
            row['chimera_score'] = result.get('score')
        db.upsert_sequence_provenance(marker_provenance)
        write_marker_qc_bridge(Path(outdir) / 'marker_qc_provenance.tsv', marker_provenance)
        # Register the exact sequence currently stored in the project DB, not a
        # declined incoming replacement. This is what future novelty runs compare.
        from branchmanager.marker_evidence import upsert_active_marker_evidence
        source_type = str(getattr(args, 'marker_source_type', 'auto') or 'auto').upper()
        if source_type == 'AUTO':
            source_type = 'PARTNER_FASTA'
        preserve_ids = set(declined_resequencing_ids) | set(failed_resequencing_ids)
        evidence_records = []
        for sequence_id, _incoming_sequence in mapped_records:
            if sequence_id in preserve_ids:
                continue
            with db.connect() as conn:
                stored = conn.execute('SELECT sequence FROM sequences WHERE id=?', (sequence_id,)).fetchone()
            if stored and str(stored[0] or '').strip():
                evidence_records.append((sequence_id, str(stored[0])))
        # A cohort FASTA may mix source types; a per-row source in the QC sidecar wins.
        row_source = {
            row['sequence_id']: str(row.get('marker_source_type') or '').upper()
            for row in marker_provenance
        }
        evidence_groups: dict[str, list] = {}
        for record in evidence_records:
            evidence_groups.setdefault(row_source.get(record[0]) or source_type, []).append(record)
        inserted_evidence = []
        for group_source, group_records in evidence_groups.items():
            inserted_evidence.extend(upsert_active_marker_evidence(
                db, group_records, marker_provenance,
                source_type=group_source, dataset=run_dataset, source_file=str(Path(args.input).resolve()),
            ) or [])
        logging.getLogger(__name__).info(
            '[PERFORMANCE REVIEW] Registered %d marker evidence record(s) in the v1 evidence registry',
            len(inserted_evidence),
        )
    except Exception as exc:
        raise SystemExit(f'[PERFORMANCE REVIEW] Reference-based chimera screening failed: {exc}')

    # Resolve alternate references once for every control-flow path. Novelty
    # metrics consume this list even when classification is skipped because
    # no sequences survived QC/domain filtering or assignments were supplied.
    alt_databases = _build_alt_databases(args)

    # classify (or use external taxa assignments if provided)
    class_out = None
    all_class_results: dict = {}
    run_main_db_name: str = (
        getattr(args, 'main_ref', None)
        or getattr(args, 'ref_name', None)
        or _classify_derive_db_name(effective_ref)
        or 'main'
    )
    if not mapped_records:
        class_out = str(Path(outdir) / 'taxonomy.tsv')
        with open(class_out, 'w') as handle:
            handle.write('\t'.join(classify._classification_header(run_main_db_name)) + '\n')
        empty_databases = {run_main_db_name: {}}
        for _ref, _taxa, name in alt_databases:
            empty_databases.setdefault(name, {})
        primary_name = getattr(args, 'ref_name', None) or _classify_derive_db_name(effective_ref)
        empty_databases.setdefault(primary_name, {})
        if len(empty_databases) > 1:
            classify._write_merged_taxonomy(
                outdir, empty_databases, run_main_db_name, [],
            )
    elif assignment_tsv:
        taxa_file = assignment_tsv
        logging.getLogger(__name__).info("[PERFORMANCE REVIEW] Using taxa assignments from %s (skipping classifier)", taxa_file)
        # parse provided taxonomy assignment table and insert taxonomy for mapped runtime ids
        try:
            tax_entries_local = load_taxonomy_entries_from_assignments(
                taxa_file,
                orig_to_short,
                db,
                run_dataset,
                source_fasta_path=args.input,
            )
        except Exception as e:
            raise RuntimeError(f'Failed to read taxonomy assignments {taxa_file}: {e}') from e
        if not tax_entries_local:
            raise RuntimeError('No taxonomy assignment rows mapped to the current Performance Review IDs.')

        # persist taxonomy for mapped ids
        if tax_entries_local:
            try:
                db.insert_taxonomy(tax_entries_local)
                logging.getLogger(__name__).info("[DB] Inserted/updated taxonomy for %d ids from taxa_assignments", len(tax_entries_local))
            except Exception as e:
                raise RuntimeError(f'Failed to store provided taxonomy assignments: {e}') from e

        # create a synthetic classification-like file so downstream code that
        # expects `class_out` can run unchanged. Format: qid\tbest\tidentity\ttaxon\tconfidence
        try:
            class_out_path = Path(outdir) / 'taxa_assignments_classout.tsv'
            with open(class_out_path, 'w') as cf:
                cf.write('id\tbest\tidentity\ttaxon\tconfidence\n')
                try:
                    for row in iter_assignment_rows(taxa_file, source_fasta_path=args.input):
                        qid = row.get('qid', 'NA')
                        tax = row.get('tax') if row.get('tax') is not None else 'NA'
                        conf = row.get('confidence') if row.get('confidence') is not None else 'NA'
                        cf.write(f"{qid}\t\tNA\t{tax}\t{conf}\n")
                except Exception as exc:
                    raise RuntimeError(f'Failed to normalise provided taxonomy assignments: {exc}') from exc
            class_out = str(class_out_path)
        except Exception as exc:
            raise RuntimeError(f'Failed to create normalised taxonomy evidence: {exc}') from exc
    else:
        # Only run classification here if we did not already run an early
        # pre-insert classification for kingdom filtering and the user did
        # not supply taxa_assignments.
        ref_name = getattr(args, 'ref_name', None) or _classify_derive_db_name(effective_ref)
        run_main_db_name = getattr(args, 'main_ref', None) or ref_name

        if alt_databases:
            if early_class_out:
                logging.getLogger(__name__).info(
                    "[PERFORMANCE REVIEW] Domain-filter preclassification will not be reused for final output; "
                    "running multi-database classification so alternate taxonomy is persisted."
                )
            logging.getLogger(__name__).info(
                "[PERFORMANCE REVIEW] Multi-database classification: primary=%s, alt=%s, main=%s",
                ref_name, [n for _, _, n in alt_databases], run_main_db_name,
            )
            logging.getLogger(__name__).info(
                '[PERFORMANCE REVIEW][DIAG] CLASSIFICATION_ORCHESTRATOR_START input=%s primary=%s authoritative=%s alternatives=%s threads=%s',
                current_dataset_fasta, ref_name, run_main_db_name,
                [name for _ref, _taxa, name in alt_databases], threads,
            )
            class_out, all_class_results = classify.run_all_classifications(
                str(current_dataset_fasta), outdir,
                primary_ref=effective_ref,
                primary_taxa=effective_taxa_tsv,
                primary_name=ref_name,
                alt_refs=alt_databases,
                threads=threads,
                main_db=run_main_db_name,
            )
            logging.getLogger(__name__).info(
                '[PERFORMANCE REVIEW][DIAG] CLASSIFICATION_ORCHESTRATOR_COMPLETE primary_output=%s databases=%s result_sizes=%s',
                class_out, list(all_class_results), {name: len(rows) for name, rows in all_class_results.items()},
            )
            logging.getLogger(__name__).info('[PERFORMANCE REVIEW][DIAG] ALT_TAXONOMY_DB_STORE_START')
            _store_alt_taxonomy_in_db(db, all_class_results, run_main_db_name)
            logging.getLogger(__name__).info('[PERFORMANCE REVIEW][DIAG] ALT_TAXONOMY_DB_STORE_COMPLETE')
        elif early_class_out:
            class_out = early_class_out
        else:
            class_out = classify.run_classification(
                str(current_dataset_fasta), outdir, ref_fasta=effective_ref,
                taxa_tsv=effective_taxa_tsv, threads=threads, db_name=ref_name,
            )

    if (
        getattr(args, 'command', None) == 'performance-review'
        and domain_filter in ('bacteria', 'archaea')
        and 'gtdb' not in str(run_main_db_name).lower()
    ):
        raise SystemExit(
            '[PERFORMANCE REVIEW] bacterial/archaeal marker classification requires GTDB as the primary sequence-reference database. '
            'Name the GTDB reference with --ref-name GTDB or select it with --main-ref GTDB. '
            'Other databases remain available through --alt-ref for comparative reporting; partner-catalogue taxon controls pathway routing only, while GTDB/reference taxonomy controls pangenome grouping.'
        )

    # Persist primary classifier evidence before rolling novelty/genome coverage.
    # Partner-catalogue taxonomy is separately retained as the isolate identity standard.
    try:
        logging.getLogger(__name__).info(
            '[PERFORMANCE REVIEW][DIAG] AUTHORITATIVE_CLASSIFICATION_IMPORT_START path=%s database=%s',
            class_out, run_main_db_name,
        )
        tax_entries, dist_entries = load_classification_results_for_dataset(
            class_out or '',
            orig_to_short,
            db,
            run_dataset,
        )
    except Exception as e:
        raise SystemExit(f'[PERFORMANCE REVIEW] Authoritative classification output could not be parsed: {e}')
    logging.getLogger(__name__).info(
        '[PERFORMANCE REVIEW][DIAG] AUTHORITATIVE_CLASSIFICATION_IMPORT_COMPLETE taxonomy_rows=%d distance_rows=%d',
        len(tax_entries), len(dist_entries),
    )
    if tax_entries:
        db.insert_taxonomy(tax_entries)
        logging.getLogger(__name__).info("[DB] Inserted/updated taxonomy for %d ids", len(tax_entries))
    if dist_entries:
        db.insert_distances(dist_entries)
        logging.getLogger(__name__).info("[DB] Inserted/updated distances for %d ids", len(dist_entries))
    classification_evidence = []
    for row in iter_classification_rows(class_out or ''):
        query_id = str(row.get('qid') or '')
        mapped = orig_to_short.get(query_id) or (query_id if query_id in {item[0] for item in mapped_records} else None)
        if not mapped:
            continue
        classification_evidence.append({
            'sequence_id': mapped, 'ref_db': run_main_db_name,
            'best_hit': row.get('best'), 'identity': row.get('identity'),
            'query_coverage': row.get('query_coverage'), 'target_coverage': row.get('target_coverage'),
            'alignment_length': row.get('alignment_length'), 'query_length': row.get('query_length'),
            'target_length': row.get('target_length'), 'mismatches': row.get('mismatches'),
            'gaps': row.get('gaps'), 'source_classifier_confidence': row.get('confidence'),
            'assignment_status': row.get('assignment_status'),
            'taxonomy_support_pct': row.get('taxonomy_support_pct'),
            'taxonomy_candidate_count': row.get('taxonomy_candidate_count'),
            'evidence_tier': row.get('evidence_tier'),
            'alignment_evidence_status': row.get('alignment_evidence_status'),
            'assignment_evidence_level': row.get('assignment_evidence_level'),
            'evidence_interpretation': row.get('evidence_interpretation'),
        })
    db.upsert_classification_evidence(classification_evidence)

    # Register explicit cultured baselines, then migrate unlabelled historical
    # run datasets into the rolling candidate collection.
    novelty_baseline_datasets = _resolve_novelty_baseline_datasets(db, args)
    for baseline_name in novelty_baseline_datasets:
        db.upsert_dataset_role(
            baseline_name, 'baseline', genomes_available=True,
            baseline_tier=_normalise_cli_baseline_tier(baseline_name),
        )
    try:
        registered_roles = db.get_dataset_roles()
        with db.connect() as conn:
            datasets_in_db = [row[0] for row in conn.execute(
                "SELECT DISTINCT dataset FROM sequences WHERE dataset IS NOT NULL AND dataset != ''"
            ).fetchall()]
        for dataset_name in datasets_in_db:
            if dataset_name not in registered_roles:
                db.upsert_dataset_role(dataset_name, 'candidate', genomes_available=False)
    except Exception as e:
        raise RuntimeError(f'Could not establish rolling dataset roles: {e}') from e

    # Assess only this submission. Novelty searches still use every registered
    # candidate and baseline in the database as rolling comparison context.
    metrics_query_fasta = str(current_dataset_fasta)

    # Novelty metrics establish the authoritative project nearest-hit evidence.
    # Legacy novelty.tsv/novelty_matches.tsv are derived from that evidence below,
    # avoiding a second, differently configured search against another project pool.
    target_fasta = getattr(args, 'target', None)
    try:
        novelty_metrics_out = novelty.build_reference_novelty_metrics(
            metrics_query_fasta,
            effective_ref,
            outdir,
            threads=threads,
            db=db,
            run_dataset=run_dataset,
            target_fasta=target_fasta,
            baseline_datasets=novelty_baseline_datasets,
            supplementary_references=alt_databases,
            pangenome_target=getattr(
                args, 'pangenome_target', _selection_sets.DEFAULT_PANGENOME_TARGET,
            ),
        )
        logging.getLogger(__name__).info("[NOVELTY] Wrote novelty metrics to %s", novelty_metrics_out)
        novelty_out = novelty.run_novelty(
            str(current_dataset_fasta), effective_ref, outdir,
            db=db, run_dataset=run_dataset, threads=threads,
            target_fasta=target_fasta,
        )
    except Exception as e:
        raise SystemExit(f'[PERFORMANCE REVIEW] Required novelty metrics failed: {e}')

    # Safety-net: remove any run sequences with explicit non-matching kingdom assignments.
    try:
        _prune_dataset_by_kingdom(db, run_dataset, str(kingdom) if kingdom else None, '[PERFORMANCE REVIEW]')
    except Exception as e:
        raise RuntimeError(f'Final Performance Review domain filtering failed: {e}') from e

    # Initialise cluster-tracking variables used by both the tree section and
    # the assessment section below.  They will be populated if --collapse is
    # active; otherwise they stay empty and the assessment omits cluster columns.
    tree_fasta = current_dataset_fasta
    run_member_to_rep: dict = {}
    run_rep_to_members: dict = {}

    # update tree/alignment
    try:
        # report DB sequence count for diagnostics so users can verify that
        # registered baseline sequences exist and will be used to seed the backbone
        try:
            db_ids = db.get_all_ids()
            logging.getLogger(__name__).info("[TREE] DB contains %d sequences; these will be considered for backbone construction", len(db_ids))
        except Exception:
            logging.getLogger(__name__).info("[TREE] Could not determine DB sequence count before tree build")

        # Optionally collapse highly similar run sequences (for tree readability)
        if getattr(args, 'collapse', False):
            try:
                threshold = float(getattr(args, 'collapse_threshold', 99.9))
            except Exception:
                threshold = 99.9
            # build qid -> tax mapping from classification output so we only
            # cluster sequences that share the same taxonomic assignment
            qid_to_tax = {}
            try:
                if class_out:
                    with open(class_out) as cf:
                        next(cf, None)
                        for line in cf:
                            parts = line.strip().split('\t')
                            if not parts:
                                continue
                            q = parts[0]
                            tax = parts[3] if len(parts) > 3 else None
                            qid_to_tax[q] = tax
            except Exception:
                qid_to_tax = {}

            # group short ids by tax
            taxa_groups = {}
            try:
                for h, s in read_fasta(str(current_dataset_fasta)):
                    tax = qid_to_tax.get(h)
                    taxa_groups.setdefault(tax, []).append((h, s))
            except Exception:
                taxa_groups = {None: []}
            try:
                artefacts = collapse_fasta_within_taxa(
                    taxa_groups,
                    outdir,
                    'current_dataset_collapsed.fasta',
                    'collapsed_map.tsv',
                    'collapsed_members.tsv',
                    threshold=threshold,
                    threads=threads,
                    log_prefix='[COLLAPSE]',
                    strict=True,
                )
                tree_fasta = artefacts.collapsed_path
                run_member_to_rep = artefacts.member_to_rep or {}
                # Build rep -> members mapping
                for mem, rep in run_member_to_rep.items():
                    run_rep_to_members.setdefault(rep, []).append(mem)
                logging.getLogger(__name__).info("[COLLAPSE] Wrote collapsed fasta %s (reps=%d)", tree_fasta, len(artefacts.collapsed_records))
            except Exception as e:
                raise RuntimeError(f'Requested Performance Review collapse failed: {e}') from e

        # Optionally filter by phylum before tree building
        phylum_filter = getattr(args, 'phylum', None)
        if phylum_filter:
            try:
                # Parse classification to extract phylum-level assignments
                phylum_ids = set()
                if class_out:
                    with open(class_out) as cf:
                        next(cf, None)  # skip header
                        for line in cf:
                            parts = line.strip().split('\t')
                            if len(parts) >= 4:
                                qid = parts[0]
                                tax = parts[3]  # Full taxonomy string
                                # Extract phylum from taxonomy (first level after kingdom)
                                # Format is usually: k__Kingdom; p__Phylum; c__Class; ...
                                try:
                                    tax_parts = [t.strip() for t in tax.split(';')]
                                    for part in tax_parts:
                                        if part.startswith('p__'):
                                            phylum_name = part[3:].strip()
                                            # Normalise for comparison (case-insensitive, strip underscores)
                                            if phylum_name.lower().replace('_', '') == phylum_filter.lower().replace('_', ''):
                                                phylum_ids.add(qid)
                                            break
                                except Exception:
                                    pass  # could not parse taxonomy for phylum filter

                # Filter tree_fasta to only include sequences in phylum_ids
                if phylum_ids:
                    filtered_fasta = Path(outdir) / f'tree_sequences_phylum_{phylum_filter}.fasta'
                    filtered_records = []
                    total_checked = 0
                    try:
                        for h, s in read_fasta(str(tree_fasta)):
                            total_checked += 1
                            if h in phylum_ids:
                                filtered_records.append((h, s))
                    except Exception:
                        pass  # could not read tree FASTA for phylum filtering

                    if filtered_records:
                        write_fasta(filtered_records, str(filtered_fasta))
                        logging.getLogger(__name__).info(
                            "[PHYLUM] Filtered tree sequences from %d to %d for phylum: %s",
                            total_checked, len(filtered_records), phylum_filter,
                        )
                        tree_fasta = filtered_fasta
                    else:
                        raise RuntimeError(
                            f"No tree sequences remained for requested phylum '{phylum_filter}'"
                        )
                else:
                    raise RuntimeError(
                        f"Could not find requested phylum '{phylum_filter}' in classification output"
                    )
            except Exception as e:
                raise RuntimeError(f'Phylum filtering failed: {e}') from e

        # pass the Database object so the tree builder can export existing
        # registered baseline sequences from the DB to form the backbone alignment
        _force_rebuild = bool(getattr(args, 'force_rebuild', False))
        try:
            tree_seq_count = sum(1 for line in open(str(tree_fasta)) if line.startswith(">"))
            logging.getLogger(__name__).info("[TREE] About to build/update tree with %d sequences from: %s", tree_seq_count, tree_fasta)
        except Exception as e:
            logging.getLogger(__name__).warning("[TREE] Could not count sequences in tree_fasta: %s", e)
        if _force_rebuild:
            logging.getLogger(__name__).info("[TREE] Force rebuild requested (--force-rebuild/--rebuild-tree); tree will be rebuilt from scratch including all datasets")
        tree.initialise_or_update_tree(
            ref_fasta=effective_ref,
            user_fasta=str(tree_fasta),
            outdir=outdir,
            db=db,
            db_dataset=None,
            threads=threads,
            previous_review=None if _force_rebuild else getattr(args, 'previous_review', None),
            force_rebuild=_force_rebuild,
            anchor_file=getattr(args, 'anchors', None),
            tree_method=getattr(args, 'tree_method', 'fasttree'),
        )
        try:
            warning_rows = tree.collect_tree_build_warnings(user_fasta=str(tree_fasta), anchor_file=getattr(args, 'anchors', None), db=db, db_dataset=None)
            warning_rows.extend(tree.summarise_alignment_quality(str(Path(outdir) / 'current_alignment.fasta')))
            if warning_rows:
                warn_path = tree.write_tree_warning_tsv(outdir, warning_rows)
                logging.getLogger(__name__).warning("[TREE] Tree/alignment warnings written to %s", warn_path)
        except Exception as e:
            logging.getLogger(__name__).warning("[TREE] Failed to summarise tree/alignment quality: %s", e)
    except Exception as e:
        raise SystemExit(f'[PERFORMANCE REVIEW] Required MSA/tree update failed: {e}')

    # Build combined taxonomy for iTOL from DB sequences and current run results.
    combined_path = Path(outdir) / 'combined_taxonomy.tsv'
    merged = {}
    order = []

    try:
        previous_review = getattr(args, 'previous_review', None)
        previous_review_ids = None
        if previous_review:
            try:
                cand = next(
                    (path for path in _combined_taxonomy_candidates(previous_review) if path.exists()),
                    None,
                )
                if cand is not None:
                    previous_review_ids = read_combined_taxonomy_ids(cand)
            except Exception:
                pass  # could not read previous review combined taxonomy

        base_rows = collect_db_taxonomy_rows(db, previous_review_ids if previous_review_ids else None)
        merged_rows = merge_combined_taxonomy_rows(base_rows, class_out or '', orig_to_short, db)
        for rid, tax, conf in merged_rows:
            merged[rid] = (tax if tax is not None else 'NA', conf if conf is not None else 'NA')
            order.append(rid)
    except Exception:
        # fallback: leave merged empty and continue
        pass  # could not build combined taxonomy from DB and classification

    # write merged combined taxonomy preserving order (DB-order first, then new)
    try:
        write_combined_taxonomy_tsv(
            combined_path,
            [(iid, merged.get(iid, ('NA', 'NA'))[0], merged.get(iid, ('NA', 'NA'))[1]) for iid in order],
        )
        logging.getLogger(__name__).info("[ITOL] Wrote combined taxonomy to %s", combined_path)
    except Exception as e:
        logging.getLogger(__name__).warning("[ITOL] Failed to write combined taxonomy: %s", e)

    # Load previous-review ID map for iTOL when --previous-review is provided.
    id_map_for_itol = None
    previous_review = getattr(args, 'previous_review', None)
    if previous_review:
        try:
            id_map = {}
            cand_map = _find_preferred_id_map(str(previous_review))
            if cand_map:
                id_map = _load_id_map_from_tsv(cand_map, db=db)
            if id_map:
                id_map_for_itol = id_map
        except Exception:
            id_map_for_itol = None

        try:
            tree_path = Path(outdir) / 'current_tree.nwk'
            tfile = str(tree_path) if tree_path.exists() else _find_tree_file_in_dir(outdir)
            itol.generate_itol_colours(str(combined_path), outdir, user_colour_csv=getattr(args, 'user_colours', None), id_map=id_map_for_itol, tree_file=tfile, phylum_groups=getattr(args, 'group_phyla', None))
            logging.getLogger(__name__).info("[ITOL] Generated iTOL colour files in %s", outdir)
        except Exception as e:
            logging.getLogger(__name__).warning("[ITOL] Failed to generate iTOL files: %s", e)
    # Generate iTOL colour files for the Performance Review output directory.
    try:
        tree_path = Path(outdir) / 'current_tree.nwk'
        tfile = str(tree_path) if tree_path.exists() else _find_tree_file_in_dir(outdir)
        itol.generate_itol_colours(str(combined_path), outdir, user_colour_csv=getattr(args, 'user_colours', None), id_map=id_map_for_itol, tree_file=tfile, phylum_groups=getattr(args, 'group_phyla', None))
        logging.getLogger(__name__).info("[ITOL] Generated iTOL colour files in %s", outdir)
    except Exception as e:
        logging.getLogger(__name__).warning("[ITOL] Failed to generate iTOL files: %s", e)
    # Optional: write functional annotation iTOL datasets when provided
    try:
        func_tsv = getattr(args, 'functional', None)
        if func_tsv:
            try:
                written = itol.write_functional_annotations(str(func_tsv), outdir, id_map=id_map_for_itol)
                logging.getLogger(__name__).info("[ITOL] Wrote functional annotation iTOL files: %s", ','.join(written) if written else '(none)')
            except Exception as e:
                logging.getLogger(__name__).warning("[ITOL] Functional annotations generation failed: %s", e)
    except Exception:
        pass  # functional annotation argument parsing failed
    # Draft rumen functional groups (auto-generated from output taxonomy)
    if getattr(args, 'draft_rumen_functions', False):
        try:
            tsv_out, itol_out = itol.generate_rumen_function_draft(
                str(combined_path), outdir, id_map=id_map_for_itol
            )
            if tsv_out:
                logging.getLogger(__name__).info("[ITOL] Draft rumen functional annotation: %s", tsv_out)
            if itol_out:
                logging.getLogger(__name__).info("[ITOL] Rumen functional iTOL file: %s", itol_out)
        except Exception as e:
            logging.getLogger(__name__).warning("[ITOL] Draft rumen functions generation failed: %s", e)
    # produce dataset membership band (Filing Cabinet versus review batch)
    try:
        combined_tax = combined_path
        ids_in_order = []
        if combined_tax.exists():
            with open(combined_tax) as ct:
                next(ct, None)
                for line in ct:
                    iid = line.strip().split('\t')[0]
                    if iid:
                        ids_in_order.append(iid)
        # Query DB for exact dataset membership per id to support arbitrary dataset names
        ds_map = {}
        try:
            with db.connect() as conn:
                cur = conn.cursor()
                placeholders = ','.join('?' for _ in ids_in_order) if ids_in_order else ''
                if placeholders:
                    cur.execute(f"SELECT id, dataset FROM sequences WHERE id IN ({placeholders})", tuple(ids_in_order))
                    for iid_row, ds in cur.fetchall():
                        ds_map[iid_row] = ds or ''
        except Exception:
            ds_map = {}

        membership_path = Path(outdir) / 'itol_dataset_membership.itol'
        itol.write_dataset_membership_strip(str(membership_path), ids_in_order, ds_map)
        logging.getLogger(__name__).info("[ITOL] Wrote dataset membership ITOL to %s", membership_path)
        try:
            tree_ids = {
                header for header, _sequence in read_fasta(
                    str(Path(outdir) / 'current_alignment.fasta')
                )
            }
            unresolved_replacements = db.get_unresolved_marker_replacement_failures()
            flagged_tree_ids = set(unresolved_replacements) & tree_ids
            if flagged_tree_ids:
                marker_status_path = Path(outdir) / 'itol_marker_replacement_status.itol'
                itol.write_marker_replacement_strip(
                    str(marker_status_path), ids_in_order, flagged_tree_ids,
                )
                logging.getLogger(__name__).warning(
                    "[ITOL] %d tree leaf/leaves retain an old marker after resequencing QC failed; "
                    "status strip written to %s",
                    len(flagged_tree_ids), marker_status_path,
                )
        except Exception as exc:
            raise RuntimeError(f'Could not annotate retained markers after failed resequencing: {exc}') from exc
        tier_map = {}
        try:
            with db.connect() as conn:
                cur = conn.cursor()
                placeholders = ','.join('?' for _ in ids_in_order) if ids_in_order else ''
                if placeholders:
                    cur.execute(
                        "SELECT s.id, COALESCE(NULLIF(s.baseline_tier, ''), r.baseline_tier), COALESCE(r.role, '') "
                        "FROM sequences s LEFT JOIN dataset_roles r ON r.dataset = s.dataset "
                        f"WHERE s.id IN ({placeholders})",
                        tuple(ids_in_order),
                    )
                    for iid_row, tier, role in cur.fetchall():
                        if role == 'baseline' and str(tier or '').strip().lower() in {'priority', 'secondary'}:
                            tier_map[iid_row] = str(tier).strip().lower()
        except Exception:
            tier_map = {}
        if tier_map:
            tier_path = Path(outdir) / 'itol_baseline_tier.itol'
            itol.write_baseline_tier_strip(str(tier_path), ids_in_order, tier_map)
            logging.getLogger(__name__).info("[ITOL] Wrote baseline-tier ITOL to %s", tier_path)
    except Exception as e:
        logging.getLogger(__name__).warning("[ITOL] Failed to write dataset membership ITOL: %s", e)

    # Add a novelty-gradient colour strip for new sequences.
    try:
        try:
            run_ids = [r[0] for r in mapped_records]
        except Exception:
            run_ids = []
        if not run_ids:
            try:
                with db.connect() as conn:
                    cur = conn.cursor()
                    cur.execute("SELECT id FROM sequences WHERE dataset = ?", (run_dataset,))
                    run_ids = [r[0] for r in cur.fetchall()]
            except Exception:
                run_ids = []

        if run_ids:
            novelty.build_run_novelty_itol(
                outdir,
                run_ids,
                str(current_dataset_fasta),
                db,
                run_dataset,
                orig_to_short,
                threads=getattr(args, 'threads', None),
            )
    except Exception as e:
        logging.getLogger(__name__).warning("[ITOL] Failed to write novelty ITOL: %s", e)

    # placement warning summary combining classification + novelty support
    try:
        if class_out:
            warning_rows = build_placement_warning_rows(class_out, novelty_out, orig_to_short, db)
            if warning_rows:
                warn_path = write_placement_warning_tsv(Path(outdir) / 'placement_warnings.tsv', warning_rows)
                logging.getLogger(__name__).warning("[PERFORMANCE REVIEW] Placement warnings written to %s", warn_path)
            try:
                run_ids_for_assessment = [record[0] for record in mapped_records]
                # Determine which sequences are actually in the tree
                try:
                    _tree_ids = set()
                    for _h, _s in read_fasta(str(tree_fasta)):
                        _tree_ids.add(_h)
                except Exception:
                    _tree_ids = None

                # Collect alt-db taxonomy from DB for the run sequences
                alt_ref_dbs: list = []
                alt_taxonomies: dict = {}
                alt_classification_evidence: dict = {}
                try:
                    alt_ref_dbs = db.get_alt_ref_dbs()
                    if alt_ref_dbs:
                        alt_taxonomies = db.get_taxonomy_alt_for_ids(run_ids_for_assessment)
                        all_classification_evidence = db.get_classification_evidence(run_ids_for_assessment)
                        for (sequence_id, ref_db), evidence_row in all_classification_evidence.items():
                            alt_classification_evidence.setdefault(sequence_id, {})[ref_db] = evidence_row
                        logging.getLogger(__name__).info(
                            "[PERFORMANCE REVIEW] Found alt-db taxonomy for ref_dbs: %s", alt_ref_dbs
                        )
                except Exception as e:
                    logging.getLogger(__name__).warning("[PERFORMANCE REVIEW] Could not load alt-db taxonomy: %s", e)

                assessment_rows = build_sequence_assessment_rows(
                    run_ids_for_assessment,
                    class_out,
                    str(Path(outdir) / 'novelty_metrics.tsv'),
                    warning_rows,
                    orig_to_short,
                    db,
                    member_to_rep=run_member_to_rep if run_member_to_rep else None,
                    rep_to_members=run_rep_to_members if run_rep_to_members else None,
                    tree_ids=_tree_ids,
                    alt_taxonomies=alt_taxonomies,
                    alt_ref_dbs=alt_ref_dbs,
                    alt_classification_evidence=alt_classification_evidence,
                )
                from branchmanager.marker_provenance import marker_qc_flag
                project_provenance = db.get_sequence_provenance(run_ids_for_assessment)
                from branchmanager.marker_evidence import get_active_marker_evidence
                project_marker_evidence = get_active_marker_evidence(db, run_ids_for_assessment)
                unresolved_replacements = db.get_unresolved_marker_replacement_failures(run_ids_for_assessment)
                for assessment_row in assessment_rows:
                    sequence_id = str(assessment_row.get('id'))
                    provenance = project_provenance.get(sequence_id) or {}
                    evidence = project_marker_evidence.get(sequence_id) or {}
                    qc_flag = marker_qc_flag(provenance)
                    assessment_row.update({
                        'marker_source_type': evidence.get('marker_source_type', 'LEGACY_UNKNOWN'),
                        'marker_evidence_tier': evidence.get('marker_evidence_tier', 'LOW'),
                        'marker_review_decision': evidence.get('marker_review_decision', ''),
                        'marker_review_reviewer': evidence.get('marker_review_reviewer', ''),
                        'marker_review_notes': evidence.get('marker_review_notes', ''),
                        'marker_qc_class': provenance.get('marker_qc_class', 'QUALITY_UNVERIFIED'),
                        'marker_qc_recommendation': provenance.get('marker_qc_recommendation', 'MANUAL_REVIEW'),
                        'marker_qc_reasons': provenance.get('marker_qc_reasons', 'no_marker_qc_provenance'),
                        'marker_manual_review_status': provenance.get('manual_review_status', 'NOT_REVIEWED'),
                        'marker_qc_flag': qc_flag or 'MARKER_QC_PASSED',
                        'marker_source_manifest': provenance.get('source_manifest', ''),
                        'chimera_call': provenance.get('chimera_call', 'NOT_RUN'),
                        'chimera_score': provenance.get('chimera_score', 'NA'),
                    })
                    replacement_failure = unresolved_replacements.get(sequence_id)
                    if replacement_failure:
                        assessment_row['marker_qc_flag'] = 'MARKER_REPLACEMENT_QC_FAILED'
                        prior_reasons = str(assessment_row.get('marker_qc_reasons') or '').strip()
                        replacement_reason = (
                            'replacement QC failed; old stored marker retained'
                            + (f': {replacement_failure["detail"]}' if replacement_failure.get('detail') else '')
                        )
                        assessment_row['marker_qc_reasons'] = ';'.join(
                            value for value in (prior_reasons, replacement_reason) if value
                        )
                        qc_flag = 'MARKER_REPLACEMENT_QC_FAILED'
                    if qc_flag:
                        flags = [flag for flag in str(assessment_row.get('placement_flags') or '').split(';') if flag]
                        if qc_flag not in flags:
                            flags.append(qc_flag)
                        assessment_row['placement_flags'] = ';'.join(flags)
                    chimera_call = str(provenance.get('chimera_call') or 'NOT_RUN').upper()
                    if chimera_call != 'PASS':
                        chimera_flag = {
                            'CHIMERA': 'CHIMERA_DETECTED',
                            'INDETERMINATE': 'CHIMERA_INDETERMINATE',
                            'SKIPPED': 'CHIMERA_CHECK_SKIPPED',
                        }.get(chimera_call, 'CHIMERA_NOT_RUN')
                        flags = [flag for flag in str(assessment_row.get('placement_flags') or '').split(';') if flag]
                        if chimera_flag not in flags:
                            flags.append(chimera_flag)
                        assessment_row['placement_flags'] = ';'.join(flags)

                    # Explicitly surface manual approvals that coexist with
                    # severe automated failures. Approval is an auditable
                    # human decision, not evidence that a chimera/QC failure
                    # has been biologically resolved.
                    review_status = str(assessment_row.get('marker_manual_review_status') or '').upper()
                    qc_class = str(assessment_row.get('marker_qc_class') or '').upper()
                    qc_recommendation = str(assessment_row.get('marker_qc_recommendation') or '').upper()
                    chimera_call = str(assessment_row.get('chimera_call') or 'NOT_RUN').upper()
                    severe_reasons = []
                    if chimera_call in {'CHIMERA', 'CHIMERA_CONFIRMED'}:
                        severe_reasons.append('CHIMERA_DETECTED')
                    elif chimera_call == 'INDETERMINATE':
                        severe_reasons.append('CHIMERA_CHECK_INDETERMINATE')
                    if qc_class == 'FAIL_QC' or qc_recommendation in {'RESEQUENCE', 'FAIL_QC'}:
                        severe_reasons.append('MARKER_QC_FAIL_OR_RESEQUENCE')
                    if review_status == 'APPROVED' and severe_reasons:
                        assessment_row['manual_approval_warning'] = ';'.join(severe_reasons)
                        flags = [
                            flag for flag in str(assessment_row.get('placement_flags') or '').split(';')
                            if flag
                        ]
                        if 'MANUAL_APPROVAL_SEVERE_QC_OVERRIDE' not in flags:
                            flags.append('MANUAL_APPROVAL_SEVERE_QC_OVERRIDE')
                        assessment_row['placement_flags'] = ';'.join(flags)

                manual_override_warnings = [
                    {
                        'SequenceID': row.get('id', ''),
                        'PartnerID': row.get('partner_id', ''),
                        'ManualReviewStatus': row.get('marker_manual_review_status', ''),
                        'ManualReviewDecision': row.get('marker_review_decision', ''),
                        'MarkerQCClass': row.get('marker_qc_class', ''),
                        'MarkerQCRecommendation': row.get('marker_qc_recommendation', ''),
                        'ChimeraCall': row.get('chimera_call', ''),
                        'ChimeraScore': row.get('chimera_score', ''),
                        'WarningSeverity': (
                            'CRITICAL' if 'CHIMERA_DETECTED' in str(row.get('manual_approval_warning') or '')
                            else 'HIGH'
                        ),
                        'Warning': row.get('manual_approval_warning', ''),
                        'Reviewer': row.get('marker_review_reviewer', ''),
                        'ReviewNotes': row.get('marker_review_notes', ''),
                    }
                    for row in assessment_rows if row.get('manual_approval_warning')
                ]
                if manual_override_warnings:
                    import csv as _csv
                    warning_path = Path(outdir) / 'manual_qc_override_warnings.tsv'
                    with warning_path.open('w', newline='', encoding='utf-8') as warning_handle:
                        writer = _csv.DictWriter(
                            warning_handle,
                            fieldnames=list(manual_override_warnings[0]),
                            delimiter='\t', lineterminator='\n',
                        )
                        writer.writeheader()
                        writer.writerows(manual_override_warnings)
                    logging.getLogger(__name__).warning(
                        '[MARKER QC] %d manually approved isolate(s) retain severe automated QC/chimera failures; see %s',
                        len(manual_override_warnings), warning_path,
                    )

                mwl_path = getattr(args, 'mwl', None)
                mwl_entries = []
                if mwl_path:
                    try:
                        logging.getLogger(__name__).info(
                            "[MWL] Matching the primary Taxonomy column against MWL. For multi-db runs, set --main-ref to your GTDB reference name."
                        )
                        mwl_entries = _mwl.load_mwl_entries(
                            str(mwl_path),
                            sheet_name=getattr(args, 'mwl_sheet', 'MWL_V1') or 'MWL_V1',
                        )
                        _mwl.annotate_assessment_rows(
                            assessment_rows,
                            mwl_entries,
                            min_rank=getattr(args, 'mwl_min_rank', 'p') or 'p',
                        )
                        mwl_report = _mwl.write_mwl_matches_tsv(
                            Path(outdir) / 'mwl_matches.tsv',
                            assessment_rows,
                        )
                        logging.getLogger(__name__).info(
                            "[MWL] Annotated %d assessment rows with %d MWL entries → %s",
                            len(assessment_rows), len(mwl_entries), mwl_report,
                        )
                    except Exception as e:
                        raise RuntimeError(f'MWL annotation failed: {e}') from e
                _tree_nwk = str(Path(outdir) / 'current_tree.nwk')
                try:
                    _cluster_summary, _clusters_csv, _backup_tsv = _cluster_report.generate_cluster_reports(
                        outdir=outdir,
                        assessment_rows=assessment_rows,
                        tree_path=_tree_nwk if Path(_tree_nwk).exists() else None,
                    )
                    if mwl_entries:
                        try:
                            _mwl.add_evaluation_scores(assessment_rows)
                            _mwl.write_mwl_matches_tsv(Path(outdir) / 'mwl_matches.tsv', assessment_rows)
                        except Exception as _mwe:
                            logging.getLogger(__name__).warning("[MWL] Failed to refresh MWL evaluation scores: %s", _mwe)
                    if _cluster_summary:
                        logging.getLogger(__name__).info(
                            "[CLUSTER] Wrote cluster summary → %s",
                            _cluster_summary,
                        )
                    if _clusters_csv:
                        logging.getLogger(__name__).info(
                            "[CLUSTER] Wrote consolidated cluster membership → %s", _clusters_csv,
                        )
                    if _backup_tsv:
                        logging.getLogger(__name__).info(
                            "[CLUSTER] Wrote backup candidates table → %s", _backup_tsv,
                        )
                except Exception as _ce:
                    raise RuntimeError(f'Cluster report generation failed: {_ce}') from _ce

                # Group nearby assessed leaves into compact, labelled local-clade figures.
                try:
                    if Path(_tree_nwk).exists():
                        _neighbourhood_result = neighbourhood.generate_local_neighbourhood_visuals(
                            tree_path=_tree_nwk,
                            assessment_rows=assessment_rows,
                            db=db,
                            outdir=Path(outdir) / 'neighbourhoods',
                            alignment_path=Path(outdir) / 'current_alignment.fasta',
                            image_format=getattr(args, 'neighbourhood_format', 'png'),
                        )
                        logging.getLogger(__name__).info(
                            "[NEIGHBOURHOOD] %d figure(s), %d/%d assessed sequences resolved",
                            len(_neighbourhood_result['figures']),
                            _neighbourhood_result['resolved'],
                            len(assessment_rows),
                        )
                        if _neighbourhood_result['resolved'] != len(assessment_rows):
                            logging.getLogger(__name__).warning(
                                "[NEIGHBOURHOOD] Local tree context resolved %d/%d assessed "
                                "sequences; unresolved rows remain in neighbourhood_manifest.tsv "
                                "with LocalNeighbourhoodFigure=NA",
                                _neighbourhood_result['resolved'],
                                len(assessment_rows),
                            )
                    else:
                        logging.getLogger(__name__).warning(
                            "[NEIGHBOURHOOD] No current_tree.nwk was produced; local-clade figures were skipped"
                        )
                except Exception as _ne:
                    raise RuntimeError(f'Local-clade figure generation failed: {_ne}') from _ne

                # Propose rolling primary + backup isolates within the most-specific
                # GTDB rank available, partitioned by local tree clade.
                try:
                    sequencing_sets_path = _selection_sets.build_sequencing_sets(
                        assessment_rows,
                        Path(outdir) / 'sequencing_sets.tsv',
                        tree_path=_tree_nwk if Path(_tree_nwk).exists() else None,
                        db=db,
                        pangenome_target=getattr(
                            args, 'pangenome_target', _selection_sets.DEFAULT_PANGENOME_TARGET,
                        ),
                        candidate_set_size=getattr(
                            args, 'candidate_set_size', _selection_sets.DEFAULT_CANDIDATE_SET_SIZE,
                        ),
                        baseline_redundancy_identity=getattr(
                            args, 'baseline_redundancy_identity',
                            _selection_sets.DEFAULT_BASELINE_REDUNDANCY_IDENTITY,
                        ),
                        baseline_redundancy_min_query_coverage=getattr(
                            args, 'baseline_redundancy_min_query_coverage',
                            _selection_sets.DEFAULT_BASELINE_REDUNDANCY_QUERY_COVERAGE,
                        ),
                        baseline_extension_min_identity=getattr(
                            args, 'baseline_extension_min_identity',
                            _selection_sets.DEFAULT_BASELINE_EXTENSION_MIN_IDENTITY,
                        ),
                        baseline_extension_min_query_coverage=getattr(
                            args, 'baseline_extension_min_query_coverage',
                            _selection_sets.DEFAULT_BASELINE_EXTENSION_MIN_QUERY_COVERAGE,
                        ),
                    )
                    hierarchy_path = _selection_sets.write_pangenome_hierarchy_tsv(
                        assessment_rows,
                        Path(outdir) / 'pangenome_hierarchy.tsv',
                        tree_path=_tree_nwk if Path(_tree_nwk).exists() else None,
                        db=db,
                        pangenome_target=getattr(
                            args, 'pangenome_target', _selection_sets.DEFAULT_PANGENOME_TARGET,
                        ),
                    )
                    logging.getLogger(__name__).info(
                        "[PANGENOME] Wrote rank-stratified catalogue-based analysis to %s",
                        hierarchy_path,
                    )
                    logging.getLogger(__name__).info(
                        "[SELECTION SETS] Wrote rolling clade candidate sets to %s",
                        sequencing_sets_path,
                    )
                    current_statuses = db.get_isolate_statuses(
                        [row.get('id') for row in assessment_rows if row.get('id')]
                    )
                    protected_states = {
                        'SAB_APPROVED', 'DNA_EXTRACTION_PENDING', 'DNA_EXTRACTION_FAILED',
                        'LIBRARY_PENDING', 'SEQUENCED', 'GENOME_QC_FAILED',
                        'GENOME_QC_PASSED', 'WITHDRAWN',
                    }
                    for row in assessment_rows:
                        if str(row.get('sequencing_set_role') or '').upper() not in {'PRIMARY', 'BACKUP'}:
                            continue
                        sequence_id = str(row.get('id') or '')
                        if current_statuses.get(sequence_id, {}).get('status') in protected_states:
                            continue
                        db.update_isolate_status(
                            sequence_id, 'PROPOSED',
                            detail=str(row.get('sequencing_set_reason') or ''),
                            source_file=str(sequencing_sets_path),
                        )
                except Exception as _se:
                    raise RuntimeError(f'Candidate-set generation failed: {_se}') from _se

                # Candidate-set roles/ranks are assigned after local contexts
                # are known. Re-render the same deterministic figure groups so
                # PRIMARY/BACKUP/ALTERNATE labels appear in the final images.
                try:
                    if Path(_tree_nwk).exists():
                        neighbourhood.generate_local_neighbourhood_visuals(
                            tree_path=_tree_nwk,
                            assessment_rows=assessment_rows,
                            db=db,
                            outdir=Path(outdir) / 'neighbourhoods',
                            alignment_path=Path(outdir) / 'current_alignment.fasta',
                            image_format=getattr(args, 'neighbourhood_format', 'png'),
                        )
                except Exception as _ne:
                    raise RuntimeError(f'Ranked local-clade figure refresh failed: {_ne}') from _ne

                # Carry through every partner-submitted isolate ID, even those
                # that never produced a marker sequence (e.g. all reads failed
                # QC before assembly). Without this, such isolates silently
                # disappear from the final reports instead of surfacing as an
                # explicit resequence-required item for the hiring panel.
                try:
                    rejected_rows = _load_rejected_before_assembly(
                        Path(outdir) / 'rejected_before_assembly.tsv'
                    )
                    known_ids = {str(row.get('id')) for row in assessment_rows}
                    added = 0
                    for rejected in rejected_rows:
                        rid = str(rejected['id'])
                        if rid in known_ids:
                            continue
                        assessment_rows.append({
                            'id': rid,
                            'partner_id': rejected.get('partner_id') or rid,
                            'dataset': run_dataset,
                            'selected_for_genome_sequencing': bool(rejected.get('selected_for_sequencing')),
                            'already_sequenced': bool(rejected.get('already_sequenced')),
                            'marker_qc_class': 'FAIL_QC',
                            'marker_qc_recommendation': 'RESEQUENCE',
                            'marker_qc_reasons': rejected.get('reason') or 'no_output_sequence',
                            'marker_manual_review_status': 'NOT_REVIEWED',
                            'marker_qc_flag': 'MARKER_QC_FAILED',
                            'placement_flags': 'NO_MARKER_SEQUENCE',
                            'no_marker_sequence': True,
                            'chimera_call': 'NOT_RUN',
                            'chimera_score': 'NA',
                        })
                        known_ids.add(rid)
                        added += 1
                    if added:
                        logging.getLogger(__name__).info(
                            "[PERFORMANCE REVIEW] Added %d resequence-required placeholder row(s) "
                            "for isolates with no assembled marker sequence", added,
                        )
                    report_ids = {row['sequence_id'] for row in replacement_report_rows}
                    for rejected in rejected_rows:
                        source_id = str(rejected.get('id') or '')
                        resolved_id = db.resolve_sequence_id(source_id) if source_id else None
                        if not resolved_id or resolved_id in report_ids:
                            continue
                        with db.connect() as conn:
                            existing_row = conn.execute(
                                'SELECT length(sequence) FROM sequences WHERE id=?',
                                (resolved_id,),
                            ).fetchone()
                        previous_provenance = db.get_sequence_provenance([resolved_id]).get(resolved_id, {})
                        replacement_report_rows.append(_marker_replacement_report_row(
                            resolved_id,
                            {
                                'status': 'NO_INCOMING_CONSENSUS',
                                'detail': rejected.get('reason') or 'incoming reads did not produce a QC-passing marker consensus',
                                'retained_marker': 'EXISTING' if existing_row else 'NONE',
                                'previous_length': existing_row[0] if existing_row else '',
                                'previous_qc_class': previous_provenance.get('marker_qc_class') or '',
                                'previous_qc_reasons': previous_provenance.get('marker_qc_reasons') or '',
                                'previous_read_ids': previous_provenance.get('source_read_ids') or '',
                                'previous_read_files': previous_provenance.get('source_read_files') or '',
                                'incoming_qc_class': 'FAIL_QC',
                                'incoming_read_files': ';'.join(_incoming_ab1_files(resolved_id)),
                            },
                        ))
                        report_ids.add(resolved_id)
                except Exception as e:
                    logging.getLogger(__name__).warning(
                        "[PERFORMANCE REVIEW] Could not merge rejected-before-assembly isolates into final reports: %s", e
                    )

                try:
                    replacement_report_path = Path(outdir) / 'marker_replacement_report.tsv'
                    _write_marker_replacement_report(
                        replacement_report_path,
                        replacement_report_rows,
                    )
                    history = db.get_marker_replacement_history()
                    history_path = Path(outdir) / 'marker_replacement_history.tsv'
                    fields = [
                        'attempt_id', 'sequence_id', 'status', 'detail', 'old_sequence_sha256',
                        'incoming_sequence_sha256', 'incoming_length', 'incoming_qc_class',
                        'incoming_qc_reasons', 'source_manifest', 'source_fasta', 'attempted_at',
                    ]
                    with open(history_path, 'w', newline='') as handle:
                        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
                        writer.writeheader()
                        writer.writerows(history)
                    logging.getLogger(__name__).info(
                        "[PERFORMANCE REVIEW] Wrote %d marker replacement attempt(s) to %s",
                        len(history), history_path,
                    )
                    logging.getLogger(__name__).info(
                        "[PERFORMANCE REVIEW] Wrote %d per-isolate replacement outcome(s) to %s",
                        len(replacement_report_rows), replacement_report_path,
                    )
                except Exception as e:
                    raise RuntimeError(f'Marker replacement history could not be written: {e}') from e

                # Final report write after tree-derived fields have been attached.
                assess_path = write_sequence_assessment_tsv(
                    Path(outdir) / 'sequence_assessment.tsv', assessment_rows,
                    assessment_db_name=run_main_db_name,
                )
                selection_summary_path = write_selection_summary_tsv(
                    Path(outdir) / 'selection_summary.tsv', assessment_rows,
                    assessment_db_name=run_main_db_name,
                )
                baseline_hits_path = write_baseline_hits_tsv(Path(outdir) / 'baseline_hits.tsv', assessment_rows)
                try:
                    snapshot_id = f'{run_dataset}:{Path(outdir).resolve()}'
                    snapshot_source = Path(outdir).resolve() / 'assessment' / 'sequence_assessment.tsv'
                    saved = db.save_assessment_snapshot(
                        snapshot_id,
                        assessment_rows,
                        dataset=run_dataset,
                        source_path=str(snapshot_source),
                    )
                    logging.getLogger(__name__).info(
                        '[ASSESSMENT] Saved %d rows to project snapshot %s', saved, snapshot_id,
                    )
                except Exception as e:
                    raise RuntimeError(f'Could not persist project assessment snapshot: {e}') from e
                logging.getLogger(__name__).info("[PERFORMANCE REVIEW] Wrote sequence assessment to %s", assess_path)
                logging.getLogger(__name__).info("[PERFORMANCE REVIEW] Wrote SAB selection summary to %s", selection_summary_path)
                logging.getLogger(__name__).info("[PERFORMANCE REVIEW] Wrote nearest baseline hit report to %s", baseline_hits_path)

                # Emit a concise summary using the same decisions as selection_summary.tsv.
                try:
                    decisions = [build_selection_decision(row)['decision'] for row in assessment_rows]
                    logging.getLogger(__name__).info(
                        "[ASSESSMENT SUMMARY] %d assessed: %d PRIMARY, %d BACKUP, "
                        "%d DIVERSITY/SECONDARY, %d BASELINE REDUNDANT, %d REVIEW, "
                        "%d TARGET MET, %d ALREADY SEQUENCED. "
                        "See assessment/selection_summary.tsv for the decision table.",
                        len(assessment_rows),
                        sum(d == 'PRIORITISE - SET PRIMARY' for d in decisions),
                        sum(d == 'RESERVE - SET BACKUP' for d in decisions),
                        sum(d in ('STRONG CANDIDATE', 'SECONDARY CANDIDATE', 'SECONDARY - STRAIN DIVERSITY') for d in decisions),
                        sum(d == 'EXCLUDE - BASELINE REDUNDANT' for d in decisions),
                        sum(d.startswith('REVIEW') for d in decisions),
                        sum(d == 'LOWER PRIORITY - TARGET MET' for d in decisions),
                        sum(d == 'ALREADY SEQUENCED' for d in decisions),
                    )
                except Exception:
                    pass  # could not compute assessment summary counts
            except Exception as e:
                raise RuntimeError(f'Required sequence assessment failed: {e}') from e
    except Exception as e:
        raise SystemExit(f'[PERFORMANCE REVIEW] Required assessment/placement reporting failed: {e}')

    # At end of run, keep the externally useful deliverables and organise them
    # into high-level folders.
    try:
        _organise_run_outputs(outdir, primary_db_name=run_main_db_name)
    except Exception as e:
        raise SystemExit(f'[PERFORMANCE REVIEW] Failed to organise required outputs: {e}')

    # Write a manifest after files have been organised.
    try:
        _write_output_explanations(outdir)
    except Exception:
        pass  # could not write output explanations


def cmd_performance_review(args):
    """Run Performance Review transactionally and publish only after validation."""
    import fcntl
    import sqlite3
    import uuid

    # Internal callers predating the CLI command namespace are retained for
    # focused unit tests; every user-facing CLI invocation takes the atomic path.
    if not hasattr(args, 'command') or not hasattr(args, 'skip_chimera_check'):
        return _cmd_performance_review_impl(args)
    original_db = Path(args.db).expanduser().resolve()
    if str(args.db) == ':memory:':
        raise SystemExit('[PERFORMANCE REVIEW] A persistent --db path is required for the rolling workflow.')
    tree.preflight_tree_tools(getattr(args, 'tree_method', 'fasttree'))
    outdir = Path(args.out).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    original_db.parent.mkdir(parents=True, exist_ok=True)
    staged_db = outdir / '.performance_review_project.sqlite'
    lock_path = original_db.with_name(original_db.name + '.lock')
    manifest = RunManifest(outdir, 'performance_review')
    for role, source in (
        ('marker_fasta', args.input), ('partner_metadata', getattr(args, 'partner_metadata', None)),
        ('primary_reference', getattr(args, 'ref', None)), ('reference_taxonomy', getattr(args, 'taxa', None)),
        ('baseline_fasta', getattr(args, 'baseline_fasta', None)), ('mwl', getattr(args, 'mwl', None)),
        ('marker_qc', getattr(args, 'marker_qc', None)), ('marker_review', getattr(args, 'marker_review', None)),
        ('baseline_taxonomy', getattr(args, 'baseline_taxa_assignments', None)),
        ('chimera_reference', getattr(args, 'chimera_ref', None)),
        ('target_collection', getattr(args, 'target', None)),
        ('anchors', getattr(args, 'anchors', None)),
    ):
        if source:
            manifest.add_input(source, role=role, required=True)
    for source in getattr(args, 'alt_ref', None) or []:
        manifest.add_input(source, role='alternate_reference')
    for source in getattr(args, 'alt_taxa', None) or []:
        if source:
            manifest.add_input(source, role='alternate_reference_taxonomy')
    run_id = f"performance-review:{getattr(args, 'dataset', 'dataset')}:{utc_now()}:{uuid.uuid4().hex[:8]}"

    with open(lock_path, 'a+') as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        try:
            if staged_db.exists():
                staged_db.unlink()
            if original_db.is_file():
                with sqlite3.connect(str(original_db)) as source, sqlite3.connect(str(staged_db)) as target:
                    source.backup(target)
                manifest.add_input(original_db, role='project_database_before')
            staged_args = argparse.Namespace(**vars(args))
            staged_args.db = str(staged_db)
            manifest.add_stage('performance_review', 'RUNNING', detail='taxonomy, novelty, MSA/tree, and assessment')
            _cmd_performance_review_impl(staged_args)

            safe_main = ''.join(
                char if char.isalnum() or char in ('_', '-') else '_'
                for char in str(getattr(args, 'main_ref', None) or getattr(args, 'ref_name', None) or 'GTDB')
            ).strip('_') or 'GTDB'
            required = {
                'full_assessment': outdir / 'assessment' / 'sequence_assessment.tsv',
                'board_summary': outdir / 'assessment' / 'selection_summary.tsv',
                'sequencing_sets': outdir / 'assessment' / 'sequencing_sets.tsv',
                'pangenome_hierarchy': outdir / 'assessment' / 'pangenome_hierarchy.tsv',
                'novelty_metrics': outdir / 'assessment' / 'novelty_metrics.tsv',
                'tree': outdir / 'tree' / 'current_tree.nwk',
                'alignment': outdir / 'tree' / 'current_alignment.fasta',
                'taxonomy': outdir / 'taxonomy' / f'{safe_main}.tsv',
            }
            if not required['taxonomy'].is_file():
                taxonomy_candidates = sorted((outdir / 'taxonomy').glob('*.tsv'))
                if taxonomy_candidates:
                    required['taxonomy'] = taxonomy_candidates[0]
            for role, path in required.items():
                manifest.add_output(path, role=role)
            # Optional per-sequence QC audit outputs are useful even when no
            # manual override was necessary; include whichever were produced.
            for role, name in (
                ('qc_summary', 'qc.stats'),
                ('qc_rejections', 'qc_rejections.tsv'),
                ('qc_manual_overrides', 'qc_manual_overrides.tsv'),
                ('expected_taxonomy_validation', 'expected_taxonomy_validation.tsv'),
            ):
                path = outdir / 'assessment' / name
                if path.is_file():
                    manifest.add_output(path, role=role)
            missing = manifest.verify_required_outputs()
            if missing:
                raise RuntimeError('required Performance Review outputs are missing: ' + ', '.join(missing))

            from branchmanager.reporting import (
                write_decision_changes, write_genome_collection_assessment,
                write_performance_review_dashboard,
            )
            staged = Database(str(staged_db))
            changes = write_decision_changes(staged, outdir / 'assessment' / 'decision_changes.tsv')
            manifest.add_output(changes, role='decision_changes')
            genome_assessment = write_genome_collection_assessment(
                staged, outdir / 'assessment' / 'genome_collection_assessment.tsv',
            )
            manifest.add_output(genome_assessment, role='genome_collection_assessment')
            manifest.add_stage('performance_review', 'COMPLETE')
            manifest.add_stage('hiring_panel', 'COMPLETE', detail='candidate sets and board summary')
            manifest.finish('COMPLETE')
            staged.record_project_run(
                run_id, 'performance_review', 'COMPLETE', dataset=getattr(args, 'dataset', ''),
                manifest_path=str(manifest.json_path),
                started_at=manifest.data['started_at'], completed_at=manifest.data['completed_at'],
            )

            publish_tmp = original_db.with_name(original_db.name + '.publishing')
            if publish_tmp.exists():
                publish_tmp.unlink()
            with sqlite3.connect(str(staged_db)) as source, sqlite3.connect(str(publish_tmp)) as target:
                source.backup(target)
                if target.execute('PRAGMA integrity_check').fetchone()[0] != 'ok':
                    raise RuntimeError('staged project database failed SQLite integrity_check')
            os.replace(publish_tmp, original_db)
            staged_db.unlink(missing_ok=True)
            # Capture the published database (post-assessment) rather than
            # leaving the dashboard labelled with only the pre-run DB hash.
            manifest.add_output(original_db, role='project_database_after')
            dashboard = write_performance_review_dashboard(
                outdir / 'assessment' / 'selection_summary.tsv', outdir,
                sequence_assessment=outdir / 'assessment' / 'sequence_assessment.tsv',
                pangenome_hierarchy=outdir / 'assessment' / 'pangenome_hierarchy.tsv',
                novelty_metrics=outdir / 'assessment' / 'novelty_metrics.tsv',
                run_manifest=outdir / 'run_manifest.json',
                genome_assessment=genome_assessment,
            )
            manifest.add_output(dashboard, role='performance_review_dashboard')
        except BaseException as exc:
            manifest.finish('FAILED', error=exc)
            raise
        finally:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)


def _detect_taxon_rank(taxon_query: str, rank_arg: str) -> str:
    """Return the single-letter rank key to filter on.

    Priority:
      1. Explicit ``--rank`` argument (if not 'auto')
      2. GTDB-style prefix embedded in the query (``p__``, ``f__``, ``g__`` …)
      3. Known domain keywords: ``archaea`` / ``bacteria`` → 'd'
      4. Fallback: 'p' (phylum)
    """
    from branchmanager.taxonomy import RANK_ALIASES
    if rank_arg and rank_arg.lower() != 'auto':
        return RANK_ALIASES.get(rank_arg.lower(), rank_arg.lower()[:1])
    if '__' in taxon_query:
        prefix = taxon_query.split('__')[0].lower().strip()
        return RANK_ALIASES.get(prefix, prefix[:1] if prefix else 'p')
    # Bare domain-level keywords — map to 'd' so we filter at domain rank
    if taxon_query.strip().lower() in ('archaea', 'bacteria', 'eukarya', 'eukaryota'):
        return 'd'
    return 'p'


def _strip_rank_prefix(taxon: str) -> str:
    """Strip a GTDB-style rank prefix: 'p__Bacillota' → 'Bacillota'."""
    if '__' in taxon:
        return taxon.split('__', 1)[1].strip()
    return taxon.strip()


def _taxon_name_matches(val: str, query: str) -> bool:
    """Case-insensitive, underscore/space-flexible exact match."""
    import re as _re
    def _n(s): return _re.sub(r'[\s_]+', '_', str(s).strip().lower())
    return bool(val) and _n(val) == _n(query)


def cmd_org_chart(args):
    """Build a focused phylogenetic tree for a specific taxon from an existing DB.

    Fast path: if a ``current_alignment.fasta`` already exists in
    ``--from-dir`` (or the output directory), the matching sequences are
    extracted from that *pre-built* alignment and FastTree is run on the
    subset — no re-alignment needed.

    Slow path (fallback): if no existing alignment is found, the matching
    sequences are exported from the DB and a full MAFFT + FastTree build
    is performed (same as a normal run).
    """
    from branchmanager.pipeline import itol as itol_mod
    from branchmanager.pipeline import tree as tree_mod  # noqa: F401 (used in _build_org_chart)
    from branchmanager.pipeline.workflow_helpers import write_combined_taxonomy_tsv
    from branchmanager.taxonomy import parse_taxon_string
    from branchmanager.utils.fasta import read_fasta, write_fasta  # noqa: F401

    outdir = args.out
    os.makedirs(outdir, exist_ok=True)
    _configure_logging(outdir)
    log = logging.getLogger(__name__)

    taxon_query = args.taxon
    rank_key = _detect_taxon_rank(taxon_query, getattr(args, 'rank', 'auto'))
    taxon_clean = _strip_rank_prefix(taxon_query)
    threads = getattr(args, 'threads', 4)
    min_seqs = getattr(args, 'min_seqs', 3)
    ref_fasta = getattr(args, 'ref', None)
    anchor_file = getattr(args, 'anchors', None)
    from_dir = getattr(args, 'from_dir', None) or outdir
    no_tree = getattr(args, 'no_tree', False)

    log.info(
        "[ORG-CHART] Taxon query: '%s'  rank: '%s'  clean name: '%s'  from-dir: %s",
        taxon_query, rank_key, taxon_clean, from_dir,
    )

    db = Database(args.db)
    db.initialise()

    with db.connect() as conn:
        cur = conn.cursor()
        cur.execute(
            "SELECT s.id, s.sequence, t.taxonomy, t.confidence, s.dataset "
            "FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id"
        )
        all_rows = cur.fetchall()

    matched = []
    for sid, seq, tax, conf, dataset in all_rows:
        if not tax:
            continue
        parsed = parse_taxon_string(tax)
        val = parsed.get(rank_key, '')
        if _taxon_name_matches(val, taxon_clean):
            matched.append((sid, seq or '', tax, conf or 'NA', dataset or ''))

    log.info("[ORG-CHART] DB has %d sequences total; %d match taxon '%s' at rank '%s'",
             len(all_rows), len(matched), taxon_clean, rank_key)

    if len(matched) < min_seqs:
        log.warning(
            "[ORG-CHART] Only %d sequences matched (minimum required: %d). "
            "Check taxon spelling and rank. Available phyla in DB: %s",
            len(matched), min_seqs,
            sorted({parse_taxon_string(r[2]).get('p', 'unknown')
                    for r in all_rows if r[2]}),
        )
        return

    # Collect matched IDs as a set for fast lookup
    matched_ids = {sid for sid, *_ in matched}

    combined_tax_path = Path(outdir) / 'org_chart_combined_taxonomy.tsv'
    write_combined_taxonomy_tsv(
        combined_tax_path,
        [(sid, tax, conf) for sid, _, tax, conf, _ in matched],
    )
    log.info("[ORG-CHART] Wrote taxonomy for %d sequences → %s", len(matched), combined_tax_path)

    summary_path = Path(outdir) / 'org_chart_sequence_list.tsv'
    with open(summary_path, 'w') as sf:
        sf.write('ID\tTaxonomy\tConfidence\tDataset\n')
        for sid, _, tax, conf, dataset in matched:
            sf.write(f"{sid}\t{tax}\t{conf}\t{dataset}\n")
    log.info("[ORG-CHART] Wrote sequence list → %s", summary_path)

    tree_path = Path(outdir) / 'org_chart_tree.nwk'

    if not no_tree:
        _build_org_chart(
            matched=matched,
            matched_ids=matched_ids,
            from_dir=from_dir,
            outdir=outdir,
            tree_path=tree_path,
            ref_fasta=ref_fasta,
            anchor_file=anchor_file,
            threads=threads,
            taxon_clean=taxon_clean,
            log=log,
        )
    else:
        log.info("[ORG-CHART] Tree build skipped (--no-tree).")

    try:
        tfile = str(tree_path) if tree_path.exists() else None
        itol_mod.generate_itol_colours(
            str(combined_tax_path),
            outdir,
            tree_file=tfile,
            phylum_groups=getattr(args, 'group_phyla', None),
        )
        log.info("[ORG-CHART] Generated iTOL colour files in %s", outdir)
    except Exception as e:
        log.warning("[ORG-CHART] iTOL colour generation failed: %s", e)

    # Optional: write functional annotation datasets when provided
    try:
        func_tsv = getattr(args, 'functional', None)
        if func_tsv:
            try:
                written = itol_mod.write_functional_annotations(str(func_tsv), outdir, id_map=None)
                log.info("[ORG-CHART] Wrote functional annotation iTOL files: %s", ','.join(written) if written else '(none)')
            except Exception as e:
                log.warning("[ORG-CHART] Functional annotations generation failed: %s", e)
    except Exception:
        pass  # functional annotation argument parsing failed
    # Draft rumen functional groups
    if getattr(args, 'draft_rumen_functions', False):
        try:
            tsv_out, itol_out = itol_mod.generate_rumen_function_draft(
                str(combined_tax_path), outdir, id_map=None
            )
            if tsv_out:
                log.info("[ORG-CHART] Draft rumen functional annotation: %s", tsv_out)
            if itol_out:
                log.info("[ORG-CHART] Rumen functional iTOL file: %s", itol_out)
        except Exception as e:
            log.warning("[ORG-CHART] Draft rumen functions generation failed: %s", e)

    # One colour per dataset label so users can see which sequences came from
    # which dataset (Filing Cabinet, review batch, etc.) in the same tree view.
    try:
        ids_in_order = [sid for sid, *_ in matched]
        ds_map = {sid: (dataset or 'unknown') for sid, _, _tax, _conf, dataset in matched}
        membership_path = Path(outdir) / 'itol_dataset_membership.itol'
        itol_mod.write_dataset_membership_strip(
            str(membership_path), ids_in_order, ds_map,
            dataset_label='Dataset membership',
        )
        log.info("[ORG-CHART] Wrote dataset membership strip → %s", membership_path)
    except Exception as e:
        log.warning("[ORG-CHART] Dataset membership strip failed: %s", e)

    log.info(
        "[ORG-CHART] Complete. %d sequences, taxon='%s', output=%s",
        len(matched), taxon_query, outdir,
    )
    print(
        f"[org-chart] Done: {len(matched)} sequences for '{taxon_query}' "
        f"→ {outdir}"
    )


def _build_org_chart(
    matched, matched_ids, from_dir, outdir, tree_path,
    ref_fasta, anchor_file, threads, taxon_clean, log,
):
    """Build a FastTree for the given matched sequences.

    Fast path: filter an existing ``current_alignment.fasta`` from
    ``from_dir`` to the matched IDs + any reference anchors already
    embedded in that alignment, then run FastTree directly — no
    re-alignment required.

    Slow path: write the raw sequences as a FASTA and run the full
    ``initialise_or_update_tree`` pipeline (MAFFT → FastTree).
    """
    from branchmanager.pipeline.tree import (
        _run_fasttree, _make_unique_fasta,
        is_ref_anchor,
    )
    from branchmanager.utils.fasta import read_fasta, write_fasta
    from branchmanager.pipeline import tree as tree_mod
    import re

    out = Path(outdir)
    tree.preflight_tree_tools('fasttree', require_mafft=False)
    # Reuse an existing alignment when available.
    aln_candidates = [
        Path(from_dir) / 'current_alignment.fasta',
        Path(from_dir) / 'tree' / 'current_alignment.fasta',
        Path(outdir) / 'current_alignment.fasta',
        Path(outdir) / 'tree' / 'current_alignment.fasta',
    ]
    existing_aln = next((p for p in aln_candidates if p.exists()), None)

    if existing_aln:
        log.info("[ORG-CHART] Fast path: filtering existing alignment %s", existing_aln)
        # Build a normalised ID lookup for the matched IDs
        def _norm_id(x):
            x = str(x).split()[0]
            if '|' in x:
                x = x.split('|')[-1]
            x = re.sub(r'_[0-9]+$', '', x)
            return x.lower()

        matched_norm = {_norm_id(sid) for sid in matched_ids}

        kept = []
        anchors_included = 0
        for header, seq in read_fasta(str(existing_aln)):
            hid = header.split()[0]
            if is_ref_anchor(hid):
                kept.append((header, seq))
                anchors_included += 1
            elif hid in matched_ids or _norm_id(hid) in matched_norm:
                kept.append((header, seq))

        n_data = len(kept) - anchors_included
        log.info(
            "[ORG-CHART] Fast path: %d data sequences + %d anchors extracted from alignment",
            n_data, anchors_included,
        )

        if n_data < 3:
            log.warning(
                "[ORG-CHART] Only %d data sequences found in existing alignment "
                "(IDs may not match). Falling back to slow path.", n_data,
            )
            existing_aln = None   # trigger slow path below
        else:
            org_chart_alignment = out / 'org_chart_alignment.fasta'
            write_fasta(kept, str(org_chart_alignment))
            fasta_for_tree, id_map = _make_unique_fasta(str(org_chart_alignment), out)
            if _run_fasttree(Path(fasta_for_tree), tree_path):
                _finalise_org_chart_tree(out, id_map, tree_path, log)
            else:
                log.warning("[ORG-CHART] FastTree failed on fast path")
            return

    log.info("[ORG-CHART] Slow path: full MAFFT + FastTree build")
    tree.preflight_tree_tools('fasttree', require_mafft=True)
    seqs_fasta = out / 'org_chart_input_sequences.fasta'
    write_fasta([(sid, seq) for sid, seq, *_ in matched], str(seqs_fasta))

    # Empty FASTA as user_fasta; all sequences go in via the FASTA directly
    # by passing them as user_fasta and disabling DB pull (db=None)
    try:
        tree_mod.initialise_or_update_tree(
            ref_fasta=ref_fasta or '',
            user_fasta=str(seqs_fasta),
            outdir=outdir,
            db=None,
            threads=threads,
            anchor_file=anchor_file,
            force_rebuild=True,
        )
        # The standard pipeline writes to current_tree.nwk; copy/rename
        default_tree = out / 'current_tree.nwk'
        if default_tree.exists():
            import shutil
            shutil.copy2(str(default_tree), str(tree_path))
            log.info("[ORG-CHART] Slow path tree written → %s", tree_path)
    except Exception as e:
        log.warning("[ORG-CHART] Slow path tree build failed: %s", e)


def _finalise_org_chart_tree(out: Path, id_map: dict, tree_path: Path, log) -> None:
    """Remap IDs, prune anchors, and label internal nodes for Org Chart output."""
    from branchmanager.pipeline.tree import (
        _repair_internal_node_label_delimiters, _label_internal_nodes,
        _prune_anchor_leaves, REF_ANCHOR_PREFIX,
    )
    if not tree_path.exists():
        return
    newick = tree_path.read_text()
    if id_map:
        for new_id, orig_id in id_map.items():
            newick = newick.replace(new_id, orig_id)
    before = newick.count(REF_ANCHOR_PREFIX)
    newick = _prune_anchor_leaves(newick)
    after = newick.count(REF_ANCHOR_PREFIX)
    if before:
        log.info("[ORG-CHART] Pruned %d anchor leaves from focused-tree newick", before - after)
    newick = _repair_internal_node_label_delimiters(newick)
    newick = _label_internal_nodes(newick)
    tree_path.write_text(newick)
    log.info("[ORG-CHART] Focused tree finalised → %s", tree_path)


def cmd_label_maker(args):
    db = Database(args.db)
    db.initialise()
    outdir = args.out
    os.makedirs(outdir, exist_ok=True)
    _configure_logging(outdir)
    log = logging.getLogger(__name__)

    try:
        with db.connect() as conn:
            cur = conn.cursor()
            # If the outdir contains a Filing Cabinet combined taxonomy file, prefer
            # to use only those IDs so the regenerated iTOL matches the previous review
            # statistics. This mirrors the behaviour used during `run` when a
            # --previous-review is supplied.
            try:
                p = Path(outdir)
                # Prefer an existing combined taxonomy file in the outdir, but
                # only if it contains data (header + >=1 data row). Fall back
                # to the DB-wide query otherwise to avoid regenerating empty
                # iTOL outputs when a stub file exists.
                previous_review_file = None
                for cand in _combined_taxonomy_candidates(p):
                    if not cand.exists():
                        continue
                    try:
                        # count lines cheaply
                        with open(cand) as pf:
                            cnt = sum(1 for _ in pf)
                        if cnt > 1:
                            previous_review_file = cand
                            break
                        else:
                            # file exists but only header or empty -> ignore
                            continue
                    except Exception:
                        continue  # could not count lines in candidate taxonomy file
                previous_review_ids = None
                if previous_review_file is not None:
                    previous_review_ids = read_combined_taxonomy_ids(previous_review_file)
                else:
                    previous_review_ids = None
            except Exception:
                pass  # could not discover previous review taxonomy file

            previous_review_id_list = previous_review_ids if isinstance(previous_review_ids, list) else []

            # fetch id, taxonomy, confidence and dataset for all sequences (or filter)
            if getattr(args, 'include_datasets', None):
                ds_list = [d.strip() for d in args.include_datasets.split(',') if d.strip()]
                placeholders = ','.join('?' for _ in ds_list)
                cur.execute(f"SELECT s.id, t.taxonomy, t.confidence, s.dataset FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id WHERE s.dataset IN ({placeholders})", tuple(ds_list))
            elif previous_review_id_list:
                placeholders = ','.join('?' for _ in previous_review_id_list)
                cur.execute(f"SELECT s.id, t.taxonomy, t.confidence, s.dataset FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id WHERE s.id IN ({placeholders})", tuple(previous_review_id_list))
            else:
                cur.execute("SELECT s.id, t.taxonomy, t.confidence, s.dataset FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id")
            rows = cur.fetchall()
    except Exception as e:
        log.warning("[LABEL-MAKER] Failed to query DB: %s", e)
        return

    try:
        kingdom = _sequence_domain_filter(args, default=None)
        if kingdom:
            kingdom_text = str(kingdom)
            rows = [row for row in rows if row[1] and taxonomy_matches_kingdom(str(row[1]), kingdom_text)]
            log.info("[LABEL-MAKER] Retained %d rows after kingdom filter=%s", len(rows), kingdom)
    except Exception as e:
        log.warning("[LABEL-MAKER] Kingdom filter failed: %s", e)

    combined_path = Path(outdir) / 'combined_taxonomy.tsv'
    try:
        write_combined_taxonomy_tsv(combined_path, [(rid, tax, conf) for rid, tax, conf, ds in rows])
        log.info("[LABEL-MAKER] Wrote combined taxonomy for %d ids to %s", len(rows), combined_path)
    except Exception as e:
        log.warning("[LABEL-MAKER] Failed to write combined taxonomy: %s", e)
        return

    # call itol generator
    try:
        tree_path = Path(outdir) / 'current_tree.nwk'
        tfile = str(tree_path) if tree_path.exists() else _find_tree_file_in_dir(outdir)
        itol.generate_itol_colours(str(combined_path), outdir, tree_file=tfile, phylum_groups=getattr(args, 'group_phyla', None))
        log.info("[LABEL-MAKER] Generated iTOL colour files in %s", outdir)
    except Exception as e:
        log.warning("[LABEL-MAKER] Colour generation failed: %s", e)

    # Optional: write functional annotation datasets when provided
    try:
        func_tsv = getattr(args, 'functional', None)
        if func_tsv:
            try:
                written = itol.write_functional_annotations(str(func_tsv), outdir, id_map=None)
                log.info("[LABEL-MAKER] Wrote functional annotation iTOL files: %s", ','.join(written) if written else '(none)')
            except Exception as e:
                log.warning("[LABEL-MAKER] Functional annotations generation failed: %s", e)
    except Exception:
        pass  # functional annotation argument parsing failed
    # Draft rumen functional groups
    if getattr(args, 'draft_rumen_functions', False):
        try:
            tsv_out, itol_out = itol.generate_rumen_function_draft(str(combined_path), outdir, id_map=None)
            if tsv_out:
                log.info("[LABEL-MAKER] Draft rumen functional annotation: %s", tsv_out)
            if itol_out:
                log.info("[LABEL-MAKER] Rumen functional iTOL file: %s", itol_out)
        except Exception as e:
            log.warning("[LABEL-MAKER] Draft rumen functions generation failed: %s", e)

    # build dataset membership strip
    try:
        ids_in_order = [r[0] for r in rows]
        ds_map = {r[0]: (r[3] or '') for r in rows}
        membership_path = Path(outdir) / 'itol_dataset_membership.itol'
        itol.write_dataset_membership_strip(str(membership_path), ids_in_order, ds_map)
        log.info("[LABEL-MAKER] Wrote dataset membership iTOL file to %s", membership_path)
    except Exception as e:
        log.warning("[LABEL-MAKER] Failed to build/write dataset membership iTOL file: %s", e)

    # write explanations for regenerated outputs
    try:
        _write_output_explanations(outdir)
    except Exception:
        pass  # could not write output explanations


# build_parser() moved to branchmanager.cli.parser (1,217 lines extracted)
# Imported at top of file.

def cmd_background_check(args):
    """Run the reusable Background Check classification workflow."""
    from branchmanager.pipeline import background_check as _background_check_module

    outdir = args.out
    os.makedirs(outdir, exist_ok=True)
    _configure_logging(outdir)
    log = logging.getLogger(__name__)

    # Parse datasets — supports two input styles:
    #   Style A (flat pairs):  --dataset Name1 /path1 Name2 /path2 ...
    #   Style B (NAME=FASTA):  --dataset Name1=/path1 --dataset Name2=/path2
    #   Style C (mixed):       --dataset Name1=/path1 Name2 /path2 ...
    raw = args.datasets or []
    # Flatten in case we ever switch back to action='append'
    flat: list = []
    for item in raw:
        if isinstance(item, list):
            flat.extend(item)
        else:
            flat.append(item)

    datasets = []
    i = 0
    while i < len(flat):
        token = flat[i]
        if '=' in token:
            # NAME=FASTA form
            name, _, fasta = token.partition('=')
            datasets.append((name.strip(), fasta.strip()))
            i += 1
        elif i + 1 < len(flat) and not flat[i + 1].startswith('-'):
            # NAME FASTA pair (next token is not a flag)
            datasets.append((token.strip(), flat[i + 1].strip()))
            i += 2
        else:
            raise SystemExit(
                f"[BACKGROUND CHECK] Cannot parse dataset spec at position {i}: '{token}'. "
                "Expected 'NAME FASTA' pairs or 'NAME=FASTA' pairs."
            )

    if not datasets:
        raise SystemExit("[BACKGROUND CHECK] At least one dataset is required. "
                         "Use: --dataset NAME /path/to/file.fasta")

    # Validate FASTA paths exist
    for name, fasta in datasets:
        if not os.path.exists(fasta):
            raise SystemExit(
                f"[BACKGROUND CHECK] FASTA file not found for dataset '{name}': {fasta}"
            )

    log.info("[BACKGROUND CHECK] Starting pre-classification for %d dataset(s)", len(datasets))

    pipeline_tsv = _background_check_module.run_background_check(
        datasets=datasets,
        ref_fasta=args.ref,
        outdir=outdir,
        taxa_tsv=getattr(args, 'taxa', None),
        threads=int(getattr(args, 'threads', 4) or 4),
        min_identity=float(getattr(args, 'min_identity', 0.80) or 0.80),
        low_confidence_threshold=float(
            getattr(args, 'low_confidence_threshold', 0.97) or 0.97
        ),
        max_hits=int(getattr(args, 'max_hits', 10) or 10),
        max_rejects=int(getattr(args, 'max_rejects', 256) or 256),
    )

    log.info("[BACKGROUND CHECK] Done. Pipeline-ready taxonomy → %s", pipeline_tsv)
    print(
        f"[background-check] Done.\n"
        f"  Pipeline taxonomy : {pipeline_tsv}\n"
        f"  Summary           : {os.path.join(outdir, 'background_check_summary.txt')}\n\n"
        f"Use in the Filing Cabinet:\n"
        f"  branchmanager filing-cabinet --fasta <fasta> --taxa-assignments {pipeline_tsv} "
        f"--db project.db --dataset <name> -o filing_cabinet_out/"
    )


def _resolve_mailroom_interview_input(
    source: str | os.PathLike,
    *,
    allow_review: bool = False,
) -> tuple[str, str, str]:
    """Resolve and validate the immutable Mailroom hand-off used by Interview."""
    supplied = Path(source).expanduser().resolve()
    if supplied.is_dir():
        map_path = supplied / 'ab1_map.tsv'
        summary_path = supplied / 'mailroom_summary.json'
    elif supplied.is_file():
        map_path = supplied
        summary_path = supplied.parent / 'mailroom_summary.json'
    else:
        raise SystemExit(f'[interview] Mailroom input does not exist: {supplied}')

    if map_path.name != 'ab1_map.tsv' or not map_path.is_file():
        raise SystemExit(
            f'[interview] Expected a Mailroom ab1_map.tsv, but it was not found at {map_path}'
        )
    if not summary_path.is_file():
        raise SystemExit(
            f'[interview] Expected the accompanying Mailroom summary at {summary_path}'
        )

    try:
        summary = json.loads(summary_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise SystemExit(f'[interview] Could not read {summary_path}: {exc}') from exc
    status = str(summary.get('status', '')).strip().upper()
    if status not in {'PASS', 'REVIEW_REQUIRED', 'FAIL'}:
        raise SystemExit(
            f'[interview] Mailroom summary has an invalid or missing status: {status or "<blank>"}'
        )
    if status == 'FAIL':
        raise SystemExit(
            '[interview] Mailroom status is FAIL. Correct the Mailroom report before AB1 QC.'
        )
    if status == 'REVIEW_REQUIRED' and not allow_review:
        raise SystemExit(
            '[interview] Mailroom status is REVIEW_REQUIRED. Resolve its findings or repeat with '
            '--allow-mailroom-review after documenting the review.'
        )

    with open(map_path, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = {str(field).strip() for field in (reader.fieldnames or [])}
        required = {
            'sequence_id', 'dataset', 'read_file', 'primer', 'direction',
            'processing_mode', 'primer_assignment',
        }
        missing = sorted(required - fields)
        if missing:
            raise SystemExit(
                '[interview] Mailroom ab1_map.tsv is missing required column(s): '
                + ', '.join(missing)
            )
        if next(reader, None) is None:
            raise SystemExit('[interview] Mailroom ab1_map.tsv contains no reads.')

    return str(map_path), str(summary_path), status



def _infer_qc_dataset_from_path(value) -> str:
    raw = str(value or '').strip()
    if not raw:
        return ''
    try:
        parts = [part for part in Path(raw).parts if part not in {'', '/'}]
    except TypeError:
        return ''
    for index, part in enumerate(parts):
        if part == 'mailroom' and index > 0:
            return parts[index - 1]
    for part in reversed(parts):
        if part.startswith('01_interview_') and len(part) > len('01_interview_'):
            return part[len('01_interview_'):]
    for index, part in enumerate(parts):
        if part in {'01_onboarding', '02_paper_trail_merge_meeting'} and index > 0:
            return parts[index - 1]
    return ''

def cmd_interview(args):
    """Run standalone AB1 QC from a validated Mailroom batch."""
    sample_map, summary, status = _resolve_mailroom_interview_input(
        args.mailroom,
        allow_review=bool(getattr(args, 'allow_mailroom_review', False)),
    )
    args.sample_map = sample_map
    args.mailroom_summary = summary
    args.mailroom_status = status
    args.input = []
    cmd_paper_trail(args)


def cmd_paper_trail(args):
    """Handler for standalone or Assistant Paper Trail chromatogram QC."""
    is_interview = getattr(args, 'command', None) == 'interview'
    command_name = 'interview' if is_interview else 'paper-trail'
    workflow_name = 'interview' if is_interview else 'paper_trail'
    outdir = args.out
    os.makedirs(outdir, exist_ok=True)
    _configure_logging(outdir)
    primers = getattr(args, 'primers', None) or _paper_trail_module.DEFAULT_PRIMERS
    inputs = getattr(args, 'input', None) or []
    sample_map = getattr(args, 'sample_map', None)
    if not inputs and not sample_map:
        raise SystemExit(f'[{command_name}] Provide --input and/or --sample-map.')
    primer_sequences = dict(_paper_trail_module.DEFAULT_PRIMER_SEQUENCES)
    for specification in getattr(args, 'primer_sequences', None) or []:
        if '=' not in str(specification):
            raise SystemExit(f'[{command_name}] --primer-sequence must use NAME=SEQUENCE.')
        name, sequence = str(specification).split('=', 1)
        sequence = ''.join(sequence.split()).upper()
        if not name.strip() or len(sequence) < 12:
            raise SystemExit(
                f'[{command_name}] --primer-sequence requires a name and at least 12 bases.'
            )
        primer_sequences[name.strip().upper()] = sequence
    manifest = RunManifest(outdir, workflow_name)
    dataset_hint = (
        _infer_qc_dataset_from_path(sample_map)
        or _infer_qc_dataset_from_path(getattr(args, 'mailroom_summary', None))
        or _infer_qc_dataset_from_path(outdir)
    )
    if dataset_hint:
        manifest.data['dataset'] = dataset_hint
        manifest.write()
    for source, role in (
        (sample_map, 'sample_map'),
        (getattr(args, 'mailroom_summary', None), 'mailroom_summary'),
        (getattr(args, 'screen_ref', None), 'screen_reference'),
        (getattr(args, 'screen_taxa', None), 'screen_taxonomy'),
    ):
        if source:
            manifest.add_input(source, role=role)
    if getattr(args, 'mailroom_status', 'PASS') != 'PASS':
        manifest.warn(
            'Interview proceeded after explicit acceptance of a REVIEW_REQUIRED Mailroom batch.'
        )
    for source in inputs:
        manifest.add_input(
            source,
            role='raw_read' if Path(source).is_file() else 'raw_read_directory',
        )
    manifest.add_stage('paper_trail', 'RUNNING', detail='base calls, chromatogram evidence, and trimming')
    policy = _paper_trail_module.DEFAULT_QC_POLICY

    def setting(name, default):
        value = getattr(args, name, None)
        return default if value is None else value

    try:
        outputs = _paper_trail_module.run_paper_trail(
            inputs,
            outdir,
            sample_map=sample_map,
            primers=primers,
            primer_sequences=primer_sequences,
            trim_primers=bool(getattr(args, 'trim_primers', True)),
            secondary_peak_ratio=float(setting('secondary_peak_ratio', policy.secondary_peak_ratio)),
            max_mixed_peak_percent=float(setting('max_mixed_peak_percent', policy.max_mixed_peak_percent)),
            mixed_peak_min_quality=int(setting('mixed_peak_min_quality', policy.mixed_peak_min_quality)),
            screen_ref=getattr(args, 'screen_ref', None),
            screen_taxa=getattr(args, 'screen_taxa', None),
            selected_read_ids=getattr(args, 'selected_read_ids', None),
            threads=int(setting('threads', 4)),
            min_quality=int(setting('min_quality', policy.min_quality)),
            min_length=int(setting('min_length', policy.min_final_length)),
            min_read_length=int(setting('min_read_length', policy.min_read_length)),
            min_mean_quality=float(setting('min_mean_quality', policy.min_mean_quality)),
            mask_quality=int(setting('mask_quality', policy.mask_quality)),
            max_read_expected_errors=float(setting('max_read_expected_errors', policy.max_read_expected_errors)),
            max_output_expected_errors=float(setting('max_output_expected_errors', policy.max_output_expected_errors)),
            warn_n_percent=float(setting('warn_n_percent', policy.warn_n_percent)),
            max_n_percent=float(setting('max_n_percent', policy.max_n_percent)),
            warn_internal_low_quality_run=int(setting('warn_internal_low_quality_run', policy.warn_internal_low_quality_run)),
            max_internal_low_quality_run=int(setting('max_internal_low_quality_run', policy.max_internal_low_quality_run)),
            max_conflict_density=float(setting('max_conflict_density', policy.max_conflict_density)),
            max_ambiguous_overlap_conflicts_without_review=int(setting('max_ambiguous_overlap_conflicts_without_review', policy.max_ambiguous_overlap_conflicts_without_review)),
            quality_difference=int(setting('quality_difference', policy.quality_difference)),
            allow_missing_quality=bool(getattr(args, 'allow_missing_quality', False)),
            min_overlap=int(setting('min_overlap', policy.min_overlap)),
            min_overlap_identity=float(setting('min_overlap_identity', policy.min_overlap_identity)),
            assemble=bool(getattr(args, 'assemble', True)),
            recursive=bool(getattr(args, 'recursive', True)),
            max_report_image_height=int(getattr(args, 'max_report_image_height', 2400) or 2400),
        )
        manifest.add_stage('paper_trail', 'COMPLETE')
        manifest.add_stage('merge_meeting', 'COMPLETE', detail='per-isolate assembly or best-read selection')
        for key in (
            'raw_fasta', 'trimmed_fasta', 'assembled_fasta', 'read_qc_tsv',
            'per_base_error_tsv', 'assembly_tsv', 'assembly_placements_tsv',
            'recommendations_tsv', 'marker_review_template_tsv', 'qc_policy_tsv',
            'failed_final_fasta',
            'failed_read_fasta', 'failed_manifest_tsv', 'failed_read_manifest_tsv',
            'failed_qc_guide', 'visual_manifest_tsv', 'summary',
        ):
            manifest.add_output(outputs[key], role=key)
        for key in ('read_error_pngs', 'chromatogram_pngs', 'assembly_pngs'):
            for path in outputs[key]:
                manifest.add_output(path, role=key)
        manifest.finish('COMPLETE')
    except Exception as exc:
        manifest.finish('FAILED', error=exc)
        raise
    logging.getLogger(__name__).info(
        "[%s] Final assembled FASTA: %s", command_name.upper(), outputs['assembled_fasta'],
    )
    completion = (
        'Standalone AB1 QC complete.' if is_interview
        else 'Paper Trail + Merge Meeting complete.'
    )
    next_step = (
        'No project database was read or changed. After resolving manual-review decisions, '
        'provide assembled.fasta and assembly_report.tsv to Assistant as a reviewed FASTA submission.'
        if is_interview else
        'Use in Performance Review:\n'
        f"  branchmanager performance-review --input {outputs['assembled_fasta']} --partner-metadata <metadata.tsv> ..."
    )
    print(
        f"[{command_name}] {completion}\n"
        f"  Assembled FASTA : {outputs['assembled_fasta']}\n"
        f"  Trimmed reads   : {outputs['trimmed_fasta']}\n"
        f"  Read QC         : {outputs['read_qc_tsv']}\n"
        f"  Per-base errors : {outputs['per_base_error_tsv']}\n"
        f"  Assembly report : {outputs['assembly_tsv']}\n"
        f"  Read placements : {outputs['assembly_placements_tsv']}\n\n"
        f"  Resequence list : {outputs['recommendations_tsv']}\n"
        f"  Review template : {outputs['marker_review_template_tsv']}\n"
        f"  QC policy       : {outputs['qc_policy_tsv']}\n"
        f"  Failed QC seqs  : {outputs['failed_qc_dir']}\n\n"
        f"  Failed isolates : {outputs['failed_manifest_tsv']}\n"
        f"  Failed reads    : {outputs['failed_read_manifest_tsv']}\n\n"
        f"  Failure guide   : {outputs['failed_qc_guide']}\n\n"
        f"  Read visuals    : {len(outputs['read_error_pngs'])} page(s)\n"
        f"  Chromatograms   : {len(outputs['chromatogram_pngs'])} page(s)\n"
        f"  Assembly visuals: {len(outputs['assembly_pngs'])} page(s)\n"
        f"  Visual manifest : {outputs['visual_manifest_tsv']}\n\n"
        f"{next_step}"
    )


def cmd_mailroom(args):
    from branchmanager.mailroom import prepare_ab1_map, prepare_pacbio_map

    manifest = RunManifest(args.out, 'mailroom')
    manifest.add_input(args.read_dir, role=f'{args.technology}_directory')
    manifest.add_input(args.manifest, role='supplier_metadata')
    try:
        if args.technology == 'sanger':
            result = prepare_ab1_map(
                args.read_dir, args.manifest, args.out, dataset=args.dataset,
                forward_primer=args.forward_primer, reverse_primer=args.reverse_primer,
                processing_mode=args.processing_mode, recursive=args.recursive,
            )
            for role in ('ab1_map', 'inventory', 'report', 'summary'):
                manifest.add_output(result[role], role=role)
        else:
            result = prepare_pacbio_map(
                args.read_dir, args.manifest, args.out, dataset=args.dataset,
                recursive=args.recursive, library_prep=args.pacbio_library_prep,
            )
            for role in ('pacbio_map', 'inventory', 'report', 'summary'):
                manifest.add_output(result[role], role=role)
            if result['status'] == 'PASS':
                from branchmanager.pipeline.kinnex import run_kinnex_import
                pacbio_outputs = run_kinnex_import(
                    result['pacbio_map'], Path(args.out) / 'pacbio_16s_import',
                    library_prep=args.pacbio_library_prep,
                    threads=args.pacbio_threads,
                    strand=args.pacbio_vsearch_strand,
                    query_cov=args.pacbio_vsearch_query_cov,
                    maxaccepts=args.pacbio_vsearch_maxaccepts,
                    maxrejects=args.pacbio_vsearch_maxrejects,
                )
                for role in ('fasta', 'marker_qc', 'read_qc', 'report'):
                    manifest.add_output(pacbio_outputs[role], role=f'pacbio_{role}')
                result['pacbio_outputs'] = pacbio_outputs
                if pacbio_outputs['accepted'] != pacbio_outputs['total']:
                    result['status'] = 'REVIEW_REQUIRED'
        manifest.add_stage(
            'mailroom', result['status'],
            detail=f"{result['mapped_reads']} mapped reads; {result['errors']} errors",
        )
        manifest.finish('COMPLETE' if result['status'] in {'PASS', 'REVIEW_REQUIRED'} else 'FAILED')
    except Exception as exc:
        manifest.finish('FAILED', error=exc)
        raise
    print(
        f"[mailroom] {result['status']}: {result.get('physical_ab1_files', result.get('physical_fastq_files', 0))} {args.technology} file(s), "
        f"{result['mapped_reads']} mapped read(s), {result['isolates']} metadata isolate(s), "
        f"{result['mapped_isolates']} with mapped reads, "
        f"{result['errors']} error(s), {result.get('unresolved_primers', 0)} unresolved primer(s).\n"
        f"  Batch map : {result.get('ab1_map', result.get('pacbio_map'))}\n"
        f"  Inventory : {result['inventory']}\n"
        f"  Report    : {result['report']}\n"
        f"  Summary   : {result['summary']}"
        + (f"\n  Representatives : {result['pacbio_outputs']['fasta']}\n  Marker QC       : {result['pacbio_outputs']['marker_qc']}" if result.get('pacbio_outputs') else '')
    )
    if result['status'] != 'PASS':
        raise SystemExit(2)


def cmd_replace_selected_markers(args):
    """QC new Sanger reads and atomically replace eligible selected markers."""
    from branchmanager.marker_provenance import load_marker_provenance, write_marker_qc_bridge
    from branchmanager.utils.fasta import read_fasta, write_fasta

    root = Path(args.out).expanduser().resolve()
    sanger_out = root / 'sanger_qc'
    manifest = RunManifest(root, 'replace-selected-markers')
    manifest.add_input(args.sample_map, role='replacement_sample_map')
    manifest.add_input(args.db, role='project_database')
    policy = _paper_trail_module.DEFAULT_QC_POLICY
    primer_sequences = dict(_paper_trail_module.DEFAULT_PRIMER_SEQUENCES)
    for specification in args.primer_sequences or []:
        if '=' not in str(specification):
            raise SystemExit('[replace-selected-markers] --primer-sequence must use NAME=SEQUENCE.')
        name, sequence = str(specification).split('=', 1)
        sequence = ''.join(sequence.split()).upper()
        if not name.strip() or len(sequence) < 12:
            raise SystemExit('[replace-selected-markers] --primer-sequence requires a name and at least 12 bases.')
        primer_sequences[name.strip().upper()] = sequence
    try:
        primers = args.primers or _paper_trail_module.DEFAULT_PRIMERS
        db = Database(args.db)
        db.initialise()
        read_meta, listed_files = _paper_trail_module._load_sample_map(args.sample_map, primers)
        previous = db.get_sequence_provenance()
        candidate_files = {}
        for path in listed_files:
            path = Path(path).expanduser().resolve()
            if path.suffix.lower() not in {'.ab1', '.abi'}:
                continue
            entry = next((read_meta.get(key) for key in _paper_trail_module._metadata_key_variants(path)
                          if read_meta.get(key)), {})
            sequence_id = str(entry.get('sequence_id') or '').strip()
            if sequence_id:
                candidate_files.setdefault(sequence_id, []).append((path, entry))

        changed = {}
        for sequence_id, files in candidate_files.items():
            old = previous.get(sequence_id, {})
            old_paths = set()
            for value in str(old.get('source_read_files') or '').split(';'):
                if value.strip():
                    try:
                        old_paths.add(str(Path(value.strip()).expanduser().resolve()))
                    except OSError:
                        old_paths.add(value.strip())
            try:
                old_hashes = json.loads(old.get('source_read_checksums_json') or '{}')
            except (TypeError, json.JSONDecodeError):
                old_hashes = {}
            for path, _entry in files:
                path_text = str(path)
                current_hash = sha256_file(path) if path.is_file() else ''
                old_hash = old_hashes.get(path_text) or old_hashes.get(path.name) or ''
                if not old_hash:
                    old_hash = next((value for old_path, value in old_hashes.items()
                                     if str(Path(old_path).expanduser().resolve()) == path_text), '')
                if path_text not in old_paths or (old_hash and current_hash and old_hash != current_hash) or (path_text in old_paths and not old_hash):
                    changed.setdefault(sequence_id, []).append((path, current_hash, 'new path' if path_text not in old_paths else ('content changed' if old_hash else 'legacy provenance; content not fingerprinted')))

        review_path = root / 'resequencing_file_review.tsv'
        reviewed_map = root / 'reviewed_replacement_sample_map.csv'
        accepted = set()
        review_rows = []
        print('[replace-selected-markers] Incoming AB1 files requiring review:')
        for sequence_id in sorted(changed):
            print(f'  {sequence_id}')
            for path, _fingerprint, reason in changed[sequence_id]:
                print(f'    {path} ({reason})')
            try:
                response = input(f'Re-examine isolate {sequence_id}? [y/N] ').strip().lower()
            except EOFError:
                response = ''
            answer = response in {'y', 'yes'}
            for path, fingerprint, reason in changed[sequence_id]:
                review_rows.append({'sequence_id': sequence_id, 'file': str(path), 'sha256': fingerprint,
                                    'change': reason, 'decision': 'YES' if answer else 'NO'})
            if answer:
                accepted.add(sequence_id)

        root.mkdir(parents=True, exist_ok=True)
        with open(review_path, 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=['sequence_id', 'file', 'sha256', 'change', 'decision'], delimiter='\t')
            writer.writeheader(); writer.writerows(review_rows)
        with open(reviewed_map, 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=['sequence_id', 'primer', 'direction', 'ab1_file'])
            writer.writeheader()
            for sequence_id in sorted(accepted):
                for path, entry in candidate_files[sequence_id]:
                    writer.writerow({'sequence_id': sequence_id, 'primer': entry.get('primer', ''),
                                     'direction': entry.get('direction', ''), 'ab1_file': str(path)})
        manifest.add_output(review_path, role='resequencing_file_review')
        manifest.add_output(reviewed_map, role='reviewed_replacement_sample_map')
        if not accepted:
            manifest.add_stage('resequencing_review', 'COMPLETE', detail='No isolates approved for re-examination')
            manifest.finish('COMPLETE')
            print(f'[replace-selected-markers] No isolates selected for re-examination. Review: {review_path}')
            return
        qc_outputs = _paper_trail_module.run_paper_trail(
            [], sanger_out, sample_map=str(reviewed_map),
            primers=primers,
            primer_sequences=primer_sequences, min_length=800,
            min_read_length=policy.min_read_length, min_quality=policy.min_quality,
            min_mean_quality=policy.min_mean_quality, mask_quality=policy.mask_quality,
            max_read_expected_errors=policy.max_read_expected_errors,
            max_output_expected_errors=policy.max_output_expected_errors,
            warn_n_percent=policy.warn_n_percent, max_n_percent=policy.max_n_percent,
            warn_internal_low_quality_run=policy.warn_internal_low_quality_run,
            max_internal_low_quality_run=policy.max_internal_low_quality_run,
            secondary_peak_ratio=policy.secondary_peak_ratio,
            max_mixed_peak_percent=policy.max_mixed_peak_percent,
            mixed_peak_min_quality=policy.mixed_peak_min_quality,
            max_conflict_density=policy.max_conflict_density,
            max_ambiguous_overlap_conflicts_without_review=policy.max_ambiguous_overlap_conflicts_without_review,
            quality_difference=policy.quality_difference, allow_missing_quality=args.allow_missing_quality,
            min_overlap=args.min_overlap, min_overlap_identity=args.min_overlap_identity,
            threads=args.threads,
        )
        assembled = Path(qc_outputs['assembled_fasta'])
        provenance, _qc_source = load_marker_provenance(
            assembled, qc_path=qc_outputs['assembly_tsv'], accept_unverified=False,
        )
        sequence_by_id = {sequence_id: sequence for sequence_id, sequence in read_fasta(assembled)}
        eligible, report_rows = [], []
        for row in provenance:
            sequence_id = str(row['sequence_id'])
            qc_class = str(row.get('marker_qc_class') or '').upper()
            if sequence_id not in sequence_by_id:
                db.record_marker_replacement_attempt(
                    sequence_id, 'FAILED_QC',
                    row.get('marker_qc_reasons') or qc_class or 'replacement sequence was not assembled',
                    source_manifest=str(root / 'run_manifest.json'),
                    source_fasta=str(assembled),
                    incoming_qc_class=qc_class,
                    incoming_qc_reasons=row.get('marker_qc_reasons') or '',
                    incoming_sequence=sequence_by_id.get(sequence_id, ''),
                )
                with db.connect() as conn:
                    old = conn.execute(
                        'SELECT length(sequence) FROM sequences WHERE id=?', (sequence_id,),
                    ).fetchone()
                old_provenance = db.get_sequence_provenance([sequence_id]).get(sequence_id, {})
                report_rows.append(_marker_replacement_report_row(
                    sequence_id,
                    {
                        'status': 'FAILED_QC',
                        'detail': row.get('marker_qc_reasons') or qc_class or 'replacement sequence was not assembled',
                        'retained_marker': 'EXISTING' if old else 'NONE',
                        'previous_length': old[0] if old else '',
                        'previous_qc_class': old_provenance.get('marker_qc_class') or '',
                        'previous_qc_reasons': old_provenance.get('marker_qc_reasons') or '',
                        'previous_read_ids': old_provenance.get('source_read_ids') or '',
                        'incoming_qc_class': qc_class,
                        'incoming_qc_reasons': row.get('marker_qc_reasons') or '',
                        'incoming_length': len(sequence_by_id.get(sequence_id, '')),
                        'incoming_read_ids': row.get('source_read_ids') or '',
                        'incoming_read_files': row.get('source_read_files') or '',
                    },
                ))
                continue
            eligible.append({**row, 'sequence': sequence_by_id[sequence_id]})
        applied = db.replace_selected_markers(
            eligible, source_manifest=str(root / 'run_manifest.json'), source_fasta=str(assembled),
        )
        applied_by_id = {row['sequence_id']: row for row in applied if row.get('sequence_id')}
        replaced_ids = {row['sequence_id'] for row in applied if row['status'] == 'REPLACED'}
        if replaced_ids:
            db.upsert_sequence_provenance([row for row in eligible if row['sequence_id'] in replaced_ids])
            for row in applied:
                if row['status'] == 'REPLACED':
                    db.update_isolate_status(
                        row['sequence_id'], 'MARKER_REPLACED', detail=row['detail'],
                        source_file=str(assembled),
                    )
        for row in eligible:
            result = applied_by_id.get(row['sequence_id'], {
                'sequence_id': row['sequence_id'], 'status': 'NOT_PROCESSED', 'detail': '',
            })
            report_rows.append(_marker_replacement_report_row(row['sequence_id'], result))
        report_path = root / 'marker_replacement_report.tsv'
        root.mkdir(parents=True, exist_ok=True)
        _write_marker_replacement_report(report_path, report_rows)
        replacement_fasta = root / 'replaced_markers.fasta'
        write_fasta(
            [(sequence_id, sequence_by_id[sequence_id]) for sequence_id in sorted(replaced_ids)],
            replacement_fasta,
        )
        provenance_path = write_marker_qc_bridge(root / 'marker_qc_provenance.tsv', provenance)
        for path, role in (
            (qc_outputs['assembled_fasta'], 'sanger_assembled_fasta'),
            (qc_outputs['assembly_tsv'], 'sanger_assembly_report'),
            (report_path, 'marker_replacement_report'), (replacement_fasta, 'replaced_markers'),
            (provenance_path, 'marker_qc_provenance'),
            (review_path, 'resequencing_file_review'), (reviewed_map, 'reviewed_replacement_sample_map'),
        ):
            manifest.add_output(path, role=role)
        manifest.add_stage('paper_trail', 'COMPLETE', detail='800 bp minimum Sanger replacement QC')
        manifest.add_stage('marker_replacement', 'COMPLETE', detail=f'{len(replaced_ids)} markers replaced')
        manifest.finish('COMPLETE')
    except Exception as exc:
        manifest.finish('FAILED', error=exc)
        raise
    print(
        f"[replace-selected-markers] {len(replaced_ids)} marker(s) replaced.\n"
        f"  Replacement report : {report_path}\n"
        f"  Replacement FASTA  : {replacement_fasta}\n"
        '  Re-run Performance Review before using taxonomy, novelty, or tree results for these isolates.'
    )


def cmd_onboarding(args):
    from branchmanager.onboarding import validate_submission, write_onboarding_outputs
    from branchmanager.pipeline.paper_trail import DEFAULT_PRIMERS

    manifest = RunManifest(args.out, 'onboarding')
    sample_map = getattr(args, 'sample_map', None)
    fasta = getattr(args, 'fasta', None)
    manifest.add_input(sample_map or fasta, role='sample_map' if sample_map else 'marker_fasta')
    if args.partner_metadata:
        manifest.add_input(args.partner_metadata, role='partner_metadata')
    try:
        result = validate_submission(
            sample_map,
            fasta=fasta,
            partner_metadata=args.partner_metadata,
            read_dir=args.read_dir,
            primers=args.primers or DEFAULT_PRIMERS,
            expected_partner_id=getattr(args, 'partner_id', None),
            dataset=getattr(args, 'dataset', None),
        )
        outputs = write_onboarding_outputs(args.out, result)
        for role, path in outputs.items():
            manifest.add_output(path, role=role)
        manifest.add_stage('onboarding', result['status'], detail=f"{result['isolates']} isolates; {result['errors']} errors")
        manifest.finish('COMPLETE' if result['status'] == 'PASS' else 'FAILED')
    except Exception as exc:
        manifest.finish('FAILED', error=exc)
        raise
    print(
        f"[onboarding] {result['status']}: {result['isolates']} isolate(s), "
        f"input={result['input_type']}, {result['read_files']} read file(s), {result['errors']} error(s).\n"
        f"  Normalised map : {outputs['normalised']}\n"
        f"  Validation     : {outputs['report']}"
        + (f"\n  Normalised FASTA: {outputs['fasta']}" if outputs.get('fasta') else '')
    )
    if result['status'] != 'PASS':
        raise SystemExit(2)


def cmd_pacbio_16s_import(args):
    from branchmanager.pipeline.kinnex import run_kinnex_import

    manifest = RunManifest(args.out, 'pacbio-16s-import')
    manifest.add_input(args.pacbio_map, role='pacbio_fastq_map')
    try:
        outputs = run_kinnex_import(
            args.pacbio_map, args.out, min_read_length=args.min_read_length,
            max_read_length=args.max_read_length, min_mean_quality=args.min_mean_quality,
            cluster_identity=args.cluster_identity, min_reads=args.min_reads,
            min_dominant_fraction=args.min_dominant_fraction, library_prep=args.library_prep,
            threads=args.threads, strand=args.vsearch_strand,
            query_cov=args.vsearch_query_cov, maxaccepts=args.vsearch_maxaccepts,
            maxrejects=args.vsearch_maxrejects,
        )
        for role in ('fasta', 'marker_qc', 'report', 'read_qc'):
            manifest.add_output(outputs[role], role=f'pacbio_{role}')
        manifest.add_stage('pacbio_16s_import', 'COMPLETE', detail=f"{outputs['accepted']}/{outputs['total']} isolates accepted")
        manifest.finish('COMPLETE')
    except Exception as exc:
        manifest.finish('FAILED', error=exc)
        raise
    print(f"[pacbio-16s-import] {outputs['accepted']}/{outputs['total']} isolate(s) accepted.\n"
          f"  Representatives : {outputs['fasta']}\n  Marker QC       : {outputs['marker_qc']}\n  QC report       : {outputs['report']}")


def _cmd_project_import(args, *, genomes: bool):
    import fcntl
    import sqlite3
    from branchmanager.project_state import import_genome_results, import_status_updates, write_import_report

    workflow = 'records_update' if genomes else 'status_meeting'
    manifest = RunManifest(args.out, workflow)
    manifest.add_input(args.input, role='genome_results' if genomes else 'status_updates')
    original = Path(args.db).expanduser().resolve()
    if not original.is_file():
        raise SystemExit(f'[records-update] Project database does not exist: {original}')
    manifest.add_input(original, role='project_database_before')
    staged_path = Path(args.out).resolve() / f'.{workflow}_project.sqlite'
    lock_path = original.with_name(original.name + '.lock')
    with open(lock_path, 'a+') as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        try:
            staged_path.unlink(missing_ok=True)
            with sqlite3.connect(str(original)) as source, sqlite3.connect(str(staged_path)) as target:
                source.backup(target)
            db = Database(str(staged_path))
            db.initialise()
            rows = import_genome_results(
                db, args.input,
                min_completeness=float(getattr(args, 'min_completeness', 90.0)),
                max_contamination=float(getattr(args, 'max_contamination', 5.0)),
            ) if genomes else import_status_updates(db, args.input)
            report_name = 'records_update_report.tsv' if genomes else 'status_meeting_report.tsv'
            report = write_import_report(Path(args.out) / report_name, rows)
            manifest.add_output(report, role='import_report')
            rejected = sum(row.get('result') == 'REJECTED' for row in rows)
            manifest.add_stage(workflow, 'COMPLETE' if rejected == 0 else 'COMPLETE_WITH_WARNINGS', detail=f'{len(rows) - rejected} imported; {rejected} rejected')
            manifest.finish('COMPLETE' if rejected == 0 else 'COMPLETE_WITH_WARNINGS')
            db.record_project_run(
                f'{workflow.replace("_", "-")}:{manifest.data["started_at"]}', workflow, manifest.data['status'],
                manifest_path=str(manifest.json_path), started_at=manifest.data['started_at'],
                completed_at=manifest.data['completed_at'],
            )
            publish_tmp = original.with_name(original.name + '.publishing')
            publish_tmp.unlink(missing_ok=True)
            with sqlite3.connect(str(staged_path)) as source, sqlite3.connect(str(publish_tmp)) as target:
                source.backup(target)
                if target.execute('PRAGMA integrity_check').fetchone()[0] != 'ok':
                    raise RuntimeError('staged Records Update database failed integrity_check')
            os.replace(publish_tmp, original)
            staged_path.unlink(missing_ok=True)
        except BaseException as exc:
            manifest.finish('FAILED', error=exc)
            raise
        finally:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
    command_label = 'records-update' if genomes else 'status-meeting'
    print(f'[{command_label}] {len(rows) - rejected} row(s) imported; {rejected} rejected.\n  Report: {report}')


def cmd_add_external_genomes(args):
    """Add genome evidence acquired outside BranchManager's selection cycle."""
    import fcntl
    import sqlite3
    from branchmanager.project_state import import_external_genomes, write_import_report

    root = Path(args.out).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    manifest = RunManifest(root, 'external_genome_import')
    manifest.add_input(args.input, role='external_genomes')
    if args.marker_fasta:
        manifest.add_input(args.marker_fasta, role='optional_marker_sequences')
    original = Path(args.db).expanduser().resolve()
    if not original.is_file():
        raise SystemExit(f'[external-genomes] Project database does not exist: {original}')
    staged = root / '.external_genomes_project.sqlite'
    with open(original.with_name(original.name + '.lock'), 'a+') as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        try:
            staged.unlink(missing_ok=True)
            with sqlite3.connect(str(original)) as source, sqlite3.connect(str(staged)) as target:
                source.backup(target)
            db = Database(str(staged))
            db.initialise()
            results = import_external_genomes(
                db, args.input, marker_fasta=args.marker_fasta,
                default_dataset=args.dataset, default_tier=args.baseline_tier,
            )
            report = write_import_report(root / 'external_genome_import_report.tsv', results)
            manifest.add_output(report, role='import_report')
            rejected = sum(row.get('result') == 'REJECTED' for row in results)
            manifest.add_stage('external_genome_import', 'COMPLETE' if rejected == 0 else 'COMPLETE_WITH_WARNINGS',
                               detail=f'{len(results) - rejected} imported; {rejected} rejected')
            manifest.finish('COMPLETE' if rejected == 0 else 'COMPLETE_WITH_WARNINGS')
            db.record_project_run(
                f'external-genome-import:{manifest.data["started_at"]}', 'external_genome_import',
                manifest.data['status'], manifest_path=str(manifest.json_path),
                started_at=manifest.data['started_at'], completed_at=manifest.data['completed_at'],
            )
            publish = original.with_name(original.name + '.publishing')
            publish.unlink(missing_ok=True)
            with sqlite3.connect(str(staged)) as source, sqlite3.connect(str(publish)) as target:
                source.backup(target)
            os.replace(publish, original)
            staged.unlink(missing_ok=True)
            manifest.add_output(original, role='project_database')
        except BaseException as exc:
            manifest.finish('FAILED', error=exc)
            raise
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
    print(f'[external-genomes] {len(results) - rejected} genome(s) imported; {rejected} rejected.\n  Report: {report}')


def cmd_add_maldi_identifications(args):
    """Register markerless partner isolates with source-labeled MALDI context."""
    import fcntl
    import sqlite3
    from branchmanager.project_state import import_maldi_identifications, write_maldi_import_report

    root = Path(args.out).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    manifest = RunManifest(root, 'maldi_identification_import')
    manifest.add_input(args.input, role='maldi_identifications')
    original = Path(args.db).expanduser().resolve()
    if not original.is_file():
        raise SystemExit(f'[maldi-identifications] Project database does not exist: {original}')
    staged = root / '.maldi_identifications_project.sqlite'
    with open(original.with_name(original.name + '.lock'), 'a+') as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        try:
            staged.unlink(missing_ok=True)
            with sqlite3.connect(str(original)) as source, sqlite3.connect(str(staged)) as target:
                source.backup(target)
            db = Database(str(staged))
            db.initialise()
            results = import_maldi_identifications(
                db, args.input, default_dataset=args.dataset,
            )
            report = write_maldi_import_report(root / 'maldi_identification_import_report.tsv', results)
            manifest.add_output(report, role='import_report')
            rejected = sum(row.get('result') == 'REJECTED' for row in results)
            imported = sum(row.get('result') == 'IMPORTED' for row in results)
            manifest.add_stage(
                'maldi_identification_import',
                'COMPLETE' if rejected == 0 else 'COMPLETE_WITH_WARNINGS',
                detail=f'{imported} imported; {rejected} rejected',
            )
            manifest.finish('COMPLETE' if rejected == 0 else 'COMPLETE_WITH_WARNINGS')
            db.record_project_run(
                f'maldi-identification-import:{manifest.data["started_at"]}',
                'maldi_identification_import', manifest.data['status'],
                manifest_path=str(manifest.json_path), started_at=manifest.data['started_at'],
                completed_at=manifest.data['completed_at'],
            )
            publish = original.with_name(original.name + '.publishing')
            publish.unlink(missing_ok=True)
            with sqlite3.connect(str(staged)) as source, sqlite3.connect(str(publish)) as target:
                source.backup(target)
            os.replace(publish, original)
            staged.unlink(missing_ok=True)
            manifest.add_output(original, role='project_database')
        except BaseException as exc:
            manifest.finish('FAILED', error=exc)
            raise
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
    print(f'[maldi-identifications] {imported} identification(s) imported; {rejected} rejected.\n  Report: {report}')


def cmd_import_partner_workbook(args):
    """Import isolate metadata and MALDI observations from a partner XLSM."""
    import fcntl
    import sqlite3
    from branchmanager.project_state import (
        import_partner_workbook, write_partner_workbook_import_report,
        import_genome_only_partner_workbook, write_genome_only_partner_import_report,
    )

    root = Path(args.out).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    manifest = RunManifest(root, 'partner_workbook_import')
    manifest.add_input(args.input, role='partner_metadata_workbook')
    original = Path(args.db).expanduser().resolve()
    if not original.is_file():
        raise SystemExit(f'[partner-workbook] Project database does not exist: {original}')
    staged = root / '.partner_workbook_project.sqlite'
    with open(original.with_name(original.name + '.lock'), 'a+') as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        try:
            staged.unlink(missing_ok=True)
            with sqlite3.connect(str(original)) as source, sqlite3.connect(str(staged)) as target:
                source.backup(target)
            db = Database(str(staged))
            db.initialise()
            if args.genome_only:
                results = import_genome_only_partner_workbook(
                    db, args.input, sheet_name=args.sheet, dataset=args.dataset,
                )
                report = write_genome_only_partner_import_report(
                    root / 'genome_only_partner_import_report.tsv', results,
                )
                inserted = sum(row.get('result') == 'INSERTED' for row in results)
                conflicts = sum(row.get('result') == 'CONFLICT' for row in results)
                rejected = sum(row.get('result') == 'REJECTED' for row in results)
                mode_label = 'genome-only partner registration'
                detail = (f'{len(results)} workbook row(s); {inserted} new isolate(s); '
                          f'{conflicts} conflict(s); {rejected} rejected')
            else:
                results = import_partner_workbook(
                    db, args.input, sheet_name=args.sheet, dataset=args.dataset,
                )
                report = write_partner_workbook_import_report(
                    root / 'partner_workbook_import_report.tsv', results,
                )
                inserted = sum(row.get('metadata_result') == 'INSERTED_ISOLATE' for row in results)
                maldi_imported = sum(row.get('maldi_result') == 'SOFT_TAXON_CONTEXT_IMPORTED' for row in results)
                needs_review = sum(row.get('maldi_result') in {
                    'MALDI_RESULT_NEEDS_REVIEW', 'MALDI_SELECTED_NO_RESULT',
                    'MALDI_FIELDS_WITHOUT_SELECTION',
                } for row in results)
                conflicts = 0
                rejected = 0
                mode_label = 'partner workbook import'
                detail = (f'{len(results)} isolate rows; {inserted} inserted; '
                          f'{maldi_imported} MALDI context records; {needs_review} MALDI row(s) need review')
            manifest.add_output(report, role='partner_workbook_import_report')
            status = 'COMPLETE_WITH_WARNINGS' if conflicts or rejected or (not args.genome_only and needs_review) else 'COMPLETE'
            manifest.add_stage(mode_label, status, detail=detail)
            manifest.finish(status)
            db.record_project_run(
                f'partner-workbook-import:{manifest.data["started_at"]}',
                'partner_workbook_import', status,
                manifest_path=str(manifest.json_path), started_at=manifest.data['started_at'],
                completed_at=manifest.data['completed_at'],
            )
            publish = original.with_name(original.name + '.publishing')
            publish.unlink(missing_ok=True)
            with sqlite3.connect(str(staged)) as source, sqlite3.connect(str(publish)) as target:
                source.backup(target)
                if target.execute('PRAGMA integrity_check').fetchone()[0] != 'ok':
                    raise RuntimeError('partner-workbook staged database failed integrity_check')
            os.replace(publish, original)
            staged.unlink(missing_ok=True)
            manifest.add_output(original, role='project_database')
        except BaseException as exc:
            manifest.finish('FAILED', error=exc)
            raise
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
    if args.genome_only:
        print(
            f'[partner-workbook] {len(results)} row(s); {inserted} new genome-only isolate(s); '
            f'{conflicts} conflict(s); {rejected} rejected.\n  Report: {report}'
        )
    else:
        print(
            f'[partner-workbook] {len(results)} isolate row(s); {inserted} inserted; '
            f'{maldi_imported} soft MALDI context record(s); {needs_review} need review.\n'
            f'  Report: {report}'
        )


def cmd_exit_interview(args):
    """Preview or atomically apply an audited sequence withdrawal."""
    import fcntl
    import sqlite3
    from branchmanager.personnel import (
        load_exit_requests, write_departing_fasta, write_exit_interview_report,
    )

    output = Path(args.out).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    original = Path(args.db).expanduser().resolve()
    if not original.is_file():
        raise SystemExit(f'[exit-interview] Project database does not exist: {original}')
    requests = load_exit_requests(
        args.sequence_id, input_path=args.input, default_reason=args.reason,
    )
    manifest = RunManifest(output, 'exit_interview')
    manifest.add_input(original, role='project_database_before')
    if args.input:
        manifest.add_input(args.input, role='exit_interview_requests')

    if not args.apply:
        db = Database(str(original))
        rows = db.plan_sequence_removals(
            requests,
            allow_baseline=bool(args.allow_baseline),
            allow_genome_records=bool(args.allow_genome_records),
        )
        report = write_exit_interview_report(output / 'exit_interview_plan.tsv', rows)
        manifest.add_output(report, role='exit_interview_plan')
        ready = sum(row.get('status') == 'READY' for row in rows)
        blocked = len(rows) - ready
        status = 'COMPLETE' if blocked == 0 else 'COMPLETE_WITH_WARNINGS'
        manifest.add_stage('exit_interview', status, detail=f'{ready} ready; {blocked} blocked')
        manifest.finish(status)
        print(
            f'[exit-interview] Preview only: {ready} sequence(s) ready; {blocked} blocked.\n'
            f'  Plan: {report}\n'
            '  Review the plan, then repeat the command with --apply.'
        )
        return

    staged_path = output / '.exit_interview_project.sqlite'
    lock_path = original.with_name(original.name + '.lock')
    blocked_rows = []
    report = ''
    removed_count = 0
    with open(lock_path, 'a+') as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        try:
            staged_path.unlink(missing_ok=True)
            with sqlite3.connect(str(original)) as source, sqlite3.connect(str(staged_path)) as target:
                source.backup(target)
            db = Database(str(staged_path))
            db.initialise()
            rows = db.plan_sequence_removals(
                requests,
                allow_baseline=bool(args.allow_baseline),
                allow_genome_records=bool(args.allow_genome_records),
            )
            blocked_rows = [row for row in rows if row.get('status') != 'READY']
            if blocked_rows:
                report = write_exit_interview_report(output / 'exit_interview_plan.tsv', rows)
                manifest.add_output(report, role='exit_interview_plan')
                manifest.add_stage(
                    'exit_interview', 'COMPLETE_WITH_WARNINGS',
                    detail=f'0 removed; {len(blocked_rows)} blocked',
                )
                manifest.finish('COMPLETE_WITH_WARNINGS')
            else:
                archive = write_departing_fasta(output / 'departing_sequences.fasta', rows)
                manifest.add_output(archive, role='departing_sequences')
                if args.backup:
                    backup = output / 'project_before_exit_interview.sqlite'
                    backup.unlink(missing_ok=True)
                    with sqlite3.connect(str(original)) as source, sqlite3.connect(str(backup)) as target:
                        source.backup(target)
                    manifest.add_output(backup, role='project_database_backup')
                db.apply_sequence_removals(
                    rows, source_request=str(Path(args.input).resolve()) if args.input else 'command_line',
                )
                report = write_exit_interview_report(output / 'exit_interview_report.tsv', rows)
                manifest.add_output(report, role='exit_interview_report')
                with db.connect() as conn:
                    integrity = conn.execute('PRAGMA integrity_check').fetchone()[0]
                    remaining = conn.execute(
                        'SELECT COUNT(*) FROM sequences WHERE id IN ('
                        + ','.join('?' for _ in rows) + ')',
                        tuple(row['sequence_id'] for row in rows),
                    ).fetchone()[0]
                if integrity != 'ok':
                    raise RuntimeError(f'staged Exit Interview database failed integrity_check: {integrity}')
                if remaining:
                    raise RuntimeError(f'{remaining} planned sequence(s) remained after Exit Interview')
                removed_count = len(rows)
                manifest.add_stage('exit_interview', 'COMPLETE', detail=f'{removed_count} removed')
                manifest.finish('COMPLETE')
                db.record_project_run(
                    f'exit-interview:{manifest.data["started_at"]}', 'exit_interview', 'COMPLETE',
                    manifest_path=str(manifest.json_path), started_at=manifest.data['started_at'],
                    completed_at=manifest.data['completed_at'],
                )
                publish_tmp = original.with_name(original.name + '.publishing')
                publish_tmp.unlink(missing_ok=True)
                with sqlite3.connect(str(staged_path)) as source, sqlite3.connect(str(publish_tmp)) as target:
                    source.backup(target)
                    if target.execute('PRAGMA integrity_check').fetchone()[0] != 'ok':
                        raise RuntimeError('published Exit Interview database failed integrity_check')
                os.replace(publish_tmp, original)
                staged_path.unlink(missing_ok=True)
        except BaseException as exc:
            if not blocked_rows:
                manifest.finish('FAILED', error=exc)
            raise
        finally:
            staged_path.unlink(missing_ok=True)
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
    if blocked_rows:
        print(
            f'[exit-interview] No changes applied; {len(blocked_rows)} request(s) are blocked.\n'
            f'  Plan: {report}'
        )
        raise SystemExit(2)
    print(
        f'[exit-interview] Removed {removed_count} sequence(s) from active project state.\n'
        f'  Audit report: {report}\n'
        f'  Recovery FASTA: {output / "departing_sequences.fasta"}'
    )


def cmd_marker_evidence_migrate(args):
    """Create the v1 marker-evidence registry, conservatively backfill legacy rows, and report."""
    import json
    import shutil
    from branchmanager.marker_evidence import (
        ensure_marker_evidence_schema, backfill_legacy_marker_evidence,
        sync_registry_with_sequences, marker_evidence_migration_report,
    )
    source = Path(args.db).expanduser().resolve()
    output = Path(args.out).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    if source.is_file():
        backup = output / (source.name + '.pre_marker_evidence_migration.bak')
        if not backup.exists():
            shutil.copy2(source, backup)
    db = Database(args.db)
    db.initialise()
    ensure_marker_evidence_schema(db)
    counts = backfill_legacy_marker_evidence(db, force=bool(args.force))
    counts['synced_from_sequences'] = sync_registry_with_sequences(db)
    report = output / 'marker_evidence_backfill.tsv'
    with open(report, 'w') as handle:
        handle.write('Category\tCount\n')
        for key, value in counts.items():
            handle.write(f'{key}\t{value}\n')
    summary = marker_evidence_migration_report(db)
    (output / 'marker_evidence_migration_report.json').write_text(json.dumps(summary, indent=2, sort_keys=True, default=str))
    with open(output / 'marker_evidence_migration_report.tsv', 'w') as handle:
        handle.write('Category\tCount\n')
        for key, value in summary.items():
            if isinstance(value, dict):
                for k2, v2 in value.items():
                    handle.write(f'{key}:{k2}\t{v2}\n')
            else:
                handle.write(f'{key}\t{value}\n')
    print(f'[marker-evidence-migrate] Evidence registry ready. Report: {report}')


def _run_on_staged_database(args, label, work):
    """Run ``work(db, root, manifest)`` on a staged DB copy; publish only after integrity_check."""
    import fcntl
    import sqlite3

    root = Path(args.out).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    original = Path(args.db).expanduser().resolve()
    if not original.is_file():
        raise SystemExit(f'[{label}] Project database does not exist: {original}')
    manifest = RunManifest(root, label.replace('-', '_'))
    staged = root / f'.{label}_project.sqlite'
    with open(original.with_name(original.name + '.lock'), 'a+') as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        try:
            staged.unlink(missing_ok=True)
            with sqlite3.connect(str(original)) as src, sqlite3.connect(str(staged)) as dst:
                src.backup(dst)
            db = Database(str(staged))
            db.initialise()
            result = work(db, root, manifest)
            with sqlite3.connect(str(staged)) as conn:
                check = conn.execute('PRAGMA integrity_check').fetchone()[0]
            if check != 'ok':
                raise RuntimeError(f'Staged database failed integrity_check: {check}')
            if not getattr(args, 'dry_run', False):
                publish = original.with_name(original.name + '.publishing')
                publish.unlink(missing_ok=True)
                with sqlite3.connect(str(staged)) as src, sqlite3.connect(str(publish)) as dst:
                    src.backup(dst)
                os.replace(publish, original)
            staged.unlink(missing_ok=True)
            manifest.finish('COMPLETE')
            return result
        except BaseException as exc:
            manifest.finish('FAILED', error=exc)
            raise
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


def cmd_import_genome_results(args):
    from branchmanager.genome_evidence_import import import_precomputed_genome_results

    def work(db, root, manifest):
        manifest.add_input(args.results_tsv, role='precomputed_genome_results')
        result = import_precomputed_genome_results(
            db, args.results_tsv, dry_run=args.dry_run,
        )
        report_path = root / 'genome_results_import_summary.json'
        report_path.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
        manifest.add_output(report_path, role='genome_results_import_summary')
        manifest.add_stage(
            'precomputed_genome_results_import', 'COMPLETE',
            detail=f"{result['records']} genome record(s) imported across {len(result['datasets'])} dataset(s)",
        )
        return result

    result = _run_on_staged_database(args, 'import-genome-results', work)
    print(json.dumps(result, indent=2))


def cmd_import_genome_evidence(args):
    from branchmanager.genome_evidence_import import (
        import_partner_genome_evidence, write_genome_evidence_report,
    )

    def work(db, root, manifest):
        for path, role in (
            (args.workbook, 'partner_workbook'),
            (args.gtdbtk, 'gtdbtk_summary'), (args.checkm2, 'checkm2_quality_report'),
        ):
            manifest.add_input(path, role=role)
        genome_dir = Path(args.genome_dir).expanduser().resolve()
        for fasta in sorted(genome_dir.iterdir()):
            if fasta.is_file() and re.search(r'\.(?:fa|fna|fasta)(?:\.gz)?$', fasta.name, re.I):
                manifest.add_input(fasta, role='genome_fasta')
        rows = import_partner_genome_evidence(
            db, workbook=args.workbook, genome_dir=args.genome_dir,
            gtdbtk=args.gtdbtk, checkm2=args.checkm2, dataset=args.dataset,
            sheet_name=args.sheet, min_completeness=args.min_completeness,
            max_contamination=args.max_contamination,
            gtdbtk_version=args.gtdbtk_version, checkm2_version=args.checkm2_version,
            dry_run=args.dry_run,
        )
        report = write_genome_evidence_report(root / 'genome_evidence_import_report.tsv', rows)
        manifest.add_output(report, role='genome_evidence_import_report')
        imported = sum(r.get('result') in {'IMPORTED', 'DRY_RUN_READY'} for r in rows)
        incomplete = sum(r.get('result') == 'INCOMPLETE' for r in rows)
        manifest.add_stage('partner_genome_evidence_import', 'COMPLETE_WITH_WARNINGS' if incomplete else 'COMPLETE',
                           detail=f'{imported} genome(s) ready/imported; {incomplete} incomplete')
        return rows

    rows = _run_on_staged_database(args, 'import-genome-evidence', work)
    imported = sum(r.get('result') in {'IMPORTED', 'DRY_RUN_READY'} for r in rows)
    incomplete = sum(r.get('result') == 'INCOMPLETE' for r in rows)
    print(f'[import-genome-evidence] {imported} genome(s) ready/imported; {incomplete} incomplete. '
          f"Report: {Path(args.out).expanduser().resolve() / 'genome_evidence_import_report.tsv'}")


def cmd_import_genome_16s(args):
    from branchmanager.genome_16s_import import import_genome_16s

    def work(db, root, manifest):
        gdir = Path(args.genomic_dir).expanduser().resolve()
        for name in ('genome_16s_summary.tsv', 'best_16s_per_genome.fna', 'run_manifest.json'):
            if (gdir / name).is_file():
                manifest.add_input(gdir / name, role=f'genomic_16s:{name}')
        gm = gdir / 'run_manifest.json'
        if gm.is_file():
            try:
                gdata = json.loads(gm.read_text())
                manifest.set_reference_identity(
                    genomic_16s_manifest=str(gm),
                    reference=json.dumps({k: v for k, v in gdata.items() if 'sha256' in k.lower() or 'reference' in k.lower() or 'gtdb' in k.lower()}, default=str)[:2000],
                )
            except (OSError, ValueError):
                pass
        if args.genome_metadata:
            manifest.add_input(args.genome_metadata, role='genome_metadata')
        return import_genome_16s(
            db, gdir, root, genome_metadata=args.genome_metadata,
            dataset=args.dataset, partner_id=args.partner_id,
            genome_id_isolate_id=args.genome_id_isolate_id,
        )

    results = _run_on_staged_database(args, 'import-genome-16s', work)
    warnings = sum(bool(r.get('warning')) for r in results)
    print(f'[import-genome-16s] {len(results)} genome(s) processed; {warnings} GTDB-Tk/16S disagreement warning(s).')


def cmd_marker_review_apply(args):
    from branchmanager.marker_evidence import apply_marker_review_decisions

    def work(db, root, manifest):
        manifest.add_input(args.decisions, role='marker_review_decisions')
        results = apply_marker_review_decisions(db, args.decisions, dry_run=bool(args.dry_run))
        with open(root / 'marker_review_apply_report.tsv', 'w', newline='') as handle:
            fields = sorted({k for r in results for k in r}) or ['sequence_id', 'result']
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
            writer.writeheader()
            writer.writerows(results)
        errors = [r for r in results if r.get('result') == 'ERROR']
        if errors:
            raise SystemExit(
                f'[marker-review-apply] {len(errors)} invalid decision row(s); nothing was applied. '
                f'See {root / "marker_review_apply_report.tsv"}'
            )
        return results

    results = _run_on_staged_database(args, 'marker-review-apply', work)
    applied = sum(r.get('result') == 'APPLIED' for r in results)
    print(f'[marker-review-apply] {applied} decision(s) applied of {len(results)} row(s).')


def cmd_register_marker_evidence(args):
    """Register supplementary marker observations without replacing the active isolate marker."""
    from branchmanager.marker_evidence import (
        ensure_marker_evidence_schema, backfill_legacy_marker_evidence,
        register_supplemental_marker_fasta,
    )
    db = Database(args.db)
    db.initialise()
    ensure_marker_evidence_schema(db)
    backfill_legacy_marker_evidence(db)
    results = register_supplemental_marker_fasta(
        db, args.fasta, args.metadata,
        default_source_type=args.source_type,
        default_evidence_tier=args.evidence_tier,
    )
    output = Path(args.out).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    report = output / 'marker_evidence_import.tsv'
    fields = sorted({key for row in results for key in row}) or ['sequence_id', 'result']
    with open(report, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(results)
    imported = sum(row.get('result') == 'IMPORTED' for row in results)
    existing = sum(row.get('result') == 'ALREADY_PRESENT' for row in results)
    print(f'[register-marker-evidence] {imported} imported; {existing} already present. Report: {report}')


def cmd_annual_report(args):
    from branchmanager.reporting import write_annual_report

    db = Database(args.db)
    db.initialise()
    manifest = RunManifest(args.out, 'annual_report')
    manifest.add_input(args.db, role='project_database')
    try:
        outputs = write_annual_report(db, args.out)
        for role, path in outputs.items():
            manifest.add_output(path, role=role)
        manifest.add_stage('annual_report', 'COMPLETE')
        manifest.finish('COMPLETE')
    except Exception as exc:
        manifest.finish('FAILED', error=exc)
        raise
    print(f"[annual-report] Cumulative project report: {outputs['html']}")


def cmd_it_desk(args):
    from branchmanager.it_desk import run_it_desk_checks, write_it_desk_report

    rows = run_it_desk_checks(
        db_path=args.db, references=args.references, output_dir=args.out,
        tree_method=args.tree_method,
    )
    outputs = write_it_desk_report(args.out, rows)
    for row in rows:
        print(f"[{row['status']}] {row['check']}: {row['detail']}")
    print(f"[it-desk] {outputs['status']}; report: {outputs['tsv']}")
    if args.strict and outputs['status'] != 'PASS':
        raise SystemExit(2)


_COHORT_INPUT_TYPES = {'sanger', 'partner_fasta', 'pacbio_hifi', 'genome_16s'}
_COHORT_INPUT_ALIASES = {
    'sanger_ab1': 'sanger', 'ab1': 'sanger', 'fasta': 'partner_fasta', 'marker_fasta': 'partner_fasta',
    'pacbio': 'pacbio_hifi', 'genome': 'genome_16s', 'genomic_16s': 'genome_16s',
    'barrnap': 'genome_16s', 'genome_barrnap': 'genome_16s',
}
_COHORT_MARKER_SOURCE = {'sanger': 'SANGER_AB1', 'partner_fasta': 'PARTNER_FASTA', 'pacbio_hifi': 'PACBIO_HIFI'}
_COHORT_REVIEW_FIELDS = [
    'sequence_id', 'decision', 'reviewer', 'reviewed_at', 'notes', 'review_type',
    'eligible_for_manual_acceptance', 'candidate_fasta', 'candidate_source',
    'candidate_read_id', 'candidate_length', 'candidate_n_percent', 'qc_class', 'qc_reasons',
    'taxonomy_concordance', 'read_taxonomy_assignments', 'output_length', 'n_percent',
    'read_level_evidence', 'quality_status', 'manual_acceptance_basis', 'manual_acceptance_reason',
    'minimum_acceptable_length', 'length_deficit_bp', 'allowed_decisions',
    'fallback_read_options', 'recommended_read_id', 'selected_read_id',
    'read_selection_source', 'fallback_reason', 'processing_mode', 'assembly_status', 'output_type',
    'partner_id', 'dataset', 'marker_source_type',
]
_COHORT_ALLOWED = {
    'WARNING_REVIEW': 'APPROVE|REJECT|RESEQUENCE',
    'FAILED_QC_REVIEW': 'APPROVE_FAIL_QC|REJECT|RESEQUENCE',
    'UNVERIFIED_QC_REVIEW': 'APPROVE_UNVERIFIED|REJECT|RESEQUENCE',
}
_COHORT_STATE_NAME = 'cohort_prepare_state.json'
# Bump this when review-template inclusion rules change; it invalidates cached Paper Trail output.
_COHORT_REVIEW_TEMPLATE_POLICY_VERSION = '10_consistent_read_selection_and_output_provenance'


def _read_assistant_cohort_manifest(path: str | Path) -> list[dict[str, str]]:
    """Read a typed cohort manifest; relative paths resolve against the manifest directory."""
    manifest_path = Path(path).expanduser().resolve()
    with open(manifest_path, newline='') as handle:
        sample = handle.read(4096)
        handle.seek(0)
        suffix = manifest_path.suffix.lower()
        if suffix in {'.tsv', '.txt'}:
            delimiter = '\t'
        elif suffix == '.csv':
            delimiter = ','
        else:
            try:
                delimiter = csv.Sniffer().sniff(sample, delimiters=',\t').delimiter
            except csv.Error:
                delimiter = '\t'
        raw_rows = list(csv.DictReader(handle, delimiter=delimiter))

    def resolve(value: str) -> str:
        if not value:
            return ''
        candidate = Path(value).expanduser()
        return str((candidate if candidate.is_absolute() else manifest_path.parent / candidate).resolve())

    rows = []
    seen_datasets = set()
    for line_number, raw in enumerate(raw_rows, start=2):
        row = {
            str(key or '').strip().lower().replace('-', '_').replace(' ', '_'):
                str(value or '').strip()
            for key, value in raw.items() if key is not None
        }
        if not any(row.values()):
            continue
        partner = row.get('partner_id') or row.get('partner') or row.get('partner_acronym')
        dataset = row.get('dataset') or row.get('batch_id')
        sample_map = row.get('sample_map') or row.get('ab1_map') or row.get('mailroom_map')
        marker_fasta = row.get('marker_fasta') or row.get('fasta')
        genome_dir = row.get('genome_dir') or row.get('genomic_16s_dir')
        maldi_file = row.get('maldi_file') or row.get('maldi')
        input_type = (row.get('input_type') or '').lower().replace('-', '_')
        input_type = _COHORT_INPUT_ALIASES.get(input_type, input_type)
        if not input_type:
            input_type = (
                'sanger' if sample_map else 'partner_fasta' if marker_fasta
                else 'genome_16s' if genome_dir else ('maldi_only' if maldi_file else '')
            )
        missing = [name for name, value in (('partner_id', partner), ('dataset', dataset)) if not value]
        if input_type and input_type not in _COHORT_INPUT_TYPES | {'maldi_only'}:
            raise ValueError(
                f'Cohort manifest row {line_number} has unsupported input_type {input_type!r}; '
                f'use one of {", ".join(sorted(_COHORT_INPUT_TYPES))}'
            )
        required = {
            'sanger': ('sample_map', sample_map), 'partner_fasta': ('marker_fasta', marker_fasta),
            'pacbio_hifi': ('marker_fasta', marker_fasta), 'genome_16s': ('genome_dir', genome_dir),
            'maldi_only': ('maldi_file', maldi_file),
        }
        if not input_type:
            missing.append('input_type (or a sample_map/marker_fasta/genome_dir/maldi_file)')
        else:
            name, value = required[input_type]
            if not value:
                missing.append(f'{name} (required for input_type {input_type})')
        if missing:
            raise ValueError(
                f'Cohort manifest row {line_number} is missing required field(s): {", ".join(missing)}'
            )
        if dataset.casefold() in seen_datasets:
            raise ValueError(f'Cohort manifest repeats dataset {dataset!r} on row {line_number}')
        seen_datasets.add(dataset.casefold())
        workbook = row.get('partner_workbook') or row.get('metadata_workbook')
        entry = {
            'partner_id': partner, 'dataset': dataset, 'input_type': input_type,
            'sample_map': resolve(sample_map) if input_type == 'sanger' else '',
            'read_dir': resolve(row.get('read_dir')) if input_type == 'sanger' else '',
            'marker_fasta': resolve(marker_fasta) if input_type in {'partner_fasta', 'pacbio_hifi'} else '',
            'marker_qc': resolve(row.get('marker_qc')) if input_type in {'partner_fasta', 'pacbio_hifi'} else '',
            'maldi_file': resolve(maldi_file),
            'genome_metadata': resolve(row.get('genome_metadata')) if input_type == 'genome_16s' else '',
            'genome_id_isolate_id': str(row.get('genome_id_isolate_id') or '').strip().lower() in {'1', 'true', 'yes'},
            'gtdb_summary': resolve(row.get('gtdb_summary')) if input_type == 'genome_16s' else '',
            'genome_dir': resolve(genome_dir) if input_type == 'genome_16s' else '',
            'partner_workbook': resolve(workbook),
            'partner_workbook_sheet': row.get('partner_workbook_sheet') or 'Bulk Upload Template',
        }
        for key in ('marker_fasta', 'marker_qc', 'maldi_file', 'genome_metadata', 'gtdb_summary'):
            if entry[key] and not Path(entry[key]).is_file():
                raise ValueError(f'Cohort manifest row {line_number}: {key} not found: {entry[key]}')
        if entry['genome_dir'] and not Path(entry['genome_dir']).is_dir():
            raise ValueError(f'Cohort manifest row {line_number}: genome_dir not found: {entry["genome_dir"]}')
        rows.append(entry)
    if not rows:
        raise ValueError(f'Cohort manifest contains no partner batches: {manifest_path}')
    return rows


def _merge_tsv_files(paths: list[str | Path], output: str | Path) -> list[dict]:
    rows = []
    fields = []
    for path in paths:
        with open(path, newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            for field in reader.fieldnames or []:
                if field not in fields:
                    fields.append(field)
            rows.extend(dict(row) for row in reader)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    return rows


def _cohort_write_rows(path: Path, fields: list[str], rows: list[dict]) -> Path:
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    return path


def _cohort_read_rows(path: Path) -> list[dict]:
    with open(path, newline='', encoding='utf-8-sig') as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter='\t')]


def _cohort_paper_trail_outputs_complete(trace_dir: Path) -> tuple[bool, str]:
    """Check that cached Paper Trail results include a complete review hand-off.

    An assembly report alone is insufficient: older/incomplete runs could contain
    all isolate rows and candidate FASTAs while marker_review_template.tsv was
    header-only.  PREPARE must rebuild those outputs rather than reuse them.
    """
    trace_dir = Path(trace_dir)
    required = {
        'assembly_report.tsv': trace_dir / 'assembly_report.tsv',
        'marker_review_template.tsv': trace_dir / 'marker_review_template.tsv',
        'failed_review_candidates.fasta': trace_dir / 'failed_qc_sequences' / 'failed_review_candidates.fasta',
        'failed_qc_manifest.tsv': trace_dir / 'failed_qc_sequences' / 'failed_qc_manifest.tsv',
        'failed_read_manifest.tsv': trace_dir / 'failed_qc_sequences' / 'failed_read_manifest.tsv',
    }
    missing_files = [name for name, path in required.items() if not path.is_file()]
    if missing_files:
        return False, 'missing required output(s): ' + ', '.join(missing_files)

    assembly_rows = _cohort_read_rows(required['assembly_report.tsv'])
    if not assembly_rows:
        return False, 'assembly_report.tsv contains no isolate records'
    review_rows = _cohort_read_rows(required['marker_review_template.tsv'])
    expected_review_ids = {
        _cohort_row_get(row, 'SequenceID', 'sequence_id')
        for row in assembly_rows
        if _cohort_row_get(row, 'QCClass', 'qc_class').upper() == 'FAIL_QC'
        or _cohort_row_get(row, 'Recommendation', 'recommendation').upper() == 'MANUAL_REVIEW'
    }
    observed_review_ids = {
        _cohort_row_get(row, 'sequence_id')
        for row in review_rows if _cohort_row_get(row, 'sequence_id')
    }
    missing_review_ids = sorted(expected_review_ids - observed_review_ids)
    if missing_review_ids:
        return False, (
            f'marker_review_template.tsv is incomplete ({len(missing_review_ids)} of '
            f'{len(expected_review_ids)} required review records missing; examples: '
            + ', '.join(missing_review_ids[:8]) + ')'
        )

    from branchmanager.utils.fasta import read_fasta
    candidate_records = {
        str(header).split()[0]: ''.join(str(sequence).split()).upper()
        for header, sequence in read_fasta(str(required['failed_review_candidates.fasta']))
    }
    missing_candidate_ids = []
    for row in review_rows:
        if _cohort_row_get(row, 'qc_class').upper() != 'FAIL_QC':
            continue
        eligible = _cohort_row_get(row, 'eligible_for_manual_acceptance').lower() in {'yes', 'true', '1'}
        if not eligible:
            continue
        sid = _cohort_row_get(row, 'sequence_id')
        candidate_sequence = candidate_records.get(sid, '') if sid else ''
        if not candidate_sequence:
            missing_candidate_ids.append(sid or '(blank sequence_id)')
            continue
        declared_length = _cohort_row_get(row, 'candidate_length')
        if declared_length:
            try:
                if int(float(declared_length)) != len(candidate_sequence):
                    missing_candidate_ids.append(f'{sid} (candidate length does not match FASTA)')
            except (TypeError, ValueError):
                missing_candidate_ids.append(f'{sid} (invalid candidate_length)')
    if missing_candidate_ids:
        return False, (
            'manual-approval eligibility is set but the failed-review candidate FASTA '
            'does not contain sequence(s): ' + ', '.join(missing_candidate_ids[:8])
        )
    return True, ''


def _cohort_row_get(row: dict, *names: str) -> str:
    lowered = {str(k).strip().lower().replace('_', ''): v for k, v in row.items() if k is not None}
    for name in names:
        value = lowered.get(name.lower().replace('_', ''))
        if value not in (None, ''):
            return str(value).strip()
    return ''


def _cohort_read_fallback_option_ids(row: dict) -> set[str]:
    raw = str(row.get('fallback_read_options') or '').strip()
    if not raw:
        return set()
    try:
        options = json.loads(raw)
    except (TypeError, ValueError, json.JSONDecodeError):
        return set()
    if not isinstance(options, list):
        return set()
    return {
        str(option.get('read_id') or '').strip()
        for option in options
        if isinstance(option, dict) and str(option.get('read_id') or '').strip()
    }


def _cohort_manual_read_selections(decisions_path: str | Path | None) -> dict[str, str]:
    """Read manual fallback read choices from a saved Studio review TSV."""
    if not decisions_path:
        return {}
    source = Path(decisions_path).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(f'Previous cohort review decisions file not found: {source}')
    selections: dict[str, str] = {}
    for row in _cohort_read_rows(source):
        sequence_id = str(row.get('sequence_id') or '').strip()
        selected_read_id = str(row.get('selected_read_id') or '').strip()
        source_type = str(row.get('read_selection_source') or '').strip().upper()
        if not sequence_id or not selected_read_id or source_type != 'MANUAL':
            continue
        previous = selections.get(sequence_id)
        if previous and previous != selected_read_id:
            raise ValueError(
                f'Conflicting manual fallback reads for {sequence_id!r} in {source}: '
                f'{previous!r} and {selected_read_id!r}'
            )
        selections[sequence_id] = selected_read_id
    return selections


def _cohort_restore_review_state(
    review_rows: list[dict],
    template_path: Path,
    decisions_path: str | Path | None = None,
) -> int:
    """Carry reviewer decisions and still-valid manual read choices forward.

    New QC evidence always wins. A manual selected_read_id is restored only if
    that read remains among the newly generated fallback candidates for the
    isolate; stale choices are discarded in favour of the current recommendation.
    """
    previous = {
        row.get('sequence_id', ''): row
        for row in _cohort_read_rows(template_path)
        if row.get('sequence_id')
    } if template_path.is_file() else {}

    saved = {}
    if decisions_path:
        source = Path(decisions_path).expanduser().resolve()
        if not source.is_file():
            raise FileNotFoundError(f'Previous cohort review decisions file not found: {source}')
        saved = {
            row.get('sequence_id', ''): row
            for row in _cohort_read_rows(source)
            if row.get('sequence_id')
        }

    restored = 0
    for row in review_rows:
        sequence_id = str(row.get('sequence_id') or '')
        source_row = saved.get(sequence_id) if sequence_id in saved else previous.get(sequence_id)
        if not source_row:
            continue
        has_decision = bool(str(source_row.get('decision') or '').strip())
        for key in ('decision', 'reviewer', 'reviewed_at', 'notes'):
            if key in source_row:
                row[key] = source_row.get(key, '')

        current_options = _cohort_read_fallback_option_ids(row)
        previous_selected = str(source_row.get('selected_read_id') or '').strip()
        previous_source = str(source_row.get('read_selection_source') or '').strip().upper()
        if previous_source == 'MANUAL' and previous_selected in current_options:
            row['selected_read_id'] = previous_selected
            row['read_selection_source'] = 'MANUAL'
        else:
            row['selected_read_id'] = str(row.get('recommended_read_id') or '') if current_options else ''
            row['read_selection_source'] = 'AUTO_LONGEST_USABLE_READ' if current_options else ''
        restored += int(has_decision)
    return restored


def _cohort_failed_qc_candidate(sequence: str) -> tuple[bool, float]:
    """Return whether a failed marker has a non-empty candidate for human review.

    Automated length and N-content thresholds determine the technical QC class,
    but do not prevent a human reviewer from considering a failed marker.
    Approval remains explicit and is validated against this prepared candidate
    sequence at cohort-resume time.
    """
    sequence = ''.join(str(sequence or '').split()).upper()
    if not sequence:
        return False, 100.0
    n_percent = 100.0 * sequence.count('N') / len(sequence)
    return True, n_percent


def _cohort_batch_fingerprint(args, batch: dict, selected_read_ids: dict[str, str] | None = None) -> str:
    import hashlib
    digest = hashlib.sha256()
    digest.update(_COHORT_REVIEW_TEMPLATE_POLICY_VERSION.encode())
    digest.update(sha256_file(batch['sample_map']).encode())
    digest.update(json.dumps(selected_read_ids or {}, sort_keys=True, separators=(',', ':')).encode())
    if batch['read_dir'] and Path(batch['read_dir']).is_dir():
        base = Path(batch['read_dir'])
        for path in sorted(p for p in base.rglob('*') if p.is_file()):
            info = path.stat()
            digest.update(f'{path.relative_to(base)}|{info.st_size}|{info.st_mtime_ns}'.encode())
    skip = {'out', 'db', 'cohort_manifest', 'cohort_resume', 'cohort_prepare_only', 'marker_review', 'config', 'command'}
    digest.update(json.dumps({k: v for k, v in sorted(vars(args).items()) if k not in skip}, default=str, sort_keys=True).encode())
    return digest.hexdigest()


def _cohort_register_and_import(project_db, batches, batch_ids, report_dir: Path) -> dict:
    """Register isolates and replay MALDI/workbook/genome sidecars on a project DB copy."""
    from branchmanager.genome_16s_import import import_genome_16s
    from branchmanager.project_state import (
        import_maldi_identifications, import_partner_workbook,
        write_maldi_import_report, write_partner_workbook_import_report,
    )
    reports: dict[str, list[str]] = {'maldi': [], 'genome': [], 'workbook': []}
    for index, batch in enumerate(batches, start=1):
        if batch['partner_workbook']:
            imported = import_partner_workbook(
                project_db, batch['partner_workbook'], sheet_name=batch['partner_workbook_sheet'],
                dataset=batch['dataset'],
            )
            reports['workbook'].append(str(write_partner_workbook_import_report(
                report_dir / f'{index:02d}_{batch["dataset"]}_workbook_import.tsv', imported)))
    by_dataset: dict[str, list[tuple[str, str]]] = {}
    for sequence_id, dataset in batch_ids.items():
        existing = project_db.resolve_sequence_id(sequence_id)
        if existing:
            with project_db.connect() as conn:
                stored = conn.execute('SELECT sequence, dataset FROM sequences WHERE id=?', (existing,)).fetchone()
            if not stored:
                continue
            if str(stored[1] or '') in project_db.get_dataset_names_by_role('baseline'):
                raise ValueError(f'Partner isolate {sequence_id!r} resolves to a registered baseline')
            if not str(stored[0] or '').strip():
                with project_db.connect() as conn:
                    conn.execute('UPDATE sequences SET dataset=? WHERE id=?', (dataset, existing))
            continue
        by_dataset.setdefault(dataset, []).append((sequence_id, ''))
    for dataset, records in by_dataset.items():
        project_db.upsert_dataset_role(dataset, 'candidate', genomes_available=False)
        project_db.insert_sequences(records, dataset=dataset)
    for index, batch in enumerate(batches, start=1):
        if batch['maldi_file']:
            results = import_maldi_identifications(project_db, batch['maldi_file'], default_dataset=batch['dataset'])
            reports['maldi'].append(str(write_maldi_import_report(
                report_dir / f'{index:02d}_{batch["dataset"]}_maldi_import.tsv', results)))
    for index, batch in enumerate(batches, start=1):
        if batch['input_type'] == 'genome_16s':
            gdir = report_dir / f'{index:02d}_{batch["dataset"]}_genome_16s'
            import_genome_16s(
                project_db, batch['genome_dir'], gdir, genome_metadata=batch['genome_metadata'] or None,
                dataset=batch['dataset'], partner_id=batch['partner_id'],
                genome_id_isolate_id=batch['genome_id_isolate_id'],
            )
            reports['genome'].append(str(gdir / 'genome_16s_import_report.tsv'))
    return reports


def _cohort_copy_db(original_db: Path, target: Path) -> None:
    import sqlite3
    target.unlink(missing_ok=True)
    with sqlite3.connect(str(original_db)) as source, sqlite3.connect(str(target)) as dest:
        source.backup(dest)


def _cohort_prepare(args, root: Path, batches: list[dict]) -> Path:
    """Stage 1: onboarding, QC, review template and evidence manifest. No Hiring Panel, no publish."""
    import fcntl
    from branchmanager.onboarding import validate_submission, write_onboarding_outputs
    from branchmanager.utils.fasta import read_fasta, write_fasta

    prep_root = root / '00_cohort_preparation'
    prep_root.mkdir(parents=True, exist_ok=True)
    manual_read_selections = _cohort_manual_read_selections(getattr(args, 'marker_review', None))
    if manual_read_selections:
        print(
            f'[assistant-cohort] Applying {len(manual_read_selections)} reviewer-selected fallback read(s) '
            'where a valid consensus cannot be formed.'
        )
    manifest = RunManifest(prep_root, 'assistant_cohort_prepare')
    manifest.add_input(args.cohort_manifest, role='cohort_manifest')
    manifest.add_input(args.partner_metadata, role='cumulative_partner_catalogue')
    for batch in batches:
        for key in ('sample_map', 'marker_fasta', 'marker_qc', 'maldi_file', 'genome_metadata', 'gtdb_summary', 'partner_workbook'):
            if batch[key]:
                manifest.add_input(batch[key], role=f'{key}:{batch["dataset"]}')
        if batch['read_dir']:
            manifest.add_input(batch['read_dir'], role=f'raw_reads:{batch["dataset"]}')
        if batch['genome_dir']:
            manifest.add_input(batch['genome_dir'], role=f'genomic_16s_output:{batch["dataset"]}')
        if batch['gtdb_summary']:
            manifest.set_reference_identity(**{f'gtdb_summary_sha256:{batch["dataset"]}': sha256_file(batch['gtdb_summary'])})
    try:
        batch_ids: dict[str, str] = {}
        results = []
        for index, batch in enumerate(batches, start=1):
            res = {'batch': batch, 'qc_rows': [], 'review_rows': [], 'passing_fasta': '', 'failed_fasta': '',
                   'failed_manifest': ''}
            source = _COHORT_MARKER_SOURCE.get(batch['input_type'], '')
            if batch['input_type'] == 'sanger':
                onboard_dir = prep_root / f'{index:02d}_{batch["dataset"]}_onboarding'
                result = validate_submission(
                    batch['sample_map'], partner_metadata=args.partner_metadata,
                    read_dir=batch['read_dir'] or None, primers=args.primers or (),
                    expected_partner_id=batch['partner_id'], dataset=batch['dataset'],
                )
                outputs = write_onboarding_outputs(onboard_dir, result)
                manifest.add_output(outputs['report'], role=f"onboarding_report:{batch['dataset']}")
                if result['status'] != 'PASS':
                    raise RuntimeError(
                        f"Onboarding failed for {batch['partner_id']} / {batch['dataset']}; inspect {outputs['report']}"
                    )
                for row in result['normalised']:
                    sequence_id = str(row['sequence_id'])
                    if sequence_id in batch_ids:
                        raise ValueError(
                            f"Isolate {sequence_id!r} appears in multiple cohort batches "
                            f"({batch_ids[sequence_id]!r} and {batch['dataset']!r})"
                        )
                    batch_ids[sequence_id] = batch['dataset']
                batch_sequence_ids = {str(row['sequence_id']) for row in result['normalised']}
                batch_manual_read_selections = {
                    sid: read_id for sid, read_id in manual_read_selections.items()
                    if sid in batch_sequence_ids
                }
                trace_dir = prep_root / f'{index:02d}_{batch["dataset"]}_paper_trail'
                fingerprint = _cohort_batch_fingerprint(args, batch, batch_manual_read_selections)
                state_file = trace_dir / '.prepared.json'
                outputs_complete, outputs_problem = _cohort_paper_trail_outputs_complete(trace_dir)
                try:
                    fingerprint_matches = (
                        state_file.is_file()
                        and json.loads(state_file.read_text()).get('fingerprint') == fingerprint
                    )
                except (OSError, ValueError, json.JSONDecodeError):
                    fingerprint_matches = False
                reusable = (
                    fingerprint_matches
                    and (trace_dir / 'assembled.fasta').is_file()
                    and outputs_complete
                )
                if reusable:
                    print(f'[assistant-cohort] Reusing unchanged, complete Paper Trail output for {batch["dataset"]}')
                else:
                    if state_file.is_file() and not outputs_complete:
                        print(
                            f'[assistant-cohort] Rebuilding incomplete Paper Trail output for {batch["dataset"]}: '
                            f'{outputs_problem}'
                        )
                    trace_args = argparse.Namespace(**vars(args))
                    trace_args.command = 'paper-trail'
                    trace_args.input = []
                    trace_args.sample_map = outputs['read_map']
                    trace_args.out = str(trace_dir)
                    trace_args.read_dir = batch['read_dir'] or None
                    trace_args.screen_ref = args.ref
                    trace_args.screen_taxa = args.taxa
                    trace_args.mailroom_summary = None
                    trace_args.mailroom_status = 'PASS'
                    trace_args.assemble = True
                    trace_args.recursive = True
                    trace_args.selected_read_ids = batch_manual_read_selections
                    cmd_paper_trail(trace_args)
                    outputs_complete, outputs_problem = _cohort_paper_trail_outputs_complete(trace_dir)
                    if not outputs_complete:
                        raise RuntimeError(
                            f'Paper Trail did not produce a complete review hand-off for {batch["dataset"]}: '
                            f'{outputs_problem}. The cohort review template has not been generated.'
                        )
                    state_file.write_text(json.dumps({'fingerprint': fingerprint}))
                res['passing_fasta'] = str(trace_dir / 'assembled.fasta')
                res['failed_fasta'] = str(trace_dir / 'failed_qc_sequences' / 'failed_review_candidates.fasta')
                res['qc_rows'] = _cohort_read_rows(trace_dir / 'assembly_report.tsv')
                template = trace_dir / 'marker_review_template.tsv'
                res['review_rows'] = _cohort_read_rows(template) if template.is_file() else []
                failed_manifest = trace_dir / 'failed_qc_sequences' / 'failed_qc_manifest.tsv'
                res['failed_manifest'] = str(failed_manifest) if failed_manifest.is_file() else ''
            elif source:
                bdir = prep_root / f'{index:02d}_{batch["dataset"]}_marker_input'
                bdir.mkdir(parents=True, exist_ok=True)
                records = [(str(h).split()[0], str(s).upper()) for h, s in read_fasta(batch['marker_fasta'])]
                seen = set()
                for sequence_id, sequence in records:
                    if sequence_id in seen or sequence_id in batch_ids:
                        raise ValueError(f'Duplicate marker ID {sequence_id!r} in cohort batch {batch["dataset"]!r}')
                    if not sequence:
                        raise ValueError(f'Marker {sequence_id!r} in {batch["dataset"]!r} has an empty sequence')
                    seen.add(sequence_id)
                    batch_ids[sequence_id] = batch['dataset']
                qc_by_id = {}
                if batch['marker_qc']:
                    for row in _cohort_read_rows(Path(batch['marker_qc'])):
                        sid = _cohort_row_get(row, 'SequenceID', 'sequence_id', 'id')
                        if sid:
                            qc_by_id[sid] = row
                passing, failed = [], []
                for sequence_id, sequence in records:
                    row = qc_by_id.get(sequence_id)
                    if row is None:
                        qc_class, recommendation, reasons = 'QUALITY_UNVERIFIED', 'MANUAL_REVIEW', 'no_marker_qc_sidecar'
                        row = {}
                    else:
                        qc_class = (_cohort_row_get(row, 'QCClass', 'qc_class') or 'QUALITY_UNVERIFIED').upper()
                        recommendation = _cohort_row_get(row, 'Recommendation') or 'MANUAL_REVIEW'
                        reasons = _cohort_row_get(row, 'Reasons')
                    qc_row = dict(row)
                    qc_row.update({'SequenceID': sequence_id, 'QCClass': qc_class,
                                   'Recommendation': recommendation, 'Reasons': reasons,
                                   'OutputLength': len(sequence)})
                    res['qc_rows'].append(qc_row)
                    (failed if qc_class == 'FAIL_QC' else passing).append((sequence_id, sequence))
                    if qc_class == 'PASS_HIGH_CONFIDENCE' and recommendation != 'MANUAL_REVIEW':
                        continue

                    min_candidate_length = int(getattr(args, 'min_length', 800))
                    failed_candidate, n_percent = _cohort_failed_qc_candidate(sequence)
                    review_type = (
                        'FAILED_QC_REVIEW' if qc_class == 'FAIL_QC'
                        else 'UNVERIFIED_QC_REVIEW' if qc_class == 'QUALITY_UNVERIFIED' else 'WARNING_REVIEW'
                    )
                    # Every non-passing marker is listed for a common manual review.
                    # Failed markers with any non-empty candidate sequence can be
                    # explicitly approved by the reviewer, regardless of the QC
                    # threshold(s) failed. No candidate means no approval option.
                    eligible = failed_candidate if qc_class == 'FAIL_QC' else bool(sequence)
                    manual_acceptance_basis = (
                        'HUMAN_REVIEW_OVERRIDE'
                        if qc_class == 'FAIL_QC' and eligible else ''
                    )
                    manual_acceptance_reason = (
                        f'Candidate sequence is available ({len(sequence)} bp; {n_percent:.3f}% Ns). '
                        f'Automated QC class remains FAIL_QC. Recorded reason(s): {reasons}. '
                        'The reviewer may explicitly approve this marker regardless of which automated QC '
                        'threshold(s) failed. Inspect all read, chromatogram, sequence and taxonomy evidence; '
                        'a named reviewer and explanatory notes are mandatory. This is a human override, '
                        'not an automated pass.'
                        if manual_acceptance_basis else ''
                    )
                    res['review_rows'].append({
                        'sequence_id': sequence_id, 'review_type': review_type,
                        'eligible_for_manual_acceptance': 'yes' if eligible else 'no', 'qc_class': qc_class,
                        'qc_reasons': reasons, 'output_length': len(sequence), 'n_percent': f'{n_percent:.3f}',
                        'quality_status': 'PARTNER_SUPPLIED_QC' if row else 'NO_QC_PROVENANCE',
                        'manual_acceptance_basis': manual_acceptance_basis,
                        'manual_acceptance_reason': manual_acceptance_reason,
                        'minimum_acceptable_length': '',
                        'length_deficit_bp': max(0, min_candidate_length - len(sequence)),
                        'candidate_fasta': str(bdir / ('failed_markers.fasta' if qc_class == 'FAIL_QC' else 'passing_markers.fasta')),
                    })
                write_fasta(passing, str(bdir / 'passing_markers.fasta'))
                write_fasta(failed, str(bdir / 'failed_markers.fasta'))
                res['passing_fasta'] = str(bdir / 'passing_markers.fasta')
                res['failed_fasta'] = str(bdir / 'failed_markers.fasta')
            for qc_row in res['qc_rows']:
                qc_row['marker_source_type'] = source
                qc_row['partner_id'] = batch['partner_id']
                qc_row['cohort_dataset'] = batch['dataset']
            for review_row in res['review_rows']:
                review_row['marker_source_type'] = source
                review_row['partner_id'] = batch['partner_id']
                review_row['dataset'] = batch['dataset']
                allowed = _COHORT_ALLOWED.get(review_row.get('review_type', ''), '')
                if str(review_row.get('eligible_for_manual_acceptance', '')).lower() != 'yes':
                    allowed = 'REJECT|RESEQUENCE'
                review_row['allowed_decisions'] = allowed
            results.append(res)

        combined_records, used_ids, failed_records, source_by_id = [], set(), {}, {}
        for res in results:
            for fasta_key, bucket in (('passing_fasta', None), ('failed_fasta', failed_records)):
                path = res[fasta_key]
                if not path or not Path(path).is_file():
                    continue
                for sequence_id, sequence in read_fasta(path):
                    sequence_id = str(sequence_id).split()[0]
                    if bucket is None:
                        if sequence_id in used_ids:
                            raise ValueError(f'Cohort contains duplicate passing marker ID {sequence_id!r}')
                        used_ids.add(sequence_id)
                        combined_records.append((sequence_id, sequence))
                    else:
                        bucket[sequence_id] = sequence
        combined_fasta = prep_root / 'cohort_passing_markers.fasta'
        write_fasta(combined_records, str(combined_fasta))
        qc_rows = [row for res in results for row in res['qc_rows']]
        qc_fields = []
        for row in qc_rows:
            qc_fields.extend(k for k in row if k not in qc_fields and not str(k).startswith('_'))
        combined_qc = _cohort_write_rows(prep_root / 'cohort_marker_qc.tsv', qc_fields or ['SequenceID', 'QCClass'], qc_rows)

        review_rows = [row for res in results for row in res['review_rows']]
        for row in review_rows:
            allowed = _COHORT_ALLOWED.get(row.get('review_type', ''), '')
            if str(row.get('eligible_for_manual_acceptance', '')).lower() != 'yes':
                allowed = 'REJECT|RESEQUENCE'
            row['allowed_decisions'] = allowed
        template_path = prep_root / 'cohort_marker_review_template.tsv'
        restored_decisions = _cohort_restore_review_state(
            review_rows, template_path, getattr(args, 'marker_review', None),
        )
        restored_read_selections = sum(
            str(row.get('read_selection_source') or '').upper() == 'MANUAL'
            and bool(str(row.get('selected_read_id') or '').strip())
            for row in review_rows
        )
        if getattr(args, 'marker_review', None):
            manifest.add_input(args.marker_review, role='previous_cohort_review_decisions')
        _cohort_write_rows(template_path, _COHORT_REVIEW_FIELDS, review_rows)
        if restored_decisions:
            print(f'[assistant-cohort] Restored {restored_decisions} existing manual review decision(s) into the regenerated template.')
        if restored_read_selections:
            print(
                f'[assistant-cohort] Preserved {restored_read_selections} manual fallback read selection(s) '
                'in the combined cohort review template.'
            )

        recoverable = [(sid, failed_records[sid]) for sid in sorted(failed_records)
                       if any(r['sequence_id'] == sid
                              and r.get('review_type') == 'FAILED_QC_REVIEW'
                              and str(r.get('eligible_for_manual_acceptance', '')).lower() == 'yes'
                              for r in review_rows)]
        recoverable_fasta = prep_root / 'cohort_recoverable_failed_markers.fasta'
        write_fasta(recoverable, str(recoverable_fasta))
        failed_manifests = [res['failed_manifest'] for res in results if res['failed_manifest']]
        failed_rows = _merge_tsv_files(failed_manifests, prep_root / 'cohort_failed_qc_manifest.tsv') if failed_manifests else []
        if not failed_manifests:
            (prep_root / 'cohort_failed_qc_manifest.tsv').write_text('SequenceID\tQCClass\tRecommendation\tStatus\tReasons\n')
        failure_reasons = {}
        for row in failed_rows:
            sid = str(row.get('SequenceID') or '').strip()
            if sid:
                failure_reasons[sid] = str(row.get('Reasons') or row.get('Status') or 'failed_qc')
        for row in qc_rows:
            if str(row.get('QCClass') or '').upper() == 'FAIL_QC' and row.get('SequenceID'):
                failure_reasons.setdefault(row['SequenceID'], str(row.get('Reasons') or 'failed_qc'))

        review_ids = {r['sequence_id'] for r in review_rows}
        recoverable_ids = {sid for sid, _ in recoverable}
        evidence_rows = []
        for res in results:
            for row in res['qc_rows']:
                sid, qc_class = row['SequenceID'], str(row.get('QCClass') or '').upper()
                if qc_class == 'PASS_HIGH_CONFIDENCE':
                    status = 'INCLUDED_PASSING'
                elif qc_class == 'FAIL_QC':
                    if sid in recoverable_ids:
                        status = 'FAILED_RECOVERABLE_PENDING_REVIEW'
                    elif sid in review_ids:
                        status = 'FAILED_PENDING_REVIEW'
                    else:
                        status = 'FAILED_UNRECOVERABLE'
                elif sid in review_ids:
                    status = 'PENDING_REVIEW'
                else:
                    status = 'INCLUDED_UNREVIEWED'
                evidence_rows.append({
                    'sequence_id': sid, 'partner_id': res['batch']['partner_id'], 'dataset': res['batch']['dataset'],
                    'source_type': row['marker_source_type'], 'qc_class': qc_class, 'status': status,
                    'in_vsearch_fasta': 'yes' if sid in used_ids else 'no',
                })
        for res in results:
            if res['batch']['input_type'] == 'genome_16s':
                evidence_rows.append({
                    'sequence_id': '*', 'partner_id': res['batch']['partner_id'], 'dataset': res['batch']['dataset'],
                    'source_type': 'GENOME_BARRNAP', 'qc_class': 'GENOME_DERIVED', 'status': 'SUPPLEMENTARY_SECONDARY_POOL',
                    'in_vsearch_fasta': 'no',
                })
        evidence_manifest = _cohort_write_rows(
            prep_root / 'cohort_marker_evidence_manifest.tsv',
            ['sequence_id', 'partner_id', 'dataset', 'source_type', 'qc_class', 'status', 'in_vsearch_fasta'], evidence_rows)
        failure_path = prep_root / 'cohort_qc_failures.tsv'
        _cohort_write_rows(failure_path, ['ID', 'Reasons'],
                           [{'ID': k, 'Reasons': v} for k, v in sorted(failure_reasons.items())])

        # Import sidecars on a scratch copy only to produce reports; the live DB is never touched here.
        original_db = Path(args.db).expanduser().resolve()
        sidecar_reports: dict[str, list[str]] = {'maldi': [], 'genome': [], 'workbook': []}
        if original_db.is_file() and any(b['maldi_file'] or b['genome_dir'] or b['partner_workbook'] for b in batches):
            scratch = prep_root / '.prepare_scratch.sqlite'
            with open(original_db.with_name(original_db.name + '.lock'), 'a+') as lock:
                fcntl.flock(lock.fileno(), fcntl.LOCK_SH)
                try:
                    _cohort_copy_db(original_db, scratch)
                finally:
                    fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
            try:
                scratch_db = Database(str(scratch))
                scratch_db.initialise()
                sidecar_reports = _cohort_register_and_import(scratch_db, batches, batch_ids, prep_root)
            finally:
                scratch.unlink(missing_ok=True)
        for kind, name in (('maldi', 'cohort_maldi_import_report.tsv'), ('genome', 'cohort_genome_import_report.tsv')):
            if sidecar_reports[kind]:
                _merge_tsv_files(sidecar_reports[kind], prep_root / name)
            else:
                (prep_root / name).write_text('result\tdetail\nNOT_APPLICABLE\tno input of this kind in the cohort manifest\n')

        outputs = {
            'passing_markers': combined_fasta, 'marker_qc': combined_qc, 'review_template': template_path,
            'recoverable_failed': recoverable_fasta, 'evidence_manifest': evidence_manifest,
            'qc_failures': failure_path, 'maldi_report': prep_root / 'cohort_maldi_import_report.tsv',
            'genome_report': prep_root / 'cohort_genome_import_report.tsv',
        }
        for role, path in outputs.items():
            manifest.add_output(path, role=f'cohort_{role}')
        checksums = {name: sha256_file(path) for name, path in outputs.items() if name != 'review_template'}
        (prep_root / _COHORT_STATE_NAME).write_text(json.dumps({
            'batches': batches,
            'batch_ids': batch_ids,
            'checksums': {k: [str(outputs[k]), v] for k, v in checksums.items()},
            'marker_source_by_id': {r['SequenceID']: r['marker_source_type'] for r in qc_rows},
            'review_template_policy_version': _COHORT_REVIEW_TEMPLATE_POLICY_VERSION,
            'failed_qc_review_policy': {
                'manual_approval': 'REVIEWER_OVERRIDE_FOR_ANY_NONEMPTY_CANDIDATE',
                'minimum_length': int(getattr(args, 'min_length', 800)),
                'maximum_n_percent': float(getattr(args, 'max_n_percent', 5.0)),
                'minimum_mean_quality': float(getattr(args, 'min_mean_quality', 20.0)),
                'maximum_output_expected_errors': float(getattr(args, 'max_output_expected_errors', 5.0)),
                'maximum_conflict_density': float(getattr(args, 'max_conflict_density', 1.0)),
                'note': 'These are automated QC thresholds, not manual-approval eligibility gates. Any non-empty prepared candidate can be approved with reviewer identity and notes.',
            },
        }, indent=2, sort_keys=True))
        status = 'COMPLETE_WITH_WARNINGS' if review_rows or failure_reasons else 'COMPLETE'
        manifest.add_stage('cohort_prepare', status,
                           detail=f'{len(batches)} batch(es); {len(used_ids)} marker(s) for assessment; '
                                  f'{len(review_rows)} awaiting review; {len(failure_reasons)} QC failure(s)')
        manifest.finish(status)
        print(
            f'[assistant-cohort] PREPARE {status}: {len(batches)} batch(es), {len(batch_ids)} isolate(s), '
            f'{len(review_rows)} manual-review row(s).\n'
            f'  Review template: {template_path}\n'
            f'  Next: fill decisions, then run `branchmanager assistant --cohort-resume {root} --marker-review <decisions.tsv> ...`'
        )
    except BaseException as exc:
        manifest.finish('FAILED', error=exc)
        raise
    return prep_root




def _cohort_validate_decisions(prep_root: Path, decisions_path: Path | None):
    """Validate manual-review decisions against the prepared cohort evidence."""
    from branchmanager.marker_evidence import REVIEW_DECISIONS, validate_review_decision
    from branchmanager.utils.fasta import read_fasta

    template = {
        r['sequence_id']: r
        for r in _cohort_read_rows(prep_root / 'cohort_marker_review_template.tsv')
        if r.get('sequence_id')
    }
    if not decisions_path:
        return [], []

    recoverable_path = prep_root / 'cohort_recoverable_failed_markers.fasta'
    recoverable = {
        str(header).split()[0]: ''.join(str(sequence).split()).upper()
        for header, sequence in read_fasta(str(recoverable_path))
    } if recoverable_path.is_file() else {}

    qc_by_id = {
        _cohort_row_get(row, 'SequenceID', 'sequence_id'): row
        for row in _cohort_read_rows(prep_root / 'cohort_marker_qc.tsv')
        if _cohort_row_get(row, 'SequenceID', 'sequence_id')
    }

    # Manual approval is deliberately reviewer-controlled. The parent resume
    # routine has already verified the preparation checksums; this validator
    # verifies that a candidate exists, but does not impose a second automated
    # QC pass gate.

    errors, accepted, seen = [], [], set()
    for line, row in enumerate(_cohort_read_rows(decisions_path), start=2):
        sid = _cohort_row_get(row, 'sequence_id')
        decision = _cohort_row_get(row, 'decision', 'review_decision').upper()
        if not sid or not decision:
            continue
        if sid in seen:
            errors.append(f'line {line}: duplicate decision for {sid}')
            continue
        seen.add(sid)

        if sid not in template:
            errors.append(f'line {line}: {sid} is not in the cohort review template')
            continue
        if decision not in REVIEW_DECISIONS:
            errors.append(f'line {line}: {sid} has unsupported decision {decision!r}')
            continue

        info = template[sid]
        review_type = str(info.get('review_type') or '').upper()
        qc_class = str(info.get('qc_class') or '').upper()
        qc_row = qc_by_id.get(sid, {})
        prepared_qc_class = _cohort_row_get(qc_row, 'QCClass', 'qc_class').upper()
        if prepared_qc_class and qc_class != prepared_qc_class:
            errors.append(
                f'line {line}: {sid} QC class was changed in the review template '
                f'({qc_class!r}; prepared value is {prepared_qc_class!r})'
            )
            continue

        expected_review_type = {
            'APPROVE': {'WARNING_REVIEW'},
            'APPROVE_FAIL_QC': {'FAILED_QC_REVIEW'},
            'APPROVE_UNVERIFIED': {'UNVERIFIED_QC_REVIEW'},
        }
        if decision in expected_review_type and review_type not in expected_review_type[decision]:
            errors.append(
                f'line {line}: {decision} is not valid for review_type {review_type or "(blank)"!r}'
            )
            continue

        expected_type_for_qc = {
            'FAIL_QC': 'FAILED_QC_REVIEW',
            'QUALITY_UNVERIFIED': 'UNVERIFIED_QC_REVIEW',
            'PASS_WITH_WARNINGS': 'WARNING_REVIEW',
        }
        expected_type = expected_type_for_qc.get(qc_class)
        if expected_type and review_type != expected_type:
            errors.append(
                f'line {line}: {sid} review_type {review_type!r} does not match prepared QC class {qc_class!r}'
            )
            continue

        has_candidate = str(info.get('eligible_for_manual_acceptance', '')).lower() == 'yes'
        if review_type == 'FAILED_QC_REVIEW':
            if qc_class != 'FAIL_QC':
                errors.append(f'line {line}: {sid} is FAILED_QC_REVIEW but its QC class is {qc_class!r}')
                continue
            # Failed records can always be rejected or marked for resequencing.
            # Explicit approval is permitted for any non-empty candidate in the
            # checksummed recoverable FASTA, regardless of the failed QC metric(s).
            if decision == 'APPROVE_FAIL_QC':
                sequence = recoverable.get(sid, '')
                if not sequence:
                    errors.append(
                        f'line {line}: {sid} has no non-empty candidate in '
                        'cohort_recoverable_failed_markers.fasta; it cannot be approved'
                    )
                    continue
                has_candidate = True

        error = validate_review_decision(
            decision, qc_class, reviewer=_cohort_row_get(row, 'reviewer'),
            notes=_cohort_row_get(row, 'notes', 'review_notes'),
            has_consensus=has_candidate,
        )
        if error:
            errors.append(f'line {line}: {sid}: {error}')
            continue

        accepted.append({
            'sequence_id': sid,
            'decision': decision,
            'reviewer': _cohort_row_get(row, 'reviewer'),
            'notes': _cohort_row_get(row, 'notes', 'review_notes'),
            'reviewed_at': _cohort_row_get(row, 'reviewed_at'),
            'qc_class': qc_class,
            'review_type': review_type,
            'candidate_fasta': info.get('candidate_fasta', ''),
        })
    return accepted, errors



def _cohort_review_application_audit_rows(
    template_rows: list[dict],
    decisions: list[dict],
    passing_sequences: dict[str, str],
    failed_sequences: dict[str, str],
    included_ids: set[str],
    decisions_source: str = '',
) -> list[dict]:
    """Build a per-isolate audit of how validated human decisions were applied."""
    decided = {str(row.get('sequence_id') or ''): row for row in decisions}
    output = []
    for info in template_rows:
        sid = str(info.get('sequence_id') or '').strip()
        if not sid:
            continue
        decision = decided.get(sid, {})
        value = str(decision.get('decision') or '').upper()
        sequence = passing_sequences.get(sid) or failed_sequences.get(sid) or ''
        included = sid in included_ids
        if value == 'APPROVE_FAIL_QC' and included:
            status = 'INCLUDED_BY_MANUAL_FAILED_QC_OVERRIDE'
            reason = 'Non-empty failed-QC candidate included following validated reviewer approval; technical QC remains FAIL_QC'
        elif value == 'APPROVE' and included:
            status = 'INCLUDED_AFTER_MANUAL_APPROVAL'
            reason = 'Warning-level marker included following validated reviewer approval'
        elif value == 'APPROVE_UNVERIFIED' and included:
            status = 'INCLUDED_AFTER_MANUAL_UNVERIFIED_APPROVAL'
            reason = 'Marker without an automated QC sidecar included following validated reviewer approval'
        elif value == 'RESEQUENCE':
            status = 'EXCLUDED_RESEQUENCE'
            reason = 'Reviewer requested resequencing; marker excluded from joint assessment'
        elif value == 'REJECT':
            status = 'EXCLUDED_REJECTED'
            reason = 'Reviewer rejected the marker; marker excluded from joint assessment'
        elif included:
            status = 'INCLUDED_WITHOUT_MANUAL_REVIEW'
            reason = 'Marker entered the joint set without an individual manual-review decision'
        else:
            status = 'EXCLUDED_PENDING_REVIEW'
            reason = 'No valid approval decision was supplied for this review-required marker'
        output.append({
            'sequence_id': sid,
            'partner_id': info.get('partner_id', ''),
            'dataset': info.get('dataset', ''),
            'review_type': info.get('review_type', ''),
            'prepared_qc_class': info.get('qc_class', ''),
            'qc_reasons': info.get('qc_reasons', ''),
            'decision': value,
            'reviewer': decision.get('reviewer', ''),
            'reviewed_at': decision.get('reviewed_at', ''),
            'notes': decision.get('notes', ''),
            'candidate_available': 'YES' if sequence else 'NO',
            'candidate_length': len(sequence),
            'candidate_sha256': hashlib.sha256(sequence.encode('ascii', errors='replace')).hexdigest() if sequence else '',
            'included_in_joint_assessment': 'YES' if included else 'NO',
            'application_status': status,
            'application_reason': reason,
            'decision_source_file': decisions_source,
        })
    return output


def _cohort_assess(args, root: Path) -> None:
    """Stage 2: validate decisions, build the joint marker set, assess once, publish after integrity checks."""
    import fcntl
    import hashlib
    import sqlite3
    from branchmanager.marker_evidence import record_review_decisions
    from branchmanager.partner_metadata import load_partner_sequencing_metadata
    from branchmanager.utils.fasta import read_fasta, write_fasta

    prep_root = root / '00_cohort_preparation'
    joint_root = root / '01_joint_assistant'
    state_path = prep_root / _COHORT_STATE_NAME
    if not state_path.is_file():
        raise SystemExit(f'[assistant-cohort] No completed preparation found in {prep_root}; run the prepare stage first.')
    state = json.loads(state_path.read_text())
    for name, (path, expected) in state['checksums'].items():
        if not Path(path).is_file() or sha256_file(path) != expected:
            raise SystemExit(
                f'[assistant-cohort] Preparation output {name} changed since it was prepared ({path}); '
                'rerun the prepare stage rather than resuming.'
            )
    batches, batch_ids = state['batches'], state['batch_ids']
    decisions_path = Path(args.marker_review).expanduser().resolve() if getattr(args, 'marker_review', None) else None
    decisions, errors = _cohort_validate_decisions(prep_root, decisions_path)
    if errors:
        raise SystemExit('[assistant-cohort] Invalid manual-review decision(s); nothing was changed:\n  ' + '\n  '.join(errors))

    manifest = RunManifest(root, 'assistant_cohort')
    manifest.add_input(prep_root / _COHORT_STATE_NAME, role='cohort_prepare_state')
    manifest.add_input(prep_root / 'cohort_marker_review_template.tsv', role='cohort_marker_review_template')
    if decisions_path:
        manifest.add_input(decisions_path, role='marker_review_decisions')
        history = root / 'review_history'
        history.mkdir(parents=True, exist_ok=True)
        digest = hashlib.sha256(decisions_path.read_bytes()).hexdigest()[:12]
        kept = history / f'{digest}_{decisions_path.name}'
        if not kept.exists():
            kept.write_bytes(decisions_path.read_bytes())
    approved_failed = {d['sequence_id'] for d in decisions if d['decision'] == 'APPROVE_FAIL_QC'}
    decided = {d['sequence_id']: d for d in decisions}

    # Warning-level and quality-unverified markers need an explicit approval before they are assessed.
    template_rows = {r['sequence_id']: r for r in _cohort_read_rows(prep_root / 'cohort_marker_review_template.tsv')
                     if r.get('sequence_id')}
    held_back = {}
    for sid, info in template_rows.items():
        if info.get('review_type') == 'FAILED_QC_REVIEW':
            continue
        decision = decided.get(sid, {}).get('decision', '')
        if decision not in {'APPROVE', 'APPROVE_UNVERIFIED'}:
            held_back[sid] = f'manual_review_{decision.lower() or "pending"}'
    records = [(str(h).split()[0], s) for h, s in read_fasta(str(prep_root / 'cohort_passing_markers.fasta'))
               if str(h).split()[0] not in held_back]
    failed_pool = {str(h).split()[0]: s for h, s in read_fasta(str(prep_root / 'cohort_recoverable_failed_markers.fasta'))}
    missing_approved_failed = sorted(approved_failed - set(failed_pool))
    if missing_approved_failed:
        raise SystemExit(
            '[assistant-cohort] Approved failed-QC marker(s) are missing from the prepared recoverable FASTA: '
            + ', '.join(missing_approved_failed)
            + '. No database changes were made; rerun PREPARE and review the regenerated template.'
        )
    present = {sid for sid, _ in records}
    for sid in sorted(approved_failed):
        if sid not in present:
            records.append((sid, failed_pool[sid]))
    joint_inputs = root / '01_joint_inputs'
    joint_inputs.mkdir(parents=True, exist_ok=True)
    joint_fasta = joint_inputs / 'cohort_joint_markers.fasta'
    write_fasta(records, str(joint_fasta))
    included = {sid for sid, _ in records}

    failure_rows = _cohort_read_rows(prep_root / 'cohort_qc_failures.tsv')
    failure_reasons = {r['ID']: r.get('Reasons', '') for r in failure_rows if r.get('ID') and r['ID'] not in included}
    failure_reasons.update({sid: reason for sid, reason in held_back.items() if sid not in included})
    failure_path = joint_inputs / 'cohort_qc_failures.tsv'
    _cohort_write_rows(failure_path, ['ID', 'Reasons'], [{'ID': k, 'Reasons': v} for k, v in sorted(failure_reasons.items())])
    qc_rows = _cohort_read_rows(prep_root / 'cohort_marker_qc.tsv')
    qc_fields = list(qc_rows[0].keys()) if qc_rows else ['SequenceID', 'QCClass']
    joint_qc = _cohort_write_rows(joint_inputs / 'cohort_joint_marker_qc.tsv', qc_fields,
                                  [r for r in qc_rows if r.get('SequenceID') in included])

    passing_pool = {
        str(h).split()[0]: str(sequence).upper()
        for h, sequence in read_fasta(str(prep_root / 'cohort_passing_markers.fasta'))
    }
    review_template_rows = _cohort_read_rows(prep_root / 'cohort_marker_review_template.tsv')
    review_audit_rows = _cohort_review_application_audit_rows(
        review_template_rows, decisions, passing_pool, failed_pool, included,
        str(decisions_path or ''),
    )
    review_audit_path = joint_inputs / 'cohort_review_application_audit.tsv'
    audit_fields = [
        'sequence_id', 'partner_id', 'dataset', 'review_type', 'prepared_qc_class',
        'qc_reasons', 'decision', 'reviewer', 'reviewed_at', 'notes',
        'candidate_available', 'candidate_length', 'candidate_sha256',
        'included_in_joint_assessment', 'application_status', 'application_reason',
        'decision_source_file',
    ]
    _cohort_write_rows(review_audit_path, audit_fields, review_audit_rows)

    original_db = Path(args.db).expanduser().resolve()
    if not original_db.is_file():
        raise SystemExit(f'[assistant-cohort] Project database does not exist: {original_db}')
    working_db = prep_root / '.assistant_cohort_project.sqlite'
    try:
        with open(original_db.with_name(original_db.name + '.lock'), 'a+') as lock:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
            try:
                _cohort_copy_db(original_db, working_db)
                project_db = Database(str(working_db))
                project_db.initialise()
                _cohort_register_and_import(project_db, batches, batch_ids, joint_inputs)
                if records:
                    joint_args = argparse.Namespace(**vars(args))
                    joint_args.cohort_manifest = None
                    joint_args.cohort_resume = None
                    joint_args.cohort_prepare_only = False
                    joint_args.command = 'assistant'
                    joint_args.sample_map = None
                    joint_args.pacbio_map = None
                    joint_args.fasta = str(joint_fasta)
                    joint_args.marker_qc = str(joint_qc)
                    joint_args.marker_source_type = 'SANGER_AB1'
                    joint_args.marker_review = str(decisions_path) if decisions_path else None
                    joint_args.partner_id = None
                    joint_args.cohort_qc_failures = str(failure_path)
                    joint_args.db = str(working_db)
                    joint_args.out = str(joint_root)
                    cmd_assistant(joint_args)
                else:
                    metadata_by_id = {r['source_id']: r for r in load_partner_sequencing_metadata(args.partner_metadata)}
                    for sequence_id, reason in failure_reasons.items():
                        row = metadata_by_id.get(sequence_id, {})
                        project_db.upsert_sequencing_metadata([{
                            'id': sequence_id, 'partner_id': row.get('partner_id') or sequence_id.split('_', 1)[0],
                            'dataset': batch_ids.get(sequence_id, args.dataset),
                            'selected_for_sequencing': row.get('selected_for_sequencing', False),
                            'selected_for_wgs': row.get('selected_for_wgs', False),
                            'source_id': sequence_id, 'source_file': args.partner_metadata,
                            'raw_selected_value': row.get('raw_selected_value', ''),
                            'raw_commitment_value': row.get('raw_commitment_value', ''),
                            'project_metadata': {
                                **(row.get('project_metadata') or {}),
                                **({'expected_taxon': row.get('expected_taxon')} if row.get('expected_taxon') else {}),
                                'expected_domain': row.get('expected_domain') or 'Bacteria',
                            },
                        }])
                        project_db.upsert_sequence_provenance([{
                            'sequence_id': sequence_id, 'marker_qc_class': 'FAIL_QC',
                            'marker_qc_recommendation': 'RESEQUENCE', 'marker_qc_reasons': reason,
                            'manual_review_status': 'NOT_REVIEWED', 'chimera_call': 'NOT_RUN',
                        }])
                        project_db.update_isolate_status(sequence_id, 'MARKER_QC_FAILED', detail=reason, source_file=str(failure_path))
                # Persist decisions for records that never entered the registry (rejected/unrescued failures).
                unrescued = [d for d in decisions if d['sequence_id'] not in included]
                if unrescued:
                    project_db.upsert_sequence_provenance([{
                        'sequence_id': d['sequence_id'], 'marker_qc_class': d['qc_class'] or 'FAIL_QC',
                        'marker_qc_recommendation': 'RESEQUENCE', 'marker_qc_reasons': failure_reasons.get(d['sequence_id'], ''),
                        'manual_review_status': 'REJECTED' if d['decision'] in {'REJECT', 'RESEQUENCE'} else 'NOT_REVIEWED',
                        'review_decision': d['decision'], 'reviewer': d['reviewer'], 'review_notes': d['notes'],
                        'reviewed_at': d['reviewed_at'], 'chimera_call': 'NOT_RUN',
                    } for d in unrescued])
                    record_review_decisions(project_db, unrescued, source_file=str(decisions_path or ''))
                with sqlite3.connect(str(working_db)) as conn:
                    check = conn.execute('PRAGMA integrity_check').fetchone()[0]
                if check != 'ok':
                    raise RuntimeError(f'cohort project database failed SQLite integrity_check: {check}')
                publish = original_db.with_name(original_db.name + '.publishing')
                publish.unlink(missing_ok=True)
                with sqlite3.connect(str(working_db)) as source, sqlite3.connect(str(publish)) as target:
                    source.backup(target)
                os.replace(publish, original_db)
                working_db.unlink(missing_ok=True)
            finally:
                fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
        manifest.add_output(joint_fasta, role='cohort_joint_marker_fasta')
        manifest.add_output(joint_qc, role='cohort_joint_marker_qc')
        manifest.add_output(review_audit_path, role='cohort_review_application_audit')
        if records:
            manifest.add_output(joint_root / 'run_manifest.json', role='joint_assistant_manifest')
            manifest.add_output(joint_root / '03_performance_review_hiring_panel' / 'assessment' / 'selection_summary.tsv', role='joint_selection_summary')
        status = 'COMPLETE_WITH_WARNINGS' if failure_reasons else 'COMPLETE'
        manifest.add_stage('assistant_cohort', status,
                           detail=f'{len(batch_ids)} isolate(s), {len(records)} marker(s) assessed, '
                                  f'{len(approved_failed & included)} approved failed marker(s), {len(failure_reasons)} excluded')
        manifest.finish(status)
        print(f'[assistant-cohort] ASSESS {status}: {len(records)} marker(s) assessed jointly; {len(failure_reasons)} excluded.\n'
              f'  Review application audit: {review_audit_path}\n'
              f'  Cohort report: {manifest.json_path}')
    except BaseException as exc:
        manifest.finish('FAILED', error=exc)
        raise


def cmd_assistant_cohort(args):
    """Mixed-partner cohort Assistant: PREPARE, then REVIEW AND ASSESS (resume)."""
    # Emit Python's native traceback if a C extension or native dependency
    # terminates the process with SIGSEGV/SIGABRT. This is diagnostic only.
    import faulthandler
    import sys
    if not faulthandler.is_enabled():
        faulthandler.enable(file=sys.stderr, all_threads=True)
    logging.getLogger(__name__).info('[ASSISTANT][DIAG] faulthandler enabled for native-crash diagnostics')
    if getattr(args, 'cohort_resume', None):
        root = Path(args.cohort_resume).expanduser().resolve()
        return _cohort_assess(args, root)
    root = Path(args.out).expanduser().resolve()
    original_db = Path(args.db).expanduser().resolve()
    if not original_db.is_file():
        raise SystemExit(f'[assistant-cohort] Project database does not exist: {original_db}')
    batches = _read_assistant_cohort_manifest(args.cohort_manifest)
    _cohort_prepare(args, root, batches)
    if getattr(args, 'cohort_prepare_only', False):
        return None
    return _cohort_assess(args, root)


def cmd_assistant(args):
    """Run the Assistant workflow with explicit scientific stage boundaries."""
    if getattr(args, 'cohort_manifest', None) or getattr(args, 'cohort_resume', None):
        return cmd_assistant_cohort(args)
    from branchmanager.it_desk import run_it_desk_checks, write_it_desk_report

    root = Path(args.out).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    manifest = RunManifest(root, 'assistant')
    for role, source in (
        ('sample_map', args.sample_map), ('marker_fasta', args.fasta), ('pacbio_map', args.pacbio_map),
        ('partner_metadata', args.partner_metadata), ('marker_qc', args.marker_qc),
        ('marker_review', args.marker_review),
        ('primary_reference', args.ref), ('reference_taxonomy', args.taxa),
        ('baseline_fasta', args.baseline_fasta), ('mwl', args.mwl),
        ('baseline_taxonomy', args.baseline_taxa_assignments),
        ('chimera_reference', args.chimera_ref),
    ):
        if source:
            manifest.add_input(source, role=role)
    for source in args.alt_ref or []:
        manifest.add_input(source, role='alternate_reference')
    for source in args.alt_taxa or []:
        if source:
            manifest.add_input(source, role='alternate_reference_taxonomy')
    try:
        it_desk_dir = root / '00_it_desk'
        checks = run_it_desk_checks(
            db_path=args.db,
            references=[args.ref, *(args.alt_ref or [])],
            output_dir=root,
            tree_method=args.tree_method,
        )
        it_desk_report = write_it_desk_report(it_desk_dir, checks)
        manifest.add_output(it_desk_report['tsv'], role='it_desk_report')
        if it_desk_report['status'] != 'PASS':
            raise RuntimeError(f'IT Desk checks failed; inspect {it_desk_report["tsv"]}')

        kinnex_dir = root / '01_pacbio_16s_import'
        kinnex_outputs = None
        if args.pacbio_map:
            from branchmanager.pipeline.kinnex import run_kinnex_import
            kinnex_outputs = run_kinnex_import(
                args.pacbio_map, kinnex_dir,
                min_read_length=args.pacbio_min_read_length, max_read_length=args.pacbio_max_read_length,
                min_mean_quality=args.pacbio_min_mean_quality,
                cluster_identity=args.pacbio_cluster_identity, min_reads=args.pacbio_min_reads,
                min_dominant_fraction=args.pacbio_min_dominant_fraction,
                library_prep=args.pacbio_library_prep,
                threads=args.pacbio_threads,
                strand=args.pacbio_vsearch_strand,
                query_cov=args.pacbio_vsearch_query_cov,
                maxaccepts=args.pacbio_vsearch_maxaccepts,
                maxrejects=args.pacbio_vsearch_maxrejects,
            )
            manifest.add_stage('pacbio_16s_import', 'COMPLETE', detail=f"{kinnex_outputs['accepted']}/{kinnex_outputs['total']} isolates accepted")

        onboarding_dir = root / '01_onboarding'
        onboarding_args = argparse.Namespace(
            sample_map=args.sample_map, fasta=kinnex_outputs['fasta'] if kinnex_outputs else args.fasta, partner_metadata=args.partner_metadata,
            read_dir=args.read_dir, primers=args.primers, partner_id=args.partner_id,
            dataset=args.dataset, out=str(onboarding_dir),
        )
        cmd_onboarding(onboarding_args)
        manifest.add_stage('onboarding', 'COMPLETE')

        paper_trail_dir = root / '02_paper_trail_merge_meeting'
        if args.sample_map:
            trace_args = argparse.Namespace(
                command='paper-trail', input=[], out=str(paper_trail_dir),
                sample_map=str(onboarding_dir / 'normalised_read_map.tsv'), primers=args.primers,
                primer_sequences=args.primer_sequences, trim_primers=args.trim_primers,
                min_quality=args.min_quality, min_length=args.min_length,
                min_read_length=args.min_read_length, min_mean_quality=args.min_mean_quality,
                mask_quality=args.mask_quality,
                max_read_expected_errors=args.max_read_expected_errors,
                max_output_expected_errors=args.max_output_expected_errors,
                warn_n_percent=args.warn_n_percent, max_n_percent=args.max_n_percent,
                warn_internal_low_quality_run=args.warn_internal_low_quality_run,
                max_internal_low_quality_run=args.max_internal_low_quality_run,
                secondary_peak_ratio=args.secondary_peak_ratio,
                max_mixed_peak_percent=args.max_mixed_peak_percent,
                mixed_peak_min_quality=args.mixed_peak_min_quality,
                max_conflict_density=args.max_conflict_density,
                max_ambiguous_overlap_conflicts_without_review=args.max_ambiguous_overlap_conflicts_without_review,
                quality_difference=args.quality_difference,
                allow_missing_quality=args.allow_missing_quality,
                min_overlap=args.min_overlap, min_overlap_identity=args.min_overlap_identity,
                screen_ref=args.ref, screen_taxa=args.taxa,
                threads=args.threads,
                max_report_image_height=args.max_report_image_height,
            )
            cmd_paper_trail(trace_args)
            marker_input = paper_trail_dir / 'assembled.fasta'
            marker_qc = paper_trail_dir / 'assembly_report.tsv'
            manifest.add_stage('paper_trail', 'COMPLETE')
            manifest.add_stage('merge_meeting', 'COMPLETE')
        elif kinnex_outputs:
            marker_input = onboarding_dir / 'normalised_input.fasta'
            marker_qc = Path(kinnex_outputs['marker_qc'])
            manifest.add_stage('paper_trail', 'SKIPPED', detail='Kinnex clustering replaced Sanger trace processing')
            manifest.add_stage('merge_meeting', 'SKIPPED', detail='Kinnex representative selected from dominant read cluster')
        else:
            if not args.marker_qc and not args.accept_unverified_marker_qc:
                raise RuntimeError(
                    'FASTA Assistant submissions require --marker-qc or the explicit '
                    '--accept-unverified-marker-qc audit flag'
                )
            marker_input = onboarding_dir / 'normalised_input.fasta'
            marker_qc = Path(args.marker_qc) if args.marker_qc else None
            manifest.add_stage('paper_trail', 'SKIPPED', detail='partner supplied FASTA')
            manifest.add_stage('merge_meeting', 'SKIPPED', detail='partner supplied FASTA')

        performance_review_dir = root / '03_performance_review_hiring_panel'
        performance_review_args = argparse.Namespace(
            command='performance-review', input=str(marker_input),
            partner_metadata=args.partner_metadata, db=args.db, dataset=args.dataset,
            out=str(performance_review_dir), ref=args.ref, taxa=args.taxa, ref_name=args.ref_name,
            alt_ref=args.alt_ref, alt_taxa=args.alt_taxa, alt_ref_name=args.alt_ref_name,
            main_ref=args.ref_name, baseline_fasta=args.baseline_fasta,
            baseline_dataset=args.baseline_dataset,
            baseline_tier=args.baseline_tier,
            baseline_taxa_assignments=args.baseline_taxa_assignments,
            baseline_skip_classify=False, baseline_colours=None, baseline_shorten_ids=False,
            novelty_baseline_datasets=[], mwl=args.mwl, mwl_sheet='MWL_V1', mwl_min_rank='p',
            taxa_assignments=None, previous_review=None, shorten_ids=False,
            min_len=800, max_n_percent=args.max_n_percent,
            sequence_domain=args.sequence_domain,
            phylum=None, target=None, force_rebuild=False, anchors=args.anchors,
            threads=args.threads, tree_method=args.tree_method,
            neighbourhood_format='png', pangenome_target=args.pangenome_target,
            candidate_set_size=args.candidate_set_size, user_colours=None,
            baseline_redundancy_identity=args.baseline_redundancy_identity,
            baseline_redundancy_min_query_coverage=args.baseline_redundancy_min_query_coverage,
            baseline_extension_min_identity=args.baseline_extension_min_identity,
            baseline_extension_min_query_coverage=args.baseline_extension_min_query_coverage,
            group_phyla=None, functional=None, draft_rumen_functions=False,
            marker_qc=str(marker_qc) if marker_qc else None,
            marker_review=args.marker_review,
            marker_source_type=(
                str(getattr(args, 'marker_source_type', 'auto') or 'auto').upper()
                if str(getattr(args, 'marker_source_type', 'auto') or 'auto').lower() != 'auto'
                else ('SANGER_AB1' if args.sample_map else 'PACBIO_HIFI' if args.pacbio_map else 'PARTNER_FASTA')
            ),
            accept_unverified_marker_qc=bool(args.accept_unverified_marker_qc),
            unverified_override_authorizer=getattr(args, 'unverified_override_authorizer', ''),
            unverified_override_reason=getattr(args, 'unverified_override_reason', ''),
            chimera_ref=args.chimera_ref or args.ref,
            skip_chimera_check=args.skip_chimera_check, collapse=False,
        )
        cmd_performance_review(performance_review_args)
        manifest.add_stage('performance_review', 'COMPLETE')
        manifest.add_stage('hiring_panel', 'COMPLETE')
        for role, path in (
            ('marker_fasta', marker_input),
            ('performance_review_dashboard', performance_review_dir / 'performance_review_dashboard.html'),
            ('selection_summary', performance_review_dir / 'assessment' / 'selection_summary.tsv'),
            ('sequencing_sets', performance_review_dir / 'assessment' / 'sequencing_sets.tsv'),
            ('pangenome_hierarchy', performance_review_dir / 'assessment' / 'pangenome_hierarchy.tsv'),
        ):
            manifest.add_output(path, role=role)
        missing = manifest.verify_required_outputs()
        if missing:
            raise RuntimeError('Assistant workflow outputs missing: ' + ', '.join(missing))
        manifest.finish('COMPLETE')
    except BaseException as exc:
        manifest.finish('FAILED', error=exc)
        raise
    print(
        '[assistant] Onboarding through Hiring Panel complete.\n'
        f'  Decision dashboard: {performance_review_dir / "performance_review_dashboard.html"}\n'
        f'  Sequencing sets  : {performance_review_dir / "assessment" / "sequencing_sets.tsv"}\n'
        f'  Pangenome ranks  : {performance_review_dir / "assessment" / "pangenome_hierarchy.tsv"}'
    )


def _cmd_quarterly_review_impl(args):
    """Run a project-wide post-coverage genome-selection round."""
    from datetime import datetime, timezone

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    _configure_logging(str(outdir))
    log = logging.getLogger(__name__)
    db = Database(args.db)
    db.initialise()

    # A final cumulative metadata file may be supplied without re-running the
    # expensive classification/tree workflow. It updates factual genome status
    # only; Quarterly Review recommendations never write back to that status.
    if getattr(args, 'partner_metadata', None):
        with db.connect() as conn:
            sequence_rows = conn.execute('SELECT id, dataset FROM sequences').fetchall()
        roles = db.get_dataset_roles()
        candidate_ids = [
            str(sequence_id) for sequence_id, dataset in sequence_rows
            if roles.get(str(dataset or ''), {}).get('role') != 'baseline'
        ]
        setattr(args, 'dataset', 'QuarterlyReview')
        _load_partner_metadata_for_run(args, db, str(outdir), {}, candidate_ids)

    round_id = getattr(args, 'round_id', None) or datetime.now(timezone.utc).strftime(
        'quarterly_review_%Y%m%dT%H%M%SZ'
    )
    assessment_paths = list(getattr(args, 'assessment', None) or [])
    if assessment_paths:
        try:
            imported_rows = quarterly_review.read_assessment_tables(assessment_paths)
        except Exception as exc:
            raise SystemExit(f'[QUARTERLY REVIEW] Failed to import assessment table: {exc}')
        if not imported_rows:
            raise SystemExit('[QUARTERLY REVIEW] The supplied assessment table(s) contained no sequence rows.')
        snapshot_id = f'quarterly-review-import:{round_id}'
        db.save_assessment_snapshot(
            snapshot_id,
            imported_rows,
            dataset='quarterly_review_import',
            source_path=';'.join(str(Path(path).resolve()) for path in assessment_paths),
        )
        log.info('[QUARTERLY REVIEW] Imported %d latest assessment rows into %s', len(imported_rows), snapshot_id)

    latest = db.get_latest_assessment_rows()
    if not latest:
        raise SystemExit(
            '[QUARTERLY REVIEW] No assessment snapshots are stored. Run Performance Review once with the current '
            'BranchManager version or provide --assessment path/to/sequence_assessment.tsv.'
        )

    tree_path = getattr(args, 'tree', None)
    if not tree_path and getattr(args, 'from_dir', None):
        tree_path = _find_tree_file_in_dir(args.from_dir)
    if tree_path and not Path(tree_path).exists():
        raise SystemExit(f'[QUARTERLY REVIEW] Tree file not found: {tree_path}')
    if not tree_path:
        log.warning(
            '[QUARTERLY REVIEW] No cumulative tree supplied; ranking will use marker identities and other '
            'assessment evidence without marginal patristic distance.'
        )

    alignment_path = getattr(args, 'alignment', None)
    if not alignment_path and getattr(args, 'from_dir', None):
        for candidate in (
            Path(args.from_dir) / 'tree' / 'current_alignment.fasta',
            Path(args.from_dir) / 'current_alignment.fasta',
        ):
            if candidate.exists():
                alignment_path = str(candidate)
                break
    if not alignment_path and tree_path:
        candidate = Path(tree_path).parent / 'current_alignment.fasta'
        if candidate.exists():
            alignment_path = str(candidate)
    if alignment_path and not Path(alignment_path).exists():
        raise SystemExit(f'[QUARTERLY REVIEW] Alignment file not found: {alignment_path}')

    parameters = {
        'GenomeBudget': int(args.genome_budget),
        'BackupsPerPrimary': int(args.backups_per_primary),
        'PangenomeTarget': int(args.pangenome_target),
        'BaselineRedundancyIdentity': float(args.baseline_redundancy_identity),
        'BaselineRedundancyMinQueryCoverage': float(args.baseline_redundancy_min_query_coverage),
        'BaselineExtensionMinIdentity': float(args.baseline_extension_min_identity),
        'BaselineExtensionMinQueryCoverage': float(args.baseline_extension_min_query_coverage),
        'IncludeModerateEvidence': bool(args.include_moderate_evidence),
        'Tree': str(Path(tree_path).resolve()) if tree_path else 'None',
        'Alignment': str(Path(alignment_path).resolve()) if alignment_path else 'None',
        'AssessmentRows': len(latest),
    }
    recommendations = quarterly_review.build_quarterly_review(
        list(latest.values()),
        db=db,
        tree_path=tree_path,
        alignment_path=alignment_path,
        genome_budget=args.genome_budget,
        backups_per_primary=args.backups_per_primary,
        pangenome_target=args.pangenome_target,
        baseline_redundancy_identity=args.baseline_redundancy_identity,
        baseline_redundancy_min_query_coverage=args.baseline_redundancy_min_query_coverage,
        baseline_extension_min_identity=args.baseline_extension_min_identity,
        baseline_extension_min_query_coverage=args.baseline_extension_min_query_coverage,
        include_moderate=args.include_moderate_evidence,
    )

    # Quarterly Review recommendations need round-specific figures. Reusing the most
    # recent Performance Review image would display stale PRIMARY/BACKUP stars.
    if tree_path:
        try:
            recommendation_by_id = {row['sequence_id']: row for row in recommendations}
            visual_rows = [quarterly_review.normalise_assessment_row(row) for row in latest.values()]
            for row in visual_rows:
                recommendation = recommendation_by_id.get(row['id'], {})
                role = str(recommendation.get('role') or '')
                row['sequencing_set_role'] = {
                    'PRIMARY': 'PRIMARY',
                    'BACKUP': 'BACKUP',
                    'ALREADY_SELECTED': 'COMMITTED',
                    'ALREADY_SEQUENCED': 'SEQUENCED',
                }.get(role, '')
                row['sequencing_set_rank'] = (
                    recommendation.get('round_rank', 'NA') if role in {'PRIMARY', 'BACKUP'} else 'NA'
                )
                if role == 'PRIMARY':
                    row['sequencing_set_badge'] = f'P{recommendation.get("round_rank", "")}'
                elif role == 'BACKUP':
                    row['sequencing_set_badge'] = (
                        f'B{recommendation.get("round_rank", "")}.{recommendation.get("backup_rank", "")}'
                    )
                row['selected_for_genome_sequencing'] = str(role == 'ALREADY_SELECTED')
                row['already_sequenced'] = str(role == 'ALREADY_SEQUENCED')
            visual_result = neighbourhood.generate_local_neighbourhood_visuals(
                tree_path=tree_path,
                assessment_rows=visual_rows,
                db=db,
                outdir=outdir / 'neighbourhoods',
                alignment_path=alignment_path,
                image_format=args.neighbourhood_format,
            )
            for row in visual_rows:
                recommendation = recommendation_by_id.get(row['id'])
                if recommendation is not None:
                    recommendation['local_tree_figure'] = row.get('local_neighbourhood_figure', 'NA')
            parameters['NeighbourhoodFigures'] = len(visual_result['figures'])
            parameters['NeighbourhoodFormat'] = args.neighbourhood_format
        except Exception as exc:
            log.warning('[QUARTERLY REVIEW] Could not generate round-specific neighbourhood figures: %s', exc)
    outputs = quarterly_review.write_quarterly_review_reports(outdir, round_id, recommendations, parameters)
    db.save_selection_round(round_id, recommendations, mode='quarterly_review', parameters=parameters)

    primary_count = sum(row.get('role') == 'PRIMARY' for row in recommendations)
    backup_count = sum(row.get('role') == 'BACKUP' for row in recommendations)
    log.info(
        '[QUARTERLY REVIEW] Round %s nominated %d primary and %d backup isolate(s)',
        round_id, primary_count, backup_count,
    )
    print(
        '[quarterly-review] Done.\n'
        f'  Round             : {round_id}\n'
        f'  Primary nominees  : {primary_count}\n'
        f'  Backup nominees   : {backup_count}\n'
        f'  Next genome set   : {outputs["selected"]}\n'
        f'  Full audit        : {outputs["summary"]}\n'
        f'  Round manifest    : {outputs["manifest"]}'
    )


def cmd_quarterly_review(args):
    import fcntl
    import sqlite3

    manifest = RunManifest(args.out, 'quarterly_review')
    original = Path(args.db).expanduser().resolve()
    if not original.is_file():
        raise SystemExit(f'[QUARTERLY REVIEW] Project database does not exist: {original}')
    manifest.add_input(original, role='project_database')
    for source in getattr(args, 'assessment', None) or []:
        manifest.add_input(source, role='assessment_import')
    for role, source in (
        ('partner_metadata', getattr(args, 'partner_metadata', None)),
        ('tree', getattr(args, 'tree', None)), ('alignment', getattr(args, 'alignment', None)),
    ):
        if source:
            manifest.add_input(source, role=role)
    staged_path = Path(args.out).resolve() / '.quarterly_review_project.sqlite'
    lock_path = original.with_name(original.name + '.lock')
    with open(lock_path, 'a+') as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        try:
            staged_path.unlink(missing_ok=True)
            with sqlite3.connect(str(original)) as source, sqlite3.connect(str(staged_path)) as target:
                source.backup(target)
            staged_args = argparse.Namespace(**vars(args))
            staged_args.db = str(staged_path)
            manifest.add_stage('quarterly_review', 'RUNNING')
            _cmd_quarterly_review_impl(staged_args)
            for role, name in (
                ('quarterly_review_summary', 'quarterly_review_summary.tsv'),
                ('next_genome_set', 'next_genome_set.tsv'),
                ('quarterly_review_parameters', 'quarterly_review_manifest.tsv'),
            ):
                manifest.add_output(Path(args.out) / name, role=role)
            missing = manifest.verify_required_outputs()
            if missing:
                raise RuntimeError('required Quarterly Review outputs are missing: ' + ', '.join(missing))
            manifest.add_stage('quarterly_review', 'COMPLETE')
            manifest.finish('COMPLETE')
            db = Database(str(staged_path))
            db.initialise()
            db.record_project_run(
                str(getattr(args, 'round_id', None) or f'quarterly-review:{manifest.data["started_at"]}'),
                'quarterly_review', 'COMPLETE', manifest_path=str(manifest.json_path),
                started_at=manifest.data['started_at'], completed_at=manifest.data['completed_at'],
            )
            publish_tmp = original.with_name(original.name + '.publishing')
            publish_tmp.unlink(missing_ok=True)
            with sqlite3.connect(str(staged_path)) as source, sqlite3.connect(str(publish_tmp)) as target:
                source.backup(target)
                if target.execute('PRAGMA integrity_check').fetchone()[0] != 'ok':
                    raise RuntimeError('staged Quarterly Review database failed integrity_check')
            os.replace(publish_tmp, original)
            staged_path.unlink(missing_ok=True)
        except BaseException as exc:
            manifest.finish('FAILED', error=exc)
            raise
        finally:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)


def cmd_dashboard(args):
    """Rebuild the interactive HTML dashboard from an existing run directory.

    This is deliberately a reporting-only operation: it reads existing TSV/JSON
    outputs and does not open or mutate the project database.
    """
    supplied_dir = Path(args.output_dir).expanduser().resolve()
    if not supplied_dir.is_dir():
        raise SystemExit(f"[DASHBOARD] Output directory does not exist: {supplied_dir}")

    # Accept either the run root or its assessment/ directory.
    if (supplied_dir / "selection_summary.tsv").is_file():
        assessment_dir = supplied_dir
        run_dir = supplied_dir.parent if supplied_dir.name.lower() == "assessment" else supplied_dir
    elif (supplied_dir / "assessment" / "selection_summary.tsv").is_file():
        run_dir = supplied_dir
        assessment_dir = supplied_dir / "assessment"
    else:
        raise SystemExit(
            "[DASHBOARD] Could not find selection_summary.tsv. Expected either "
            f"{supplied_dir / 'selection_summary.tsv'} or "
            f"{supplied_dir / 'assessment' / 'selection_summary.tsv'}."
        )

    selection_summary = assessment_dir / "selection_summary.tsv"
    sequence_assessment = assessment_dir / "sequence_assessment.tsv"
    pangenome_hierarchy = assessment_dir / "pangenome_hierarchy.tsv"
    novelty_metrics = assessment_dir / "novelty_metrics.tsv"
    run_manifest = run_dir / "run_manifest.json"
    genome_assessment = assessment_dir / "genome_collection_assessment.tsv"

    # The dashboard generator treats optional report inputs as empty when they
    # are absent. Log the omissions so a partial dashboard is not mistaken for
    # a complete report.
    for label, path in (
        ("sequence assessment", sequence_assessment),
        ("pangenome hierarchy", pangenome_hierarchy),
        ("novelty metrics", novelty_metrics),
        ("run manifest", run_manifest),
        ("genome collection assessment", genome_assessment),
    ):
        if not path.is_file():
            logging.getLogger(__name__).warning(
                "[DASHBOARD] Optional input is missing (%s): %s", label, path
            )

    from branchmanager.reporting import write_performance_review_dashboard

    dashboard = write_performance_review_dashboard(
        selection_summary,
        run_dir,
        sequence_assessment=sequence_assessment if sequence_assessment.is_file() else None,
        pangenome_hierarchy=pangenome_hierarchy if pangenome_hierarchy.is_file() else None,
        novelty_metrics=novelty_metrics if novelty_metrics.is_file() else None,
        run_manifest=run_manifest if run_manifest.is_file() else None,
        genome_assessment=genome_assessment if genome_assessment.is_file() else None,
        database=getattr(args, 'db', None),
    )
    logging.getLogger(__name__).info("[DASHBOARD] Rebuilt dashboard: %s", dashboard)
    print(f"Dashboard regenerated: {dashboard}")
    return dashboard


def cmd_regenerate_assessment(args):
    """Rebuild selection/pangenome reports from an existing assessment without rerunning biology."""
    assessment_dir = Path(args.assessment_dir).expanduser().resolve()
    source = assessment_dir / 'sequence_assessment.tsv'
    if not source.is_file():
        raise FileNotFoundError(f'Existing assessment TSV not found: {source}')
    db_path = Path(args.db).expanduser().resolve()
    if not db_path.is_file():
        raise FileNotFoundError(f'Existing project database not found: {db_path}')
    tree_path = Path(args.tree).expanduser().resolve() if args.tree else None
    if tree_path is None:
        for candidate in (assessment_dir / 'current_tree.nwk', assessment_dir.parent / 'current_tree.nwk',
                          assessment_dir.parent / 'tree' / 'current_tree.nwk'):
            if candidate.is_file():
                tree_path = candidate
                break
    if tree_path is not None and not tree_path.is_file():
        raise FileNotFoundError(f'Tree file not found: {tree_path}')

    def _snake_header(value):
        value = str(value or '').strip()
        # The assessment writer uses presentation headers (ID, PartnerID,
        # LocalNeighbourhoodFigure); selection functions consume snake_case.
        value = re.sub(r'([A-Z]+)([A-Z][a-z])', r'\1_\2', value)
        value = re.sub(r'([a-z0-9])([A-Z])', r'\1_\2', value)
        return value.replace('-', '_').replace(' ', '_').lower()

    with source.open('r', encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        rows = [
            {_snake_header(key): value for key, value in row.items() if key is not None}
            for row in reader
        ]
    if not rows:
        raise ValueError(f'Assessment TSV contains no data rows: {source}')

    # The original workflow creates candidate sets before appending isolates
    # that never produced a marker consensus. Preserve those placeholders in
    # final reports, but do not feed them into marker-candidate selection.
    marker_rows = [row for row in rows if str(row.get('no_marker_sequence', '')).strip().lower() not in {'1', 'true', 'yes'}]
    if not marker_rows:
        raise ValueError('No marker-bearing assessment rows are available for selection regeneration.')

    db = Database(str(db_path))
    # These reporting functions query the database; this command deliberately
    # does not update isolate statuses or write an assessment snapshot.
    sequencing_path = assessment_dir / 'sequencing_sets.tsv'
    hierarchy_path = assessment_dir / 'pangenome_hierarchy.tsv'
    _selection_sets.build_sequencing_sets(
        marker_rows, sequencing_path, tree_path=tree_path, db=db,
        pangenome_target=args.pangenome_target, candidate_set_size=args.candidate_set_size,
        baseline_redundancy_identity=args.baseline_redundancy_identity,
        baseline_redundancy_min_query_coverage=args.baseline_redundancy_min_query_coverage,
        baseline_extension_min_identity=args.baseline_extension_min_identity,
        baseline_extension_min_query_coverage=args.baseline_extension_min_query_coverage,
    )
    _selection_sets.write_pangenome_hierarchy_tsv(
        marker_rows, hierarchy_path, tree_path=tree_path, db=db,
        pangenome_target=args.pangenome_target,
    )

    # Carry the newly calculated selection fields back to the complete report
    # rows while retaining no-marker placeholders unchanged.
    updated = {str(r.get('id') or r.get('sequence_id') or ''): r for r in marker_rows}
    final_rows = [updated.get(str(r.get('id') or r.get('sequence_id') or ''), r) for r in rows]
    write_sequence_assessment_tsv(source, final_rows)
    write_selection_summary_tsv(assessment_dir / 'selection_summary.tsv', final_rows)
    write_baseline_hits_tsv(assessment_dir / 'baseline_hits.tsv', final_rows)

    logging.getLogger(__name__).info('[REGENERATE] Read %d assessment rows (%d marker-bearing).', len(rows), len(marker_rows))
    logging.getLogger(__name__).info('[REGENERATE] Wrote %s and %s', sequencing_path, hierarchy_path)
    logging.getLogger(__name__).info('[REGENERATE] Database statuses and assessment snapshots were not modified.')

    if not args.no_dashboard:
        from branchmanager.reporting import write_genome_collection_assessment, write_performance_review_dashboard
        run_dir = assessment_dir.parent
        genome_assessment = write_genome_collection_assessment(db, assessment_dir / 'genome_collection_assessment.tsv')
        dashboard = write_performance_review_dashboard(
            assessment_dir / 'selection_summary.tsv', run_dir,
            sequence_assessment=source,
            pangenome_hierarchy=hierarchy_path,
            novelty_metrics=assessment_dir / 'novelty_metrics.tsv',
            run_manifest=run_dir / 'run_manifest.json',
            genome_assessment=genome_assessment,
        )
        logging.getLogger(__name__).info('[REGENERATE] Rebuilt dashboard: %s', dashboard)

def main(argv=None):
    from branchmanager.cli.parser import build_parser
    parser = build_parser()
    selected_config = config_path_from_argv(argv)
    if selected_config:
        project_config = load_project_config(selected_config)
        apply_project_config_defaults(parser, project_config)
    args = parser.parse_args(argv)
    if args.command in {'ground-truth', 'validate'}:
        from branchmanager.validation import run_cli as _run_validation_cli
        result = _run_validation_cli(args)
        if result:
            sys.exit(result)
        return
    if args.command == 'assistant':
        # A cohort manifest is itself the complete input declaration. Ignore any
        # per-batch input path inherited from a reusable config profile.
        if getattr(args, 'cohort_manifest', None) or getattr(args, 'cohort_resume', None):
            args.sample_map = args.fasta = args.pacbio_map = None
        input_count = sum(bool(getattr(args, key, None)) for key in ('sample_map', 'fasta', 'pacbio_map', 'cohort_manifest', 'cohort_resume'))
        if input_count != 1:
            parser.error('assistant requires exactly one of --sample-map, --fasta, --pacbio-map, --cohort-manifest, or --cohort-resume, supplied on the command line or in --config')
        if getattr(args, 'cohort_prepare_only', False) and not args.cohort_manifest:
            parser.error('--cohort-prepare-only requires --cohort-manifest')
        if getattr(args, 'cohort_resume', None) and not args.out:
            args.out = args.cohort_resume
        if not args.dataset:
            parser.error('assistant requires --dataset or a dataset in --config')
        if not args.out:
            parser.error('assistant requires --out or an output_directory in --config')
    if args.command == 'quarterly-review' and args.genome_budget is None:
        parser.error('quarterly-review requires --genome-budget or a project selection_budget in --config')
    if selected_config and args.command == 'filing-cabinet':
        args.configured_baselines = project_config.get('baselines', [])
        if not args.configured_baselines and (not args.fasta or not args.dataset):
            parser.error('filing-cabinet requires --fasta and --dataset, or a project config with a baselines list')
    if args.command == 'filing-cabinet':
        if getattr(args, 'configured_baselines', None):
            cmd_filing_cabinet_configured_baselines(args)
        else:
            cmd_filing_cabinet(args)
    elif args.command == 'performance-review':
        cmd_performance_review(args)
    elif args.command == 'quarterly-review':
        cmd_quarterly_review(args)
    elif args.command == 'label-maker':
        cmd_label_maker(args)
    elif args.command == 'org-chart':
        cmd_org_chart(args)
    elif args.command == 'background-check':
        cmd_background_check(args)
    elif args.command == 'paper-trail':
        cmd_paper_trail(args)
    elif args.command == 'sanger-studio':
        from branchmanager.sanger_studio import serve
        serve(port=args.port, open_browser=not args.no_browser)
    elif args.command == 'mailroom':
        cmd_mailroom(args)
    elif args.command == 'interview':
        cmd_interview(args)
    elif args.command == 'onboarding':
        cmd_onboarding(args)
    elif args.command == 'pacbio-16s-import':
        cmd_pacbio_16s_import(args)
    elif args.command == 'status-meeting':
        _cmd_project_import(args, genomes=False)
    elif args.command == 'records-update':
        _cmd_project_import(args, genomes=True)
    elif args.command == 'add-external-genomes':
        cmd_add_external_genomes(args)
    elif args.command == 'add-maldi-identifications':
        cmd_add_maldi_identifications(args)
    elif args.command == 'marker-evidence-migrate':
        cmd_marker_evidence_migrate(args)
    elif args.command == 'register-marker-evidence':
        cmd_register_marker_evidence(args)
    elif args.command == 'import-genome-16s':
        cmd_import_genome_16s(args)
    elif args.command == 'import-genome-evidence':
        cmd_import_genome_evidence(args)
    elif args.command == 'import-genome-results':
        cmd_import_genome_results(args)
    elif args.command == 'marker-review-apply':
        cmd_marker_review_apply(args)
    elif args.command == 'import-partner-workbook':
        cmd_import_partner_workbook(args)
    elif args.command == 'exit-interview':
        cmd_exit_interview(args)
    elif args.command == 'annual-report':
        cmd_annual_report(args)
    elif args.command == 'it-desk':
        cmd_it_desk(args)
    elif args.command == 'dashboard':
        cmd_dashboard(args)
    elif args.command == 'regenerate-assessment':
        cmd_regenerate_assessment(args)
    elif args.command == 'assistant':
        cmd_assistant(args)
    elif getattr(args, '_tool', None) is not None:
        tool = args._tool
        if tool.__name__.endswith('id_sanitation'):
            sys.exit(tool.run(args))
        elif tool.__name__.endswith('benchmark_selection'):
            tool.run_cli(args, args._tool_parser)
        else:
            result = tool.run_cli(args)
            if isinstance(result, int) and result:
                sys.exit(result)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()