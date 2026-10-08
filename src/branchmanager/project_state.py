"""Operational isolate and genome state for BranchManager's rolling project."""

from __future__ import annotations

import csv
import gzip
import hashlib
from pathlib import Path
from typing import Iterable, Optional


ISOLATE_STATES = (
    'RECEIVED',
    'TRACE_REVIEW',
    'MARKER_QC_PASSED',
    'PROPOSED',
    'SAB_APPROVED',
    'DNA_EXTRACTION_PENDING',
    'DNA_EXTRACTION_FAILED',
    'LIBRARY_PENDING',
    'SEQUENCED',
    'GENOME_QC_FAILED',
    'GENOME_QC_PASSED',
    'WITHDRAWN',
)

TRUE_VALUES = {'1', 'true', 't', 'yes', 'y', 'pass', 'passed'}
FALSE_VALUES = {'0', 'false', 'f', 'no', 'n', 'fail', 'failed'}


def _normalise(name: object) -> str:
    return str(name or '').strip().lstrip('\ufeff').lower().replace('-', '_').replace(' ', '_')


def _open(path: str | Path):
    return gzip.open(path, 'rt', newline='') if str(path).lower().endswith('.gz') else open(path, newline='')


def _read_rows(path: str | Path) -> list[dict[str, str]]:
    with _open(path) as handle:
        sample = handle.read(4096)
        handle.seek(0)
        delimiter = '\t' if str(path).lower().endswith(('.tsv', '.tsv.gz', '.txt', '.txt.gz')) else ','
        try:
            delimiter = csv.Sniffer().sniff(sample, delimiters=',\t').delimiter
        except Exception:
            pass
        return [
            {_normalise(key): str(value or '').strip() for key, value in row.items() if key is not None}
            for row in csv.DictReader(handle, delimiter=delimiter)
        ]


def _value(row: dict[str, str], *names: str) -> str:
    for name in names:
        value = row.get(_normalise(name), '')
        if value:
            return value
    return ''


def _optional_float(value: object) -> Optional[float]:
    if value in (None, '', 'NA', 'N/A', 'None'):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _optional_bool(value: object) -> Optional[bool]:
    text = str(value or '').strip().lower()
    if not text:
        return None
    if text in TRUE_VALUES:
        return True
    if text in FALSE_VALUES:
        return False
    return None


def _sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def normalise_state(value: object, *, default: str = 'RECEIVED') -> str:
    state = str(value or default).strip().upper().replace('-', '_').replace(' ', '_')
    aliases = {
        'APPROVED': 'SAB_APPROVED',
        'DNA_FAILED': 'DNA_EXTRACTION_FAILED',
        'DNA_PENDING': 'DNA_EXTRACTION_PENDING',
        'QC_PASSED': 'GENOME_QC_PASSED',
        'QC_FAILED': 'GENOME_QC_FAILED',
        'GENOME_AVAILABLE': 'GENOME_QC_PASSED',
    }
    state = aliases.get(state, state)
    if state not in ISOLATE_STATES:
        raise ValueError(f'Unknown isolate status {value!r}; expected one of {", ".join(ISOLATE_STATES)}')
    return state


def import_status_updates(db, path: str | Path) -> list[dict]:
    rows = _read_rows(path)
    known = set(db.get_all_ids())
    results = []
    for line_number, row in enumerate(rows, start=2):
        sequence_id = _value(row, 'sequence_id', 'isolate_id', 'sample_id', 'id')
        state_text = _value(row, 'status', 'operational_status', 'isolate_status')
        detail = _value(row, 'detail', 'status_detail', 'notes')
        if not sequence_id:
            results.append({'line': line_number, 'sequence_id': '', 'result': 'REJECTED', 'detail': 'missing sequence_id'})
            continue
        if sequence_id not in known:
            results.append({'line': line_number, 'sequence_id': sequence_id, 'result': 'REJECTED', 'detail': 'sequence ID is not in the project database'})
            continue
        try:
            state = normalise_state(state_text)
        except ValueError as exc:
            results.append({'line': line_number, 'sequence_id': sequence_id, 'result': 'REJECTED', 'detail': str(exc)})
            continue
        db.update_isolate_status(sequence_id, state, detail=detail, source_file=str(Path(path).resolve()))
        results.append({'line': line_number, 'sequence_id': sequence_id, 'result': 'IMPORTED', 'detail': state})
    return results


def import_genome_results(
    db, path: str | Path, *, min_completeness: float = 90.0,
    max_contamination: float = 5.0,
) -> list[dict]:
    """Import completed genome evidence without equating sequencing with QC success."""
    rows = _read_rows(path)
    known = set(db.get_all_ids())
    accepted = []
    results = []
    for line_number, row in enumerate(rows, start=2):
        sequence_id = _value(row, 'sequence_id', 'isolate_id', 'sample_id', 'id')
        genome_id = _value(row, 'genome_id', 'assembly_id', 'accession', 'genome_accession')
        if not sequence_id:
            results.append({'line': line_number, 'sequence_id': '', 'genome_id': genome_id, 'result': 'REJECTED', 'detail': 'missing sequence_id'})
            continue
        if sequence_id not in known:
            results.append({'line': line_number, 'sequence_id': sequence_id, 'genome_id': genome_id, 'result': 'REJECTED', 'detail': 'sequence ID is not in the project database'})
            continue
        if not genome_id:
            results.append({'line': line_number, 'sequence_id': sequence_id, 'genome_id': '', 'result': 'REJECTED', 'detail': 'missing genome_id/accession'})
            continue

        status_text = _value(row, 'genome_status', 'status') or 'SEQUENCED'
        try:
            status = normalise_state(status_text, default='SEQUENCED')
        except ValueError as exc:
            results.append({'line': line_number, 'sequence_id': sequence_id, 'genome_id': genome_id, 'result': 'REJECTED', 'detail': str(exc)})
            continue
        if status not in {'SEQUENCED', 'GENOME_QC_FAILED', 'GENOME_QC_PASSED'}:
            results.append({'line': line_number, 'sequence_id': sequence_id, 'genome_id': genome_id, 'result': 'REJECTED', 'detail': 'genome_status must be SEQUENCED, GENOME_QC_FAILED, or GENOME_QC_PASSED'})
            continue

        completeness = _optional_float(_value(row, 'completeness', 'genome_completeness'))
        contamination = _optional_float(_value(row, 'contamination', 'genome_contamination'))
        explicit_pass = _optional_bool(_value(row, 'genome_qc_pass', 'qc_pass', 'passed_qc'))
        metrics_available = completeness is not None and contamination is not None
        metrics_pass = bool(
            metrics_available
            and completeness >= float(min_completeness)
            and contamination <= float(max_contamination)
        )
        if explicit_pass is None and status == 'SEQUENCED' and metrics_available:
            qc_pass = metrics_pass
            status = 'GENOME_QC_PASSED' if metrics_pass else 'GENOME_QC_FAILED'
            qc_basis = f'automatic:{min_completeness:g}% completeness/{max_contamination:g}% contamination'
        else:
            qc_pass = status == 'GENOME_QC_PASSED' if explicit_pass is None else explicit_pass
            qc_basis = 'explicit_status_or_qc_flag'
        if status == 'GENOME_QC_FAILED' and qc_pass:
            results.append({'line': line_number, 'sequence_id': sequence_id, 'genome_id': genome_id, 'result': 'REJECTED', 'detail': 'GENOME_QC_FAILED conflicts with genome_qc_pass=yes'})
            continue
        if status == 'GENOME_QC_PASSED' and not qc_pass:
            results.append({'line': line_number, 'sequence_id': sequence_id, 'genome_id': genome_id, 'result': 'REJECTED', 'detail': 'GENOME_QC_PASSED conflicts with genome_qc_pass=no'})
            continue
        if qc_pass and metrics_available and not metrics_pass:
            results.append({
                'line': line_number, 'sequence_id': sequence_id, 'genome_id': genome_id,
                'result': 'REJECTED',
                'detail': (
                    f'QC pass conflicts with completeness={completeness}, contamination={contamination}; '
                    f'thresholds are >= {min_completeness} and <= {max_contamination}'
                ),
            })
            continue

        genome_path = _value(row, 'genome_path', 'assembly_path', 'fasta_path')
        genome_sha256 = ''
        if genome_path:
            candidate = Path(genome_path).expanduser()
            if not candidate.is_absolute():
                candidate = Path(path).resolve().parent / candidate
            if not candidate.is_file():
                results.append({'line': line_number, 'sequence_id': sequence_id, 'genome_id': genome_id, 'result': 'REJECTED', 'detail': f'genome file not found: {candidate}'})
                continue
            genome_path = str(candidate.resolve())
            genome_sha256 = _sha256(candidate)

        accepted.append({
            'genome_id': genome_id,
            'sequence_id': sequence_id,
            'accession': _value(row, 'accession', 'genome_accession'),
            'genome_status': status,
            'genome_qc_pass': bool(qc_pass),
            'completeness': completeness,
            'contamination': contamination,
            'gtdb_taxonomy': _value(row, 'gtdb_taxonomy', 'gtdbtk_taxonomy', 'taxonomy'),
            'ani_cluster': _value(row, 'ani_cluster', 'species_cluster', 'ani_species'),
            'genome_path': genome_path,
            'genome_sha256': genome_sha256,
            'notes': '; '.join(filter(None, (_value(row, 'notes', 'detail'), f'qc_basis={qc_basis}'))),
            'source_file': str(Path(path).resolve()),
        })
        results.append({'line': line_number, 'sequence_id': sequence_id, 'genome_id': genome_id, 'result': 'IMPORTED', 'detail': status})

    db.upsert_genome_records(accepted)
    for row in accepted:
        sequence_id = row['sequence_id']
        state = 'GENOME_QC_PASSED' if row['genome_qc_pass'] else row['genome_status']
        db.update_isolate_status(sequence_id, state, detail=row.get('notes', ''), source_file=str(Path(path).resolve()))
        if row['genome_qc_pass']:
            current = db.get_sequencing_metadata_for_ids([sequence_id]).get(sequence_id, {})
            db.upsert_sequencing_metadata([{
                'id': sequence_id,
                'partner_id': current.get('partner_id', ''),
                'dataset': current.get('dataset', ''),
                'selected_for_sequencing': current.get('selected_for_sequencing', False),
                'selected_for_wgs': True,
                'source_id': sequence_id,
                'source_file': str(Path(path).resolve()),
                'raw_selected_value': 'genome_qc_pass',
                'raw_commitment_value': current.get('raw_commitment_value', ''),
                'project_metadata': current.get('project_metadata', {}),
            }])
    return results


def import_external_genomes(
    db, path: str | Path, *, marker_fasta: str | Path | None = None,
    default_dataset: str = 'ExternalGenomeBaseline', default_tier: str = 'secondary',
) -> list[dict]:
    """Register genomes acquired outside the candidate-selection workflow.

    A marker sequence is optional. Genome-only records receive an empty marker
    slot, remain available for genome-coverage/taxonomy context, and are
    excluded from marker similarity searches and trees.
    """
    from branchmanager.taxonomy import canonicalise_sequence_id
    from branchmanager.utils.fasta import read_fasta

    rows = _read_rows(path)
    markers = {}
    if marker_fasta:
        for header, sequence in read_fasta(str(marker_fasta)):
            markers[canonicalise_sequence_id(header)] = sequence
    results = []
    sequences_by_dataset = {}
    tax_entries = []
    genome_entries = []
    metadata_entries = []
    role_entries = {}
    known = db.get_all_ids()
    existing_metadata = db.get_sequencing_metadata_for_ids()
    for line_number, row in enumerate(rows, start=2):
        sequence_id = _value(row, 'sequence_id', 'isolate_id', 'sample_id', 'id')
        genome_id = _value(row, 'genome_id', 'assembly_id', 'accession', 'genome_accession')
        if not sequence_id or not genome_id:
            results.append({'line': line_number, 'sequence_id': sequence_id, 'genome_id': genome_id,
                            'result': 'REJECTED', 'detail': 'sequence_id and genome_id/accession are required'})
            continue
        dataset = _value(row, 'dataset', 'collection', 'baseline_dataset') or default_dataset
        tier = (_value(row, 'baseline_tier', 'tier') or default_tier).strip().lower()
        if tier not in {'priority', 'secondary'}:
            results.append({'line': line_number, 'sequence_id': sequence_id, 'genome_id': genome_id,
                            'result': 'REJECTED', 'detail': 'baseline_tier must be priority or secondary'})
            continue
        marker = markers.get(canonicalise_sequence_id(sequence_id), '')
        if sequence_id not in known:
            sequences_by_dataset.setdefault(dataset, []).append((sequence_id, marker))
            role_entries[dataset] = tier
            known.add(sequence_id)
        elif marker:
            with db.connect() as conn:
                conn.execute(
                    "UPDATE sequences SET sequence=?, length=? WHERE id=? AND (sequence IS NULL OR sequence='')",
                    (marker, len(marker), sequence_id),
                )
        taxonomy = _value(row, 'gtdb_taxonomy', 'gtdbtk_taxonomy', 'taxonomy')
        old = existing_metadata.get(sequence_id, {})
        if taxonomy:
            confidence = _optional_float(_value(row, 'taxonomy_confidence', 'confidence'))
            tax_entries.append((sequence_id, taxonomy, confidence, old.get('dataset') or dataset))
        metadata_entries.append({
            'id': sequence_id, 'partner_id': _value(row, 'partner_id', 'partner') or old.get('partner_id', 'external'),
            'dataset': old.get('dataset') or dataset,
            'selected_for_sequencing': bool(old.get('selected_for_sequencing')),
            'selected_for_wgs': True, 'source_id': sequence_id, 'source_file': str(Path(path).resolve()),
            'raw_selected_value': 'external_genome_available',
            'raw_commitment_value': old.get('raw_commitment_value', ''),
            'project_metadata': {
                key: value for key, value in row.items()
                if key not in {
                    'sequence_id', 'isolate_id', 'sample_id', 'id', 'genome_id',
                    'assembly_id', 'accession', 'genome_accession', 'dataset',
                    'collection', 'baseline_dataset', 'baseline_tier', 'tier',
                    'gtdb_taxonomy', 'gtdbtk_taxonomy', 'taxonomy',
                    'taxonomy_confidence', 'confidence', 'genome_status',
                    'genome_qc_pass', 'qc_pass', 'completeness', 'genome_completeness',
                    'contamination', 'genome_contamination', 'ani_cluster',
                    'species_cluster', 'genome_path', 'assembly_path', 'notes', 'detail',
                }
            },
        })
        genome_status = _value(row, 'genome_status') or 'SEQUENCED'
        genome_qc = _optional_bool(_value(row, 'genome_qc_pass', 'qc_pass'))
        if genome_qc is None:
            genome_qc = genome_status.strip().upper() == 'GENOME_QC_PASSED'
        genome_entries.append({
            'genome_id': genome_id, 'sequence_id': sequence_id,
            'accession': _value(row, 'accession', 'genome_accession') or genome_id,
            'genome_status': genome_status,
            'genome_qc_pass': genome_qc,
            'completeness': _optional_float(_value(row, 'completeness', 'genome_completeness')),
            'contamination': _optional_float(_value(row, 'contamination', 'genome_contamination')),
            'gtdb_taxonomy': taxonomy,
            'ani_cluster': _value(row, 'ani_cluster', 'species_cluster'),
            'genome_path': _value(row, 'genome_path', 'assembly_path'),
            'notes': _value(row, 'notes', 'detail'), 'source_file': str(Path(path).resolve()),
        })
        results.append({'line': line_number, 'sequence_id': sequence_id, 'genome_id': genome_id,
                        'result': 'IMPORTED', 'detail': 'genome available; marker ' + ('attached' if marker else 'not supplied')})

    for dataset, tier in role_entries.items():
        db.upsert_dataset_role(dataset, 'baseline', genomes_available=True, baseline_tier=tier)
    if sequences_by_dataset:
        # Empty sequence values are intentional genome-only records.
        for dataset, subset in sequences_by_dataset.items():
            db.insert_sequences(subset, dataset=dataset, baseline_tier=role_entries[dataset])
    if tax_entries:
        db.insert_taxonomy(tax_entries)
    if metadata_entries:
        db.upsert_sequencing_metadata(metadata_entries)
    if genome_entries:
        db.upsert_genome_records(genome_entries)
    return results


def write_import_report(path: str | Path, rows: Iterable[dict]) -> str:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = ['line', 'sequence_id', 'genome_id', 'result', 'detail']
    with open(output, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    return str(output)
