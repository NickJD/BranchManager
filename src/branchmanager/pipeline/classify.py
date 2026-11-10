"""
classify.py — Taxonomic classification for BranchManager.

Uses VSEARCH --usearch_global against a reference database (GreenGenes2,
SILVA, or custom) to assign taxonomy to query sequences.

Key design decisions
--------------------
- Reference anchors (BRANCHMANAGER_REF_* headers) are excluded from output.
- The taxonomy lookup map is built once at load time for performance.
- Taxonomy classification uses a high-stringency production identity/coverage
  ladder (98/95 down to 90/85), followed by a separately labelled diagnostic
  ladder ending at the legacy 80/70 sensitivity; diagnostic-only assignments
  remain explicitly tentative.
- Taxonomy is assigned by rank-wise dominant consensus across near-best hits;
  conflicts and assignment support are retained as explicit review evidence.
  Novelty detection (novelty.py) runs separate nearest-hit and density searches.
- Output columns are database-labeled. Identity is percent identity across the
  aligned region. Query/target coverage report how much of each sequence aligned.
  SourceClassifierConfidence is populated only when a source taxonomy table
  explicitly supplies a confidence value; it is never inferred from identity.
"""

from __future__ import annotations

import logging
import os
import re
import shlex
from pathlib import Path
from typing import Optional

from branchmanager.pipeline.evidence_policy import classification_evidence_interpretation
from branchmanager.utils.reference import ensure_uncompressed_reference
from branchmanager.taxonomy import parse_reference_header_taxonomy, reference_lookup_keys
from branchmanager.taxonomy_io import iter_taxonomy_assignment_rows
from branchmanager.utils.fasta import read_fasta, write_fasta
from branchmanager.utils.subprocess import run_cmd
from branchmanager.pipeline.vsearch_audit import run_vsearch_audited
from branchmanager.pipeline.tree import is_ref_anchor

logger = logging.getLogger(__name__)


VSEARCH_USERFIELDS = 'query+target+id+alnlen+mism+gaps+qlo+qhi+tlo+thi+ql+tl+qcov+tcov'


def _vsearch_safe_query(input_fasta, outdir) -> tuple[str, dict[str, str]]:
    """Write a copy of *input_fasta* with whitespace-free IDs.

    vsearch truncates query labels at the first whitespace, which breaks IDs
    such as 'JJ A11'. Returns the temp FASTA path and a safe-ID -> original-ID map.
    """
    safe_path = os.path.join(str(outdir), 'vsearch_query_safe.fasta')
    id_map: dict[str, str] = {}
    records = []
    for i, (header, seq) in enumerate(read_fasta(input_fasta)):
        safe = f'q{i}'
        id_map[safe] = header
        records.append((safe, seq))
    write_fasta(records, safe_path)
    return safe_path, id_map


def _parse_vsearch_match(parts: list[str]) -> Optional[dict]:
    if len(parts) < 3:
        return None
    try:
        identity = float(parts[2])
    except (TypeError, ValueError):
        return None
    result = {
        'query': parts[0], 'target': parts[1], 'identity': identity,
        'alignment_length': None, 'mismatches': None, 'gaps': None,
        'query_length': None, 'target_length': None,
        'query_coverage': None, 'target_coverage': None,
    }
    if len(parts) >= 12:
        try:
            result.update({
                'alignment_length': int(float(parts[3])),
                'mismatches': int(float(parts[4])),
                'gaps': int(float(parts[5])),
                'query_length': int(float(parts[10])),
                'target_length': int(float(parts[11])),
            })
            if len(parts) >= 14:
                result['query_coverage'] = max(0.0, min(100.0, float(parts[12])))
                result['target_coverage'] = max(0.0, min(100.0, float(parts[13])))
            else:
                # Compatibility for older cached 12-column matches. Alignment length
                # includes gap columns, so it is only a fallback and must never yield
                # biologically impossible coverage above 100%.
                if result['query_length'] > 0:
                    result['query_coverage'] = min(
                        100.0, 100.0 * result['alignment_length'] / result['query_length'],
                    )
                if result['target_length'] > 0:
                    result['target_coverage'] = min(
                        100.0, 100.0 * result['alignment_length'] / result['target_length'],
                    )
        except (TypeError, ValueError, ZeroDivisionError):
            pass
    return result


def _alignment_columns(hit: Optional[dict]) -> list[str]:
    if not hit:
        return ['NA'] * 7
    def number(key, digits=2):
        value = hit.get(key)
        return 'NA' if value is None else f'{float(value):.{digits}f}'
    def integer(key):
        value = hit.get(key)
        return 'NA' if value is None else str(int(value))
    return [
        number('query_coverage'), number('target_coverage'), integer('alignment_length'),
        integer('query_length'), integer('target_length'), integer('mismatches'), integer('gaps'),
    ]


# Re-export consolidated ID normalisation helpers for backward compatibility
from branchmanager.utils.text import norm_id as _norm_id
from branchmanager.utils.text import canon_id as _canon_id


# Taxonomy map loading

def _load_taxa_map(taxa_tsv: str, wanted_ids: Optional[set[str]] = None) -> dict[str, tuple[str, Optional[float]]]:
    """
    Load a taxonomy TSV/CSV into a lookup dict.

    Accepts two formats:
      FeatureID  Taxon  Confidence       (QIIME2 / GreenGenes2 style)
      ID  Taxon  Confidence              (plain two-column)
      and equivalent comma-delimited CSV files, including .gz variants.

    Builds THREE index keys per entry for resilient lookup:
      - raw ID
      - _norm_id(ID)
      - _canon_id(ID)

    Returns {key: (taxon_string, confidence_float_or_None)}
    """
    taxa_map: dict[str, tuple[str, Optional[float]]] = {}
    if wanted_ids is not None and not wanted_ids:
        logger.info('[CLASSIFY] No matched reference IDs; skipping taxonomy-table scan: %s', taxa_tsv)
        return taxa_map

    try:
        # Log start of taxonomy loading
        try:
            size_mb = os.path.getsize(taxa_tsv) / (1024 * 1024)
            logger.info("[CLASSIFY] Loading taxonomy map from %s (%.1f MB)...", taxa_tsv, size_mb)
        except Exception:
            logger.info("[CLASSIFY] Loading taxonomy map from %s...", taxa_tsv)

        rows_loaded = 0
        wanted_keys = set()
        if wanted_ids is not None:
            for wanted in wanted_ids:
                wanted_keys.update(reference_lookup_keys(wanted))
                wanted_keys.update((_norm_id(wanted), _canon_id(wanted)))
        for row in iter_taxonomy_assignment_rows(taxa_tsv):
            fid = str(row.get('id', '')).strip()
            if wanted_ids is not None:
                candidate_keys = set(reference_lookup_keys(fid)) | {_norm_id(fid), _canon_id(fid)}
                if not candidate_keys.intersection(wanted_keys):
                    continue
            tax = str(row.get('taxonomy', '')).strip()
            conf = row.get('confidence')
            _add_taxa_entry(fid, tax, conf, taxa_map)
            rows_loaded += 1

        logger.info("[CLASSIFY] Loaded %d taxonomy assignment rows from %s",
                    rows_loaded, taxa_tsv)
    except Exception as e:
        logger.warning("[CLASSIFY] Failed to load taxa file %s: %s", taxa_tsv, e)

    return taxa_map


def _load_taxa_map_from_reference_fasta(ref_fasta: str, wanted_ids: Optional[set[str]] = None) -> tuple[dict[str, tuple[str, Optional[float]]], list[dict[str, str]]]:
    """Build taxonomy lookups directly from reference FASTA headers.

    Supports GTDB 16S headers where the accession is followed by a lineage.
    Returns (taxa_map, warning_rows).
    """
    taxa_map: dict[str, tuple[str, Optional[float]]] = {}
    warning_rows = []
    if wanted_ids is not None and not wanted_ids:
        return taxa_map, warning_rows
    total = 0
    parsed = 0
    wanted_keys = set()
    if wanted_ids is not None:
        for wanted in wanted_ids:
            wanted_keys.update(reference_lookup_keys(wanted))
            wanted_keys.update((_norm_id(wanted), _canon_id(wanted)))
    for header, _ in read_fasta(ref_fasta):
        ref_id, taxonomy = parse_reference_header_taxonomy(header)
        if wanted_ids is not None:
            header_keys = set(reference_lookup_keys(header)) | set(reference_lookup_keys(ref_id or ''))
            if not header_keys.intersection(wanted_keys):
                continue
        total += 1
        if not ref_id or not taxonomy:
            warning_rows.append({
                'category': 'REFERENCE_HEADER_WITHOUT_TAXONOMY',
                'subject_id': ref_id or str(header).split()[0],
                'detail': 'Reference FASTA header did not contain a parseable lineage string.',
            })
            continue
        parsed += 1
        for key in list(reference_lookup_keys(header)) + list(reference_lookup_keys(ref_id)):
            if key and key not in taxa_map:
                taxa_map[key] = (taxonomy, None)
    if total and parsed / float(total) < 0.8:
        warning_rows.insert(0, {
            'category': 'LOW_HEADER_TAXONOMY_PARSE_RATE',
            'subject_id': f'{parsed}/{total}',
            'detail': 'Fewer than 80% of reference FASTA headers contained parseable taxonomy. Check header format or provide --taxa explicitly.',
        })
    return taxa_map, warning_rows


def _add_taxa_entry(fid: str, tax: str, conf: object, taxa_map: dict) -> None:
    fid = str(fid or '').strip()
    tax = str(tax or '').strip()
    if not fid or not tax:
        return
    try:
        conf_value: Optional[float] = float(conf) if conf is not None else None
    except (ValueError, TypeError):
        conf_value = None

    for key in list(reference_lookup_keys(fid)) + [_norm_id(fid), _canon_id(fid)]:
        if key and key not in taxa_map:
            taxa_map[key] = (tax, conf_value)


def _lookup_tax(sid: str, taxa_map: dict) -> tuple[Optional[str], Optional[float]]:
    """
    Look up a VSEARCH hit ID in the taxa_map.
    Tries: raw → norm → canon → token-level splitting.
    Returns (taxon, confidence) or (None, None).
    """
    if not taxa_map:
        return None, None

    for key in list(reference_lookup_keys(sid)) + [_norm_id(sid), _canon_id(sid)]:
        if key and key in taxa_map:
            return taxa_map[key]

    # Token-level splitting for complex IDs like RS_GCF_000001405.40
    tokens = re.split(r"[|_:\-.]+", sid)
    for tok in tokens:
        if not tok:
            continue
        for key in (tok, _norm_id(tok), _canon_id(tok)):
            if key and key in taxa_map:
                return taxa_map[key]

    return None, None


# Reference decompression

def _ensure_uncompressed(ref_fasta: str, outdir: str, out_name: str = "ref_uncompressed.fasta") -> str:
    """Return a validated, source-specific uncompressed reference cache."""
    return ensure_uncompressed_reference(ref_fasta, outdir, out_name=out_name)


def _iter_taxa_file_rows(taxa_tsv: str):
    for row in iter_taxonomy_assignment_rows(taxa_tsv):
        values = [str(row.get('id', '')), str(row.get('taxonomy', ''))]
        conf = row.get('confidence')
        if conf is not None:
            values.append(str(conf))
        yield values


def validate_reference_taxonomy_consistency(ref_fasta: Optional[str], taxa_tsv: Optional[str]):
    """Return warning rows about mismatches between reference FASTA and taxonomy table.

    Skips validation for large reference files (>500 MB) as it's too expensive.
    """
    warnings = []
    if not ref_fasta or not taxa_tsv:
        return warnings

    # Skip validation for large reference files - too expensive to parse
    try:
        ref_size = os.path.getsize(ref_fasta)
        if ref_size > 500 * 1024 * 1024:  # 500 MB threshold
            logger.info(
                '[CLASSIFY] Skipping reference/taxonomy consistency validation for large reference '
                'file (%.1f MB). This check is too expensive for production databases.',
                ref_size / (1024 * 1024)
            )
            return warnings
    except Exception:
        logger.debug('Could not check reference file size; proceeding with validation')

    ref_ids = []
    ref_key_to_id = {}
    try:
        logger.info('[CLASSIFY] Validating reference/taxonomy consistency (this may take a minute for large files)...')
        for header, _ in read_fasta(ref_fasta):
            ref_ids.append(header)
            for key in reference_lookup_keys(header):
                ref_key_to_id.setdefault(key, header)
    except Exception as e:
        return [{'category': 'REFERENCE_READ_FAILED', 'subject_id': str(ref_fasta), 'detail': str(e)}]

    taxa_ids = []
    taxa_key_to_id = {}
    try:
        for parts in _iter_taxa_file_rows(taxa_tsv):
            if len(parts) < 2:
                continue
            tid = parts[0].strip()
            if not tid:
                continue
            taxa_ids.append(tid)
            for key in reference_lookup_keys(tid):
                taxa_key_to_id.setdefault(key, tid)
    except Exception as e:
        return [{'category': 'TAXONOMY_READ_FAILED', 'subject_id': str(taxa_tsv), 'detail': str(e)}]

    matched_ref = [rid for rid in ref_ids if any(key in taxa_key_to_id for key in reference_lookup_keys(rid))]
    unmatched_ref = [rid for rid in ref_ids if rid not in matched_ref]
    unmatched_taxa = [tid for tid in taxa_ids if not any(key in ref_key_to_id for key in reference_lookup_keys(tid))]

    if ref_ids:
        overlap = len(matched_ref) / float(len(ref_ids))
        if overlap < 0.8:
            warnings.append({
                'category': 'LOW_REFERENCE_TAXONOMY_OVERLAP',
                'subject_id': f'{len(matched_ref)}/{len(ref_ids)}',
                'detail': f'Only {overlap:.1%} of reference FASTA IDs matched taxonomy table IDs; taxonomy source may be inconsistent with the classifier reference.',
            })
    for rid in unmatched_ref[:25]:
        warnings.append({'category': 'REFERENCE_ID_MISSING_TAXONOMY', 'subject_id': rid, 'detail': 'Reference FASTA entry had no matching taxonomy row.'})
    for tid in unmatched_taxa[:25]:
        warnings.append({'category': 'TAXONOMY_ID_MISSING_REFERENCE', 'subject_id': tid, 'detail': 'Taxonomy table entry had no matching reference FASTA record.'})
    return warnings


def write_reference_taxonomy_warnings(outdir: str, warning_rows):
    path = Path(outdir) / 'taxonomy_input_warnings.tsv'
    with open(path, 'w') as fh:
        fh.write('Category\tSubjectID\tDetail\n')
        for row in warning_rows:
            fh.write(f"{row.get('category', 'WARNING')}\t{row.get('subject_id', '')}\t{row.get('detail', '')}\n")
    return str(path)


# Main classification function

def run_classification(
    input_fasta: str,
    outdir: str,
    ref_fasta: Optional[str] = None,
    taxa_tsv: Optional[str] = None,
    threads: Optional[int] = None,
    db_name: Optional[str] = None,
) -> str:
    """Compatibility entry point using the same adaptive classifier as multi-DB runs."""
    if not ref_fasta:
        raise ValueError("ref_fasta must be provided to run_classification")
    output, _results = run_classification_single(
        input_fasta=input_fasta, outdir=outdir, ref_fasta=ref_fasta,
        taxa_tsv=taxa_tsv, threads=threads, db_name='main',
        column_db_name=db_name or 'primary',
    )
    return output


# Multi-database classification

def _derive_db_name(fasta_path: str) -> str:
    """Derive a short, filesystem-safe name from a reference FASTA path.

    Strips common suffixes (_reps, _ssu, _16s, _rna, _seqs, _NR<digits>) before
    sanitising and truncating to 24 characters.
    """
    import re as _re
    base = Path(fasta_path).stem
    base = _re.sub(r'_(reps?|ssu_?reps?|16[sS]|rna|nr\d*|seqs?)$', '', base, flags=_re.IGNORECASE)
    safe = _re.sub(r'[^A-Za-z0-9]+', '_', base).strip('_')[:24]
    return safe if safe else Path(fasta_path).stem[:16]


def _classification_header(db_name: str) -> list[str]:
    """Build database-labeled classifier columns for a taxonomy TSV."""
    import re as _re
    safe = _re.sub(r'[^A-Za-z0-9]+', '_', str(db_name)).strip('_') or 'primary'
    metrics = [
        'QueryCoverage', 'TargetCoverage', 'AlignmentLength', 'QueryLength',
        'TargetLength', 'Mismatches', 'Gaps',
    ]
    return [
        'ID', f'BestHit_{safe}', f'Identity_{safe}', f'Taxon_{safe}',
        f'SourceClassifierConfidence_{safe}', *[f'{metric}_{safe}' for metric in metrics],
        f'AssignmentStatus_{safe}', f'TaxonomySupportPct_{safe}',
        f'TaxonomyLeadPct_{safe}', f'AssignedRank_{safe}',
        f'EvidenceTier_{safe}', f'TaxonomyCandidateCount_{safe}',
        f'AlignmentEvidenceStatus_{safe}', f'AssignmentEvidenceLevel_{safe}',
        f'EvidenceInterpretation_{safe}',
    ]



# High-stringency adaptive 16S classification. Each tier is a minimum identity
# and query-coverage requirement. Only queries with fewer than the required
# number of distinct candidate references advance to the next tier.
CLASSIFICATION_SEARCH_TIERS = (
    (98.0, 95.0),
    (97.0, 95.0),
    (96.0, 95.0),
    (95.0, 95.0),
    (93.0, 90.0),
    (90.0, 85.0),
)
CLASSIFICATION_MIN_CANDIDATES = 3
CLASSIFICATION_MAX_ACCEPTS = 20
CLASSIFICATION_MAX_REJECTS = 256
CLASSIFICATION_MAX_HITS = 20
# Diagnostic searches are run only for queries not sufficiently represented by
# the production ladder. The final tier reproduces the previous background
# check's sensitivity (80% identity / 70% query coverage) and hit limits.
DIAGNOSTIC_SEARCH_TIERS = (
    (90.0, 80.0),
    (90.0, 75.0),
    (85.0, 80.0),
    (85.0, 75.0),
    (80.0, 70.0),
)
LEGACY_DIAGNOSTIC_MAX_ACCEPTS = 10
LEGACY_DIAGNOSTIC_MAX_HITS = 10

# Conservative rank-wise consensus across near-best 16S reference hits.
TAXONOMY_CONSENSUS_MIN_SUPPORT = 0.70
TAXONOMY_CONSENSUS_MIN_LEAD = 0.20
TAXONOMY_CONSENSUS_MIN_REFERENCES = 3


def _run_adaptive_classification_search(
    safe_query: str, query_id_map: dict[str, str], ref_fasta: str,
    matches_path: str, outdir: str, db_name: str, threads: Optional[int],
) -> dict[str, list[dict]]:
    """Run production and diagnostic VSEARCH tiers for unresolved queries.

    The production ladder is evaluated first. Queries with fewer than
    CLASSIFICATION_MIN_CANDIDATES distinct references then enter the diagnostic
    ladder, which ends with the legacy 80% identity / 70% query-coverage search.
    Diagnostic-only hits remain explicitly labelled and cannot be mistaken for
    production-supported classifications.
    """
    import csv
    import tempfile

    records = list(read_fasta(safe_query))
    sequences = dict(records)
    hits_by_query: dict[str, dict[str, dict]] = {}
    summary = []
    thread_flag = f" --threads {int(threads)}" if threads and int(threads) > 0 else ""
    out_path = Path(matches_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=f'.{out_path.stem}_tiers_', dir=str(out_path.parent)) as temp_dir:
        root = Path(temp_dir)
        unresolved = set(sequences)
        for tier_index, (identity, coverage) in enumerate(CLASSIFICATION_SEARCH_TIERS, 1):
            if not unresolved:
                break
            queries_this_tier = len(unresolved)
            query_path = root / f'queries_tier_{tier_index}.fasta'
            raw_path = root / f'hits_tier_{tier_index}.tsv'
            write_fasta([(qid, sequences[qid]) for qid in sorted(unresolved)], str(query_path))
            cmd = (
                f"vsearch --usearch_global {shlex.quote(str(query_path))} "
                f"--db {shlex.quote(str(ref_fasta))} --id {identity / 100.0:.2f} "
                f"--query_cov {coverage / 100.0:.2f} --strand both "
                f"--userout {shlex.quote(str(raw_path))} --userfields {VSEARCH_USERFIELDS} "
                f"--maxaccepts {CLASSIFICATION_MAX_ACCEPTS} "
                f"--maxrejects {CLASSIFICATION_MAX_REJECTS} "
                f"--maxhits {CLASSIFICATION_MAX_HITS}{thread_flag}"
            )
            logger.info(
                '[CLASSIFY][%s][LADDER] tier=%d identity=%.1f%% qcov=%.1f%% queries=%d',
                db_name, tier_index, identity, coverage, len(unresolved),
            )
            run_vsearch_audited(
                cmd, audit_path=out_path.parent / 'vsearch_search_audit.tsv',
                stage='A', purpose=f'taxonomic_classification:{db_name}:tier{tier_index}',
                query_fasta=query_path, reference_fasta=ref_fasta, runner=run_cmd,
            )
            tier_hits = 0
            if raw_path.exists():
                with raw_path.open(encoding='utf-8', errors='replace') as handle:
                    for line in handle:
                        hit = _parse_vsearch_match(line.rstrip('\n').split('\t'))
                        if hit is not None:
                            hit['_raw_line'] = line.rstrip('\n')
                        if not hit or hit['query'] not in unresolved:
                            continue
                        qid = hit['query']
                        target = str(hit.get('target') or '')
                        if not target:
                            continue
                        # Keep the strongest evidence for a repeated target.
                        previous = hits_by_query.setdefault(qid, {}).get(target)
                        rank = (hit.get('identity') or 0, hit.get('query_coverage') or 0, hit.get('alignment_length') or 0)
                        old_rank = (previous.get('identity') or 0, previous.get('query_coverage') or 0, previous.get('alignment_length') or 0) if previous else None
                        if old_rank is None or rank > old_rank:
                            hits_by_query[qid][target] = hit
                        tier_hits += 1
            resolved_now = 0
            for qid in list(unresolved):
                if len(hits_by_query.get(qid, {})) >= CLASSIFICATION_MIN_CANDIDATES:
                    unresolved.remove(qid)
                    resolved_now += 1
            summary.append({
                'database': db_name, 'tier': tier_index, 'minimum_identity_pct': identity,
                'minimum_query_coverage_pct': coverage, 'queries_searched': queries_this_tier,
                'raw_hits': tier_hits, 'queries_resolved_at_tier': resolved_now,
                'queries_remaining': len(unresolved),
            })
            logger.info(
                '[CLASSIFY][%s][LADDER] tier=%d raw_hits=%d newly_resolved=%d remaining=%d',
                db_name, tier_index, tier_hits, resolved_now, len(unresolved),
            )

        # Diagnostic rescue ladder for queries without a qualifying production
        # hit. These lower thresholds are explicitly diagnostic, not part of the
        # production assignment policy. Queries stop once three distinct rescue
        # references have been collected, matching the production candidate rule.
        diagnostic_unresolved = set(unresolved)
        for diagnostic_index, (identity, coverage) in enumerate(DIAGNOSTIC_SEARCH_TIERS, 1):
            if not diagnostic_unresolved:
                break
            query_path = root / f'queries_diagnostic_{diagnostic_index}.fasta'
            raw_path = root / f'hits_diagnostic_{diagnostic_index}.tsv'
            write_fasta(
                [(qid, sequences[qid]) for qid in sorted(diagnostic_unresolved)],
                str(query_path),
            )
            # The final 80/70 tier deliberately matches the previous background
            # check's maxaccepts/maxhits settings as well as its identity/qcov.
            is_legacy_tier = (identity, coverage) == (80.0, 70.0)
            max_accepts = LEGACY_DIAGNOSTIC_MAX_ACCEPTS if is_legacy_tier else CLASSIFICATION_MAX_ACCEPTS
            max_hits = LEGACY_DIAGNOSTIC_MAX_HITS if is_legacy_tier else CLASSIFICATION_MAX_HITS
            cmd = (
                f"vsearch --usearch_global {shlex.quote(str(query_path))} "
                f"--db {shlex.quote(str(ref_fasta))} --id {identity / 100.0:.2f} "
                f"--query_cov {coverage / 100.0:.2f} --strand both "
                f"--userout {shlex.quote(str(raw_path))} --userfields {VSEARCH_USERFIELDS} "
                f"--maxaccepts {max_accepts} --maxrejects {CLASSIFICATION_MAX_REJECTS} "
                f"--maxhits {max_hits}{thread_flag}"
            )
            tier_label = f'DIAGNOSTIC_{int(identity)}_{int(coverage)}'
            logger.info(
                '[CLASSIFY][%s][DIAGNOSTIC] tier=%d identity=%.1f%% qcov=%.1f%% queries=%d',
                db_name, diagnostic_index, identity, coverage, len(diagnostic_unresolved),
            )
            run_vsearch_audited(
                cmd, audit_path=out_path.parent / 'vsearch_search_audit.tsv',
                stage='A_DIAGNOSTIC',
                purpose=f'taxonomic_classification_diagnostic:{db_name}:{int(identity)}_{int(coverage)}',
                query_fasta=query_path, reference_fasta=ref_fasta, runner=run_cmd,
            )
            tier_hits = 0
            if raw_path.exists():
                with raw_path.open(encoding='utf-8', errors='replace') as handle:
                    for line in handle:
                        hit = _parse_vsearch_match(line.rstrip('\n').split('\t'))
                        if not hit or hit['query'] not in diagnostic_unresolved:
                            continue
                        hit['_raw_line'] = line.rstrip('\n')
                        hit['_evidence_tier'] = tier_label
                        qid, target = hit['query'], str(hit.get('target') or '')
                        if not target:
                            continue
                        previous = hits_by_query.setdefault(qid, {}).get(target)
                        rank = (
                            hit.get('identity') or 0,
                            hit.get('query_coverage') or 0,
                            hit.get('alignment_length') or 0,
                        )
                        old_rank = (
                            previous.get('identity') or 0,
                            previous.get('query_coverage') or 0,
                            previous.get('alignment_length') or 0,
                        ) if previous else None
                        if old_rank is None or rank > old_rank:
                            hits_by_query[qid][target] = hit
                        tier_hits += 1
            resolved_now = 0
            for qid in list(diagnostic_unresolved):
                diagnostic_count = sum(
                    1 for hit in hits_by_query.get(qid, {}).values()
                    if str(hit.get('_evidence_tier', '')).startswith('DIAGNOSTIC_')
                )
                if diagnostic_count >= CLASSIFICATION_MIN_CANDIDATES:
                    diagnostic_unresolved.remove(qid)
                    resolved_now += 1
            summary.append({
                'database': db_name, 'tier': tier_label,
                'minimum_identity_pct': identity,
                'minimum_query_coverage_pct': coverage,
                'queries_searched': len(diagnostic_unresolved) + resolved_now,
                'raw_hits': tier_hits,
                'queries_resolved_at_tier': resolved_now,
                'queries_remaining': len(diagnostic_unresolved),
            })
            logger.info(
                '[CLASSIFY][%s][DIAGNOSTIC] tier=%d raw_hits=%d newly_resolved=%d remaining=%d',
                db_name, diagnostic_index, tier_hits, resolved_now, len(diagnostic_unresolved),
            )

    # Write one combined, headerless VSEARCH userout for compatibility with
    # existing downstream parsers. Restore original query IDs.
    with out_path.open('w', encoding='utf-8') as handle:
        for qsafe, _seq in records:
            for hit in hits_by_query.get(qsafe, {}).values():
                raw_line = hit.get('_raw_line')
                if raw_line:
                    handle.write(raw_line + '\n')

    summary_path = out_path.with_name(f'classify_search_summary_{re.sub(r"[^A-Za-z0-9]+", "_", db_name)}.tsv')
    fields = ['database', 'tier', 'minimum_identity_pct', 'minimum_query_coverage_pct', 'queries_searched', 'raw_hits', 'queries_resolved_at_tier', 'queries_remaining']
    with summary_path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(summary)
    standard_with_hits = sum(
        any(not str(hit.get('_evidence_tier', '')).startswith('DIAGNOSTIC_') for hit in hits_by_query.get(q, {}).values())
        for q in sequences
    )
    diagnostic_only = sum(
        bool(hits_by_query.get(q)) and
        not any(not str(hit.get('_evidence_tier', '')).startswith('DIAGNOSTIC_') for hit in hits_by_query.get(q, {}).values())
        for q in sequences
    )
    no_hits_after_diagnostic = sum(not bool(hits_by_query.get(q)) for q in sequences)
    logger.info(
        '[CLASSIFY][%s][LADDER] complete queries=%d standard_hits=%d diagnostic_only=%d no_hits_after_diagnostic=%d summary=%s',
        db_name, len(records), standard_with_hits, diagnostic_only,
        no_hits_after_diagnostic, summary_path,
    )

    # Convert safe query IDs to the original IDs for the caller.
    return {query_id_map.get(qid, qid): list(targets.values()) for qid, targets in hits_by_query.items()}

def _taxonomy_parts(taxonomy: object) -> list[str]:
    if taxonomy is None:
        return []
    return [part.strip() for part in str(taxonomy).split(';') if part.strip()]


def _dominant_taxonomy_consensus(candidates: list[dict], taxa_map: dict) -> dict:
    """Return the deepest defensible rank and retain conflicts at every rank.

    Conflicts below the reported rank do not invalidate a broader assignment.
    Strict confidence still requires the existing support, lead and reference
    diversity thresholds at every rank through species. A plurality supported
    lineage is retained at the deepest rank with >50% support and a positive
    lead; unresolved lower ranks are recorded in the conflict ledger.
    """
    ranked = []
    for hit in candidates:
        tax, _confidence = _lookup_tax(hit.get('target', ''), taxa_map)
        parts = _taxonomy_parts(tax)
        if parts:
            ranked.append((hit, parts))
    if not ranked:
        return {'taxonomy': None, 'status': 'UNRESOLVED', 'support': 0.0,
                'lead': 0.0, 'rank': '', 'conflicts': [], 'near_best_count': 0,
                'strict_rank': '', 'candidate_taxonomy': None}

    selected, conflicts, rank_evidence = [], [], []
    parent_subset = ranked
    strict_rank = ''
    rank_names = ('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species')
    max_depth = max(len(parts) for _hit, parts in ranked)

    for index in range(max_depth):
        available = [(hit, parts) for hit, parts in parent_subset if len(parts) > index]
        if not available:
            break
        counts, refs = {}, {}
        for hit, parts in available:
            lineage = parts[index]
            counts[lineage] = counts.get(lineage, 0) + 1
            refs.setdefault(lineage, set()).add(str(hit.get('target') or ''))

        top_lineage = next((parts[index] for _hit, parts in available), None)
        winner_count = max(counts.values())
        tied = {taxon for taxon, count in counts.items() if count == winner_count}
        winner = top_lineage if top_lineage in tied else sorted(tied)[0]
        ordered = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
        runner_count = max((count for taxon, count in counts.items() if taxon != winner), default=0)
        denominator = len(available)
        support = winner_count / denominator
        lead = (winner_count - runner_count) / denominator
        distinct_refs = len({ref for ref in refs[winner] if ref})
        rank_name = rank_names[index] if index < len(rank_names) else f'rank_{index + 1}'
        strict = (support >= TAXONOMY_CONSENSUS_MIN_SUPPORT
                  and lead >= TAXONOMY_CONSENSUS_MIN_LEAD
                  and distinct_refs >= TAXONOMY_CONSENSUS_MIN_REFERENCES)
        rank_evidence.append({
            'rank': rank_name, 'taxon': winner, 'support': support,
            'lead': lead, 'distinct_refs': distinct_refs, 'strict': strict,
            'contested': len(ordered) > 1,
        })
        if strict:
            strict_rank = rank_name
        if len(ordered) > 1:
            runner = max(
                (taxon for taxon, count in counts.items() if taxon != winner and count == runner_count),
                key=lambda taxon: next((i for i, (_hit, parts) in enumerate(available)
                                        if parts[index] == taxon), len(available)),
            )
            conflicts.append({
                'rank': rank_name, 'dominant_taxon': winner,
                'dominant_support_pct': 100.0 * support,
                'runner_up_taxon': runner,
                'runner_up_support_pct': 100.0 * runner_count / denominator,
                'lead_pct_points': 100.0 * lead,
                'distinct_supporting_references': distinct_refs,
                'candidate_count': denominator,
            })
        selected.append(winner)
        parent_subset = [(hit, parts) for hit, parts in available if parts[index] == winner]

    # Report the deepest rank with an actual majority and positive lead.
    supported_depth = 0
    for index, evidence in enumerate(rank_evidence):
        if evidence['support'] >= 0.5 and evidence['lead'] > 0:
            supported_depth = index + 1
        else:
            break
    if supported_depth == 0:
        supported_depth = 1  # retain the top-hit domain as a broad provisional call
    reported = selected[:supported_depth]
    reported_evidence = rank_evidence[supported_depth - 1]
    reported_rank = reported_evidence['rank']
    # Use the actual rank-evidence order rather than the fixed standard
    # taxonomy tuple. This also supports additional ranks represented as
    # rank_8, rank_9, etc.
    rank_positions = {
        evidence['rank']: index
        for index, evidence in enumerate(rank_evidence)
    }

    conflicts_at_or_above = [
        item for item in conflicts
        if rank_positions.get(item['rank'], len(rank_evidence)) < supported_depth
        or (
            supported_depth == 1
            and rank_positions.get(item['rank'], len(rank_evidence)) == 1
        )
    ]
    if reported_rank == 'species' and all(item['strict'] for item in rank_evidence[:supported_depth]):
        status = 'CONFIDENT'
    elif conflicts_at_or_above:
        status = 'CONTESTED'
    else:
        status = 'PROVISIONAL'
    return {
        'taxonomy': ';'.join(reported) if reported else None,
        'candidate_taxonomy': ';'.join(selected) if selected else None,
        'status': status,
        'support': reported_evidence['support'],
        'lead': reported_evidence['lead'],
        'rank': reported_rank,
        'conflicts': conflicts,
        'near_best_count': len(ranked),
        'strict_rank': strict_rank,
    }


def _diagnostic_hit_reason(candidates: list[dict]) -> str:
    """Explain why no diagnostic candidate meets the production 90%/85% floor."""
    if not candidates:
        return 'NO_HIT_AT_DIAGNOSTIC_FLOOR'
    has_identity = any(float(h.get('identity') or 0.0) >= 90.0 for h in candidates)
    has_coverage = any(float(h.get('query_coverage') or 0.0) >= 85.0 for h in candidates)
    has_joint = any(
        float(h.get('identity') or 0.0) >= 90.0
        and float(h.get('query_coverage') or 0.0) >= 85.0
        for h in candidates
    )
    if has_joint:
        return 'MEETS_PRODUCTION_FLOOR_DIAGNOSTIC_ONLY'
    if has_identity and not has_coverage:
        return 'QUERY_COVERAGE_LIMITED'
    if has_coverage and not has_identity:
        return 'IDENTITY_LIMITED'
    return 'IDENTITY_AND_COVERAGE_LIMITED'


def run_classification_single(
    input_fasta: str,
    outdir: str,
    ref_fasta: str,
    taxa_tsv: Optional[str] = None,
    threads: Optional[int] = None,
    db_name: str = 'main',
    column_db_name: Optional[str] = None,
) -> tuple[str, dict]:
    """Classify *input_fasta* against a single reference database.

    Parameters
    ----------
    db_name : Explicit label for this database. Named databases write
              ``taxonomy_<db_name>.tsv`` / ``matches_<db_name>.tsv``.
              The special ``'main'`` label retains legacy filenames for direct callers.

    Returns
    -------
    (tsv_path, results_dict) where
    results_dict = {qid: (hit_id, pct_identity, taxon, source_classifier_confidence, alignment metrics, taxonomy status and evidence interpretation)}
    """
    import re as _re

    safe_db_name = _re.sub(r'[^A-Za-z0-9]+', '_', str(db_name)).strip('_') or 'primary'
    if db_name == 'main':
        # Legacy direct-call compatibility; multi-database runs use explicit names.
        tsv_name = 'taxonomy.tsv'
        match_name = 'matches.tsv'
    else:
        tsv_name = f'taxonomy_{safe_db_name}.tsv'
        match_name = f'matches_{safe_db_name}.tsv'

    output = os.path.join(outdir, tsv_name)
    matches_path = os.path.join(outdir, match_name)

    # Use a db-specific decompressed filename so multiple gzipped references
    # used in the same run do not overwrite each other's cached copy.
    if db_name == 'main':
        decomp_name = 'ref_uncompressed.fasta'
    else:
        decomp_name = f'ref_uncompressed_{safe_db_name}.fasta'
    ref_to_use = _ensure_uncompressed(ref_fasta, outdir, out_name=decomp_name)

    import time
    import resource
    started = time.monotonic()
    try:
        input_bytes = os.path.getsize(input_fasta)
    except OSError:
        input_bytes = -1
    try:
        ref_bytes = os.path.getsize(ref_to_use)
    except OSError:
        ref_bytes = -1
    logger.info('[CLASSIFY][%s][DIAG] START pid=%d input=%s input_bytes=%d ref=%s ref_bytes=%d taxa=%s threads=%s',
                db_name, os.getpid(), input_fasta, input_bytes, ref_to_use, ref_bytes,
                taxa_tsv or 'FASTA_HEADERS', threads)

    # Search first at high stringency. Taxonomy is loaded afterwards and can be
    # restricted to the reference IDs actually returned by VSEARCH.
    safe_query, query_id_map = _vsearch_safe_query(input_fasta, outdir)
    logger.info('[CLASSIFY][%s][DIAG] SAFE_QUERY_READY records=%d path=%s', db_name, len(query_id_map), safe_query)
    logger.info('[CLASSIFY][%s][DIAG] VSEARCH_LADDER_START output=%s', db_name, matches_path)
    hits_by_query = _run_adaptive_classification_search(
        safe_query, query_id_map, ref_to_use, matches_path, outdir, db_name, threads,
    )
    logger.info('[CLASSIFY][%s][DIAG] VSEARCH_LADDER_COMPLETE elapsed=%.2fs output_exists=%s', db_name, time.monotonic() - started, os.path.exists(matches_path))
    all_query_ids = [h for h, _ in read_fasta(input_fasta) if not is_ref_anchor(h)]
    matched_target_ids = {hit['target'] for candidates in hits_by_query.values() for hit in candidates}
    taxa_map: dict = {}
    if taxa_tsv:
        logger.info('[CLASSIFY][%s] Loading taxonomy only for %d matched reference IDs from %s', db_name, len(matched_target_ids), taxa_tsv)
        taxa_map = _load_taxa_map(taxa_tsv, wanted_ids=matched_target_ids)
    else:
        try:
            taxa_map, warning_rows = _load_taxa_map_from_reference_fasta(ref_to_use, wanted_ids=matched_target_ids)
            if warning_rows:
                warn_path = write_reference_taxonomy_warnings(outdir, warning_rows)
                logger.warning('[CLASSIFY][%s] Reference-header taxonomy warnings → %s', db_name, warn_path)
        except Exception as e:
            logger.warning('[CLASSIFY][%s] Could not parse matched reference taxonomy: %s', db_name, e)
    logger.info('[CLASSIFY][%s][DIAG] TAXA_MAP_READY entries=%d elapsed=%.2fs', db_name, len(taxa_map), time.monotonic() - started)

    # Retain multiple candidate hits, rank them deterministically, and preserve
    # supporting evidence. VSEARCH's first accepted hit is not definitive.

    def _rank_hit(hit):
        return (
            float(hit.get('identity') or 0.0),
            float(hit.get('query_coverage') or 0.0),
            int(hit.get('alignment_length') or 0),
            str(hit.get('target') or ''),
        )

    best_hits: dict[str, dict] = {}
    consensus_by_query: dict[str, dict] = {}
    support_path = os.path.join(outdir, f"classification_supporting_hits_{safe_db_name}.tsv")
    conflicts_path = os.path.join(outdir, f"taxonomy_conflicts_{safe_db_name}.tsv")
    with open(support_path, "w", encoding='utf-8') as support_fh:
        support_fh.write("Query\tRank\tTarget\tIdentity\tQueryCoverage\tTargetCoverage\tAlignmentLength\tTaxonomy\tNearBest\tEvidenceTier\n")
        for qid, candidates in hits_by_query.items():
            candidates.sort(key=_rank_hit, reverse=True)
            standard_candidates = [
                h for h in candidates
                if not str(h.get('_evidence_tier', '')).startswith('DIAGNOSTIC_')
            ]
            # Diagnostic hits may be retained as supporting/audit evidence,
            # but they must never displace a production-ladder hit or alter
            # the production-supported taxonomy consensus for that query.
            assignment_candidates = standard_candidates or candidates
            best_hits[qid] = assignment_candidates[0]
            best_identity = float(assignment_candidates[0].get('identity') or 0.0)
            near_best = [
                h for h in assignment_candidates
                if float(h.get('identity') or 0.0) >= best_identity - 1.0
            ]
            consensus_by_query[qid] = _dominant_taxonomy_consensus(near_best, taxa_map)
            if not standard_candidates:
                # A tentative lineage can be useful, but it must never be
                # presented as equivalent to a hit meeting the production floor.
                consensus_by_query[qid]['status'] = 'LOW_SIMILARITY'
            near_targets = {str(h.get('target') or '') for h in near_best}
            for rank, hit in enumerate(candidates, 1):
                candidate_tax, _tax_conf = _lookup_tax(hit['target'], taxa_map)
                support_fh.write(
                    f"{qid}\t{rank}\t{hit['target']}\t{hit['identity']:.2f}\t"
                    f"{hit.get('query_coverage') if hit.get('query_coverage') is not None else 'NA'}\t"
                    f"{hit.get('target_coverage') if hit.get('target_coverage') is not None else 'NA'}\t"
                    f"{hit.get('alignment_length') if hit.get('alignment_length') is not None else 'NA'}\t"
                    f"{candidate_tax or 'NA'}\t{'YES' if str(hit.get('target') or '') in near_targets else 'NO'}\t{hit.get('_evidence_tier', 'STANDARD_LADDER')}\n"
                )

    with open(conflicts_path, 'w', encoding='utf-8') as conflict_fh:
        conflict_fh.write('Query\tAssignmentStatus\tReportedTaxonomy\tAssignedRank\tConflictRank\tDominantTaxon\tDominantSupportPct\tRunnerUpTaxon\tRunnerUpSupportPct\tLeadPctPoints\tDistinctSupportingReferences\tCandidateCount\n')
        for qid, consensus in consensus_by_query.items():
            for conflict in consensus['conflicts']:
                conflict_fh.write(
                    f"{qid}\t{consensus['status']}\t{consensus['taxonomy'] or 'NA'}\t"
                    f"{consensus['rank'] or 'NA'}\t{conflict['rank']}\t{conflict['dominant_taxon']}\t"
                    f"{conflict['dominant_support_pct']:.2f}\t{conflict['runner_up_taxon']}\t"
                    f"{conflict['runner_up_support_pct']:.2f}\t{conflict['lead_pct_points']:.2f}\t"
                    f"{conflict['distinct_supporting_references']}\t{conflict['candidate_count']}\n"
                )
    logger.info('[CLASSIFY][%s] Supporting hits → %s; taxonomy conflicts/review → %s',
                db_name, support_path, conflicts_path)

    logger.info('[CLASSIFY][%s][DIAG] MATCH_PARSE_COMPLETE queries_with_hits=%d candidate_hits=%d', db_name, len(best_hits), sum(map(len, hits_by_query.values())))

    # A successful VSEARCH match is not a successful taxonomic assignment.
    # Fail explicitly if every matched reference lacks a taxonomy lookup; this
    # prevents malformed/unsupported header formats (as previously happened
    # with SILVA) from silently producing an all-NA taxonomy table.
    # Count a query as taxonomically supported when any retained candidate
    # maps to taxonomy; the top-scoring accession itself may lack a taxonomy
    # row while another near-best accession is correctly mapped.
    taxonomic_hits = sum(
        1 for candidates in hits_by_query.values()
        if any(_lookup_tax(hit.get('target', ''), taxa_map)[0] is not None for hit in candidates)
    )
    logger.info(
        '[CLASSIFY][%s][DIAG] TAXONOMY_LOOKUP_SUMMARY matched_hits=%d taxonomic_hits=%d missing_taxonomy=%d taxa_map_entries=%d',
        db_name, len(best_hits), taxonomic_hits,
        max(0, len(best_hits) - taxonomic_hits), len(taxa_map),
    )
    if best_hits and taxonomic_hits == 0:
        message = (
            f'{db_name}: VSEARCH returned {len(best_hits)} hit(s), but none could '
            'be mapped to taxonomy. Check the reference-header parser or taxa file; '
            'refusing to report this database as successfully classified.'
        )
        logger.error('[CLASSIFY][%s][DIAG] TAXONOMY_LOOKUP_FAILED %s', db_name, message)
        raise RuntimeError(message)
    if best_hits and taxonomic_hits < len(best_hits):
        logger.warning(
            '[CLASSIFY][%s] Taxonomy missing for %d/%d VSEARCH hits (%.1f%%).',
            db_name, len(best_hits) - taxonomic_hits, len(best_hits),
            100.0 * (len(best_hits) - taxonomic_hits) / len(best_hits),
        )

    diagnostic_path = os.path.join(outdir, f"classification_diagnostics_{safe_db_name}.tsv")
    with open(diagnostic_path, 'w', encoding='utf-8') as diagnostic_fh:
        diagnostic_fh.write(
            'Query\tSearchDiagnosticStatus\tAssignmentDiagnosticStatus\t'
            'BestHit\tBestIdentityPct\tBestQueryCoveragePct\tBestTargetCoveragePct\t'
            'AlignmentLength\tQueryLength\tTargetLength\tAssignedTaxonomy\tAssignedRank\t'
            'CandidateCount\tEvidenceTier\tAssignmentStatus\tTaxonomySupportPct\t'
            'AlignmentEvidenceStatus\tAssignmentEvidenceLevel\tEvidenceInterpretation\n'
        )
        for qid in all_query_ids:
            candidates = hits_by_query.get(qid, [])
            diagnostic_candidates = [h for h in candidates if str(h.get('_evidence_tier', '')).startswith('DIAGNOSTIC_')]
            standard_candidates = [h for h in candidates if not str(h.get('_evidence_tier', '')).startswith('DIAGNOSTIC_')]
            if standard_candidates:
                best = sorted(standard_candidates, key=_rank_hit, reverse=True)[0]
                reason, evidence = 'STANDARD_LADDER_HIT', 'STANDARD_LADDER'
            elif diagnostic_candidates:
                best = sorted(diagnostic_candidates, key=_rank_hit, reverse=True)[0]
                reason = _diagnostic_hit_reason(diagnostic_candidates)
                evidence = str(best.get('_evidence_tier') or 'DIAGNOSTIC_RESCUE')
            else:
                best, reason, evidence = None, 'NO_HIT_AT_DIAGNOSTIC_FLOOR', 'NONE'
            consensus = consensus_by_query.get(qid, {})
            assigned_taxonomy = consensus.get('taxonomy') or 'NA'
            assigned_rank = consensus.get('rank') or 'NA'
            if reason == 'NO_HIT_AT_DIAGNOSTIC_FLOOR':
                assignment_reason = 'NO_TAXONOMIC_ASSIGNMENT_NO_REFERENCE_HIT'
            elif reason != 'STANDARD_LADDER_HIT':
                if not best or _lookup_tax(best.get('target', ''), taxa_map)[0] is None:
                    assignment_reason = 'DIAGNOSTIC_HIT_TAXONOMY_UNMAPPED'
                else:
                    assignment_reason = 'LOW_SIMILARITY_TENTATIVE_ASSIGNMENT'
            elif consensus.get('status') == 'UNRESOLVED':
                mapped_taxonomies = [
                    _lookup_tax(hit.get('target', ''), taxa_map)[0]
                    for hit in candidates
                ]
                assignment_reason = (
                    'TAXONOMY_LOOKUP_FAILED'
                    if not any(value is not None for value in mapped_taxonomies)
                    else 'TAXONOMY_LINEAGE_UNRESOLVED'
                )
            elif consensus.get('status') == 'CONTESTED':
                assignment_reason = 'ADEQUATE_ALIGNMENT_CONSENSUS_CONTESTED'
            else:
                assignment_reason = 'ADEQUATE_ALIGNMENT_TAXONOMIC_ASSIGNMENT'
            diagnostic_fh.write(
                f"{qid}\t{reason}\t{assignment_reason}\t"
                f"{best.get('target') if best else 'NA'}\t"
                f"{best.get('identity', 'NA') if best else 'NA'}\t"
                f"{best.get('query_coverage', 'NA') if best else 'NA'}\t"
                f"{best.get('target_coverage', 'NA') if best else 'NA'}\t"
                f"{best.get('alignment_length', 'NA') if best else 'NA'}\t"
                f"{best.get('query_length', 'NA') if best else 'NA'}\t"
                f"{best.get('target_length', 'NA') if best else 'NA'}\t"
                f"{assigned_taxonomy}\t{assigned_rank}\t"
                f"{consensus.get('near_best_count', len(candidates))}\t{evidence}\t"
                f"{consensus.get('status', 'UNRESOLVED')}\t"
                f"{float(consensus.get('support') or 0.0) * 100.0:.2f}\t"
                f"{'NOT_ASSESSED' if not best else classification_evidence_interpretation(best.get('identity'), best.get('query_coverage'), best.get('target_coverage'), evidence_tier=evidence, assignment_status=consensus.get('status'), candidate_count=consensus.get('near_best_count'), alignment_length=best.get('alignment_length'), query_length=best.get('query_length'), target_length=best.get('target_length'))['alignment_evidence_status']}\t"
                f"{'UNRESOLVED' if not best else classification_evidence_interpretation(best.get('identity'), best.get('query_coverage'), best.get('target_coverage'), evidence_tier=evidence, assignment_status=consensus.get('status'), candidate_count=consensus.get('near_best_count'), alignment_length=best.get('alignment_length'), query_length=best.get('query_length'), target_length=best.get('target_length'))['assignment_evidence_level']}\t"
                f"{'insufficient alignment evidence' if not best else classification_evidence_interpretation(best.get('identity'), best.get('query_coverage'), best.get('target_coverage'), evidence_tier=evidence, assignment_status=consensus.get('status'), candidate_count=consensus.get('near_best_count'), alignment_length=best.get('alignment_length'), query_length=best.get('query_length'), target_length=best.get('target_length'))['evidence_interpretation']}\n"
            )
    logger.info('[CLASSIFY][%s] Per-query diagnostic summary → %s', db_name, diagnostic_path)

    results: dict = {}
    written = 0
    no_hit = 0
    with open(output, "w") as out_fh:
        out_fh.write('\t'.join(_classification_header(column_db_name or (db_name if db_name != 'main' else 'primary'))) + '\n')
        for qid in all_query_ids:
            if qid not in best_hits:
                out_fh.write(f"{qid}\tNA\tNA\tNA\tNA\t" + '\t'.join(_alignment_columns(None)) + "\tUNRESOLVED\t0.00\t0.00\tNA\tNONE\t0\tNOT_ASSESSED\tUNRESOLVED\tinsufficient alignment evidence\n")
                results[qid] = ('NA', 0.0, 'NA', None, None, None, None, None, None, None, None, 'UNRESOLVED', 0.0, 0.0, 'NA', 'NONE', 0, 'NOT_ASSESSED', 'UNRESOLVED', 'insufficient alignment evidence')
                no_hit += 1
                continue
            hit = best_hits[qid]
            sid, pct = hit['target'], hit['identity']
            _top_tax, conf = _lookup_tax(sid, taxa_map)
            consensus = consensus_by_query.get(qid, {'taxonomy': None, 'status': 'UNRESOLVED', 'support': 0.0, 'lead': 0.0, 'rank': ''})
            tax = consensus.get('taxonomy')
            # A missing source classifier confidence stays missing. Never
            # substitute percent identity: identity is not a probability.
            tax_str = tax if tax is not None else 'NA'
            status = consensus.get('status', 'UNRESOLVED')
            support_pct = float(consensus.get('support') or 0.0) * 100.0
            lead_pct = float(consensus.get('lead') or 0.0) * 100.0
            assigned_rank = consensus.get('rank') or 'NA'
            if status != 'CONFIDENT':
                logger.warning(
                    '[CLASSIFY][%s] Query %s assignment=%s taxonomy=%s support=%.1f%% lead=%.1fpp rank=%s; review=%s',
                    db_name, qid, status, tax_str, support_pct, lead_pct, assigned_rank, conflicts_path,
                )
            evidence_tier = str(hit.get('_evidence_tier') or 'STANDARD_LADDER')
            evidence_summary = classification_evidence_interpretation(
                pct, hit.get('query_coverage'), hit.get('target_coverage'),
                evidence_tier=evidence_tier, assignment_status=status,
                candidate_count=consensus.get('near_best_count', len(consensus_by_query.get(qid, {}).get('conflicts', []))),
                alignment_length=hit.get('alignment_length'),
                query_length=hit.get('query_length'),
                target_length=hit.get('target_length'),
            )
            candidate_count = int(consensus.get('near_best_count', 0) or 0)
            source_confidence = f'{float(conf):.4f}' if conf is not None else 'NA'
            interpretation = evidence_summary['evidence_interpretation'].replace('\t', ' ').replace('\n', ' ')
            out_fh.write(
                f"{qid}\t{sid}\t{pct:.2f}\t{tax_str}\t{source_confidence}\t"
                + '\t'.join(_alignment_columns(hit))
                + f"\t{status}\t{support_pct:.2f}\t{lead_pct:.2f}\t{assigned_rank}"
                + f"\t{evidence_tier}\t{candidate_count}\t{evidence_summary['alignment_evidence_status']}"
                + f"\t{evidence_summary['assignment_evidence_level']}\t{interpretation}\n"
            )
            results[qid] = (
                sid, pct, tax_str, conf, hit.get('query_coverage'), hit.get('target_coverage'),
                hit.get('alignment_length'), hit.get('query_length'), hit.get('target_length'),
                hit.get('mismatches'), hit.get('gaps'), status, support_pct, lead_pct, assigned_rank,
                evidence_tier, candidate_count, evidence_summary['alignment_evidence_status'],
                evidence_summary['assignment_evidence_level'], interpretation,
            )
            written += 1

    usage = resource.getrusage(resource.RUSAGE_SELF)
    logger.info('[CLASSIFY][%s][DIAG] OUTPUT_COMPLETE written=%d no_hit=%d results=%d elapsed=%.2fs maxrss=%s',
                db_name, written, no_hit, len(results), time.monotonic() - started, usage.ru_maxrss)
    logger.info(
        "[CLASSIFY][%s] %d/%d sequences classified, %d with no hit → %s",
        db_name, written, len(all_query_ids), no_hit, output,
    )
    return output, results


def run_all_classifications(
    input_fasta: str,
    outdir: str,
    primary_ref: str,
    primary_taxa: Optional[str] = None,
    primary_name: str = 'main',
    alt_refs: Optional[list] = None,
    threads: Optional[int] = None,
    main_db: Optional[str] = None,
) -> tuple[str, dict]:
    """Classify *input_fasta* against the primary and all alternative databases.

    Parameters
    ----------
    primary_ref   : Path to the primary reference FASTA.
    primary_taxa  : Optional taxa TSV/CSV for the primary reference.
    primary_name  : Display name for the primary database (default ``'main'``).
    alt_refs      : List of ``(ref_fasta, taxa_table_or_None, db_name)`` tuples.
    main_db       : Authoritative database name used to order the merged table.
                    Every database retains its explicit name in its output filename.

    Returns
    -------
    (primary_taxonomy_tsv, all_results)
    where ``all_results`` retains hit identity, source-provided classifier confidence (if any), alignment metrics, taxonomy consensus, and evidence interpretation.

    Side effects
    ------------
    * ``taxonomy_<name>.tsv``            – one file per database, including primary
    * ``matches_<name>.tsv``              – raw VSEARCH matches per database
    * ``taxonomy_all_dbs.tsv``            – wide merged table (multi-DB only)
    """
    # Older callers may omit primary_name. Resolve that to a meaningful
    # reference-derived label rather than reintroducing generic 'main' output.
    if primary_name == 'main':
        primary_name = _derive_db_name(primary_ref)
    effective_main = main_db or primary_name

    all_dbs = [(primary_ref, primary_taxa, primary_name)] + (alt_refs or [])
    all_results: dict = {}
    logger.info('[CLASSIFY][MULTI][DIAG] START pid=%d primary=%s authoritative=%s databases=%s input=%s threads=%s',
                os.getpid(), primary_name, effective_main, [n for _r, _t, n in all_dbs], input_fasta, threads)
    for index, (ref, taxa, name) in enumerate(all_dbs, start=1):
        logger.info('[CLASSIFY][MULTI][DIAG] DATABASE_BEGIN index=%d/%d name=%s ref=%s taxa=%s',
                    index, len(all_dbs), name, ref, taxa or 'FASTA_HEADERS')
        output_path, results = run_classification_single(
            input_fasta, outdir, ref, taxa, threads, name, column_db_name=name,
        )
        all_results[name] = results
        logger.info('[CLASSIFY][MULTI][DIAG] DATABASE_END index=%d/%d name=%s output=%s result_count=%d',
                    index, len(all_dbs), name, output_path, len(results))

    # Write wide merged taxonomy when more than one database was used
    if len(all_dbs) > 1:
        all_query_ids = [h for h, _ in read_fasta(input_fasta) if not is_ref_anchor(h)]
        _write_merged_taxonomy(outdir, all_results, effective_main, all_query_ids)

    safe_primary = re.sub(r'[^A-Za-z0-9]+', '_', str(effective_main)).strip('_') or 'primary'
    primary_tsv = os.path.join(outdir, f'taxonomy_{safe_primary}.tsv')
    if not os.path.isfile(primary_tsv):
        raise RuntimeError(f'Primary taxonomy output was not created for {effective_main}: {primary_tsv}')
    logger.info('[CLASSIFY][MULTI][DIAG] COMPLETE primary_output=%s databases=%d', primary_tsv, len(all_results))
    return primary_tsv, all_results


def _write_merged_taxonomy(
    outdir: str,
    all_results: dict,
    primary_name: str,
    all_query_ids: list,
) -> str:
    """Write ``taxonomy_all_dbs.tsv`` — a wide table with one set of columns per DB.

    Every database, including the primary reference, has a suffixed set of
    columns (``BestHit_<DB>``, ``Identity_<DB>``, ``Taxon_<DB>``,
    ``SourceClassifierConfidence_<DB>``) so assignments remain identifiable in merged tables.
    """
    import re as _re
    output = os.path.join(outdir, 'taxonomy_all_dbs.tsv')

    db_names = [primary_name] + [n for n in all_results if n != primary_name]

    with open(output, 'w') as fh:
        metric_names = ['QueryCoverage', 'TargetCoverage', 'AlignmentLength', 'QueryLength', 'TargetLength', 'Mismatches', 'Gaps']
        header = ['ID']
        for name in db_names:
            safe = _re.sub(r'[^A-Za-z0-9]+', '_', name).strip('_')
            header += [f'BestHit_{safe}', f'Identity_{safe}', f'Taxon_{safe}', f'SourceClassifierConfidence_{safe}']
            header += [f'{metric}_{safe}' for metric in metric_names]
            header += [f'AssignmentStatus_{safe}', f'TaxonomySupportPct_{safe}',
                       f'TaxonomyLeadPct_{safe}', f'AssignedRank_{safe}',
                       f'EvidenceTier_{safe}', f'TaxonomyCandidateCount_{safe}',
                       f'AlignmentEvidenceStatus_{safe}', f'AssignmentEvidenceLevel_{safe}',
                       f'EvidenceInterpretation_{safe}']
        fh.write('\t'.join(header) + '\n')

        primary_res = all_results.get(primary_name, {})
        for qid in all_query_ids:
            result = primary_res.get(qid, ('NA', 0.0, 'NA', 'NA', None, None, None, None, None, None, None))
            hit, pct, tax, conf = result[:4]
            row = [
                qid,
                str(hit),
                f"{pct:.2f}" if isinstance(pct, float) else str(pct),
                str(tax),
                f"{conf:.4f}" if isinstance(conf, float) else str(conf),
            ]
            row += ['NA' if value is None else (f'{value:.2f}' if isinstance(value, float) else str(value)) for value in result[4:11]]
            row += [str(result[11]) if len(result) > 11 else 'UNRESOLVED',
                    f'{float(result[12]):.2f}' if len(result) > 12 else '0.00',
                    f'{float(result[13]):.2f}' if len(result) > 13 else '0.00',
                    str(result[14]) if len(result) > 14 else 'NA',
                    str(result[15]) if len(result) > 15 else 'NA',
                    str(result[16]) if len(result) > 16 else 'NA',
                    str(result[17]) if len(result) > 17 else 'NA',
                    str(result[18]) if len(result) > 18 else 'NA',
                    str(result[19]) if len(result) > 19 else 'NA']
            for name in db_names[1:]:
                alt_result = all_results.get(name, {}).get(qid, ('NA', 0.0, 'NA', 'NA', None, None, None, None, None, None, None))
                a_hit, a_pct, a_tax, a_conf = alt_result[:4]
                row += [
                    str(a_hit),
                    f"{a_pct:.2f}" if isinstance(a_pct, float) else str(a_pct),
                    str(a_tax),
                    f"{a_conf:.4f}" if isinstance(a_conf, float) else str(a_conf),
                ]
                row += ['NA' if value is None else (f'{value:.2f}' if isinstance(value, float) else str(value)) for value in alt_result[4:11]]
                row += [str(alt_result[11]) if len(alt_result) > 11 else 'UNRESOLVED',
                        f'{float(alt_result[12]):.2f}' if len(alt_result) > 12 else '0.00',
                        f'{float(alt_result[13]):.2f}' if len(alt_result) > 13 else '0.00',
                        str(alt_result[14]) if len(alt_result) > 14 else 'NA',
                        str(alt_result[15]) if len(alt_result) > 15 else 'NA',
                        str(alt_result[16]) if len(alt_result) > 16 else 'NA',
                        str(alt_result[17]) if len(alt_result) > 17 else 'NA',
                        str(alt_result[18]) if len(alt_result) > 18 else 'NA',
                        str(alt_result[19]) if len(alt_result) > 19 else 'NA']
            fh.write('\t'.join(row) + '\n')

    logger.info("[CLASSIFY] Wrote merged multi-db taxonomy → %s", output)
    return output