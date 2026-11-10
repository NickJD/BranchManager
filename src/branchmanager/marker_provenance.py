"""Bridge marker preparation QC into downstream candidate assessment."""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Optional

from branchmanager.marker_evidence import (
    APPROVAL_DECISIONS, REVIEW_DECISIONS, validate_review_decision,
)
from branchmanager.run_manifest import sha256_file
from branchmanager.utils.io import open_text
from branchmanager.utils.text import normalise_header


REJECTED_REVIEW_VALUES = {'rejected', 'reject', 'failed', 'fail', 'no'}


def _read_table(path: str | Path) -> list[dict[str, str]]:
    """Read a small CSV/TSV sidecar without relying on heuristic sniffing.

    File extensions are authoritative for the formats BranchManager writes.
    Sniffing is retained only for files with an unfamiliar extension. In
    particular, csv.Sniffer raises on valid but short/header-only samples;
    that must never turn a valid TSV into a fatal exception.
    """
    path_text = str(path).lower()
    if path_text.endswith(('.tsv', '.tsv.gz', '.txt', '.txt.gz')):
        delimiter = '\t'
    elif path_text.endswith(('.csv', '.csv.gz')):
        delimiter = ','
    else:
        delimiter = None

    with open_text(path) as handle:
        sample = handle.read(4096)
        handle.seek(0)
        if delimiter is None:
            try:
                delimiter = csv.Sniffer().sniff(sample, delimiters=',\t').delimiter
            except csv.Error:
                # Deterministic fallback for extensionless/unknown sidecars.
                # Prefer TSV on ties because most BranchManager sidecars are TSV.
                delimiter = '\t' if sample.count('\t') >= sample.count(',') else ','
                logging.getLogger(__name__).debug(
                    'Could not sniff CSV delimiter for %r; using %r', path, delimiter,
                )
        return [
            {normalise_header(key): str(value or '').strip() for key, value in row.items() if key is not None}
            for row in csv.DictReader(handle, delimiter=delimiter)
        ]


def discover_marker_qc(input_fasta: str | Path, explicit: Optional[str | Path] = None) -> Optional[Path]:
    if explicit:
        candidate = Path(explicit).expanduser()
        if not candidate.is_file():
            raise ValueError(f'Marker-QC sidecar not found: {candidate}')
        return candidate.resolve()
    parent = Path(input_fasta).expanduser().resolve().parent
    for name in ('assembly_report.tsv', 'marker_qc.tsv'):
        candidate = parent / name
        if candidate.is_file():
            return candidate
    return None


def _review_decisions(path: Optional[str | Path]) -> dict[str, dict[str, str]]:
    """Read explicit v1 review decisions; safe legacy rejection aliases remain accepted."""
    if not path:
        return {}
    rows = _read_table(path)
    if not rows:
        raise ValueError(
            f'Marker review decisions file is empty or contains no data rows: {path}'
        )
    headers = set(rows[0])
    if not headers.intersection({'sequence_id', 'isolate_id', 'id'}):
        raise ValueError(
            f'Marker review decisions file has no sequence identifier column '
            f'(expected sequence_id, isolate_id or id): {path}'
        )
    if not headers.intersection({'decision', 'review_status', 'manual_review_status'}):
        raise ValueError(
            f'Marker review decisions file has no decision column '
            f'(expected decision, review_status or manual_review_status): {path}'
        )
    decisions = {}
    for row in rows:
        sequence_id = row.get('sequence_id') or row.get('isolate_id') or row.get('id') or ''
        raw = row.get('decision') or row.get('review_status') or row.get('manual_review_status') or ''
        key = raw.strip().lower()
        explicit = raw.strip().upper() if raw.strip().upper() in REVIEW_DECISIONS else ''
        if explicit in APPROVAL_DECISIONS:
            status = 'APPROVED'
        elif explicit in {'REJECT', 'RESEQUENCE'}:
            status = 'REJECTED'
        elif key in REJECTED_REVIEW_VALUES:
            explicit = 'RESEQUENCE' if key in {'failed', 'fail'} else 'REJECT'
            status = 'REJECTED'
        elif key in {'approved', 'approve', 'accepted', 'accept', 'pass', 'yes'}:
            raise ValueError(
                f"Manual approval for {sequence_id or '<unknown>'} must use an explicit v1 decision "
                "and include reviewer and explanatory notes"
            )
        else:
            status = 'NOT_REVIEWED'
        if sequence_id:
            decisions[sequence_id] = {
                'status': status,
                'decision': explicit,
                'reviewer': row.get('reviewer', ''),
                'notes': row.get('notes') or row.get('detail') or '',
                'reviewed_at': row.get('reviewed_at', ''),
                'review_type': row.get('review_type', ''),
            }
    return decisions


def _has_consensus(fasta_sequences: dict, *ids: str) -> bool:
    return any(fasta_sequences.get(i) for i in ids if i)


def load_marker_provenance(
    input_fasta: str | Path,
    *,
    qc_path: Optional[str | Path] = None,
    review_path: Optional[str | Path] = None,
    id_map: Optional[dict[str, str]] = None,
    accept_unverified: bool = False,
    unverified_override_authorizer: str = '',
    unverified_override_reason: str = '',
) -> tuple[list[dict], Optional[str]]:
    """Return per-sequence provenance and the resolved QC sidecar path."""
    if accept_unverified and not (
        str(unverified_override_authorizer or '').strip()
        and str(unverified_override_reason or '').strip()
    ):
        raise ValueError(
            "--accept-unverified-marker-qc requires an override authorizer and reason; "
            "it permits processing but does not approve the marker for selection"
        )
    resolved = discover_marker_qc(input_fasta, qc_path)
    decisions = _review_decisions(review_path)
    mapping = id_map or {}
    input_hash = sha256_file(input_fasta)
    manifest = Path(input_fasta).resolve().parent / 'run_manifest.json'
    read_qc = Path(input_fasta).resolve().parent / 'read_qc.tsv'
    read_files_by_sequence: dict[str, list[str]] = {}
    read_checksums_by_sequence: dict[str, dict[str, str]] = {}
    if read_qc.is_file():
        for row in _read_table(read_qc):
            sequence_id = row.get('sequenceid') or row.get('sequence_id') or ''
            source = row.get('sourcefile') or row.get('source_file') or ''
            if sequence_id and source:
                read_files_by_sequence.setdefault(sequence_id, []).append(source)
                source_path = Path(source).expanduser()
                if source_path.is_file():
                    try:
                        read_checksums_by_sequence.setdefault(sequence_id, {})[source] = sha256_file(source_path)
                    except OSError:
                        pass

    from branchmanager.utils.fasta import read_fasta as _read_fasta
    fasta_sequences = {str(h).split()[0]: str(q) for h, q in _read_fasta(str(input_fasta))}
    records = []
    if resolved:
        for row in _read_table(resolved):
            source_id = row.get('sequenceid') or row.get('sequence_id') or row.get('id') or ''
            sequence_id = mapping.get(source_id, source_id)
            if not sequence_id:
                continue
            qc_class = row.get('qcclass') or row.get('qc_class') or 'QUALITY_UNVERIFIED'
            recommendation = row.get('recommendation') or 'MANUAL_REVIEW'
            reasons = row.get('reasons') or ''
            review = decisions.get(source_id) or decisions.get(sequence_id) or {}
            if review.get('decision'):
                error = validate_review_decision(
                    review['decision'], qc_class, reviewer=review.get('reviewer', ''),
                    notes=review.get('notes', ''),
                    has_consensus=_has_consensus(fasta_sequences, source_id, sequence_id),
                )
                if error:
                    raise ValueError(f'Manual review decision for {sequence_id}: {error}')
            records.append({
                'sequence_id': sequence_id,
                'source_manifest': str(manifest) if manifest.is_file() else '',
                'source_sequence_file': str(Path(input_fasta).resolve()),
                'source_sequence_sha256': input_hash,
                'source_read_ids': row.get('readids') or row.get('read_ids') or '',
                'source_used_read_ids': row.get('usedreadids') or row.get('used_read_ids') or '',
                'source_read_files': ';'.join(sorted(set(read_files_by_sequence.get(source_id, [])))),
                'source_read_checksums_json': json.dumps(
                    read_checksums_by_sequence.get(source_id, {}), sort_keys=True,
                ),
                'marker_qc_class': qc_class,
                'marker_qc_recommendation': recommendation,
                'marker_qc_reasons': reasons,
                'manual_review_status': review.get('status', 'NOT_REVIEWED'),
                'review_decision': review.get('decision', ''),
                'reviewer': review.get('reviewer', ''),
                'review_notes': review.get('notes', ''),
                'reviewed_at': review.get('reviewed_at', ''),
                'review_type': review.get('review_type', ''),
                'marker_source_type': str(row.get('marker_source_type') or row.get('markersourcetype') or '').upper(),
            })
    else:
        from branchmanager.utils.fasta import read_fasta
        for source_id, _sequence in read_fasta(input_fasta):
            sequence_id = mapping.get(source_id, source_id)
            review = decisions.get(source_id) or decisions.get(sequence_id) or {}
            if review.get('decision'):
                error = validate_review_decision(
                    review['decision'], 'QUALITY_UNVERIFIED', reviewer=review.get('reviewer', ''),
                    notes=review.get('notes', ''), has_consensus=bool(_sequence),
                )
                if error:
                    raise ValueError(f'Manual review decision for {sequence_id}: {error}')
            records.append({
                'sequence_id': sequence_id,
                'source_manifest': '',
                'source_sequence_file': str(Path(input_fasta).resolve()),
                'source_sequence_sha256': input_hash,
                'source_read_ids': '',
                'source_used_read_ids': '',
                'source_read_files': '',
                'source_read_checksums_json': '{}',
                'marker_qc_class': 'QUALITY_UNVERIFIED',
                'marker_qc_recommendation': 'MANUAL_REVIEW',
                'marker_qc_reasons': (
                    'no_marker_qc_sidecar; '
                    f'batch_processing_override_authorized_by={unverified_override_authorizer.strip()}; '
                    f'batch_processing_override_reason={unverified_override_reason.strip()}'
                    if accept_unverified else 'no_marker_qc_sidecar'
                ),
                'manual_review_status': review.get('status', 'NOT_REVIEWED'),
                'review_decision': review.get('decision', ''),
                'reviewer': review.get('reviewer', ''),
                'review_notes': review.get('notes', ''),
                'reviewed_at': review.get('reviewed_at', ''),
            })
    return records, str(resolved) if resolved else None


def marker_qc_flag(record: Optional[dict]) -> str:
    if not record:
        return 'MARKER_QC_UNVERIFIED'
    qc_class = str(record.get('marker_qc_class') or '').upper()
    review = str(record.get('manual_review_status') or '').upper()
    decision = str(record.get('review_decision') or '').upper()
    if review == 'REJECTED' or (qc_class == 'FAIL_QC' and not (review == 'APPROVED' and decision == 'APPROVE_FAIL_QC')):
        return 'MARKER_QC_FAILED'
    if qc_class == 'PASS_HIGH_CONFIDENCE':
        return ''
    if review == 'APPROVED':
        return 'MARKER_QC_REVIEW_APPROVED'
    return 'MARKER_QC_REVIEW_REQUIRED'


def write_marker_qc_bridge(path: str | Path, records: list[dict]) -> str:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        'sequence_id', 'marker_qc_class', 'marker_qc_recommendation', 'manual_review_status',
        'review_decision', 'reviewer', 'review_notes', 'marker_qc_reasons', 'source_sequence_file', 'source_sequence_sha256',
        'source_read_ids', 'source_used_read_ids', 'source_read_files',
        'source_read_checksums_json', 'source_manifest', 'chimera_call', 'chimera_score',
    ]
    with open(output, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(records)
    return str(output)
