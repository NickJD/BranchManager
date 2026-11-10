"""Generic marker-sequence QC with an auditable manual-review override."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable

from branchmanager.utils.fasta import read_fasta, write_fasta

logger = logging.getLogger(__name__)


def _sequence_id(header: str) -> str:
    """Return the first FASTA token, excluding optional pipe-delimited metadata."""
    token = str(header or '').strip().split()[0] if str(header or '').strip() else ''
    return token.split('|', 1)[0]


def run_qc(
    input_fasta,
    outdir,
    min_len: int = 800,
    max_n_percent: float = 5.0,
    manual_approved_ids: Iterable[str] | None = None,
):
    """Filter sequences by length/N content, respecting explicit failed-QC approvals.

    ``manual_approved_ids`` contains sequence IDs with an explicit
    ``APPROVE_FAIL_QC`` decision in a validated review file. Such sequences are
    included in ``qc.fasta`` even when they fail this generic length/N filter.
    They are not reclassified as QC-passing: the manual override is recorded in
    ``qc_manual_overrides.tsv`` and, when relevant, ``qc_rejections.tsv``.

    The decision itself is validated elsewhere against the marker QC sidecar,
    reviewer identity, notes, and prepared candidate FASTA. This function only
    applies the already-recorded decision to avoid silently filtering the
    sequence a second time during downstream assessment.
    """
    max_n_percent = float(max_n_percent)
    min_len = int(min_len)
    if not 0.0 <= max_n_percent <= 100.0:
        raise ValueError('max_n_percent must be between 0 and 100')
    if min_len < 0:
        raise ValueError('min_len cannot be negative')

    output_dir = Path(outdir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / 'qc.fasta'
    override_ids = {str(value).strip() for value in (manual_approved_ids or []) if str(value).strip()}

    filtered = []
    rejected_too_short = []
    rejected_too_many_n = []
    rejected_details = []
    override_rows = []
    matched_override_ids = set()
    manual_override_applied = 0
    passed_generic_qc = 0
    total = 0

    for header, seq in read_fasta(str(input_fasta)):
        total += 1
        seq = str(seq or '').upper()
        sequence_id = _sequence_id(header)
        seq_len = len(seq)
        n_count = seq.count('N')
        n_percent = (100.0 * n_count / seq_len) if seq_len else 0.0
        reasons = []
        if seq_len < min_len:
            reasons.append('too_short')
        if n_percent > max_n_percent:
            reasons.append('too_many_n')

        manually_approved = sequence_id in override_ids or str(header).strip() in override_ids
        if manually_approved:
            matched_override_ids.add(sequence_id)
            override_rows.append({
                'sequence_id': sequence_id,
                'fasta_header': str(header),
                'length': seq_len,
                'n_count': n_count,
                'n_percent': f'{n_percent:.6f}',
                'automatic_qc_reasons': ';'.join(reasons),
                'manual_review_decision': 'APPROVE_FAIL_QC',
                'disposition': (
                    'INCLUDED_BY_MANUAL_REVIEW' if reasons
                    else 'APPROVED_MARKER_PASSED_GENERIC_PREFILTER'
                ),
            })

        if reasons:
            # These lists and their counters describe sequences actually
            # excluded by the generic filter. Manually approved records are
            # instead reported in qc_manual_overrides.tsv and marked as
            # INCLUDED_BY_MANUAL_REVIEW in qc_rejections.tsv.
            if not manually_approved:
                if 'too_short' in reasons:
                    rejected_too_short.append(header)
                if 'too_many_n' in reasons:
                    rejected_too_many_n.append(header)
            disposition = 'INCLUDED_BY_MANUAL_REVIEW' if manually_approved else 'EXCLUDED'
            rejected_details.append((
                header, seq_len, n_count, n_percent, ';'.join(reasons),
                'APPROVE_FAIL_QC' if manually_approved else '', disposition,
            ))
            if not manually_approved:
                continue
            manual_override_applied += 1
        else:
            passed_generic_qc += 1
        filtered.append((header, seq))

    write_fasta(filtered, str(output))

    rejection_path = output_dir / 'qc_rejections.tsv'
    with open(rejection_path, 'w', encoding='utf-8', newline='') as rf:
        rf.write(
            'ID\tLength\tNCount\tNPercent\tReasons\tMinLength\tMaxNPercent\t'
            'ManualReviewDecision\tDisposition\n'
        )
        for header, seq_len, n_count, n_percent, reasons, decision, disposition in rejected_details:
            rf.write(
                f'{header}\t{seq_len}\t{n_count}\t{n_percent:.6f}\t{reasons}\t'
                f'{min_len}\t{max_n_percent:g}\t{decision}\t{disposition}\n'
            )

    override_path = output_dir / 'qc_manual_overrides.tsv'
    with open(override_path, 'w', encoding='utf-8', newline='') as handle:
        handle.write(
            'sequence_id\tfasta_header\tlength\tn_count\tn_percent\t'
            'automatic_qc_reasons\tmanual_review_decision\tdisposition\n'
        )
        for row in override_rows:
            handle.write(
                f"{row['sequence_id']}\t{row['fasta_header']}\t{row['length']}\t"
                f"{row['n_count']}\t{row['n_percent']}\t{row['automatic_qc_reasons']}\t"
                f"{row['manual_review_decision']}\t{row['disposition']}\n"
            )
        for missing_id in sorted(override_ids - matched_override_ids):
            handle.write(
                f'{missing_id}\t\t\t\t\t\tAPPROVE_FAIL_QC\tNOT_FOUND_IN_INPUT\n'
            )

    stats_path = output_dir / 'qc.stats'
    with open(stats_path, 'w', encoding='utf-8') as sf:
        sf.write(f'total_input\t{total}\n')
        sf.write(f'kept\t{len(filtered)}\n')
        sf.write(f'passed_generic_qc\t{passed_generic_qc}\n')
        sf.write(f'retained_by_manual_override\t{manual_override_applied}\n')
        sf.write(f'manual_approved_records_present\t{len(matched_override_ids)}\n')
        sf.write(f'manual_override_requested\t{len(override_ids)}\n')
        sf.write(f'rejected_total\t{sum(1 for row in rejected_details if row[-1] == "EXCLUDED")}\n')
        sf.write(f'rejected_too_short\t{len(rejected_too_short)}\n')
        sf.write(f'rejected_too_many_n\t{len(rejected_too_many_n)}\n')
        sf.write(f'min_len\t{min_len}\n')
        sf.write(f'max_n_percent\t{max_n_percent:g}\n')
        sf.write(f'rejection_details\t{rejection_path}\n')
        sf.write(f'manual_override_details\t{override_path}\n')
        sf.write('sample_ids_kept:\n')
        for header, _ in filtered[:10]:
            sf.write(f'\t{header}\n')
        if rejected_too_short:
            sf.write('sample_ids_too_short:\n')
            for header in rejected_too_short[:5]:
                sf.write(f'\t{header}\n')
        if rejected_too_many_n:
            sf.write('sample_ids_too_many_n:\n')
            for header in rejected_too_many_n[:5]:
                sf.write(f'\t{header}\n')

    rejected_total = sum(1 for row in rejected_details if row[-1] == 'EXCLUDED')
    logger.info(
        '[QC] %d input → %d retained (%d by automated QC; %d by explicit manual override), '
        '%d excluded; %d too short (<%dbp), %d high-N (>%.2f%%). Stats: %s; details: %s',
        total, len(filtered),
        passed_generic_qc, manual_override_applied,
        rejected_total, len(rejected_too_short), min_len,
        len(rejected_too_many_n), max_n_percent, stats_path, rejection_path,
    )
    if rejected_too_short:
        logger.debug('[QC] Too-short input records: %s', ', '.join(rejected_too_short[:10]))
    if rejected_too_many_n:
        logger.debug('[QC] High-N input records: %s', ', '.join(rejected_too_many_n[:10]))
    return str(output)