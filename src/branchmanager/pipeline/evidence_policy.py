"""Single evidence-quality policy shared by selection, selection sets and Quarterly Review.

HIGH      eligible for normal assessment and ranking.
MODERATE  eligible for assessment but kept as a qualified/review-level candidate.
LOW       excluded from automated candidate selection until the evidence is resolved.

Genome availability is deliberately NOT part of this policy.
"""
from __future__ import annotations

from typing import Optional

BLOCKING_FLAGS = frozenset({
    'NO_REFERENCE_HIT', 'NO_CLASSIFICATION',
    'CHIMERA', 'CHIMERA_CONFIRMED', 'CHIMERA_DETECTED',
    'CHIMERA_INDETERMINATE', 'CHIMERA_CHECK_SKIPPED', 'CHIMERA_NOT_RUN',
    'VERY_SHORT', 'HIGH_N_CONTENT', 'MARKER_QC_FAILED',
    'MARKER_QC_REVIEW_REQUIRED', 'MARKER_QC_UNVERIFIED',
    'MARKER_REPLACEMENT_QC_FAILED', 'MANUAL_APPROVAL_SEVERE_QC_OVERRIDE',
})
_ORDER = {'LOW': 0, 'MODERATE': 1, 'HIGH': 2}
_APPROVALS = {'APPROVE', 'APPROVE_FAIL_QC', 'APPROVE_UNVERIFIED'}


def _num(value) -> Optional[float]:
    try:
        return None if value in (None, '', 'NA', 'None') else float(value)
    except (TypeError, ValueError):
        return None


def _provenance_tier(row: dict) -> Optional[str]:
    """Tier implied by explicit marker provenance, or None when the row carries none."""
    source = str(row.get('marker_source_type') or '').strip().upper()
    tier = str(row.get('marker_evidence_tier') or '').strip().upper()
    qc_class = str(row.get('marker_qc_class') or '').strip().upper()
    review = str(row.get('marker_manual_review_status') or '').strip().upper()
    decision = str(row.get('marker_review_decision') or '').strip().upper()
    if tier not in _ORDER:
        return None
    if tier == 'LOW' or review == 'REJECTED':
        return 'LOW'
    if qc_class in {'FAIL_QC', 'QUALITY_UNVERIFIED'}:
        if review == 'APPROVED' and decision in _APPROVALS:
            return 'MODERATE'
        return 'LOW'
    if source == 'GENOME_BARRNAP' and tier == 'HIGH':
        return 'MODERATE'
    return tier


def _heuristic_tier(row: dict, *, require_taxonomy: bool) -> str:
    flags = {f.strip().upper() for f in str(row.get('placement_flags') or '').split(';') if f.strip()}
    identity = _num(row.get('classification_identity'))
    coverage = _num(row.get('classification_query_coverage'))
    taxonomy = str(row.get('taxonomy') or '').strip().lower()
    marker_flag = str(row.get('marker_qc_flag') or '').upper()
    if marker_flag in {'MARKER_QC_FAILED', 'MARKER_QC_REVIEW_REQUIRED', 'MARKER_QC_UNVERIFIED'}:
        return 'LOW'
    if flags & BLOCKING_FLAGS or (identity is not None and identity < 90.0):
        return 'LOW'
    if coverage is not None and coverage < 80.0:
        return 'LOW'
    if require_taxonomy and identity is None and taxonomy in ('', 'na', 'none'):
        return 'LOW'
    if (
        flags & {'LOW_CLASSIFICATION', 'DISAGREEMENT', 'CONFLICT'}
        or (identity is not None and identity < 95.0)
        or (coverage is not None and coverage < 90.0)
    ):
        return 'MODERATE'
    return 'HIGH'


def evidence_quality(row: dict) -> str:
    """Return HIGH/MODERATE/LOW. Provenance sets a ceiling; classification checks can only lower it."""
    explicit = _provenance_tier(row)
    if explicit is not None:
        return min(explicit, _heuristic_tier(row, require_taxonomy=False), key=_ORDER.get)
    declared = str(row.get('evidence_quality') or '').strip().upper()
    if declared in _ORDER:
        return declared
    return _heuristic_tier(row, require_taxonomy=True)



def classification_evidence_interpretation(
    identity, query_coverage, target_coverage, evidence_tier='STANDARD_LADDER',
    assignment_status='PROVISIONAL', candidate_count=None,
    alignment_length=None, query_length=None, target_length=None,
):
    """Describe alignment and taxonomy evidence without treating identity as probability.

    Thresholds are reporting bands, not biological species cut-offs:
    HIGH alignment: >=98% identity, >=95% query coverage, >=90% target coverage.
    MODERATE alignment: >=95% identity, >=90% query coverage, >=80% target coverage.
    Lower values are LOW. Diagnostic-only evidence is always TENTATIVE.
    """
    def number(value):
        try:
            if value in (None, '', 'NA', 'None'):
                return None
            return float(value)
        except (TypeError, ValueError):
            return None

    ident, qcov, tcov = number(identity), number(query_coverage), number(target_coverage)
    tier = str(evidence_tier or '').upper()
    status = str(assignment_status or '').upper()
    count = number(candidate_count)

    if ident is None or qcov is None:
        alignment_status = 'NOT_ASSESSED'
    elif qcov < 90.0 and tcov is not None and tcov < 80.0:
        alignment_status = 'BOTH_PARTIAL'
    elif qcov < 90.0:
        alignment_status = 'QUERY_PARTIAL'
    elif tcov is not None and tcov < 80.0:
        alignment_status = 'TARGET_PARTIAL'
    elif qcov < 95.0 or (tcov is not None and tcov < 90.0):
        alignment_status = 'PARTIAL'
    else:
        alignment_status = 'NEAR_FULL_LENGTH'

    if tier.startswith('DIAGNOSTIC'):
        level = 'TENTATIVE'
    elif ident is None or qcov is None:
        level = 'UNRESOLVED'
    elif ident < 90.0 or qcov < 80.0 or (tcov is not None and tcov < 70.0):
        level = 'LOW'
    elif ident < 95.0 or qcov < 90.0 or (tcov is not None and tcov < 80.0):
        level = 'LOW'
    elif ident < 98.0 or qcov < 95.0 or (tcov is not None and tcov < 90.0):
        level = 'MODERATE'
    elif status in {'CONTESTED', 'UNRESOLVED'}:
        level = 'LOW'
    elif status != 'CONFIDENT' or (count is not None and count < 3):
        level = 'MODERATE'
    else:
        level = 'HIGH'

    reasons = []
    if ident is not None:
        reasons.append(f'{ident:.1f}% identity across the aligned region')
    if qcov is not None:
        reasons.append(f'{qcov:.1f}% query coverage')
    if tcov is not None:
        reasons.append(f'{tcov:.1f}% reference coverage')
    try:
        aln = int(float(alignment_length)) if alignment_length not in (None, '', 'NA') else None
        qlen = int(float(query_length)) if query_length not in (None, '', 'NA') else None
        tlen = int(float(target_length)) if target_length not in (None, '', 'NA') else None
    except (TypeError, ValueError):
        aln, qlen, tlen = None, None, None
    if aln is not None:
        lengths = f'{aln} bp aligned'
        if qlen is not None:
            lengths += f' of {qlen} bp query'
        if tlen is not None:
            lengths += f' and {tlen} bp reference'
        reasons.append(lengths)
    if tier.startswith('DIAGNOSTIC'):
        reasons.append('diagnostic-only search evidence; not production-supported')
    if status == 'CONTESTED':
        reasons.append('supporting references disagree taxonomically')
    elif status == 'PROVISIONAL':
        reasons.append('taxonomic consensus is provisional')
    if count is not None:
        reasons.append(f'{int(count)} near-best candidate reference(s)')
        if count < 2:
            reasons.append('single-reference support is not independent taxonomic confirmation')
    interpretation = '; '.join(reasons) if reasons else 'insufficient alignment evidence'
    return {
        'alignment_evidence_status': alignment_status,
        'assignment_evidence_level': level,
        'evidence_interpretation': interpretation,
    }
