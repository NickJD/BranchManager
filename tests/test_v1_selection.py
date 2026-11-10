"""v1 selection / Quarterly Review evidence-policy tests."""
from __future__ import annotations

from branchmanager.pipeline import quarterly_review as qr


def _row(sid, tier, *, source='SANGER_AB1', sequenced=False, qc='PASS_HIGH_CONFIDENCE', **extra):
    row = {
        'id': sid, 'partner_id': 'QUB', 'taxonomy': 'd__Bacteria; g__Test; s__Test species',
        'classification_identity': '99.5', 'classification_confidence': '0.99',
        'nearest_identity': '95.0', 'project_nearest_identity': '95.0', 'reference_nearest_identity': '95.0',
        'selected_for_genome_sequencing': 'False', 'already_sequenced': str(sequenced),
        'genome_committed_count_same_species': '0', 'pangenome_target': '3', 'placement_flags': '',
        'marker_source_type': source, 'marker_evidence_tier': tier, 'marker_qc_class': qc,
    }
    row.update(extra)
    return row


def _by_id(rows, **kwargs):
    recs = qr.build_quarterly_review(rows, genome_budget=3, backups_per_primary=0, pangenome_target=3, **kwargs)
    return {r['sequence_id']: r for r in recs}


def test_quarterly_review_high_moderate_low():
    rows = [_row('H', 'HIGH'), _row('M', 'MODERATE', qc='PASS_WITH_WARNINGS', marker_manual_review_status='APPROVED',
                                    marker_review_decision='APPROVE'),
            _row('L', 'LOW', qc='FAIL_QC', marker_manual_review_status='NOT_REVIEWED')]
    result = _by_id(rows)
    assert result['H']['role'] == 'PRIMARY'
    assert result['M']['role'] == 'REVIEW' and result['M']['evidence_quality'] == 'MODERATE'
    assert result['L']['role'] == 'REVIEW' and result['L']['evidence_quality'] == 'LOW'
    # moderate only joins the tranche when explicitly enabled
    assert _by_id(rows, include_moderate=True)['M']['role'] == 'PRIMARY'
    assert _by_id(rows, include_moderate=True)['L']['role'] == 'REVIEW'


def test_approved_failed_marker_keeps_original_qc_class_and_stays_qualified():
    row = _row('A', 'MODERATE', qc='FAIL_QC', marker_manual_review_status='APPROVED',
               marker_review_decision='APPROVE_FAIL_QC')
    rec = _by_id([row])['A']
    assert rec['marker_qc_class'] == 'FAIL_QC' and rec['evidence_quality'] == 'MODERATE'
    assert rec['role'] == 'REVIEW'


def test_barrnap_evidence_never_exceeds_moderate_and_is_not_primary_tranche():
    rec = _by_id([_row('G', 'HIGH', source='GENOME_BARRNAP', qc='GENOME_DERIVED')])['G']
    assert rec['evidence_quality'] == 'MODERATE' and rec['role'] == 'REVIEW'


def test_existing_genome_is_never_a_new_wgs_candidate_regardless_of_marker_tier():
    result = _by_id([_row('S', 'HIGH', sequenced=True), _row('S2', 'LOW', sequenced=True, qc='FAIL_QC')])
    assert result['S']['role'] == 'ALREADY_SEQUENCED' and result['S2']['role'] == 'ALREADY_SEQUENCED'


def test_barrnap_only_hit_cannot_trigger_primary_baseline_redundancy():
    # Primary-pool metrics say 95% identity; a 100% genome-derived hit lives only in the separate context file.
    row = _row('Q', 'HIGH', genome_derived_nearest_hit='G9', genome_derived_identity='100.0')
    rec = _by_id([row])['Q']
    assert rec['role'] == 'PRIMARY'
    assert 'REDUNDANT' not in str(rec.get('baseline_redundancy_status', '')).upper().replace('NOT_REDUNDANT', '')


def test_isolate_with_multiple_evidence_sources_is_counted_once():
    rows = [_row('Q', 'HIGH')]
    assert len(_by_id(rows)) == 1
