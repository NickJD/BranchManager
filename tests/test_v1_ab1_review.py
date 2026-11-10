"""v1 AB1 quality handling and manual-review workflow tests (types A and C)."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
from test_taxonomy_and_db import _write_minimal_ab1  # noqa: E402

from branchmanager.db.interface import Database
from branchmanager.marker_evidence import (
    apply_marker_review_decisions, validate_review_decision, upsert_active_marker_evidence,
)
from branchmanager.pipeline import paper_trail as pt
from branchmanager.utils.fasta import reverse_complement

FWD = 'AAAACCCCGGGGTTTTAAAA'
REV = 'CCCCGGGGTTTTAAAACCCC'


def _run(tmp_path, quals_f, quals_r, **kwargs):
    reads = tmp_path / 'reads'
    reads.mkdir(exist_ok=True)
    _write_minimal_ab1(reads / 'a.ab1', FWD, quals_f)
    _write_minimal_ab1(reads / 'b.ab1', reverse_complement(REV), quals_r)
    sample_map = tmp_path / 'map.tsv'
    sample_map.write_text('isolate_id\t27F\t907R\nIsoX\ta.ab1\tb.ab1\n')
    return pt.run_paper_trail(
        [str(reads)], tmp_path / 'out', sample_map=str(sample_map), min_length=8,
        min_overlap=8, min_overlap_identity=0.9, **kwargs,
    )


def _rows(path):
    with open(path, newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def test_high_confidence_sanger_passes_with_verified_quality(tmp_path):
    out = _run(tmp_path, [40] * 20, [40] * 20)
    row = _rows(out['assembly_tsv'])[0]
    # the synthetic AB1 has no trace channels, which is the only (warning) reason
    assert row['QCClass'] in {'PASS_HIGH_CONFIDENCE', 'PASS_WITH_WARNINGS'}
    assert row['QualityStatus'] == 'PHRED_VERIFIED'
    assert 'quality' not in row['Reasons']


@pytest.mark.parametrize('label,quals,expected', [
    ('missing', [], 'MISSING'),
    ('zero', [0] * 20, 'ZERO_FILLED'),
])
def test_unusable_pcon_is_never_high_confidence(tmp_path, label, quals, expected):
    out = _run(tmp_path, quals, quals)
    row = _rows(out['assembly_tsv'])[0]
    assert row['QCClass'] != 'PASS_HIGH_CONFIDENCE'
    assert expected in row['QualityStatus']
    read_qc = Path(out['read_qc_tsv']).read_text()
    assert 'ab1_zero_filled_pcon' in read_qc or 'ab1_missing_pcon' in read_qc
    if label == 'zero':
        assert 'zero_filled_pcon' in row['Reasons'] or 'zero_filled' in read_qc


def test_allow_missing_quality_keeps_zero_filled_reads_unverified(tmp_path):
    out = _run(tmp_path, [0] * 20, [0] * 20, allow_missing_quality=True)
    row = _rows(out['assembly_tsv'])[0]
    assert row['QCClass'] != 'PASS_HIGH_CONFIDENCE'
    assert 'output_quality_unverified' in row['Reasons']
    assert 'ZERO_FILLED' in row['QualityStatus']
    template = _rows(Path(out['assembly_tsv']).parent / 'marker_review_template.tsv')
    assert template and template[0]['quality_status'].startswith('PHRED_UNVERIFIED')


def test_decision_validation_rules():
    ok = dict(reviewer='NJ', notes='trace inspected', has_consensus=True)
    assert validate_review_decision('APPROVE_FAIL_QC', 'FAIL_QC', **ok) == ''
    assert validate_review_decision('APPROVE_FAIL_QC', 'FAIL_QC', reviewer='NJ', notes='x', has_consensus=False)
    assert validate_review_decision('APPROVE_FAIL_QC', 'PASS_HIGH_CONFIDENCE', **ok)
    assert validate_review_decision('APPROVE_FAIL_QC', 'FAIL_QC', reviewer='', notes='x', has_consensus=True)
    assert validate_review_decision('APPROVE_FAIL_QC', 'FAIL_QC', reviewer='NJ', notes='', has_consensus=True)
    assert validate_review_decision('REJECT', 'FAIL_QC', reviewer='NJ', notes='bad', has_consensus=False) == ''


def _db_with_failed(tmp_path, seq='ACGT' * 300):
    db = Database(str(tmp_path / 'p.sqlite'))
    db.initialise()
    db.upsert_dataset_role('P', 'candidate', genomes_available=False)
    db.insert_sequences([('F1', seq)], dataset='P')
    prov = [{'sequence_id': 'F1', 'marker_qc_class': 'FAIL_QC', 'marker_qc_recommendation': 'RESEQUENCE',
             'marker_qc_reasons': 'low_quality', 'manual_review_status': 'NOT_REVIEWED'}]
    db.upsert_sequence_provenance(prov)
    upsert_active_marker_evidence(db, [('F1', seq)], prov, source_type='SANGER_AB1', dataset='P')
    return db


def _evidence(db, sid='F1'):
    with db.connect() as conn:
        return conn.execute(
            'SELECT evidence_tier,qc_class,manual_review_status,review_decision,reviewer,active_for_selection,eligible_primary '
            'FROM marker_evidence WHERE sequence_id=?', (sid,)).fetchone()


def test_unreviewed_failed_marker_is_excluded(tmp_path):
    db = _db_with_failed(tmp_path)
    tier, qc, _status, _dec, _rev, active, primary = _evidence(db)
    assert (tier, qc) == ('LOW', 'FAIL_QC') and not (active and primary)


def test_approved_failed_marker_is_moderate_and_keeps_fail_qc(tmp_path):
    db = _db_with_failed(tmp_path)
    tsv = tmp_path / 'd.tsv'
    tsv.write_text('sequence_id\tdecision\treviewer\tnotes\nF1\tAPPROVE_FAIL_QC\tNJ\tchromatogram clean\n')
    results = apply_marker_review_decisions(db, tsv)
    assert results[0]['result'] == 'APPLIED'
    tier, qc, status, decision, reviewer, active, primary = _evidence(db)
    assert (tier, qc, status, decision, reviewer) == ('MODERATE', 'FAIL_QC', 'APPROVED', 'APPROVE_FAIL_QC', 'NJ')
    assert active and primary
    with db.connect() as conn:
        assert conn.execute('SELECT COUNT(*) FROM marker_review_decisions').fetchone()[0] == 1
    apply_marker_review_decisions(db, tsv)  # idempotent history
    with db.connect() as conn:
        assert conn.execute('SELECT COUNT(*) FROM marker_review_decisions').fetchone()[0] == 1


def test_rejected_failed_marker_stays_excluded(tmp_path):
    db = _db_with_failed(tmp_path)
    tsv = tmp_path / 'd.tsv'
    tsv.write_text('sequence_id\tdecision\treviewer\tnotes\nF1\tREJECT\tNJ\tunusable\n')
    apply_marker_review_decisions(db, tsv)
    tier, _qc, status, _d, _r, active, primary = _evidence(db)
    assert (tier, status) == ('LOW', 'REJECTED') and not (active and primary)


def test_failed_marker_without_consensus_cannot_be_rescued(tmp_path):
    db = Database(str(tmp_path / 'p.sqlite'))
    db.initialise()
    db.insert_sequences([('N1', '')], dataset='P')
    tsv = tmp_path / 'd.tsv'
    tsv.write_text('sequence_id\tdecision\treviewer\tnotes\nN1\tAPPROVE_FAIL_QC\tNJ\tplease\n')
    results = apply_marker_review_decisions(db, tsv)
    assert results[0]['result'] == 'ERROR'
