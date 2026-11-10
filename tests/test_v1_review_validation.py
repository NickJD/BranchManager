"""Manual-review aliases and batch quality override audit tests."""
from __future__ import annotations

import csv

import pytest

from branchmanager.marker_evidence import upsert_active_marker_evidence
from branchmanager.marker_provenance import _review_decisions, load_marker_provenance
from branchmanager.db.interface import Database


def _tsv(path, rows):
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    return path


def test_legacy_rejection_alias_is_mapped_but_approval_alias_is_rejected(tmp_path):
    reject = _tsv(tmp_path / 'reject.tsv', [{
        'sequence_id': 'I1', 'decision': 'reject', 'reviewer': 'R', 'notes': 'not supported',
    }])
    assert _review_decisions(reject)['I1']['decision'] == 'REJECT'
    approval = _tsv(tmp_path / 'approve.tsv', [{
        'sequence_id': 'I1', 'decision': 'approved', 'reviewer': 'R', 'notes': 'looked good',
    }])
    with pytest.raises(ValueError, match='explicit v1 decision'):
        _review_decisions(approval)


def test_explicit_approval_alias_path_requires_reviewer_notes_and_consensus(tmp_path):
    fasta = tmp_path / 'markers.fa'
    fasta.write_text('>I1\nACGTACGT\n')
    review = _tsv(tmp_path / 'review.tsv', [{
        'sequence_id': 'I1', 'decision': 'APPROVE_UNVERIFIED', 'reviewer': 'R', 'notes': '',
    }])
    with pytest.raises(ValueError, match='explanatory notes'):
        load_marker_provenance(fasta, review_path=review)


def test_accept_unverified_is_audited_batch_override_not_individual_approval(tmp_path):
    fasta = tmp_path / 'markers.fa'
    fasta.write_text('>I1\nACGTACGT\n')
    with pytest.raises(ValueError, match='authorizer and reason'):
        load_marker_provenance(fasta, accept_unverified=True)

    provenance, qc_path = load_marker_provenance(
        fasta, accept_unverified=True, unverified_override_authorizer='release-manager',
        unverified_override_reason='partner supplied sequence independently verified against source file',
    )
    row = provenance[0]
    assert qc_path is None
    assert row['manual_review_status'] == 'NOT_REVIEWED'
    assert row['review_decision'] == ''
    assert 'release-manager' in row['marker_qc_reasons']
    assert 'independently verified' in row['marker_qc_reasons']

    db = Database(str(tmp_path / 'p.sqlite'))
    db.initialise()
    db.insert_sequences([('I1', 'ACGTACGT')], dataset='P')
    upsert_active_marker_evidence(db, [('I1', 'ACGTACGT')], provenance, source_type='PARTNER_FASTA', dataset='P')
    with db.connect() as conn:
        tier, status, decision = conn.execute(
            "SELECT evidence_tier,manual_review_status,review_decision FROM marker_evidence WHERE sequence_id='I1'"
        ).fetchone()
    assert (tier, status, decision) == ('LOW', 'NOT_REVIEWED', '')
