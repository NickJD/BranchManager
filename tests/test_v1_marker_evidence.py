"""Regression tests for BranchManager v1 marker-evidence policy.

Run from the BranchManager repository root with:
    python -m pytest tests/test_v1_marker_evidence.py -q
"""
from __future__ import annotations

import csv
import sqlite3
from pathlib import Path

import pytest

from branchmanager.marker_evidence import (
    backfill_legacy_marker_evidence,
    ensure_marker_evidence_schema,
    query_marker_evidence_pool,
    register_supplemental_marker_fasta,
)
from branchmanager.pipeline.paper_trail import _classify_ab1_quality_values
from branchmanager.pipeline.quarterly_review import _evidence_quality


class SQLiteProject:
    """Minimal Database-compatible wrapper for isolated evidence-registry tests."""
    def __init__(self, path: Path):
        self.path = str(path)

    def connect(self):
        return sqlite3.connect(self.path)


def make_db(path: Path) -> SQLiteProject:
    db = SQLiteProject(path)
    with db.connect() as conn:
        conn.executescript(
            """
            CREATE TABLE sequences(id TEXT PRIMARY KEY, sequence TEXT, dataset TEXT);
            CREATE TABLE dataset_roles(dataset TEXT PRIMARY KEY, role TEXT);
            CREATE TABLE sequence_provenance(
                sequence_id TEXT PRIMARY KEY,
                marker_qc_class TEXT,
                manual_review_status TEXT,
                review_decision TEXT,
                marker_qc_reasons TEXT
            );
            INSERT INTO sequences VALUES
                ('BASE_1', 'ACGTACGT', 'Hungate'),
                ('SANGER_1', 'ACGTACGA', 'PartnerA'),
                ('SANGER_2', 'ACGTACGC', 'PartnerA'),
                ('SANGER_FAIL', 'ACGTACGG', 'PartnerA'),
                ('MALDI_ONLY', '', 'PartnerA');
            INSERT INTO dataset_roles VALUES ('Hungate', 'baseline');
            INSERT INTO dataset_roles VALUES ('PartnerA', 'candidate');
            INSERT INTO sequence_provenance VALUES
                ('SANGER_1', 'PASS_HIGH_CONFIDENCE', 'NOT_REQUIRED', '', ''),
                ('SANGER_2', 'PASS_WITH_WARNINGS', 'APPROVED', 'APPROVE', 'minor trace warning'),
                ('SANGER_FAIL', 'FAIL_QC', 'NOT_REVIEWED', '', 'poor trace');
            """
        )
    return db


def test_zero_filled_ab1_pcon_is_not_quality_available():
    assert _classify_ab1_quality_values([0] * 1416, 1416) == (False, 'AB1_ZERO_FILLED_PCON')
    assert _classify_ab1_quality_values([], 1416) == (False, 'AB1_MISSING_PCON')
    assert _classify_ab1_quality_values([30] * 800, 1416) == (False, 'AB1_PARTIAL_PCON')
    assert _classify_ab1_quality_values([0, 30, 20], 3) == (True, 'AB1_PCON_OK')


def test_legacy_backfill_is_conservative_and_primary_pool_excludes_failed_and_markerless(tmp_path):
    db = make_db(tmp_path / 'project.sqlite')
    ensure_marker_evidence_schema(db)
    counts = backfill_legacy_marker_evidence(db)
    assert counts['baseline_high'] == 1
    assert counts['partner_high'] == 1
    assert counts['partner_moderate'] == 1
    assert counts['partner_low'] == 1

    with db.connect() as conn:
        records, datasets = query_marker_evidence_pool(conn, pool='primary')
    ids = {sequence_id for sequence_id, _ in records}
    assert ids == {'BASE_1', 'SANGER_1', 'SANGER_2'}
    assert 'SANGER_FAIL' not in ids
    assert 'MALDI_ONLY' not in ids


def test_genome_barrnap_is_secondary_only_and_stays_separate_from_active_marker(tmp_path):
    db = make_db(tmp_path / 'project.sqlite')
    ensure_marker_evidence_schema(db)
    backfill_legacy_marker_evidence(db)

    fasta = tmp_path / 'genomic_16s.fna'
    fasta.write_text('>GCF_000001.1\n' + 'ACGT' * 375 + '\n')
    metadata = tmp_path / 'genomic_16s_metadata.tsv'
    metadata.write_text(
        'sequence_id\tfasta_id\tgenome_id\tsource_type\tevidence_tier\tqc_class\t'
        'genome_taxonomy\ttaxonomy\tgenome_16s_concordance\n'
        'SANGER_1\tGCF_000001.1\tGCF_000001.1\tGENOME_BARRNAP\tMODERATE\tGENOME_DERIVED\t'
        'd__Bacteria;p__Bacillota\td__Bacteria;p__Bacillota\tPHYLUM_CONCORDANT\n'
    )
    result = register_supplemental_marker_fasta(db, fasta, metadata)
    assert result[0]['result'] == 'IMPORTED'
    assert result[0]['evidence_tier'] == 'MODERATE'

    with db.connect() as conn:
        primary, _ = query_marker_evidence_pool(conn, pool='primary')
        genome, _ = query_marker_evidence_pool(conn, pool='genome_derived')
        observations = conn.execute(
            "SELECT source_type,evidence_tier,active_for_selection,eligible_primary,genome_16s_concordance "
            "FROM marker_evidence WHERE sequence_id='SANGER_1' AND source_type='GENOME_BARRNAP'"
        ).fetchone()
    assert 'SANGER_1' in {sid for sid, _ in primary}  # existing Sanger observation remains active
    assert 'SANGER_1' in {sid for sid, _ in genome}  # predicted 16S is separately visible
    assert observations == ('GENOME_BARRNAP', 'MODERATE', 0, 0, 'PHYLUM_CONCORDANT')


def test_quarterly_review_honours_source_tier_and_manual_approval():
    assert _evidence_quality({
        'marker_source_type': 'GENOME_BARRNAP', 'marker_evidence_tier': 'HIGH',
        'marker_qc_class': 'GENOME_DERIVED',
    }) == 'MODERATE'
    assert _evidence_quality({
        'marker_source_type': 'SANGER_AB1', 'marker_evidence_tier': 'LOW',
        'marker_qc_class': 'FAIL_QC', 'marker_manual_review_status': 'NOT_REVIEWED',
    }) == 'LOW'
    assert _evidence_quality({
        'marker_source_type': 'SANGER_AB1', 'marker_evidence_tier': 'MODERATE',
        'marker_qc_class': 'FAIL_QC', 'marker_manual_review_status': 'APPROVED',
        'marker_review_decision': 'APPROVE_FAIL_QC',
    }) == 'MODERATE'
