"""v1 tests for partner FASTA, MALDI, genomic 16S and evidence pools (types B and D)."""
from __future__ import annotations

import csv
import shutil
from pathlib import Path

import pytest

from branchmanager.db.interface import Database
from branchmanager.genome_16s_import import import_genome_16s
from branchmanager.marker_evidence import (
    apply_marker_review_decisions,
    query_marker_evidence_pool, register_supplemental_marker_fasta, upsert_active_marker_evidence,
    ensure_marker_evidence_schema, marker_evidence_migration_report, backfill_legacy_marker_evidence,
    _derive_legacy_tier,
)
from branchmanager.pipeline.evidence_policy import evidence_quality
from branchmanager.pipeline.novelty import _write_pool_fasta
from branchmanager.project_state import import_maldi_identifications

SEQ = 'ACGT' * 350


def _db(tmp_path):
    db = Database(str(tmp_path / 'p.sqlite'))
    db.initialise()
    db.upsert_dataset_role('Base', 'baseline', genomes_available=False)
    db.upsert_dataset_role('PartnerA', 'candidate', genomes_available=False)
    return db


def _tsv(path, rows):
    with open(path, 'w', newline='') as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    return path


def _genome_dir(tmp_path, concordance='CONCORDANT_GENUS', quality='FULL_LENGTH_LIKE', gid='G1'):
    d = tmp_path / 'genomic'
    d.mkdir(exist_ok=True)
    _tsv(d / 'genome_16s_summary.tsv', [{
        'genome_id': gid, 'original_genome_id': 'GCA_000001.1', 'input_assembly': '/data/g.fa',
        'gtdbtk_taxonomy': 'd__Bacteria;p__Bacillota;g__Bacillus;s__Bacillus subtilis', 'barrnap_status': 'OK',
        'best_candidate_id': gid, 'best_quality_class': quality, 'best_vsearch_identity_pct': '99.1',
        'best_vsearch_query_coverage_pct': '100', 'gtdbtk_concordance': concordance,
    }])
    _tsv(d / 'taxonomy_assignment_audit.tsv', [{
        'genome_id': gid, 'vsearch_best_reference_taxonomy': 'd__Bacteria;p__Bacillota;g__Bacillus'}])
    (d / 'best_16s_per_genome.fna').write_text(f'>{gid}\n{SEQ}\n')
    return d


def _import_genome(db, genomic_dir, out_dir, **kwargs):
    rows = list(csv.DictReader(open(Path(genomic_dir) / 'genome_16s_summary.tsv'), delimiter='\t'))
    mapping = _tsv(Path(out_dir).parent / 'genome_map.tsv', [
        {'genome_id': row['genome_id'], 'sequence_id': row['genome_id']} for row in rows
    ])
    return import_genome_16s(db, genomic_dir, out_dir, genome_metadata=mapping, **kwargs)


def _registry(db, sid):
    with db.connect() as conn:
        return conn.execute(
            'SELECT source_type,evidence_tier,qc_class,active_for_selection,eligible_primary,genome_16s_concordance '
            'FROM marker_evidence WHERE sequence_id=? ORDER BY source_type', (sid,)).fetchall()


def _sanger(db, sid='S1', seq='GGGG' * 300, qc='PASS_HIGH_CONFIDENCE'):
    db.insert_sequences([(sid, seq)], dataset='PartnerA')
    prov = [{'sequence_id': sid, 'marker_qc_class': qc, 'manual_review_status': 'NOT_REQUIRED'}]
    db.upsert_sequence_provenance(prov)
    upsert_active_marker_evidence(db, [(sid, seq)], prov, source_type='SANGER_AB1', dataset='PartnerA')


def test_high_sanger_is_primary_eligible(tmp_path):
    db = _db(tmp_path)
    _sanger(db)
    assert _registry(db, 'S1') == [('SANGER_AB1', 'HIGH', 'PASS_HIGH_CONFIDENCE', 1, 1, '')]


def test_partner_fasta_with_and_without_qc_provenance(tmp_path):
    db = _db(tmp_path)
    db.insert_sequences([('P1', SEQ), ('P2', SEQ[::-1])], dataset='PartnerA')
    upsert_active_marker_evidence(db, [('P1', SEQ)], [{
        'sequence_id': 'P1', 'marker_qc_class': 'PASS_HIGH_CONFIDENCE', 'manual_review_status': 'NOT_REQUIRED'}],
        source_type='PARTNER_FASTA', dataset='PartnerA')
    upsert_active_marker_evidence(db, [('P2', SEQ[::-1])], [], source_type='PARTNER_FASTA', dataset='PartnerA')
    p1, p2 = _registry(db, 'P1')[0], _registry(db, 'P2')[0]
    assert p1[1] in {'HIGH', 'MODERATE'} and p1[4] == 1
    assert p2[1] == 'LOW' and p2[4] == 0  # unverified partner sequence stays out of the primary pool


def test_maldi_only_isolate_has_no_sequence_and_never_enters_vsearch(tmp_path):
    db = _db(tmp_path)
    m = _tsv(tmp_path / 'maldi.tsv', [{'sequence_id': 'M1', 'partner_id': 'PX',
                                       'taxon_name': 'Bacillus subtilis', 'taxonomic_rank': 'species'}])
    res = import_maldi_identifications(db, m, default_dataset='PartnerM')
    assert res[0]['result'] != 'REJECTED'
    with db.connect() as conn:
        assert conn.execute("SELECT sequence FROM sequences WHERE id='M1'").fetchone()[0] in ('', None)
        assert conn.execute('SELECT COUNT(*) FROM maldi_identifications').fetchone()[0] == 1
        records, _ = query_marker_evidence_pool(conn, pool='primary')
        genome, _ = query_marker_evidence_pool(conn, pool='genome_derived')
    assert 'M1' not in {s for s, _ in records} | {s for s, _ in genome}
    fasta = _write_pool_fasta(records, tmp_path / 'pool.fa')
    assert fasta is None or 'M1' not in Path(fasta).read_text()


def test_sanger_isolate_with_maldi_keeps_marker_and_counts_once(tmp_path):
    db = _db(tmp_path)
    _sanger(db, 'S1')
    m = _tsv(tmp_path / 'maldi.tsv', [{'sequence_id': 'S1', 'partner_id': 'PX',
                                       'taxon_name': 'Clostridium butyricum', 'taxonomic_rank': 'species'}])
    import_maldi_identifications(db, m, default_dataset='PartnerA')
    with db.connect() as conn:
        assert conn.execute("SELECT COUNT(*) FROM sequences WHERE id='S1'").fetchone()[0] == 1
        assert conn.execute("SELECT sequence FROM sequences WHERE id='S1'").fetchone()[0] == 'GGGG' * 300
        assert conn.execute("SELECT COUNT(*) FROM maldi_identifications WHERE sequence_id='S1'").fetchone()[0] == 1
    assert _registry(db, 'S1')[0][:2] == ('SANGER_AB1', 'HIGH')


def test_barrnap_marker_registers_supplementary_and_marks_genome_sequenced(tmp_path):
    db = _db(tmp_path)
    res = _import_genome(db, _genome_dir(tmp_path), tmp_path / 'rep', dataset='PartnerG', partner_id='PG')
    assert res[0]['evidence_tier'] == 'MODERATE'
    rows = _registry(db, 'G1')
    assert rows[0][0] == 'GENOME_BARRNAP' and rows[0][4] == 0 and rows[0][3] == 0
    assert rows[0][5] == 'CONCORDANT_GENUS'
    assert db.get_genome_records()[0]['accession'] == 'GCA_000001.1'
    assert db.get_sequencing_metadata_for_ids()['G1']['selected_for_wgs']
    with db.connect() as conn:
        primary, _ = query_marker_evidence_pool(conn, pool='primary')
        genome, _ = query_marker_evidence_pool(conn, pool='genome_derived')
        # GTDB-Tk taxonomy and VSEARCH marker taxonomy stay in separate columns
        gt, tax = conn.execute("SELECT genome_taxonomy,taxonomy FROM marker_evidence WHERE sequence_id='G1'").fetchone()
    assert 'G1' not in {s for s, _ in primary} and 'G1' in {s for s, _ in genome}
    assert gt and tax


def test_barrnap_discordance_warns_without_rejecting(tmp_path):
    db = _db(tmp_path)
    res = _import_genome(db, _genome_dir(tmp_path, concordance='DISCORDANT_AT_PHYLUM'), tmp_path / 'rep',
                         dataset='PartnerG')
    assert res[0]['warning'] == 'GTDBTK_16S_DISCORDANT'
    assert res[0]['evidence_tier'] == 'MODERATE'
    assert _registry(db, 'G1')[0][5] == 'DISCORDANT_AT_PHYLUM'


def test_truncated_barrnap_prediction_is_low(tmp_path):
    db = _db(tmp_path)
    res = _import_genome(db, _genome_dir(tmp_path, quality='REVIEW_TRUNCATED'), tmp_path / 'rep', dataset='PartnerG')
    assert res[0]['evidence_tier'] == 'LOW'
    with db.connect() as conn:
        genome, _ = query_marker_evidence_pool(conn, pool='genome_derived')
    assert 'G1' not in {s for s, _ in genome}


def test_barrnap_cannot_replace_accepted_sanger_marker(tmp_path):
    db = _db(tmp_path)
    _sanger(db, 'G1')
    _import_genome(db, _genome_dir(tmp_path), tmp_path / 'rep', dataset='PartnerG')
    by_source = {r[0]: r for r in _registry(db, 'G1')}
    assert by_source['SANGER_AB1'][3:5] == (1, 1)
    assert by_source['GENOME_BARRNAP'][3:5] == (0, 0)
    with db.connect() as conn:
        assert conn.execute("SELECT sequence FROM sequences WHERE id='G1'").fetchone()[0] == 'GGGG' * 300
        assert conn.execute("SELECT COUNT(*) FROM sequences WHERE id='G1'").fetchone()[0] == 1
        primary, _ = query_marker_evidence_pool(conn, pool='primary')
    assert dict(primary)['G1'] == 'GGGG' * 300


def test_duplicate_import_does_not_duplicate_evidence(tmp_path):
    db = _db(tmp_path)
    gd = _genome_dir(tmp_path)
    _import_genome(db, gd, tmp_path / 'r1', dataset='PartnerG')
    _import_genome(db, gd, tmp_path / 'r2', dataset='PartnerG')
    _sanger(db, 'S1')
    _sanger(db, 'S1')
    with db.connect() as conn:
        assert conn.execute("SELECT COUNT(*) FROM marker_evidence WHERE sequence_id='G1'").fetchone()[0] == 1
        assert conn.execute("SELECT COUNT(*) FROM marker_evidence WHERE sequence_id='S1'").fetchone()[0] == 1
        assert conn.execute('SELECT COUNT(*) FROM genome_records').fetchone()[0] == 1


def test_primary_pool_excludes_unreviewed_failed_maldi_and_genome(tmp_path):
    db = _db(tmp_path)
    db.insert_sequences([('B1', SEQ)], dataset='Base')
    _sanger(db, 'S1')
    _sanger(db, 'F1', seq='TTTT' * 300, qc='FAIL_QC')
    _import_genome(db, _genome_dir(tmp_path, gid='G9'), tmp_path / 'rep', dataset='PartnerG')
    ensure_marker_evidence_schema(db)
    backfill_legacy_marker_evidence(db)
    with db.connect() as conn:
        primary = {s for s, _ in query_marker_evidence_pool(conn, pool='primary')[0]}
        secondary = {s for s, _ in query_marker_evidence_pool(conn, pool='genome_derived')[0]}
    assert primary == {'B1', 'S1'}
    assert secondary == {'G9'}


def test_novelty_primary_pool_is_registry_only(tmp_path):
    db = _db(tmp_path)
    db.insert_sequences([('B1', SEQ)], dataset='Base')
    _sanger(db, 'F1', seq='TTTT' * 300, qc='FAIL_QC')
    _import_genome(db, _genome_dir(tmp_path, gid='G9'), tmp_path / 'rep', dataset='PartnerG')
    ensure_marker_evidence_schema(db)
    backfill_legacy_marker_evidence(db)
    with db.connect() as conn:
        ids = {s for s, _ in query_marker_evidence_pool(conn, pool='primary')[0]}
    assert 'F1' not in ids and 'G9' not in ids


def test_quarterly_policy_high_moderate_low():
    assert evidence_quality({'marker_source_type': 'SANGER_AB1', 'marker_qc_class': 'PASS_HIGH_CONFIDENCE',
                             'marker_evidence_tier': 'HIGH'}) == 'HIGH'
    assert evidence_quality({'marker_source_type': 'GENOME_BARRNAP', 'marker_evidence_tier': 'HIGH'}) == 'MODERATE'
    assert evidence_quality({'marker_source_type': 'SANGER_AB1', 'marker_qc_class': 'FAIL_QC',
                             'marker_manual_review_status': 'NOT_REVIEWED'}) == 'LOW'
    assert evidence_quality({'marker_source_type': 'SANGER_AB1', 'marker_qc_class': 'FAIL_QC',
                             'marker_manual_review_status': 'APPROVED',
                             'marker_review_decision': 'APPROVE_FAIL_QC',
                             'marker_evidence_tier': 'MODERATE'}) == 'MODERATE'


def test_migration_on_db_copy_is_idempotent_and_preserves_ids(tmp_path):
    db = _db(tmp_path)
    _sanger(db, 'S1')
    db.insert_sequences([('L1', SEQ)], dataset='PartnerA')
    live = tmp_path / 'p.sqlite'
    copy = tmp_path / 'copy.sqlite'
    shutil.copy(live, copy)
    cdb = Database(str(copy))
    ensure_marker_evidence_schema(cdb)
    first = backfill_legacy_marker_evidence(cdb)
    second = backfill_legacy_marker_evidence(cdb)
    with cdb.connect() as conn:
        ids = {r[0] for r in conn.execute('SELECT id FROM sequences')}
        n = conn.execute('SELECT COUNT(*) FROM marker_evidence').fetchone()[0]
    assert ids == {'S1', 'L1'} and n == 2
    assert marker_evidence_migration_report(cdb) is not None
    assert first is not None and second is not None


def test_reject_then_identical_reimport_preserves_review_and_pool_exclusion(tmp_path):
    db = _db(tmp_path)
    _sanger(db, 'R1')
    decisions = _tsv(tmp_path / 'reject.tsv', [{
        'sequence_id': 'R1', 'decision': 'REJECT', 'reviewer': 'reviewer-1', 'notes': 'trace failed inspection',
    }])
    assert apply_marker_review_decisions(db, decisions)[0]['result'] == 'APPLIED'

    # A routine re-import writes empty/NOT_REQUIRED review metadata, not a new decision.
    _sanger(db, 'R1')
    with db.connect() as conn:
        evidence = conn.execute(
            "SELECT manual_review_status,review_decision,reviewer,review_notes,active_for_selection,eligible_primary "
            "FROM marker_evidence WHERE sequence_id='R1' AND source_type='SANGER_AB1'"
        ).fetchone()
        provenance = conn.execute(
            "SELECT manual_review_status,review_decision,reviewer,review_notes FROM sequence_provenance WHERE sequence_id='R1'"
        ).fetchone()
        primary, _ = query_marker_evidence_pool(conn, pool='primary')
    assert evidence == ('REJECTED', 'REJECT', 'reviewer-1', 'trace failed inspection', 0, 0)
    assert provenance == ('REJECTED', 'REJECT', 'reviewer-1', 'trace failed inspection')
    assert 'R1' not in {sid for sid, _ in primary}


def test_legacy_backfill_rejection_precedes_high_and_warning_qc():
    rejected_high = {
        'marker_qc_class': 'PASS_HIGH_CONFIDENCE',
        'manual_review_status': 'REJECTED',
        'review_decision': 'REJECT',
    }
    rejected_warning = {
        'marker_qc_class': 'PASS_WITH_WARNINGS',
        'manual_review_status': 'REJECTED',
        'review_decision': 'REJECT',
    }
    assert _derive_legacy_tier('candidate', rejected_high) == ('LEGACY_UNKNOWN', 'LOW', 'PASS_HIGH_CONFIDENCE', 0)
    assert _derive_legacy_tier('candidate', rejected_warning) == ('LEGACY_UNKNOWN', 'LOW', 'PASS_WITH_WARNINGS', 0)
    assert _derive_legacy_tier('baseline', rejected_high) == ('LEGACY_BASELINE', 'LOW', 'CURATED_BASELINE', 0)


def test_rejected_legacy_high_and_warning_records_backfill_low(tmp_path):
    db = _db(tmp_path)
    db.insert_sequences([('RH', SEQ), ('RW', SEQ[::-1])], dataset='PartnerA')
    db.upsert_sequence_provenance([
        {'sequence_id': 'RH', 'marker_qc_class': 'PASS_HIGH_CONFIDENCE', 'manual_review_status': 'REJECTED',
         'review_decision': 'REJECT', 'reviewer': 'reviewer-1', 'review_notes': 'rejected high'},
        {'sequence_id': 'RW', 'marker_qc_class': 'PASS_WITH_WARNINGS', 'manual_review_status': 'REJECTED',
         'review_decision': 'REJECT', 'reviewer': 'reviewer-1', 'review_notes': 'rejected warning'},
    ])
    backfill_legacy_marker_evidence(db)
    with db.connect() as conn:
        rows = conn.execute(
            "SELECT sequence_id,evidence_tier,eligible_primary FROM marker_evidence ORDER BY sequence_id"
        ).fetchall()
    assert rows == [('RH', 'LOW', 0), ('RW', 'LOW', 0)]


def test_genome_qc_unknown_is_nullable_and_mapping_is_explicit(tmp_path):
    db = _db(tmp_path)
    genomic = _genome_dir(tmp_path, gid='G_UNKNOWN')
    with pytest.raises(ValueError, match='Genome-to-isolate mapping is required'):
        import_genome_16s(db, genomic, tmp_path / 'must_fail', dataset='PartnerG')
    mapping = _tsv(tmp_path / 'mapping.tsv', [{'genome_id': 'G_UNKNOWN', 'sequence_id': 'I_UNKNOWN'}])
    import_genome_16s(db, genomic, tmp_path / 'unknown', dataset='PartnerG', genome_metadata=mapping)
    record = db.get_genome_records()[0]
    assert record['genome_status'] == 'SEQUENCED' and record['genome_qc_pass'] is None


def test_quarterly_review_uses_passed_gtdb_genomes_as_species_anchors(tmp_path):
    from branchmanager.pipeline.quarterly_review import _available_ids_by_species
    db = _db(tmp_path)
    genomic = _genome_dir(tmp_path, gid='G_PASS')
    mapping = _tsv(tmp_path / 'mapping.tsv', [{
        'genome_id': 'G_PASS', 'sequence_id': 'I_PASS', 'genome_qc_pass': 'true',
    }])
    import_genome_16s(db, genomic, tmp_path / 'passed', dataset='PartnerG', genome_metadata=mapping)
    anchors = _available_ids_by_species(db)
    assert anchors.get('bacillussubtilis') == ['I_PASS']


def test_failed_genome_is_available_but_candidate_can_be_reselected(tmp_path):
    from branchmanager.pipeline.quarterly_review import build_quarterly_review
    db = _db(tmp_path)
    _sanger(db, 'I_FAIL')
    genomic = _genome_dir(tmp_path, gid='G_FAIL')
    mapping = _tsv(tmp_path / 'mapping.tsv', [{
        'genome_id': 'G_FAIL', 'sequence_id': 'I_FAIL', 'genome_qc_pass': 'false',
    }])
    import_genome_16s(db, genomic, tmp_path / 'failed', dataset='PartnerG', genome_metadata=mapping)
    row = {
        'id': 'I_FAIL', 'partner_id': 'PG', 'taxonomy': 'd__Bacteria; g__Bacillus; s__Bacillus subtilis',
        'classification_identity': '99.5', 'classification_confidence': '0.99',
        'nearest_identity': '95.0', 'project_nearest_identity': '95.0', 'reference_nearest_identity': '95.0',
        'selected_for_genome_sequencing': 'False', 'already_sequenced': 'False',
        'genome_committed_count_same_species': '0', 'pangenome_target': '3', 'placement_flags': '',
        'marker_source_type': 'SANGER_AB1', 'marker_evidence_tier': 'HIGH',
        'marker_qc_class': 'PASS_HIGH_CONFIDENCE',
    }
    recommendation = build_quarterly_review([row], db=db, genome_budget=1, pangenome_target=3)[0]
    assert recommendation['genome_available'] == 'True'
    assert recommendation['genome_qc_status'] == 'GENOME_QC_FAILED'
    assert recommendation['genome_resequencing_eligible'] == 'True'
    assert recommendation['role'] == 'PRIMARY'


def test_multiple_genomes_for_one_isolate_do_not_inflate_species_anchors(tmp_path):
    from branchmanager.pipeline.quarterly_review import _available_ids_by_species
    db = _db(tmp_path)
    first = _genome_dir(tmp_path, gid='G_A')
    mapping = _tsv(tmp_path / 'mapping.tsv', [
        {'genome_id': 'G_A', 'sequence_id': 'I_SHARED', 'genome_qc_pass': 'true'},
    ])
    import_genome_16s(db, first, tmp_path / 'a', dataset='PartnerG', genome_metadata=mapping)
    second = _genome_dir(tmp_path, gid='G_B')
    mapping = _tsv(tmp_path / 'mapping.tsv', [
        {'genome_id': 'G_A', 'sequence_id': 'I_SHARED', 'genome_qc_pass': 'true'},
        {'genome_id': 'G_B', 'sequence_id': 'I_SHARED', 'genome_qc_pass': 'true'},
    ])
    # Re-import with both mappings records two assemblies attached to one isolate.
    import_genome_16s(db, second, tmp_path / 'b', dataset='PartnerG', genome_metadata=mapping)
    assert _available_ids_by_species(db)['bacillussubtilis'] == ['I_SHARED']


def test_secondary_genome_pool_is_reported_separately_with_vsearch(tmp_path):
    import shutil as sh
    import pytest
    if not sh.which('vsearch'):
        pytest.skip('vsearch not installed')
    from branchmanager.pipeline.novelty import _write_genome_derived_context
    db = _db(tmp_path)
    _import_genome(db, _genome_dir(tmp_path, gid='G9'), tmp_path / 'rep', dataset='PartnerG')
    query = tmp_path / 'q.fa'
    query.write_text(f'>Q1\n{SEQ}\n')
    out = tmp_path / 'ctx'
    out.mkdir()
    summary = _write_genome_derived_context(str(query), str(out), db, ['Q1'])
    row = next(csv.DictReader(open(summary), delimiter='\t'))
    assert row['genome_derived_nearest_hit'] == 'G9' and float(row['identity']) > 99


def test_ledgers_count_each_isolate_once_and_flag_maldi_disagreement(tmp_path):
    from branchmanager.v1_ledgers import v1_summary_counts, write_v1_ledgers
    db = _db(tmp_path)
    _sanger(db, 'G1')
    _import_genome(db, _genome_dir(tmp_path, concordance='DISCORDANT_AT_PHYLUM'), tmp_path / 'rep', dataset='PartnerG')
    import_maldi_identifications(db, _tsv(tmp_path / 'm.tsv', [{
        'sequence_id': 'G1', 'partner_id': 'PX', 'taxon_name': 'Escherichia coli', 'taxonomic_rank': 'species'}]),
        default_dataset='PartnerA')
    counts = v1_summary_counts(db)
    assert [n for s, label, n in counts if label == 'Total isolates'] == [1]
    paths = write_v1_ledgers(db, tmp_path / 'ledgers')
    for name in ('project_taxonomy_concordance.tsv', 'project_unreviewed_markers.tsv', 'project_genome_ledger.tsv'):
        assert (tmp_path / 'ledgers' / name).is_file()
    assert paths


def test_maldi_counts_distinguish_direct_and_genomic_marker_evidence(tmp_path):
    from branchmanager.v1_ledgers import v1_summary_counts
    db = _db(tmp_path)
    for sid, sequence in [('M_ONLY', ''), ('M_DIRECT', 'GGGG' * 300), ('M_GENOME', ''), ('M_BOTH', 'CCCC' * 300)]:
        db.insert_sequences([(sid, sequence)], dataset='PartnerA')
    _sanger(db, 'M_DIRECT', seq='GGGG' * 300)
    _sanger(db, 'M_BOTH', seq='CCCC' * 300)
    genome = _genome_dir(tmp_path, gid='G_M')
    mapping = _tsv(tmp_path / 'genome_map.tsv', [{'genome_id': 'G_M', 'sequence_id': 'M_GENOME'}])
    import_genome_16s(db, genome, tmp_path / 'genome-only', dataset='PartnerG', genome_metadata=mapping)
    genome = _genome_dir(tmp_path, gid='G_B')
    mapping = _tsv(tmp_path / 'genome_map.tsv', [{'genome_id': 'G_B', 'sequence_id': 'M_BOTH'}])
    import_genome_16s(db, genome, tmp_path / 'genome-both', dataset='PartnerG', genome_metadata=mapping)
    maldi = _tsv(tmp_path / 'maldi.tsv', [
        {'sequence_id': sid, 'partner_id': 'PM', 'taxon_name': 'Bacillus subtilis', 'taxonomic_rank': 'species'}
        for sid in ('M_ONLY', 'M_DIRECT', 'M_GENOME', 'M_BOTH')
    ])
    import_maldi_identifications(db, maldi)
    labels = {label: n for section, label, n in v1_summary_counts(db) if section == 'MALDI'}
    assert labels['MALDI only (no direct or genomic marker)'] == 1
    assert labels['MALDI plus direct marker only'] == 1
    assert labels['MALDI plus genomic 16S only'] == 1
    assert labels['MALDI plus direct and genomic marker evidence'] == 1


def test_per_observation_concordance_ledger_preserves_multiplicity(tmp_path):
    from branchmanager.v1_ledgers import write_v1_ledgers
    db = _db(tmp_path)
    _sanger(db, 'I1')
    for gid in ('G1', 'G2'):
        genomic = _genome_dir(tmp_path, gid=gid, concordance='DISCORDANT_AT_PHYLUM')
        mapping = _tsv(tmp_path / 'mapping.tsv', [{'genome_id': gid, 'sequence_id': 'I1'}])
        import_genome_16s(db, genomic, tmp_path / f'import_{gid}', dataset='PartnerG', genome_metadata=mapping)
    maldi = _tsv(tmp_path / 'maldi.tsv', [
        {'sequence_id': 'I1', 'partner_id': 'PM', 'taxon_name': name, 'taxonomic_rank': 'species'}
        for name in ('Bacillus subtilis', 'Escherichia coli')
    ])
    import_maldi_identifications(db, maldi)
    paths = write_v1_ledgers(db, tmp_path / 'ledgers')
    observations = list(csv.DictReader(open(paths['taxonomy_observations']), delimiter='\t'))
    i1 = [row for row in observations if row['sequence_id'] == 'I1']
    assert sum(row['observation_type'] == 'GENOME_RECORD' for row in i1) == 2
    assert sum(row['observation_type'] == 'MALDI_IDENTIFICATION' for row in i1) == 2
    assert sum(row['source_type'] == 'GENOME_BARRNAP' for row in i1) == 2


def test_stored_status_only_approval_is_demoted_on_reimport(tmp_path):
    db = _db(tmp_path)
    _sanger(db, 'W1', qc='PASS_WITH_WARNINGS')
    with db.connect() as conn:
        conn.execute(
            "UPDATE marker_evidence SET manual_review_status='APPROVED',review_decision='',"
            "evidence_tier='MODERATE',eligible_primary=1 WHERE sequence_id='W1'")
    _sanger(db, 'W1', qc='PASS_WITH_WARNINGS')
    with db.connect() as conn:
        row = conn.execute(
            "SELECT manual_review_status,review_decision,evidence_tier,eligible_primary "
            "FROM marker_evidence WHERE sequence_id='W1'").fetchone()
    assert row == ('NOT_REVIEWED', '', 'LOW', 0)


def test_failed_genome_commitment_does_not_reduce_pangenome_gap(tmp_path):
    from branchmanager.pipeline.quarterly_review import build_quarterly_review
    db = _db(tmp_path)
    _sanger(db, 'I_FAIL')
    mapping = _tsv(tmp_path / 'mapping.tsv', [{
        'genome_id': 'G_FAIL', 'sequence_id': 'I_FAIL', 'genome_qc_pass': 'false'}])
    import_genome_16s(db, _genome_dir(tmp_path, gid='G_FAIL'), tmp_path / 'f', dataset='PartnerG',
                      genome_metadata=mapping)
    row = {
        'id': 'I_FAIL', 'partner_id': 'PG', 'taxonomy': 'd__Bacteria; g__Bacillus; s__Bacillus subtilis',
        'classification_identity': '99.5', 'classification_confidence': '0.99',
        'nearest_identity': '95.0', 'project_nearest_identity': '95.0', 'reference_nearest_identity': '95.0',
        'selected_for_genome_sequencing': 'True', 'already_sequenced': 'True',
        'genome_committed_count_same_species': '0', 'pangenome_target': '3', 'placement_flags': '',
        'marker_source_type': 'SANGER_AB1', 'marker_evidence_tier': 'HIGH',
        'marker_qc_class': 'PASS_HIGH_CONFIDENCE',
    }
    rec = build_quarterly_review([row], db=db, genome_budget=1, pangenome_target=3)[0]
    assert str(rec['available_genomes_before']) == '0'
