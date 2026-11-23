import json
import sqlite3
from pathlib import Path

from branchmanager.continuation import (
    backup_sqlite,
    create_version_dir,
    export_candidate_marker_fasta,
    latest_performance_review_manifest,
)
from branchmanager.db.interface import Database
from branchmanager.utils.fasta import read_fasta


def test_backup_sqlite_creates_valid_immutable_copy(tmp_path):
    source = tmp_path / 'project.sqlite'
    with sqlite3.connect(source) as conn:
        conn.execute('CREATE TABLE sample (value TEXT)')
        conn.execute('INSERT INTO sample VALUES (?)', ('before',))
    snapshot = backup_sqlite(source, tmp_path / 'versions' / 'previous.sqlite')
    with sqlite3.connect(source) as conn:
        conn.execute('UPDATE sample SET value=?', ('after',))
    with sqlite3.connect(snapshot) as conn:
        assert conn.execute('SELECT value FROM sample').fetchone()[0] == 'before'
        assert conn.execute('PRAGMA integrity_check').fetchone()[0] == 'ok'


def test_create_version_dir_never_overwrites_existing_version(tmp_path):
    first = create_version_dir(tmp_path, stamp='20261120_120000')
    second = create_version_dir(tmp_path, stamp='20261120_120000')
    assert first.name == '20261120_120000'
    assert second.name == '20261120_120000_01'
    assert first != second


def test_export_candidate_fasta_excludes_baselines_and_empty_markers(tmp_path):
    db_path = tmp_path / 'test.sqlite'
    db = Database(str(db_path))
    db.initialise()
    db.insert_sequences([('BASE_01', 'AAAA')], dataset='Baseline')
    db.upsert_dataset_role('Baseline', 'baseline', baseline_tier='priority')
    db.insert_sequences([('ISO_01', 'ACGT'), ('ISO_02', '')], dataset='Partner_01')
    db.upsert_dataset_role('Partner_01', 'candidate')
    fasta_path = export_candidate_marker_fasta(db, tmp_path / 'candidate.fasta')
    records = list(read_fasta(fasta_path))
    assert records == [('ISO_01', 'ACGT')]


def test_latest_performance_review_manifest_is_found(tmp_path):
    db_path = tmp_path / 'project.sqlite'
    manifest = tmp_path / 'review' / 'run_manifest.json'
    manifest.parent.mkdir()
    manifest.write_text(json.dumps({'workflow': 'performance_review'}))
    with sqlite3.connect(db_path) as conn:
        conn.execute('CREATE TABLE project_runs (run_id TEXT, workflow TEXT, dataset TEXT, status TEXT, manifest_path TEXT, started_at TEXT, completed_at TEXT, error TEXT)')
        conn.execute('INSERT INTO project_runs VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
                    ('r1', 'performance_review', 'batch', 'COMPLETE', str(manifest), '', '', ''))
    assert latest_performance_review_manifest(db_path) == manifest.resolve()


def test_inherit_review_args_prefers_absolute_manifest_inputs(tmp_path):
    from branchmanager.continuation import inherit_performance_review_args

    ref = tmp_path / 'ref.fasta'
    ref.write_text('>ref\nACGT\n')
    seqs = tmp_path / 'input.fasta'
    seqs.write_text('>iso\nACGT\n')
    db_path = tmp_path / 'project.sqlite'
    out = tmp_path / 'review'
    manifest = tmp_path / 'run_manifest.json'
    manifest.write_text(json.dumps({
        'command': ['branchmanager', 'performance-review', '--input', 'relative.fasta', '--db', 'relative.db', '--out', 'relative_out', '--dataset', 'Batch', '--ref', 'relative_ref.fasta'],
        'inputs': [{'role': 'primary_reference', 'path': str(ref), 'exists': True}],
    }))
    recovered, data = inherit_performance_review_args(manifest)
    assert recovered.ref == str(ref)
    assert recovered.dataset == 'Batch'
    assert data['command'][1] == 'performance-review'


def test_continue_project_publishes_new_version_only_after_mock_assessment(tmp_path, monkeypatch):
    import csv
    import argparse
    import sys
    from branchmanager.cli_core import cmd_continue_project

    root = tmp_path
    db_path = root / 'project.sqlite'
    db = Database(str(db_path))
    db.initialise()
    db.insert_sequences([('AU_0001', 'ACGT')], dataset='AU_01')
    db.upsert_dataset_role('AU_01', 'candidate')
    db.upsert_sequencing_metadata([{
        'id': 'AU_0001', 'partner_id': 'AU', 'dataset': 'AU_01',
        'selected_for_sequencing': False, 'selected_for_wgs': False,
        'source_id': 'AU_0001', 'project_metadata': {'partner_own_id': 'partner-1'},
    }])
    ref = root / 'reference.fasta'
    ref.write_text('>ref\nACGT\n')
    prior_out = root / 'previous_review'
    prior_out.mkdir()
    prior_manifest = prior_out / 'run_manifest.json'
    prior_manifest.write_text(json.dumps({
        'workflow': 'performance_review',
        'command': ['branchmanager', 'performance-review', '--input', 'old.fasta', '--db', 'old.sqlite', '--out', str(prior_out), '--dataset', 'Batch_01', '--ref', str(ref)],
        'inputs': [{'role': 'primary_reference', 'path': str(ref), 'exists': True}],
        'outputs': [{'role': 'full_assessment', 'path': str(prior_out / 'assessment' / 'sequence_assessment.tsv')}],
    }))
    with sqlite3.connect(db_path) as conn:
        conn.execute('INSERT INTO project_runs VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
                    ('review-old', 'performance_review', 'Batch_01', 'COMPLETE', str(prior_manifest), '2026-11-19', '2026-11-19', ''))
    old_assessment = prior_out / 'assessment'
    old_assessment.mkdir()
    (old_assessment / 'sequence_assessment.tsv').write_text('sequence_id\nAU_0001\n')

    catalogue = root / 'catalogue.csv'
    with catalogue.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=['sequence_id', 'partner_id', 'selected_for_genome_sequencing', 'already_sequenced', 'addition_method', 'taxon', 'partner_own_ids'])
        writer.writeheader()
        writer.writerow({'sequence_id': 'AU_0001', 'partner_id': 'AU', 'selected_for_genome_sequencing': 'no', 'already_sequenced': 'no', 'addition_method': 'marker_sequence', 'taxon': 'Bacteria', 'partner_own_ids': 'partner-1'})
        writer.writerow({'sequence_id': 'AU_0002', 'partner_id': 'AU', 'selected_for_genome_sequencing': 'yes', 'already_sequenced': 'yes', 'addition_method': '', 'taxon': 'Bacteria', 'partner_own_ids': 'partner-2'})
    genomes = root / 'genomes.tsv'
    genomes.write_text('sequence_id\tgenome_id\tgenome_status\tgenome_qc_pass\tcompleteness\tcontamination\tgtdb_taxonomy\nAU_0002\tGCA_TEST.1\tGENOME_QC_PASSED\ttrue\t98\t1\td__Bacteria;p__Bacillota\n')

    def fake_review(review_args):
        review_root = Path(review_args.out)
        assessment = review_root / 'assessment'
        assessment.mkdir(parents=True, exist_ok=True)
        for name in ('sequence_assessment.tsv', 'selection_summary.tsv', 'novelty_metrics.tsv', 'pangenome_hierarchy.tsv', 'sequencing_sets.tsv', 'genome_collection_assessment.tsv'):
            (assessment / name).write_text('sequence_id\nAU_0001\n')
        (review_root / 'performance_review_dashboard.html').write_text('<html>mock</html>')
        (review_root / 'run_manifest.json').write_text('{}')
        return None

    monkeypatch.setattr('branchmanager.cli_core.cmd_performance_review', fake_review)
    args = argparse.Namespace(
        db=str(db_path), partner_catalogue=str(catalogue), genome_results=str(genomes),
        partner_workbook=None, partner_workbook_sheet='Bulk Upload Template', partner_workbook_dataset='AU_01',
        marker_fasta=None, marker_qc=None, marker_review=None, replacement_policy='prompt', marker_dataset='Continuation_Markers',
        ref=None, taxa=None, main_ref=None, mwl=None, mwl_sheet=None, threads=None,
        baseline_fasta=None, baseline_taxa_assignments=None, chimera_ref=None, anchors=None,
        alt_ref=None, alt_taxa=None, alt_ref_name=None, sequence_domain=None,
        accept_unverified_marker_qc=False, unverified_override_authorizer='', unverified_override_reason='',
        version_root=str(root / 'versions'), dry_run=False, out=str(root / 'continue_output'),
    )
    cmd_continue_project(args)
    assert 'AU_0002' in db.get_all_ids()
    with sqlite3.connect(db_path) as conn:
        assert conn.execute('PRAGMA integrity_check').fetchone()[0] == 'ok'
        genome_count = conn.execute('SELECT COUNT(*) FROM genome_records WHERE sequence_id="AU_0002"').fetchone()[0]
    assert genome_count == 1
    versions = sorted(path for path in (root / 'versions').iterdir() if path.is_dir())
    assert len(versions) == 1
    version = versions[0]
    assert (version / 'previous_project.sqlite').is_file()
    assert (version / 'project.sqlite').is_file()
    assert (version / 'assessment' / 'performance_review_dashboard.html').is_file()
    assert json.loads((version / 'version_manifest.json').read_text())['published'] is True
    assert (root / 'versions' / 'latest.json').is_file()


def test_continue_project_failed_assessment_does_not_publish(tmp_path, monkeypatch):
    import csv
    import argparse
    from branchmanager.cli_core import cmd_continue_project

    db_path = tmp_path / 'project.sqlite'
    db = Database(str(db_path))
    db.initialise()
    db.insert_sequences([('AU_0001', 'ACGT')], dataset='AU_01')
    db.upsert_dataset_role('AU_01', 'candidate')
    ref = tmp_path / 'reference.fasta'
    ref.write_text('>ref\nACGT\n')
    prior_dir = tmp_path / 'previous'
    prior_dir.mkdir()
    prior_manifest = prior_dir / 'run_manifest.json'
    prior_manifest.write_text(json.dumps({
        'command': ['branchmanager', 'performance-review', '--input', 'old.fasta', '--db', 'old.db', '--out', str(prior_dir), '--dataset', 'Batch', '--ref', str(ref)],
        'inputs': [{'role': 'primary_reference', 'path': str(ref), 'exists': True}],
        'outputs': [],
    }))
    with sqlite3.connect(db_path) as conn:
        conn.execute('INSERT INTO project_runs VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
                    ('review', 'performance_review', 'Batch', 'COMPLETE', str(prior_manifest), '', '', ''))
    catalogue = tmp_path / 'catalogue.csv'
    with catalogue.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=['sequence_id', 'partner_id', 'selected_for_genome_sequencing', 'already_sequenced', 'addition_method', 'taxon'])
        writer.writeheader()
        writer.writerow({'sequence_id': 'AU_0001', 'partner_id': 'AU', 'selected_for_genome_sequencing': 'no', 'already_sequenced': 'no', 'addition_method': 'marker_sequence', 'taxon': 'Bacteria'})
        writer.writerow({'sequence_id': 'AU_0002', 'partner_id': 'AU', 'selected_for_genome_sequencing': 'yes', 'already_sequenced': 'yes', 'addition_method': '', 'taxon': 'Bacteria'})
    genomes = tmp_path / 'genomes.tsv'
    genomes.write_text('sequence_id\tgenome_id\tgenome_status\tgenome_qc_pass\nAU_0002\tGCA_TEST.2\tSEQUENCED\ttrue\n')
    def failing_review(_review_args):
        raise RuntimeError('simulated assessment failure')
    monkeypatch.setattr('branchmanager.cli_core.cmd_performance_review', failing_review)
    args = argparse.Namespace(
        db=str(db_path), partner_catalogue=str(catalogue), genome_results=str(genomes), partner_workbook=None,
        partner_workbook_sheet=None, partner_workbook_dataset='AU_01', marker_fasta=None, marker_qc=None,
        marker_review=None, taxa_assignments=None, replacement_policy='prompt', marker_dataset='Continuation_Markers',
        ref=None, taxa=None, main_ref=None, mwl=None, mwl_sheet=None, threads=None, baseline_fasta=None,
        baseline_taxa_assignments=None, chimera_ref=None, anchors=None, alt_ref=None, alt_taxa=None, alt_ref_name=None,
        sequence_domain=None, accept_unverified_marker_qc=False, unverified_override_authorizer='',
        unverified_override_reason='', version_root=str(tmp_path / 'versions'), dry_run=False,
        out=str(tmp_path / 'outputs'),
    )
    import pytest
    with pytest.raises(RuntimeError, match='simulated assessment failure'):
        cmd_continue_project(args)
    assert 'AU_0002' not in db.get_all_ids()
    with sqlite3.connect(db_path) as conn:
        assert conn.execute('SELECT COUNT(*) FROM genome_records').fetchone()[0] == 0
    assert not (tmp_path / 'versions' / 'latest.json').exists()
    version_dirs = [path for path in (tmp_path / 'versions').iterdir() if path.is_dir()]
    assert len(version_dirs) == 1
    assert json.loads((version_dirs[0] / 'version_manifest.json').read_text())['status'] == 'FAILED'


def test_continue_parser_exposes_marker_replacement_and_assessment_overrides():
    from branchmanager.cli.parser import build_parser
    parsed = build_parser().parse_args([
        'continue-project', '--db', 'project.sqlite', '--partner-catalogue', 'catalogue.csv',
        '--genome-results', 'genomes.tsv', '--marker-fasta', 'incoming.fasta',
        '--replacement-policy', 'approve_qc_pass', '--ref', 'reference.fasta',
        '--taxa-assignments', 'incoming_taxonomy.tsv', '--version-root', 'versions', '-o', 'out',
    ])
    assert parsed.replacement_policy == 'approve_qc_pass'
    assert parsed.marker_fasta == 'incoming.fasta'
    assert parsed.taxa_assignments == 'incoming_taxonomy.tsv'
    assert parsed.version_root == 'versions'
