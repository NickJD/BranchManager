import csv
import json
from pathlib import Path

import pytest
from openpyxl import Workbook

from branchmanager import cli_core
from branchmanager.db.interface import Database
from branchmanager.genome_evidence_import import import_precomputed_genome_results
from branchmanager.partner_metadata import load_genome_only_registration_table
from branchmanager.project_state import import_genome_only_partner_workbook
from branchmanager.release_manager import ReleaseManager


def _write_table(path, rows, *, delimiter='\t'):
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def _genome_row(sequence_id, partner_id, genome_id):
    return {
        'sequence_id': sequence_id, 'partner_id': partner_id, 'genome_id': genome_id,
        'genome_status': 'SEQUENCED', 'completeness': '98', 'contamination': '1',
        'genome_qc_pass': 'yes', 'gtdb_taxonomy': 'd__Bacteria',
        'genome_path': '/example/genome.fasta', 'notes': '',
    }


def test_cli_registers_table_then_imports_genome_results(tmp_path):
    source = tmp_path / 'Genome_Results.tsv'
    _write_table(source, [
        _genome_row('AAFC_0001', 'AAFC', 'LRC159'),
        _genome_row('IAS_CAAS_0001', 'IAS_CAAS', 'GENOME1'),
    ])
    db_path = tmp_path / 'project.sqlite'
    db = Database(str(db_path))
    db.initialise()

    cli_core.main([
        'import-partner-workbook', '--db', str(db_path), '--input', str(source),
        '--genome-only', '--dataset', 'Genome_Only', '--out', str(tmp_path / 'registration'),
    ])

    db = Database(str(db_path))
    assert set(db.get_all_ids()) == {'AAFC_0001', 'IAS_CAAS_0001'}
    metadata = db.get_sequencing_metadata_for_ids()
    assert metadata['AAFC_0001']['partner_id'] == 'AAFC'
    assert metadata['AAFC_0001']['selected_for_wgs'] is False
    assert metadata['AAFC_0001']['project_metadata']['genome_registration_source'] == str(source)
    with db.connect() as conn:
        assert conn.execute('SELECT COUNT(*) FROM genome_records').fetchone()[0] == 0
        assert conn.execute("SELECT sequence FROM sequences WHERE id='AAFC_0001'").fetchone()[0] in ('', None)
    result = import_precomputed_genome_results(db, source)
    assert result['records'] == 2
    with db.connect() as conn:
        assert conn.execute('SELECT COUNT(*) FROM genome_records').fetchone()[0] == 2


def test_csv_registration_is_idempotent_and_preserves_existing_marker(tmp_path):
    source = tmp_path / 'results.csv'
    _write_table(source, [
        _genome_row('UoM_0001', 'UoM', 'G1'),
        _genome_row('UoM_0001', 'UoM', 'G2'),
    ], delimiter=',')
    db = Database(str(tmp_path / 'project.sqlite'))
    db.initialise()
    db.insert_sequences([('UoM_0001', 'ACGT')], dataset='Candidates')
    db.upsert_sequencing_metadata([{
        'id': 'UoM_0001', 'partner_id': 'UoM', 'dataset': 'Candidates',
        'selected_for_wgs': True, 'selected_for_sequencing': False,
        'source_id': 'UoM_0001', 'source_file': 'catalogue',
        'project_metadata': {'partner_culture_id': 'original'},
    }])
    for _ in range(2):
        report = import_genome_only_partner_workbook(db, source, dataset='Candidates')
        assert [row['result'] for row in report] == ['ALREADY_REGISTERED_UPDATED']
    metadata = db.get_sequencing_metadata_for_ids()['UoM_0001']
    assert metadata['selected_for_wgs'] is True
    assert metadata['project_metadata']['partner_culture_id'] == 'original'
    assert 'marker_sequence_status' not in metadata['project_metadata']
    with db.connect() as conn:
        assert conn.execute("SELECT sequence FROM sequences WHERE id='UoM_0001'").fetchone()[0] == 'ACGT'


def test_release_registration_uses_release_database(tmp_path):
    root = tmp_path / 'project'
    root.mkdir()
    config = root / 'project_config.json'
    config.write_text(json.dumps({
        'schema_version': 1, 'project_directory': str(root), 'data': {}, 'defaults': {},
    }))
    manager = ReleaseManager(root)
    manager.init_project('Test', config_source=config, overwrite=True)
    release = manager.create_release('r1')
    release_dir = Path(release['release_dir'])
    db = Database(str(release_dir / 'project.sqlite'))
    db.initialise()
    source = release_dir / 'genome_results.tsv'
    _write_table(source, [_genome_row('AAFC_0001', 'AAFC', 'LRC159')])

    cli_core.main([
        'import-partner-workbook', '--release', str(release_dir), '--input', str(source),
        '--genome-only', '--dataset', 'Genome_Only', '-o', str(release_dir / 'registration'),
    ])

    assert Database(str(release_dir / 'project.sqlite')).get_all_ids() == {'AAFC_0001'}
    assert manager.get_release_manifest(release_dir)['workflow_status']['registration'] == 'COMPLETE'


def test_table_registration_does_not_reassign_existing_dataset(tmp_path):
    source = tmp_path / 'results.tsv'
    _write_table(source, [_genome_row('AAFC_0001', 'AAFC', 'G1')])
    db = Database(str(tmp_path / 'project.sqlite'))
    db.initialise()
    db.insert_sequences([('AAFC_0001', 'ACGT')], dataset='Existing')
    report = import_genome_only_partner_workbook(db, source, dataset='New')
    assert report[0]['result'] == 'CONFLICT'
    with db.connect() as conn:
        assert conn.execute("SELECT dataset, sequence FROM sequences WHERE id='AAFC_0001'").fetchone() == ('Existing', 'ACGT')


@pytest.mark.parametrize('content, message', [
    ('genome_id\tpartner_id\nG1\tAAFC\n', 'requires columns: sequence_id'),
    ('sequence_id\tpartner_id\nAAFC_0001\t\n', 'line 2'),
    ('sequence_id\tpartner_id\nAAFC_0001\tAAFC\nAAFC_0001\tOther\n', 'Conflicting partner_id'),
    ('sequence_id\tpartner_id\n', 'No genome registration rows'),
])
def test_invalid_tables_fail_explicitly(tmp_path, content, message):
    source = tmp_path / 'invalid.tsv'
    source.write_text(content)
    with pytest.raises(ValueError, match=message):
        load_genome_only_registration_table(source)


def test_existing_workbook_registration_remains_supported(tmp_path):
    source = tmp_path / 'registration.xlsx'
    workbook = Workbook()
    sheet = workbook.active
    sheet.title = 'Bulk Upload Template'
    sheet.append(['record_id', 'partner_id', 'partner_culture_id'])
    sheet.append([1, 'Agriculture and Agri-Food Canada (AAFC)', 'LRC159'])
    workbook.save(source)
    db = Database(str(tmp_path / 'project.sqlite'))
    db.initialise()
    report = import_genome_only_partner_workbook(db, source, dataset='AAFC')
    assert report[0]['sequence_id'] == 'AAFC_0001'
    assert report[0]['result'] == 'INSERTED'
    assert db.get_sequencing_metadata_for_ids()['AAFC_0001']['project_metadata']['partner_culture_id'] == 'LRC159'
