import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'src'
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from branchmanager.partner_metadata import load_partner_workbook_rows, load_partner_sequencing_metadata
from branchmanager.db.interface import Database
from branchmanager.project_state import import_partner_catalogue


def test_partner_survey_csv_headers_are_supported(tmp_path):
    csv_path = tmp_path / 'survey_export.csv'
    with csv_path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            'Sequential identifier:', 'Partner ID:', 'Partner culture collection isolate ID:',
            'Microorganism type:',
            'Taxonomic classification method: (choice=Amplicon)',
            'Taxonomic classification method: (choice=MALDI-TOF)',
            'NCBI scientific name of the isolate from amplicon taxonomic identification:',
            'Percentage identity of the amplicon assignment (NCBI):',
            'Taxon name/ID returned by MALDI-TOF MS:',
            'What sequencing platform was used for amplicon sequencing?',
        ])
        writer.writeheader()
        writer.writerow({
            'Sequential identifier:': 'AU_0001',
            'Partner ID:': 'Aarhus University (AU)',
            'Partner culture collection isolate ID:': 'ss172',
            'Microorganism type:': 'Bacteria',
            'Taxonomic classification method: (choice=Amplicon)': 'Checked',
            'Taxonomic classification method: (choice=MALDI-TOF)': 'Unchecked',
            'NCBI scientific name of the isolate from amplicon taxonomic identification:': 'Segatella paludivivens',
            'Percentage identity of the amplicon assignment (NCBI):': '90.73',
            'Taxon name/ID returned by MALDI-TOF MS:': '',
            'What sequencing platform was used for amplicon sequencing?': 'Sanger',
        })
    rows = load_partner_workbook_rows(csv_path, sheet_name='ignored_for_csv')
    assert len(rows) == 1
    row = rows[0]
    assert row['source_id'] == 'AU_0001'
    assert row['partner_id'] == 'AU'
    assert row['amplicon_method_reported'] is True
    assert row['maldi_method_reported'] is False
    assert row['project_metadata']['partner_culture_id'] == 'ss172'
    assert row['project_metadata']['amplicon_ncbi_name_reported'] == 'Segatella paludivivens'
    assert row['project_metadata']['amplicon_ncbi_identity_reported'] == '90.73'
    assert row['project_metadata']['partner_csv_raw_fields']['What sequencing platform was used for amplicon sequencing?'] == 'Sanger'
    assert 'Add the full or partial 16S/18S/ITS sequence.' not in row['project_metadata']['partner_csv_raw_fields']


def test_partner_catalogue_loader_retains_partner_own_id(tmp_path):
    catalogue = tmp_path / 'Partner_Isolate_Catalogue_November.csv'
    with catalogue.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            'sequence_id', 'partner_id', 'selected_for_genome_sequencing',
            'already_sequenced', 'addition_method', 'taxon', 'partner_own_ids',
        ])
        writer.writeheader()
        writer.writerow({
            'sequence_id': 'AU_0001', 'partner_id': 'AU',
            'selected_for_genome_sequencing': 'no', 'already_sequenced': 'yes',
            'addition_method': 'marker_sequence', 'taxon': 'Bacteria',
            'partner_own_ids': 'ss172',
        })
    rows = load_partner_sequencing_metadata(catalogue)
    assert rows[0]['project_metadata']['partner_own_id'] == 'ss172'


def test_catalogue_import_registers_new_records_and_is_idempotent(tmp_path):
    catalogue = tmp_path / 'catalogue.csv'
    with catalogue.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            'sequence_id', 'partner_id', 'selected_for_genome_sequencing',
            'already_sequenced', 'addition_method', 'taxon', 'partner_own_ids',
        ])
        writer.writeheader()
        for sid, own_id in [('AU_0001', 'ss172'), ('AU_0002', 'ss174')]:
            writer.writerow({
                'sequence_id': sid, 'partner_id': 'AU',
                'selected_for_genome_sequencing': 'no', 'already_sequenced': 'yes',
                'addition_method': 'marker_sequence', 'taxon': 'Bacteria',
                'partner_own_ids': own_id,
            })
    db = Database(str(tmp_path / 'project.sqlite'))
    db.initialise()
    first = import_partner_catalogue(db, catalogue)
    second = import_partner_catalogue(db, catalogue)
    assert sum(r['result'] == 'INSERTED' for r in first) == 2
    assert sum(r['result'] == 'UPDATED' for r in second) == 2
    assert set(db.get_all_ids()) == {'AU_0001', 'AU_0002'}
    metadata = db.get_sequencing_metadata_for_ids()
    assert metadata['AU_0001']['project_metadata']['partner_own_id'] == 'ss172'
    assert metadata['AU_0001']['selected_for_wgs'] is True


def test_annual_report_genome_ledger_displays_both_ids(tmp_path):
    from branchmanager.reporting import write_annual_report

    db = Database(str(tmp_path / 'project.sqlite'))
    db.initialise()
    db.insert_sequences([('AU_0001', '')], dataset='AU_01')
    db.upsert_sequencing_metadata([{
        'id': 'AU_0001', 'partner_id': 'AU', 'dataset': 'AU_01',
        'selected_for_sequencing': False, 'selected_for_wgs': True,
        'source_id': 'AU_0001', 'source_file': 'test catalogue',
        'project_metadata': {'partner_own_id': 'ss172'},
    }])
    db.upsert_genome_records([{
        'genome_id': 'GCA_000001.1', 'sequence_id': 'AU_0001',
        'genome_status': 'GENOME_QC_PASSED', 'genome_qc_pass': True,
        'completeness': 99.0, 'contamination': 0.5,
        'gtdb_taxonomy': 'd__Bacteria;p__Bacillota;g__Example',
    }])
    output_dir = tmp_path / 'annual'
    write_annual_report(db, output_dir)
    ledger = output_dir / 'project_genome_ledger.tsv'
    with ledger.open(encoding='utf-8') as handle:
        header = handle.readline().rstrip('\n').split('\t')
        row = handle.readline().rstrip('\n').split('\t')
    assert 'partner_own_id' in header
    assert row[header.index('sequence_id')] == 'AU_0001'
    assert row[header.index('partner_own_id')] == 'ss172'


def test_dashboard_database_reader_exposes_partner_own_id(tmp_path):
    from branchmanager.reporting import _read_dashboard_database

    db = Database(str(tmp_path / 'project.sqlite'))
    db.initialise()
    db.insert_sequences([('AU_0001', '')], dataset='AU_01')
    db.upsert_sequencing_metadata([{
        'id': 'AU_0001', 'partner_id': 'AU', 'dataset': 'AU_01',
        'selected_for_sequencing': False, 'selected_for_wgs': False,
        'source_id': 'AU_0001', 'source_file': 'test catalogue',
        'project_metadata': {'partner_own_id': 'ss172'},
    }])
    states, genomes, identity, registered = _read_dashboard_database(db.path)
    assert states['AU_0001']['partner_own_id'] == 'ss172'
    assert registered[0]['partner_own_id'] == 'ss172'
    assert genomes == []
