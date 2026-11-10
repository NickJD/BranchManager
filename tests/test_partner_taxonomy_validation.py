import tempfile
import unittest
from pathlib import Path

from branchmanager.cli import _validate_partner_taxonomy
from branchmanager.partner_metadata import load_partner_sequencing_metadata


class PartnerTaxonomyValidationTests(unittest.TestCase):
    def test_queries_without_classifier_rows_are_retained_as_unclassified(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / 'input.fasta'
            fasta.write_text('>QUB_001\nACGTACGT\n>QUB_002\nACGTACGA\n', encoding='utf-8')
            metadata = root / 'partner_metadata.tsv'
            metadata.write_text(
                'sequence_id\tpartner_id\talready_sequenced\ttaxon\n'
                'QUB_001\tQUB\tno\t\n'
                'QUB_002\tQUB\tno\t\n',
                encoding='utf-8',
            )
            # A classifier may omit queries for which it found no accepted hit.
            classification = root / 'taxonomy.tsv'
            classification.write_text(
                'id\tbest\tidentity\ttaxon\tconfidence\n', encoding='utf-8',
            )

            retained = _validate_partner_taxonomy(
                str(classification), str(metadata), 'bacteria', str(root),
                input_fasta=str(fasta),
            )

            self.assertEqual(retained, {'QUB_001', 'QUB_002'})
            report = (root / 'expected_taxonomy_validation.tsv').read_text(encoding='utf-8')
            self.assertEqual(report.count('UNCLASSIFIED_EXPECTED_BACTERIAL'), 2)
            self.assertEqual(report.count('\tYES\t'), 2)


    def test_generic_bacterial_taxon_is_treated_as_domain_expectation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            metadata = root / 'partner_metadata.tsv'
            metadata.write_text(
                'sequence_id\tpartner_id\talready_sequenced\ttaxon\n'
                'QUB_001\tQUB\tno\tbacterial\n',
                encoding='utf-8',
            )
            row = load_partner_sequencing_metadata(metadata)[0]
            self.assertEqual(row['expected_domain'], 'Bacteria')
            self.assertEqual(row['expected_taxon'], '')
            self.assertEqual(row['project_metadata']['partner_taxon_raw'], 'bacterial')

    def test_generic_bacterial_taxon_does_not_reject_specific_bacterial_assignment(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / 'input.fasta'
            fasta.write_text('>QUB_001\nACGTACGT\n', encoding='utf-8')
            metadata = root / 'partner_metadata.tsv'
            metadata.write_text(
                'sequence_id\tpartner_id\talready_sequenced\ttaxon\n'
                'QUB_001\tQUB\tno\tbacterial\n',
                encoding='utf-8',
            )
            classification = root / 'taxonomy.tsv'
            classification.write_text(
                'id\tbest\tidentity\ttaxon\tconfidence\n'
                'QUB_001\tref1\t99.9\td__Bacteria;p__Firmicutes;g__Winkia;s__Winkia anitrata\t1.0\n',
                encoding='utf-8',
            )
            retained = _validate_partner_taxonomy(
                str(classification), str(metadata), 'bacteria', str(root),
                input_fasta=str(fasta),
            )
            self.assertEqual(retained, {'QUB_001'})
            report = (root / 'expected_taxonomy_validation.tsv').read_text(encoding='utf-8')
            self.assertIn('ASSIGNED_REQUESTED_DOMAIN', report)
            self.assertIn('\tBacteria\tNA\t', report)

    def test_specific_partner_taxon_disagreement_is_reported_but_catalogue_record_is_retained(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / 'input.fasta'
            fasta.write_text('>QUB_001\nACGTACGT\n', encoding='utf-8')
            metadata = root / 'partner_metadata.tsv'
            metadata.write_text(
                'sequence_id\tpartner_id\talready_sequenced\ttaxon\n'
                'QUB_001\tQUB\tno\tEscherichia coli\n',
                encoding='utf-8',
            )
            classification = root / 'taxonomy.tsv'
            classification.write_text(
                'id\tbest\tidentity\ttaxon\tconfidence\n'
                'QUB_001\tref1\t99.9\td__Bacteria;p__Firmicutes;g__Winkia;s__Winkia anitrata\t1.0\n',
                encoding='utf-8',
            )
            retained = _validate_partner_taxonomy(
                str(classification), str(metadata), 'bacteria', str(root),
                input_fasta=str(fasta),
            )
            self.assertEqual(retained, {'QUB_001'})
            report = (root / 'expected_taxonomy_validation.tsv').read_text(encoding='utf-8')
            self.assertIn('TAXON_DISAGREEMENT', report)
            self.assertIn('catalogue value is retained as the gold-standard identity', report)



if __name__ == '__main__':
    unittest.main()


def test_classifier_domain_disagreement_is_reported_but_not_used_to_remove_catalogue_isolate(tmp_path):
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        fasta = root / 'input.fasta'
        fasta.write_text('>QUB_001\nACGTACGT\n', encoding='utf-8')
        metadata = root / 'Partner_Isolate_Catalogue_November.csv'
        metadata.write_text(
            'sequence_id,partner_id,already_sequenced,taxon\n'
            'QUB_001,QUB,no,Escherichia coli\n', encoding='utf-8',
        )
        classification = root / 'taxonomy.tsv'
        classification.write_text(
            'id\tbest\tidentity\ttaxon\tconfidence\n'
            'QUB_001\tref1\t99.9\td__Archaea;p__Euryarchaeota;g__Methanobrevibacter;s__Methanobrevibacter smithii\t1.0\n',
            encoding='utf-8',
        )
        retained = _validate_partner_taxonomy(
            str(classification), str(metadata), 'bacteria', str(root), input_fasta=str(fasta),
        )
        assert retained == {'QUB_001'}
        report = (root / 'expected_taxonomy_validation.tsv').read_text(encoding='utf-8')
        assert 'DOMAIN_DISAGREEMENT' in report
        assert 'Escherichia coli' in report


def test_partner_catalogue_taxon_is_preserved_as_raw_authoritative_metadata(tmp_path):
    catalogue = tmp_path / 'Partner_Isolate_Catalogue_November.csv'
    catalogue.write_text(
        'sequence_id,partner_id,already_sequenced,taxon,taxon_rank\n'
        'QUB_001,QUB,no,Escherichia coli,species\n', encoding='utf-8',
    )
    row = load_partner_sequencing_metadata(catalogue)[0]
    assert row['expected_taxon'] == 'Escherichia coli'
    assert row['project_metadata']['partner_taxon_raw'] == 'Escherichia coli'
    assert row['project_metadata']['partner_taxon_rank'] == 'species'
    assert row['project_metadata']['partner_taxonomy_source_file'] == 'Partner_Isolate_Catalogue_November.csv'


def test_combined_tree_taxonomy_uses_partner_catalogue_value_not_classifier(tmp_path):
    import json
    import sqlite3
    from branchmanager.pipeline.workflow_helpers import merge_combined_taxonomy_rows

    connection = sqlite3.connect(':memory:')
    connection.execute(
        'CREATE TABLE sequencing_metadata (id TEXT, project_metadata_json TEXT)'
    )
    connection.execute(
        'INSERT INTO sequencing_metadata VALUES (?, ?)',
        ('QUB_001', json.dumps({
            'partner_taxon_raw': 'Escherichia coli',
            'expected_taxon': 'Escherichia coli',
            'partner_taxonomy_source_file': 'Partner_Isolate_Catalogue_November.csv',
        })),
    )

    class FakeDB:
        def connect(self):
            return connection

    merged = merge_combined_taxonomy_rows(
        [('QUB_001', 'd__Bacteria;g__Wrong;s__Wrong species', '0.99')],
        '', {}, FakeDB(),
    )
    assert merged[0][1] == 'g__Escherichia;s__Escherichia coli'
    assert merged[0][2] == 'PARTNER_CATALOGUE_GOLD_STANDARD'


def test_november_catalogue_taxon_cannot_be_overwritten_by_other_metadata_imports():
    from branchmanager.partner_metadata import merge_partner_project_metadata

    existing = {
        'partner_taxon_raw': 'Escherichia coli',
        'expected_taxon': 'Escherichia coli',
        'partner_taxon_rank': 'species',
        'partner_taxonomy_source_file': 'Partner_Isolate_Catalogue_November.csv',
        'host_species': 'Bos taurus',
    }
    incoming = {
        'partner_taxon_raw': 'Lactobacillus acidophilus',
        'expected_taxon': 'Lactobacillus acidophilus',
        'partner_taxonomy_source_file': 'Partner_Metadata_October.csv',
        'host_breed': 'Holstein',
    }
    merged, protected = merge_partner_project_metadata(existing, incoming)
    assert protected is True
    assert merged['partner_taxon_raw'] == 'Escherichia coli'
    assert merged['expected_taxon'] == 'Escherichia coli'
    assert merged['partner_taxon_rank'] == 'species'
    assert merged['host_breed'] == 'Holstein'


def test_same_november_catalogue_can_explicitly_update_its_taxon():
    from branchmanager.partner_metadata import merge_partner_project_metadata

    existing = {
        'partner_taxon_raw': 'Escherichia coli',
        'partner_taxonomy_source_file': 'Partner_Isolate_Catalogue_November.csv',
    }
    incoming = {
        'partner_taxon_raw': 'Escherichia albertii',
        'expected_taxon': 'Escherichia albertii',
        'partner_taxon_rank': 'species',
        'partner_taxonomy_source_file': 'Partner_Isolate_Catalogue_November.csv',
    }
    merged, protected = merge_partner_project_metadata(existing, incoming)
    assert protected is True
    assert merged['partner_taxon_raw'] == 'Escherichia albertii'


def test_common_higher_rank_free_text_taxa_are_ranked_conservatively():
    from branchmanager.taxonomy import infer_taxon_rank

    assert infer_taxon_rank('Firmicutes') == 'p'
    assert infer_taxon_rank('Bacillota') == 'p'
    assert infer_taxon_rank('Lachnospiraceae') == 'f'
    assert infer_taxon_rank('Clostridiales') == 'o'
    assert infer_taxon_rank('Clostridia') == 'c'
    assert infer_taxon_rank('Streptococcus') == 'g'
    assert infer_taxon_rank('Streptococcus equinus') == 's'
    assert infer_taxon_rank('Firmicutes', 'genus') == 'g'


def test_bacterial_taxon_catalogue_column_alias_is_recognised(tmp_path):
    catalogue = tmp_path / 'Partner_Isolate_Catalogue_November.csv'
    catalogue.write_text(
        'sequence_id,partner_id,already_sequenced,taxon (bacteria)\n'
        'QUB_001,QUB,no,Escherichia coli\n', encoding='utf-8',
    )
    row = load_partner_sequencing_metadata(catalogue)[0]
    assert row['expected_taxon'] == 'Escherichia coli'
    assert row['project_metadata']['partner_taxon_raw'] == 'Escherichia coli'
