import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from branchmanager.tools.primers_to_manifest import process_files, _infer_primer_from_filename


class PrimersToManifestTests(unittest.TestCase):
    def test_warns_when_ab1_file_is_assigned_to_different_isolates(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / 'reads'
            source.mkdir()
            (root / 'mapping.xlsx').touch()
            (source / 'shared.ab1').write_bytes(b'ABIF-placeholder')
            rows = pd.DataFrame({
                'Partner culture collection isolate ID': ['ISO1', 'ISO2'],
                'Filename of primer sequence (27F)': ['shared.ab1', 'shared.ab1'],
            })
            output = io.StringIO()

            with patch('pandas.read_excel', return_value=rows), contextlib.redirect_stdout(output):
                process_files(root / 'mapping.xlsx', source, root / 'output', 'manifest.csv')

            self.assertIn('AB1 file', output.getvalue())
            self.assertIn('assigned to multiple isolates: ISO1, ISO2', output.getvalue())
            manifest = (root / 'output' / 'manifest.csv').read_text()
            self.assertIn('ISO1,ISO1_27F.ab1,forward,27F', manifest)
            self.assertIn('ISO2,ISO2_27F.ab1,forward,27F', manifest)


    def test_preserves_resequencing_reads_from_other_columns_as_unique_attempts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / 'reads'
            (source / 'submission_1').mkdir(parents=True)
            (source / 'submission_2').mkdir(parents=True)
            (source / 'submission_1' / 'ss172_907R.ab1').write_bytes(b'ABIF-one')
            (source / 'submission_1' / 'ss172_785F.ab1').write_bytes(b'ABIF-two')
            (source / 'submission_2' / 'ss172_907R.ab1').write_bytes(b'ABIF-three')
            mapping = root / 'mapping.xlsx'
            mapping.touch()
            rows = pd.DataFrame({
                'Partner culture collection isolate ID': ['ss172'],
                'Filename of primer sequence (Other #1)': ['submission_1/ss172_907R'],
                'Filename of primer sequence (Other #2)': ['submission_1/ss172_785F'],
                'Filename of primer sequence (Other #4)': ['submission_2/ss172_907R'],
            })
            with patch('pandas.read_excel', return_value=rows):
                result = process_files(mapping, source, root / 'out', 'manifest.csv')
            self.assertEqual(result, 0)

            import csv
            with (root / 'out' / 'manifest.csv').open(newline='') as handle:
                manifest = list(csv.DictReader(handle))
            self.assertEqual(len(manifest), 3)
            reverse_reads = [row for row in manifest if row['primer'] == '907R']
            self.assertEqual(len(reverse_reads), 2)
            self.assertEqual(
                {row['read_file'] for row in reverse_reads},
                {'ss172_907R_attempt01.ab1', 'ss172_907R_attempt02.ab1'},
            )
            self.assertEqual({row['direction'] for row in reverse_reads}, {'reverse'})
            self.assertEqual({row['sequence_id'] for row in manifest}, {'ss172'})
            self.assertTrue(all(row['source_sha256'] for row in manifest))
            self.assertEqual(
                {row['attempt_group'] for row in reverse_reads},
                {'submission_1', 'submission_2'},
            )
            self.assertTrue((root / 'out' / 'ss172_907R_attempt01.ab1').is_file())
            self.assertTrue((root / 'out' / 'ss172_907R_attempt02.ab1').is_file())

    def test_resolves_workbook_relative_path_without_ab1_extension(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / 'reads'
            (source / 'EC00170083').mkdir(parents=True)
            trace = source / 'EC00170083' / 'ss172_907R.ab1'
            trace.write_bytes(b'ABIF-placeholder')
            mapping = root / 'mapping.xlsx'
            mapping.touch()
            rows = pd.DataFrame({
                'Partner culture collection isolate ID': ['ss172'],
                'Filename of primer sequence (Other #1)': ['EC00170083/ss172_907R'],
            })
            with patch('pandas.read_excel', return_value=rows):
                result = process_files(mapping, source, root / 'out', 'manifest.csv')
            self.assertEqual(result, 0)
            self.assertTrue((root / 'out' / 'ss172_907R.ab1').is_file())



    def test_resolves_submission_folder_suffix_when_export_adds_date_prefix(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / 'reads'
            actual_folder = source / '2026-11-12 - Samantha Joan Noel EC00170083'
            actual_folder.mkdir(parents=True)
            (actual_folder / 'ss172_907R.ab1').write_bytes(b'ABIF-placeholder')
            mapping = root / 'mapping.xlsx'
            mapping.touch()
            rows = pd.DataFrame({
                'Partner culture collection isolate ID': ['ss172'],
                'Filename of primer sequence (Other #1)': ['EC00170083/ss172_907R'],
            })
            with patch('pandas.read_excel', return_value=rows):
                result = process_files(mapping, source, root / 'out', 'manifest.csv')
            self.assertEqual(result, 0)
            import csv
            with (root / 'out' / 'manifest.csv').open(newline='') as handle:
                row = next(csv.DictReader(handle))
            self.assertEqual(row['primer'], '907R')
            self.assertEqual(row['attempt_group'], actual_folder.name)



    def test_missing_qualified_reference_does_not_fall_back_to_other_submission(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / 'reads'
            (source / 'EC00124257_last_ss').mkdir(parents=True)
            actual = source / 'EC00124257_last_ss' / 'ss194_907R.ab1'
            actual.write_bytes(b'ABIF-placeholder')
            mapping = root / 'mapping.xlsx'
            mapping.touch()
            rows = pd.DataFrame({
                'Partner culture collection isolate ID': ['ss194'],
                'Filename of primer sequence (Other #1)': ['EC00170083/ss194_907R'],
                'Filename of primer sequence (Other #2)': ['EC00124257_last_ss/ss194_907R'],
            })
            output = io.StringIO()
            with patch('pandas.read_excel', return_value=rows), contextlib.redirect_stdout(output):
                result = process_files(mapping, source, root / 'out', 'manifest.csv')
            self.assertEqual(result, 0)
            self.assertIn('EC00170083/ss194_907R', output.getvalue())
            self.assertNotIn('duplicate mapping reference', output.getvalue())
            import csv
            with (root / 'out' / 'manifest.csv').open(newline='') as handle:
                manifest = list(csv.DictReader(handle))
            self.assertEqual(len(manifest), 1)
            self.assertEqual(manifest[0]['attempt_group'], 'EC00124257_last_ss')

    def test_writes_console_output_to_default_log_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / 'reads'
            source.mkdir()
            (source / 'sample_907-R.ab1').write_bytes(b'ABIF-placeholder')
            mapping = root / 'mapping.xlsx'
            mapping.touch()
            rows = pd.DataFrame({
                'Partner culture collection isolate ID': ['ISO1'],
                'Filename of primer sequence (Other #1)': ['sample_907-R'],
            })
            with patch('pandas.read_excel', return_value=rows):
                result = process_files(mapping, source, root / 'out', 'manifest.csv')
            self.assertEqual(result, 0)
            log = (root / 'out' / 'primers_to_manifest.log').read_text()
            self.assertIn('Reading Excel:', log)
            self.assertIn('Manifest created at:', log)
            self.assertIn('Run finished with exit status: 0', log)
            self.assertTrue((root / 'out' / 'ISO1_907R.ab1').is_file())

    def test_infers_additional_primer_notations_from_filenames(self):
        cases = {
            'ss1_907R.ab1': '907R',
            'ss1_785F.ab1': '785F',
            'ss1_1492R.ab1': '1492R',
            'ss1_907-R.ab1': '907R',
            'ss1_1492_R.ab1': '1492R',
            'ss1_ITS-1.ab1': 'ITS1',
            'ss1_18S.ab1': '18S',
            'ss1_DS-2.ab1': 'DS2',
        }
        for filename, expected in cases.items():
            with self.subTest(filename=filename):
                self.assertEqual(_infer_primer_from_filename(Path(filename)), expected)

    def test_unqualified_basename_is_ambiguous_across_submissions(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / 'reads'
            for folder in ('submission_1', 'submission_2'):
                (source / folder).mkdir(parents=True)
                (source / folder / 'ss172_907R.ab1').write_bytes(folder.encode())
            mapping = root / 'mapping.xlsx'
            mapping.touch()
            rows = pd.DataFrame({
                'Partner culture collection isolate ID': ['ss172'],
                'Filename of primer sequence (Other #1)': ['ss172_907R'],
            })
            output = io.StringIO()
            with patch('pandas.read_excel', return_value=rows), contextlib.redirect_stdout(output):
                process_files(mapping, source, root / 'out', 'manifest.csv')
            self.assertIn('multiple source files match unqualified reference', output.getvalue())
            self.assertFalse((root / 'out' / 'manifest.csv').exists())


if __name__ == '__main__':
    unittest.main()
