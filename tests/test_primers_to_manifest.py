import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from branchmanager.tools.primers_to_manifest import process_files


class PrimersToManifestTests(unittest.TestCase):
    def test_warns_when_ab1_file_is_assigned_to_different_isolates(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / 'reads'
            source.mkdir()
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
            self.assertIn('ISO1,shared.ab1,forward,27F', manifest)
            self.assertIn('ISO2,shared.ab1,forward,27F', manifest)


if __name__ == '__main__':
    unittest.main()
