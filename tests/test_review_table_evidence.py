"""Focused tests for the cohort-review read-level evidence display."""
from __future__ import annotations

import ast
import unittest
from pathlib import Path
from types import SimpleNamespace
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
PAPER_TRAIL = ROOT / "src/branchmanager/pipeline/paper_trail.py"


def load_review_evidence_helpers():
    tree = ast.parse(PAPER_TRAIL.read_text(encoding="utf-8"))
    wanted = {"_n_percent", "_read_level_review_evidence"}
    nodes = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in wanted]
    assert {node.name for node in nodes} == wanted
    module = ast.Module(body=nodes, type_ignores=[])
    namespace = {"Iterable": Iterable}
    exec(compile(module, str(PAPER_TRAIL), "exec"), namespace)
    return namespace["_read_level_review_evidence"]


class ReviewEvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.format_evidence = staticmethod(load_review_evidence_helpers())

    def test_empty_consensus_can_show_failed_read_metrics(self):
        read = SimpleNamespace(
            read_id="JDP771_61777713_61777713",
            primer="27F",
            direction="forward",
            qc_class="FAIL_QC",
            status="filtered",
            trimmed_sequence="A" * 987,
            masked_sequence=("N" * 86) + ("A" * (987 - 86)),
            warning="n_percent_gt_5;mixed_peak_percent_8.32",
        )
        evidence = self.format_evidence([read])
        self.assertIn("trimmed 987 bp", evidence)
        self.assertIn("masked N 8.713%", evidence)
        self.assertIn("FAIL_QC", evidence)
        self.assertIn("n_percent_gt_5", evidence)

    def test_multiple_reads_are_all_preserved(self):
        reads = [
            SimpleNamespace(read_id="FWD", primer="27F", direction="forward", qc_class="FAIL_QC",
                            status="filtered", trimmed_sequence="A" * 987,
                            masked_sequence="N" * 86 + "A" * 901, warning="n_percent_gt_5"),
            SimpleNamespace(read_id="REV", primer="1492R", direction="reverse", qc_class="PASS_HIGH_CONFIDENCE",
                            status="kept", trimmed_sequence="C" * 961,
                            masked_sequence="C" * 961, warning=""),
        ]
        evidence = self.format_evidence(reads)
        self.assertIn("FWD (27F/forward)", evidence)
        self.assertIn("REV (1492R/reverse)", evidence)
        self.assertIn("masked N 0.000%", evidence)

    def test_no_read_group_returns_empty_evidence(self):
        self.assertEqual(self.format_evidence([]), "")


if __name__ == "__main__":
    unittest.main()
