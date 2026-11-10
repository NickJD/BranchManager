"""Tests for preserving saved decisions during Assistant cohort PREPARE."""
from __future__ import annotations

import ast
import csv
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "src/branchmanager/cli_core.py"


def load_restore_helper():
    tree = ast.parse(CLI.read_text(encoding="utf-8"))
    wanted = {"_cohort_read_rows", "_cohort_read_fallback_option_ids", "_cohort_manual_read_selections", "_cohort_restore_review_state"}
    nodes = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in wanted]
    assert {node.name for node in nodes} == wanted
    namespace = {"Path": Path, "csv": csv, "json": __import__("json")}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(CLI), "exec"), namespace)
    return namespace["_cohort_restore_review_state"]


class RestoreDecisionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.restore = staticmethod(load_restore_helper())

    def _write_tsv(self, path: Path, row: dict):
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
            writer.writeheader()
            writer.writerow(row)

    def test_saved_decisions_file_takes_priority_over_old_template(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            template = root / "cohort_marker_review_template.tsv"
            decisions = root / "cohort_marker_review_decisions.tsv"
            self._write_tsv(template, {
                "sequence_id": "QUB_0327", "decision": "REJECT", "reviewer": "old", "reviewed_at": "2026-11-01", "notes": "old note"
            })
            self._write_tsv(decisions, {
                "sequence_id": "QUB_0327", "decision": "RESEQUENCE", "reviewer": "nickd", "reviewed_at": "2026-11-02", "notes": "re-run trace"
            })
            rows = [{"sequence_id": "QUB_0327", "qc_class": "FAIL_QC", "read_level_evidence": "new evidence"}]
            restored = self.restore(rows, template, decisions)
            self.assertEqual(restored, 1)
            self.assertEqual(rows[0]["decision"], "RESEQUENCE")
            self.assertEqual(rows[0]["reviewer"], "nickd")
            self.assertEqual(rows[0]["notes"], "re-run trace")
            self.assertEqual(rows[0]["read_level_evidence"], "new evidence")

    def test_decisions_file_can_explicitly_clear_an_old_decision(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            template = root / "cohort_marker_review_template.tsv"
            decisions = root / "cohort_marker_review_decisions.tsv"
            self._write_tsv(template, {
                "sequence_id": "QUB_0327", "decision": "REJECT", "reviewer": "old", "reviewed_at": "2026-11-01", "notes": "old note"
            })
            self._write_tsv(decisions, {
                "sequence_id": "QUB_0327", "decision": "", "reviewer": "nickd", "reviewed_at": "", "notes": ""
            })
            rows = [{"sequence_id": "QUB_0327"}]
            self.restore(rows, template, decisions)
            self.assertEqual(rows[0]["decision"], "")
            self.assertEqual(rows[0]["reviewer"], "nickd")

    def test_previous_template_is_fallback_when_no_decisions_file_given(self):
        with tempfile.TemporaryDirectory() as directory:
            template = Path(directory) / "cohort_marker_review_template.tsv"
            self._write_tsv(template, {
                "sequence_id": "QUB_0327", "decision": "RESEQUENCE", "reviewer": "nickd", "reviewed_at": "2026-11-02", "notes": "repeat"
            })
            rows = [{"sequence_id": "QUB_0327"}]
            self.restore(rows, template)
            self.assertEqual(rows[0]["decision"], "RESEQUENCE")


if __name__ == "__main__":
    unittest.main()

    def test_cli_extracts_manual_fallback_choices_from_studio_decisions(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "cohort_marker_review_decisions.tsv"
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=[
                    "sequence_id", "selected_read_id", "read_selection_source", "decision"
                ], delimiter="\t")
                writer.writeheader()
                writer.writerows([
                    {"sequence_id": "BSI_0016", "selected_read_id": "BSI_0016_1492R", "read_selection_source": "MANUAL", "decision": ""},
                    {"sequence_id": "BSI_0017", "selected_read_id": "BSI_0017_27F", "read_selection_source": "AUTO_LONGEST_USABLE_READ", "decision": ""},
                    {"sequence_id": "BSI_0018", "selected_read_id": "", "read_selection_source": "", "decision": "RESEQUENCE"},
                ])
            tree = ast.parse(CLI.read_text(encoding="utf-8"))
            wanted = {"_cohort_read_rows", "_cohort_manual_read_selections"}
            nodes = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in wanted]
            assert {node.name for node in nodes} == wanted
            local = {"Path": Path, "csv": csv}
            exec(compile(ast.Module(body=nodes, type_ignores=[]), str(CLI), "exec"), local)
            self.assertEqual(local["_cohort_manual_read_selections"](path), {
                "BSI_0016": "BSI_0016_1492R"
            })

    def test_manual_read_selection_restored_only_if_candidate_still_exists(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            decisions = root / 'cohort_marker_review_decisions.tsv'
            self._write_tsv(decisions, {
                'sequence_id': 'BSI_0016', 'decision': '', 'reviewer': 'nickd',
                'reviewed_at': '2026-11-03', 'notes': '',
                'selected_read_id': 'BSI_0016_27F', 'read_selection_source': 'MANUAL',
            })
            rows = [{
                'sequence_id': 'BSI_0016', 'fallback_read_options': '[{"read_id":"BSI_0016_27F"},{"read_id":"BSI_0016_1492R"}]',
                'recommended_read_id': 'BSI_0016_1492R',
            }]
            self.restore(rows, root / 'no_previous_template.tsv', decisions)
            self.assertEqual(rows[0]['selected_read_id'], 'BSI_0016_27F')
            self.assertEqual(rows[0]['read_selection_source'], 'MANUAL')

            rows = [{
                'sequence_id': 'BSI_0016', 'fallback_read_options': '[{"read_id":"BSI_0016_1492R"}]',
                'recommended_read_id': 'BSI_0016_1492R',
            }]
            self.restore(rows, root / 'no_previous_template.tsv', decisions)
            self.assertEqual(rows[0]['selected_read_id'], 'BSI_0016_1492R')
            self.assertEqual(rows[0]['read_selection_source'], 'AUTO_LONGEST_USABLE_READ')
