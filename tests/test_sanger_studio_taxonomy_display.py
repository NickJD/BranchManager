"""Focused tests for taxonomy evidence exposed by Sanger Studio inspections."""
from __future__ import annotations

import ast
import json
import tempfile
import unittest
import urllib.parse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BACKEND = ROOT / "src/branchmanager/sanger_studio.py"
HTML = ROOT / "src/branchmanager/data/sanger_studio.html"


def load_backend_helpers(names):
    tree = ast.parse(BACKEND.read_text(encoding="utf-8"))
    nodes = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in names]
    assert {node.name for node in nodes} == set(names)
    namespace = {"Path": Path, "json": json, "urllib": urllib, "re": __import__("re"), "csv": __import__("csv")}
    module = ast.Module(body=nodes, type_ignores=[])
    exec(compile(module, str(BACKEND), "exec"), namespace)
    return namespace


class StudioTaxonomyDisplayTests(unittest.TestCase):
    def test_inspection_serialises_per_read_summary_and_candidate_hits(self):
        names = {"_is_within", "_normalise_key", "_table_delimiter", "_read_tsv", "_row_get", "_read_report_row", "_parse_fallback_read_options", "_serialise_inspection"}
        ns = load_backend_helpers(names)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            job = root / "inspections" / "job1"
            screen = job / "read_taxonomy_screen"
            screen.mkdir(parents=True)
            (job / "assembly.tsv").write_text(
                "SequenceID\tStatus\tQCClass\tRecommendation\tReasons\tOutputLength\tOverlapLength\tOverlapIdentity\tUsedReads\tReadCount\tTaxonomyConcordance\tReadTaxonomyAssignments\tSuggestedAction\n"
                "QUB_0289\tassembled\tPASS_WITH_WARNINGS\tMANUAL_REVIEW\tprimer_read_taxonomic_conflict_at_genus\t1400\t120\t99.2\t2\t2\tPOTENTIAL_CONFLICT\tJTQ765=g__Sarcina | JTQ766=g__Veillonella\tReview taxonomy\n",
                encoding="utf-8",
            )
            (job / "read_qc.tsv").write_text(
                "ReadID\tSequenceID\tPrimer\tDirection\tRawLength\tTrimmedLength\tMeanRawQuality\tMeanTrimmedQuality\tRawExpectedErrors\tTrimmedExpectedErrors\tMaskedNPercent\tQCClass\tStatus\tReasons\n"
                "JTQ765\tQUB_0289\t27F\tforward\t1450\t1400\t40\t42\t2\t1\t0.1\tPASS_HIGH_CONFIDENCE\tkept\t\n"
                "JTQ766\tQUB_0289\t1492R\treverse\t1450\t1390\t40\t41\t2\t1\t0.2\tPASS_HIGH_CONFIDENCE\tkept\t\n"
                "OTHER\tQUB_0001\t27F\tforward\t1000\t900\t40\t40\t2\t1\t0.0\tPASS_HIGH_CONFIDENCE\tkept\t\n",
                encoding="utf-8",
            )
            (job / "assembly_read_placements.tsv").write_text("SequenceID\tReadID\nQUB_0289\tJTQ765\n", encoding="utf-8")
            (screen / "taxonomy.tsv").write_text(
                "ID\tBestHit\tIdentity\tQueryCoverage\tTargetCoverage\tQualifiedHitCount\tNearBestHitCount\tAssignedRank\tTaxon\tReadClassificationStatus\tTaxonomySupport\n"
                "JTQ765\trefA\t99.8\t99.0\t99.0\t3\t2\tg\tg__Sarcina\tRESOLVED\tMULTIPLE_NEAR_BEST_HITS_AGREE_AT_GENUS\n"
                "JTQ766\trefB\t99.7\t98.0\t98.0\t2\t1\tg\tg__Veillonella\tRESOLVED\tSINGLE_NEAR_BEST_REFERENCE_HIT\n"
                "OTHER\trefC\t99.0\t99.0\t99.0\t1\t1\tg\tg__Other\tRESOLVED\tSINGLE_NEAR_BEST_REFERENCE_HIT\n",
                encoding="utf-8",
            )
            (screen / "read_taxonomy_hits.tsv").write_text(
                "ReadID\tTargetID\tIdentity\tQueryCoverage\tTargetCoverage\tAlignmentLength\tTaxonomy\tQualified\tNearBest\tIdentityDelta\n"
                "JTQ765\trefA\t99.8\t99.0\t99.0\t1400\tg__Sarcina\tyes\tyes\t0.0\n"
                "JTQ765\trefA2\t99.6\t98.0\t98.0\t1390\tg__Sarcina\tyes\tyes\t0.2\n"
                "JTQ766\trefB\t99.7\t98.0\t98.0\t1380\tg__Veillonella\tyes\tyes\t0.0\n"
                "OTHER\trefC\t99.0\t99.0\t99.0\t900\tg__Other\tyes\tyes\t0.0\n",
                encoding="utf-8",
            )
            output = {
                "assembly_tsv": str(job / "assembly.tsv"),
                "read_qc_tsv": str(job / "read_qc.tsv"),
                "assembly_placements_tsv": str(job / "assembly_read_placements.tsv"),
                "assembly_pngs": [], "chromatogram_pngs": [], "read_error_pngs": [],
                "taxonomy_screen_tsv": str(job / "read_taxonomy_concordance.tsv"),
                "summary": str(job / "summary.txt"),
            }
            (job / "read_taxonomy_concordance.tsv").write_text("SequenceID\tConcordanceStatus\nQUB_0289\tPOTENTIAL_CONFLICT\n", encoding="utf-8")
            (job / "summary.txt").write_text("summary\n", encoding="utf-8")
            result = ns["_serialise_inspection"](output, "QUB_0289", job, root)
            self.assertTrue(result["taxonomy_screen_requested"])
            self.assertEqual([row["ID"] for row in result["taxonomy_reads"]], ["JTQ765", "JTQ766"])
            self.assertEqual(len(result["taxonomy_hits"]), 3)
            self.assertNotIn("OTHER", {row["ReadID"] for row in result["taxonomy_hits"]})
            labels = {item["label"] for item in result["files"]}
            self.assertIn("Per-read taxonomy reference hits", labels)

    def test_tsv_with_json_commas_keeps_tab_delimiter_and_serialises_assembly_row(self):
        ns = load_backend_helpers({"_table_delimiter", "_read_tsv", "_normalise_key", "_row_get", "_read_report_row"})
        with tempfile.TemporaryDirectory() as directory:
            report = Path(directory) / "assembly_report.tsv"
            report.write_text(
                "SequenceID\tStatus\tQCClass\tRecommendation\tOutputLength\tUsedReads\tReadCount\tSelectedReadID\tFallbackReadOptions\n"
                'BSI_0079\tpartial_no_overlap\tFAIL_QC\tRESEQUENCE\t861\t1\t2\tBSI_0079_1492R\t'
                '[{"read_id":"BSI_0079_1492R","length":861,"mean_quality":40.67},'
                '{"read_id":"BSI_0079_27F","length":316,"mean_quality":49.20}]\n',
                encoding="utf-8",
            )
            rows = ns["_read_tsv"](report)
            self.assertEqual(len(rows), 1)
            self.assertIn("FallbackReadOptions", rows[0])
            row = ns["_read_report_row"](report, "BSI_0079")
            self.assertEqual(row["Status"], "partial_no_overlap")
            self.assertEqual(row["QCClass"], "FAIL_QC")
            self.assertEqual(row["OutputLength"], "861")
            self.assertEqual(row["SelectedReadID"], "BSI_0079_1492R")
            self.assertEqual(len(json.loads(row["FallbackReadOptions"])), 2)

    def test_paper_trail_report_writers_keep_selection_columns_aligned(self):
        source = (ROOT / "src/branchmanager/pipeline/paper_trail.py").read_text(encoding="utf-8")
        assert r"SelectionSource\tFallbackReason\tSelectionWarning\tRecommendedReadID\tFallbackReadOptions" in source
        assert 'stats.get("selection_warning", "")' in source
        assert r"SelectionWarning\tRecommendedReadID\tFallbackReadOptions\tAssemblyStatus\tMergeMethod" in source

    def test_legacy_assembly_report_with_missing_selection_warning_is_repaired(self):
        ns = load_backend_helpers({"_table_delimiter", "_read_tsv", "_normalise_key", "_row_get", "_read_report_row"})
        with tempfile.TemporaryDirectory() as directory:
            report = Path(directory) / "assembly_report.tsv"
            header = (
                "SequenceID\tStatus\tReadCount\tUsedReads\tOutputLength\tMeanQuality\t"
                "OutputExpectedErrors\tOutputNPercent\tConflictDensity\tQCClass\tRecommendation\t"
                "TaxonomyConcordance\tTaxonomyComparedRank\tReadTaxonomyAssignments\tOverlapIdentity\t"
                "OverlapLength\tConflicts\tAmbiguousConflicts\tUnmergedReads\tReadIDs\tKeptReadIDs\t"
                "UsedReadIDs\tMergeMethod\tProcessingMode\tSelectedReadID\tSelectionSource\tFallbackReason\t"
                "SelectionWarning\tRecommendedReadID\tFallbackReadOptions\tAssemblyStatus\tPassLengthQC\t"
                "Reasons\tSuggestedAction\tFailedReadIDs\tModeWarning\tQualityStatus\n"
            )
            # The legacy row omits SelectionWarning, so fields after FallbackReason are one column early.
            row = (
                "BSI_0079\tpartial_no_overlap\t2\t1\t861\t40.67\t0.4532\t2.207\t0.0\tFAIL_QC\tRESEQUENCE\t"
                "NOT_SCREENED\tNA\tNA\tNA\tNA\t0\t0\tBSI_0079_27F\tBSI_0079_1492R;BSI_0079_27F\t"
                "BSI_0079_1492R;BSI_0079_27F\tBSI_0079_1492R\tlongest_read_no_overlap_fallback\tassemble\t"
                "BSI_0079_1492R\tAUTO_LONGEST_USABLE_READ\tno_qualifying_overlap\tBSI_0079_1492R\t"
                '[{"read_id":"BSI_0079_1492R","length":861},{"read_id":"BSI_0079_27F","length":316}]\t'
                "partial_no_overlap\tyes\tpartial_no_overlap;unmerged_primer_reads_1\t"
                "Repeat Sanger sequencing for this isolate.\t\t\tPHRED_VERIFIED\n"
            )
            report.write_text(header + row, encoding="utf-8")
            repaired = ns["_read_report_row"](report, "BSI_0079")
            assert repaired['SelectionWarning'] == ''
            assert repaired['RecommendedReadID'] == 'BSI_0079_1492R'
            assert json.loads(repaired['FallbackReadOptions'])[0]['length'] == 861
            assert repaired['AssemblyStatus'] == 'partial_no_overlap'
            assert repaired['PassLengthQC'] == 'yes'
            assert repaired['Reasons'] == 'partial_no_overlap;unmerged_primer_reads_1'
            assert repaired['SuggestedAction'] == 'Repeat Sanger sequencing for this isolate.'

    def test_report_row_missing_is_a_clear_error_not_an_empty_result(self):
        ns = load_backend_helpers({"_table_delimiter", "_read_tsv", "_normalise_key", "_row_get", "_read_report_row"})
        with tempfile.TemporaryDirectory() as directory:
            report = Path(directory) / "assembly_report.tsv"
            report.write_text("SequenceID\tStatus\nQUB_0001\tassembled\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "no row for isolate 'QUB_9999'"):
                ns["_read_report_row"](report, "QUB_9999")

    def test_project_assistant_config_supplies_default_reference(self):
        ns = load_backend_helpers({"_resolve_config_path", "_discover_taxonomy_screen_settings"})
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "Taxa_Databases" / "GTDB").mkdir(parents=True)
            reference = root / "Taxa_Databases" / "GTDB" / "gtdb.fna"
            reference.write_text(">ref\nACGT\n", encoding="utf-8")
            config = {
                "project_directory": str(root),
                "defaults": {"assistant": {"ref": "Taxa_Databases/GTDB/gtdb.fna", "taxa": None}},
            }
            (root / "Rumen_Gateway_Assistant.json").write_text(json.dumps(config), encoding="utf-8")
            found = ns["_discover_taxonomy_screen_settings"](root)
            self.assertEqual(found["screen_ref"], str(reference.resolve()))
            self.assertEqual(found["screen_taxa"], "")
            self.assertTrue(found["config_path"].endswith("Rumen_Gateway_Assistant.json"))

    def test_frontend_renders_candidate_hits_and_uses_taxonomy_settings(self):
        html = HTML.read_text(encoding="utf-8")
        self.assertIn('id="taxonomyResults"', html)
        self.assertIn('function renderTaxonomyResults(r)', html)
        self.assertIn("x.NearBest", html)
        self.assertIn("/api/pick-taxonomy-reference", html)
        self.assertIn("/api/set-taxonomy-screen", html)
        self.assertIn("id=\"runBtnSidebar\"", html)
        self.assertIn('id="chooseInputFolderBtn"', html)
        self.assertIn('id="shutdownBtn"', html)
        self.assertNotIn('id="chooseOnboardingResultsBtn"', html)
        self.assertIn('automatically detects Assistant cohort outputs, standalone sample maps, or raw AB1/ABI chromatograms', html)


if __name__ == "__main__":
    unittest.main()
