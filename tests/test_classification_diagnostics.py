import unittest
import csv
import tempfile
import os

from branchmanager.pipeline.evidence_policy import classification_evidence_interpretation
from branchmanager.pipeline.background_check import _load_candidate_hits
from branchmanager.pipeline.workflow_helpers import iter_classification_rows

from branchmanager.pipeline.classify import (
    _diagnostic_hit_reason,
    _dominant_taxonomy_consensus,
    _parse_vsearch_match,
    DIAGNOSTIC_SEARCH_TIERS,
    LEGACY_DIAGNOSTIC_MAX_ACCEPTS,
    LEGACY_DIAGNOSTIC_MAX_HITS,
)


class ClassificationConsensusTests(unittest.TestCase):
    def test_contested_species_retains_best_guess(self):
        taxa = {
            "r1": ("Bacteria;Firmicutes;Bacilli;Order;Family;Genus;Genus_a", None),
            "r2": ("Bacteria;Firmicutes;Bacilli;Order;Family;Genus;Genus_a", None),
            "r3": ("Bacteria;Firmicutes;Bacilli;Order;Family;Genus;Genus_b", None),
            "r4": ("Bacteria;Firmicutes;Bacilli;Order;Family;Genus;Genus_c", None),
        }
        hits = [{"target": f"r{i}", "identity": 99.0 - i / 10} for i in range(1, 5)]
        result = _dominant_taxonomy_consensus(hits, taxa)
        self.assertEqual(result["status"], "CONTESTED")
        self.assertTrue(result["taxonomy"].endswith("Genus_a"))
        self.assertEqual(result["rank"], "species")
        self.assertEqual(result["conflicts"][0]["rank"], "species")

    def test_single_reference_is_provisional_not_confident(self):
        taxa = {"r1": ("Bacteria;Firmicutes;Bacilli", None)}
        result = _dominant_taxonomy_consensus([{"target": "r1"}], taxa)
        self.assertEqual(result["taxonomy"], "Bacteria;Firmicutes;Bacilli")
        self.assertEqual(result["status"], "PROVISIONAL")


    def test_plurality_tie_uses_top_ranked_reference(self):
        taxa = {
            "top": ("Bacteria;Firmicutes;Bacilli;OrderA", None),
            "other": ("Bacteria;Proteobacteria;Gammaproteobacteria;OrderB", None),
        }
        result = _dominant_taxonomy_consensus(
            [{"target": "top", "identity": 99.9}, {"target": "other", "identity": 99.8}],
            taxa,
        )
        self.assertEqual(result["status"], "CONTESTED")
        self.assertEqual(result["taxonomy"], "Bacteria")
        self.assertEqual(result["rank"], "domain")

    def test_no_taxonomy_remains_unresolved(self):
        result = _dominant_taxonomy_consensus([{"target": "r1"}], {})
        self.assertIsNone(result["taxonomy"])
        self.assertEqual(result["status"], "UNRESOLVED")

    def test_identity_and_coverage_limited(self):
        self.assertEqual(
            _diagnostic_hit_reason([
                {"identity": 92.0, "query_coverage": 70.0},
                {"identity": 88.0, "query_coverage": 90.0},
            ]),
            "IDENTITY_AND_COVERAGE_LIMITED",
        )

    def test_coverage_limited(self):
        self.assertEqual(
            _diagnostic_hit_reason([{"identity": 92.0, "query_coverage": 80.0}]),
            "QUERY_COVERAGE_LIMITED",
        )

    def test_identity_limited(self):
        self.assertEqual(
            _diagnostic_hit_reason([{"identity": 88.0, "query_coverage": 90.0}]),
            "IDENTITY_LIMITED",
        )


    def test_vsearch_tabular_row_is_parsed(self):
        row = "q1\tref1\t99.5\t1450\t7\t0\t1\t1450\t1\t1450\t1450\t1500\t100.0\t96.7"
        parsed = _parse_vsearch_match(row.split("\t"))
        self.assertIsNotNone(parsed)
        self.assertEqual(parsed["query"], "q1")
        self.assertEqual(parsed["target"], "ref1")
        self.assertAlmostEqual(parsed["query_coverage"], 100.0)
        self.assertAlmostEqual(parsed["target_coverage"], 96.7)

    def test_fully_supported_species_consensus_is_confident(self):
        taxa = {
            "r1": ("Bacteria;Firmicutes;Bacilli;Order;Family;Genus;Genus_a", None),
            "r2": ("Bacteria;Firmicutes;Bacilli;Order;Family;Genus;Genus_a", None),
            "r3": ("Bacteria;Firmicutes;Bacilli;Order;Family;Genus;Genus_a", None),
        }
        result = _dominant_taxonomy_consensus(
            [{"target": "r1"}, {"target": "r2"}, {"target": "r3"}], taxa
        )
        self.assertEqual(result["status"], "CONFIDENT")
        self.assertEqual(result["rank"], "species")

    def test_no_hit_at_diagnostic_floor(self):
        self.assertEqual(_diagnostic_hit_reason([]), "NO_HIT_AT_DIAGNOSTIC_FLOOR")

    def test_diagnostic_ladder_includes_legacy_background_search(self):
        self.assertEqual(
            DIAGNOSTIC_SEARCH_TIERS[-1],
            (80.0, 70.0),
        )
        self.assertEqual(LEGACY_DIAGNOSTIC_MAX_ACCEPTS, 10)
        self.assertEqual(LEGACY_DIAGNOSTIC_MAX_HITS, 10)


    def test_high_identity_partial_alignment_is_not_high_evidence(self):
        result = classification_evidence_interpretation(
            99.9, 82.2, 69.0, evidence_tier="DIAGNOSTIC_D1",
            assignment_status="LOW_SIMILARITY", candidate_count=1,
        )
        self.assertEqual(result["alignment_evidence_status"], "BOTH_PARTIAL")
        self.assertEqual(result["assignment_evidence_level"], "TENTATIVE")
        self.assertIn("99.9% identity across the aligned region", result["evidence_interpretation"])
        self.assertIn("82.2% query coverage", result["evidence_interpretation"])
        self.assertIn("69.0% reference coverage", result["evidence_interpretation"])
        self.assertIn("diagnostic-only", result["evidence_interpretation"])

    def test_full_length_high_identity_requires_taxonomy_support_for_high_level(self):
        result = classification_evidence_interpretation(
            99.9, 99.8, 99.5, evidence_tier="STANDARD_LADDER",
            assignment_status="CONFIDENT", candidate_count=4,
        )
        self.assertEqual(result["alignment_evidence_status"], "NEAR_FULL_LENGTH")
        self.assertEqual(result["assignment_evidence_level"], "HIGH")

    def test_partial_production_alignment_is_not_high_level(self):
        result = classification_evidence_interpretation(
            99.9, 82.2, 69.0, evidence_tier="STANDARD_LADDER",
            assignment_status="PROVISIONAL", candidate_count=1,
        )
        self.assertEqual(result["assignment_evidence_level"], "LOW")


    def test_background_check_retains_both_alignment_coverages(self):
        with tempfile.NamedTemporaryFile(mode="w", newline="", delete=False) as handle:
            path = handle.name
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow([
                "q1", "ref1", "99.9", "853", "1", "0", "1", "853",
                "1", "853", "1037", "1235", "82.2", "69.0",
            ])
        try:
            hits = _load_candidate_hits(path, "round1")
            hit = hits["q1"][0]
            self.assertEqual(hit[0], "ref1")
            self.assertAlmostEqual(hit[1], 99.9)
            self.assertAlmostEqual(hit[2], 82.2)
            self.assertAlmostEqual(hit[3], 69.0)
            self.assertEqual(hit[4], 853)
        finally:
            os.unlink(path)


    def test_legacy_confidence_column_is_not_reinterpreted_as_source_confidence(self):
        with tempfile.NamedTemporaryFile(mode="w", newline="", delete=False) as handle:
            path = handle.name
            handle.write("ID\tBestHit\tIdentity\tTaxon\tConfidence\n")
            handle.write("q1\tr1\t99.9\td__Bacteria\t0.999\n")
        try:
            row = list(iter_classification_rows(path))[0]
            self.assertIsNone(row["confidence"])
        finally:
            os.unlink(path)


    def test_source_classifier_confidence_is_retained_when_explicitly_named(self):
        with tempfile.NamedTemporaryFile(mode="w", newline="", delete=False) as handle:
            path = handle.name
            handle.write("ID\tBestHit\tIdentity\tTaxon\tSourceClassifierConfidence\n")
            handle.write("q1\tr1\t99.9\td__Bacteria\t0.87\n")
        try:
            row = list(iter_classification_rows(path))[0]
            self.assertAlmostEqual(row["confidence"], 0.87)
        finally:
            os.unlink(path)


if __name__ == "__main__":
    unittest.main()


def test_species_conflict_does_not_discard_supported_genus():
    from branchmanager.pipeline.classify import _dominant_taxonomy_consensus

    taxa = {
        "a": ("Bacteria;Firmicutes;Clostridia;OrderA;FamilyA;GenusA;GenusA species1", None),
        "b": ("Bacteria;Firmicutes;Clostridia;OrderA;FamilyA;GenusA;GenusA species2", None),
    }
    result = _dominant_taxonomy_consensus(
        [{"target": key, "identity": 99.9} for key in ("a", "b")], taxa
    )
    assert result["rank"] == "genus"
    assert result["taxonomy"].endswith("GenusA")
    assert result["status"] == "PROVISIONAL"
    assert any(row["rank"] == "species" for row in result["conflicts"])
