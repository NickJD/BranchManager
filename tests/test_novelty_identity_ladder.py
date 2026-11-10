import csv
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import shlex

from branchmanager.pipeline import novelty


class NoveltyIdentityLadderTests(unittest.TestCase):
    def test_ladder_retries_only_unresolved_and_keeps_partial_coverage(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            query = root / "queries.fasta"
            ref = root / "reference.fasta"
            output = root / "novelty_reference_nearest_matches.tsv"
            query.write_text(">q_full\n" + "A" * 100 + "\n>q_mid\n" + "C" * 100
                             + "\n>q_partial\n" + "G" * 100 + "\n")
            ref.write_text(">reference\n" + "A" * 100 + "\n")

            calls = []

            def fake_run_cmd(command):
                args = shlex.split(command)
                query_path = Path(args[args.index("--usearch_global") + 1])
                userout = Path(args[args.index("--userout") + 1])
                identity = float(args[args.index("--id") + 1])
                coverage = float(args[args.index("--query_cov") + 1])
                qids = [h for h, _ in novelty.read_fasta(str(query_path))]
                calls.append((identity, coverage, tuple(qids)))
                rows = []
                for qid in qids:
                    if qid == "q_full" and identity == 0.90 and coverage == 0.95:
                        pct, qcov = 97.5, 98.0
                    elif qid == "q_mid" and identity == 0.80 and coverage == 0.95:
                        pct, qcov = 84.0, 96.0
                    elif qid == "q_partial" and identity == 0.90 and coverage == 0.80:
                        pct, qcov = 92.0, 85.0
                    else:
                        continue
                    rows.append([
                        qid, "ref_" + qid, str(pct), "90", "0", "0",
                        "1", "90", "1", "90", "100", "100",
                        str(qcov), "90",
                    ])
                userout.write_text("".join("\t".join(row) + "\n" for row in rows))

            with patch.object(novelty, "run_cmd", side_effect=fake_run_cmd):
                best = novelty._run_nearest_search(
                    str(query), str(ref), output, threads=1,
                )

            self.assertEqual(best["q_full"][0], 97.5)
            self.assertEqual(best["q_full"][6], "FULL")
            self.assertEqual(best["q_mid"][0], 84.0)
            self.assertEqual(best["q_mid"][6], "FULL")
            self.assertEqual(best["q_partial"][6], "PARTIAL")
            self.assertEqual(best["q_partial"][2], 85.0)

            # Resolved full-coverage queries are not searched in later identity rounds.
            self.assertEqual(calls[0][2], ("q_full", "q_mid", "q_partial"))
            self.assertNotIn("q_full", calls[1][2])
            self.assertNotIn("q_mid", calls[-1][2])

            with output.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 3)
            self.assertEqual(rows[2]["match_status"], "PARTIAL_COVERAGE")
            self.assertEqual(rows[1]["match_status"], "MATCH")
            summary = root / "novelty_reference_search_summary.tsv"
            self.assertTrue(summary.exists())


if __name__ == "__main__":
    unittest.main()
