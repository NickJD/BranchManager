"""Focused tests for multi-hit near-best taxonomy resolution."""
from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path

SOURCE = Path(__file__).parents[1] / "src/branchmanager/pipeline/paper_trail.py"

# paper_trail imports FASTA helpers at module import time. Stub these here so
# this focused test can run independently of the complete BranchManager tree.
branchmanager = types.ModuleType("branchmanager")
utils = types.ModuleType("branchmanager.utils")
fasta = types.ModuleType("branchmanager.utils.fasta")
fasta.read_fasta = lambda *args, **kwargs: []
fasta.reverse_complement = lambda seq: seq[::-1]
fasta.write_fasta = lambda *args, **kwargs: None
sys.modules.setdefault("branchmanager", branchmanager)
sys.modules.setdefault("branchmanager.utils", utils)
sys.modules.setdefault("branchmanager.utils.fasta", fasta)

spec = importlib.util.spec_from_file_location("paper_trail_taxonomy_test", SOURCE)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
assert spec.loader is not None
spec.loader.exec_module(module)


def parse_taxonomy(lineage: str) -> dict[str, str]:
    result = {}
    for item in (lineage or "").split(";"):
        if "__" in item:
            rank, value = item.split("__", 1)
            result[rank] = value
    return result


def hit(target: str, identity: float, taxonomy: str) -> dict:
    return {"target": target, "identity": identity, "taxonomy": taxonomy}


def test_near_best_hits_agreeing_at_genus_resolve_genus():
    result = module._resolve_near_best_taxonomy([
        hit("ref1", 99.90, "d__Bacteria;p__Bacillota;c__Clostridia;o__Eubacteriales;f__Clostridiaceae;g__Sarcina;s__Sarcina perfringens"),
        hit("ref2", 99.75, "d__Bacteria;p__Bacillota;c__Clostridia;o__Eubacteriales;f__Clostridiaceae;g__Sarcina;s__Sarcina sp."),
    ], parse_taxonomy)
    assert result["assigned_rank"] in {"g", "s"}
    assert "g__Sarcina" in result["taxonomy"]
    assert result["status"] == "RESOLVED"
    assert len(result["near_best_hits"]) == 2


def test_near_best_hits_disagreeing_at_genus_fall_back_to_family():
    result = module._resolve_near_best_taxonomy([
        hit("ref1", 99.90, "d__Bacteria;p__Bacillota;c__Clostridia;o__Eubacteriales;f__Clostridiaceae;g__Sarcina"),
        hit("ref2", 99.70, "d__Bacteria;p__Bacillota;c__Clostridia;o__Eubacteriales;f__Clostridiaceae;g__OtherGenus"),
    ], parse_taxonomy)
    assert result["assigned_rank"] == "f"
    assert result["taxonomy"].endswith("f__Clostridiaceae")
    assert result["status"] == "RESOLVED"


def test_near_best_hits_without_shared_family_do_not_force_a_genus():
    result = module._resolve_near_best_taxonomy([
        hit("ref1", 99.90, "d__Bacteria;p__Bacillota;c__Negativicutes;o__Veillonellales;f__Veillonellaceae;g__Veillonella"),
        hit("ref2", 99.60, "d__Bacteria;p__Bacillota;c__Clostridia;o__Eubacteriales;f__Clostridiaceae;g__Sarcina"),
    ], parse_taxonomy)
    assert result["assigned_rank"] in {"o", "c", "p", "d", "NA"}
    assert result["status"] != "RESOLVED"


def test_distant_lower_identity_hit_does_not_override_best_genus():
    result = module._resolve_near_best_taxonomy([
        hit("ref1", 99.90, "d__Bacteria;p__Bacillota;c__Clostridia;o__Eubacteriales;f__Clostridiaceae;g__Sarcina"),
        hit("ref2", 99.20, "d__Bacteria;p__Bacillota;c__Negativicutes;o__Veillonellales;f__Veillonellaceae;g__Veillonella"),
    ], parse_taxonomy)
    assert len(result["near_best_hits"]) == 1
    assert "g__Sarcina" in result["taxonomy"]


def test_unmapped_near_best_taxonomy_is_not_reported_as_genus():
    result = module._resolve_near_best_taxonomy([
        hit("ref1", 99.90, ""),
    ], parse_taxonomy)
    assert result["assigned_rank"] == "NA"
    assert result["status"] == "TAXONOMY_UNMAPPED"
    assert result["taxonomy"] == ""


def make_read(read_id: str, sequence_id: str):
    return module.SangerRead(
        read_id=read_id, sequence_id=sequence_id, primer="27F", direction="forward",
        source_file="", raw_sequence="A" * 1000, raw_qualities=[40] * 1000,
    )


def test_single_reference_hit_disagreement_is_potential_conflict():
    group = [make_read("read_a", "isolate_1"), make_read("read_b", "isolate_1")]
    summaries = {
        "read_a": {
            "taxonomy": "d__Bacteria;p__Bacillota;c__Clostridia;o__Eubacteriales;f__Clostridiaceae;g__Sarcina",
            "best_identity": 99.8, "best_hit": "ref_a", "best_query_coverage": 99.0,
            "near_best_hit_count": 1, "assigned_rank": "g", "support": "SINGLE_NEAR_BEST_REFERENCE_HIT",
        },
        "read_b": {
            "taxonomy": "d__Bacteria;p__Bacillota;c__Negativicutes;o__Veillonellales;f__Veillonellaceae;g__Veillonella",
            "best_identity": 99.7, "best_hit": "ref_b", "best_query_coverage": 99.0,
            "near_best_hit_count": 1, "assigned_rank": "g", "support": "SINGLE_NEAR_BEST_REFERENCE_HIT",
        },
    }
    result = module._compare_primer_read_taxonomy_group(group, summaries, parse_taxonomy)
    assert result["conflict"] is True
    assert result["status"] == "POTENTIAL_CONFLICT"
    assert result["compared_rank"] == "genus"
    assert "single near-best reference hit" in result["reason"]


def test_multi_reference_disagreement_is_reported_as_conflict():
    group = [make_read("read_a", "isolate_1"), make_read("read_b", "isolate_1")]
    summaries = {
        "read_a": {
            "taxonomy": "d__Bacteria;p__Bacillota;c__Clostridia;o__Eubacteriales;f__Clostridiaceae;g__Sarcina",
            "best_identity": 99.8, "best_hit": "ref_a", "best_query_coverage": 99.0,
            "near_best_hit_count": 3, "assigned_rank": "g", "support": "MULTIPLE_NEAR_BEST_HITS_AGREE_AT_GENUS",
        },
        "read_b": {
            "taxonomy": "d__Bacteria;p__Bacillota;c__Negativicutes;o__Veillonellales;f__Veillonellaceae;g__Veillonella",
            "best_identity": 99.7, "best_hit": "ref_b", "best_query_coverage": 99.0,
            "near_best_hit_count": 2, "assigned_rank": "g", "support": "MULTIPLE_NEAR_BEST_HITS_AGREE_AT_GENUS",
        },
    }
    result = module._compare_primer_read_taxonomy_group(group, summaries, parse_taxonomy)
    assert result["conflict"] is True
    assert result["status"] == "CONFLICT"
    assert result["compared_rank"] == "genus"


def test_genus_disagreement_with_shared_family_is_reported_at_genus():
    group = [make_read("read_a", "isolate_1"), make_read("read_b", "isolate_1")]
    shared_family = "d__Bacteria;p__Bacillota;c__Clostridia;o__Eubacteriales;f__Clostridiaceae;"
    summaries = {
        "read_a": {
            "taxonomy": shared_family + "g__Sarcina", "near_best_hit_count": 2,
            "best_identity": 99.8, "best_hit": "ref_a", "best_query_coverage": 99.0,
        },
        "read_b": {
            "taxonomy": shared_family + "g__OtherGenus", "near_best_hit_count": 2,
            "best_identity": 99.8, "best_hit": "ref_b", "best_query_coverage": 99.0,
        },
    }
    result = module._compare_primer_read_taxonomy_group(group, summaries, parse_taxonomy)
    assert result["status"] == "CONFLICT"
    assert result["compared_rank"] == "genus"


def test_one_unclassified_read_does_not_hide_two_read_agreement():
    group = [
        make_read("read_a", "isolate_1"), make_read("read_b", "isolate_1"),
        make_read("read_c", "isolate_1"),
    ]
    lineage = "d__Bacteria;p__Bacillota;c__Clostridia;o__Eubacteriales;f__Clostridiaceae;g__Sarcina"
    summaries = {
        "read_a": {"taxonomy": lineage, "near_best_hit_count": 2, "best_identity": 99.8,
                   "best_hit": "ref_a", "best_query_coverage": 99.0},
        "read_b": {"taxonomy": lineage, "near_best_hit_count": 2, "best_identity": 99.8,
                   "best_hit": "ref_b", "best_query_coverage": 99.0},
        "read_c": {"taxonomy": "NA", "near_best_hit_count": 0, "best_identity": None,
                   "best_hit": "", "best_query_coverage": None},
    }
    result = module._compare_primer_read_taxonomy_group(group, summaries, parse_taxonomy)
    assert result["status"] == "CONCORDANT"
    assert result["compared_rank"] == "genus"


def test_screen_function_runs_multi_hit_search_and_writes_audit_files(tmp_path, monkeypatch):
    import csv
    from types import SimpleNamespace

    pipeline = types.ModuleType("branchmanager.pipeline")
    taxonomy_module = types.ModuleType("branchmanager.taxonomy")
    taxonomy_module.parse_taxon_string = parse_taxonomy

    taxonomy_map = {
        "ref_a": "d__Bacteria;p__Bacillota;c__Clostridia;o__Eubacteriales;f__Clostridiaceae;g__Sarcina",
        "ref_a2": "d__Bacteria;p__Bacillota;c__Clostridia;o__Eubacteriales;f__Clostridiaceae;g__Sarcina",
        "ref_b": "d__Bacteria;p__Bacillota;c__Negativicutes;o__Veillonellales;f__Veillonellaceae;g__Veillonella",
        "ref_b2": "d__Bacteria;p__Bacillota;c__Negativicutes;o__Veillonellales;f__Veillonellaceae;g__Veillonella",
    }

    def parse_match(parts):
        if len(parts) < 14:
            return None
        return {
            "query": parts[0], "target": parts[1], "identity": float(parts[2]),
            "alignment_length": int(parts[3]), "query_coverage": float(parts[12]),
            "target_coverage": float(parts[13]),
        }

    classify = types.ModuleType("branchmanager.pipeline.classify")
    classify.VSEARCH_USERFIELDS = "query+target+id+alnlen+mism+gaps+qlo+qhi+tlo+thi+ql+tl+qcov+tcov"
    classify._ensure_uncompressed = lambda path, outdir, out_name: path
    classify._load_taxa_map = lambda path: taxonomy_map
    classify._load_taxa_map_from_reference_fasta = lambda path: (taxonomy_map, [])
    classify._vsearch_safe_query = lambda path, outdir: (str(path), {"q0": "read_a", "q1": "read_b"})
    classify._parse_vsearch_match = parse_match
    classify._lookup_tax = lambda target, mapping: (mapping.get(target), None)
    classify.validate_reference_taxonomy_consistency = lambda ref, taxa: []
    classify.write_reference_taxonomy_warnings = lambda outdir, rows: None
    pipeline.classify = classify
    branchmanager.pipeline = pipeline
    sys.modules["branchmanager.pipeline"] = pipeline
    sys.modules["branchmanager.pipeline.classify"] = classify
    sys.modules["branchmanager.taxonomy"] = taxonomy_module

    commands = []
    raw_lines = [
        "q0\tref_a\t99.80\t1000\t2\t0\t1\t1000\t1\t1500\t1000\t1500\t100.00\t66.67",
        "q0\tref_a2\t99.70\t1000\t3\t0\t1\t1000\t1\t1500\t1000\t1500\t100.00\t66.67",
        "q1\tref_b\t99.80\t1000\t2\t0\t1\t1000\t1\t1500\t1000\t1500\t100.00\t66.67",
        "q1\tref_b2\t99.70\t1000\t3\t0\t1\t1000\t1\t1500\t1000\t1500\t100.00\t66.67",
    ]

    def fake_run(command, capture_output, text, check):
        commands.append(command)
        if "--usearch_global" in command:
            outpath = Path(command[command.index("--userout") + 1])
            outpath.write_text("\n".join(raw_lines) + "\n")
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        return SimpleNamespace(returncode=0, stdout="vsearch v2.30.0\n", stderr="")

    monkeypatch.setattr(module.subprocess, "run", fake_run)
    output_dir = tmp_path / "paper_trail"
    output_dir.mkdir()
    trimmed_fasta = output_dir / "trimmed.fasta"
    trimmed_fasta.write_text(">read_a\nAAAA\n>read_b\nCCCC\n")
    reads = [make_read("read_a", "isolate_1"), make_read("read_b", "isolate_1")]

    results, report_path = module._screen_primer_read_taxonomy(
        reads, trimmed_fasta, output_dir,
        ref_fasta="reference.fasta", taxa_tsv="taxonomy.tsv", threads=2,
    )

    search_command = commands[0]
    assert search_command[search_command.index("--maxaccepts") + 1] == "20"
    assert search_command[search_command.index("--maxrejects") + 1] == "32"
    assert search_command[search_command.index("--maxhits") + 1] == "20"
    assert results["isolate_1"]["status"] == "CONFLICT"
    assert Path(report_path).is_file()

    screen_dir = output_dir / "read_taxonomy_screen"
    assert (screen_dir / "read_taxonomy_hits.tsv").is_file()
    assert (screen_dir / "taxonomy_screen_policy.tsv").is_file()
    with open(screen_dir / "taxonomy_screen_policy.tsv") as handle:
        policy = next(csv.DictReader(handle, delimiter="\t"))
    assert policy["MaxAccepts"] == "20"
    assert policy["MaxRejects"] == "32"
    assert policy["MaxHits"] == "20"
    with open(screen_dir / "taxonomy.tsv") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 2
    assert rows[0]["TaxonomySupport"] == "MULTIPLE_NEAR_BEST_HITS_AGREE_AT_GENUS"
