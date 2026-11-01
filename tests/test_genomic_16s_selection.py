from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from branchmanager.tools import genomic_16s_selection as tool  # noqa: E402


def test_parse_lineage_and_taxonomic_concordance():
    gtdb = "d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__Lactobacillus acidophilus"
    ref = "d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__Lactobacillus acidophilus"
    status, common, disagreement = tool.compare_taxonomies(
        gtdb, ref, ref_identity_pct=99.8, query_coverage_pct=99.0, target_coverage_pct=99.0,
        minimum_identity_pct=94.5, minimum_query_coverage_pct=95.0, minimum_target_coverage_pct=80.0,
    )
    assert status == "CONCORDANT_SPECIES"
    assert common == "s"
    assert disagreement == ""

    discordant = ref.replace("p__Bacillota", "p__Pseudomonadota")
    status, common, disagreement = tool.compare_taxonomies(
        gtdb, discordant, ref_identity_pct=99.8, query_coverage_pct=99.0, target_coverage_pct=99.0,
        minimum_identity_pct=94.5, minimum_query_coverage_pct=95.0, minimum_target_coverage_pct=80.0,
    )
    assert status == "DISCORDANT_AT_PHYLUM"
    assert disagreement == "p"

    status, common, disagreement = tool.compare_taxonomies(
        gtdb, ref, ref_identity_pct=97.0, query_coverage_pct=99.0, target_coverage_pct=99.0,
        minimum_identity_pct=94.5, minimum_query_coverage_pct=95.0, minimum_target_coverage_pct=80.0,
        species_identity_threshold_pct=98.65,
    )
    assert status == "CONCORDANT_GENUS"
    assert common == "g"
    assert disagreement == ""


def test_reference_taxonomy_header():
    description = (
        "RS_GCF_031457235.1 d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
        "o__Burkholderiales;f__Burkholderiaceae_C;g__Hydrogenophaga;"
        "s__Hydrogenophaga laconesensis [locus_tag=NZ_X-#1] [location=5..1531]"
    )
    tax = tool.reference_taxonomy_from_description(description)
    assert tax.startswith("d__Bacteria;p__Pseudomonadota")
    assert "locus_tag" not in tax


def _write_executable(path: Path, contents: str) -> Path:
    path.write_text(contents, encoding="utf-8")
    path.chmod(0o755)
    return path


def _read_tsv(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_end_to_end_fake_barrnap_vsearch(tmp_path: Path):
    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    outdir = tmp_path / "out"
    # 1,800 bp contig; fake Barrnap reports a 1,500-bp 16S feature away from the edge.
    genome_seq = ("ACGT" * 450)
    (genome_dir / "G001.fna").write_text(f">ctg1\n{genome_seq}\n", encoding="utf-8")

    barrnap = _write_executable(tmp_path / "barrnap", '''#!/usr/bin/env python3
import sys
if "--version" in sys.argv:
    print("barrnap 1.10.6")
    raise SystemExit(0)
print("##gff-version 3")
print("ctg1\\tinfernal\\trRNA\\t101\\t1600\\t1.7e-49\\t+\\t.\\tName=16S_rRNA;product=16S ribosomal RNA")
''')
    vsearch = _write_executable(tmp_path / "vsearch", '''#!/usr/bin/env python3
import sys
from pathlib import Path
if "--version" in sys.argv:
    print("vsearch v2.test")
    raise SystemExit(0)
def value(flag):
    return sys.argv[sys.argv.index(flag) + 1]
query = Path(value("--usearch_global"))
userout = Path(value("--userout"))
notmatched = Path(value("--notmatched"))
queries = []
for line in query.read_text().splitlines():
    if line.startswith(">"):
        queries.append(line[1:].split()[0])
with userout.open("w") as out:
    for q in queries:
        out.write(f"{q}\\tREF1\\t99.5\\t1490\\t100.0\\t98.8\\t3\\t0\\t+\\n")
notmatched.write_text("")
''')
    reference = tmp_path / "gtdb_ssu.fna"
    reference.write_text(
        ">REF1 d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__Lactobacillus testis [locus_tag=x]\nACGT\n",
        encoding="utf-8",
    )
    summary = tmp_path / "gtdbtk.bac120.summary.tsv"
    summary.write_text(
        "user_genome\tclassification\n"
        "G001\td__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__Lactobacillus testis\n",
        encoding="utf-8",
    )
    args = argparse.Namespace(
        genome_dir=genome_dir, out=outdir, reference=reference,
        gtdb_summary=[summary], kingdom="bac", barrnap=str(barrnap), vsearch=str(vsearch),
        jobs=1, threads_per_genome=1, vsearch_threads=1, minimum_search_identity=70.0,
        minimum_search_query_coverage=70.0, max_hits_per_candidate=20,
        minimum_taxonomy_identity=94.5, minimum_taxonomy_query_coverage=95.0,
        minimum_taxonomy_target_coverage=80.0, species_identity_threshold=98.65,
        near_top_identity_delta=0.5, min_full_length=1400, max_full_length=1700,
        preferred_length=1500, max_ambiguity_percent=1.0,
    )
    assert tool.run_cli(args) == 0

    fasta = (outdir / "best_16s_per_genome.fna").read_text(encoding="utf-8")
    assert fasta.startswith(">G001\n")
    candidate_rows = _read_tsv(outdir / "16s_candidates.tsv")
    assert len(candidate_rows) == 1
    assert candidate_rows[0]["selected_best_for_genome"] == "yes"
    assert candidate_rows[0]["gtdbtk_concordance"] == "CONCORDANT_SPECIES"
    summary_rows = _read_tsv(outdir / "genome_16s_summary.tsv")
    assert summary_rows[0]["selection_status"] == "SELECTED_GOOD_CANDIDATE"
    assignments = _read_tsv(outdir / "branchmanager_taxa_assignments.tsv")
    assert assignments[0]["ID"] == "G001"
    assert assignments[0]["Taxon"].startswith("d__Bacteria")
    assert (outdir / "gtdb16s_hits.tsv").is_file()
    assert (outdir / "run_manifest.json").is_file()


def test_vsearch_supplies_fallback_taxonomy_without_gtdbtk(tmp_path: Path):
    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    outdir = tmp_path / "out"
    genome_seq = "ACGT" * 450
    (genome_dir / "G002.fna").write_text(f">ctg1\n{genome_seq}\n", encoding="utf-8")

    barrnap = _write_executable(tmp_path / "barrnap", r'''#!/usr/bin/env python3
import sys
if "--version" in sys.argv:
    print("barrnap test")
    raise SystemExit(0)
print("##gff-version 3")
print("ctg1\tinfernal\trRNA\t101\t1600\t1.7e-49\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA")
''')
    vsearch = _write_executable(tmp_path / "vsearch", r'''#!/usr/bin/env python3
import sys
from pathlib import Path
if "--version" in sys.argv:
    print("vsearch test")
    raise SystemExit(0)
def value(flag):
    return sys.argv[sys.argv.index(flag) + 1]
queries = [line[1:].split()[0] for line in Path(value("--usearch_global")).read_text().splitlines() if line.startswith(">")]
with Path(value("--userout")).open("w") as out:
    for q in queries:
        out.write(f"{q}\tREF1\t99.5\t1490\t100.0\t98.8\t3\t0\t+\n")
Path(value("--notmatched")).write_text("")
''')
    reference = tmp_path / "reference.fna"
    reference.write_text(
        ">REF1 d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__Lactobacillus testis [locus_tag=x]\nACGT\n",
        encoding="utf-8",
    )
    args = argparse.Namespace(
        genome_dir=genome_dir, out=outdir, reference=reference, gtdb_summary=[], kingdom="bac",
        barrnap=str(barrnap), vsearch=str(vsearch), jobs=1, threads_per_genome=1, vsearch_threads=1,
        minimum_search_identity=70.0, minimum_search_query_coverage=70.0, max_hits_per_candidate=20,
        minimum_taxonomy_identity=94.5, minimum_taxonomy_query_coverage=95.0,
        minimum_taxonomy_target_coverage=80.0, species_identity_threshold=98.65,
        near_top_identity_delta=0.5, min_full_length=1400, max_full_length=1700,
        preferred_length=1500, max_ambiguity_percent=1.0,
    )
    assert tool.run_cli(args) == 0
    assignments = _read_tsv(outdir / "branchmanager_taxa_assignments.tsv")
    assert assignments[0]["ID"] == "G002"
    assert assignments[0]["Taxon"].endswith("s__Lactobacillus testis")
    audit = _read_tsv(outdir / "taxonomy_assignment_audit.tsv")
    assert audit[0]["taxonomy_source"] == "VSEARCH_GTDB_SSU"
