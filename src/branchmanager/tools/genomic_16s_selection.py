#!/usr/bin/env python3
"""BranchManager genomic 16S extraction, GTDB SSU matching and representative selection.

Runs Barrnap on each assembled genome, extracts every 16S rRNA copy, searches all
candidates against a GTDB SSU FASTA with VSEARCH, compares the best supported
reference taxonomy with optional GTDB-Tk summaries, and writes exactly one
best-available genome-identifiable marker per genome with at least one detected
16S candidate. All candidates and decisions remain auditable.

This is a standalone BranchManager tool module: CLI registration is handled in
branchmanager.cli through the standard branchmanager.tools tool registry.
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

from branchmanager.utils.io import open_text, sha256sum, write_tsv

FASTA_SUFFIXES = (".fna.gz", ".fasta.gz", ".fas.gz", ".fa.gz", ".ffn.gz",
                  ".fna", ".fasta", ".fas", ".fa", ".ffn")
DNA_COMPLEMENT = str.maketrans(
    "ACGTRYMKBDHVNacgtrymkbdhvn",
    "TGCAYRKMVHDBNtgcayrkmvhdbn",
)
RANKS = ("d", "p", "c", "o", "f", "g", "s")
RANK_NAMES = {
    "d": "domain", "p": "phylum", "c": "class", "o": "order",
    "f": "family", "g": "genus", "s": "species",
}


@dataclass
class Candidate:
    genome_id: str
    original_genome_id: str
    assembly_path: str
    contig: str
    start: int
    end: int
    strand: str
    sequence: str
    barrnap_evalue: float
    copy_number: int = 0
    contig_length: int = 0
    touches_contig_end: bool = False
    ambiguous_bases: int = 0
    ambiguous_pct: float = 0.0
    quality_class: str = ""
    candidate_id: str = ""
    selected: bool = False

    @property
    def length(self) -> int:
        return len(self.sequence)


def strip_fasta_suffix(name: str) -> str:
    low = name.lower()
    for suffix in sorted(FASTA_SUFFIXES, key=len, reverse=True):
        if low.endswith(suffix):
            return name[:-len(suffix)]
    return Path(name).stem


def safe_token(value: str) -> str:
    """Create a deterministic, simple token; original values are stored separately."""
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    token = re.sub(r"_+", "_", token).strip("_.-")
    if not token:
        token = "genome"
    if not token[0].isalnum():
        token = "G_" + token
    return token


def genome_aliases(value: str) -> set[str]:
    """Aliases for joining assembly stems to GTDB-Tk user_genome labels."""
    raw = str(value or "").strip()
    if not raw:
        return set()
    base = strip_fasta_suffix(Path(raw).name)
    values = {raw, Path(raw).name, base, safe_token(raw), safe_token(base)}
    # GTDB-Tk summaries commonly omit the NCBI download suffix `_genomic`.
    for item in list(values):
        values.add(re.sub(r"(?:[._-])genomic$", "", item, flags=re.IGNORECASE))
    return {v for v in values if v}


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current: Optional[str] = None
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                if current in records:
                    raise ValueError(f"Duplicate FASTA sequence ID {current!r} in {path}")
                records[current] = []
            elif current is None:
                raise ValueError(f"Sequence text before FASTA header in {path}")
            else:
                records[current].append(line)
    return {key: "".join(chunks).upper() for key, chunks in records.items()}


def iter_fasta_headers(path: Path) -> Iterable[tuple[str, str]]:
    """Yield (first-token ID, full description after '>') without loading reference bases."""
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                description = line[1:].strip()
                if description:
                    yield description.split()[0], description


def read_gtdb_summaries(paths: list[Path]) -> tuple[dict[str, str], dict[str, str]]:
    """Load GTDB-Tk taxonomy and retain the user_genome ID used in the summary."""
    taxonomy: dict[str, str] = {}
    matched_label: dict[str, str] = {}
    for path in paths:
        with path.open("rt", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames:
                raise ValueError(f"GTDB-Tk summary has no header: {path}")
            fields = {f.lower().strip(): f for f in reader.fieldnames if f}
            id_field = next((fields[key] for key in ("user_genome", "genome_id", "genome", "name", "user genome") if key in fields), None)
            tax_field = next((fields[key] for key in ("classification", "gtdb_taxonomy", "taxonomy") if key in fields), None)
            if id_field is None or tax_field is None:
                raise ValueError(
                    f"GTDB-Tk summary {path} must contain a user_genome (or genome_id) and classification (or taxonomy) column; "
                    f"found: {reader.fieldnames}"
                )
            for row in reader:
                raw_id = str(row.get(id_field) or "").strip()
                tax = str(row.get(tax_field) or "").strip()
                if not raw_id:
                    continue
                for alias in genome_aliases(raw_id):
                    # First match wins deterministically if multiple summary files include aliases.
                    taxonomy.setdefault(alias, tax)
                    matched_label.setdefault(alias, raw_id)
    return taxonomy, matched_label


def parse_lineage(taxonomy: str) -> dict[str, str]:
    """Parse a GTDB lineage such as d__Bacteria;p__Bacillota;... into rank values."""
    out: dict[str, str] = {}
    for token in str(taxonomy or "").strip().split(";"):
        token = token.strip()
        match = re.match(r"^([dpcofgs])__(.*)$", token, flags=re.IGNORECASE)
        if match:
            value = match.group(2).strip()
            if value and value.lower() not in {"na", "none", "unclassified", "unassigned"}:
                out[match.group(1).lower()] = value
    return out


def reference_taxonomy_from_description(description: str) -> str:
    """Read the lineage embedded in a GTDB SSU representative FASTA header."""
    # GTDB SSU headers are `ACCESSION d__...;... [locus_tag=...] [location=...]`.
    body = description.split(maxsplit=1)[1] if " " in description else ""
    body = re.split(r"\s+\[", body, maxsplit=1)[0].strip()
    return body if re.search(r"(?:^|;)\s*d__", body) else ""


def compare_taxonomies(gtdbtk_taxonomy: str, ref_taxonomy: str, *,
                       ref_identity_pct: float | None, query_coverage_pct: float | None,
                       target_coverage_pct: float | None, minimum_identity_pct: float,
                       minimum_query_coverage_pct: float, minimum_target_coverage_pct: float,
                       species_identity_threshold_pct: float = 98.65) -> tuple[str, str, str]:
    """Return concordance status, deepest common rank, and first disagreement rank.

    Comparison is reported only when the VSEARCH hit passes explicit identity and
    query-coverage cut-offs. This avoids treating an arbitrary low-similarity top hit
    as meaningful taxonomic evidence.
    """
    if not str(gtdbtk_taxonomy or "").strip():
        return "GTDB_TK_NOT_AVAILABLE", "", ""
    if not str(ref_taxonomy or "").strip() or ref_identity_pct is None:
        return "NO_REFERENCE_TAXONOMY", "", ""
    if query_coverage_pct is None or query_coverage_pct < minimum_query_coverage_pct:
        return "INSUFFICIENT_QUERY_COVERAGE", "", ""
    if target_coverage_pct is None or target_coverage_pct < minimum_target_coverage_pct:
        return "INSUFFICIENT_TARGET_COVERAGE", "", ""
    if ref_identity_pct < minimum_identity_pct:
        return "REFERENCE_IDENTITY_BELOW_ASSIGNMENT_THRESHOLD", "", ""
    gtdb_lineage = parse_lineage(gtdbtk_taxonomy)
    ref_lineage = parse_lineage(ref_taxonomy)
    if not gtdb_lineage:
        return "GTDB_TK_TAXONOMY_UNPARSEABLE", "", ""
    if not ref_lineage:
        return "REFERENCE_TAXONOMY_UNPARSEABLE", "", ""
    common_rank = ""
    for rank in RANKS:
        # Do not report species-level concordance/disagreement below the
        # commonly used 98.65% 16S species-screening heuristic. It is not a
        # species boundary, so genus-level agreement may still be reported.
        if rank == "s" and ref_identity_pct < species_identity_threshold_pct:
            continue
        gv, rv = gtdb_lineage.get(rank), ref_lineage.get(rank)
        if not gv or not rv:
            continue
        if gv.casefold() != rv.casefold():
            return f"DISCORDANT_AT_{RANK_NAMES[rank].upper()}", common_rank, rank
        common_rank = rank
    if common_rank == "s":
        return "CONCORDANT_SPECIES", common_rank, ""
    if common_rank == "g":
        return "CONCORDANT_GENUS", common_rank, ""
    if common_rank == "f":
        return "CONCORDANT_FAMILY", common_rank, ""
    if common_rank:
        return f"CONCORDANT_{RANK_NAMES[common_rank].upper()}", common_rank, ""
    return "NO_COMPARABLE_RANK", "", ""


def parse_gff_for_16s(gff_text: str) -> list[tuple[str, int, int, str, float, str]]:
    hits = []
    for line in gff_text.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            continue
        seqid, _source, feature, start, end, score, strand, _phase, attrs = fields
        attrs_lower = attrs.lower()
        is_16s = (
            "16s ribosomal rna" in attrs_lower
            or re.search(r"(?:^|[;=])name=16s[_ ]?rrna(?:;|$)", attrs_lower) is not None
            or re.search(r"(?:^|[;=])name=16s_rna(?:;|$)", attrs_lower) is not None
        )
        if feature.lower() != "rrna" or not is_16s:
            continue
        try:
            start_i, end_i = int(start), int(end)
            evalue = float(score) if score not in (".", "") else float("inf")
        except ValueError:
            continue
        if strand not in {"+", "-"} or start_i < 1 or end_i < start_i:
            continue
        hits.append((seqid, start_i, end_i, strand, evalue, attrs))
    return hits


def quality_class(candidate: Candidate, *, min_full_length: int, max_full_length: int,
                  max_ambiguity_pct: float) -> str:
    if (min_full_length <= candidate.length <= max_full_length
            and not candidate.touches_contig_end
            and candidate.ambiguous_pct <= max_ambiguity_pct):
        return "FULL_LENGTH_LIKE"
    if min_full_length <= candidate.length <= max_full_length:
        return "FULL_LENGTH_REVIEW_EDGE_OR_AMBIGUITY"
    if min_full_length - 200 <= candidate.length <= max_full_length + 100:
        return "PARTIAL_OR_UNUSUAL_LENGTH_REVIEW"
    return "SHORT_OR_UNUSUAL_LENGTH_REVIEW"


def quality_rank(quality: str) -> int:
    return {
        "FULL_LENGTH_LIKE": 0,
        "FULL_LENGTH_REVIEW_EDGE_OR_AMBIGUITY": 1,
        "PARTIAL_OR_UNUSUAL_LENGTH_REVIEW": 2,
        "SHORT_OR_UNUSUAL_LENGTH_REVIEW": 3,
    }.get(quality, 4)




def write_fasta(path: Path, records: Iterable[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def _load_reference_taxonomy(path: Path) -> dict[str, str]:
    output: dict[str, str] = {}
    for identifier, description in iter_fasta_headers(path):
        output[identifier] = reference_taxonomy_from_description(description)
    return output


def _version(executable: str, args: list[str]) -> str:
    try:
        result = subprocess.run([executable, *args], text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, check=False, timeout=30)
        return re.sub(r"\s+", " ", (result.stdout + " " + result.stderr).strip()) or f"exit={result.returncode}"
    except Exception as exc:
        return f"unavailable: {exc}"


def _run_one(payload: tuple) -> tuple[str, str, str, str, list[Candidate], str, str, str]:
    (assembly_path, genome_id, raw_genome_id, outdir, kingdom, barrnap_exe,
     threads_per_genome, min_full_length, max_full_length, max_ambiguity_pct) = payload
    gff_path = outdir / "gff" / f"{genome_id}.gff3"
    log_path = outdir / "logs" / f"{genome_id}.log"
    candidates: list[Candidate] = []
    try:
        assembly = read_fasta(assembly_path)
        # Some older Barrnap builds infer input format from the filename extension.
        # Normalise compressed/non-.fna inputs into a temporary .fna while preserving contig IDs.
        barrnap_input = assembly_path
        temp_dir: Optional[tempfile.TemporaryDirectory] = None
        if not assembly_path.name.lower().endswith(".fna"):
            temp_dir = tempfile.TemporaryDirectory(prefix=f"bm16s_{genome_id}_")
            barrnap_input = Path(temp_dir.name) / f"{genome_id}.fna"
            write_fasta(barrnap_input, assembly.items())
        cmd = [barrnap_exe, "--kingdom", kingdom, "--threads", str(threads_per_genome), str(barrnap_input)]
        try:
            proc = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        finally:
            if temp_dir is not None:
                temp_dir.cleanup()
        gff_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        gff_path.write_text(proc.stdout, encoding="utf-8")
        log_path.write_text("COMMAND: " + " ".join(cmd) + "\nRETURN_CODE: " + str(proc.returncode)
                            + "\n\n" + proc.stderr, encoding="utf-8")
        if proc.returncode != 0:
            return genome_id, raw_genome_id, str(assembly_path), "BARRNAP_FAILED", [], str(gff_path), str(log_path), f"Barrnap exited {proc.returncode}"
        for contig, start, end, strand, evalue, _attrs in parse_gff_for_16s(proc.stdout):
            if contig not in assembly:
                continue
            contig_seq = assembly[contig]
            if end > len(contig_seq):
                continue
            seq = contig_seq[start - 1:end]
            if strand == "-":
                seq = seq.translate(DNA_COMPLEMENT)[::-1]
            ambiguous = sum(base not in "ACGT" for base in seq)
            ambiguity_pct = 100.0 * ambiguous / len(seq) if seq else 100.0
            candidates.append(Candidate(
                genome_id=genome_id, original_genome_id=raw_genome_id,
                assembly_path=str(assembly_path), contig=contig, start=start, end=end,
                strand=strand, sequence=seq, barrnap_evalue=evalue,
                contig_length=len(contig_seq), touches_contig_end=(start <= 1 or end >= len(contig_seq)),
                ambiguous_bases=ambiguous, ambiguous_pct=ambiguity_pct,
            ))
        # Coordinate-order copy IDs are stable across repeated runs on the same assembly.
        candidates.sort(key=lambda c: (c.contig, c.start, c.end, c.strand))
        for number, candidate in enumerate(candidates, 1):
            candidate.copy_number = number
            candidate.quality_class = quality_class(candidate, min_full_length=min_full_length,
                                                    max_full_length=max_full_length,
                                                    max_ambiguity_pct=max_ambiguity_pct)
            candidate.candidate_id = (
                f"{genome_id}__16S_copy{number:02d}__contig={safe_token(candidate.contig)}"
                f"__start={candidate.start}__end={candidate.end}__strand={'plus' if candidate.strand == '+' else 'minus'}"
            )
        return genome_id, raw_genome_id, str(assembly_path), ("OK" if candidates else "NO_16S_HIT"), candidates, str(gff_path), str(log_path), ""
    except Exception as exc:
        try:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            log_path.write_text(f"Input/Barrnap error: {type(exc).__name__}: {exc}\n", encoding="utf-8")
        except Exception:
            logging.getLogger(__name__).debug('Failed to write error log for %s', genome_id)
        return genome_id, raw_genome_id, str(assembly_path), "INPUT_OR_BARRNAP_ERROR", [], str(gff_path), str(log_path), f"{type(exc).__name__}: {exc}"


def _read_vsearch_hits(path: Path, reference_taxonomy: dict[str, str]) -> tuple[dict[str, list[dict]], list[dict]]:
    fieldnames = ["candidate_id", "target_id", "identity_pct", "alignment_length", "query_coverage_pct",
                  "target_coverage_pct", "mismatches", "gap_openings", "query_strand", "reference_taxonomy"]
    by_query: dict[str, list[dict]] = {}
    raw_rows: list[dict] = []
    if not path.exists():
        return by_query, raw_rows
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for vals in reader:
            if len(vals) < 9:
                continue
            query, target, identity, alnlen, qcov, tcov, mism, opens, qstrand = vals[:9]
            try:
                identity_value = float(identity)
                if identity_value <= 1.0:
                    identity_value *= 100.0
            except (ValueError, TypeError):
                identity_value = None
            try: aln_value = int(float(alnlen))
            except (ValueError, TypeError): aln_value = 0
            try: qcov_value = float(qcov)
            except (ValueError, TypeError): qcov_value = 0.0
            try: tcov_value = float(tcov)
            except (ValueError, TypeError): tcov_value = 0.0
            row = {
                "candidate_id": query,
                "target_id": "" if target == "*" else target,
                "identity_pct": identity_value if identity_value is not None else "",
                "alignment_length": aln_value,
                "query_coverage_pct": qcov_value,
                "target_coverage_pct": tcov_value,
                "mismatches": mism,
                "gap_openings": opens,
                "query_strand": qstrand,
                "reference_taxonomy": reference_taxonomy.get(target, "") if target != "*" else "",
            }
            raw_rows.append(row)
            if target != "*" and identity_value is not None:
                by_query.setdefault(query, []).append(row)
    for query in by_query:
        by_query[query].sort(key=lambda row: (
            -float(row.get("identity_pct") or 0),
            -float(row.get("query_coverage_pct") or 0),
            -float(row.get("target_coverage_pct") or 0),
            -int(row.get("alignment_length") or 0),
            str(row.get("target_id") or ""),
        ))
    return by_query, raw_rows


def _common_taxonomy(taxonomies: list[str]) -> tuple[str, str]:
    lineages = [parse_lineage(tax) for tax in taxonomies if tax]
    if not lineages:
        return "", ""
    common = {}
    for rank in RANKS:
        values = [lineage.get(rank) for lineage in lineages]
        if any(value is None for value in values) or len({str(v).casefold() for v in values}) != 1:
            break
        common[rank] = values[0]
    rank = max(common, key=lambda r: RANKS.index(r)) if common else ""
    if not rank:
        return "", ""
    return rank, ";".join(f"{r}__{common[r]}" for r in RANKS if r in common)


def _candidate_reference_evidence(candidate: Candidate, hits: list[dict], *, gtdb_taxonomy: str,
                                  minimum_identity_pct: float, minimum_query_coverage_pct: float,
                                  minimum_target_coverage_pct: float, species_identity_threshold_pct: float,
                                  near_top_delta_pct: float) -> dict:
    hit_rows = [h for h in hits if h.get("target_id")]
    best = hit_rows[0] if hit_rows else None
    if best is None:
        return {
            "vsearch_best_hit": "", "vsearch_best_identity_pct": "", "vsearch_best_query_coverage_pct": "",
            "vsearch_best_target_coverage_pct": "", "vsearch_best_alignment_length": "", "vsearch_best_reference_taxonomy": "",
            "vsearch_n_hits": 0, "near_top_hit_count": 0, "near_top_taxonomic_consensus_rank": "",
            "near_top_taxonomic_consensus": "", "gtdbtk_concordance": "NO_VSEARCH_HIT",
            "gtdbtk_common_rank": "", "gtdbtk_disagreement_rank": "",
        }
    best_identity = float(best["identity_pct"])
    best_qcov = float(best["query_coverage_pct"])
    best_tax = str(best.get("reference_taxonomy") or "")
    near_hits = [h for h in hit_rows if float(h["identity_pct"]) >= best_identity - near_top_delta_pct
                 and float(h["query_coverage_pct"]) >= minimum_query_coverage_pct
                 and float(h["target_coverage_pct"]) >= minimum_target_coverage_pct]
    consensus_rank, consensus_tax = _common_taxonomy([str(h.get("reference_taxonomy") or "") for h in near_hits])
    # Limit the exported common-lineage assignment to genus below the species
    # similarity heuristic even if the reference header contains species labels.
    if best_identity < species_identity_threshold_pct and consensus_rank == "s":
        consensus_parts = [part for part in consensus_tax.split(";") if not part.startswith("s__")]
        consensus_tax = ";".join(consensus_parts)
        consensus_rank = "g" if any(part.startswith("g__") for part in consensus_parts) else ""
    concordance, common_rank, disagreement_rank = compare_taxonomies(
        gtdb_taxonomy, best_tax, ref_identity_pct=best_identity,
        query_coverage_pct=best_qcov,
        target_coverage_pct=float(best["target_coverage_pct"]),
        minimum_identity_pct=minimum_identity_pct,
        minimum_query_coverage_pct=minimum_query_coverage_pct,
        minimum_target_coverage_pct=minimum_target_coverage_pct,
        species_identity_threshold_pct=species_identity_threshold_pct,
    )
    return {
        "vsearch_best_hit": best["target_id"],
        "vsearch_best_identity_pct": f"{best_identity:.4f}",
        "vsearch_best_query_coverage_pct": f"{best_qcov:.4f}",
        "vsearch_best_target_coverage_pct": f"{float(best['target_coverage_pct']):.4f}",
        "vsearch_best_alignment_length": best["alignment_length"],
        "vsearch_best_reference_taxonomy": best_tax,
        "vsearch_n_hits": len(hit_rows),
        "near_top_hit_count": len(near_hits),
        "near_top_taxonomic_consensus_rank": consensus_rank,
        "near_top_taxonomic_consensus": consensus_tax,
        "gtdbtk_concordance": concordance,
        "gtdbtk_common_rank": common_rank,
        "gtdbtk_disagreement_rank": disagreement_rank,
    }


def _concordance_rank(status: str) -> int:
    if status == "CONCORDANT_SPECIES": return 0
    if status == "CONCORDANT_GENUS": return 1
    if status == "CONCORDANT_FAMILY": return 2
    if status.startswith("CONCORDANT_"): return 3
    if status.startswith("DISCORDANT_AT_"): return 5
    return 4


def _selection_key(candidate: Candidate, evidence: dict, preferred_length: int) -> tuple:
    # Physical sequence quality takes precedence. Taxonomic concordance and match
    # strength break ties after full/partial-length, contig-edge and ambiguity evidence.
    try: identity = float(evidence.get("vsearch_best_identity_pct") or -1)
    except (TypeError, ValueError): identity = -1
    try: qcov = float(evidence.get("vsearch_best_query_coverage_pct") or -1)
    except (TypeError, ValueError): qcov = -1
    try: tcov = float(evidence.get("vsearch_best_target_coverage_pct") or -1)
    except (TypeError, ValueError): tcov = -1
    status = str(evidence.get("gtdbtk_concordance") or "")
    return (
        quality_rank(candidate.quality_class),
        int(candidate.touches_contig_end),
        candidate.ambiguous_pct,
        abs(candidate.length - preferred_length),
        _concordance_rank(status),
        -identity,
        -qcov,
        -tcov,
        candidate.barrnap_evalue,
        candidate.contig,
        candidate.start,
    )


def _run_vsearch(candidate_fasta: Path, reference_fasta: Path, raw_hits_path: Path, *,
                 vsearch_exe: str, threads: int, minimum_search_identity_pct: float,
                 minimum_search_query_coverage_pct: float, max_hits: int, no_hit_fasta: Path) -> list[str]:
    """Run a high-stringency adaptive search, advancing only unresolved candidates.

    The default path is 98/95, 97/95, 96/95, 95/95. User-specified lower
    identity/coverage floors are appended as rescue tiers, never used first.
    """
    import tempfile

    from branchmanager.utils.fasta import read_fasta, write_fasta

    identity_floor = float(minimum_search_identity_pct)
    coverage_floor = float(minimum_search_query_coverage_pct)
    identity_steps = [98.0, 97.0, 96.0, 95.0, 93.0, 90.0, 85.0, 80.0, 75.0]
    identity_steps = [value for value in identity_steps if value >= identity_floor]
    if not identity_steps or identity_steps[-1] > identity_floor:
        identity_steps.append(identity_floor)
    coverage_steps = [max(95.0, coverage_floor)]
    for value in (90.0, 85.0, 80.0, 75.0, 70.0):
        if value < coverage_steps[-1] and value >= coverage_floor:
            coverage_steps.append(value)
    if coverage_steps[-1] > coverage_floor:
        coverage_steps.append(coverage_floor)
    tiers = [(identity, coverage_steps[0]) for identity in identity_steps]
    tiers.extend((identity_floor, coverage) for coverage in coverage_steps[1:])

    records = list(read_fasta(str(candidate_fasta)))
    sequence_by_id = dict(records)
    unresolved = set(sequence_by_id)
    hits_by_query: dict[str, dict[str, tuple[tuple, str]]] = {}
    commands: list[str] = []
    raw_hits_path.parent.mkdir(parents=True, exist_ok=True)
    no_hit_fasta.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="branchmanager_16s_ladder_", dir=str(raw_hits_path.parent)) as temp_name:
        temp_dir = Path(temp_name)
        for tier_index, (identity, coverage) in enumerate(tiers, 1):
            if not unresolved:
                break
            tier_query = temp_dir / f"queries_tier_{tier_index}.fna"
            tier_hits = temp_dir / f"hits_tier_{tier_index}.tsv"
            write_fasta([(qid, sequence_by_id[qid]) for qid in sorted(unresolved)], str(tier_query))
            cmd = [
                vsearch_exe, "--usearch_global", str(tier_query), "--db", str(reference_fasta),
                "--id", str(identity / 100.0), "--query_cov", str(coverage / 100.0),
                "--maxaccepts", str(max_hits), "--maxrejects", "256", "--maxhits", str(max_hits),
                "--strand", "both", "--userout", str(tier_hits),
                "--userfields", "query+target+id+alnlen+qcov+tcov+mism+opens+qstrand",
                "--threads", str(threads),
            ]
            commands.append(" ".join(cmd))
            proc = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
            tier_hits.with_suffix(tier_hits.suffix + ".log").write_text(
                "COMMAND: " + " ".join(cmd) + f"\nRETURN_CODE: {proc.returncode}\n\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}",
                encoding="utf-8",
            )
            if proc.returncode != 0:
                raise RuntimeError(f"VSEARCH failed at adaptive tier {tier_index} ({identity:.1f}% identity, {coverage:.1f}% coverage); see {tier_hits}.log")

            for line in tier_hits.read_text(encoding="utf-8", errors="replace").splitlines():
                fields = line.split("\t")
                if len(fields) < 9 or fields[0] not in unresolved or fields[1] in ("", "*"):
                    continue
                try:
                    identity_value = float(fields[2])
                    qcov_value = float(fields[4])
                    aln_value = int(float(fields[3]))
                except (ValueError, TypeError):
                    continue
                target = fields[1]
                rank = (identity_value, qcov_value, aln_value, target)
                previous = hits_by_query.setdefault(fields[0], {}).get(target)
                if previous is None or rank > previous[0]:
                    hits_by_query[fields[0]][target] = (rank, line)

            # Three distinct candidates provide a useful near-top comparison set.
            resolved = [qid for qid in unresolved if len(hits_by_query.get(qid, {})) >= min(3, max_hits)]
            unresolved.difference_update(resolved)

    with raw_hits_path.open("w", encoding="utf-8") as handle:
        for qid, _sequence in records:
            candidates = hits_by_query.get(qid, {})
            for _rank, line in sorted(candidates.values(), key=lambda item: item[0], reverse=True):
                handle.write(line + "\n")
    write_fasta([(qid, sequence_by_id[qid]) for qid in sorted(unresolved)], str(no_hit_fasta))
    return commands


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Register CLI options for `branchmanager genomic-16s-selection`."""
    parser.add_argument("--genome-dir", required=True, type=Path,
                        help="Directory containing one assembled genome FASTA per isolate; searched recursively.")
    parser.add_argument("--reference", "--ref", dest="reference", required=True, type=Path,
                        help="GTDB SSU representative FASTA (.fna or .fna.gz), e.g. gtdb_r232_bac_arch_ssu_reps.fna.gz.")
    parser.add_argument("--gtdb-summary", "--gtdbtk-summary", dest="gtdb_summary", action="append", default=[], type=Path,
                        help="Optional GTDB-Tk summary TSV. Repeat for bacterial and archaeal summaries if needed.")
    parser.add_argument("--kingdom", choices=("bac", "arc"), default="bac",
                        help="Barrnap model: bac=bacteria or arc=Archaea (default: bac). This subtool selects prokaryotic 16S/SSU genes.")
    parser.add_argument("--barrnap", default="barrnap", help="Barrnap executable or path (default: barrnap on PATH).")
    parser.add_argument("--vsearch", default="vsearch", help="VSEARCH executable or path (default: vsearch on PATH).")
    parser.add_argument("--jobs", type=int, default=8,
                        help="Number of genomes processed concurrently by Barrnap (default: 8).")
    parser.add_argument("--threads-per-genome", type=int, default=1,
                        help="Threads per Barrnap process; jobs × this value should fit the CPU allocation (default: 1).")
    parser.add_argument("--vsearch-threads", type=int, default=8,
                        help="CPU threads for the VSEARCH phase (default: 8).")
    parser.add_argument("--minimum-search-identity", type=float, default=95.0,
                        help="Minimum identity percentage passed to VSEARCH to collect candidate hits (default: 95).")
    parser.add_argument("--minimum-search-query-coverage", type=float, default=95.0,
                        help="Minimum query coverage percentage for VSEARCH search hits (default: 95).")
    parser.add_argument("--max-hits-per-candidate", type=int, default=20,
                        help="Maximum accepted GTDB reference hits retained per 16S candidate (default: 20).")
    parser.add_argument("--minimum-taxonomy-identity", type=float, default=94.5,
                        help="Minimum identity percentage before a best-hit taxonomy comparison is considered interpretable (default: 94.5; heuristic).")
    parser.add_argument("--minimum-taxonomy-query-coverage", type=float, default=95.0,
                        help="Minimum query coverage percentage for the GTDB-Tk concordance assessment (default: 95).")
    parser.add_argument("--minimum-taxonomy-target-coverage", type=float, default=80.0,
                        help="Minimum GTDB reference-target coverage percentage before taxonomy concordance is assessed (default: 80).")
    parser.add_argument("--species-identity-threshold", type=float, default=98.65,
                        help="16S identity heuristic required before species-level concordance is reported (default: 98.65%).")
    parser.add_argument("--near-top-identity-delta", type=float, default=0.5,
                        help="Treat reference hits within this many identity percentage points of the best as near-ties (default: 0.5).")
    parser.add_argument("--min-full-length", type=int, default=1400,
                        help="Lower bound for full-length-like 16S candidates (default: 1400 bp).")
    parser.add_argument("--max-full-length", type=int, default=1700,
                        help="Upper bound for full-length-like 16S candidates (default: 1700 bp).")
    parser.add_argument("--preferred-length", type=int, default=1500,
                        help="Preferred sequence length used to break quality-equivalent ties (default: 1500 bp).")
    parser.add_argument("--max-ambiguity-percent", type=float, default=1.0,
                        help="Maximum ambiguous-base percentage for FULL_LENGTH_LIKE status (default: 1).")
    parser.add_argument("-o", "--out", required=True, type=Path,
                        help="Output directory for candidate FASTAs, VSEARCH results, summaries, GFFs and logs.")


def run_cli(args: argparse.Namespace) -> int:
    log = logging.getLogger("branchmanager.genomic_16s_selection")
    genome_root = args.genome_dir.expanduser().resolve()
    outdir = args.out.expanduser().resolve()
    reference = args.reference.expanduser().resolve()
    if not genome_root.is_dir():
        raise SystemExit(f"[genomic-16s-selection] Genome directory does not exist: {genome_root}")
    if not reference.is_file():
        raise SystemExit(f"[genomic-16s-selection] GTDB SSU reference does not exist: {reference}")
    for summary in args.gtdb_summary:
        if not summary.is_file():
            raise SystemExit(f"[genomic-16s-selection] GTDB-Tk summary does not exist: {summary}")
    for name, value in (("jobs", args.jobs), ("threads-per-genome", args.threads_per_genome),
                        ("vsearch-threads", args.vsearch_threads), ("max-hits-per-candidate", args.max_hits_per_candidate)):
        if value < 1:
            raise SystemExit(f"[genomic-16s-selection] --{name} must be >= 1")
    for name, value in (("minimum-search-identity", args.minimum_search_identity),
                        ("minimum-search-query-coverage", args.minimum_search_query_coverage),
                        ("minimum-taxonomy-identity", args.minimum_taxonomy_identity),
                        ("minimum-taxonomy-query-coverage", args.minimum_taxonomy_query_coverage),
                        ("minimum-taxonomy-target-coverage", args.minimum_taxonomy_target_coverage),
                        ("species-identity-threshold", args.species_identity_threshold),
                        ("max-ambiguity-percent", args.max_ambiguity_percent)):
        if not 0 <= value <= 100:
            raise SystemExit(f"[genomic-16s-selection] --{name} must be between 0 and 100")
    if args.min_full_length < 1 or args.max_full_length < args.min_full_length:
        raise SystemExit("[genomic-16s-selection] --min-full-length must be positive and <= --max-full-length")
    if not 0 <= args.near_top_identity_delta <= 100:
        raise SystemExit("[genomic-16s-selection] --near-top-identity-delta must be between 0 and 100")
    if shutil.which(args.barrnap) is None and not Path(args.barrnap).is_file():
        raise SystemExit(f"[genomic-16s-selection] Could not find Barrnap executable: {args.barrnap}")
    if shutil.which(args.vsearch) is None and not Path(args.vsearch).is_file():
        raise SystemExit(f"[genomic-16s-selection] Could not find VSEARCH executable: {args.vsearch}")
    try:
        outdir.relative_to(genome_root)
        raise SystemExit("[genomic-16s-selection] --out must be outside --genome-dir to prevent outputs being rediscovered as input genomes.")
    except ValueError:
        pass

    genome_paths = sorted(path for path in genome_root.rglob("*") if path.is_file()
                          and any(path.name.lower().endswith(suffix) for suffix in FASTA_SUFFIXES))
    if not genome_paths:
        raise SystemExit(f"[genomic-16s-selection] No supported genome FASTA files found under {genome_root}")
    genomes = []
    seen_ids: dict[str, Path] = {}
    for path in genome_paths:
        original_id = strip_fasta_suffix(path.name)
        genome_id = safe_token(original_id)
        if genome_id in seen_ids:
            raise SystemExit(
                f"[genomic-16s-selection] Genome identifiers collide after normalisation: {seen_ids[genome_id]} and {path} -> {genome_id!r}. "
                "Rename files or process them separately."
            )
        seen_ids[genome_id] = path
        genomes.append((path, genome_id, original_id))

    outdir.mkdir(parents=True, exist_ok=True)
    for subdir in ("gff", "logs", "intermediate"):
        (outdir / subdir).mkdir(parents=True, exist_ok=True)

    gtdb_taxonomy, gtdb_labels = read_gtdb_summaries(args.gtdb_summary)
    def taxonomy_for(genome_id: str, original_id: str) -> tuple[str, str]:
        for alias in list(genome_aliases(original_id)) + list(genome_aliases(genome_id)):
            if alias in gtdb_taxonomy:
                return gtdb_taxonomy[alias], gtdb_labels.get(alias, "")
        return "", ""

    barrnap_version = _version(args.barrnap, ["--version"])
    vsearch_version = _version(args.vsearch, ["--version"])
    log.info("Processing %d genome assemblies; Barrnap=%s; VSEARCH=%s", len(genomes), barrnap_version, vsearch_version)
    tasks = [
        (path, genome_id, original, outdir, args.kingdom, args.barrnap, args.threads_per_genome,
         args.min_full_length, args.max_full_length, args.max_ambiguity_percent)
        for path, genome_id, original in genomes
    ]
    results = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(_run_one, task) for task in tasks]
        for i, future in enumerate(as_completed(futures), start=1):
            result = future.result()
            results.append(result)
            if i % 25 == 0 or i == len(futures):
                print(f"[genomic-16s-selection] Barrnap: processed {i}/{len(futures)} genomes", file=sys.stderr)
    results.sort(key=lambda row: row[0])

    candidate_by_id: dict[str, Candidate] = {}
    all_candidates: list[Candidate] = []
    genome_rows: dict[str, dict] = {}
    all_fastas: list[tuple[str, str]] = []
    run_rows = []
    for genome_id, original_id, assembly_path, run_status, candidates, gff_path, log_path, error in results:
        gtdbtk_taxon, matched_gtdb_id = taxonomy_for(genome_id, original_id)
        genome_rows[genome_id] = {
            "genome_id": genome_id, "original_genome_id": original_id, "input_assembly": assembly_path,
            "gtdbtk_user_genome": matched_gtdb_id, "gtdbtk_taxonomy": gtdbtk_taxon,
            "barrnap_status": run_status, "n_16s_hits": len(candidates), "gff_path": gff_path,
            "log_path": log_path, "error": error,
        }
        run_rows.append({**genome_rows[genome_id], "barrnap_version": barrnap_version,
                         "vsearch_version": vsearch_version})
        for candidate in candidates:
            candidate_by_id[candidate.candidate_id] = candidate
            all_candidates.append(candidate)
            all_fastas.append((candidate.candidate_id, candidate.sequence))
    all_fasta = outdir / "all_16s_candidates.fna"
    write_fasta(all_fasta, all_fastas)

    hits_by_candidate: dict[str, list[dict]] = {}
    raw_hit_rows: list[dict] = []
    ref_tax: dict[str, str] = {}
    vsearch_command: list[str] = []
    with tempfile.TemporaryDirectory(prefix="branchmanager_gtdb_ssu_") as temp_name:
        temp_dir = Path(temp_name)
        ref_for_vsearch = reference
        if reference.name.lower().endswith(".gz"):
            ref_for_vsearch = temp_dir / "gtdb_ssu_reference.fna"
            with gzip.open(reference, "rb") as src, ref_for_vsearch.open("wb") as dst:
                shutil.copyfileobj(src, dst, length=1024 * 1024)
        ref_tax = _load_reference_taxonomy(reference)
        # If compressed, parsing through gzip was already handled. VSEARCH target IDs are first header tokens.
        raw_hits_path = outdir / "gtdb16s_hits.raw.tsv"
        no_hit_fasta = outdir / "intermediate" / "vsearch_no_hits.fna"
        if all_candidates:
            vsearch_command = _run_vsearch(
                all_fasta, ref_for_vsearch, raw_hits_path, vsearch_exe=args.vsearch,
                threads=args.vsearch_threads, minimum_search_identity_pct=args.minimum_search_identity,
                minimum_search_query_coverage_pct=args.minimum_search_query_coverage,
                max_hits=args.max_hits_per_candidate, no_hit_fasta=no_hit_fasta,
            )
            hits_by_candidate, raw_hit_rows = _read_vsearch_hits(raw_hits_path, ref_tax)
        else:
            raw_hits_path.write_text("", encoding="utf-8")
            no_hit_fasta.write_text("", encoding="utf-8")
        # Preserve original GTDB reference path in the manifest, not temporary decompressed filename.

    hits_header = ["candidate_id", "target_id", "identity_pct", "alignment_length", "query_coverage_pct",
                   "target_coverage_pct", "mismatches", "gap_openings", "query_strand", "reference_taxonomy"]
    write_tsv(outdir / "gtdb16s_hits.tsv", hits_header, raw_hit_rows)
    try:
        (outdir / "gtdb16s_hits.raw.tsv").unlink(missing_ok=True)
    except OSError:
        pass

    candidate_rows = []
    representative_records = []
    assignment_rows = []
    taxonomy_assessment_rows = []
    best_candidate_by_genome: dict[str, tuple[Candidate, dict]] = {}
    for candidate in all_candidates:
        gtdbtk_taxon = genome_rows[candidate.genome_id]["gtdbtk_taxonomy"]
        evidence = _candidate_reference_evidence(
            candidate, hits_by_candidate.get(candidate.candidate_id, []), gtdb_taxonomy=gtdbtk_taxon,
            minimum_identity_pct=args.minimum_taxonomy_identity,
            minimum_query_coverage_pct=args.minimum_taxonomy_query_coverage,
            minimum_target_coverage_pct=args.minimum_taxonomy_target_coverage,
            species_identity_threshold_pct=args.species_identity_threshold,
            near_top_delta_pct=args.near_top_identity_delta,
        )
        candidate_rows.append({
            "genome_id": candidate.genome_id, "original_genome_id": candidate.original_genome_id,
            "gtdbtk_user_genome": genome_rows[candidate.genome_id]["gtdbtk_user_genome"],
            "gtdbtk_taxonomy": gtdbtk_taxon, "candidate_id": candidate.candidate_id,
            "copy_number_coordinate_order": candidate.copy_number, "contig_id": candidate.contig,
            "start_1based": candidate.start, "end_1based_inclusive": candidate.end,
            "strand": candidate.strand, "contig_length_bp": candidate.contig_length,
            "length_bp": candidate.length, "ambiguous_bases": candidate.ambiguous_bases,
            "ambiguous_pct": f"{candidate.ambiguous_pct:.4f}",
            "touches_contig_end": str(candidate.touches_contig_end).lower(),
            "barrnap_evalue": "inf" if math.isinf(candidate.barrnap_evalue) else f"{candidate.barrnap_evalue:.6g}",
            "quality_class": candidate.quality_class, **evidence, "selected_best_for_genome": "no",
        })
        taxonomy_assessment_rows.append({
            "genome_id": candidate.genome_id, "candidate_id": candidate.candidate_id,
            "gtdbtk_taxonomy": gtdbtk_taxon, **evidence,
        })
    # Choose exactly one sequence per genome with a Barrnap hit; quality is primary,
    # taxonomy concordance and VSEARCH evidence are tie-breakers, never hard exclusions.
    candidate_lookup = {row["candidate_id"]: row for row in candidate_rows}
    for genome_id, genome_row in genome_rows.items():
        candidates = [c for c in all_candidates if c.genome_id == genome_id]
        if not candidates:
            genome_row.update({
                "best_candidate_id": "", "best_length_bp": "", "best_quality_class": "",
                "best_vsearch_hit": "", "best_vsearch_identity_pct": "", "best_vsearch_query_coverage_pct": "",
                "gtdbtk_concordance": "NO_16S_CANDIDATE", "selected_taxonomy": genome_row["gtdbtk_taxonomy"],
                "selected_taxonomy_source": "GTDB-Tk" if genome_row["gtdbtk_taxonomy"] else "",
                "selection_status": genome_row["barrnap_status"],
            })
            continue
        def evidence_for(c):
            return _candidate_reference_evidence(
                c, hits_by_candidate.get(c.candidate_id, []), gtdb_taxonomy=genome_row["gtdbtk_taxonomy"],
                minimum_identity_pct=args.minimum_taxonomy_identity,
                minimum_query_coverage_pct=args.minimum_taxonomy_query_coverage,
                minimum_target_coverage_pct=args.minimum_taxonomy_target_coverage,
                species_identity_threshold_pct=args.species_identity_threshold,
                near_top_delta_pct=args.near_top_identity_delta,
            )
        ranked = sorted(candidates, key=lambda c: _selection_key(c, evidence_for(c), args.preferred_length))
        best = ranked[0]
        best.selected = True
        best_evidence = evidence_for(best)
        candidate_lookup[best.candidate_id]["selected_best_for_genome"] = "yes"
        selection_status = "SELECTED_GOOD_CANDIDATE" if best.quality_class == "FULL_LENGTH_LIKE" else "SELECTED_BEST_AVAILABLE_REVIEW"
        # Prefer the genome-wide GTDB-Tk assignment when available. Without it,
        # use the conservative common lineage of near-top SSU matches only when
        # both identity and query coverage meet the explicit assignment thresholds.
        selected_taxonomy = str(genome_row["gtdbtk_taxonomy"] or "")
        selected_taxonomy_source = "GTDB-Tk" if selected_taxonomy else ""
        if not selected_taxonomy:
            try:
                best_identity_for_assignment = float(best_evidence.get("vsearch_best_identity_pct") or -1)
                best_qcov_for_assignment = float(best_evidence.get("vsearch_best_query_coverage_pct") or -1)
            except (TypeError, ValueError):
                best_identity_for_assignment = best_qcov_for_assignment = -1
            if (best_identity_for_assignment >= args.minimum_taxonomy_identity
                    and best_qcov_for_assignment >= args.minimum_taxonomy_query_coverage
                    and float(best_evidence.get("vsearch_best_target_coverage_pct") or -1) >= args.minimum_taxonomy_target_coverage):
                selected_taxonomy = str(best_evidence.get("near_top_taxonomic_consensus")
                                        or best_evidence.get("vsearch_best_reference_taxonomy") or "")
                if selected_taxonomy:
                    selected_taxonomy_source = "VSEARCH_GTDB_SSU"
        genome_row.update({
            "best_candidate_id": best.candidate_id,
            "best_length_bp": best.length,
            "best_quality_class": best.quality_class,
            "best_vsearch_hit": best_evidence["vsearch_best_hit"],
            "best_vsearch_identity_pct": best_evidence["vsearch_best_identity_pct"],
            "best_vsearch_query_coverage_pct": best_evidence["vsearch_best_query_coverage_pct"],
            "gtdbtk_concordance": best_evidence["gtdbtk_concordance"],
            "selected_taxonomy": selected_taxonomy,
            "selected_taxonomy_source": selected_taxonomy_source,
            "selection_status": selection_status,
        })
        # Keep the FASTA ID equal to the stable genome ID, for direct use as a marker FASTA in BranchManager.
        representative_records.append((genome_id, best.sequence))
        if selected_taxonomy:
            assignment_rows.append({
                "ID": genome_id,
                "Taxon": selected_taxonomy,
                # GTDB-Tk summaries do not provide a scalar confidence in the
                # standard classification column; leave this blank rather than inventing one.
                "Confidence": "",
            })
        best_candidate_by_genome[genome_id] = (best, best_evidence)

    # Selection and per-genome summary are emitted in stable sorted order.
    for row in candidate_rows:
        row["selected_best_for_genome"] = "yes" if row["candidate_id"] in {c.candidate_id for c, _e in best_candidate_by_genome.values()} else "no"
    write_tsv(outdir / "16s_candidates.tsv", list(candidate_rows[0].keys()) if candidate_rows else [
        "genome_id", "original_genome_id", "candidate_id", "contig_id", "start_1based", "end_1based_inclusive",
        "strand", "length_bp", "quality_class", "vsearch_best_hit", "gtdbtk_concordance", "selected_best_for_genome",
    ], candidate_rows)
    write_tsv(outdir / "16s_taxonomy_assessment.tsv", list(taxonomy_assessment_rows[0].keys()) if taxonomy_assessment_rows else [
        "genome_id", "candidate_id", "gtdbtk_taxonomy", "vsearch_best_hit", "vsearch_best_identity_pct",
        "vsearch_best_query_coverage_pct", "gtdbtk_concordance",
    ], taxonomy_assessment_rows)
    write_tsv(outdir / "genome_16s_summary.tsv", list(genome_rows[next(iter(genome_rows))].keys()) if genome_rows else [],
              (genome_rows[key] for key in sorted(genome_rows)))
    write_fasta(outdir / "best_16s_per_genome.fna", representative_records)
    write_tsv(outdir / "branchmanager_taxa_assignments.tsv", ["ID", "Taxon", "Confidence"], assignment_rows)
    write_tsv(outdir / "taxonomy_assignment_audit.tsv", [
        "genome_id", "taxonomy", "taxonomy_source", "gtdbtk_taxonomy", "gtdbtk_concordance",
        "selected_16s_candidate", "vsearch_best_hit", "vsearch_best_identity_pct",
        "vsearch_best_query_coverage_pct", "vsearch_best_reference_taxonomy",
        "near_top_taxonomic_consensus_rank", "near_top_taxonomic_consensus",
    ], ({
        "genome_id": row["genome_id"], "taxonomy": row.get("selected_taxonomy", ""),
        "taxonomy_source": row.get("selected_taxonomy_source", ""),
        "gtdbtk_taxonomy": row.get("gtdbtk_taxonomy", ""),
        "gtdbtk_concordance": row.get("gtdbtk_concordance", ""),
        "selected_16s_candidate": row.get("best_candidate_id", ""),
        "vsearch_best_hit": row.get("best_vsearch_hit", ""),
        "vsearch_best_identity_pct": row.get("best_vsearch_identity_pct", ""),
        "vsearch_best_query_coverage_pct": row.get("best_vsearch_query_coverage_pct", ""),
        "vsearch_best_reference_taxonomy": best_candidate_by_genome.get(row["genome_id"], (None, {}))[1].get("vsearch_best_reference_taxonomy", ""),
        "near_top_taxonomic_consensus_rank": best_candidate_by_genome.get(row["genome_id"], (None, {}))[1].get("near_top_taxonomic_consensus_rank", ""),
        "near_top_taxonomic_consensus": best_candidate_by_genome.get(row["genome_id"], (None, {}))[1].get("near_top_taxonomic_consensus", ""),
    } for row in (genome_rows[key] for key in sorted(genome_rows))))
    no_hits_rows = [row for row in genome_rows.values() if not row.get("best_candidate_id")]
    write_tsv(outdir / "no_16s_or_failed.tsv", list(genome_rows[next(iter(genome_rows))].keys()) if genome_rows else [],
              sorted(no_hits_rows, key=lambda row: row["genome_id"]))
    write_tsv(outdir / "run_manifest.tsv", list(run_rows[0].keys()) if run_rows else [], run_rows)

    manifest = {
        "tool": "branchmanager genomic-16s-selection",
        "tool_version": "1.0.0",
        "started_with": {
            "genome_dir": str(genome_root), "reference": str(reference),
            "reference_sha256": sha256sum(reference), "gtdb_summary": [str(p.resolve()) for p in args.gtdb_summary],
            "gtdb_summary_sha256": {str(p.resolve()): sha256sum(p) for p in args.gtdb_summary},
            "kingdom": args.kingdom, "jobs": args.jobs, "threads_per_genome": args.threads_per_genome,
            "vsearch_threads": args.vsearch_threads, "minimum_search_identity_pct": args.minimum_search_identity,
            "minimum_search_query_coverage_pct": args.minimum_search_query_coverage,
            "minimum_taxonomy_identity_pct": args.minimum_taxonomy_identity,
            "minimum_taxonomy_query_coverage_pct": args.minimum_taxonomy_query_coverage,
            "minimum_taxonomy_target_coverage_pct": args.minimum_taxonomy_target_coverage,
            "species_identity_threshold_pct": args.species_identity_threshold,
            "max_hits_per_candidate": args.max_hits_per_candidate,
            "near_top_identity_delta_pct": args.near_top_identity_delta,
            "min_full_length": args.min_full_length, "max_full_length": args.max_full_length,
            "preferred_length": args.preferred_length, "max_ambiguity_percent": args.max_ambiguity_percent,
        },
        "versions": {"barrnap": barrnap_version, "vsearch": vsearch_version},
        "counts": {
            "genomes": len(genomes), "genomes_with_16s_candidate": len(representative_records),
            "candidate_16s_sequences": len(all_candidates), "genomes_without_candidate": len(no_hits_rows),
        },
        "outputs": {
            "representatives": "best_16s_per_genome.fna", "all_candidates": "all_16s_candidates.fna",
            "candidate_table": "16s_candidates.tsv", "vsearch_hits": "gtdb16s_hits.tsv",
            "taxonomy_assessment": "16s_taxonomy_assessment.tsv", "genome_summary": "genome_16s_summary.tsv",
            "branchmanager_taxa_assignments": "branchmanager_taxa_assignments.tsv",
            "taxonomy_assignment_audit": "taxonomy_assignment_audit.tsv",
            "no_16s_or_failed": "no_16s_or_failed.tsv", "run_manifest_table": "run_manifest.tsv",
        },
        "vsearch_command": vsearch_command,
        "selection_policy": (
            "Sequence quality takes precedence (full-length-like status, contig-edge status, ambiguous bases, "
            "length closeness); GTDB SSU hit evidence and GTDB-Tk concordance only break ties. No candidate is "
            "discarded solely because its taxonomy disagrees with GTDB-Tk."
        ),
    }
    (outdir / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    status_counts: dict[str, int] = {}
    for row in genome_rows.values():
        status = str(row.get("selection_status") or "UNKNOWN")
        status_counts[status] = status_counts.get(status, 0) + 1
    print(
        f"[genomic-16s-selection] Complete: {len(genomes)} genomes; {len(all_candidates)} 16S candidates; "
        f"{len(representative_records)} representative sequence(s).",
        file=sys.stderr,
    )
    print("[genomic-16s-selection] Status: " + ", ".join(f"{k}={v}" for k, v in sorted(status_counts.items())), file=sys.stderr)
    print(f"[genomic-16s-selection] Outputs: {outdir}", file=sys.stderr)
    return 0
