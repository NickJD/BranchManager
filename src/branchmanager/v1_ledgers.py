"""v1 audit ledgers and summary counts for the Annual Report.

One row per isolate (or per marker observation where stated); an isolate with several
evidence sources is counted once in isolate-level tallies.
"""
from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

from branchmanager.marker_evidence import ensure_marker_evidence_schema
from branchmanager.taxonomy import normalise_taxon_name, parse_taxon_string


def _taxon_at(taxonomy: str, ranks=("s", "g")) -> tuple[str, str]:
    parsed = parse_taxon_string(str(taxonomy or ""))
    for rank in ranks:
        value = str(parsed.get(rank) or "").strip()
        if value and value.lower() not in {"na", "none", "unclassified"}:
            return rank, normalise_taxon_name(value)
    return "", ""


def _genus(taxonomy: str) -> str:
    """Genus of a GTDB-style string or a plain MALDI name such as 'Prevotella ruminicola'."""
    text = str(taxonomy or "").strip()
    if not text:
        return ""
    if "__" in text:
        return _taxon_at(text, ("g",))[1]
    return normalise_taxon_name(text.split()[0])


def _compare(a: str, b: str) -> str:
    ga, gb = _genus(a), _genus(b)
    if not ga or not gb:
        return "NOT_COMPARABLE"
    return "CONCORDANT_GENUS" if ga == gb else "DISCORDANT_GENUS"


def _write(path: Path, fields: list[str], rows: list[dict]) -> str:
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return str(path)


def collect_taxonomy_concordance(db) -> list[dict]:
    """Compare marker (VSEARCH), genome (GTDB-Tk) and MALDI taxonomy without merging them."""
    ensure_marker_evidence_schema(db)
    with db.connect() as conn:
        conn.row_factory = None
        primary = {}
        for sid, tax, tier in conn.execute(
            "SELECT sequence_id, COALESCE(taxonomy,''), evidence_tier FROM marker_evidence "
            "WHERE eligible_primary=1 AND active_for_selection=1 ORDER BY evidence_id"
        ):
            primary[sid] = (tax, tier)
        taxonomy = {}
        for sid, tax in conn.execute("SELECT id, taxonomy FROM taxonomy"):
            taxonomy.setdefault(sid, tax or "")
        genomic = {}
        for sid, tax, conc in conn.execute(
            "SELECT sequence_id, COALESCE(taxonomy,''), COALESCE(genome_16s_concordance,'') "
            "FROM marker_evidence WHERE source_type='GENOME_BARRNAP' ORDER BY evidence_id"
        ):
            genomic[sid] = (tax, conc)
        genomes = {}
        for sid, tax in conn.execute("SELECT sequence_id, COALESCE(gtdb_taxonomy,'') FROM genome_records ORDER BY genome_id"):
            genomes.setdefault(sid, tax)
        maldi = {}
        for sid, name in conn.execute(
            "SELECT sequence_id, COALESCE(taxon_name,'') FROM maldi_identifications ORDER BY identification_id"
        ):
            maldi[sid] = name
        ids = sorted(set(primary) | set(genomic) | set(genomes) | set(maldi))
    rows = []
    for sid in ids:
        marker_tax = ""
        if sid in primary:
            marker_tax = primary[sid][0] or taxonomy.get(sid, "")
        genome_tax = genomes.get(sid, "")
        maldi_tax = maldi.get(sid, "")
        g16_tax, g16_conc = genomic.get(sid, ("", ""))
        marker_vs_genome = _compare(marker_tax, genome_tax) if marker_tax and genome_tax else "NOT_COMPARABLE"
        maldi_vs_marker = _compare(maldi_tax, marker_tax) if maldi_tax and marker_tax else "NOT_COMPARABLE"
        genome16_vs_gtdb = g16_conc or ("NOT_COMPARABLE" if not g16_tax else _compare(g16_tax, genome_tax))
        warnings = []
        if marker_vs_genome.startswith("DISCORDANT"):
            warnings.append("GTDBTK_MARKER_DISAGREEMENT")
        if maldi_vs_marker.startswith("DISCORDANT"):
            warnings.append("MALDI_MARKER_DISAGREEMENT")
        if genome16_vs_gtdb.startswith("DISCORDANT"):
            warnings.append("GTDBTK_GENOMIC16S_DISAGREEMENT")
        rows.append({
            "sequence_id": sid, "marker_taxonomy_vsearch": marker_tax, "genomic_16s_taxonomy_vsearch": g16_tax,
            "genome_taxonomy_gtdbtk": genome_tax, "maldi_taxon": maldi_tax,
            "marker_vs_genome": marker_vs_genome, "maldi_vs_marker": maldi_vs_marker,
            "genomic_16s_vs_gtdbtk": genome16_vs_gtdb, "warnings": ";".join(warnings),
        })
    return rows


def write_v1_ledgers(db, output_dir: str | Path) -> dict:
    """Write taxonomy-concordance, unreviewed-marker and extended genome ledgers."""
    ensure_marker_evidence_schema(db)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    concordance = collect_taxonomy_concordance(db)
    tc_path = _write(output / "project_taxonomy_concordance.tsv", [
        "sequence_id", "marker_taxonomy_vsearch", "genomic_16s_taxonomy_vsearch", "genome_taxonomy_gtdbtk",
        "maldi_taxon", "marker_vs_genome", "maldi_vs_marker", "genomic_16s_vs_gtdbtk", "warnings",
    ], concordance)
    with db.connect() as conn:
        unreviewed = [dict(zip(
            ("sequence_id", "dataset", "source_type", "evidence_tier", "qc_class", "manual_review_status",
             "review_decision", "reviewer", "qc_detail"), r))
            for r in conn.execute(
                "SELECT sequence_id, dataset, source_type, evidence_tier, qc_class, manual_review_status, "
                "review_decision, reviewer, COALESCE(qc_detail,'') FROM marker_evidence "
                "WHERE source_type != 'GENOME_BARRNAP' AND (manual_review_status IN ('PENDING','NOT_REVIEWED','REJECTED') "
                "OR evidence_tier='LOW') ORDER BY sequence_id"
            )]
        genome_rows = [dict(zip(
            ("genome_id", "sequence_id", "accession", "genome_status", "genome_qc_pass", "completeness",
             "contamination", "gtdb_taxonomy", "ani_cluster", "genome_path", "genome_sha256", "source_file"), r))
            for r in conn.execute(
                "SELECT genome_id, sequence_id, accession, genome_status, genome_qc_pass, completeness, contamination, "
                "gtdb_taxonomy, ani_cluster, genome_path, genome_sha256, source_file FROM genome_records "
                "ORDER BY sequence_id, genome_id")]
        g16 = {}
        for sid, gid, tier, qc, conc, tax, ident, cov in conn.execute(
            "SELECT sequence_id, genome_id, evidence_tier, qc_class, COALESCE(genome_16s_concordance,''), "
            "COALESCE(taxonomy,''), taxonomy_identity, query_coverage FROM marker_evidence "
            "WHERE source_type='GENOME_BARRNAP' ORDER BY evidence_id"
        ):
            g16[(sid, gid)] = (tier, qc, conc, tax, ident, cov)
        observation_rows = []
        for eid, sid, source, tier, qc, tax, gid, conc, digest in conn.execute(
            "SELECT evidence_id, sequence_id, source_type, evidence_tier, qc_class, COALESCE(taxonomy,''), "
            "COALESCE(genome_id,''), COALESCE(genome_16s_concordance,''), sequence_sha256 "
            "FROM marker_evidence ORDER BY evidence_id"
        ):
            observation_rows.append({
                "sequence_id": sid, "observation_type": "MARKER_EVIDENCE", "observation_id": eid,
                "source_type": source, "genome_id": gid, "taxonomy_source": "VSEARCH_MARKER",
                "taxonomy": tax, "evidence_tier": tier, "qc_class": qc,
                "genome_16s_concordance": conc, "accession": "", "genome_status": "",
                "genome_qc_pass": "", "sequence_sha256": digest or "",
            })
        for gid, sid, accession, status, qc_pass, taxonomy in conn.execute(
            "SELECT genome_id, sequence_id, accession, genome_status, genome_qc_pass, "
            "COALESCE(gtdb_taxonomy,'') FROM genome_records ORDER BY genome_id"
        ):
            observation_rows.append({
                "sequence_id": sid, "observation_type": "GENOME_RECORD", "observation_id": gid,
                "source_type": "GENOME", "genome_id": gid, "taxonomy_source": "GTDB_TK",
                "taxonomy": taxonomy, "evidence_tier": "", "qc_class": "", "genome_16s_concordance": "",
                "accession": accession or "", "genome_status": status, "genome_qc_pass": qc_pass,
                "sequence_sha256": "",
            })
        for ident, sid, taxon_id, name, rank, instrument, library, date, source_file in conn.execute(
            "SELECT identification_id, sequence_id, COALESCE(taxon_id,''), COALESCE(taxon_name,''), "
            "COALESCE(taxonomic_rank,''), COALESCE(instrument,''), COALESCE(library_name,''), "
            "COALESCE(identification_date,''), COALESCE(source_file,'') "
            "FROM maldi_identifications ORDER BY identification_id"
        ):
            observation_rows.append({
                "sequence_id": sid, "observation_type": "MALDI_IDENTIFICATION", "observation_id": ident,
                "source_type": "MALDI", "genome_id": "", "taxonomy_source": "MALDI_TOF",
                "taxonomy": name, "evidence_tier": "", "qc_class": "", "genome_16s_concordance": "",
                "accession": "", "genome_status": "", "genome_qc_pass": "", "sequence_sha256": "",
                "taxon_id": taxon_id, "taxonomic_rank": rank, "instrument": instrument,
                "library_name": library, "identification_date": date, "source_file": source_file,
            })
    for row in genome_rows:
        tier, qc, conc, tax, ident, cov = g16.get((row["sequence_id"], row["genome_id"]), ("", "", "", "", "", ""))
        row.update({
            "genomic_16s_tier": tier, "genomic_16s_qc_class": qc, "gtdbtk_16s_concordance": conc,
            "genomic_16s_taxonomy": tax, "genomic_16s_identity": ident, "genomic_16s_query_coverage": cov,
        })
    um_path = _write(output / "project_unreviewed_markers.tsv", [
        "sequence_id", "dataset", "source_type", "evidence_tier", "qc_class", "manual_review_status",
        "review_decision", "reviewer", "qc_detail"], unreviewed)
    genome_path = _write(output / "project_genome_ledger.tsv", [
        "genome_id", "sequence_id", "accession", "genome_status", "genome_qc_pass", "completeness",
        "contamination", "gtdb_taxonomy", "ani_cluster", "genome_path", "genome_sha256", "source_file",
        "genomic_16s_tier", "genomic_16s_qc_class", "gtdbtk_16s_concordance", "genomic_16s_taxonomy",
        "genomic_16s_identity", "genomic_16s_query_coverage"], genome_rows)
    observation_path = _write(output / "project_taxonomy_concordance_observations.tsv", [
        "sequence_id", "observation_type", "observation_id", "source_type", "genome_id",
        "taxonomy_source", "taxonomy", "evidence_tier", "qc_class", "genome_16s_concordance",
        "accession", "genome_status", "genome_qc_pass", "sequence_sha256", "taxon_id",
        "taxonomic_rank", "instrument", "library_name", "identification_date", "source_file",
    ], observation_rows)
    return {"taxonomy_concordance": tc_path, "unreviewed_markers": um_path, "genomes_v1": genome_path,
            "taxonomy_observations": observation_path, "concordance_rows": concordance}


def v1_summary_counts(db, concordance_rows: list[dict] | None = None) -> list[tuple[str, str, int]]:
    """Return (section, label, count) rows; each isolate is counted once per section."""
    ensure_marker_evidence_schema(db)
    out: list[tuple[str, str, int]] = []
    with db.connect() as conn:
        def one(sql):
            return conn.execute(sql).fetchone()[0]
        out.append(("Totals", "Total isolates", one("SELECT COUNT(*) FROM sequences")))
        for partner, n in conn.execute(
            "SELECT COALESCE(NULLIF(m.partner_id,''),'UNKNOWN'), COUNT(DISTINCT s.id) FROM sequences s "
            "LEFT JOIN sequencing_metadata m ON m.id=s.id GROUP BY 1 ORDER BY 1"):
            out.append(("By partner", partner, n))
        for src, n in conn.execute(
            "SELECT source_type, COUNT(DISTINCT sequence_id) FROM marker_evidence GROUP BY 1 ORDER BY 1"):
            out.append(("By marker source", src, n))
        for qc, n in conn.execute(
            "SELECT qc_class, COUNT(*) FROM marker_evidence GROUP BY 1 ORDER BY 1"):
            out.append(("Marker QC", qc, n))
        for tier, n in conn.execute(
            "SELECT evidence_tier, COUNT(DISTINCT sequence_id) FROM marker_evidence WHERE source_type != 'GENOME_BARRNAP' GROUP BY 1 ORDER BY 1"):
            out.append(("Direct-marker evidence tier", tier, n))
        for state, n in conn.execute(
            "SELECT manual_review_status, COUNT(*) FROM marker_evidence WHERE source_type != 'GENOME_BARRNAP' GROUP BY 1 ORDER BY 1"):
            out.append(("Manual review", state, n))
        marker_ids = "(SELECT DISTINCT sequence_id FROM marker_evidence WHERE source_type != 'GENOME_BARRNAP')"
        genomic_ids = "(SELECT DISTINCT sequence_id FROM marker_evidence WHERE source_type='GENOME_BARRNAP')"
        out.append(("MALDI", "MALDI only (no direct or genomic marker)", one(
            f"SELECT COUNT(DISTINCT m.sequence_id) FROM maldi_identifications m "
            f"WHERE m.sequence_id NOT IN {marker_ids} AND m.sequence_id NOT IN {genomic_ids}")))
        out.append(("MALDI", "MALDI plus direct marker only", one(
            f"SELECT COUNT(DISTINCT m.sequence_id) FROM maldi_identifications m "
            f"WHERE m.sequence_id IN {marker_ids} AND m.sequence_id NOT IN {genomic_ids}")))
        out.append(("MALDI", "MALDI plus genomic 16S only", one(
            f"SELECT COUNT(DISTINCT m.sequence_id) FROM maldi_identifications m "
            f"WHERE m.sequence_id NOT IN {marker_ids} AND m.sequence_id IN {genomic_ids}")))
        out.append(("MALDI", "MALDI plus direct and genomic marker evidence", one(
            f"SELECT COUNT(DISTINCT m.sequence_id) FROM maldi_identifications m "
            f"WHERE m.sequence_id IN {marker_ids} AND m.sequence_id IN {genomic_ids}")))
        out.append(("Genomes", "Genomic 16S records", one(
            "SELECT COUNT(*) FROM marker_evidence WHERE source_type='GENOME_BARRNAP'")))
        out.append(("Genomes", "Isolates with existing genome", one("SELECT COUNT(DISTINCT sequence_id) FROM genome_records")))
        for status, n in conn.execute("SELECT genome_status, COUNT(*) FROM genome_records GROUP BY 1 ORDER BY 1"):
            out.append(("Genome QC", status, n))
    rows = concordance_rows if concordance_rows is not None else collect_taxonomy_concordance(db)
    counts = Counter(w for r in rows for w in r["warnings"].split(";") if w)
    for key in ("GTDBTK_MARKER_DISAGREEMENT", "MALDI_MARKER_DISAGREEMENT", "GTDBTK_GENOMIC16S_DISAGREEMENT"):
        out.append(("Taxonomic disagreements", key, counts.get(key, 0)))
    return out
