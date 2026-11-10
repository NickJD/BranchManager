"""Integrate genomic-16S-selection outputs (Barrnap + VSEARCH + GTDB-Tk) into the project database.

Genome-derived 16S is supplementary GENOME_BARRNAP evidence. It never replaces an
accepted Sanger marker, GTDB-Tk taxonomy and 16S taxonomy stay in separate fields,
and a genome record marks the isolate as already sequenced.
"""
from __future__ import annotations

import csv
from pathlib import Path

from branchmanager.marker_evidence import (
    ensure_marker_evidence_schema, register_supplemental_marker_fasta, sync_registry_with_sequences,
)
from branchmanager.utils.fasta import read_fasta

GOOD_QUALITY_CLASS = "FULL_LENGTH_LIKE"
MIN_COMPLETENESS = 90.0
MAX_CONTAMINATION = 5.0


def _read_tsv(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _float(value):
    try:
        return None if value in (None, "", "NA") else float(value)
    except (TypeError, ValueError):
        return None


def _bool(value):
    text = str(value or "").strip().lower()
    if text in {"1", "yes", "y", "true", "pass", "passed"}:
        return True
    if text in {"0", "no", "n", "false", "fail", "failed"}:
        return False
    return None


def classify_genome_16s(quality_class: str, genome_qc: bool | None, has_sequence: bool) -> tuple[str, str, str]:
    """Return (tier, qc_class, detail). Unknown genome metrics do not count as a failure."""
    if not has_sequence:
        return "LOW", "GENOME_DERIVED_FAIL", "no_16s_sequence"
    if quality_class != GOOD_QUALITY_CLASS:
        return "LOW", "GENOME_DERIVED_FAIL", f"16s_quality_class={quality_class or 'UNKNOWN'}"
    if genome_qc is False:
        return "LOW", "GENOME_DERIVED_FAIL", "genome_qc_failed"
    return "MODERATE", "GENOME_DERIVED", f"16s_quality_class={quality_class}"


def import_genome_16s(
    db, genomic_dir, out_dir, *, genome_metadata=None, dataset="Genomes", partner_id="",
    genome_id_isolate_id: bool = False,
) -> list[dict]:
    """Register genomes and their selected Barrnap 16S as supplementary marker evidence.

    genome_metadata TSV maps genome_id to sequence_id and may include partner_id, dataset,
    accession, completeness, contamination, genome_qc_pass, and genome_path. Set
    genome_id_isolate_id only when genome identifiers are deliberately stable isolate IDs.
    """
    genomic_dir = Path(genomic_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ensure_marker_evidence_schema(db)
    summary = {r["genome_id"]: r for r in _read_tsv(genomic_dir / "genome_16s_summary.tsv")}
    audit_path = genomic_dir / "taxonomy_assignment_audit.tsv"
    audit = {r["genome_id"]: r for r in _read_tsv(audit_path)} if audit_path.is_file() else {}
    fasta_path = genomic_dir / "best_16s_per_genome.fna"
    sequences = {str(h).split()[0]: s for h, s in read_fasta(str(fasta_path))} if fasta_path.is_file() else {}
    meta = {}
    if genome_metadata:
        meta = {r.get("genome_id", ""): r for r in _read_tsv(genome_metadata)}

    unmapped = [
        genome_id for genome_id in summary
        if not (meta.get(genome_id, {}).get("sequence_id") or meta.get(genome_id, {}).get("isolate_id"))
        and not genome_id_isolate_id
    ]
    if unmapped:
        raise ValueError(
            "Genome-to-isolate mapping is required for every genome; add sequence_id to genome_metadata "
            "or explicitly set genome_id_isolate_id=True. Unmapped genome IDs: "
            + ", ".join(sorted(unmapped)[:10])
        )

    known = set(db.get_all_ids())
    existing_meta = db.get_sequencing_metadata_for_ids()
    results, genome_rows, new_by_dataset, meta_rows, marker_rows = [], [], {}, [], []
    keep_fasta = []
    for genome_id, row in sorted(summary.items()):
        m = meta.get(genome_id, {})
        sid = (m.get("sequence_id") or m.get("isolate_id") or genome_id).strip()
        ds = m.get("dataset") or dataset
        pid = m.get("partner_id") or partner_id or existing_meta.get(sid, {}).get("partner_id", "") or sid
        if sid not in known:
            new_by_dataset.setdefault(ds, []).append((sid, ""))
            known.add(sid)
        completeness = _float(m.get("completeness"))
        contamination = _float(m.get("contamination"))
        explicit = _bool(m.get("genome_qc_pass"))
        if explicit is None and completeness is not None and contamination is not None:
            explicit = completeness >= MIN_COMPLETENESS and contamination <= MAX_CONTAMINATION
        status = "GENOME_QC_PASSED" if explicit else ("GENOME_QC_FAILED" if explicit is False else "SEQUENCED")
        genome_rows.append({
            "genome_id": genome_id, "sequence_id": sid, "accession": m.get("accession") or row.get("original_genome_id") or genome_id,
            "genome_status": status, "genome_qc_pass": explicit, "completeness": completeness,
            "contamination": contamination, "gtdb_taxonomy": row.get("gtdbtk_taxonomy", ""),
            "genome_path": m.get("genome_path") or row.get("input_assembly", ""),
            "notes": f"barrnap_status={row.get('barrnap_status', '')}; gtdb16s_concordance={row.get('gtdbtk_concordance', '')}",
            "source_file": str(genomic_dir / "genome_16s_summary.tsv"),
        })
        old = existing_meta.get(sid, {})
        meta_rows.append({
            "id": sid, "partner_id": pid, "dataset": old.get("dataset") or ds,
            "selected_for_sequencing": bool(old.get("selected_for_sequencing")),
            "selected_for_wgs": True, "source_id": sid, "source_file": str(genomic_dir),
            "raw_selected_value": "genome_available",
            "raw_commitment_value": old.get("raw_commitment_value", ""),
            "project_metadata": {
                **(old.get("project_metadata") or {}),
                "genome_available": True,
                "genome_qc_status": status,
                "genome_resequencing_eligible": explicit is False,
            },
        })
        seq = sequences.get(row.get("best_candidate_id", "")) or sequences.get(genome_id, "")
        tier, qc, detail = classify_genome_16s(row.get("best_quality_class", ""), explicit, bool(seq))
        concordance = row.get("gtdbtk_concordance", "")
        warning = "GTDBTK_16S_DISCORDANT" if concordance.startswith("DISCORDANT") else ""
        if seq:
            keep_fasta.append((genome_id, seq))
            marker_rows.append({
                "sequence_id": sid, "fasta_id": genome_id, "dataset": old.get("dataset") or ds,
                "evidence_tier": tier, "qc_class": qc, "qc_detail": detail, "genome_id": genome_id,
                "genome_taxonomy": row.get("gtdbtk_taxonomy", ""), "genome_16s_concordance": concordance,
                "taxonomy": audit.get(genome_id, {}).get("vsearch_best_reference_taxonomy", ""),
                "taxonomy_source": "VSEARCH_16S_GTDB",
                "taxonomy_identity": row.get("best_vsearch_identity_pct", ""),
                "query_coverage": row.get("best_vsearch_query_coverage_pct", ""),
                "source_file": str(fasta_path),
            })
        results.append({
            "genome_id": genome_id, "sequence_id": sid, "evidence_tier": tier if seq else "NONE",
            "qc_class": qc if seq else "NO_16S", "concordance": concordance, "warning": warning,
            "detail": detail,
        })

    for ds, records in new_by_dataset.items():
        db.insert_sequences(records, dataset=ds)
    db.upsert_genome_records(genome_rows)
    db.upsert_sequencing_metadata(meta_rows)
    for r in genome_rows:
        db.update_isolate_status(r["sequence_id"], r["genome_status"], detail=r["notes"], source_file=r["source_file"])

    if marker_rows:
        meta_tsv = out_dir / "genome_16s_marker_metadata.tsv"
        fields = sorted({k for r in marker_rows for k in r})
        with open(meta_tsv, "w", newline="") as handle:
            w = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            w.writeheader()
            w.writerows(marker_rows)
        fasta_out = out_dir / "genome_16s_markers.fna"
        with open(fasta_out, "w") as handle:
            for gid, seq in keep_fasta:
                handle.write(f">{gid}\n{seq}\n")
        status = {
            (r["sequence_id"], r.get("genome_id", "")): r["result"]
            for r in register_supplemental_marker_fasta(
                db, fasta_out, meta_tsv, default_source_type="GENOME_BARRNAP")
        }
        for r in results:
            r["import_result"] = status.get((r["sequence_id"], r["genome_id"]), "")
            if not r["import_result"] and r["evidence_tier"] != "NONE":
                r["import_result"] = "IMPORTED"
    sync_registry_with_sequences(db)
    report = out_dir / "genome_16s_import_report.tsv"
    fields = ["genome_id", "sequence_id", "evidence_tier", "qc_class", "concordance", "warning", "detail", "import_result"]
    with open(report, "w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(results)
    return results
