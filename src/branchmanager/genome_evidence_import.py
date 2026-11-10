"""Import partner genome assemblies with GTDB-Tk and CheckM2 evidence.

Partner culture IDs are the join key across workbook metadata, FASTA filenames,
GTDB-Tk user_genome and CheckM2 Name. The workbook resolves those IDs to the
canonical BranchManager sequence IDs.
"""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

from branchmanager.partner_metadata import load_genome_only_partner_workbook_rows
from branchmanager.utils.io import sha256sum


def _read_table(path: str | Path) -> list[dict[str, str]]:
    path = Path(path)
    with open(path, newline="", encoding="utf-8-sig") as handle:
        sample = handle.read(8192)
        handle.seek(0)
        try:
            delimiter = csv.Sniffer().sniff(sample, delimiters="\t,").delimiter
        except csv.Error:
            delimiter = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"No header row found in {path}")
        return [
            {str(k or "").strip().lstrip("\ufeff"): str(v or "").strip()
             for k, v in row.items() if k is not None}
            for row in reader
        ]


def _get(row: dict, *names: str) -> str:
    wanted = {n.casefold().replace(" ", "_") for n in names}
    for key, value in row.items():
        if key.casefold().replace(" ", "_") in wanted and value not in (None, ""):
            return str(value).strip()
    return ""


def _number(value):
    text = str(value or "").strip()
    if not text or text.upper() in {"N/A", "NA", "NONE", "-"}:
        return None
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def _normalise_partner_id(value: str) -> str:
    """Normalise harmless zero-padding differences in trailing numeric IDs."""
    value = Path(str(value or "").strip()).name
    value = re.sub(r"\.(?:fasta|fa|fna)(?:\.gz)?$", "", value, flags=re.I)
    match = re.fullmatch(r"(.*?)(\d+)", value)
    if match:
        return match.group(1).casefold() + str(int(match.group(2)))
    return value.casefold()


def _sha256(path: Path) -> str:
    return sha256sum(path)


def _version_from_manifest(path: Path | None, keys: tuple[str, ...]) -> str:
    if not path or not path.is_file():
        return "not_recorded"
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return "not_recorded"
    flat = json.dumps(data, sort_keys=True).casefold()
    for key in keys:
        match = re.search(rf'"{re.escape(key.casefold())}"\s*:\s*"?([\w.+-]+)', flat)
        if match:
            return match.group(1)
    return "not_recorded"


def import_partner_genome_evidence(
    db, *, workbook, genome_dir, gtdbtk, checkm2, dataset,
    sheet_name="Bulk Upload Template", min_completeness=90.0,
    max_contamination=5.0, gtdbtk_version="not_recorded",
    checkm2_version="not_recorded", dry_run=False,
) -> list[dict]:
    workbook = Path(workbook).expanduser().resolve()
    genome_dir = Path(genome_dir).expanduser().resolve()
    gtdbtk = Path(gtdbtk).expanduser().resolve()
    checkm2 = Path(checkm2).expanduser().resolve()
    for label, path in (("workbook", workbook), ("genome directory", genome_dir),
                        ("GTDB-Tk summary", gtdbtk), ("CheckM2 report", checkm2)):
        if not path.exists():
            raise FileNotFoundError(f"{label} not found: {path}")

    workbook_rows = load_genome_only_partner_workbook_rows(workbook, sheet_name=sheet_name)
    metadata_by_partner_id = {}
    metadata_by_sequence_id = {}
    for row in workbook_rows:
        culture_id = row.get("partner_culture_id", "").strip()
        if not culture_id:
            continue
        key = _normalise_partner_id(culture_id)
        if key in metadata_by_partner_id:
            raise ValueError(
                f"Partner culture IDs collide after zero-padding normalisation: "
                f"{metadata_by_partner_id[key]['partner_culture_id']!r} and {culture_id!r}"
            )
        row = dict(row)
        row["dataset"] = dataset
        metadata_by_partner_id[key] = row
        metadata_by_sequence_id[row["source_id"]] = row

    if not metadata_by_partner_id:
        raise ValueError("The partner workbook contains no populated partner_culture_id values")

    known = set(db.get_all_ids())
    missing_registered = [r["source_id"] for r in workbook_rows if r["source_id"] not in known]
    if missing_registered:
        raise ValueError(
            "Some workbook isolate IDs are not registered in this project database. "
            "Run import-partner-workbook --genome-only first. Examples: "
            + ", ".join(missing_registered[:10])
        )

    # Build FASTA lookup using partner culture IDs; support .fa/.fna/.fasta and gzip.
    fasta_by_partner_id = {}
    for path in sorted(genome_dir.iterdir()):
        if not path.is_file() or not re.search(r"\.(?:fa|fna|fasta)(?:\.gz)?$", path.name, re.I):
            continue
        key = _normalise_partner_id(path.name)
        if key in fasta_by_partner_id:
            raise ValueError(f"Multiple FASTA files resolve to the same partner culture ID: {path.name}")
        fasta_by_partner_id[key] = path

    gtdb_rows = _read_table(gtdbtk)
    checkm_rows = _read_table(checkm2)
    gtdb_by_key, checkm_by_key = {}, {}
    culture_keys = set(metadata_by_partner_id)

    def resolve_annotation_key(raw_id: str) -> str:
        key = _normalise_partner_id(raw_id)
        if key in culture_keys:
            return key
        # Some partner GTDB-Tk/CheckM2 runs append a stable suffix (e.g.
        # RG471RumEpi) to the submitted culture ID (RG471). Permit this only
        # when exactly one workbook culture ID is an unambiguous prefix.
        candidates = [candidate for candidate in culture_keys
                      if key.startswith(candidate) and key[len(candidate):].isalpha()]
        if len(candidates) == 1:
            return candidates[0]
        return key

    for label, rows, target, id_col in (
        ("GTDB-Tk", gtdb_rows, gtdb_by_key, ("user_genome",)),
        ("CheckM2", checkm_rows, checkm_by_key, ("Name", "name")),
    ):
        for row in rows:
            raw_id = _get(row, *id_col)
            if not raw_id:
                continue
            key = resolve_annotation_key(raw_id)
            if key in target:
                raise ValueError(f"Duplicate {label} identifier after partner-ID resolution: {raw_id}")
            target[key] = row

    metadata_existing = db.get_sequencing_metadata_for_ids()
    for row in workbook_rows:
        registered = metadata_existing.get(row["source_id"], {})
        registered_dataset = str(registered.get("dataset") or "").strip()
        registered_partner = str(registered.get("partner_id") or "").strip()
        if registered_dataset and registered_dataset != dataset:
            raise ValueError(
                f"Dataset mismatch for {row['source_id']}: database has {registered_dataset!r}, "
                f"but --dataset is {dataset!r}"
            )
        if registered_partner and registered_partner.casefold() != row["partner_id"].casefold():
            raise ValueError(
                f"Partner mismatch for {row['source_id']}: database has {registered_partner!r}, "
                f"workbook has {row['partner_id']!r}"
            )
    results, genome_rows, evidence_rows = [], [], []
    source_files = {
        "workbook": str(workbook), "gtdbtk": str(gtdbtk),
        "checkm2": str(checkm2), "genome_dir": str(genome_dir),
    }
    source_hashes = {name: _sha256(Path(path)) for name, path in source_files.items() if name != "genome_dir"}
    gtdbtk_version = gtdbtk_version or "not_recorded"
    checkm2_version = checkm2_version or "not_recorded"

    for key, meta in sorted(metadata_by_partner_id.items(), key=lambda item: item[1]["source_id"]):
        culture_id = meta["partner_culture_id"]
        sequence_id = meta["source_id"]
        fasta = fasta_by_partner_id.get(key)
        gtdb = gtdb_by_key.get(key)
        checkm = checkm_by_key.get(key)
        missing = []
        if fasta is None: missing.append("FASTA_MISSING")
        if gtdb is None: missing.append("GTDBTK_MISSING")
        if checkm is None: missing.append("CHECKM2_MISSING")
        if fasta is None:
            results.append({"sequence_id": sequence_id, "partner_culture_id": culture_id,
                            "genome_id": culture_id, "result": "INCOMPLETE", "detail": ";".join(missing)})
            continue

        completeness = _number(_get(checkm or {}, "Completeness"))
        contamination = _number(_get(checkm or {}, "Contamination"))
        qc_pass = None
        if completeness is not None and contamination is not None:
            qc_pass = completeness >= float(min_completeness) and contamination <= float(max_contamination)
        status = "SEQUENCED" if qc_pass is None else ("GENOME_QC_PASSED" if qc_pass else "GENOME_QC_FAILED")
        taxonomy = _get(gtdb or {}, "classification")
        method = _get(gtdb or {}, "classification_method")
        ani = _number(_get(gtdb or {}, "closest_genome_ani"))
        af = _number(_get(gtdb or {}, "closest_genome_af"))
        closest_ref = _get(gtdb or {}, "closest_genome_reference")
        notes = "; ".join([
            f"partner_culture_id={culture_id}",
            f"gtdbtk_classification_method={method or 'missing'}",
            f"gtdbtk_closest_reference={closest_ref or 'missing'}",
            f"gtdbtk_closest_ani={ani if ani is not None else 'NA'}",
            f"gtdbtk_closest_af={af if af is not None else 'NA'}",
            f"missing_evidence={','.join(missing) if missing else 'none'}",
        ])
        fasta_hash = _sha256(fasta)
        genome_rows.append({
            "genome_id": culture_id, "sequence_id": sequence_id, "accession": "",
            "genome_status": status, "genome_qc_pass": qc_pass,
            "completeness": completeness, "contamination": contamination,
            "gtdb_taxonomy": taxonomy, "ani_cluster": "",
            "genome_path": str(fasta.resolve()), "genome_sha256": fasta_hash,
            "notes": notes, "source_file": str(workbook),
        })
        existing = metadata_existing.get(sequence_id, {})
        project_metadata = dict(existing.get("project_metadata") or {})
        project_metadata.update({
            "genome_available": True, "genome_fasta": str(fasta.resolve()),
            "genome_sha256": fasta_hash, "genome_qc_status": status,
            "genome_evidence_source": "partner_workbook_gtdbtk_checkm2",
            "genome_resequencing_eligible": qc_pass is False,
            "gtdbtk_classification_method": method,
            "gtdbtk_closest_genome_reference": closest_ref,
            "gtdbtk_closest_genome_ani": ani,
            "gtdbtk_closest_genome_af": af,
            "gtdbtk_version": gtdbtk_version,
            "checkm2_version": checkm2_version,
            "partner_culture_id": culture_id,
        })
        if not dry_run:
            db.upsert_sequencing_metadata([{
                "id": sequence_id, "partner_id": existing.get("partner_id", meta["partner_id"]),
                "dataset": existing.get("dataset", dataset),
                "selected_for_sequencing": existing.get("selected_for_sequencing", False),
                "selected_for_wgs": True, "source_id": existing.get("source_id", sequence_id),
                "source_file": str(fasta.resolve()), "raw_selected_value": "genome_available",
                "raw_commitment_value": existing.get("raw_commitment_value", ""),
                "project_metadata": project_metadata,
            }])
            db.update_isolate_status(sequence_id, status, detail=notes, source_file=str(fasta.resolve()))

        for evidence_type, source_path, source_hash, tool, version, raw in (
            ("ASSEMBLY_FASTA", fasta, fasta_hash, "Partner submission", "not_applicable",
             {"path": str(fasta.resolve()), "sha256": fasta_hash, "partner_culture_id": culture_id}),
            ("GTDB_TK", gtdbtk, source_hashes["gtdbtk"], "GTDB-Tk", gtdbtk_version, gtdb or {}),
            ("CHECKM2", checkm2, source_hashes["checkm2"], "CheckM2", checkm2_version, checkm or {}),
        ):
            if raw:
                evidence_rows.append({
                    "sequence_id": sequence_id, "genome_id": culture_id,
                    "evidence_type": evidence_type, "tool_name": tool,
                    "tool_version": version, "source_file": str(source_path),
                    "source_sha256": source_hash, "raw_json": json.dumps(raw, sort_keys=True),
                })
        results.append({
            "sequence_id": sequence_id, "partner_culture_id": culture_id,
            "genome_id": culture_id, "genome_path": str(fasta.resolve()),
            "completeness": completeness, "contamination": contamination,
            "genome_qc_pass": "" if qc_pass is None else str(qc_pass).lower(),
            "gtdbtk_taxonomy": taxonomy, "gtdbtk_method": method,
            "result": "DRY_RUN_READY" if dry_run else "IMPORTED",
            "detail": ";".join(missing) if missing else status,
        })

    matched_keys = set(metadata_by_partner_id)
    for label, lookup, result_name in (
        ("FASTA", fasta_by_partner_id, "UNMATCHED_FASTA"),
        ("GTDB-Tk", gtdb_by_key, "UNMATCHED_GTDBTK"),
        ("CheckM2", checkm_by_key, "UNMATCHED_CHECKM2"),
    ):
        for key in sorted(set(lookup) - matched_keys):
            results.append({"sequence_id": "", "partner_culture_id": "", "genome_id": key,
                            "result": result_name, "detail": f"{label} identifier has no workbook culture-ID mapping"})

    if not dry_run:
        db.upsert_genome_records(genome_rows)
        with db.connect() as conn:
            conn.executemany(
                "INSERT OR IGNORE INTO genome_annotation_evidence "
                "(sequence_id, genome_id, evidence_type, tool_name, tool_version, source_file, source_sha256, raw_json) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                [(r["sequence_id"], r["genome_id"], r["evidence_type"], r["tool_name"], r["tool_version"],
                  r["source_file"], r["source_sha256"], r["raw_json"]) for r in evidence_rows],
            )
        # All registered partner genomes are available, irrespective of QC outcome.
        datasets = sorted({metadata_existing.get(r["sequence_id"], {}).get("dataset") or dataset for r in genome_rows})
        for ds in datasets:
            db.upsert_dataset_role(ds, "candidate", genomes_available=True)
    return results


def write_genome_evidence_report(path: str | Path, rows: list[dict]) -> str:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row}) or ["sequence_id", "result", "detail"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return str(path)

def _parse_qc_pass(value: str, *, row_number: int) -> bool | None:
    """Parse a genome QC flag without silently treating unknown values as failure."""
    text = str(value or "").strip().casefold()
    if text in {"", "na", "n/a", "unknown", "none", "not_assessed", "not assessed", "pending"}:
        return None
    if text in {"1", "true", "yes", "pass", "passed", "qc_pass"}:
        return True
    if text in {"0", "false", "no", "fail", "failed", "qc_fail"}:
        return False
    raise ValueError(
        f"Unrecognised genome_qc_pass value on TSV row {row_number}: {value!r}. "
        "Use a recognised true/false value or leave it blank for unknown."
    )


def import_precomputed_genome_results(db, results_tsv: str | Path, *, dry_run: bool = False) -> dict:
    """Import a per-isolate genome QC/GTDB summary TSV without marker evidence."""
    source = Path(results_tsv).expanduser().resolve()
    rows = _read_table(source)
    if not rows:
        raise ValueError(f"No genome result rows found in {source}")
    known = set(db.get_all_ids())
    dataset_by_id = {}
    with db.connect() as conn:
        dataset_by_id = {str(sid): str(dataset or '') for sid, dataset in conn.execute(
            "SELECT id, dataset FROM sequences"
        )}
    prepared, missing, seen = [], [], set()
    for row_number, row in enumerate(rows, start=2):
        sid = _get(row, "sequence_id")
        genome_id = _get(row, "genome_id", "accession") or sid
        if not sid or not genome_id:
            continue
        if sid not in known:
            missing.append(sid)
            continue
        if genome_id in seen:
            raise ValueError(f"Duplicate genome_id in genome results: {genome_id}")
        seen.add(genome_id)
        qc_pass = _parse_qc_pass(_get(row, "genome_qc_pass"), row_number=row_number)
        default_status = (
            "GENOME_QC_PASSED" if qc_pass is True else
            "GENOME_QC_FAILED" if qc_pass is False else
            "SEQUENCED"
        )
        prepared.append({
            "sequence_id": sid, "genome_id": genome_id,
            "genome_status": _get(row, "genome_status") or default_status,
            "genome_qc_pass": qc_pass,
            "completeness": _number(_get(row, "completeness")),
            "contamination": _number(_get(row, "contamination")),
            "gtdb_taxonomy": _get(row, "gtdb_taxonomy"),
            "ani_cluster": _get(row, "ani_cluster"),
            "genome_path": _get(row, "genome_path"),
            "notes": _get(row, "notes"),
            "source_file": str(source),
        })
    if missing:
        raise ValueError("Genome results contain unregistered sequence IDs: " + ", ".join(missing[:15]))
    if not prepared:
        raise ValueError(f"No valid genome rows found in {source}")
    if not dry_run:
        db.upsert_genome_records(prepared)
        for dataset in sorted({dataset_by_id.get(row["sequence_id"], "") for row in prepared} - {""}):
            db.upsert_dataset_role(dataset, "candidate", genomes_available=True)
    return {
        "source": str(source), "records": len(prepared), "dry_run": bool(dry_run),
        "datasets": sorted({dataset_by_id.get(row["sequence_id"], "") for row in prepared} - {""}),
        "qc_pass": sum(row["genome_qc_pass"] is True for row in prepared),
        "qc_fail": sum(row["genome_qc_pass"] is False for row in prepared),
        "qc_unknown": sum(row["genome_qc_pass"] is None for row in prepared),
    }
