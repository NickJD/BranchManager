"""Evidence provenance and eligibility policy for BranchManager marker sequences.

This module deliberately keeps marker evidence separate from isolate identity,
MALDI observations, and genome taxonomy.  It uses a small SQLite registry so
multiple marker observations (for example Sanger and genome-predicted 16S) can
be retained for one isolate without replacing the active marker in ``sequences``.

The registry is additive: existing ``sequences`` and ``sequence_provenance``
tables remain the operational source for the current active marker.  This
module records a traceable evidence copy and determines which observations may
enter primary or secondary sequence-comparison pools.
"""
from __future__ import annotations

import csv
import hashlib
import json
import logging
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from branchmanager.utils.fasta import read_fasta

logger = logging.getLogger(__name__)

SOURCE_TYPES = {
    "SANGER_AB1",
    "PARTNER_FASTA",
    "PACBIO_HIFI",
    "GENOME_BARRNAP",
    "LEGACY_BASELINE",
    "LEGACY_UNKNOWN",
}
EVIDENCE_TIERS = {"HIGH", "MODERATE", "LOW"}
MANUAL_STATES = {"NOT_REQUIRED", "PENDING", "APPROVED", "REJECTED", "NOT_REVIEWED"}

SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS marker_evidence (
    evidence_id INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id TEXT NOT NULL,
    dataset TEXT NOT NULL DEFAULT '',
    sequence TEXT NOT NULL,
    sequence_sha256 TEXT NOT NULL,
    source_type TEXT NOT NULL,
    evidence_tier TEXT NOT NULL,
    qc_class TEXT NOT NULL DEFAULT 'QUALITY_UNVERIFIED',
    manual_review_status TEXT NOT NULL DEFAULT 'NOT_REVIEWED',
    review_decision TEXT NOT NULL DEFAULT '',
    reviewer TEXT NOT NULL DEFAULT '',
    review_notes TEXT NOT NULL DEFAULT '',
    taxonomy_source TEXT NOT NULL DEFAULT '',
    taxonomy TEXT NOT NULL DEFAULT '',
    taxonomy_rank TEXT NOT NULL DEFAULT '',
    taxonomy_identity REAL,
    query_coverage REAL,
    target_coverage REAL,
    genome_id TEXT NOT NULL DEFAULT '',
    genome_taxonomy TEXT NOT NULL DEFAULT '',
    genome_16s_concordance TEXT NOT NULL DEFAULT '',
    source_file TEXT NOT NULL DEFAULT '',
    source_sha256 TEXT NOT NULL DEFAULT '',
    active_for_selection INTEGER NOT NULL DEFAULT 0,
    eligible_primary INTEGER NOT NULL DEFAULT 0,
    created_at TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_marker_evidence_sequence ON marker_evidence(sequence_id);
CREATE INDEX IF NOT EXISTS idx_marker_evidence_pool ON marker_evidence(eligible_primary, active_for_selection, evidence_tier);
CREATE INDEX IF NOT EXISTS idx_marker_evidence_source ON marker_evidence(source_type, evidence_tier);
CREATE TABLE IF NOT EXISTS marker_evidence_meta (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
"""


def _now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def sha256_text(sequence: str) -> str:
    return hashlib.sha256(str(sequence).upper().encode("ascii", errors="replace")).hexdigest()


def ensure_marker_evidence_schema(db_or_connection) -> None:
    """Create/extend the evidence registry without changing legacy tables."""
    def apply(conn):
        conn.executescript(SCHEMA_SQL)
        columns = _table_columns(conn, 'marker_evidence')
        if 'genome_16s_concordance' not in columns:
            conn.execute(
                "ALTER TABLE marker_evidence ADD COLUMN genome_16s_concordance TEXT NOT NULL DEFAULT ''"
            )
    if hasattr(db_or_connection, "connect"):
        with db_or_connection.connect() as conn:
            apply(conn)
    else:
        apply(db_or_connection)


def _table_columns(conn: sqlite3.Connection, table: str) -> set[str]:
    try:
        return {str(row[1]) for row in conn.execute(f"PRAGMA table_info({table})").fetchall()}
    except sqlite3.Error:
        return set()


def _registry_has_rows(conn: sqlite3.Connection) -> bool:
    try:
        return conn.execute("SELECT 1 FROM marker_evidence LIMIT 1").fetchone() is not None
    except sqlite3.Error:
        return False


def _derive_legacy_tier(role: str, provenance: dict) -> tuple[str, str, str, int]:
    """Conservative one-time mapping for pre-registry project sequences."""
    if role == "baseline":
        return "LEGACY_BASELINE", "HIGH", "CURATED_BASELINE", 1
    qc = str(provenance.get("marker_qc_class") or "QUALITY_UNVERIFIED").upper()
    review = str(provenance.get("manual_review_status") or "NOT_REVIEWED").upper()
    decision = str(provenance.get("review_decision") or "").upper()
    if qc == "PASS_HIGH_CONFIDENCE":
        return "LEGACY_UNKNOWN", "HIGH", qc, 1
    if qc == "PASS_WITH_WARNINGS" and review == "APPROVED":
        return "LEGACY_UNKNOWN", "MODERATE", qc, 1
    if review == "APPROVED" and decision == "APPROVE_FAIL_QC":
        return "LEGACY_UNKNOWN", "MODERATE", qc, 1
    if review == "APPROVED" and decision in {"APPROVE_UNVERIFIED", "APPROVE"}:
        return "LEGACY_UNKNOWN", "MODERATE", qc, 1
    return "LEGACY_UNKNOWN", "LOW", qc, 0


def backfill_legacy_marker_evidence(db, *, force: bool = False) -> dict[str, int]:
    """Backfill existing sequence rows once; unknown partner provenance defaults LOW.

    Curated registered baseline datasets are kept as HIGH baseline evidence.
    Partner sequences are admitted only when existing provenance shows a valid
    PASS_HIGH_CONFIDENCE marker or an explicit manual approval. This avoids
    silently promoting old unverified FASTA records.
    """
    ensure_marker_evidence_schema(db)
    counts = {"baseline_high": 0, "partner_high": 0, "partner_moderate": 0, "partner_low": 0, "skipped": 0}
    with db.connect() as conn:
        if not force:
            done = conn.execute(
                "SELECT value FROM marker_evidence_meta WHERE key='legacy_backfill_v1'"
            ).fetchone()
            if done and done[0] == "complete":
                return counts
        seq_rows = conn.execute(
            "SELECT id, sequence, COALESCE(dataset, '') FROM sequences "
            "WHERE sequence IS NOT NULL AND sequence != '' ORDER BY id"
        ).fetchall()
        role_map = {}
        try:
            role_map = {
                str(dataset or ""): str(role or "")
                for dataset, role in conn.execute("SELECT dataset, role FROM dataset_roles").fetchall()
            }
        except sqlite3.Error:
            pass
        provenance_map = {}
        provenance_cols = _table_columns(conn, "sequence_provenance")
        if provenance_cols and "sequence_id" in provenance_cols:
            select_cols = ["sequence_id"]
            for field in ("marker_qc_class", "manual_review_status", "review_decision", "marker_qc_reasons"):
                if field in provenance_cols:
                    select_cols.append(field)
            try:
                for values in conn.execute(
                    "SELECT " + ",".join(select_cols) + " FROM sequence_provenance"
                ).fetchall():
                    provenance_map[str(values[0])] = dict(zip(select_cols[1:], values[1:]))
            except sqlite3.Error:
                pass
        now = _now()
        for sequence_id, sequence, dataset in seq_rows:
            role = role_map.get(str(dataset or ""), "")
            source, tier, qc_class, eligible = _derive_legacy_tier(role, provenance_map.get(str(sequence_id), {}))
            if tier == "HIGH" and source == "LEGACY_BASELINE":
                counts["baseline_high"] += 1
            elif tier == "HIGH":
                counts["partner_high"] += 1
            elif tier == "MODERATE":
                counts["partner_moderate"] += 1
            else:
                counts["partner_low"] += 1
            digest = sha256_text(sequence)
            existing = conn.execute(
                "SELECT evidence_id FROM marker_evidence WHERE sequence_id=? AND sequence_sha256=? AND source_type=?",
                (str(sequence_id), digest, source),
            ).fetchone()
            if existing:
                counts["skipped"] += 1
                continue
            conn.execute(
                "INSERT INTO marker_evidence (sequence_id,dataset,sequence,sequence_sha256,source_type,"
                "evidence_tier,qc_class,manual_review_status,active_for_selection,eligible_primary,created_at) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                (str(sequence_id), str(dataset or ""), str(sequence).upper(), digest, source, tier,
                 qc_class, str(provenance_map.get(str(sequence_id), {}).get("manual_review_status") or "NOT_REVIEWED"),
                 1, int(eligible), now),
            )
        conn.execute(
            "INSERT INTO marker_evidence_meta(key,value) VALUES('legacy_backfill_v1','complete') "
            "ON CONFLICT(key) DO UPDATE SET value=excluded.value"
        )
    return counts


def _tier_for_marker(source_type: str, qc_class: str, manual_status: str, review_decision: str) -> tuple[str, int]:
    source = source_type.upper()
    qc = qc_class.upper()
    review = manual_status.upper()
    decision = review_decision.upper()
    if source == "GENOME_BARRNAP":
        # Genome-predicted 16S is supplementary context, never a HIGH Sanger-equivalent marker.
        return ("LOW", 0) if qc in {"FAIL_QC", "GENOME_DERIVED_FAIL"} else ("MODERATE", 0)
    if qc == "PASS_HIGH_CONFIDENCE":
        return "HIGH", 1
    if qc == "PASS_WITH_WARNINGS" and review == "APPROVED":
        return "MODERATE", 1
    if review == "APPROVED" and decision == "APPROVE_FAIL_QC":
        return "MODERATE", 1
    if review == "APPROVED" and decision in {"APPROVE_UNVERIFIED", "APPROVE"}:
        return "MODERATE", 1
    return "LOW", 0


def upsert_active_marker_evidence(
    db,
    records: Iterable[tuple[str, str]],
    provenance_rows: Iterable[dict],
    *,
    source_type: str,
    dataset: str,
    source_file: str = "",
) -> list[dict]:
    """Persist the marker currently used for selection for each processed isolate."""
    ensure_marker_evidence_schema(db)
    source_type = str(source_type or "PARTNER_FASTA").upper()
    allowed_active_sources = SOURCE_TYPES - {"LEGACY_BASELINE", "LEGACY_UNKNOWN", "GENOME_BARRNAP"}
    if source_type not in allowed_active_sources and source_type != "MIXED":
        raise ValueError(f"Unsupported active marker source type: {source_type}")
    provenance = {str(row.get("sequence_id")): dict(row) for row in provenance_rows if row.get("sequence_id")}
    now = _now()
    inserted = []
    with db.connect() as conn:
        for sequence_id, sequence in records:
            sequence_id = str(sequence_id)
            sequence = str(sequence).upper()
            p = provenance.get(sequence_id, {})
            record_source_type = str(p.get("marker_source_type") or ("PARTNER_FASTA" if source_type == "MIXED" else source_type)).upper()
            if record_source_type not in allowed_active_sources:
                raise ValueError(f"Unsupported active marker source type {record_source_type!r} for {sequence_id}")
            qc = str(p.get("marker_qc_class") or "QUALITY_UNVERIFIED").upper()
            manual = str(p.get("manual_review_status") or "NOT_REVIEWED").upper()
            decision = str(p.get("review_decision") or p.get("decision") or "").upper()
            tier, eligible = _tier_for_marker(record_source_type, qc, manual, decision)
            digest = sha256_text(sequence)
            conn.execute(
                "UPDATE marker_evidence SET active_for_selection=0 WHERE sequence_id=? AND active_for_selection=1",
                (sequence_id,),
            )
            conn.execute(
                "INSERT INTO marker_evidence (sequence_id,dataset,sequence,sequence_sha256,source_type,"
                "evidence_tier,qc_class,manual_review_status,review_decision,reviewer,review_notes,"
                "source_file,active_for_selection,eligible_primary,created_at) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (sequence_id, str(dataset or ""), sequence, digest, record_source_type, tier, qc, manual,
                 decision, str(p.get("reviewer") or ""), str(p.get("review_notes") or p.get("marker_qc_reasons") or ""),
                 str(source_file or p.get("source_sequence_file") or ""), 1, int(eligible), now),
            )
            inserted.append({
                "sequence_id": sequence_id, "source_type": record_source_type, "evidence_tier": tier,
                "qc_class": qc, "manual_review_status": manual, "eligible_primary": bool(eligible),
                "sequence_sha256": digest,
            })
    return inserted


def register_supplemental_marker_fasta(
    db,
    fasta_path: str | Path,
    metadata_path: str | Path,
    *,
    default_source_type: str = "GENOME_BARRNAP",
    default_evidence_tier: str = "MODERATE",
) -> list[dict]:
    """Register supplemental marker observations (normally Barrnap 16S) by isolate ID.

    Metadata TSV columns: sequence_id, optional fasta_id, dataset, source_type,
    evidence_tier, qc_class, genome_id, genome_taxonomy, genome_16s_concordance,
    taxonomy, taxonomy_source, taxonomy_identity, query_coverage, target_coverage,
    source_file, source_sha256. FASTA headers may use the genome ID as the first
    token when `fasta_id` maps it to the stable BranchManager `sequence_id`.
    """
    ensure_marker_evidence_schema(db)
    with open(metadata_path, newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_fasta_id = {}
    for row in rows:
        sid = str(row.get("sequence_id") or row.get("SequenceID") or row.get("id") or "").strip()
        fasta_id = str(row.get("fasta_id") or row.get("genome_id") or sid).strip()
        if not sid:
            continue
        if not fasta_id:
            raise ValueError(f"Supplemental marker metadata has no fasta_id/genome_id for {sid!r}")
        if fasta_id in by_fasta_id:
            raise ValueError(f"Supplemental marker metadata contains duplicate fasta_id {fasta_id!r}")
        by_fasta_id[fasta_id] = row
    fasta_records = list(read_fasta(str(fasta_path)))
    seen = set()
    results = []
    with db.connect() as conn:
        existing_ids = {str(row[0]) for row in conn.execute("SELECT id FROM sequences").fetchall()}
        for header, sequence in fasta_records:
            fasta_id = str(header).split()[0]
            meta = by_fasta_id.get(fasta_id)
            if meta is None:
                raise ValueError(f"Supplemental marker FASTA ID {fasta_id!r} has no metadata mapping row")
            sid = str(meta.get("sequence_id") or meta.get("SequenceID") or meta.get("id") or "").strip()
            if not sid:
                raise ValueError(f"Supplemental marker mapping for FASTA ID {fasta_id!r} has no sequence_id")
            if sid in seen:
                raise ValueError(f"Multiple supplemental FASTA records map to the same isolate ID {sid!r}")
            seen.add(sid)
            if sid not in existing_ids:
                raise ValueError(f"Supplemental marker isolate {sid!r} is not present in this project database")
            source = str(meta.get("source_type") or default_source_type).upper()
            tier = str(meta.get("evidence_tier") or default_evidence_tier).upper()
            if source != "GENOME_BARRNAP":
                raise ValueError(
                    f"Supplemental marker import is reserved for GENOME_BARRNAP evidence; got {source!r} for {sid}"
                )
            if tier not in EVIDENCE_TIERS:
                raise ValueError(f"Unsupported evidence_tier {tier!r} for {sid}")
            if source == "GENOME_BARRNAP" and tier == "HIGH":
                raise ValueError(f"GENOME_BARRNAP evidence cannot be marked HIGH ({sid}); use MODERATE or LOW")
            qc = str(meta.get("qc_class") or ("GENOME_DERIVED" if source == "GENOME_BARRNAP" else "QUALITY_UNVERIFIED")).upper()
            eligible = int(source != "GENOME_BARRNAP" and tier in {"HIGH", "MODERATE"})
            digest = sha256_text(sequence)
            previous = conn.execute(
                "SELECT evidence_id FROM marker_evidence WHERE sequence_id=? AND sequence_sha256=? AND source_type=?",
                (sid, digest, source),
            ).fetchone()
            if previous:
                results.append({"sequence_id": sid, "result": "ALREADY_PRESENT", "evidence_id": previous[0]})
                continue
            conn.execute(
                "INSERT INTO marker_evidence (sequence_id,dataset,sequence,sequence_sha256,source_type,"
                "evidence_tier,qc_class,manual_review_status,review_decision,taxonomy_source,taxonomy,"
                "taxonomy_rank,taxonomy_identity,query_coverage,target_coverage,genome_id,genome_taxonomy,"
                "genome_16s_concordance,source_file,source_sha256,active_for_selection,eligible_primary,created_at) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (sid, str(meta.get("dataset") or ""), str(sequence).upper(), digest, source, tier, qc,
                 "NOT_REQUIRED", "", str(meta.get("taxonomy_source") or "VSEARCH_16S_GTDB"),
                 str(meta.get("taxonomy") or ""), str(meta.get("taxonomy_rank") or ""),
                 _maybe_float(meta.get("taxonomy_identity")), _maybe_float(meta.get("query_coverage")),
                 _maybe_float(meta.get("target_coverage")), str(meta.get("genome_id") or fasta_id),
                 str(meta.get("genome_taxonomy") or ""), str(meta.get("genome_16s_concordance") or ""),
                 str(meta.get("source_file") or fasta_path), str(meta.get("source_sha256") or ""),
                 0, eligible, _now()),
            )
            results.append({"sequence_id": sid, "result": "IMPORTED", "source_type": source,
                            "evidence_tier": tier, "sequence_sha256": digest})
    return results


def _maybe_float(value):
    try:
        if value in (None, "", "NA", "None"):
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def query_marker_evidence_pool(
    conn: sqlite3.Connection,
    *,
    pool: str = "primary",
    run_dataset: str | None = None,
    include_datasets: Iterable[str] | None = None,
    include_ids: Iterable[str] | None = None,
    exclude_datasets: Iterable[str] | None = None,
) -> tuple[list[tuple[str, str]], list[str]]:
    """Return eligible sequence records from the evidence registry."""
    include_datasets = [str(x) for x in (include_datasets or []) if str(x)]
    include_ids = [str(x) for x in (include_ids or []) if str(x)]
    exclude_datasets = [str(x) for x in (exclude_datasets or []) if str(x)]
    clauses = []
    params: list[object] = []
    if pool == "primary":
        clauses.extend(["eligible_primary=1", "active_for_selection=1", "evidence_tier IN ('HIGH','MODERATE')"])
    elif pool == "genome_derived":
        clauses.extend(["source_type='GENOME_BARRNAP'", "evidence_tier='MODERATE'", "qc_class='GENOME_DERIVED'"])
    else:
        raise ValueError(f"Unknown marker evidence pool: {pool}")
    if run_dataset:
        clauses.append("(dataset != ? OR dataset IS NULL)")
        params.append(run_dataset)
    if include_datasets or include_ids:
        sub = []
        if include_datasets:
            sub.append("dataset IN (" + ",".join("?" for _ in include_datasets) + ")")
            params.extend(include_datasets)
        if include_ids:
            sub.append("sequence_id IN (" + ",".join("?" for _ in include_ids) + ")")
            params.extend(include_ids)
        clauses.append("(" + " OR ".join(sub) + ")")
    if exclude_datasets:
        clauses.append("(dataset NOT IN (" + ",".join("?" for _ in exclude_datasets) + ") OR dataset IS NULL)")
        params.extend(exclude_datasets)
    sql = (
        "SELECT sequence_id, sequence, dataset FROM marker_evidence WHERE " + " AND ".join(clauses)
        + " ORDER BY sequence_id, created_at DESC, evidence_id DESC"
    )
    # Keep one sequence per isolate in each pool. Registry rows may retain historical marker versions.
    rows = conn.execute(sql, tuple(params)).fetchall()
    records_by_id = {}
    datasets = set()
    for sid, sequence, dataset in rows:
        if sid and sequence:
            records_by_id.setdefault(str(sid), str(sequence))
            if dataset:
                datasets.add(str(dataset))
    return sorted(records_by_id.items()), sorted(datasets)


def write_marker_evidence_report(db, output_path: str | Path) -> str:
    ensure_marker_evidence_schema(db)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with db.connect() as conn:
        rows = conn.execute(
            "SELECT sequence_id,dataset,source_type,evidence_tier,qc_class,manual_review_status,"
            "review_decision,reviewer,review_notes,taxonomy_source,taxonomy,taxonomy_rank,taxonomy_identity,"
            "query_coverage,target_coverage,genome_id,genome_taxonomy,genome_16s_concordance,"
            "sequence_sha256,source_sha256,active_for_selection,eligible_primary,source_file,created_at "
            "FROM marker_evidence ORDER BY sequence_id,created_at,evidence_id"
        ).fetchall()
    fields = ["sequence_id","dataset","source_type","evidence_tier","qc_class","manual_review_status",
              "review_decision","reviewer","review_notes","taxonomy_source","taxonomy","taxonomy_rank",
              "taxonomy_identity","query_coverage","target_coverage","genome_id","genome_taxonomy",
              "genome_16s_concordance","sequence_sha256","source_sha256","active_for_selection",
              "eligible_primary","source_file","created_at"]
    with open(output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(fields)
        writer.writerows(rows)
    return str(output)


def register_baseline_marker_records(db, fasta_path: str | Path, *, dataset: str) -> int:
    """Register newly loaded Filing Cabinet records as curated HIGH baseline markers."""
    ensure_marker_evidence_schema(db)
    now = _now()
    records = list(read_fasta(str(fasta_path)))
    added = 0
    with db.connect() as conn:
        for sequence_id, sequence in records:
            sid = str(sequence_id).split()[0]
            seq = str(sequence).upper()
            digest = sha256_text(seq)
            existing = conn.execute(
                "SELECT evidence_id FROM marker_evidence WHERE sequence_id=? AND sequence_sha256=? AND source_type='LEGACY_BASELINE'",
                (sid, digest),
            ).fetchone()
            if existing:
                continue
            conn.execute(
                "UPDATE marker_evidence SET active_for_selection=0 WHERE sequence_id=? AND active_for_selection=1",
                (sid,),
            )
            conn.execute(
                "INSERT INTO marker_evidence (sequence_id,dataset,sequence,sequence_sha256,source_type,"
                "evidence_tier,qc_class,manual_review_status,active_for_selection,eligible_primary,created_at) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                (sid, str(dataset), seq, digest, 'LEGACY_BASELINE', 'HIGH', 'CURATED_BASELINE',
                 'NOT_REQUIRED', 1, 1, now),
            )
            added += 1
    return added


def get_active_marker_evidence(db, sequence_ids: Iterable[str]) -> dict[str, dict]:
    """Return the latest active marker evidence row per isolate."""
    ids = list(dict.fromkeys(str(x) for x in sequence_ids if str(x)))
    if not ids:
        return {}
    ensure_marker_evidence_schema(db)
    output = {}
    with db.connect() as conn:
        for start in range(0, len(ids), 800):
            chunk = ids[start:start + 800]
            rows = conn.execute(
                "SELECT sequence_id,source_type,evidence_tier,qc_class,manual_review_status,"
                "review_decision,reviewer,review_notes,sequence_sha256 "
                "FROM marker_evidence WHERE active_for_selection=1 AND sequence_id IN (" +
                ",".join("?" for _ in chunk) + ") ORDER BY created_at DESC,evidence_id DESC",
                tuple(chunk),
            ).fetchall()
            for row in rows:
                output.setdefault(str(row[0]), {
                    'marker_source_type': row[1], 'marker_evidence_tier': row[2],
                    'marker_qc_class': row[3], 'marker_manual_review_status': row[4],
                    'marker_review_decision': row[5], 'marker_review_reviewer': row[6],
                    'marker_review_notes': row[7], 'marker_sequence_sha256': row[8],
                })
    return output
