"""Additive, read-only ground-truth and drift validation for BranchManager.

This module deliberately operates on completed run directories. It does not
modify workflow outputs, project databases, or biological decisions.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
import shutil
import sqlite3
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

VALIDATION_SCHEMA_VERSION = 1
REFERENCE_SCHEMA_VERSION = 1

# Core decision-bearing outputs. The relative source paths are retained in the
# reference manifest; the validator can therefore handle the same workflow
# layout without hard-coding a particular run identifier.
DEFAULT_OUTPUTS = {
    "sequence_assessment.tsv": {
        "category": "assessment",
        "key": ["ID"],
        "ignore_columns": [
            "LocalNeighbourhoodFigure", "LocalPairwisePidentTable",
            "MarkerSourceManifest",
        ],
    },
    "selection_summary.tsv": {
        "category": "selection",
        "key": ["SequenceID"],
        "ignore_columns": [],
    },
    "sequencing_sets.tsv": {
        "category": "selection",
        "key": ["CandidateID"],
        "ignore_columns": ["LocalTreeFigure"],
    },
    "novelty_metrics.tsv": {
        "category": "novelty",
        "key": ["ID"],
        "ignore_columns": [],
    },
    "cohort_review_application_audit.tsv": {
        "category": "marker_review",
        "key": ["sequence_id"],
        "ignore_columns": [],
    },
    "cohort_joint_marker_qc.tsv": {
        "category": "marker_qc",
        "key": ["SequenceID"],
        "ignore_columns": [],
    },
}

DATABASE_TABLES = (
    "sequences", "genome_records", "isolate_status",
    "sequencing_metadata", "marker_evidence", "marker_review_decisions",
)
REQUIRED_OUTPUTS = {"sequence_assessment.tsv", "selection_summary.tsv"}
OPTIONAL_OUTPUTS = set(DEFAULT_OUTPUTS) - REQUIRED_OUTPUTS
APPROVAL_RECORD = "ground_truth_approval.json"


NULL_VALUES = {"", "na", "n/a", "nan", "none", "null"}


def _now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _path_sha256(path: Path) -> str:
    """Match RunManifest.file_record hashing for files and directory inputs."""
    if path.is_file():
        return _sha256(path)
    if path.is_dir():
        digest = hashlib.sha256()
        for child in sorted(p for p in path.rglob("*") if p.is_file()):
            relative = child.relative_to(path).as_posix()
            digest.update(relative.encode("utf-8"))
            digest.update(b"\0")
            digest.update(_sha256(child).encode("ascii"))
            digest.update(b"\n")
        return digest.hexdigest()
    raise FileNotFoundError(path)


def _json_write(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return value


def _find_file(root: Path, name: str) -> Path | None:
    """Resolve a named output without silently choosing between ambiguous copies."""
    direct = root / name
    if direct.is_file():
        return direct
    matches = sorted(p for p in root.rglob(name) if p.is_file())
    if not matches:
        return None
    if name in {"sequence_assessment.tsv", "selection_summary.tsv",
                "sequencing_sets.tsv", "novelty_metrics.tsv"}:
        preferred = [p for p in matches if "03_performance_review_hiring_panel" in str(p)]
    elif name in {"cohort_review_application_audit.tsv", "cohort_joint_marker_qc.tsv"}:
        preferred = [p for p in matches if "01_joint_inputs" in str(p)]
    else:
        preferred = []
    candidates = preferred or matches
    if len(candidates) > 1:
        raise ValueError(
            f"Ambiguous output {name!r}: found {len(candidates)} candidates: "
            + "; ".join(str(p) for p in candidates[:8])
        )
    return candidates[0]


def _performance_manifest(root: Path) -> Path | None:
    candidates = sorted(root.rglob("run_manifest.json"))
    for path in candidates:
        try:
            data = _read_json(path)
        except Exception:
            continue
        if str(data.get("workflow", "")).lower() == "performance_review":
            return path
    # Assistant cohort manifests may point to the nested performance review
    # through their stages; do not silently substitute unrelated stage manifests.
    return None


def _load_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"TSV has no header: {path}")
        headers = [str(x) for x in reader.fieldnames]
        if any(not header.strip() for header in headers):
            raise ValueError(f"TSV contains a blank column name: {path}")
        if len(headers) != len(set(headers)):
            raise ValueError(f"TSV contains duplicate column names: {path}")
        rows = [{k: (v if v is not None else "") for k, v in row.items()} for row in reader]
    return headers, rows


def _normalise(value: Any) -> str:
    text = "" if value is None else str(value).strip()
    if text.lower() in NULL_VALUES:
        return ""
    return text


def _equivalent(a: Any, b: Any, tolerance: float) -> bool:
    left, right = _normalise(a), _normalise(b)
    if left == right:
        return True
    if left and right:
        try:
            x, y = float(left), float(right)
            if math.isfinite(x) and math.isfinite(y):
                return abs(x - y) <= tolerance
        except (TypeError, ValueError, OverflowError):
            pass
    return False


def _table_metadata(path: Path, key: list[str]) -> dict[str, Any]:
    headers, rows = _load_tsv(path)
    missing = [k for k in key if k not in headers]
    if missing:
        raise ValueError(f"{path.name}: missing key column(s): {missing}")
    ids = [tuple(_normalise(row.get(k, "")) for k in key) for row in rows]
    counts = Counter(ids)
    return {
        "rows": len(rows),
        "columns": headers,
        "key": key,
        "duplicate_key_count": sum(n - 1 for n in counts.values() if n > 1),
        "blank_key_count": sum(1 for key_value in ids if any(not part for part in key_value)),
    }


def _resolve_database(root: Path, database: str | Path | None = None) -> Path | None:
    if database is not None:
        path = Path(database).expanduser()
        if not path.is_absolute():
            cwd_candidate = path.resolve()
            path = cwd_candidate if cwd_candidate.is_file() else (root / path).resolve()
        path = path.resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Validation database does not exist: {path}")
        return path
    candidates = []
    for pattern in ("project.sqlite", "*.sqlite", "*.db"):
        candidates.extend(p.resolve() for p in root.rglob(pattern) if p.is_file())
    unique = sorted(set(candidates))
    if not unique:
        return None
    if len(unique) > 1:
        raise ValueError(
            "Multiple SQLite databases found; specify --db explicitly: "
            + "; ".join(str(p) for p in unique[:8])
        )
    return unique[0]


def _database_state(root: Path, database: str | Path | None = None) -> dict[str, Any] | None:
    db_path = _resolve_database(root, database)
    if db_path is None:
        return None
    uri = db_path.as_uri() + "?mode=ro"
    conn = sqlite3.connect(uri, uri=True, timeout=10)
    try:
        integrity = conn.execute("PRAGMA integrity_check").fetchone()[0]
        tables = {
            row[0] for row in conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'"
            )
        }
        counts: dict[str, int] = {}
        fingerprints: dict[str, Any] = {}
        for table in DATABASE_TABLES:
            if table not in tables:
                continue
            quoted = '"' + table.replace('"', '""') + '"'
            info = conn.execute(f"PRAGMA table_info({quoted})").fetchall()
            columns = [row[1] for row in info]
            primary_key = [row[1] for row in sorted((r for r in info if r[5]), key=lambda r: r[5])]
            records: dict[str, str] = {}
            row_hashes: list[str] = []
            for row in conn.execute(f"SELECT * FROM {quoted}"):
                payload = json.dumps(list(row), ensure_ascii=False, separators=(",", ":"), default=str)
                digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()
                row_hashes.append(digest)
                if primary_key:
                    values = [row[columns.index(col)] for col in primary_key]
                    record_key = json.dumps(values, ensure_ascii=False, separators=(",", ":"), default=str)
                    if record_key in records:
                        raise ValueError(f"Duplicate primary key in database table {table}: {record_key}")
                    records[record_key] = digest
            counts[table] = len(row_hashes)
            fingerprints[table] = {
                "columns": columns,
                "primary_key": primary_key,
                "rows": len(row_hashes),
                "records": dict(sorted(records.items())) if primary_key else None,
                "unordered_row_hashes": sorted(row_hashes) if not primary_key else None,
            }
        schema_version = None
        if "schema_migrations" in tables:
            schema_version = conn.execute(
                "SELECT MAX(version) FROM schema_migrations"
            ).fetchone()[0]
        foreign_key_violations = len(conn.execute("PRAGMA foreign_key_check").fetchall())
        try:
            display_path = str(db_path.relative_to(root))
        except ValueError:
            display_path = str(db_path)
        return {
            "path": display_path,
            "sha256": _sha256(db_path),
            "integrity_check": integrity,
            "schema_version": schema_version,
            "table_counts": counts,
            "table_fingerprints": fingerprints,
            "foreign_key_violations": foreign_key_violations,
        }
    finally:
        conn.close()


def _manifest_inputs(perf_manifest: dict[str, Any] | None) -> dict[str, dict[str, Any]]:
    """Retain every input record, including repeated roles such as alternate_reference."""
    if not perf_manifest:
        return {}
    grouped: dict[str, list[dict[str, Any]]] = {}
    for item in perf_manifest.get("inputs", []) or []:
        role = str(item.get("role") or item.get("path") or "").strip()
        if not role:
            continue
        grouped.setdefault(role, []).append({
            "path": str(item.get("path") or ""),
            "sha256": str(item.get("sha256") or ""),
            "required": bool(item.get("required", True)),
            "exists": bool(item.get("exists", False)),
            "kind": str(item.get("kind") or ""),
        })
    result: dict[str, dict[str, Any]] = {}
    for role, records in sorted(grouped.items()):
        records.sort(key=lambda item: (item.get("path", ""), item.get("sha256", "")))
        for index, record in enumerate(records, start=1):
            key = role if len(records) == 1 else f"{role}#{index}"
            result[key] = record
    return result


def _historical_database_evidence(
    role: str, item: dict[str, Any], perf_manifest: dict[str, Any] | None,
) -> dict[str, Any] | None:
    """Return paired before/after manifest evidence for a removed cohort DB.

    This exception is deliberately narrow: only the temporary project database
    captured as ``project_database_before`` may be verified historically. The
    performance-review manifest must record both the input and output state,
    with matching paths, non-empty SHA256 values and ``exists=True``.
    """
    if role != "project_database_before" or not perf_manifest:
        return None
    path = str(item.get("path") or "")
    before = [
        record for record in (perf_manifest.get("inputs") or [])
        if record.get("role") == role
        and str(record.get("path") or "") == path
        and str(record.get("sha256") or "") == str(item.get("sha256") or "")
        and record.get("exists") is True
    ]
    after = [
        record for record in (perf_manifest.get("outputs") or [])
        if record.get("role") == "project_database_after"
        and str(record.get("path") or "") == path
        and bool(record.get("sha256"))
        and record.get("exists") is True
    ]
    if len(before) != 1 or len(after) != 1:
        return None
    before_sha = str(before[0].get("sha256") or "")
    after_sha = str(after[0].get("sha256") or "")
    if any(
        len(value) != 64 or any(character not in "0123456789abcdefABCDEF" for character in value)
        for value in (before_sha, after_sha)
    ):
        return None
    return {
        "before": {
            "role": role, "path": path,
            "sha256": before_sha,
            "exists_at_workflow_time": True,
        },
        "after": {
            "role": "project_database_after", "path": path,
            "sha256": after_sha,
            "exists_at_workflow_time": True,
        },
    }


def _verify_input_provenance(
    inputs: dict[str, dict[str, Any]],
    perf_manifest: dict[str, Any] | None = None,
) -> tuple[str, list[str]]:
    if not inputs:
        return "UNAVAILABLE", ["No input records were present in the Performance Review manifest"]
    problems = []
    unavailable = []
    historically_verified = []
    for role, item in inputs.items():
        path_text = item.get("path", "")
        expected = item.get("sha256", "")
        if not path_text or not expected:
            if item.get("required", True):
                unavailable.append(f"{role}: missing path or recorded SHA256")
            continue
        path = Path(path_text).expanduser()
        if not path.exists() or not (path.is_file() or path.is_dir()):
            if not item.get("required", True):
                continue
            evidence = _historical_database_evidence(role, item, perf_manifest)
            if evidence is not None:
                historically_verified.append(
                    f"{role}: file no longer exists; before/after database checksums "
                    "are recorded in the Performance Review manifest"
                )
            else:
                unavailable.append(f"{role}: input unavailable ({path})")
            continue
        observed = _path_sha256(path)
        if observed != expected:
            problems.append(f"{role}: recorded SHA256 does not match current file ({path})")
    if problems:
        return "INVALID", problems + unavailable + historically_verified
    if unavailable:
        return "UNAVAILABLE", unavailable + historically_verified
    if historically_verified:
        return "VERIFIED_FROM_RUN_MANIFEST", historically_verified
    return "VERIFIED", []


def _manifest_environment(perf_manifest: dict[str, Any] | None) -> dict[str, Any]:
    if not perf_manifest:
        return {}
    return {
        "branchmanager_version": perf_manifest.get("branchmanager_version"),
        "workflow_version": perf_manifest.get("workflow_version"),
        "environment": perf_manifest.get("environment", {}),
        "external_tools": perf_manifest.get("external_tools", {}),
        "policies": perf_manifest.get("policies", {}),
        "reference_identity": perf_manifest.get("reference_identity", {}),
    }


def create_reference(run_dir: str | Path, reference_dir: str | Path, *,
                     label: str = "", tolerance: float = 1e-6,
                     database: str | Path | None = None) -> Path:
    """Create an immutable-by-default PROVISIONAL snapshot from a completed run."""
    root = Path(run_dir).expanduser().resolve()
    destination = Path(reference_dir).expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"Run directory does not exist: {root}")
    if destination.exists():
        raise FileExistsError(f"Reference destination already exists; refusing to overwrite: {destination}")

    expected_dir = destination / "expected"
    expected_dir.mkdir(parents=True)
    outputs = {}
    for name, spec in DEFAULT_OUTPUTS.items():
        try:
            source = _find_file(root, name)
            if source is None:
                continue
            rel = source.relative_to(root)
            target_name = name
            target = expected_dir / target_name
            shutil.copy2(source, target)
            meta = _table_metadata(target, spec["key"])
        except Exception:
            shutil.rmtree(destination, ignore_errors=True)
            raise
        if meta["blank_key_count"] or meta["duplicate_key_count"]:
            shutil.rmtree(destination, ignore_errors=True)
            raise ValueError(
                f"Cannot create ground truth from {name}: "
                f"{meta['blank_key_count']} blank and {meta['duplicate_key_count']} duplicate keys"
            )
        outputs[name] = {
            "source_relative_path": str(rel),
            "expected_file": f"expected/{target_name}",
            "sha256": _sha256(target),
            "category": spec["category"],
            "key": spec["key"],
            "ignore_columns": spec["ignore_columns"],
            "rows": meta["rows"],
            "columns": meta["columns"],
            "duplicate_key_count": meta["duplicate_key_count"],
            "blank_key_count": meta["blank_key_count"],
            "comparison": "keyed_cellwise",
        }

    missing_outputs = sorted(set(DEFAULT_OUTPUTS) - set(outputs))
    missing_required = sorted(REQUIRED_OUTPUTS - set(outputs))
    if missing_required:
        shutil.rmtree(destination, ignore_errors=True)
        raise FileNotFoundError(
            "Run directory lacks required assessment outputs: " + ", ".join(missing_required)
        )

    try:
        perf_path = _performance_manifest(root)
        perf = _read_json(perf_path) if perf_path else {}
        input_records = _manifest_inputs(perf)
        provenance_evidence = None
        if perf_path is not None:
            evidence_dir = destination / "evidence"
            evidence_dir.mkdir(parents=True, exist_ok=True)
            evidence_path = evidence_dir / "performance_run_manifest.json"
            shutil.copy2(perf_path, evidence_path)
            provenance_evidence = {
                "file": "evidence/performance_run_manifest.json",
                "sha256": _sha256(evidence_path),
                "source_relative_path": str(perf_path.relative_to(root)),
            }
        provenance_status, provenance_messages = _verify_input_provenance(input_records, perf)
        db_state = _database_state(root, database=database)
    except Exception:
        shutil.rmtree(destination, ignore_errors=True)
        raise
    manifest = {
        "reference_schema_version": REFERENCE_SCHEMA_VERSION,
        "reference_id": destination.name,
        "label": label or destination.name,
        "status": "PROVISIONAL",
        "created_at": _now(),
        "created_from_run": str(root),
        "source_workflow": perf.get("workflow"),
        "branchmanager_version": perf.get("branchmanager_version"),
        "workflow_version": perf.get("workflow_version"),
        "performance_run_manifest": str(perf_path.relative_to(root)) if perf_path else None,
        "input_fingerprints": input_records,
        "input_provenance_status": provenance_status,
        "input_provenance_messages": provenance_messages,
        "provenance_evidence": provenance_evidence,
        "required_outputs": sorted(REQUIRED_OUTPUTS),
        "missing_required_outputs": missing_required,
        "missing_optional_outputs": sorted(set(missing_outputs) - REQUIRED_OUTPUTS),
        "environment": _manifest_environment(perf),
        "database_state": db_state,
        "numeric_tolerance": float(tolerance),
        "outputs": outputs,
        "approval_history": [],
        "important": (
            "This snapshot is PROVISIONAL. It is not an authoritative ground truth "
            "until explicitly reviewed and approved with `branchmanager ground-truth approve`."
        ),
    }
    _json_write(destination / "ground_truth_manifest.json", manifest)
    _write_reference_readme(destination, manifest)
    return destination


def _write_reference_readme(destination: Path, manifest: dict[str, Any]) -> None:
    lines = [
        f"# BranchManager ground truth: {manifest['label']}",
        "",
        f"- Status: **{manifest['status']}**",
        f"- Created: {manifest['created_at']}",
        f"- Source run: `{manifest['created_from_run']}`",
        f"- BranchManager version: `{manifest.get('branchmanager_version') or 'unknown'}`",
        f"- Workflow version: `{manifest.get('workflow_version') or 'unknown'}`",
        "",
        "This is a provisional snapshot. Review the source run, input fingerprints,",
        "external reference versions and expected outputs before approving it.",
        "Expected files are checksummed; do not edit them manually.",
        "",
        "Captured outputs:",
    ]
    for name, item in manifest["outputs"].items():
        lines.append(f"- `{name}` — {item['rows']} rows; {len(item['columns'])} columns")
    (destination / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _approval_payload(manifest: dict[str, Any], root: Path) -> dict[str, Any]:
    return {
        "reference_id": manifest.get("reference_id"),
        "manifest_sha256": _sha256(root / "ground_truth_manifest.json"),
        "expected_outputs": {
            name: _sha256(root / item["expected_file"])
            for name, item in sorted(manifest.get("outputs", {}).items())
        },
    }


def _verify_approval(root: Path, manifest: dict[str, Any]) -> tuple[bool, str]:
    approval_path = root / APPROVAL_RECORD
    if not approval_path.is_file():
        return False, "Approved reference has no detached approval record"
    try:
        approval = _read_json(approval_path)
        expected = _approval_payload(manifest, root)
    except Exception as exc:
        return False, f"Could not verify approval record: {exc}"
    if approval != expected:
        return False, "Approved manifest or expected outputs differ from the detached approval digest"
    return True, "Approved manifest and expected outputs match the detached approval digest"


def approve_reference(reference_dir: str | Path, *, approved_by: str, reason: str) -> Path:
    """Approve a provisional reference and create a detached tamper-evident digest."""
    root = Path(reference_dir).expanduser().resolve()
    manifest_path = root / "ground_truth_manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Ground truth manifest not found: {manifest_path}")
    manifest = _read_json(manifest_path)
    if manifest.get("status") != "PROVISIONAL":
        raise ValueError(f"Only PROVISIONAL references can be approved (status={manifest.get('status')!r})")
    if not approved_by.strip() or not reason.strip():
        raise ValueError("Both --approved-by and --reason are required")
    if manifest.get("missing_required_outputs"):
        raise ValueError("Cannot approve reference with missing required outputs: " + ", ".join(manifest["missing_required_outputs"]))
    provenance_status = manifest.get("input_provenance_status", "UNAVAILABLE")
    if provenance_status not in {"VERIFIED", "VERIFIED_FROM_RUN_MANIFEST"}:
        raise ValueError(
            "Cannot approve reference without verified input provenance; "
            f"status={provenance_status}. "
            "Resolve missing or inconsistent required input provenance first."
        )
    evidence = manifest.get("provenance_evidence")
    if provenance_status == "VERIFIED_FROM_RUN_MANIFEST" and not evidence:
        raise ValueError("Cannot approve historically verified provenance without captured workflow-manifest evidence")
    if evidence:
        evidence_path = root / evidence.get("file", "")
        if (not evidence_path.is_file()
                or _sha256(evidence_path) != evidence.get("sha256")):
            raise ValueError("Cannot approve reference: supporting workflow-manifest evidence is missing or has changed")
    if not manifest.get("outputs"):
        raise ValueError("Cannot approve a reference with no captured outputs")
    for name, item in manifest.get("outputs", {}).items():
        expected = root / item["expected_file"]
        if not expected.is_file() or _sha256(expected) != item["sha256"]:
            raise ValueError(f"Reference output checksum mismatch: {name}")
        meta = _table_metadata(expected, item["key"])
        if meta["blank_key_count"] or meta["duplicate_key_count"]:
            raise ValueError(f"Cannot approve {name}: blank or duplicate row keys detected")
    manifest["status"] = "APPROVED"
    manifest["approved_at"] = _now()
    manifest["approval_history"].append({
        "approved_at": manifest["approved_at"],
        "approved_by": approved_by.strip(),
        "reason": reason.strip(),
    })
    _json_write(manifest_path, manifest)
    _write_reference_readme(root, manifest)
    _json_write(root / APPROVAL_RECORD, _approval_payload(manifest, root))
    return root


def _check_row_invariants(name: str, headers: list[str], rows: list[dict[str, str]],
                          report: list[dict[str, Any]]) -> None:
    key = DEFAULT_OUTPUTS[name]["key"]
    keys = [tuple(_normalise(row.get(k, "")) for k in key) for row in rows]
    counts = Counter(keys)
    duplicates = [k for k, count in counts.items() if count > 1]
    blank = [k for k in keys if any(not part for part in k)]
    invalid = bool(duplicates or blank)
    report.append({
        "category": "invariants",
        "check": f"{name}:unique_keys",
        "severity": "CRITICAL" if invalid else "PASS",
        "status": "FAIL" if invalid else "PASS",
        "reference": "non-empty unique row keys",
        "observed": f"{len(duplicates)} duplicate keys; {len(blank)} blank keys",
        "details": "; ".join(
            ["duplicate=" + "|".join(k) for k in duplicates[:5]]
            + ["blank=" + "|".join(k) for k in blank[:5]]
        ),
    })


def _compare_table(name: str, current: Path, reference: Path,
                   spec: dict[str, Any], tolerance: float,
                   report: list[dict[str, Any]], mode: str = "STRICT_REGRESSION",
                   inputs_changed: bool = False) -> None:
    try:
        ref_headers, ref_rows = _load_tsv(reference)
        cur_headers, cur_rows = _load_tsv(current)
    except Exception as exc:
        report.append({
            "category": spec["category"], "check": name, "severity": "CRITICAL",
            "status": "FAIL", "reference": str(reference), "observed": str(current),
            "details": f"Could not read comparison table: {exc}",
        })
        return

    key = spec["key"]
    _check_row_invariants(name, cur_headers, cur_rows, report)
    _check_row_invariants(name, ref_headers, ref_rows, report)
    ignored = set(spec.get("ignore_columns", []))
    ref_cols = [c for c in ref_headers if c not in ignored]
    cur_cols = [c for c in cur_headers if c not in ignored]
    if set(ref_cols) != set(cur_cols):
        report.append({
            "category": spec["category"], "check": f"{name}:columns",
            "severity": "HIGH", "status": "FAIL",
            "reference": str(len(ref_cols)) + " columns",
            "observed": str(len(cur_cols)) + " columns",
            "details": f"Missing={sorted(set(ref_cols)-set(cur_cols))}; added={sorted(set(cur_cols)-set(ref_cols))}",
        })
    else:
        report.append({
            "category": spec["category"], "check": f"{name}:columns",
            "severity": "PASS", "status": "PASS",
            "reference": str(len(ref_cols)), "observed": str(len(cur_cols)),
            "details": "Comparable column sets match",
        })

    def keyed(rows):
        result = {}
        for row in rows:
            k = tuple(_normalise(row.get(field, "")) for field in key)
            result[k] = row
        return result

    ref_map, cur_map = keyed(ref_rows), keyed(cur_rows)
    missing = sorted(set(ref_map) - set(cur_map))
    added = sorted(set(cur_map) - set(ref_map))
    if mode == "PRODUCTION_COMPATIBILITY" and inputs_changed:
        row_severity = "INFO" if (missing or added) else "PASS"
    else:
        row_severity = "HIGH" if (missing or added) else "PASS"
    report.append({
        "category": spec["category"], "check": f"{name}:row_keys",
        "severity": row_severity, "status": "INFO" if row_severity == "INFO" else ("FAIL" if row_severity == "HIGH" else "PASS"),
        "reference": len(ref_map), "observed": len(cur_map),
        "details": f"Missing={len(missing)} {missing[:8]}; added={len(added)} {added[:8]}",
    })

    changed = []
    comparable_cols = [c for c in ref_cols if c in cur_cols]
    for k in sorted(set(ref_map) & set(cur_map)):
        old, new = ref_map[k], cur_map[k]
        for col in comparable_cols:
            if not _equivalent(old.get(col, ""), new.get(col, ""), tolerance):
                changed.append((k, col, _normalise(old.get(col, "")), _normalise(new.get(col, ""))))
    report.append({
        "category": spec["category"], "check": f"{name}:values",
        "severity": "HIGH" if changed else "PASS",
        "status": "FAIL" if changed else "PASS",
        "reference": "unchanged values",
        "observed": f"{len(changed)} changed cells",
        "details": "; ".join(
            f"{'|'.join(k)}:{col} [{old!r} -> {new!r}]"
            for k, col, old, new in changed[:30]
        ),
    })


def _validate_database(current_root: Path, reference_state: dict[str, Any] | None,
                       report: list[dict[str, Any]], inputs_changed: bool,
                       database: str | Path | None = None, mode: str = "STRICT_REGRESSION") -> None:
    try:
        state = _database_state(current_root, database=database)
    except Exception as exc:
        report.append({
            "category": "database", "check": "database_readability",
            "severity": "CRITICAL", "status": "FAIL", "reference": "readable SQLite database",
            "observed": "error", "details": str(exc),
        })
        return
    if state is None:
        report.append({
            "category": "database", "check": "database_presence",
            "severity": "HIGH", "status": "FAIL", "reference": "project SQLite database",
            "observed": "not found", "details": "No database found; provide --db if stored outside the run directory",
        })
        return
    integrity_ok = state["integrity_check"] == "ok" and state["foreign_key_violations"] == 0
    report.append({
        "category": "database", "check": "sqlite_integrity",
        "severity": "PASS" if integrity_ok else "CRITICAL", "status": "PASS" if integrity_ok else "FAIL",
        "reference": "integrity_check=ok; foreign_key_violations=0",
        "observed": f"{state['integrity_check']}; foreign_key_violations={state['foreign_key_violations']}",
        "details": state["path"],
    })
    if not reference_state:
        report.append({
            "category": "database", "check": "database_reference_state",
            "severity": "LOW", "status": "WARN", "reference": "database snapshot",
            "observed": "not captured", "details": "Reference was created without a readable project database",
        })
        return
    if state.get("schema_version") != reference_state.get("schema_version"):
        report.append({
            "category": "database", "check": "schema_version", "severity": "HIGH", "status": "FAIL",
            "reference": reference_state.get("schema_version"), "observed": state.get("schema_version"),
            "details": "Database schema version differs from the reference",
        })
    else:
        report.append({
            "category": "database", "check": "schema_version", "severity": "PASS", "status": "PASS",
            "reference": reference_state.get("schema_version"), "observed": state.get("schema_version"),
            "details": "Database schema version matches",
        })
    reference_tables = reference_state.get("table_fingerprints", {})
    current_tables = state.get("table_fingerprints", {})
    if not reference_tables:
        # Backward compatibility for early schema-v1 snapshots that stored counts only.
        for table, expected_count in reference_state.get("table_counts", {}).items():
            observed_count = state.get("table_counts", {}).get(table)
            if observed_count != expected_count:
                severity = "INFO" if mode == "PRODUCTION_COMPATIBILITY" and inputs_changed else "HIGH"
                report.append({
                    "category": "database", "check": f"table_count:{table}",
                    "severity": severity, "status": "INFO" if severity == "INFO" else "FAIL",
                    "reference": expected_count, "observed": observed_count,
                    "details": "Legacy reference contains counts only; row-level comparison unavailable",
                })
        return
    for table in sorted(set(reference_tables) | set(current_tables)):
        old = reference_tables.get(table)
        new = current_tables.get(table)
        if old is None or new is None:
            severity = "INFO" if mode == "PRODUCTION_COMPATIBILITY" and inputs_changed else "HIGH"
            report.append({
                "category": "database", "check": f"table_presence:{table}",
                "severity": severity, "status": "INFO" if severity == "INFO" else "FAIL",
                "reference": "present" if old else "absent", "observed": "present" if new else "absent",
                "details": "Database table presence changed",
            })
            continue
        if old.get("columns") != new.get("columns"):
            report.append({
                "category": "database", "check": f"table_columns:{table}", "severity": "HIGH", "status": "FAIL",
                "reference": old.get("columns"), "observed": new.get("columns"),
                "details": "Database table columns changed",
            })
            continue
        old_records, new_records = old.get("records"), new.get("records")
        if old_records is None or new_records is None:
            equal = old.get("unordered_row_hashes") == new.get("unordered_row_hashes")
            if not equal:
                severity = "INFO" if mode == "PRODUCTION_COMPATIBILITY" and inputs_changed else "HIGH"
                report.append({
                    "category": "database", "check": f"table_rows:{table}", "severity": severity,
                    "status": "INFO" if severity == "INFO" else "FAIL",
                    "reference": old.get("rows"), "observed": new.get("rows"),
                    "details": "Unkeyed table row content differs",
                })
            continue
        old_keys, new_keys = set(old_records), set(new_records)
        added, removed = sorted(new_keys-old_keys), sorted(old_keys-new_keys)
        changed = sorted(k for k in old_keys & new_keys if old_records[k] != new_records[k])
        if added or removed or changed:
            severity = "INFO" if mode == "PRODUCTION_COMPATIBILITY" and inputs_changed and not changed else "HIGH"
            report.append({
                "category": "database", "check": f"table_records:{table}",
                "severity": severity, "status": "INFO" if severity == "INFO" else "FAIL",
                "reference": old.get("rows"), "observed": new.get("rows"),
                "details": f"added={len(added)} {added[:5]}; removed={len(removed)} {removed[:5]}; changed={len(changed)} {changed[:5]}",
            })
        else:
            report.append({
                "category": "database", "check": f"table_records:{table}", "severity": "PASS", "status": "PASS",
                "reference": old.get("rows"), "observed": new.get("rows"),
                "details": "Database records match by primary key",
            })


def _write_reports(outdir: Path, result: dict[str, Any]) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    _json_write(outdir / "validation_summary.json", result)
    fields = ["category", "check", "severity", "status", "reference", "observed", "details"]
    with (outdir / "drift_report.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(result["checks"])
    counts = Counter(item["severity"] for item in result["checks"])
    rows = "\n".join(
        "<tr>" + "".join(f"<td>{html.escape(str(item.get(k, '')))}</td>"
                          for k in fields) + "</tr>"
        for item in result["checks"]
    )
    page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>BranchManager validation</title>
<style>body{{font:14px system-ui,sans-serif;margin:2rem;color:#20242a}}table{{border-collapse:collapse;width:100%;font-size:12px}}th,td{{border:1px solid #ddd;padding:6px;text-align:left;vertical-align:top}}th{{background:#f2f4f7}}.status{{font-weight:700}}</style></head>
<body><h1>BranchManager validation report</h1>
<p class="status">Status: {html.escape(result['status'])}</p>
<p>Reference: {html.escape(result['reference_id'])} ({html.escape(result['reference_status'])})</p>
<p>Run: {html.escape(result['run_dir'])}</p>
<p>Checks: {len(result['checks'])}; severity counts: {html.escape(json.dumps(dict(counts), sort_keys=True))}</p>
<table><thead><tr>{''.join(f'<th>{html.escape(k)}</th>' for k in fields)}</tr></thead><tbody>{rows}</tbody></table>
</body></html>"""
    (outdir / "validation_report.html").write_text(page, encoding="utf-8")


def validate_run(run_dir: str | Path, reference_dir: str | Path, *,
                 output_dir: str | Path | None = None,
                 database: str | Path | None = None,
                 mode: str = "STRICT_REGRESSION") -> dict[str, Any]:
    """Compare a completed run with a reference; never mutate the run or DB."""
    root = Path(run_dir).expanduser().resolve()
    reference = Path(reference_dir).expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"Run directory does not exist: {root}")
    manifest_path = reference / "ground_truth_manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Ground truth manifest not found: {manifest_path}")
    manifest = _read_json(manifest_path)
    if manifest.get("reference_schema_version") != REFERENCE_SCHEMA_VERSION:
        raise ValueError("Unsupported ground truth manifest schema version")

    if mode not in {"STRICT_REGRESSION", "PRODUCTION_COMPATIBILITY"}:
        raise ValueError("mode must be STRICT_REGRESSION or PRODUCTION_COMPATIBILITY")
    checks: list[dict[str, Any]] = []
    perf_path = _performance_manifest(root)
    perf = _read_json(perf_path) if perf_path else {}
    current_inputs = _manifest_inputs(perf)
    reference_inputs_raw = manifest.get("input_fingerprints", {})
    reference_inputs = {
        key: (value if isinstance(value, dict) else {"sha256": str(value), "path": "", "required": True})
        for key, value in reference_inputs_raw.items()
    }
    current_provenance, provenance_messages = _verify_input_provenance(current_inputs, perf)
    reference_provenance = manifest.get("input_provenance_status", "UNAVAILABLE")
    inputs_changed = bool(reference_inputs and current_inputs and (
        set(reference_inputs) != set(current_inputs) or {
            key: value.get("sha256") for key, value in reference_inputs.items()
        } != {
            key: value.get("sha256") for key, value in current_inputs.items()
        }
    ))
    if current_provenance == "VERIFIED" and inputs_changed:
        current_provenance = "CHANGED"
    provenance_severity = "PASS"
    provenance_result = "PASS"
    accepted_provenance = {"VERIFIED", "VERIFIED_FROM_RUN_MANIFEST"}
    if current_provenance == "INVALID" or reference_provenance == "INVALID":
        provenance_severity, provenance_result = "CRITICAL", "FAIL"
    elif current_provenance == "UNAVAILABLE" or reference_provenance not in accepted_provenance:
        provenance_severity, provenance_result = "HIGH", "FAIL"
    elif inputs_changed:
        provenance_severity, provenance_result = ("INFO", "INFO") if mode == "PRODUCTION_COMPATIBILITY" else ("HIGH", "FAIL")
    checks.append({
        "category": "provenance", "check": "input_fingerprints",
        "severity": provenance_severity, "status": provenance_result,
        "reference": reference_provenance, "observed": current_provenance,
        "details": "; ".join(provenance_messages) if provenance_messages else (
            "Input files were independently re-hashed and recorded fingerprints verified"
            if current_provenance == "VERIFIED"
            else "Inputs verified; removed temporary database state is supported by the captured workflow manifest"
            if current_provenance == "VERIFIED_FROM_RUN_MANIFEST"
            else "Input fingerprints differ from the reference" if inputs_changed else "Input provenance verified"
        ),
    })
    reference_evidence = manifest.get("provenance_evidence")
    if reference_evidence:
        evidence_path = reference / reference_evidence.get("file", "")
        evidence_ok = (
            evidence_path.is_file()
            and _sha256(evidence_path) == reference_evidence.get("sha256")
        )
        checks.append({
            "category": "provenance",
            "check": "supporting_workflow_manifest_integrity",
            "severity": "PASS" if evidence_ok else "CRITICAL",
            "status": "PASS" if evidence_ok else "FAIL",
            "reference": reference_evidence.get("sha256"),
            "observed": _sha256(evidence_path) if evidence_path.is_file() else None,
            "details": (
                "Captured supporting workflow manifest matches its recorded SHA256"
                if evidence_ok else
                "Captured supporting workflow manifest is missing or has changed"
            ),
        })
    elif manifest.get("input_provenance_status") == "VERIFIED_FROM_RUN_MANIFEST":
        checks.append({
            "category": "provenance",
            "check": "supporting_workflow_manifest_integrity",
            "severity": "CRITICAL",
            "status": "FAIL",
            "reference": None,
            "observed": None,
            "details": "Historical provenance has no captured supporting workflow manifest",
        })

    current_environment = _manifest_environment(perf)
    ref_environment = manifest.get("environment", {})
    env_changes = []
    for field in ("branchmanager_version", "workflow_version", "external_tools", "policies", "reference_identity"):
        if ref_environment.get(field) != current_environment.get(field):
            env_changes.append(field)
    checks.append({
        "category": "provenance", "check": "software_and_reference_environment",
        "severity": "INFO" if env_changes else "PASS",
        "status": "INFO" if env_changes else "PASS",
        "reference": "recorded environment", "observed": "recorded environment",
        "details": "Changed: " + ", ".join(env_changes) if env_changes else "Recorded versions and policies match",
    })

    for name, spec in manifest.get("outputs", {}).items():
        expected = reference / spec["expected_file"]
        if not expected.is_file() or _sha256(expected) != spec.get("sha256"):
            checks.append({
                "category": "ground_truth", "check": f"{name}:reference_integrity",
                "severity": "CRITICAL", "status": "FAIL",
                "reference": spec.get("sha256"), "observed": "missing or checksum mismatch",
                "details": "Reference output is missing or has been modified",
            })
            continue
        current = root / spec.get("source_relative_path", "")
        if not current.is_file():
            try:
                current = _find_file(root, name)
            except Exception as exc:
                checks.append({
                    "category": spec["category"], "check": f"{name}:resolution",
                    "severity": "CRITICAL", "status": "FAIL", "reference": "one unambiguous output",
                    "observed": "ambiguous", "details": str(exc),
                })
                continue
        if current is None or not current.is_file():
            checks.append({
                "category": spec["category"], "check": f"{name}:presence",
                "severity": "CRITICAL", "status": "FAIL",
                "reference": "present", "observed": "missing",
                "details": "Expected output was not found in current run",
            })
            continue
        _compare_table(name, current, expected, spec,
                       float(manifest.get("numeric_tolerance", 1e-6)), checks,
                       mode=mode, inputs_changed=inputs_changed)

    for name in manifest.get("missing_optional_outputs", []):
        try:
            newly_available = _find_file(root, name)
        except Exception as exc:
            checks.append({
                "category": "outputs", "check": f"{name}:optional_resolution",
                "severity": "HIGH", "status": "FAIL", "reference": "not captured in reference",
                "observed": "ambiguous", "details": str(exc),
            })
            continue
        checks.append({
            "category": "outputs", "check": f"{name}:optional_coverage",
            "severity": "INFO" if newly_available else "PASS",
            "status": "INFO" if newly_available else "PASS",
            "reference": "not captured in reference",
            "observed": "now present" if newly_available else "still absent",
            "details": "This output is not compared because it was absent from the reference; create a new reference to include it" if newly_available else "Optional output remains absent",
        })

    # Cross-output and biological/operational invariants.
    try:
        current_assessment = _find_file(root, "sequence_assessment.tsv")
        current_summary = _find_file(root, "selection_summary.tsv")
    except Exception as exc:
        current_assessment = current_summary = None
        checks.append({
            "category": "invariants", "check": "assessment_output_resolution",
            "severity": "CRITICAL", "status": "FAIL", "reference": "one unambiguous assessment/summary pair",
            "observed": "ambiguous", "details": str(exc),
        })
    if current_assessment and current_summary:
        try:
            ah, ar = _load_tsv(current_assessment)
            sh, sr = _load_tsv(current_summary)
            aid = {row.get("ID", "") for row in ar}
            sid = {row.get("SequenceID", "") for row in sr}
            ok = sid.issubset(aid) and "" not in sid
            checks.append({
                "category": "invariants", "check": "SAB_shortlist_is_assessment_subset",
                "severity": "PASS" if ok else "CRITICAL",
                "status": "PASS" if ok else "FAIL",
                "reference": "all SAB SequenceID values occur in assessment ID",
                "observed": f"{len(sid)} shortlist IDs; {len(sid - aid)} absent from assessment",
                "details": "SAB summary is a subset of the full assessment" if ok else str(sorted(sid - aid)[:20]),
            })
        except Exception as exc:
            checks.append({
                "category": "invariants", "check": "SAB_shortlist_is_assessment_subset",
                "severity": "CRITICAL", "status": "FAIL",
                "reference": "readable assessment and summary", "observed": "error",
                "details": str(exc),
            })

    try:
        audit = _find_file(root, "cohort_review_application_audit.tsv")
    except Exception as exc:
        audit = None
        checks.append({
            "category": "invariants", "check": "marker_review_audit_resolution",
            "severity": "CRITICAL", "status": "FAIL", "reference": "one unambiguous review audit",
            "observed": "ambiguous", "details": str(exc),
        })
    if audit:
        try:
            headers, rows = _load_tsv(audit)
            if {"decision", "prepared_qc_class", "application_status", "included_in_joint_assessment"}.issubset(headers):
                bad_override = [
                    r for r in rows if r["decision"].upper() == "APPROVE_FAIL_QC"
                    and (r["prepared_qc_class"] != "FAIL_QC"
                         or r["application_status"] != "INCLUDED_BY_MANUAL_FAILED_QC_OVERRIDE")
                ]
                bad_reseq = [
                    r for r in rows if r["decision"].upper() == "RESEQUENCE"
                    and r["included_in_joint_assessment"].upper() == "YES"
                ]
                ok = not bad_override and not bad_reseq
                checks.append({
                    "category": "invariants", "check": "manual_marker_review_semantics",
                    "severity": "PASS" if ok else "CRITICAL",
                    "status": "PASS" if ok else "FAIL",
                    "reference": "APPROVE_FAIL_QC retained as failed QC override; RESEQUENCE excluded",
                    "observed": f"invalid overrides={len(bad_override)}; resequence included={len(bad_reseq)}",
                    "details": "Review audit semantics are consistent" if ok else "Inspect cohort_review_application_audit.tsv",
                })
        except Exception as exc:
            checks.append({
                "category": "invariants", "check": "manual_marker_review_semantics",
                "severity": "CRITICAL", "status": "FAIL",
                "reference": "readable review audit", "observed": "error",
                "details": str(exc),
            })

    approval_ok, approval_detail = _verify_approval(reference, manifest) if manifest.get("status") == "APPROVED" else (True, "Reference is provisional")
    checks.append({
        "category": "ground_truth", "check": "approval_integrity",
        "severity": "PASS" if approval_ok else "CRITICAL",
        "status": "PASS" if approval_ok else "FAIL",
        "reference": "valid detached approval digest", "observed": approval_detail,
        "details": approval_detail,
    })
    _validate_database(root, manifest.get("database_state"), checks, inputs_changed,
                       database=database, mode=mode)

    severities = Counter(item["severity"] for item in checks)
    critical = severities.get("CRITICAL", 0)
    high = severities.get("HIGH", 0)
    reference_status = str(manifest.get("status", "UNKNOWN")).upper()
    if critical:
        status = "CRITICAL_FAILURE"
    elif reference_status != "APPROVED":
        status = "REFERENCE_PROVISIONAL"
    elif provenance_severity == "HIGH":
        status = "VALIDATION_INCOMPLETE"
    elif high:
        status = "DRIFT_DETECTED"
    else:
        status = "PASS"
    result = {
        "validation_schema_version": VALIDATION_SCHEMA_VERSION,
        "validated_at": _now(),
        "status": status,
        "run_dir": str(root),
        "reference_id": manifest.get("reference_id"),
        "reference_status": reference_status,
        "reference_created_at": manifest.get("created_at"),
        "branchmanager_version": perf.get("branchmanager_version"),
        "workflow_version": perf.get("workflow_version"),
        "validation_mode": mode,
        "provenance_status": current_provenance,
        "database_path": str(database) if database else None,
        "checks_total": len(checks),
        "severity_counts": dict(severities),
        "checks": checks,
        "read_only": True,
        "note": "Validation is observational: it does not alter run outputs or the project database.",
    }
    out = Path(output_dir).expanduser().resolve() if output_dir else root / "validation"
    _write_reports(out, result)
    result["report_directory"] = str(out)
    return result


def add_cli_arguments(subparsers) -> None:
    """Register additive validation commands on the existing parser."""
    gt = subparsers.add_parser(
        "ground-truth",
        help="Create and approve versioned reference snapshots for regression validation.",
        description="Manage immutable-by-default BranchManager ground-truth snapshots.",
    )
    gt_sub = gt.add_subparsers(dest="ground_truth_action", required=True)
    init = gt_sub.add_parser("init", help="Create a provisional reference snapshot from a completed run.")
    init.add_argument("--run-dir", required=True, help="Completed BranchManager run directory.")
    init.add_argument("--out", required=True, help="New reference directory; must not already exist.")
    init.add_argument("--label", default="", help="Human-readable reference label.")
    init.add_argument("--numeric-tolerance", type=float, default=1e-6)
    init.add_argument("--db", default=None, help="Project SQLite database; required if not unambiguously inside the run directory.")
    approve = gt_sub.add_parser("approve", help="Approve a reviewed provisional ground-truth snapshot.")
    approve.add_argument("--reference", required=True, help="Reference directory.")
    approve.add_argument("--approved-by", required=True, help="Name or role of the reviewer.")
    approve.add_argument("--reason", required=True, help="Reason for approving this reference.")
    validate = subparsers.add_parser(
        "validate",
        help="Compare a completed run against ground truth and check biological/workflow invariants.",
        description="Read-only validation; does not modify workflow outputs or the project database.",
    )
    validate.add_argument("--run-dir", required=True, help="Completed BranchManager run directory.")
    validate.add_argument("--reference", required=True, help="Ground-truth reference directory.")
    validate.add_argument("--out", default=None, help="Validation report directory (default: RUN_DIR/validation).")
    validate.add_argument("--db", default=None, help="Project SQLite database to inspect (read-only).")
    validate.add_argument("--mode", choices=("STRICT_REGRESSION", "PRODUCTION_COMPATIBILITY"),
                          default="STRICT_REGRESSION", help="Fixed-cohort regression or changing production cohort comparison.")


def run_cli(args) -> int:
    if args.command == "ground-truth":
        if args.ground_truth_action == "init":
            path = create_reference(args.run_dir, args.out, label=args.label,
                                    tolerance=args.numeric_tolerance, database=args.db)
            print(f"Created PROVISIONAL ground truth: {path}")
            print("Review the reference outputs and provenance before approving it.")
            return 0
        if args.ground_truth_action == "approve":
            path = approve_reference(args.reference, approved_by=args.approved_by, reason=args.reason)
            print(f"Approved ground truth: {path}")
            return 0
    if args.command == "validate":
        result = validate_run(args.run_dir, args.reference, output_dir=args.out,
                              database=args.db, mode=args.mode)
        print(f"Validation status: {result['status']}")
        print(f"Checks: {result['checks_total']}; severity counts: {result['severity_counts']}")
        print(f"Report: {result['report_directory']}")
        if result["status"] == "CRITICAL_FAILURE":
            return 3
        if result["status"] == "DRIFT_DETECTED":
            return 1
        if result["status"] == "REFERENCE_PROVISIONAL":
            return 2
        if result["status"] == "VALIDATION_INCOMPLETE":
            return 4
        return 0
    return 0
