#!/usr/bin/env python3
"""Audit saved APPROVE_FAIL_QC decisions against the prepared recovery FASTA."""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sequence_id", "decision", "reviewer", "notes"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Review TSV is missing required columns: {', '.join(sorted(missing))}")
        return [dict(row) for row in reader]


def read_fasta_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    current_id = None
    sequence_length = 0
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None and sequence_length > 0:
                    ids.add(current_id)
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header in {path}")
                current_id = header.split()[0].split("|", 1)[0]
                sequence_length = 0
            else:
                if current_id is None:
                    raise ValueError(f"Sequence data appears before a FASTA header in {path}")
                sequence_length += len("".join(line.split()))
    if current_id is not None and sequence_length > 0:
        ids.add(current_id)
    return ids


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--review", required=True, type=Path, help="Saved Sanger Studio cohort review decisions TSV")
    parser.add_argument("--recoverable", required=True, type=Path, help="PREPARE 00_cohort_preparation/cohort_recoverable_failed_markers.fasta")
    args = parser.parse_args()

    rows = read_tsv(args.review)
    ids = [str(row.get("sequence_id") or "").strip() for row in rows]
    if len(ids) != len(set(ids)):
        duplicates = sorted(sid for sid, count in Counter(ids).items() if sid and count > 1)
        raise ValueError("Duplicate sequence IDs in review TSV: " + ", ".join(duplicates[:20]))
    decisions = Counter(str(row.get("decision") or "").strip().upper() for row in rows)
    approved = [row for row in rows if str(row.get("decision") or "").strip().upper() == "APPROVE_FAIL_QC"]
    incomplete = [
        str(row.get("sequence_id") or "") for row in approved
        if not str(row.get("reviewer") or "").strip() or not str(row.get("notes") or "").strip()
    ]
    prepared_ids = read_fasta_ids(args.recoverable)
    missing = sorted(
        str(row.get("sequence_id") or "").strip() for row in approved
        if str(row.get("sequence_id") or "").strip() not in prepared_ids
    )

    print(f"Review rows: {len(rows)}")
    print("Decisions: " + ", ".join(f"{key or '(blank)'}={value}" for key, value in sorted(decisions.items())))
    print(f"APPROVE_FAIL_QC decisions: {len(approved)}")
    print(f"Prepared recoverable FASTA records: {len(prepared_ids)}")
    print(f"Approved failed-QC IDs present: {len(approved) - len(missing)}")
    if incomplete:
        print("ERROR: approvals without reviewer or notes: " + ", ".join(incomplete))
    if missing:
        print(f"ERROR: {len(missing)} approved failed-QC marker(s) are absent from the prepared recoverable FASTA:")
        for sid in missing:
            print(f"  {sid}")
        print("Do not run cohort RESUME. Regenerate PREPARE with the current review policy, then rerun this audit.")
        return 2
    if incomplete:
        return 2
    print("PASS: every APPROVE_FAIL_QC decision has a prepared non-empty sequence and reviewer notes.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
