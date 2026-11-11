import csv
import json
import sqlite3
from pathlib import Path

import pytest

from branchmanager.validation import (
    approve_reference,
    create_reference,
    validate_run,
    _resolve_database,
)


def _write_tsv(path: Path, headers, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _make_run(root: Path, *, recommendation="PRIMARY", include_summary_id=True):
    assessment = root / "04_assistant_cohort/01_joint_assistant/03_performance_review_hiring_panel/assessment"
    marker_input = root / "inputs" / "markers.fasta"
    marker_input.parent.mkdir(parents=True, exist_ok=True)
    marker_input.write_text(">QUB_0001\nACGT\n", encoding="utf-8")
    import hashlib
    marker_digest = hashlib.sha256(marker_input.read_bytes()).hexdigest()
    inputs = [{"role": "marker_fasta", "path": str(marker_input.resolve()),
               "sha256": marker_digest, "required": True, "exists": True, "kind": "file"}]
    manifest = {
        "workflow": "performance_review",
        "branchmanager_version": "1.2.0",
        "workflow_version": "1.2.0",
        "inputs": inputs,
        "environment": {"python": "3.12"},
        "external_tools": {"vsearch": {"version": "2.29"}},
        "policies": {"qc": "2.6"},
        "reference_identity": {"GTDB": "r232"},
    }
    assessment.mkdir(parents=True, exist_ok=True)
    (assessment.parent / "run_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    _write_tsv(assessment / "sequence_assessment.tsv",
               ["ID", "MarkerQCClass", "Recommendation"],
               [{"ID": "QUB_0001", "MarkerQCClass": "PASS_HIGH_CONFIDENCE",
                 "Recommendation": recommendation},
                {"ID": "QUB_0002", "MarkerQCClass": "PASS_WITH_WARNINGS",
                 "Recommendation": "SECONDARY"}])
    _write_tsv(assessment / "selection_summary.tsv",
               ["SequenceID", "Recommendation", "SABDecision"],
               [{"SequenceID": "QUB_0001", "Recommendation": recommendation, "SABDecision": ""}])
    _write_tsv(assessment / "sequencing_sets.tsv",
               ["CandidateID", "SetID", "SetRole"],
               [{"CandidateID": "QUB_0001", "SetID": "SET_1", "SetRole": "PRIMARY"}])
    _write_tsv(assessment / "novelty_metrics.tsv",
               ["ID", "NearestIdentity"],
               [{"ID": "QUB_0001", "NearestIdentity": "98.5"},
                {"ID": "QUB_0002", "NearestIdentity": "97.0"}])
    audit = root / "04_assistant_cohort/01_joint_inputs/cohort_review_application_audit.tsv"
    _write_tsv(audit,
               ["sequence_id", "decision", "prepared_qc_class", "application_status",
                "included_in_joint_assessment"],
               [{"sequence_id": "QUB_0001", "decision": "APPROVE",
                 "prepared_qc_class": "PASS_HIGH_CONFIDENCE",
                 "application_status": "INCLUDED_AFTER_MANUAL_APPROVAL",
                 "included_in_joint_assessment": "YES"},
                {"sequence_id": "QUB_0002", "decision": "RESEQUENCE",
                 "prepared_qc_class": "FAIL_QC",
                 "application_status": "EXCLUDED_RESEQUENCE",
                 "included_in_joint_assessment": "NO"}])
    _write_tsv(root / "04_assistant_cohort/01_joint_inputs/cohort_joint_marker_qc.tsv",
               ["SequenceID", "Status"],
               [{"SequenceID": "QUB_0001", "Status": "PASS"},
                {"SequenceID": "QUB_0002", "Status": "FAIL"}])
    db = sqlite3.connect(root / "project.sqlite")
    db.execute("CREATE TABLE schema_migrations (version INTEGER)")
    db.execute("INSERT INTO schema_migrations VALUES (6)")
    db.execute("CREATE TABLE sequences (id TEXT PRIMARY KEY)")
    db.executemany("INSERT INTO sequences VALUES (?)", [("QUB_0001",), ("QUB_0002",)])
    db.execute("CREATE TABLE genome_records (id TEXT)")
    db.commit()
    db.close()


def test_relative_database_path_resolves_from_working_directory_before_run_directory(tmp_path, monkeypatch):
    run = tmp_path / "run"
    run.mkdir()
    external_db = tmp_path / "project.sqlite"
    sqlite3.connect(external_db).close()
    monkeypatch.chdir(tmp_path)
    assert _resolve_database(run, "project.sqlite") == external_db.resolve()


def test_create_reference_is_provisional_and_does_not_touch_run(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    _make_run(run)
    before = (run / "project.sqlite").read_bytes()
    ref = create_reference(run, tmp_path / "reference", label="test cohort")
    manifest = json.loads((ref / "ground_truth_manifest.json").read_text())
    assert manifest["status"] == "PROVISIONAL"
    assert set(manifest["outputs"]) >= {"sequence_assessment.tsv", "selection_summary.tsv"}
    assert (run / "project.sqlite").read_bytes() == before


def test_approve_reference_requires_reason_and_preserves_expected_checksums(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    _make_run(run)
    ref = create_reference(run, tmp_path / "reference")
    with pytest.raises(ValueError):
        approve_reference(ref, approved_by="reviewer", reason="")
    approve_reference(ref, approved_by="reviewer", reason="Reviewed against accepted cohort")
    manifest = json.loads((ref / "ground_truth_manifest.json").read_text())
    assert manifest["status"] == "APPROVED"
    assert manifest["approval_history"][0]["approved_by"] == "reviewer"


def test_identical_approved_run_passes_and_writes_reports(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    _make_run(run)
    ref = create_reference(run, tmp_path / "reference")
    approve_reference(ref, approved_by="reviewer", reason="Validated reference")
    result = validate_run(run, ref)
    assert result["status"] == "PASS"
    assert (run / "validation/validation_summary.json").is_file()
    assert (run / "validation/drift_report.tsv").is_file()
    assert (run / "validation/validation_report.html").is_file()


def test_changed_recommendation_is_reported_as_drift(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    _make_run(run)
    ref = create_reference(run, tmp_path / "reference")
    approve_reference(ref, approved_by="reviewer", reason="Validated reference")
    # Change the recommendation in the current run only.
    path = run / "04_assistant_cohort/01_joint_assistant/03_performance_review_hiring_panel/assessment/sequence_assessment.tsv"
    rows = path.read_text().replace("QUB_0001\tPASS_HIGH_CONFIDENCE\tPRIMARY",
                                    "QUB_0001\tPASS_HIGH_CONFIDENCE\tSECONDARY")
    path.write_text(rows)
    result = validate_run(run, ref)
    assert result["status"] == "DRIFT_DETECTED"
    assert any(x["check"].endswith(":values") and "Recommendation" in x["details"]
               for x in result["checks"])


def test_provisional_reference_never_reports_pass(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    _make_run(run)
    ref = create_reference(run, tmp_path / "reference")
    result = validate_run(run, ref)
    assert result["status"] == "REFERENCE_PROVISIONAL"


def test_reference_checksum_tampering_is_critical(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    _make_run(run)
    ref = create_reference(run, tmp_path / "reference")
    approve_reference(ref, approved_by="reviewer", reason="Validated reference")
    with (ref / "expected/selection_summary.tsv").open("a") as handle:
        handle.write("\n")
    result = validate_run(run, ref)
    assert result["status"] == "CRITICAL_FAILURE"
    assert any(x["check"] == "selection_summary.tsv:reference_integrity"
               for x in result["checks"])


def test_validation_commands_are_additive_to_cli_parser():
    from branchmanager.cli.parser import build_parser

    parser = build_parser()
    init = parser.parse_args([
        "ground-truth", "init", "--run-dir", "run", "--out", "reference"
    ])
    validate = parser.parse_args([
        "validate", "--run-dir", "run", "--reference", "reference"
    ])
    assert init.command == "ground-truth"
    assert init.ground_truth_action == "init"
    assert validate.command == "validate"


def test_reference_rejects_blank_and_duplicate_keys(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    _make_run(run)
    path = run / "04_assistant_cohort/01_joint_assistant/03_performance_review_hiring_panel/assessment/sequence_assessment.tsv"
    path.write_text("ID\tMarkerQCClass\tRecommendation\n\tPASS\tPRIMARY\n", encoding="utf-8")
    with pytest.raises(ValueError, match="blank and .* duplicate keys"):
        create_reference(run, tmp_path / "reference")


def test_missing_input_provenance_is_not_a_pass(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    _make_run(run)
    ref = create_reference(run, tmp_path / "reference")
    approve_reference(ref, approved_by="reviewer", reason="Reviewed reference")
    # Remove the declared source input after the reference was created.
    (run / "inputs/markers.fasta").unlink()
    result = validate_run(run, ref)
    assert result["status"] == "VALIDATION_INCOMPLETE"
    assert result["provenance_status"] == "UNAVAILABLE"


def test_removed_temporary_project_database_is_verified_from_captured_run_manifest(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    _make_run(run)

    import hashlib
    db_path = run / "04_assistant_cohort/00_cohort_preparation/.assistant_cohort_project.sqlite"
    db_path.parent.mkdir(parents=True, exist_ok=True)
    db_path.write_bytes((run / "project.sqlite").read_bytes())
    before_sha = hashlib.sha256(db_path.read_bytes()).hexdigest()
    manifest_path = run / "04_assistant_cohort/01_joint_assistant/03_performance_review_hiring_panel/run_manifest.json"
    run_manifest = json.loads(manifest_path.read_text())
    run_manifest["inputs"].append({
        "role": "project_database_before",
        "path": str(db_path.resolve()),
        "sha256": before_sha,
        "required": True,
        "exists": True,
        "kind": "file",
    })
    run_manifest.setdefault("outputs", []).append({
        "role": "project_database_after",
        "path": str(db_path.resolve()),
        "sha256": "a" * 64,
        "required": True,
        "exists": True,
        "kind": "file",
    })
    manifest_path.write_text(json.dumps(run_manifest), encoding="utf-8")

    # The temporary database is deliberately removed after its workflow
    # before/after states have been recorded.
    db_path.unlink()

    ref = create_reference(run, tmp_path / "reference")
    ground_truth = json.loads((ref / "ground_truth_manifest.json").read_text())
    assert ground_truth["input_provenance_status"] == "VERIFIED_FROM_RUN_MANIFEST"
    assert ground_truth["provenance_evidence"]["file"] == "evidence/performance_run_manifest.json"
    assert (ref / ground_truth["provenance_evidence"]["file"]).is_file()

    approve_reference(ref, approved_by="reviewer", reason="Reviewed historical database provenance")
    result = validate_run(run, ref)
    assert result["provenance_status"] == "VERIFIED_FROM_RUN_MANIFEST"
    assert result["status"] == "PASS"
    assert any(
        check["check"] == "supporting_workflow_manifest_integrity"
        and check["status"] == "PASS"
        for check in result["checks"]
    )


def test_historical_database_provenance_requires_before_and_after_records(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    _make_run(run)
    db_path = run / "04_assistant_cohort/00_cohort_preparation/.assistant_cohort_project.sqlite"
    db_path.parent.mkdir(parents=True, exist_ok=True)
    db_path.write_bytes((run / "project.sqlite").read_bytes())
    import hashlib
    db_sha = hashlib.sha256(db_path.read_bytes()).hexdigest()
    manifest_path = run / "04_assistant_cohort/01_joint_assistant/03_performance_review_hiring_panel/run_manifest.json"
    run_manifest = json.loads(manifest_path.read_text())
    run_manifest["inputs"].append({
        "role": "project_database_before", "path": str(db_path.resolve()),
        "sha256": db_sha, "required": True, "exists": True, "kind": "file",
    })
    manifest_path.write_text(json.dumps(run_manifest), encoding="utf-8")
    db_path.unlink()

    ref = create_reference(run, tmp_path / "reference")
    ground_truth = json.loads((ref / "ground_truth_manifest.json").read_text())
    assert ground_truth["input_provenance_status"] == "UNAVAILABLE"
    assert any("project_database_before" in message
               for message in ground_truth["input_provenance_messages"])


def test_modified_input_is_critical_provenance_failure(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    _make_run(run)
    ref = create_reference(run, tmp_path / "reference")
    approve_reference(ref, approved_by="reviewer", reason="Reviewed reference")
    (run / "inputs/markers.fasta").write_text(">QUB_0001\nTGCA\n", encoding="utf-8")
    result = validate_run(run, ref)
    assert result["status"] == "CRITICAL_FAILURE"
    assert result["provenance_status"] == "INVALID"


def test_production_mode_reports_added_and_removed_isolates_without_high_drift(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    _make_run(run)
    ref = create_reference(run, tmp_path / "reference")
    approve_reference(ref, approved_by="reviewer", reason="Reviewed reference")
    path = run / "04_assistant_cohort/01_joint_assistant/03_performance_review_hiring_panel/assessment/sequence_assessment.tsv"
    path.write_text(
        "ID\tMarkerQCClass\tRecommendation\n"
        "QUB_0001\tPASS_HIGH_CONFIDENCE\tPRIMARY\n"
        "QUB_0003\tPASS_WITH_WARNINGS\tSECONDARY\n", encoding="utf-8"
    )
    # A genuinely changed cohort input is declared and checksummed in the run manifest.
    marker_input = run / "inputs/markers.fasta"
    marker_input.write_text(">QUB_0001\nACGT\n>QUB_0003\nTTTT\n", encoding="utf-8")
    manifest_path = run / "04_assistant_cohort/01_joint_assistant/03_performance_review_hiring_panel/run_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    import hashlib
    manifest["inputs"][0]["sha256"] = hashlib.sha256(marker_input.read_bytes()).hexdigest()
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    result = validate_run(run, ref, mode="PRODUCTION_COMPATIBILITY")
    assert result["status"] == "PASS"
    assert result["provenance_status"] == "CHANGED"
    row_check = next(x for x in result["checks"] if x["check"] == "sequence_assessment.tsv:row_keys")
    assert row_check["severity"] == "INFO"


def test_database_row_content_change_is_detected_even_when_count_is_unchanged(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    _make_run(run)
    ref = create_reference(run, tmp_path / "reference")
    approve_reference(ref, approved_by="reviewer", reason="Reviewed reference")
    db = sqlite3.connect(run / "project.sqlite")
    db.execute("UPDATE sequences SET id='QUB_0003' WHERE id='QUB_0002'")
    db.commit()
    db.close()
    result = validate_run(run, ref)
    assert result["status"] == "DRIFT_DETECTED"
    assert any(x["check"] == "table_records:sequences" and x["status"] == "FAIL"
               for x in result["checks"])


def test_approved_reference_manifest_tampering_is_detected(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    _make_run(run)
    ref = create_reference(run, tmp_path / "reference")
    approve_reference(ref, approved_by="reviewer", reason="Reviewed reference")
    manifest_path = ref / "ground_truth_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["label"] = "tampered"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    result = validate_run(run, ref)
    assert result["status"] == "CRITICAL_FAILURE"
    assert any(x["check"] == "approval_integrity" and x["status"] == "FAIL"
               for x in result["checks"])


def test_production_mode_does_not_excuse_cohort_changes_when_inputs_are_unchanged(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    _make_run(run)
    ref = create_reference(run, tmp_path / "reference")
    approve_reference(ref, approved_by="reviewer", reason="Reviewed reference")
    path = run / "04_assistant_cohort/01_joint_assistant/03_performance_review_hiring_panel/assessment/sequence_assessment.tsv"
    path.write_text(
        "ID\tMarkerQCClass\tRecommendation\n"
        "QUB_0001\tPASS_HIGH_CONFIDENCE\tPRIMARY\n"
        "QUB_0003\tPASS_WITH_WARNINGS\tSECONDARY\n", encoding="utf-8"
    )
    result = validate_run(run, ref, mode="PRODUCTION_COMPATIBILITY")
    assert result["status"] == "DRIFT_DETECTED"


def test_current_output_with_blank_identifier_is_critical(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    _make_run(run)
    ref = create_reference(run, tmp_path / "reference")
    approve_reference(ref, approved_by="reviewer", reason="Reviewed reference")
    path = run / "04_assistant_cohort/01_joint_assistant/03_performance_review_hiring_panel/assessment/sequence_assessment.tsv"
    path.write_text("ID\tMarkerQCClass\tRecommendation\n\tPASS\tPRIMARY\n", encoding="utf-8")
    result = validate_run(run, ref)
    assert result["status"] == "CRITICAL_FAILURE"
    assert any(x["check"] == "sequence_assessment.tsv:unique_keys" and x["status"] == "FAIL"
               for x in result["checks"])


def test_manifest_input_collection_preserves_repeated_roles():
    from branchmanager.validation import _manifest_inputs

    result = _manifest_inputs({"inputs": [
        {"role": "alternate_reference", "path": "/refs/b.fasta", "sha256": "b"},
        {"role": "alternate_reference", "path": "/refs/a.fasta", "sha256": "a"},
    ]})
    assert list(result) == ["alternate_reference#1", "alternate_reference#2"]
    assert [item["sha256"] for item in result.values()] == ["a", "b"]


def test_database_discovery_requires_explicit_path_when_ambiguous(tmp_path):
    from branchmanager.validation import _resolve_database

    (tmp_path / "project.sqlite").touch()
    (tmp_path / "stage.sqlite").touch()
    with pytest.raises(ValueError, match="specify --db explicitly"):
        _resolve_database(tmp_path)
    assert _resolve_database(tmp_path, "project.sqlite") == (tmp_path / "project.sqlite").resolve()


def test_reference_records_missing_optional_outputs(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    _make_run(run)
    optional = run / "04_assistant_cohort/01_joint_assistant/03_performance_review_hiring_panel/assessment/novelty_metrics.tsv"
    optional.unlink()
    ref = create_reference(run, tmp_path / "reference")
    manifest = json.loads((ref / "ground_truth_manifest.json").read_text())
    assert "novelty_metrics.tsv" in manifest["missing_optional_outputs"]
    assert "novelty_metrics.tsv" not in manifest["outputs"]


def test_cli_exposes_database_and_comparison_mode_options():
    from branchmanager.cli.parser import build_parser

    parser = build_parser()
    args = parser.parse_args([
        "validate", "--run-dir", "run", "--reference", "reference",
        "--db", "project.sqlite", "--mode", "PRODUCTION_COMPATIBILITY",
    ])
    assert args.db == "project.sqlite"
    assert args.mode == "PRODUCTION_COMPATIBILITY"


def test_validation_does_not_modify_project_database(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    _make_run(run)
    ref = create_reference(run, tmp_path / "reference")
    approve_reference(ref, approved_by="reviewer", reason="Reviewed reference")
    db_path = run / "project.sqlite"
    before = db_path.read_bytes()
    validate_run(run, ref, database=db_path)
    assert db_path.read_bytes() == before
