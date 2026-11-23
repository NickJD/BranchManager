from pathlib import Path
import json
from branchmanager.release_manager import ReleaseManager
from branchmanager.cli.parser import build_parser

def make_release(tmp_path):
    root = tmp_path / "project"
    root.mkdir()
    cfg = root / "project_config.json"
    cfg.write_text(json.dumps({"schema_version": 1, "project_directory": str(root), "data": {}, "defaults": {}}))
    manager = ReleaseManager(root)
    manager.init_project("Test", config_source=cfg, overwrite=True)
    release = manager.create_release("r_phase2")
    return manager, Path(release["release_dir"])

def test_release_workflow_event_records_input_hashes_and_status(tmp_path):
    manager, release = make_release(tmp_path)
    input_file = tmp_path / "input.tsv"
    input_file.write_text("sequence_id\tvalue\nS1\t1\n")
    entry = manager.record_workflow_event(release, workflow="registration", command="import-partner-workbook",
        argv=[str(input_file)], status="STARTED", event="workflow_started")
    assert entry["status"] == "STARTED"
    assert entry["input_files"][0]["sha256"]
    manifest = manager.get_release_manifest(release)
    assert manifest["workflow_runs"][0]["run_id"] == entry["run_id"]

def test_release_mode_makes_database_and_output_optional_for_workflows():
    parser = build_parser()
    args = parser.parse_args(["filing-cabinet", "--release", "r1", "--build-tree"])
    assert args.release_context == "r1"
    assert args.db is None and args.out == '.'

def test_continue_project_parser_accepts_release_context_without_db_out(tmp_path):
    parser = build_parser()
    args = parser.parse_args(["continue-project", "--release", "r1",
        "--partner-catalogue", "partners.csv", "--genome-results", "genomes.tsv"])
    assert args.release_context == "r1"
    assert args.db is None and args.out is None and args.version_root is None

def test_continue_project_release_context_resolves_db_out_and_version_root(tmp_path, monkeypatch):
    manager, release = make_release(tmp_path)
    catalogue = tmp_path / "catalogue.csv"
    genomes = tmp_path / "genomes.tsv"
    catalogue.write_text("isolate_id\nX1\n")
    genomes.write_text("sequence_id\tgenome_id\nX1\tG1\n")
    captured = {}
    import branchmanager.cli_core as core
    def fake_continue(args):
        captured.update(db=args.db, out=args.out, version_root=args.version_root)
    monkeypatch.setattr(core, "cmd_continue_project", fake_continue)
    core.main(["continue-project", "--release", str(release),
        "--partner-catalogue", str(catalogue), "--genome-results", str(genomes)])
    assert Path(captured["db"]).resolve() == (release / "project.sqlite").resolve()
    assert Path(captured["out"]).resolve() == (release / "continuation").resolve()
    assert Path(captured["version_root"]).resolve() == (release / "continuation_versions").resolve()
    manifest = manager.get_release_manifest(release)
    assert manifest["workflow_status"]["continuation"] == "COMPLETE"
    assert [r["status"] for r in manifest["workflow_runs"]] == ["STARTED", "COMPLETE"]

def test_continue_project_release_failure_records_failed_status(tmp_path, monkeypatch):
    manager, release = make_release(tmp_path)
    import branchmanager.cli_core as core
    def fail_continue(args):
        raise RuntimeError("simulated failure")
    monkeypatch.setattr(core, "cmd_continue_project", fail_continue)
    import pytest
    with pytest.raises(RuntimeError, match="simulated failure"):
        core.main(["continue-project", "--release", str(release),
            "--partner-catalogue", "partners.csv", "--genome-results", "genomes.tsv"])
    manifest = manager.get_release_manifest(release)
    assert manifest["status"] == "FAILED"
    assert manifest["workflow_status"]["continuation"] == "FAILED"
    assert manifest["workflow_runs"][-1]["status"] == "FAILED"


def test_concurrent_manifest_events_are_not_lost(tmp_path):
    from concurrent.futures import ThreadPoolExecutor
    manager, release = make_release(tmp_path)
    def write_event(i):
        return manager.record_workflow_event(release, workflow="registration", command="test-event",
            argv=[str(i)], status="STARTED", event="concurrent_test_event")
    with ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(write_event, range(20)))
    manifest = manager.get_release_manifest(release)
    assert len(manifest["workflow_runs"]) == 20
    assert len({row["run_id"] for row in manifest["workflow_runs"]}) == 20


def test_manual_workflow_status_accepts_continuation():
    parser = build_parser()
    args = parser.parse_args(["release", "status", "--project", "/tmp/project", "--release", "r1",
        "--set", "IN_PROGRESS", "--workflow", "continuation"])
    assert args.workflow == "continuation"
