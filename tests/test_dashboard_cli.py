
from pathlib import Path
from types import SimpleNamespace

import pytest

from branchmanager.cli.parser import build_parser
from branchmanager.cli_core import cmd_dashboard


def _make_run_directory(tmp_path):
    """Create a minimal existing BranchManager run directory."""
    run_dir = tmp_path / "test_run"
    assessment_dir = run_dir / "assessment"
    assessment_dir.mkdir(parents=True)

    (assessment_dir / "selection_summary.tsv").write_text(
        "isolate_id\tstatus\n"
        "ISO001\tPASS\n",
        encoding="utf-8",
    )

    # Optional inputs used by the dashboard generator.
    for filename in (
        "sequence_assessment.tsv",
        "pangenome_hierarchy.tsv",
        "novelty_metrics.tsv",
        "genome_collection_assessment.tsv",
    ):
        (assessment_dir / filename).write_text(
            "isolate_id\tstatus\nISO001\tPASS\n",
            encoding="utf-8",
        )

    (run_dir / "run_manifest.json").write_text(
        '{"status": "COMPLETE"}',
        encoding="utf-8",
    )

    return run_dir, assessment_dir


def test_dashboard_parser_accepts_run_directory():
    parser = build_parser()

    args = parser.parse_args(["dashboard", "review_output"])

    assert args.command == "dashboard"
    assert args.output_dir == "review_output"


def test_dashboard_parser_accepts_assessment_directory():
    parser = build_parser()

    args = parser.parse_args(["dashboard", "review_output/assessment"])

    assert args.command == "dashboard"
    assert args.output_dir == "review_output/assessment"


def test_dashboard_parser_accepts_read_only_database_overlay():
    parser = build_parser()
    args = parser.parse_args(["dashboard", "review_output", "--db", "project.sqlite"])
    assert args.db == "project.sqlite"


def test_dashboard_generates_from_run_directory(tmp_path, monkeypatch, capsys):
    run_dir, assessment_dir = _make_run_directory(tmp_path)
    expected_dashboard = run_dir / "performance_review_dashboard.html"
    calls = {}

    def fake_write_dashboard(selection_summary, output_dir, **kwargs):
        calls["selection_summary"] = selection_summary
        calls["output_dir"] = output_dir
        calls["kwargs"] = kwargs

        expected_dashboard.write_text(
            "<html><body>Test dashboard</body></html>",
            encoding="utf-8",
        )
        return expected_dashboard

    monkeypatch.setattr(
        "branchmanager.reporting.write_performance_review_dashboard",
        fake_write_dashboard,
    )

    result = cmd_dashboard(SimpleNamespace(output_dir=str(run_dir), db='project.sqlite'))

    assert result == expected_dashboard
    assert expected_dashboard.is_file()
    assert calls["selection_summary"] == assessment_dir / "selection_summary.tsv"
    assert calls["output_dir"] == run_dir
    assert calls["kwargs"]["database"] == 'project.sqlite'

    assert calls["kwargs"]["sequence_assessment"] == (
        assessment_dir / "sequence_assessment.tsv"
    )
    assert calls["kwargs"]["pangenome_hierarchy"] == (
        assessment_dir / "pangenome_hierarchy.tsv"
    )
    assert calls["kwargs"]["novelty_metrics"] == (
        assessment_dir / "novelty_metrics.tsv"
    )
    assert calls["kwargs"]["genome_assessment"] == (
        assessment_dir / "genome_collection_assessment.tsv"
    )
    assert calls["kwargs"]["run_manifest"] == run_dir / "run_manifest.json"

    captured = capsys.readouterr()
    assert "Dashboard regenerated:" in captured.out


def test_dashboard_generates_from_assessment_directory(tmp_path, monkeypatch):
    run_dir, assessment_dir = _make_run_directory(tmp_path)
    expected_dashboard = run_dir / "performance_review_dashboard.html"
    calls = {}

    def fake_write_dashboard(selection_summary, output_dir, **kwargs):
        calls["selection_summary"] = selection_summary
        calls["output_dir"] = output_dir
        expected_dashboard.write_text("<html></html>", encoding="utf-8")
        return expected_dashboard

    monkeypatch.setattr(
        "branchmanager.reporting.write_performance_review_dashboard",
        fake_write_dashboard,
    )

    result = cmd_dashboard(SimpleNamespace(output_dir=str(assessment_dir)))

    assert result == expected_dashboard
    assert calls["selection_summary"] == assessment_dir / "selection_summary.tsv"
    assert calls["output_dir"] == run_dir


def test_dashboard_missing_selection_summary_fails(tmp_path):
    run_dir = tmp_path / "incomplete_run"
    run_dir.mkdir()

    with pytest.raises(SystemExit, match="Could not find selection_summary.tsv"):
        cmd_dashboard(SimpleNamespace(output_dir=str(run_dir)))


def test_dashboard_nonexistent_directory_fails(tmp_path):
    missing_dir = tmp_path / "does_not_exist"

    with pytest.raises(SystemExit, match="Output directory does not exist"):
        cmd_dashboard(SimpleNamespace(output_dir=str(missing_dir)))


def test_dashboard_passes_none_for_missing_optional_inputs(
    tmp_path, monkeypatch
):
    run_dir = tmp_path / "partial_run"
    assessment_dir = run_dir / "assessment"
    assessment_dir.mkdir(parents=True)

    (assessment_dir / "selection_summary.tsv").write_text(
        "isolate_id\tstatus\nISO001\tPASS\n",
        encoding="utf-8",
    )

    calls = {}

    def fake_write_dashboard(selection_summary, output_dir, **kwargs):
        calls["kwargs"] = kwargs
        dashboard = output_dir / "performance_review_dashboard.html"
        dashboard.write_text("<html></html>", encoding="utf-8")
        return dashboard

    monkeypatch.setattr(
        "branchmanager.reporting.write_performance_review_dashboard",
        fake_write_dashboard,
    )

    cmd_dashboard(SimpleNamespace(output_dir=str(run_dir)))

    assert calls["kwargs"]["sequence_assessment"] is None
    assert calls["kwargs"]["pangenome_hierarchy"] is None
    assert calls["kwargs"]["novelty_metrics"] is None
    assert calls["kwargs"]["run_manifest"] is None
    assert calls["kwargs"]["genome_assessment"] is None

