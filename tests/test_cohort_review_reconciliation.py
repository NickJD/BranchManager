from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from branchmanager import cli_core


def _write_tsv(path: Path, rows: list[dict]) -> Path:
    fields = list(rows[0]) if rows else ['sequence_id', 'decision']
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    return path


def test_discovery_finds_unique_release_review_source(tmp_path):
    release = tmp_path / 'release'
    out = release / 'cohort' / 'assistant'
    out.mkdir(parents=True)
    review_dir = release / 'inputs' / 'reviews'
    review_dir.mkdir(parents=True)
    source = review_dir / 'cohort_marker_review_decisions_november20.tsv'
    source.write_text('sequence_id\tdecision\nA\tAPPROVE\n')
    assert cli_core._cohort_discover_review_decisions(out) == source.resolve()


def test_discovery_errors_if_multiple_sources(tmp_path):
    release = tmp_path / 'release'
    out = release / 'cohort'
    out.mkdir(parents=True)
    review_dir = release / 'inputs' / 'reviews'
    review_dir.mkdir(parents=True)
    for name in ('cohort_marker_review_decisions_a.tsv', 'cohort_marker_review_decisions_b.tsv'):
        (review_dir / name).write_text('sequence_id\tdecision\nA\tAPPROVE\n')
    with pytest.raises(SystemExit, match='Multiple review-decision TSVs'):
        cli_core._cohort_discover_review_decisions(out)


def test_reconciliation_passes_exact_ids_and_decisions(tmp_path):
    source = _write_tsv(tmp_path / 'source.tsv', [
        {'sequence_id': 'A', 'decision': 'APPROVE'},
        {'sequence_id': 'B', 'decision': 'REJECT'},
    ])
    template = _write_tsv(tmp_path / 'template.tsv', [
        {'sequence_id': 'A', 'decision': 'APPROVE'},
        {'sequence_id': 'B', 'decision': 'REJECT'},
    ])
    report = cli_core._cohort_reconcile_review_decisions(source, template, tmp_path)
    assert report['status'] == 'PASS'
    assert json.loads((tmp_path / 'cohort_review_reconciliation.json').read_text())['status'] == 'PASS'
    assert (tmp_path / 'cohort_review_reconciliation.tsv').is_file()


def test_reconciliation_fails_and_reports_mismatch_and_missing_ids(tmp_path):
    source = _write_tsv(tmp_path / 'source.tsv', [
        {'sequence_id': 'A', 'decision': 'APPROVE'},
        {'sequence_id': 'B', 'decision': 'REJECT'},
    ])
    template = _write_tsv(tmp_path / 'template.tsv', [
        {'sequence_id': 'A', 'decision': 'REJECT'},
        {'sequence_id': 'C', 'decision': 'APPROVE'},
    ])
    with pytest.raises(SystemExit, match='do not reconcile'):
        cli_core._cohort_reconcile_review_decisions(source, template, tmp_path)
    report = json.loads((tmp_path / 'cohort_review_reconciliation.json').read_text())
    assert report['status'] == 'FAIL'
    assert report['missing_from_template'] == ['B']
    assert report['unexpected_in_template'] == ['C']
    assert report['decision_mismatches'][0]['sequence_id'] == 'A'


def test_prepare_auto_restores_release_decisions_and_audits(tmp_path, monkeypatch):
    from test_v1_cohort import _cohort, _args
    _cohort(tmp_path)
    release_review_dir = tmp_path / 'inputs' / 'reviews'
    release_review_dir.mkdir(parents=True)
    source = _write_tsv(release_review_dir / 'cohort_marker_review_decisions_test.tsv', [
        {'sequence_id': 'P2', 'decision': 'REJECT', 'reviewer': 'Tester', 'notes': 'reviewed'},
        {'sequence_id': 'P3', 'decision': 'REJECT', 'reviewer': 'Tester', 'notes': 'reviewed'},
    ])
    monkeypatch.setattr(cli_core, 'cmd_assistant', lambda args: pytest.fail('assessment must not run in prepare-only'))
    cli_core.cmd_assistant_cohort(_args(tmp_path, '--cohort-manifest', str(tmp_path / 'cohort.tsv'),
                                     '--cohort-prepare-only'))
    prep = tmp_path / 'out' / '00_cohort_preparation'
    rows = {row['sequence_id']: row for row in cli_core._cohort_read_rows(prep / 'cohort_marker_review_template.tsv')}
    assert rows['P2']['decision'] == 'REJECT' and rows['P3']['decision'] == 'REJECT'
    assert json.loads((prep / 'cohort_review_reconciliation.json').read_text())['status'] == 'PASS'
    assert str(source.resolve()) in (prep / 'run_manifest.json').read_text()


def test_incremental_reconciliation_allows_new_isolates_and_blank_decisions(tmp_path):
    source = _write_tsv(tmp_path / 'source.tsv', [
        {'sequence_id': 'QUB_0001', 'decision': 'APPROVE'},
        {'sequence_id': 'BSI_0002', 'decision': 'RESEQUENCE'},
    ])
    template = _write_tsv(tmp_path / 'template.tsv', [
        {'sequence_id': 'QUB_0001', 'decision': 'APPROVE'},
        {'sequence_id': 'BSI_0002', 'decision': 'RESEQUENCE'},
        {'sequence_id': 'AU_0001', 'decision': ''},
        {'sequence_id': 'UNA_0001', 'decision': ''},
    ])
    report = cli_core._cohort_reconcile_review_decisions(source, template, tmp_path)
    assert report['status'] == 'PASS'
    assert report['carried_decisions'] == ['BSI_0002', 'QUB_0001']
    assert report['new_template_ids'] == ['AU_0001', 'UNA_0001']


def test_discovery_selects_latest_previous_populated_template(tmp_path):
    versions = tmp_path / 'versions'
    old = versions / '20261118_100000' / 'cohort' / '00_cohort_preparation'
    newer = versions / '20261119_100000' / 'cohort' / '00_cohort_preparation'
    current = versions / '20261120_100000' / 'cohort'
    old.mkdir(parents=True)
    newer.mkdir(parents=True)
    current.mkdir(parents=True)
    _write_tsv(old / 'cohort_marker_review_template.tsv', [
        {'sequence_id': 'A', 'decision': 'REJECT'}])
    _write_tsv(newer / 'cohort_marker_review_template.tsv', [
        {'sequence_id': 'A', 'decision': 'APPROVE'}])
    assert cli_core._cohort_discover_review_decisions(current) == (
        newer / 'cohort_marker_review_template.tsv').resolve()
