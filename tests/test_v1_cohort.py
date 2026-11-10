"""v1 mixed-partner cohort prepare / review / resume tests."""
from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import pytest

from branchmanager import cli
from branchmanager.db.interface import Database
from branchmanager.utils.fasta import read_fasta

SEQ = 'ACGTTGCA' * 180


def _tsv(path, rows):
    with open(path, 'w', newline='') as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    return path


def _cohort(tmp_path):
    db = Database(str(tmp_path / 'live.sqlite'))
    db.initialise()
    db.upsert_dataset_role('Base', 'baseline', genomes_available=False)
    db.insert_sequences([('B1', SEQ)], dataset='Base')
    (tmp_path / 'partner_metadata.tsv').write_text('partner_id\tsequence_id\n')
    inputs = tmp_path / 'in'
    inputs.mkdir()
    (inputs / 'p.fasta').write_text(f'>P1\n{SEQ}\n>P2\n{SEQ[::-1]}\n>P3\n{SEQ[:-8]}ACGTACGT\n')
    _tsv(inputs / 'p_qc.tsv', [
        {'SequenceID': 'P1', 'QCClass': 'PASS_HIGH_CONFIDENCE', 'Recommendation': 'ACCEPT', 'Reasons': ''},
        {'SequenceID': 'P2', 'QCClass': 'FAIL_QC', 'Recommendation': 'RESEQUENCE', 'Reasons': 'low_quality'},
    ])
    _tsv(inputs / 'maldi.tsv', [{'sequence_id': 'M1', 'partner_id': 'PM', 'taxon_name': 'Bacillus subtilis',
                                 'taxonomic_rank': 'species'}])
    (tmp_path / 'cohort.tsv').write_text(
        'partner_id\tdataset\tinput_type\tmarker_fasta\tmarker_qc\tmaldi_file\n'
        'PA\tPartnerA\tpartner_fasta\tin/p.fasta\tin/p_qc.tsv\t\n'
        'PM\tPartnerM\t\t\t\tin/maldi.tsv\n'
    )
    return db


def _args(tmp_path, *extra):
    parser = cli.build_parser()
    return parser.parse_args([
        'assistant', '--partner-metadata', str(tmp_path / 'partner_metadata.tsv'),
        '--db', str(tmp_path / 'live.sqlite'), '--dataset', 'PartnerA', '--ref', str(tmp_path / 'ref.fa'),
        '-o', str(tmp_path / 'out'), *extra,
    ])


def _digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def test_manifest_requires_only_modality_fields_and_resolves_relative_paths(tmp_path):
    _cohort(tmp_path)
    rows = cli._read_assistant_cohort_manifest(tmp_path / 'cohort.tsv')
    assert rows[0]['input_type'] == 'partner_fasta' and rows[1]['input_type'] == 'maldi_only'
    assert rows[0]['marker_fasta'] == str((tmp_path / 'in' / 'p.fasta').resolve())
    bad = tmp_path / 'bad.tsv'
    bad.write_text('partner_id\tdataset\tinput_type\nPA\tD\tsanger\n')
    with pytest.raises(ValueError, match='sample_map'):
        cli._read_assistant_cohort_manifest(bad)


def test_cohort_prepare_does_not_touch_live_db_or_run_panel(tmp_path, monkeypatch):
    _cohort(tmp_path)
    before = _digest(tmp_path / 'live.sqlite')
    monkeypatch.setattr(cli, 'cmd_assistant', lambda a: pytest.fail('assistant must not run in prepare'))
    cli.cmd_assistant_cohort(_args(tmp_path, '--cohort-manifest', str(tmp_path / 'cohort.tsv'),
                                   '--cohort-prepare-only'))
    prep = tmp_path / 'out' / '00_cohort_preparation'
    assert _digest(tmp_path / 'live.sqlite') == before
    for name in ('cohort_passing_markers.fasta', 'cohort_recoverable_failed_markers.fasta',
                 'cohort_marker_qc.tsv', 'cohort_marker_review_template.tsv',
                 'cohort_marker_evidence_manifest.tsv', 'cohort_maldi_import_report.tsv',
                 'cohort_genome_import_report.tsv', 'cohort_prepare_state.json'):
        assert (prep / name).is_file(), name
    passing = {str(h).split()[0] for h, _ in read_fasta(str(prep / 'cohort_passing_markers.fasta'))}
    assert 'P1' in passing and 'P2' not in passing and 'M1' not in passing  # MALDI never in marker FASTA
    template = {r['sequence_id']: r for r in csv.DictReader(open(prep / 'cohort_marker_review_template.tsv'), delimiter='\t')}
    assert template['P2']['review_type'] == 'FAILED_QC_REVIEW' and template['P3']['review_type'] == 'UNVERIFIED_QC_REVIEW'
    recoverable = {str(h).split()[0] for h, _ in read_fasta(str(prep / 'cohort_recoverable_failed_markers.fasta'))}
    assert recoverable == {'P2'}


def _prepare(tmp_path, monkeypatch):
    _cohort(tmp_path)
    cli.cmd_assistant_cohort(_args(tmp_path, '--cohort-manifest', str(tmp_path / 'cohort.tsv'),
                                   '--cohort-prepare-only'))


def test_cohort_resume_joint_assessment_and_publish(tmp_path, monkeypatch):
    _prepare(tmp_path, monkeypatch)
    seen = {}

    def fake_assistant(a):
        seen['ids'] = {str(h).split()[0] for h, _ in read_fasta(a.fasta)}
        seen['db'] = a.db
        seen['live_ids_during'] = set(Database(str(tmp_path / 'live.sqlite')).get_all_ids())

    monkeypatch.setattr(cli, 'cmd_assistant', fake_assistant)
    decisions = _tsv(tmp_path / 'decisions.tsv', [
        {'sequence_id': 'P2', 'decision': 'APPROVE_FAIL_QC', 'reviewer': 'NJ', 'notes': 'trace inspected',
         'reviewed_at': '2025-01-01'},
        {'sequence_id': 'P3', 'decision': 'REJECT', 'reviewer': 'NJ', 'notes': 'no provenance',
         'reviewed_at': '2025-01-01'},
    ])
    cli.cmd_assistant_cohort(_args(tmp_path, '--cohort-resume', str(tmp_path / 'out'),
                                   '--marker-review', str(decisions)))
    assert seen['ids'] == {'P1', 'P2'}  # rejected P3 and MALDI-only M1 are not assessed
    assert 'M1' not in seen['live_ids_during']  # live DB untouched until publish
    final = Database(str(tmp_path / 'live.sqlite'))
    ids = set(final.get_all_ids())
    assert {'P1', 'M1', 'B1'} <= ids
    with final.connect() as conn:
        assert conn.execute('SELECT COUNT(*) FROM maldi_identifications').fetchone()[0] == 1
        decisions_logged = conn.execute(
            "SELECT review_decision,manual_review_status FROM sequence_provenance WHERE sequence_id='P3'").fetchone()
    assert decisions_logged == ('REJECT', 'REJECTED')
    assert list((tmp_path / 'out' / 'review_history').glob('*decisions.tsv'))


def test_invalid_decision_aborts_before_any_change(tmp_path, monkeypatch):
    _prepare(tmp_path, monkeypatch)
    before = _digest(tmp_path / 'live.sqlite')
    monkeypatch.setattr(cli, 'cmd_assistant', lambda a: pytest.fail('must not run'))
    decisions = _tsv(tmp_path / 'decisions.tsv', [
        {'sequence_id': 'P2', 'decision': 'APPROVE_FAIL_QC', 'reviewer': '', 'notes': ''}])
    with pytest.raises(SystemExit, match='Invalid manual-review'):
        cli.cmd_assistant_cohort(_args(tmp_path, '--cohort-resume', str(tmp_path / 'out'),
                                       '--marker-review', str(decisions)))
    assert _digest(tmp_path / 'live.sqlite') == before


def test_failure_during_assessment_leaves_live_db_unchanged(tmp_path, monkeypatch):
    _prepare(tmp_path, monkeypatch)
    before = _digest(tmp_path / 'live.sqlite')

    def boom(a):
        Database(a.db).insert_sequences([('X1', SEQ)], dataset='PartnerA')
        raise RuntimeError('simulated hiring panel failure')

    monkeypatch.setattr(cli, 'cmd_assistant', boom)
    with pytest.raises(RuntimeError):
        cli.cmd_assistant_cohort(_args(tmp_path, '--cohort-resume', str(tmp_path / 'out')))
    assert _digest(tmp_path / 'live.sqlite') == before
    assert 'X1' not in set(Database(str(tmp_path / 'live.sqlite')).get_all_ids())


def test_resume_refuses_when_prepared_outputs_changed(tmp_path, monkeypatch):
    _prepare(tmp_path, monkeypatch)
    prep = tmp_path / 'out' / '00_cohort_preparation'
    (prep / 'cohort_passing_markers.fasta').write_text('>P1\nACGT\n')
    with pytest.raises(SystemExit, match='changed since'):
        cli.cmd_assistant_cohort(_args(tmp_path, '--cohort-resume', str(tmp_path / 'out')))


def test_prepare_preserves_existing_review_decisions_on_rerun(tmp_path, monkeypatch):
    _prepare(tmp_path, monkeypatch)
    template = tmp_path / 'out' / '00_cohort_preparation' / 'cohort_marker_review_template.tsv'
    rows = list(csv.DictReader(open(template), delimiter='\t'))
    for r in rows:
        if r['sequence_id'] == 'P2':
            r.update(decision='REJECT', reviewer='NJ', notes='bad')
    _tsv(template, rows)
    cli.cmd_assistant_cohort(_args(tmp_path, '--cohort-manifest', str(tmp_path / 'cohort.tsv'),
                                   '--cohort-prepare-only'))
    again = {r['sequence_id']: r for r in csv.DictReader(open(template), delimiter='\t')}
    assert again['P2']['decision'] == 'REJECT' and again['P2']['reviewer'] == 'NJ'
