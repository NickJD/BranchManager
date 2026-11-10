from pathlib import Path
import shlex

from branchmanager.pipeline import novelty
from branchmanager.utils.fasta import write_fasta


def _row(query, target, identity, qcov):
    # query,target,id,alnlen,mism,gaps,qlo,qhi,tlo,thi,ql,tl,qcov,tcov
    return '\t'.join(map(str, [query, target, identity, 100, 0, 0, 1, 100, 1, 100, 100, 100, qcov, 100])) + '\n'


def test_density_pair_checkpoints_batches_and_reports_bounded_lower_bound_counts(tmp_path, monkeypatch):
    query = tmp_path / 'queries.fasta'
    reference = tmp_path / 'reference.fasta'
    write_fasta([(f'q{i}', 'A' * 100) for i in range(5)], str(query))
    write_fasta([('r1', 'A' * 100)], str(reference))
    calls = []

    def fake_audited(command, **kwargs):
        calls.append(command)
        tokens = shlex.split(command)
        query_path = Path(tokens[tokens.index('--usearch_global') + 1])
        output_path = Path(tokens[tokens.index('--userout') + 1])
        query_records = list(novelty.read_fasta(str(query_path)))
        with output_path.open('w') as handle:
            for qid, _seq in query_records:
                handle.write(_row(qid, 'r95', 95.0, 95.0))
                handle.write(_row(qid, 'r97', 97.0, 90.0))
                handle.write(_row(qid, 'r99', 99.0, 80.0))
        return None

    monkeypatch.setattr(novelty, 'run_vsearch_audited', fake_audited)
    output = tmp_path / 'novelty_reference_density_matches.tsv'
    qcov95, qcov80 = novelty._run_density_search_coverage_pair(
        str(query), str(reference), output, batch_size=2,
    )
    assert len(calls) == 3
    assert qcov80['q0'] == {'ge_99': 1, 'ge_97': 2, 'ge_95': 3}
    assert qcov95['q0'] == {'ge_99': 0, 'ge_97': 0, 'ge_95': 1}
    assert sum(1 for _ in output.open()) == 5
    assert all('--maxaccepts 100' in command and '--maxrejects 256' in command and '--maxhits 100' in command for command in calls)

    # A second invocation reuses all verified batch outputs.
    qcov95_again, qcov80_again = novelty._run_density_search_coverage_pair(
        str(query), str(reference), output, batch_size=2,
    )
    assert len(calls) == 3
    assert qcov95_again == qcov95
    assert qcov80_again == qcov80

    # Corrupting one checkpoint invalidates only that batch.
    checkpoint_dir = next(tmp_path.glob('.novelty_reference_density_matches_density_checkpoints_*'))
    (checkpoint_dir / 'hits_0002.tsv').write_text('truncated\n')
    qcov95_repaired, qcov80_repaired = novelty._run_density_search_coverage_pair(
        str(query), str(reference), output, batch_size=2,
    )
    assert len(calls) == 4
    assert qcov95_repaired == qcov95
    assert qcov80_repaired == qcov80


def test_density_metrics_are_explicitly_lower_bounds():
    metrics = novelty._metrics_from_nearest_and_density(
        ['q1'], {'q1': (99.0, 'r1', 100.0, 100)},
        {'q1': {'ge_99': 3, 'ge_97': 5, 'ge_95': 8}},
        source='reference', id_threshold=0.97,
    )
    assert metrics['q1']['density_count_method'] == 'LOWER_BOUND_MAX_100_HITS'


def test_novelty_metrics_defer_pangenome_gap_until_tree_validation(tmp_path):
    output = tmp_path / 'novelty_metrics.tsv'
    empty = novelty._empty_pool_metrics(['q1'])
    novelty._write_novelty_metrics_table(
        output, ['q1'], empty, empty, empty, empty, empty, empty,
        sequencing_context=novelty._empty_sequencing_context(['q1']),
    )
    lines = output.read_text().splitlines()
    header = lines[0].split('\t')
    values = dict(zip(header, lines[1].split('\t')))
    assert len(header) == len(lines[1].split('\t'))
    assert values['PangenomeTarget'] == 'NA'
    assert values['PangenomeGap'] == 'NA'
    assert values['PangenomeEligibilityStatus'] == 'DEFERRED_TO_TREE_VALIDATED_SELECTION_SETS'


def test_combined_baseline_density_is_also_bounded(tmp_path, monkeypatch):
    query = tmp_path / 'queries.fasta'
    reference = tmp_path / 'reference.fasta'
    write_fasta([('q1', 'A' * 100)], str(query))
    write_fasta([('h1', 'A' * 100), ('s1', 'A' * 100)], str(reference))
    calls = []

    def fake_audited(command, **kwargs):
        calls.append(command)
        output_path = Path(shlex.split(command)[shlex.split(command).index('--userout') + 1])
        output_path.write_text(_row('q1', 'h1', 99.0, 100.0) + _row('q1', 's1', 97.0, 90.0))
        return None

    monkeypatch.setattr(novelty, 'run_vsearch_audited', fake_audited)
    result = novelty._run_combined_baseline_density(
        str(query), str(reference), str(tmp_path),
        {'h1': 'HUNGATE', 's1': 'SECONDARY'},
    )
    assert len(calls) == 1
    assert '--maxaccepts 100' in calls[0]
    assert '--maxrejects 256' in calls[0]
    assert '--maxhits 100' in calls[0]
    assert result[80]['HUNGATE']['q1']['ge_99'] == 1
    assert result[80]['SECONDARY']['q1']['ge_97'] == 1
