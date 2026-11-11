from __future__ import annotations

import ast
import csv
import sys
import types
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

@pytest.fixture(autouse=True)
def _restore_branchmanager_imports():
    """Keep test-local sys.modules stubs from leaking into other test modules."""
    import sys
    before = {
        name: module for name, module in sys.modules.items()
        if name == 'branchmanager' or name.startswith('branchmanager.')
    }
    yield
    for name in list(sys.modules):
        if (name == 'branchmanager' or name.startswith('branchmanager.')) and name not in before:
            sys.modules.pop(name, None)
    sys.modules.update(before)



def _write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def _read_fasta(path: str):
    records = []
    header = None
    sequence = []
    with open(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(sequence)))
                header, sequence = line[1:], []
            else:
                sequence.append(line)
    if header is not None:
        records.append((header, ''.join(sequence)))
    return records


def _load_cli_completeness_helper():
    source = (ROOT / 'src' / 'branchmanager' / 'cli_core.py').read_text()
    tree = ast.parse(source)
    names = {'_cohort_read_rows', '_cohort_row_get', '_cohort_paper_trail_outputs_complete'}
    nodes = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in names]
    assert {node.name for node in nodes} == names

    branchmanager = types.ModuleType('branchmanager')
    branchmanager.__path__ = []
    sys.modules['branchmanager'] = branchmanager
    utils = types.ModuleType('branchmanager.utils')
    utils.__path__ = []
    sys.modules['branchmanager.utils'] = utils
    fasta = types.ModuleType('branchmanager.utils.fasta')
    fasta.read_fasta = _read_fasta
    sys.modules['branchmanager.utils.fasta'] = fasta

    namespace = {'csv': csv, 'Path': Path}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), 'cli.py:review-output-helper', 'exec'), namespace)
    return namespace['_cohort_paper_trail_outputs_complete']


def test_paper_trail_keeps_rows_for_review_template_and_visual_reports():
    source = (ROOT / 'src' / 'branchmanager' / 'pipeline' / 'paper_trail.py').read_text()
    tree = ast.parse(source)
    function = next(
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == 'run_paper_trail'
    )
    append_targets = set()
    for node in ast.walk(function):
        if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Attribute):
            continue
        if node.func.attr != 'append' or not isinstance(node.func.value, ast.Name):
            continue
        if node.args and isinstance(node.args[0], ast.Name) and node.args[0].id == 'row':
            append_targets.add(node.func.value.id)
    assert {'assembly_rows', 'recommendation_rows'} <= append_targets


def test_cache_validator_rejects_header_only_template_and_accepts_complete_handoff(tmp_path):
    validate = _load_cli_completeness_helper()
    trace = tmp_path / 'paper_trail'
    failed = trace / 'failed_qc_sequences'
    failed.mkdir(parents=True)

    assembly_fields = ['SequenceID', 'QCClass', 'Recommendation', 'OutputLength']
    assembly_rows = [
        {'SequenceID': 'QUB_0327', 'QCClass': 'FAIL_QC', 'Recommendation': 'RESEQUENCE', 'OutputLength': '0'},
        {'SequenceID': 'QUB_0316', 'QCClass': 'FAIL_QC', 'Recommendation': 'RESEQUENCE', 'OutputLength': '0'},
        {'SequenceID': 'QUB_0172', 'QCClass': 'PASS_WITH_WARNINGS', 'Recommendation': 'MANUAL_REVIEW', 'OutputLength': '1123'},
        {'SequenceID': 'QUB_PASS', 'QCClass': 'PASS_HIGH_CONFIDENCE', 'Recommendation': 'ACCEPT', 'OutputLength': '1500'},
    ]
    _write_tsv(trace / 'assembly_report.tsv', assembly_fields, assembly_rows)
    _write_tsv(failed / 'failed_qc_manifest.tsv', ['SequenceID', 'QCClass'], [])
    _write_tsv(failed / 'failed_read_manifest.tsv', ['SequenceID', 'ReadID'], [])
    (failed / 'failed_review_candidates.fasta').write_text(
        '>QUB_0327\n' + 'A' * 987 + '\n' +
        '>QUB_0316\n' + 'A' * 45 + '\n'
    )

    review_fields = ['sequence_id', 'review_type', 'qc_class', 'eligible_for_manual_acceptance', 'candidate_length']
    review_rows = [
        {'sequence_id': 'QUB_0327', 'review_type': 'FAILED_QC_REVIEW', 'qc_class': 'FAIL_QC', 'eligible_for_manual_acceptance': 'yes', 'candidate_length': '987'},
        {'sequence_id': 'QUB_0316', 'review_type': 'FAILED_QC_REVIEW', 'qc_class': 'FAIL_QC', 'eligible_for_manual_acceptance': 'yes', 'candidate_length': '45'},
        {'sequence_id': 'QUB_0172', 'review_type': 'WARNING_REVIEW', 'qc_class': 'PASS_WITH_WARNINGS', 'eligible_for_manual_acceptance': 'yes', 'candidate_length': '1123'},
    ]
    template = trace / 'marker_review_template.tsv'
    _write_tsv(template, review_fields, [])

    complete, reason = validate(trace)
    assert not complete
    assert '3 of 3 required review records missing' in reason

    _write_tsv(template, review_fields, review_rows)
    complete, reason = validate(trace)
    assert complete, reason

    # An approval-eligible row must correspond to a real sequence in the FASTA.
    (failed / 'failed_review_candidates.fasta').write_text('>QUB_0327\n' + 'A' * 987 + '\n')
    complete, reason = validate(trace)
    assert not complete
    assert 'QUB_0316' in reason


def test_fallback_read_review_rows_are_eligible_even_without_consensus_length(tmp_path):
    """Mirror the real QUB_0327 case: output_length 0, candidate_length 987."""
    validate = _load_cli_completeness_helper()
    trace = tmp_path / 'paper_trail'
    failed = trace / 'failed_qc_sequences'
    failed.mkdir(parents=True)
    _write_tsv(
        trace / 'assembly_report.tsv',
        ['SequenceID', 'QCClass', 'Recommendation', 'OutputLength'],
        [{'SequenceID': 'QUB_0327', 'QCClass': 'FAIL_QC', 'Recommendation': 'RESEQUENCE', 'OutputLength': '0'}],
    )
    _write_tsv(
        trace / 'marker_review_template.tsv',
        ['sequence_id', 'review_type', 'qc_class', 'eligible_for_manual_acceptance', 'candidate_length'],
        [{'sequence_id': 'QUB_0327', 'review_type': 'FAILED_QC_REVIEW', 'qc_class': 'FAIL_QC', 'eligible_for_manual_acceptance': 'yes', 'candidate_length': '987'}],
    )
    _write_tsv(failed / 'failed_qc_manifest.tsv', ['SequenceID', 'QCClass'], [])
    _write_tsv(failed / 'failed_read_manifest.tsv', ['SequenceID', 'ReadID'], [])
    (failed / 'failed_review_candidates.fasta').write_text('>QUB_0327\n' + 'A' * 987 + '\n')
    complete, reason = validate(trace)
    assert complete, reason
