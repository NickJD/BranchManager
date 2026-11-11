from __future__ import annotations

import ast
import csv
import json
import sys
import types
import pytest
from tempfile import TemporaryDirectory
from pathlib import Path

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



def load_paper_trail_module():
    """Load Paper Trail with minimal FASTA stubs for this focused unit test."""
    package_root = ROOT / 'src' / 'branchmanager'
    branchmanager = types.ModuleType('branchmanager')
    branchmanager.__path__ = [str(package_root)]
    sys.modules['branchmanager'] = branchmanager

    pipeline = types.ModuleType('branchmanager.pipeline')
    pipeline.__path__ = [str(package_root / 'pipeline')]
    sys.modules['branchmanager.pipeline'] = pipeline

    utils = types.ModuleType('branchmanager.utils')
    utils.__path__ = []
    sys.modules['branchmanager.utils'] = utils

    fasta = types.ModuleType('branchmanager.utils.fasta')
    fasta.read_fasta = lambda *args, **kwargs: []
    fasta.write_fasta = lambda *args, **kwargs: None
    fasta.reverse_complement = lambda sequence: sequence[::-1]
    sys.modules['branchmanager.utils.fasta'] = fasta

    name = 'branchmanager.pipeline.paper_trail'
    spec = __import__('importlib.util').util.spec_from_file_location(
        name, package_root / 'pipeline' / 'paper_trail.py'
    )
    module = __import__('importlib.util').util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _read_fasta(path):
    records = []
    header = None
    sequence = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(sequence)))
                header = line[1:]
                sequence = []
            else:
                sequence.append(line)
    if header is not None:
        records.append((header, ''.join(sequence)))
    return records


def load_qc_module():
    """Load qc.py with small FASTA read/write helpers for the integration check."""
    package_root = ROOT / 'src' / 'branchmanager'
    branchmanager = types.ModuleType('branchmanager')
    branchmanager.__path__ = [str(package_root)]
    sys.modules['branchmanager'] = branchmanager
    pipeline = types.ModuleType('branchmanager.pipeline')
    pipeline.__path__ = [str(package_root / 'pipeline')]
    sys.modules['branchmanager.pipeline'] = pipeline
    utils = types.ModuleType('branchmanager.utils')
    utils.__path__ = []
    sys.modules['branchmanager.utils'] = utils
    fasta = types.ModuleType('branchmanager.utils.fasta')
    fasta.read_fasta = _read_fasta

    def write_fasta(records, path):
        with open(path, 'w') as handle:
            for header, sequence in records:
                handle.write(f'>{header}\n')
                handle.write(f'{sequence}\n')

    fasta.write_fasta = write_fasta
    sys.modules['branchmanager.utils.fasta'] = fasta
    name = 'branchmanager.pipeline.qc'
    spec = __import__('importlib.util').util.spec_from_file_location(name, ROOT / 'tests' / 'fixtures' / 'qc.py')
    module = __import__('importlib.util').util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def load_studio_module():
    """Load Sanger Studio with lightweight Paper Trail/manifest stubs."""
    package_root = ROOT / 'src' / 'branchmanager'
    branchmanager = types.ModuleType('branchmanager')
    branchmanager.__path__ = [str(package_root)]
    sys.modules['branchmanager'] = branchmanager
    pipeline = types.ModuleType('branchmanager.pipeline')
    pipeline.__path__ = [str(package_root / 'pipeline')]
    sys.modules['branchmanager.pipeline'] = pipeline
    paper_trail = types.ModuleType('branchmanager.pipeline.paper_trail')
    sys.modules['branchmanager.pipeline.paper_trail'] = paper_trail
    pipeline.paper_trail = paper_trail
    run_manifest = types.ModuleType('branchmanager.run_manifest')
    run_manifest.RunManifest = type('RunManifest', (), {})
    sys.modules['branchmanager.run_manifest'] = run_manifest
    name = 'branchmanager.sanger_studio'
    spec = __import__('importlib.util').util.spec_from_file_location(
        name, package_root / 'sanger_studio.py'
    )
    module = __import__('importlib.util').util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def test_sanger_studio_offers_and_saves_failed_qc_manual_approval(tmp_path):
    studio_module = load_studio_module()
    studio = studio_module.SangerStudio(tmp_path, tmp_path / 'studio_output')
    review_path = tmp_path / 'cohort_marker_review_template.tsv'
    _write_tsv(
        review_path,
        [
            'sequence_id', 'decision', 'reviewer', 'reviewed_at', 'notes',
            'review_type', 'eligible_for_manual_acceptance', 'allowed_decisions',
            'qc_class', 'output_length', 'n_percent',
        ],
        {
            'sequence_id': 'QUB_0047', 'decision': '', 'reviewer': '', 'reviewed_at': '',
            'notes': '', 'review_type': 'FAILED_QC_REVIEW',
            'eligible_for_manual_acceptance': 'yes',
            'allowed_decisions': 'APPROVE_FAIL_QC|REJECT|RESEQUENCE',
            'qc_class': 'FAIL_QC', 'output_length': '760', 'n_percent': '0.921',
        },
    )
    # Append a second failed row without a candidate: it remains reject/resequence only.
    with open(review_path, 'a') as handle:
        handle.write(
            'QUB_EMPTY\t\t\t\t\tFAILED_QC_REVIEW\tno\tREJECT|RESEQUENCE\t'
            'FAIL_QC\t0\t0\n'
        )

    state = studio.load_marker_review(review_path)
    rows = {row['sequence_id']: row for row in state['rows']}
    assert rows['QUB_0047']['allowed_decisions_list'] == [
        'APPROVE_FAIL_QC', 'REJECT', 'RESEQUENCE'
    ]
    assert rows['QUB_EMPTY']['allowed_decisions_list'] == ['REJECT', 'RESEQUENCE']

    saved = studio.save_marker_review({
        'reviewer': 'Nick',
        'rows': [
            {'sequence_id': 'QUB_0047', 'decision': 'APPROVE_FAIL_QC',
             'notes': 'Reviewed chromatogram; accepted for this marker analysis.'},
            {'sequence_id': 'QUB_EMPTY', 'decision': '', 'notes': ''},
        ],
    })
    saved_path = Path(saved['saved'])
    with open(saved_path, newline='') as handle:
        written = {row['sequence_id']: row for row in csv.DictReader(handle, delimiter='\t')}
    assert written['QUB_0047']['decision'] == 'APPROVE_FAIL_QC'
    assert written['QUB_0047']['reviewer'] == 'Nick'
    assert written['QUB_0047']['notes'].startswith('Reviewed chromatogram')
    assert written['QUB_EMPTY']['decision'] == ''


def load_cli_helpers():
    """Compile the relevant pure CLI helpers without importing full BranchManager."""
    source = (ROOT / 'src' / 'branchmanager' / 'cli_core.py').read_text()
    tree = ast.parse(source)
    names = {
        '_cohort_read_rows', '_cohort_row_get', '_cohort_failed_qc_candidate',
        '_cohort_validate_decisions',
    }
    nodes = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in names
    ]
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
    evidence = types.ModuleType('branchmanager.marker_evidence')
    evidence.REVIEW_DECISIONS = (
        'APPROVE', 'APPROVE_FAIL_QC', 'APPROVE_UNVERIFIED', 'REJECT', 'RESEQUENCE'
    )

    def validate_review_decision(decision, qc_class, *, reviewer='', notes='', has_consensus=True):
        if decision.startswith('APPROVE') and (not reviewer or not notes or not has_consensus):
            return 'approval requires reviewer, notes and a candidate sequence'
        if decision == 'APPROVE_FAIL_QC' and qc_class != 'FAIL_QC':
            return 'APPROVE_FAIL_QC requires FAIL_QC'
        return ''

    evidence.validate_review_decision = validate_review_decision
    sys.modules['branchmanager.marker_evidence'] = evidence

    namespace = {'csv': csv, 'json': json, 'Path': Path, '_COHORT_STATE_NAME': 'cohort_prepare_state.json'}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), 'cli.py:test-helpers', 'exec'), namespace)
    return namespace


def qub0047_sequence():
    # 7/760 Ns = 0.921%, matching the uploaded cohort QC evidence.
    return 'A' * 753 + 'N' * 7


def qub0047_qc_row():
    # The only field needed by the manual-acceptance predicate is the recorded
    # hard-failure reason. The sequence and key metrics mirror the supplied
    # QUB_0047 example (760 bp; 0.921% Ns; high post-trimming quality).
    return {
        'SequenceID': 'QUB_0047',
        'QCClass': 'FAIL_QC',
        'FailureReasons': 'output_length_lt_800',
        'OutputLength': '760',
        'OutputNPercent': '0.921',
        'MeanQuality': '52.81',
        'OutputExpectedErrors': '0.1094',
        'ConflictDensity': '0',
        'QualityStatus': 'PHRED_VERIFIED',
    }


def test_any_nonempty_failed_sanger_consensus_is_reviewable():
    module = load_paper_trail_module()
    row = qub0047_qc_row()
    options = {'final_min_length': 800}

    # The specific 760 bp QUB_0047 case is eligible.
    result = module._failed_qc_manual_acceptance(qub0047_sequence(), row, **options)
    assert result['eligible'] is True
    assert result['basis'] == 'HUMAN_REVIEW_OVERRIDE'
    assert result['length_deficit_bp'] == 40
    assert 'regardless of which automated QC threshold(s) failed' in result['reason']

    # Human approval remains possible for other hard failures, low quality and
    # excessive Ns, provided a non-empty consensus exists.
    poor_qc = dict(
        row,
        FailureReasons='output_length_lt_800;output_mean_q_lt_20;output_n_percent_gt_5',
        QualityStatus='PHRED_UNVERIFIED',
        MeanQuality='7.0',
        OutputExpectedErrors='99.0',
        ConflictDensity='5.0',
    )
    poor_sequence = 'N' * 100
    result = module._failed_qc_manual_acceptance(poor_sequence, poor_qc, **options)
    assert result['eligible'] is True
    assert result['basis'] == 'HUMAN_REVIEW_OVERRIDE'
    assert 'FAIL_QC' in result['reason']

    no_candidate = module._failed_qc_manual_acceptance('', row, **options)
    assert no_candidate['eligible'] is False
    assert no_candidate['basis'] == 'NO_CANDIDATE_SEQUENCE'


def test_partner_failed_sequences_are_reviewable_regardless_of_length_or_n_content():
    helpers = load_cli_helpers()
    candidate = helpers['_cohort_failed_qc_candidate']
    assert candidate('A' * 100) == (True, 0.0)
    assert candidate('N' * 40 + 'A' * 10) == (True, 80.0)
    assert candidate('') == (False, 100.0)


def _write_tsv(path: Path, fields, row):
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerow(row)


def test_generic_qc_keeps_explicitly_approved_failed_markers(tmp_path):
    if not (ROOT / 'tests' / 'fixtures' / 'qc.py').is_file():
        pytest.skip('tests archive did not include the generic QC fixture required by this test')
    qc = load_qc_module()
    source = tmp_path / 'input.fasta'
    source.write_text(
        '>QUB_0047\n' + 'A' * 760 + '\n'
        '>QUB_0048\n' + 'N' * 90 + 'A' * 810 + '\n'
        '>QUB_0049\n' + 'A' * 700 + '\n'
        '>QUB_0050\n' + 'A' * 900 + '\n'
    )
    output_dir = tmp_path / 'qc'
    qc_fasta = qc.run_qc(
        source, output_dir, min_len=800, max_n_percent=5.0,
        manual_approved_ids={'QUB_0047', 'QUB_0048'},
    )
    retained = {header.split()[0]: sequence for header, sequence in _read_fasta(qc_fasta)}
    assert set(retained) == {'QUB_0047', 'QUB_0048', 'QUB_0050'}
    assert len(retained['QUB_0047']) == 760
    assert len(retained['QUB_0048']) == 900

    rejection_text = (output_dir / 'qc_rejections.tsv').read_text()
    assert 'QUB_0047' in rejection_text and 'INCLUDED_BY_MANUAL_REVIEW' in rejection_text
    assert 'QUB_0049' in rejection_text and 'EXCLUDED' in rejection_text
    override_text = (output_dir / 'qc_manual_overrides.tsv').read_text()
    assert 'QUB_0047' in override_text and 'QUB_0048' in override_text
    assert 'APPROVE_FAIL_QC' in override_text

    # Rejection counters describe records actually excluded, not manually
    # approved records which happened to fail the generic prefilter.
    stats = dict(
        line.split('\t', 1)
        for line in (output_dir / 'qc.stats').read_text().splitlines()
        if '\t' in line
    )
    assert stats['retained_by_manual_override'] == '2'
    assert stats['rejected_total'] == '1'
    assert stats['rejected_too_short'] == '1'
    assert stats['rejected_too_many_n'] == '0'


def test_resume_accepts_any_nonempty_failed_candidate_but_requires_reviewer_notes(tmp_path):
    helpers = load_cli_helpers()
    prep = tmp_path / '00_cohort_preparation'
    prep.mkdir()
    sequence = 'N' * 100  # deliberately short and high-N: reviewer may still approve it.
    (prep / 'cohort_recoverable_failed_markers.fasta').write_text(f'>QUB_FAIL_1\n{sequence}\n')
    _write_tsv(
        prep / 'cohort_marker_review_template.tsv',
        ['sequence_id', 'review_type', 'qc_class', 'eligible_for_manual_acceptance', 'manual_acceptance_basis'],
        {'sequence_id': 'QUB_FAIL_1', 'review_type': 'FAILED_QC_REVIEW', 'qc_class': 'FAIL_QC',
         'eligible_for_manual_acceptance': 'yes', 'manual_acceptance_basis': 'HUMAN_REVIEW_OVERRIDE'},
    )
    _write_tsv(
        prep / 'cohort_marker_qc.tsv',
        ['SequenceID', 'QCClass', 'OutputLength', 'OutputNPercent', 'QualityStatus', 'FailureReasons'],
        {'SequenceID': 'QUB_FAIL_1', 'QCClass': 'FAIL_QC', 'OutputLength': '100',
         'OutputNPercent': '100', 'QualityStatus': 'PHRED_UNVERIFIED',
         'FailureReasons': 'output_length_lt_800;output_n_percent_gt_5'},
    )
    (prep / 'cohort_prepare_state.json').write_text(json.dumps({
        'failed_qc_review_policy': {'manual_approval': 'REVIEWER_OVERRIDE_FOR_ANY_NONEMPTY_CANDIDATE'}
    }))
    decisions = tmp_path / 'decisions.tsv'
    _write_tsv(
        decisions,
        ['sequence_id', 'decision', 'reviewer', 'notes'],
        {'sequence_id': 'QUB_FAIL_1', 'decision': 'APPROVE_FAIL_QC',
         'reviewer': 'reviewer-1', 'notes': 'Reviewed chromatogram and accepted for this marker-specific analysis.'},
    )

    accepted, errors = helpers['_cohort_validate_decisions'](prep, decisions)
    assert not errors
    assert len(accepted) == 1
    assert accepted[0]['decision'] == 'APPROVE_FAIL_QC'

    # The same decision cannot be applied when no candidate sequence is available.
    (prep / 'cohort_recoverable_failed_markers.fasta').write_text('')
    accepted, errors = helpers['_cohort_validate_decisions'](prep, decisions)
    assert not accepted
    assert any('no non-empty candidate' in error for error in errors)

    # Restore candidate; a reviewer note is still mandatory.
    (prep / 'cohort_recoverable_failed_markers.fasta').write_text(f'>QUB_FAIL_1\n{sequence}\n')
    _write_tsv(
        decisions,
        ['sequence_id', 'decision', 'reviewer', 'notes'],
        {'sequence_id': 'QUB_FAIL_1', 'decision': 'APPROVE_FAIL_QC',
         'reviewer': 'reviewer-1', 'notes': ''},
    )
    accepted, errors = helpers['_cohort_validate_decisions'](prep, decisions)
    assert not accepted
    assert any('approval requires reviewer, notes' in error for error in errors)


if __name__ == '__main__':
    test_any_nonempty_failed_sanger_consensus_is_reviewable()
    test_partner_failed_sequences_are_reviewable_regardless_of_length_or_n_content()
    from tempfile import TemporaryDirectory
    with TemporaryDirectory() as temp:
        test_generic_qc_keeps_explicitly_approved_failed_markers(Path(temp))
    with TemporaryDirectory() as temp:
        test_resume_accepts_any_nonempty_failed_candidate_but_requires_reviewer_notes(Path(temp))
    print('Manual failed-QC override tests: PASS')


def test_failed_read_without_consensus_is_available_as_explicit_review_candidate():
    from types import SimpleNamespace

    module = load_paper_trail_module()
    failed_read = SimpleNamespace(
        read_id='JDP771_61777713_61777713',
        oriented_sequence='',
        masked_sequence='A' * 900 + 'N' * 87,
        trimmed_sequence='A' * 900 + 'N' * 87,
        raw_sequence='A' * 1000,
        direction='forward',
        oriented_qualities=[],
        trimmed_qualities=[35] * 987,
    )
    candidate, source, read_id = module._select_failed_qc_review_candidate('', [failed_read])
    assert candidate == failed_read.masked_sequence
    assert source == 'BEST_RECOVERABLE_READ'
    assert read_id == failed_read.read_id

    result = module._failed_qc_manual_acceptance(
        candidate,
        {'Reasons': 'no_output_sequence;failed_no_reads;n_percent_gt_5'},
        final_min_length=800,
    )
    assert result['eligible'] is True
    assert len(candidate) == 987
    assert module._select_failed_qc_review_candidate('', [])[1] == 'NO_CANDIDATE'


def test_studio_uses_candidate_length_when_final_consensus_length_is_zero(tmp_path):
    studio_module = load_studio_module()
    studio = studio_module.SangerStudio(tmp_path, tmp_path / 'studio_output')
    candidate_fasta = tmp_path / 'failed_review_candidates.fasta'
    candidate_fasta.write_text('>QUB_0327\n' + 'A' * 987 + '\n')
    review_path = tmp_path / 'cohort_marker_review_template.tsv'
    _write_tsv(
        review_path,
        [
            'sequence_id', 'decision', 'reviewer', 'reviewed_at', 'notes', 'review_type',
            'eligible_for_manual_acceptance', 'candidate_fasta', 'candidate_source',
            'candidate_read_id', 'candidate_length', 'candidate_n_percent',
            'qc_class', 'output_length', 'n_percent', 'allowed_decisions',
        ],
        {
            'sequence_id': 'QUB_0327', 'decision': '', 'reviewer': '', 'reviewed_at': '',
            'notes': '', 'review_type': 'FAILED_QC_REVIEW',
            'eligible_for_manual_acceptance': 'no', 'candidate_fasta': str(candidate_fasta),
            'candidate_source': 'BEST_RECOVERABLE_READ',
            'candidate_read_id': 'JDP771_61777713_61777713', 'candidate_length': '987',
            'candidate_n_percent': '8.713', 'qc_class': 'FAIL_QC',
            'output_length': '0', 'n_percent': '',
            'allowed_decisions': 'REJECT|RESEQUENCE',
        },
    )
    state = studio.load_marker_review(review_path)
    row = state['rows'][0]
    assert row['allowed_decisions_list'] == [
        'APPROVE_FAIL_QC', 'REJECT', 'RESEQUENCE'
    ]
    saved = studio.save_marker_review({
        'reviewer': 'reviewer-1',
        'rows': [{
            'sequence_id': 'QUB_0327', 'decision': 'APPROVE_FAIL_QC',
            'notes': 'Reviewed the individual read and chromatogram; accepted for this marker analysis.',
        }],
    })
    with open(saved['saved'], newline='') as handle:
        saved_row = next(csv.DictReader(handle, delimiter='\t'))
    assert saved_row['decision'] == 'APPROVE_FAIL_QC'
    assert saved_row['eligible_for_manual_acceptance'] == 'yes'


def _make_sanger_read(module, read_id, *, length, mean_q, direction, base):
    sequence = base * length
    return module.SangerRead(
        read_id=read_id,
        sequence_id='BSI_0016',
        primer='27F' if direction == 'forward' else '1492R',
        direction=direction,
        source_file=f'{read_id}.ab1',
        raw_sequence=sequence,
        raw_qualities=[mean_q] * length,
        quality_available=True,
        quality_source='ab1_pcon',
        trimmed_sequence=sequence,
        trimmed_qualities=[mean_q] * length,
        masked_sequence=sequence,
        masked_qualities=[mean_q] * length,
        oriented_sequence=sequence,
        oriented_qualities=[mean_q] * length,
        status='kept',
        qc_class='PASS_HIGH_CONFIDENCE',
        processing_mode='assemble',
    )


def test_no_overlap_auto_fallback_uses_longest_usable_read_not_primer_direction():
    from unittest.mock import patch

    module = load_paper_trail_module()
    forward = _make_sanger_read(module, 'BSI_0016_27F', length=399, mean_q=30, direction='forward', base='A')
    reverse = _make_sanger_read(module, 'BSI_0016_1492R', length=787, mean_q=43, direction='reverse', base='C')
    with patch.object(module, 'merge_pair', return_value=None):
        sequence, stats = module.assemble_read_group([forward, reverse])

    assert sequence == reverse.oriented_sequence
    assert len(sequence) == 787
    assert stats['status'] == 'partial_no_overlap'
    assert stats['selected_read_id'] == 'BSI_0016_1492R'
    assert stats['selection_source'] == 'AUTO_LONGEST_USABLE_READ'
    assert stats['used_read_ids'] == ['BSI_0016_1492R']
    assert stats['unmerged_reads'] == 'BSI_0016_27F'
    assert stats['fallback_reason'] == 'no_qualifying_overlap'


def test_no_overlap_manual_fallback_selection_overrides_automatic_recommendation():
    from unittest.mock import patch

    module = load_paper_trail_module()
    forward = _make_sanger_read(module, 'BSI_0016_27F', length=399, mean_q=30, direction='forward', base='A')
    reverse = _make_sanger_read(module, 'BSI_0016_1492R', length=787, mean_q=43, direction='reverse', base='C')
    with patch.object(module, 'merge_pair', return_value=None):
        sequence, stats = module.assemble_read_group(
            [forward, reverse], fallback_read_id='BSI_0016_27F',
        )

    assert sequence == forward.oriented_sequence
    assert stats['selected_read_id'] == 'BSI_0016_27F'
    assert stats['selection_source'] == 'MANUAL'
    candidate, source, read_id = module._select_failed_qc_review_candidate(
        sequence, [forward, reverse], assembly_stats=stats,
    )
    assert candidate == forward.oriented_sequence
    assert source == 'BEST_RECOVERABLE_READ'
    assert read_id == 'BSI_0016_27F'


def test_select_best_read_uses_same_longest_first_policy():
    module = load_paper_trail_module()
    forward = _make_sanger_read(module, 'BSI_0019_27F', length=748, mean_q=45, direction='forward', base='A')
    reverse = _make_sanger_read(module, 'BSI_0019_1492R', length=403, mean_q=50, direction='reverse', base='C')
    sequence, stats = module.select_best_read_group([reverse, forward])
    assert sequence == forward.oriented_sequence
    assert stats['selected_read_id'] == 'BSI_0019_27F'
    assert stats['selection_source'] == 'AUTO_LONGEST_USABLE_READ'


def test_no_overlap_fallback_below_800bp_remains_failed_qc():
    from unittest.mock import patch

    module = load_paper_trail_module()
    forward = _make_sanger_read(module, 'BSI_0016_27F', length=399, mean_q=30, direction='forward', base='A')
    reverse = _make_sanger_read(module, 'BSI_0016_1492R', length=787, mean_q=43, direction='reverse', base='C')
    with patch.object(module, 'merge_pair', return_value=None):
        sequence, stats = module.assemble_read_group([forward, reverse])
    qc = module._classify_output_qc(
        sequence, stats['qualities'], stats, [forward, reverse],
        final_min_length=800, min_mean_quality=20, max_output_expected_errors=5,
        warn_n_percent=3, max_n_percent=5, max_conflict_density=1,
        max_ambiguous_overlap_conflicts_without_review=2,
        processing_mode='assemble', mode_warning='',
    )
    assert qc['qc_class'] == 'FAIL_QC'
    assert 'output_length_lt_800' in qc['reasons']
    assert 'partial_no_overlap' in qc['reasons']
