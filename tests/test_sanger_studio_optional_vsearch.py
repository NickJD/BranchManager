from __future__ import annotations

import csv
import importlib.util
import sys
import types
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

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

PACKAGE_ROOT = ROOT / 'src' / 'branchmanager'


def write_tsv(path: Path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def load_studio_module():
    branchmanager = types.ModuleType('branchmanager')
    branchmanager.__path__ = [str(PACKAGE_ROOT)]
    sys.modules['branchmanager'] = branchmanager
    pipeline = types.ModuleType('branchmanager.pipeline')
    pipeline.__path__ = [str(PACKAGE_ROOT / 'pipeline')]
    sys.modules['branchmanager.pipeline'] = pipeline

    paper_trail = types.ModuleType('branchmanager.pipeline.paper_trail')
    paper_trail.DEFAULT_QC_POLICY = SimpleNamespace(
        min_quality=20, min_read_length=300, min_final_length=800,
        min_mean_quality=20.0, max_output_expected_errors=5.0,
        warn_n_percent=3.0, max_n_percent=5.0,
        min_overlap=80, min_overlap_identity=0.98,
    )
    paper_trail.DEFAULT_PRIMERS = ('27F', '1492R')
    pipeline.paper_trail = paper_trail
    sys.modules['branchmanager.pipeline.paper_trail'] = paper_trail

    run_manifest = types.ModuleType('branchmanager.run_manifest')

    class FakeManifest:
        last_instance = None
        def __init__(self, *_args):
            self.data = {}
            self.stages = []
            FakeManifest.last_instance = self
        def add_input(self, *_args, **_kwargs): pass
        def add_output(self, *_args, **_kwargs): pass
        def add_stage(self, name, status, **kwargs): self.stages.append((name, status, kwargs))
        def write(self): pass
        def finish(self, *_args, **_kwargs): pass

    run_manifest.RunManifest = FakeManifest
    sys.modules['branchmanager.run_manifest'] = run_manifest

    name = 'branchmanager.sanger_studio'
    spec = importlib.util.spec_from_file_location(name, PACKAGE_ROOT / 'sanger_studio.py')
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module, paper_trail


def _fake_paper_trail_run(calls):
    def run_paper_trail(_paths, outdir, **kwargs):
        calls.append(kwargs)
        out = Path(outdir)
        assembly = out / 'assembly.tsv'
        reads = out / 'read_qc.tsv'
        placements = out / 'assembly_read_placements.tsv'
        empty_files = {
            'per_base_error_tsv': out / 'per_base_error.tsv',
            'recommendations_tsv': out / 'recommendations.tsv',
            'marker_review_template_tsv': out / 'marker_review_template.tsv',
            'failed_review_candidates_fasta': out / 'failed_review_candidates.fasta',
            'qc_policy_tsv': out / 'qc_policy.tsv',
            'summary': out / 'summary.txt',
        }
        write_tsv(assembly, ['SequenceID', 'Status', 'QCClass', 'Recommendation', 'Reasons', 'OutputLength', 'UsedReads', 'ReadCount'], [
            {'SequenceID': 'QUB_1001', 'Status': 'single_read', 'QCClass': 'PASS_HIGH_CONFIDENCE', 'Recommendation': 'KEEP', 'Reasons': '', 'OutputLength': '950', 'UsedReads': '1', 'ReadCount': '1'},
        ])
        write_tsv(reads, ['ReadID', 'SequenceID', 'Primer', 'Direction', 'RawLength', 'TrimmedLength', 'MaskedNPercent', 'QCClass', 'Status'], [
            {'ReadID': 'READ001', 'SequenceID': 'QUB_1001', 'Primer': '27F', 'Direction': 'forward', 'RawLength': '1000', 'TrimmedLength': '950', 'MaskedNPercent': '0.1', 'QCClass': 'PASS_HIGH_CONFIDENCE', 'Status': 'kept'},
        ])
        write_tsv(placements, ['SequenceID', 'ReadID'], [{'SequenceID': 'QUB_1001', 'ReadID': 'READ001'}])
        for path in empty_files.values():
            path.write_text('placeholder\n', encoding='utf-8')
        trimmed = out / 'trimmed_oriented_reads.fasta'
        trimmed.write_text('>READ001\nACGTACGT\n', encoding='utf-8')
        taxonomy_screen_tsv = None
        if kwargs.get('screen_ref'):
            screen = out / 'read_taxonomy_screen'
            screen.mkdir(parents=True, exist_ok=True)
            taxonomy = screen / 'taxonomy.tsv'
            hits = screen / 'read_taxonomy_hits.tsv'
            write_tsv(taxonomy, ['ID', 'Taxon', 'ReadClassificationStatus'], [
                {'ID': 'READ001', 'Taxon': 'g__Example', 'ReadClassificationStatus': 'RESOLVED'},
            ])
            write_tsv(hits, ['ReadID', 'TargetID', 'Qualified'], [
                {'ReadID': 'READ001', 'TargetID': 'refA', 'Qualified': 'yes'},
                {'ReadID': 'READ001', 'TargetID': 'refB', 'Qualified': 'yes'},
            ])
            taxonomy_screen_tsv = str(screen / 'read_taxonomy_concordance.tsv')
            Path(taxonomy_screen_tsv).write_text('SequenceID\tConcordanceStatus\nQUB_1001\tINSUFFICIENT_CLASSIFICATION\n', encoding='utf-8')
        return {
            'assembly_tsv': str(assembly), 'read_qc_tsv': str(reads),
            'assembly_placements_tsv': str(placements),
            'assembled_fasta': str(out / 'assembled.fasta'),
            'per_base_error_tsv': str(empty_files['per_base_error_tsv']),
            'recommendations_tsv': str(empty_files['recommendations_tsv']),
            'marker_review_template_tsv': str(empty_files['marker_review_template_tsv']),
            'failed_review_candidates_fasta': str(empty_files['failed_review_candidates_fasta']),
            'qc_policy_tsv': str(empty_files['qc_policy_tsv']),
            'summary': str(empty_files['summary']),
            'read_error_pngs': [], 'chromatogram_pngs': [], 'assembly_pngs': [],
            'trimmed_fasta': str(trimmed),
            'taxonomy_screen_tsv': taxonomy_screen_tsv,
        }
    return run_paper_trail


def test_inspection_does_not_rerun_vsearch_by_default():
    module, paper_trail = load_studio_module()
    calls = []
    paper_trail.run_paper_trail = _fake_paper_trail_run(calls)
    with TemporaryDirectory() as directory:
        root = Path(directory)
        read = root / 'READ001.ab1'
        read.write_bytes(b'ABIF placeholder')
        studio = module.SangerStudio(root, root / 'studio_out')
        studio.map_rows = [{
            'sequence_id': 'QUB_1001', 'file': str(read), 'primer': '27F',
            'direction': 'forward', 'processing_mode': 'best_read',
        }]
        studio.inspect({'sequence_id': 'QUB_1001', 'settings': {}})
    assert len(calls) == 1
    assert calls[0].get('screen_ref') is None
    assert calls[0].get('screen_taxa') is None


def test_inspection_runs_vsearch_only_when_explicitly_requested():
    module, paper_trail = load_studio_module()
    calls = []
    paper_trail.run_paper_trail = _fake_paper_trail_run(calls)
    with TemporaryDirectory() as directory:
        root = Path(directory)
        read = root / 'READ001.ab1'
        read.write_bytes(b'ABIF placeholder')
        studio = module.SangerStudio(root, root / 'studio_out')
        studio.screen_ref = '/reference/gtdb.fna.gz'
        studio.screen_taxa = '/reference/bac120_taxonomy.tsv.gz'
        studio.map_rows = [{
            'sequence_id': 'QUB_1001', 'file': str(read), 'primer': '27F',
            'direction': 'forward', 'processing_mode': 'best_read',
        }]
        result = studio.inspect({'sequence_id': 'QUB_1001', 'settings': {}, 'rerun_taxonomy': True})
        manifest = module.RunManifest.last_instance
    assert len(calls) == 1
    assert calls[0].get('screen_ref') == '/reference/gtdb.fna.gz'
    assert calls[0].get('screen_taxa') == '/reference/bac120_taxonomy.tsv.gz'
    assert result['taxonomy_freshness'] == 'FRESH_VSEARCH_RERUN'
    assert result['taxonomy_screen_message'] == (
        'Fresh VSEARCH completed against gtdb.fna.gz: 1 per-read assignment(s), 2 candidate hit(s).'
    )
    assert ('taxonomy_screen', 'RUNNING',) == manifest.stages[1][:2]
    assert ('taxonomy_screen', 'COMPLETE',) == manifest.stages[2][:2]


def test_inspection_applies_saved_manual_fallback_read_selection():
    module, paper_trail = load_studio_module()
    calls = []
    paper_trail.run_paper_trail = _fake_paper_trail_run(calls)
    with TemporaryDirectory() as directory:
        root = Path(directory)
        read_a = root / 'READ_FORWARD.ab1'
        read_b = root / 'READ_REVERSE.ab1'
        read_a.write_bytes(b'ABIF placeholder')
        read_b.write_bytes(b'ABIF placeholder')
        studio = module.SangerStudio(root, root / 'studio_out')
        studio.map_rows = [
            {'sequence_id': 'QUB_1001', 'file': str(read_a), 'primer': '27F',
             'direction': 'forward', 'processing_mode': 'assemble'},
            {'sequence_id': 'QUB_1001', 'file': str(read_b), 'primer': '1492R',
             'direction': 'reverse', 'processing_mode': 'assemble'},
        ]
        studio.marker_review_rows = [{
            'sequence_id': 'QUB_1001', 'selected_read_id': 'READ_REVERSE',
            'read_selection_source': 'MANUAL',
        }]
        studio.inspect({'sequence_id': 'QUB_1001', 'settings': {}})
    assert calls[0].get('selected_read_ids') == {'QUB_1001': 'READ_REVERSE'}


def test_frontend_shows_explicit_vsearch_running_and_completion_feedback():
    html = (ROOT / 'src' / 'branchmanager' / 'data' / 'sanger_studio.html').read_text()
    inspect = html[html.index('async function inspect()'):html.index('async function processAll()')]
    message_definition = inspect.index('const runningMessage=rerun')
    status = inspect.index('setRunFeedback(runningMessage,{loading:true})')
    request = inspect.index("await api('/api/inspect',body)")
    completion = inspect.index('r.taxonomy_screen_message')
    assert message_definition < status < request < completion
    assert 'FRESH VSEARCH taxonomy screen' in inspect
    assert 'may take several minutes' in inspect


def test_selecting_another_isolate_resets_processing_mode_to_its_map_setting():
    html = (ROOT / 'src' / 'branchmanager' / 'data' / 'sanger_studio.html').read_text()
    select_sample = html[
        html.index('function selectSample(id)'):
        html.index('\nfunction readSettings()', html.index('function selectSample(id)'))
    ]
    inspect = html[
        html.index('async function inspect()'):
        html.index('\nasync function processAll()', html.index('async function inspect()'))
    ]
    assert "$('#mode').value='map'" in select_sample
    assert 'if(inspectionRunning)' in select_sample
    assert "const requestedMode=$('#mode').value||'map'" in inspect
    assert "processing_mode:requestedMode" in inspect
    assert 'inspectionRunning=true' in inspect
    assert 'inspectionRunning=false' in inspect


def test_inspection_button_and_live_status_explicitly_announce_vsearch():
    html = (ROOT / 'src' / 'branchmanager' / 'data' / 'sanger_studio.html').read_text()
    inspect = html[
        html.index('async function inspect()'):
        html.index('\nasync function processAll()', html.index('async function inspect()'))
    ]
    assert "Inspect + VSEARCH running" in inspect
    assert 'RUNNING: fresh VSEARCH' in inspect
    assert "r.taxonomy_screen_message||'Fresh VSEARCH taxonomy screening completed" in inspect
    assert 'aria-live="polite"' in html
