from __future__ import annotations

import csv
import json
import importlib.util
import sys
import pytest
import types
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

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


def read_fasta(path):
    records = []
    header = None
    seq = []
    with open(path, encoding='utf-8') as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line)
    if header is not None:
        records.append((header, ''.join(seq)))
    return records


def write_tsv(path: Path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def load_studio_module():
    # Lightweight package stubs: tests exercise Sanger Studio's file/cache logic,
    # not Biopython, image rendering, VSEARCH, or the project database.
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
    paper_trail.read_fasta = read_fasta
    paper_trail.DEFAULT_PRIMERS = ('27F', '1492R')
    pipeline.paper_trail = paper_trail
    sys.modules['branchmanager.pipeline.paper_trail'] = paper_trail

    run_manifest = types.ModuleType('branchmanager.run_manifest')
    run_manifest.RunManifest = type('RunManifest', (), {})
    sys.modules['branchmanager.run_manifest'] = run_manifest

    name = 'branchmanager.sanger_studio'
    spec = importlib.util.spec_from_file_location(name, PACKAGE_ROOT / 'sanger_studio.py')
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class TestOnboardingReuse:
    def _build_onboarding_run(self, root: Path):
        onboarding = root / '01_onboarding'
        paper = root / '02_paper_trail_merge_meeting'
        screen = paper / 'read_taxonomy_screen'
        onboarding.mkdir(parents=True)
        screen.mkdir(parents=True)
        reads_dir = root / 'raw_reads'
        reads_dir.mkdir()
        read_a = reads_dir / 'JTQ765.ab1'
        read_b = reads_dir / 'JTQ766.ab1'
        read_a.write_bytes(b'ABIF-placeholder')
        read_b.write_bytes(b'ABIF-placeholder')
        write_tsv(onboarding / 'normalised_read_map.tsv',
                  ['sequence_id', 'file', 'primer', 'direction', 'processing_mode', 'partner_id', 'dataset'], [
            {'sequence_id': 'QUB_0289', 'file': str(read_a), 'primer': '27F', 'direction': 'forward', 'processing_mode': 'assemble', 'partner_id': 'QUB', 'dataset': 'QUB_01'},
            {'sequence_id': 'QUB_0289', 'file': str(read_b), 'primer': '1492R', 'direction': 'reverse', 'processing_mode': 'assemble', 'partner_id': 'QUB', 'dataset': 'QUB_01'},
        ])
        tax_fields = ['ID', 'BestHit', 'Identity', 'QueryCoverage', 'TargetCoverage', 'QualifiedHitCount', 'NearBestHitCount', 'AssignedRank', 'Taxon', 'ReadClassificationStatus', 'TaxonomySupport']
        write_tsv(screen / 'taxonomy.tsv', tax_fields, [
            {'ID': 'JTQ765', 'BestHit': 'refA', 'Identity': '99.8', 'QueryCoverage': '99', 'TargetCoverage': '75', 'QualifiedHitCount': '4', 'NearBestHitCount': '3', 'AssignedRank': 'g', 'Taxon': 'd__Bacteria;p__Bacillota;c__Clostridia;o__Order;f__Family;g__Sarcina', 'ReadClassificationStatus': 'RESOLVED', 'TaxonomySupport': 'MULTIPLE_NEAR_BEST_HITS_AGREE_AT_GENUS'},
            {'ID': 'JTQ766', 'BestHit': 'refB', 'Identity': '99.7', 'QueryCoverage': '98', 'TargetCoverage': '74', 'QualifiedHitCount': '3', 'NearBestHitCount': '2', 'AssignedRank': 'g', 'Taxon': 'd__Bacteria;p__Bacillota;c__Clostridia;o__Order;f__Family;g__Sarcina', 'ReadClassificationStatus': 'RESOLVED', 'TaxonomySupport': 'MULTIPLE_NEAR_BEST_HITS_AGREE_AT_GENUS'},
        ])
        write_tsv(screen / 'read_taxonomy_hits.tsv',
                  ['ReadID', 'TargetID', 'Identity', 'QueryCoverage', 'TargetCoverage', 'AlignmentLength', 'Taxonomy', 'Qualified', 'NearBest', 'IdentityDelta'], [
            {'ReadID': 'JTQ765', 'TargetID': 'refA', 'Identity': '99.8', 'QueryCoverage': '99', 'TargetCoverage': '75', 'AlignmentLength': '1400', 'Taxonomy': 'g__Sarcina', 'Qualified': 'yes', 'NearBest': 'yes', 'IdentityDelta': '0.0'},
            {'ReadID': 'JTQ766', 'TargetID': 'refB', 'Identity': '99.7', 'QueryCoverage': '98', 'TargetCoverage': '74', 'AlignmentLength': '1380', 'Taxonomy': 'g__Sarcina', 'Qualified': 'yes', 'NearBest': 'yes', 'IdentityDelta': '0.0'},
        ])
        write_tsv(screen / 'taxonomy_screen_policy.tsv',
                  ['PolicyVersion', 'VSEARCHVersion', 'ReferenceFASTA', 'TaxonomyFile', 'NearBestIdentityDeltaPct'], [
            {'PolicyVersion': '2.1', 'VSEARCHVersion': 'vsearch 2.29', 'ReferenceFASTA': '/db/gtdb.fna.gz', 'TaxonomyFile': 'reference FASTA headers', 'NearBestIdentityDeltaPct': '0.5'},
        ])
        (paper / 'trimmed_oriented_reads.fasta').write_text('>JTQ765\nACGTACGT\n>JTQ766\nTTTTCCCC\n', encoding='utf-8')
        # Deliberately no read_taxonomy_concordance.tsv: the loader must not require it.
        return root, paper, read_a, read_b

    def test_loads_onboarding_assignments_hits_and_concordance_without_concordance_tsv(self):
        studio_module = load_studio_module()
        with TemporaryDirectory() as directory:
            root, paper, read_a, read_b = self._build_onboarding_run(Path(directory))
            studio = studio_module.SangerStudio(root, root / 'studio_out')
            loaded = studio.load_onboarding_outputs(root)
            assert loaded['map_loaded'] is True
            assert loaded['taxonomy_read_count'] == 2
            assert loaded['taxonomy_hit_count'] == 2
            assert loaded['mapped_read_count'] == 2
            assert loaded['isolate_count'] == 1
            assert studio.map_path == root / '01_onboarding' / 'normalised_read_map.tsv'
            assert studio.pre_screen_taxonomy_by_isolate['QUB_0289']['TaxonomyConcordance'] == 'CONCORDANT'

            inspection = {
                'reads': [{'ReadID': 'JTQ765'}, {'ReadID': 'JTQ766'}],
                'taxonomy_reads': [], 'taxonomy_hits': [],
            }
            current_trimmed = paper / 'current_trimmed.fasta'
            current_trimmed.write_text('>JTQ765\nACGTACGT\n>JTQ766\nTTTTCCCC\n', encoding='utf-8')
            enriched = studio._attach_prescreened_taxonomy(inspection, 'QUB_0289', {'trimmed_fasta': str(current_trimmed)})
            assert enriched['taxonomy_screen_requested'] is True
            assert enriched['taxonomy_freshness'] == 'MATCHES_CURRENT_TRIMMED_READS'
            assert enriched['taxonomy_concordance'] == 'CONCORDANT'
            assert len(enriched['taxonomy_reads']) == 2
            assert len(enriched['taxonomy_hits']) == 2

    def test_onboarding_map_without_processing_mode_column_gets_safe_defaults(self):
        studio_module = load_studio_module()
        with TemporaryDirectory() as directory:
            root, _paper, read_a, read_b = self._build_onboarding_run(Path(directory))
            map_path = root / '01_onboarding' / 'normalised_read_map.tsv'
            # processing_mode is not a required column in the Assistant map.
            write_tsv(map_path,
                      ['sequence_id', 'file', 'primer', 'direction', 'partner_id', 'dataset'], [
                {'sequence_id': 'QUB_0289', 'file': str(read_a), 'primer': '27F', 'direction': 'forward', 'partner_id': 'QUB', 'dataset': 'QUB_01'},
                {'sequence_id': 'QUB_0289', 'file': str(read_b), 'primer': '1492R', 'direction': 'reverse', 'partner_id': 'QUB', 'dataset': 'QUB_01'},
            ])
            studio = studio_module.SangerStudio(root, root / 'studio_out')
            studio.load_onboarding_output(root)

            # inspect() persists edits via updated['processing_mode']; every map
            # row therefore needs a value even when the input TSV omits the column.
            assert [row['processing_mode'] for row in studio.map_rows] == ['assemble', 'assemble']
            updated = {row['file']: row for row in studio.map_rows}
            for stored in studio.map_rows:
                stored.update({key: updated[stored['file']][key]
                               for key in ('primer', 'direction', 'processing_mode')})

    def test_multi_partner_assistant_preparation_pairs_sibling_onboarding_and_paper_trail_folders(self):
        studio_module = load_studio_module()
        with TemporaryDirectory() as directory:
            root = Path(directory) / '00_cohort_preparation'
            screen_fields = [
                'ID', 'BestHit', 'Identity', 'QueryCoverage', 'TargetCoverage',
                'QualifiedHitCount', 'NearBestHitCount', 'AssignedRank', 'Taxon',
                'ReadClassificationStatus', 'TaxonomySupport',
            ]
            hit_fields = [
                'ReadID', 'TargetID', 'Identity', 'QueryCoverage', 'TargetCoverage',
                'AlignmentLength', 'Taxonomy', 'Qualified', 'NearBest', 'IdentityDelta',
            ]
            for dataset in ('QUB_01', 'UOG_02'):
                onboarding = root / f'01_{dataset}_onboarding'
                paper = root / f'01_{dataset}_paper_trail'
                screen = paper / 'read_taxonomy_screen'
                onboarding.mkdir(parents=True)
                screen.mkdir(parents=True)
                read_id = f'{dataset}_READ001'
                read_path = root / 'raw_reads' / f'{read_id}.ab1'
                read_path.parent.mkdir(parents=True, exist_ok=True)
                read_path.write_bytes(b'ABIF-placeholder')
                write_tsv(onboarding / 'normalised_read_map.tsv',
                          ['sequence_id', 'file', 'primer', 'direction', 'processing_mode', 'partner_id', 'dataset'], [
                    {'sequence_id': dataset, 'file': str(read_path), 'primer': '27F', 'direction': 'forward',
                     'processing_mode': 'best_read', 'partner_id': dataset.split('_')[0], 'dataset': dataset},
                ])
                write_tsv(screen / 'taxonomy.tsv', screen_fields, [{
                    'ID': read_id, 'BestHit': 'ref1', 'Identity': '99.9', 'QueryCoverage': '99',
                    'TargetCoverage': '80', 'QualifiedHitCount': '2', 'NearBestHitCount': '2',
                    'AssignedRank': 'g', 'Taxon': 'd__Bacteria;p__Bacillota;f__Family;g__Genus',
                    'ReadClassificationStatus': 'RESOLVED', 'TaxonomySupport': 'MULTIPLE_NEAR_BEST_HITS_AGREE_AT_GENUS',
                }])
                write_tsv(screen / 'read_taxonomy_hits.tsv', hit_fields, [{
                    'ReadID': read_id, 'TargetID': 'ref1', 'Identity': '99.9', 'QueryCoverage': '99',
                    'TargetCoverage': '80', 'AlignmentLength': '900', 'Taxonomy': 'g__Genus',
                    'Qualified': 'yes', 'NearBest': 'yes', 'IdentityDelta': '0.0',
                }])
                (paper / 'trimmed_oriented_reads.fasta').write_text(f'>{read_id}\nACGTACGT\n', encoding='utf-8')

            studio = studio_module.SangerStudio(root, root / 'studio_out')
            loaded = studio.load_onboarding_output(root)
            assert loaded['isolate_count'] == 2
            assert loaded['read_count'] == 2
            assert loaded['taxonomy_read_count'] == 2
            assert set(studio.onboarding_taxonomy['map_paths']) == {
                str((root / '01_QUB_01_onboarding' / 'normalised_read_map.tsv').resolve()),
                str((root / '01_UOG_02_onboarding' / 'normalised_read_map.tsv').resolve()),
            }
            assert studio.pre_screen_taxonomy_by_isolate['QUB_01']['TaxonomyConcordance'] == 'INSUFFICIENT_CLASSIFICATION'
            assert studio.pre_screen_taxonomy_by_isolate['UOG_02']['TaxonomyConcordance'] == 'INSUFFICIENT_CLASSIFICATION'

    def test_selecting_partner_onboarding_folder_discovers_matching_paper_trail_sibling(self):
        studio_module = load_studio_module()
        with TemporaryDirectory() as directory:
            root = Path(directory) / '00_cohort_preparation'
            onboarding = root / '01_QUB_01_onboarding'
            paper = root / '01_QUB_01_paper_trail'
            screen = paper / 'read_taxonomy_screen'
            reads_dir = root / 'raw_reads'
            screen.mkdir(parents=True)
            onboarding.mkdir(parents=True)
            reads_dir.mkdir()
            read_path = reads_dir / 'QUB_READ001.ab1'
            read_path.write_bytes(b'ABIF-placeholder')
            write_tsv(onboarding / 'normalised_read_map.tsv',
                      ['sequence_id', 'file', 'primer', 'direction', 'partner_id', 'dataset'], [{
                          'sequence_id': 'QUB_0001', 'file': str(read_path), 'primer': '27F',
                          'direction': 'forward', 'partner_id': 'QUB', 'dataset': 'QUB_01',
                      }])
            write_tsv(screen / 'taxonomy.tsv',
                      ['ID', 'BestHit', 'Identity', 'QueryCoverage', 'TargetCoverage', 'QualifiedHitCount',
                       'NearBestHitCount', 'AssignedRank', 'Taxon', 'ReadClassificationStatus', 'TaxonomySupport'], [{
                          'ID': 'QUB_READ001', 'BestHit': 'ref1', 'Identity': '99.1', 'QueryCoverage': '99',
                          'TargetCoverage': '80', 'QualifiedHitCount': '1', 'NearBestHitCount': '1',
                          'AssignedRank': 'g', 'Taxon': 'd__Bacteria;p__Bacillota;g__Genus',
                          'ReadClassificationStatus': 'RESOLVED', 'TaxonomySupport': 'SINGLE_NEAR_BEST_REFERENCE_HIT',
                      }])
            write_tsv(screen / 'read_taxonomy_hits.tsv',
                      ['ReadID', 'TargetID', 'Identity', 'QueryCoverage', 'TargetCoverage', 'AlignmentLength',
                       'Taxonomy', 'Qualified', 'NearBest', 'IdentityDelta'], [{
                          'ReadID': 'QUB_READ001', 'TargetID': 'ref1', 'Identity': '99.1', 'QueryCoverage': '99',
                          'TargetCoverage': '80', 'AlignmentLength': '900', 'Taxonomy': 'g__Genus',
                          'Qualified': 'yes', 'NearBest': 'yes', 'IdentityDelta': '0',
                      }])
            (paper / 'trimmed_oriented_reads.fasta').write_text('>QUB_READ001\nACGTACGT\n', encoding='utf-8')

            studio = studio_module.SangerStudio(root, root / 'studio_out')
            loaded = studio.load_onboarding_output(onboarding)
            assert loaded['loaded'] is True
            assert loaded['map_loaded'] is True
            assert loaded['isolate_count'] == 1
            assert loaded['taxonomy_read_count'] == 1
            assert loaded['map_paths'] == [str((onboarding / 'normalised_read_map.tsv').resolve())]
            assert studio.map_rows[0]['sequence_id'] == 'QUB_0001'

    def test_changed_qc_trimmed_sequence_marks_cached_taxonomy_stale(self):
        studio_module = load_studio_module()
        with TemporaryDirectory() as directory:
            root, paper, _read_a, _read_b = self._build_onboarding_run(Path(directory))
            studio = studio_module.SangerStudio(root, root / 'studio_out')
            studio.load_onboarding_outputs(root)
            current_trimmed = paper / 'current_trimmed.fasta'
            current_trimmed.write_text('>JTQ765\nACGTACGA\n>JTQ766\nTTTTCCCC\n', encoding='utf-8')
            result = studio._attach_prescreened_taxonomy(
                {'reads': [{'ReadID': 'JTQ765'}, {'ReadID': 'JTQ766'}]},
                'QUB_0289', {'trimmed_fasta': str(current_trimmed)},
            )
            assert result['taxonomy_freshness'] == 'STALE_TRIMMED_SEQUENCE_CHANGED'


class TestReviewPersistence:
    def test_legacy_saved_review_without_allowed_decisions_restores_fail_qc_approval(self):
        studio_module = load_studio_module()
        with TemporaryDirectory() as directory:
            root = Path(directory)
            candidate = root / 'failed_markers.fasta'
            candidate.write_text('>QUB_0047|qc=FAIL_QC\nACGTACGTACGT\n', encoding='utf-8')
            review = root / 'cohort_marker_review_decisions.tsv'
            write_tsv(review,
                      ['sequence_id', 'decision', 'reviewer', 'reviewed_at', 'notes', 'review_type', 'eligible_for_manual_acceptance', 'candidate_fasta', 'qc_class', 'qc_reasons', 'output_length', 'n_percent'], [
                {'sequence_id': 'QUB_0047', 'decision': '', 'reviewer': '', 'reviewed_at': '', 'notes': '', 'review_type': 'FAILED_QC_REVIEW', 'eligible_for_manual_acceptance': 'no', 'candidate_fasta': str(candidate), 'qc_class': 'FAIL_QC', 'qc_reasons': 'output_length_lt_800', 'output_length': '760', 'n_percent': '0.9'},
            ])
            studio = studio_module.SangerStudio(root, root / 'studio_out')
            state = studio.load_marker_review(review)
            row = state['rows'][0]
            assert row['allowed_decisions_list'] == ['APPROVE_FAIL_QC', 'REJECT', 'RESEQUENCE']
            assert row['eligible_for_manual_acceptance'] == 'yes'
            saved = studio.save_marker_review({
                'reviewer': 'Nick',
                'rows': [{'sequence_id': 'QUB_0047', 'decision': 'APPROVE_FAIL_QC', 'notes': 'Manual review accepted.'}],
            })
            saved_path = Path(saved['saved'])
            with saved_path.open(newline='', encoding='utf-8') as handle:
                saved_row = next(csv.DictReader(handle, delimiter='\t'))
            assert saved_row['allowed_decisions'] == 'APPROVE_FAIL_QC|REJECT|RESEQUENCE'
            assert saved_row['eligible_for_manual_acceptance'] == 'yes'

            reopened = studio_module.SangerStudio(root, root / 'studio_out_2')
            reloaded = reopened.load_marker_review(saved_path)
            assert reloaded['rows'][0]['allowed_decisions_list'] == ['APPROVE_FAIL_QC', 'REJECT', 'RESEQUENCE']
            assert reloaded['rows'][0]['decision'] == 'APPROVE_FAIL_QC'

    def test_relative_candidate_fasta_is_resolved_against_saved_review_file(self):
        studio_module = load_studio_module()
        with TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'failed_markers.fasta').write_text('>QUB_0047\nACGTACGT\n', encoding='utf-8')
            review = root / 'cohort_marker_review_decisions.tsv'
            write_tsv(review,
                      ['sequence_id', 'decision', 'reviewer', 'reviewed_at', 'notes', 'review_type',
                       'eligible_for_manual_acceptance', 'candidate_fasta', 'qc_class'], [
                {'sequence_id': 'QUB_0047', 'decision': '', 'reviewer': '', 'reviewed_at': '', 'notes': '',
                 'review_type': 'FAILED_QC_REVIEW', 'eligible_for_manual_acceptance': 'no',
                 'candidate_fasta': 'failed_markers.fasta', 'qc_class': 'FAIL_QC'},
            ])
            studio = studio_module.SangerStudio(root, root / 'studio_out')
            row = studio.load_marker_review(review)['rows'][0]
            assert row['eligible_for_manual_acceptance'] == 'yes'
            assert 'APPROVE_FAIL_QC' in row['allowed_decisions_list']
            assert Path(row['candidate_fasta']) == (root / 'failed_markers.fasta').resolve()

    def test_manual_fallback_read_selection_is_saved_and_reloaded(self):
        studio_module = load_studio_module()
        with TemporaryDirectory() as directory:
            root = Path(directory)
            review = root / 'cohort_marker_review_template.tsv'
            options = [
                {'read_id': 'BSI_0016_27F', 'primer': '27F', 'direction': 'forward',
                 'length': 787, 'mean_quality': 43.0, 'expected_errors': 1.1, 'n_percent': 2.8,
                 'qc_class': 'PASS_HIGH_CONFIDENCE'},
                {'read_id': 'BSI_0016_1492R', 'primer': '1492R', 'direction': 'reverse',
                 'length': 403, 'mean_quality': 50.0, 'expected_errors': 0.05, 'n_percent': 0.25,
                 'qc_class': 'PASS_HIGH_CONFIDENCE'},
            ]
            write_tsv(review, [
                'sequence_id', 'decision', 'reviewer', 'reviewed_at', 'notes', 'review_type',
                'eligible_for_manual_acceptance', 'allowed_decisions', 'qc_class', 'qc_reasons',
                'output_length', 'candidate_length', 'fallback_read_options',
                'recommended_read_id', 'selected_read_id', 'read_selection_source', 'fallback_reason',
            ], [{
                'sequence_id': 'BSI_0016', 'decision': '', 'reviewer': '', 'reviewed_at': '', 'notes': '',
                'review_type': 'FAILED_QC_REVIEW', 'eligible_for_manual_acceptance': 'yes',
                'allowed_decisions': 'APPROVE_FAIL_QC|REJECT|RESEQUENCE', 'qc_class': 'FAIL_QC',
                'qc_reasons': 'partial_no_overlap;output_length_lt_800', 'output_length': '787',
                'candidate_length': '787', 'fallback_read_options': json.dumps(options),
                'recommended_read_id': 'BSI_0016_27F', 'selected_read_id': 'BSI_0016_27F',
                'read_selection_source': 'AUTO_LONGEST_USABLE_READ', 'fallback_reason': 'no_qualifying_overlap',
            }])
            studio = studio_module.SangerStudio(root, root / 'studio_out')
            loaded = studio.load_marker_review(review)
            assert loaded['rows'][0]['selected_read_id'] == 'BSI_0016_27F'
            saved = studio.save_marker_review({
                'reviewer': 'Nick',
                'rows': [{'sequence_id': 'BSI_0016', 'decision': '', 'notes': '',
                          'selected_read_id': 'BSI_0016_1492R'}],
            })
            saved_path = Path(saved['saved'])
            with saved_path.open(newline='', encoding='utf-8') as handle:
                row = next(csv.DictReader(handle, delimiter='\t'))
            assert row['selected_read_id'] == 'BSI_0016_1492R'
            assert row['read_selection_source'] == 'MANUAL'
            assert row['reviewer'] == 'Nick'

            reopened = studio_module.SangerStudio(root, root / 'studio_out_2')
            restored = reopened.load_marker_review(saved_path)['rows'][0]
            assert restored['selected_read_id'] == 'BSI_0016_1492R'
            assert restored['read_selection_source'] == 'MANUAL'

    def test_empty_failed_marker_without_candidate_does_not_offer_approval(self):
        studio_module = load_studio_module()
        with TemporaryDirectory() as directory:
            root = Path(directory)
            candidate = root / 'failed_markers.fasta'
            candidate.write_text('>OTHER\nACGT\n', encoding='utf-8')
            review = root / 'review.tsv'
            write_tsv(review, ['sequence_id', 'decision', 'reviewer', 'reviewed_at', 'notes', 'review_type', 'eligible_for_manual_acceptance', 'candidate_fasta', 'qc_class', 'output_length'], [
                {'sequence_id': 'QUB_EMPTY', 'decision': '', 'reviewer': '', 'reviewed_at': '', 'notes': '', 'review_type': 'FAILED_QC_REVIEW', 'eligible_for_manual_acceptance': 'no', 'candidate_fasta': str(candidate), 'qc_class': 'FAIL_QC', 'output_length': '0'},
            ])
            studio = studio_module.SangerStudio(root, root / 'studio_out')
            state = studio.load_marker_review(review)
            assert state['rows'][0]['allowed_decisions_list'] == ['REJECT', 'RESEQUENCE']


if __name__ == '__main__':
    import pytest
    raise SystemExit(pytest.main([__file__, '-q']))


class TestGenericInputFolderDetection:
    def test_detects_assistant_cohort_and_loads_paired_taxonomy_once(self):
        studio_module = load_studio_module()
        with TemporaryDirectory() as directory:
            root = Path(directory)
            builder = TestOnboardingReuse()
            builder._build_onboarding_run(root)
            studio = studio_module.SangerStudio(root, root / 'studio_out')
            result = studio.load_input_directory(root)
            assert result['input_kind'] == 'assistant_cohort'
            assert result['taxonomy_loaded'] is True
            assert result['onboarding_taxonomy']['loaded'] is True
            assert result['map_loaded'] is True
            assert len(result['samples']) == 1
            assert result['samples'][0]['sequence_id'] == 'QUB_0289'
            assert studio.input_source == 'assistant_cohort'

    def test_detects_normalised_read_map_without_taxonomy_as_not_screened(self):
        studio_module = load_studio_module()
        with TemporaryDirectory() as directory:
            root = Path(directory)
            onboarding = root / 'Partner_01_onboarding'
            onboarding.mkdir()
            read_path = root / 'raw_reads' / 'ISOLATE_A_27F.ab1'
            read_path.parent.mkdir()
            read_path.write_bytes(b'ABIF-placeholder')
            write_tsv(onboarding / 'normalised_read_map.tsv',
                      ['sequence_id', 'file', 'primer', 'direction', 'processing_mode'], [{
                          'sequence_id': 'ISOLATE_A', 'file': str(read_path), 'primer': '27F',
                          'direction': 'forward', 'processing_mode': 'best_read',
                      }])
            studio = studio_module.SangerStudio(root, root / 'studio_out')
            result = studio.load_input_directory(root)
            assert result['input_kind'] == 'normalised_read_map'
            assert result['taxonomy_loaded'] is False
            assert result['taxonomy_status'] == 'NOT_SCREENED'
            assert result['onboarding_taxonomy']['loaded'] is False
            assert result['samples'][0]['sequence_id'] == 'ISOLATE_A'
            assert studio.onboarding_taxonomy is None

    def test_detects_standalone_sample_map_without_assistant_cohort(self):
        studio_module = load_studio_module()
        with TemporaryDirectory() as directory:
            root = Path(directory)
            reads = root / 'reads'
            reads.mkdir()
            read_a = reads / 'SAMPLE_A_27F.ab1'
            read_b = reads / 'SAMPLE_A_1492R.ab1'
            read_a.write_bytes(b'ABIF-placeholder')
            read_b.write_bytes(b'ABIF-placeholder')
            write_tsv(root / 'sample_map.tsv',
                      ['sequence_id', 'file', 'primer', 'direction'], [
                {'sequence_id': 'SAMPLE_A', 'file': str(read_a), 'primer': '27F', 'direction': 'forward'},
                {'sequence_id': 'SAMPLE_A', 'file': str(read_b), 'primer': '1492R', 'direction': 'reverse'},
            ])
            studio = studio_module.SangerStudio(root, root / 'studio_out')
            result = studio.load_input_directory(root)
            assert result['input_kind'] == 'standalone_sample_map'
            assert result['taxonomy_loaded'] is False
            assert result['taxonomy_status'] == 'NOT_SCREENED'
            assert result['samples'][0]['sequence_id'] == 'SAMPLE_A'
            assert result['samples'][0]['read_count'] == 2
            assert studio.onboarding_taxonomy is None

    def test_detects_raw_ab1_directory_without_cohort_or_sample_map(self):
        studio_module = load_studio_module()
        with TemporaryDirectory() as directory:
            root = Path(directory)
            reads = root / 'chromatograms'
            reads.mkdir()
            (reads / 'SAMPLE_B_27F.ab1').write_bytes(b'ABIF-placeholder')
            (reads / 'SAMPLE_B_1492R.ab1').write_bytes(b'ABIF-placeholder')
            studio_module.paper_trail.infer_read_info = lambda name, primers: (
                'SAMPLE_B', '27F' if '27F' in name else '1492R',
                'forward' if '27F' in name else 'reverse',
            )
            studio = studio_module.SangerStudio(root, root / 'studio_out')
            result = studio.load_input_directory(reads)
            assert result['input_kind'] == 'raw_chromatograms'
            assert result['taxonomy_status'] == 'NOT_SCREENED'
            assert len(result['samples']) == 1
            assert result['samples'][0]['sequence_id'] == 'SAMPLE_B'
            assert result['samples'][0]['read_count'] == 2

    def test_standalone_load_clears_previously_loaded_assistant_taxonomy(self):
        studio_module = load_studio_module()
        with TemporaryDirectory() as directory:
            root = Path(directory)
            builder = TestOnboardingReuse()
            builder._build_onboarding_run(root)
            studio = studio_module.SangerStudio(root, root / 'studio_out')
            studio.load_onboarding_output(root)
            assert studio.onboarding_taxonomy is not None

            standalone = root / 'standalone'
            standalone.mkdir()
            read_file = root / 'raw_reads' / 'OTHER_27F.ab1'
            read_file.write_bytes(b'ABIF-placeholder')
            write_tsv(standalone / 'sample_map.tsv',
                      ['sequence_id', 'file', 'primer', 'direction'], [{
                          'sequence_id': 'OTHER', 'file': str(read_file),
                          'primer': '27F', 'direction': 'forward',
                      }])
            result = studio.load_input_directory(standalone)
            assert result['input_kind'] == 'standalone_sample_map'
            assert studio.onboarding_taxonomy is None
            assert studio.pre_screen_taxonomy_by_isolate == {}

    def test_unrecognised_folder_has_actionable_error(self):
        studio_module = load_studio_module()
        with TemporaryDirectory() as directory:
            root = Path(directory)
            studio = studio_module.SangerStudio(root, root / 'studio_out')
            with pytest.raises(ValueError, match='No recognised Sanger inputs'):
                studio.load_input_directory(root)


    def test_discovers_and_pairs_multiple_partner_onboarding_and_paper_trail_directories(self):
        studio_module = load_studio_module()
        with TemporaryDirectory() as directory:
            root = Path(directory)
            taxonomy_fields = [
                'ID', 'BestHit', 'Identity', 'QueryCoverage', 'TargetCoverage',
                'QualifiedHitCount', 'NearBestHitCount', 'AssignedRank', 'Taxon',
                'ReadClassificationStatus', 'TaxonomySupport',
            ]
            partner_data = {
                'QUB': ('QUB_0001', 'QUBREAD1', 'g__Ruminococcus'),
                'BSI': ('BSI_0001', 'BSIREAD1', 'g__Prevotella'),
            }
            for partner, (sequence_id, read_id, genus) in partner_data.items():
                dataset = f'{partner}_01'
                onboarding = root / f'01_{dataset}_onboarding'
                screen = root / f'01_{dataset}_paper_trail' / 'read_taxonomy_screen'
                onboarding.mkdir(parents=True)
                screen.mkdir(parents=True)
                write_tsv(onboarding / 'normalised_read_map.tsv',
                          ['sequence_id', 'file', 'primer', 'direction', 'processing_mode', 'partner_id', 'dataset'], [{
                              'sequence_id': sequence_id,
                              'file': f'/nonexistent/{read_id}.ab1',
                              'primer': '27F', 'direction': 'forward',
                              'processing_mode': 'best_read', 'partner_id': partner,
                              'dataset': dataset,
                          }])
                write_tsv(screen / 'taxonomy.tsv', taxonomy_fields, [{
                    'ID': read_id, 'BestHit': f'ref_{partner}', 'Identity': '99.0',
                    'QueryCoverage': '100', 'TargetCoverage': '70',
                    'QualifiedHitCount': '1', 'NearBestHitCount': '1',
                    'AssignedRank': 'g',
                    'Taxon': f'd__Bacteria;p__Bacillota;c__Clostridia;o__Order;f__Family;{genus}',
                    'ReadClassificationStatus': 'RESOLVED',
                    'TaxonomySupport': 'SINGLE_NEAR_BEST_REFERENCE_HIT',
                }])
            studio = studio_module.SangerStudio(root, root / 'studio_out')
            result = studio.load_input_directory(root)
            assert result['input_kind'] == 'assistant_cohort'
            assert result['taxonomy_loaded'] is True
            assert result['partner_count'] == 2
            assert result['mapped_read_count'] == 2
            assert result['taxonomy_read_count'] == 2
            assert {row['sequence_id'] for row in result['samples']} == {'QUB_0001', 'BSI_0001'}
            assert len(result['onboarding_taxonomy']['map_paths']) == 2


def test_webpage_shutdown_endpoint_stops_only_the_local_studio_server(monkeypatch, tmp_path):
    import importlib.resources
    import json as json_module
    import threading
    import time
    import urllib.request

    studio_module = load_studio_module()
    real_server = studio_module.ThreadingHTTPServer
    captured = {}

    class CapturingServer(real_server):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            captured['server'] = self

    monkeypatch.setattr(studio_module, 'ThreadingHTTPServer', CapturingServer)
    monkeypatch.setattr(importlib.resources, 'files', lambda _package: PACKAGE_ROOT)
    monkeypatch.chdir(tmp_path)
    server_thread = threading.Thread(
        target=studio_module.serve, kwargs={'port': 0, 'open_browser': False}, daemon=True,
    )
    server_thread.start()
    deadline = time.time() + 5
    while 'server' not in captured and time.time() < deadline:
        time.sleep(0.02)
    assert 'server' in captured, 'Sanger Studio server did not start'
    server = captured['server']
    url = f'http://127.0.0.1:{server.server_port}/api/shutdown'
    request = urllib.request.Request(
        url, data=b'{}', headers={'Content-Type': 'application/json'}, method='POST',
    )
    with urllib.request.urlopen(request, timeout=5) as response:
        result = json_module.loads(response.read().decode('utf-8'))
    assert result['stopping'] is True
    assert 'remain on disk' in result['message']
    server_thread.join(timeout=5)
    assert not server_thread.is_alive(), 'Studio server did not stop after shutdown request'

