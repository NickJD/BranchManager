import csv
from pathlib import Path


def _write_tsv(path: Path, rows):
    fields = list(rows[0])
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def test_selected_assembly_candidate_persists_in_marker_review_tsv(tmp_path):
    from branchmanager.sanger_studio import SangerStudio

    review = tmp_path / 'cohort_marker_review_template.tsv'
    _write_tsv(review, [{
        'sequence_id': 'AU_0003', 'decision': '', 'reviewer': '', 'reviewed_at': '',
        'notes': '', 'review_type': 'WARNING_REVIEW', 'qc_class': 'PASS_WITH_WARNINGS', 'qc_reasons': '',
        'output_length': '', 'n_percent': '', 'allowed_decisions': 'APPROVE|REJECT|RESEQUENCE', 'eligible_for_manual_acceptance': 'no',
        'fallback_read_options': '', 'recommended_read_id': '', 'selected_read_id': '',
        'read_selection_source': '', 'fallback_reason': '', 'processing_mode': 'assemble',
        'assembly_status': 'assembled', 'output_type': '',
    }])
    studio = SangerStudio(tmp_path, tmp_path / 'out')
    studio.load_marker_review(review)
    studio.last_inspection['AU_0003'] = {
        'assembly_candidates': [{
            'CandidateID': 'AU_0003_candidate_01', 'ReadIDs': 'ss175_785F;ss175_907R_attempt01',
            'OutputLength': '1451', 'Sequence': 'ACGT' * 363 + 'A',
        }],
        'fallback_read_options': [], 'recommended_read_id': '',
        'fallback_reason': '', 'processing_mode': 'assemble', 'assembly_status': 'assembled',
    }

    saved = studio.save_marker_review({
        'reviewer': 'Nick',
        'rows': [{
            'sequence_id': 'AU_0003', 'decision': '', 'notes': '',
            'selected_read_id': '', 'selected_candidate_id': 'AU_0003_candidate_01',
        }],
    })
    saved_path = Path(saved['saved'])
    with saved_path.open(newline='', encoding='utf-8') as handle:
        row = next(csv.DictReader(handle, delimiter='\t'))
    assert row['selected_candidate_id'] == 'AU_0003_candidate_01'
    assert row['selected_candidate_type'] == 'ASSEMBLY'
    assert row['selected_candidate_length'] == '1451'
    assert row['selected_candidate_source'] == 'MANUAL'
    assert row['reviewer'] == 'Nick'


def test_candidate_selection_requires_reviewer_identity(tmp_path):
    from branchmanager.sanger_studio import SangerStudio


    review = tmp_path / 'cohort_marker_review_template.tsv'
    _write_tsv(review, [{
        'sequence_id': 'AU_0003', 'decision': '', 'reviewer': '', 'reviewed_at': '',
        'notes': '', 'review_type': 'WARNING_REVIEW', 'qc_class': 'PASS_WITH_WARNINGS', 'qc_reasons': '',
        'output_length': '', 'n_percent': '', 'allowed_decisions': 'APPROVE|REJECT|RESEQUENCE', 'eligible_for_manual_acceptance': 'no',
        'fallback_read_options': '', 'recommended_read_id': '', 'selected_read_id': '',
        'read_selection_source': '', 'fallback_reason': '', 'processing_mode': 'assemble',
        'assembly_status': 'assembled', 'output_type': '',
    }])
    studio = SangerStudio(tmp_path, tmp_path / 'out')
    studio.load_marker_review(review)
    studio.last_inspection['AU_0003'] = {
        'assembly_candidates': [{'CandidateID': 'cand1', 'ReadIDs': 'a;b', 'OutputLength': '900', 'Sequence': 'ACGT'}],
        'fallback_read_options': [], 'recommended_read_id': '', 'fallback_reason': '',
        'processing_mode': 'assemble', 'assembly_status': 'assembled',
    }
    try:
        studio.save_marker_review({'reviewer': '', 'rows': [{'sequence_id': 'AU_0003', 'selected_candidate_id': 'cand1'}]})
    except ValueError as exc:
        assert 'reviewer identity is required' in str(exc)
    else:
        raise AssertionError('Expected selecting a candidate without reviewer identity to be rejected')


def test_review_ui_syncs_inspected_candidate_and_saves_candidate_id():
    html = (Path(__file__).parents[1] / 'src' / 'branchmanager' / 'data' / 'sanger_studio.html').read_text(encoding='utf-8')
    assert "reviewRow.selected_candidate_id=String(candidate.CandidateID||'')" in html
    assert "reviewRow.selected_candidate_type=componentReads.length<=1?'INDIVIDUAL_READ':'ASSEMBLY'" in html
    assert 'selected_candidate_id:cohortReview.rows.find' in html
    assert 'Selected retained candidate:' in html
    assert 'selectedCandidateId?retainedLine:finalLine+candidateLine' in html


def test_unchanged_saved_candidate_is_preserved_without_current_inspection(tmp_path):
    from branchmanager.sanger_studio import SangerStudio

    review = tmp_path / 'cohort_marker_review_template.tsv'
    _write_tsv(review, [{
        'sequence_id': 'AU_0001', 'decision': '', 'reviewer': 'Nick', 'reviewed_at': '',
        'notes': '', 'review_type': 'FAILED_QC_REVIEW', 'qc_class': 'FAIL_QC', 'qc_reasons': '',
        'output_length': '', 'n_percent': '', 'allowed_decisions': 'RESEQUENCE|APPROVE_FAIL_QC', 'eligible_for_manual_acceptance': 'yes',
        'fallback_read_options': '', 'recommended_read_id': '', 'selected_read_id': '',
        'read_selection_source': '', 'fallback_reason': '', 'processing_mode': 'assemble',
        'assembly_status': 'partial_no_overlap', 'output_type': '',
        'selected_candidate_id': 'AU_0001_single_read_05',
        'selected_candidate_type': 'INDIVIDUAL_READ', 'selected_candidate_length': '848',
        'selected_candidate_source': 'MANUAL',
    }])
    studio = SangerStudio(tmp_path, tmp_path / 'out')
    studio.load_marker_review(review)
    saved = studio.save_marker_review({'reviewer': 'Nick', 'rows': [{
        'sequence_id': 'AU_0001', 'decision': '', 'notes': '',
        'selected_read_id': '', 'selected_candidate_id': 'AU_0001_single_read_05',
        'selected_candidate_type': 'INDIVIDUAL_READ', 'selected_candidate_length': '848',
        'selected_candidate_source': 'MANUAL',
    }]})
    with Path(saved['saved']).open(newline='', encoding='utf-8') as handle:
        row = next(csv.DictReader(handle, delimiter='\t'))
    assert row['selected_candidate_id'] == 'AU_0001_single_read_05'
    assert row['selected_candidate_type'] == 'INDIVIDUAL_READ'
    assert row['selected_candidate_length'] == '848'
    assert row['selected_candidate_source'] == 'MANUAL'


def test_changed_candidate_still_requires_current_inspection(tmp_path):
    from branchmanager.sanger_studio import SangerStudio

    review = tmp_path / 'cohort_marker_review_template.tsv'
    _write_tsv(review, [{
        'sequence_id': 'AU_0001', 'decision': '', 'reviewer': 'Nick', 'reviewed_at': '',
        'notes': '', 'review_type': 'FAILED_QC_REVIEW', 'qc_class': 'FAIL_QC', 'qc_reasons': '',
        'output_length': '', 'n_percent': '', 'allowed_decisions': 'RESEQUENCE|APPROVE_FAIL_QC', 'eligible_for_manual_acceptance': 'yes',
        'fallback_read_options': '', 'recommended_read_id': '', 'selected_read_id': '',
        'read_selection_source': '', 'fallback_reason': '', 'processing_mode': 'assemble',
        'assembly_status': 'partial_no_overlap', 'output_type': '',
        'selected_candidate_id': 'old_candidate', 'selected_candidate_type': 'ASSEMBLY',
        'selected_candidate_length': '900', 'selected_candidate_source': 'MANUAL',
    }])
    studio = SangerStudio(tmp_path, tmp_path / 'out')
    studio.load_marker_review(review)
    try:
        studio.save_marker_review({'reviewer': 'Nick', 'rows': [{
            'sequence_id': 'AU_0001', 'decision': '', 'notes': '',
            'selected_candidate_id': 'new_candidate',
        }]})
    except ValueError as exc:
        assert 'newly selected candidate' in str(exc)
    else:
        raise AssertionError('Expected newly selected candidate without inspection to be rejected')


def test_marker_review_distinguishes_retained_candidate_from_inactive_fallback():
    html = (Path(__file__).parents[1] / 'src' / 'branchmanager' / 'data' / 'sanger_studio.html').read_text(encoding='utf-8')
    assert 'RETAINED SEQUENCE' in html
    assert 'Fallback recommendation (inactive)' in html
    assert "if(retained){" in html


def test_newly_selected_fallback_read_is_validated_against_latest_inspection(tmp_path):
    from branchmanager.sanger_studio import SangerStudio


    review = tmp_path / 'cohort_marker_review_template.tsv'
    _write_tsv(review, [{
        'sequence_id': 'AU_0003', 'decision': '', 'reviewer': 'Nick', 'reviewed_at': '',
        'notes': '', 'review_type': 'WARNING_REVIEW', 'qc_class': 'PASS_WITH_WARNINGS', 'qc_reasons': '',
        'output_length': '', 'n_percent': '', 'allowed_decisions': 'APPROVE|REJECT|RESEQUENCE', 'eligible_for_manual_acceptance': 'no',
        'fallback_read_options': '', 'recommended_read_id': '', 'selected_read_id': '',
        'read_selection_source': '', 'fallback_reason': '', 'processing_mode': 'assemble',
        'assembly_status': 'assembled', 'output_type': '',
    }])
    studio = SangerStudio(tmp_path, tmp_path / 'out')
    studio.load_marker_review(review)
    studio.last_inspection['AU_0003'] = {
        'assembly_candidates': [],
        'fallback_read_options': [{'read_id': 'ss175_907R_attempt01', 'length': 834}],
        'recommended_read_id': 'ss175_907R_attempt01',
        'fallback_reason': '', 'processing_mode': 'assemble', 'assembly_status': 'assembled',
    }
    saved = studio.save_marker_review({'reviewer': 'Nick', 'rows': [{
        'sequence_id': 'AU_0003', 'decision': '', 'notes': '',
        'selected_candidate_id': 'ss175_907R_attempt01', 'selected_read_id': '',
    }]})
    with Path(saved['saved']).open(newline='', encoding='utf-8') as handle:
        row = next(csv.DictReader(handle, delimiter='\t'))
    assert row['selected_candidate_id'] == 'ss175_907R_attempt01'
    assert row['selected_candidate_type'] == 'INDIVIDUAL_READ'
    assert row['selected_candidate_length'] == '834'
    assert row['selected_candidate_source'] == 'MANUAL'
