from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
HTML = ROOT / 'src' / 'branchmanager' / 'data' / 'sanger_studio.html'


def test_attach_review_first_prompts_for_chromatograms_when_none_loaded():
    source = HTML.read_text(encoding='utf-8')
    assert 'async function attachReview()' in source
    attach = source[source.index('async function attachReview()'):source.index('function renderCohortReview(')]
    assert 'await chooseReads()' in attach
    assert 'reconcileReviewReads()' in attach


def test_review_read_reconciliation_reports_ids_and_unmatched_files():
    source = HTML.read_text(encoding='utf-8')
    assert 'function reconcileReviewReads()' in source
    for text in ('Exact canonical ID', 'Partner alias match', 'Ambiguous alias', 'Unmatched'):
        assert text in source
    assert 'AB1 inferred ID' in source
    assert 'AB1 filename(s)' in source
    assert 'Unmatched or ambiguous IDs are not silently renamed.' in source


def test_alias_reconciliation_only_uses_aliases_explicitly_present_in_review_rows():
    source = HTML.read_text(encoding='utf-8')
    assert "keys=['partner_own_ids','partner_alias','partner_aliases'" in source
    assert 'canonical.set(sid,sid)' in source
