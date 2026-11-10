from pathlib import Path

HTML = Path(__file__).parents[1] / "src/branchmanager/data/sanger_studio.html"

def test_cohort_review_header_is_sticky_inside_scroll_container():
    html = HTML.read_text(encoding="utf-8")
    assert ".reviewtablewrap{max-height:68vh;overflow:auto;position:relative;" in html
    assert ".reviewtable thead th{position:sticky;top:0;z-index:4;" in html
    assert ".reviewtable{min-width:1660px;border-collapse:separate;border-spacing:0}" in html

def test_manual_failed_qc_override_is_explained_without_relabelling_qc():
    html = HTML.read_text(encoding="utf-8")
    assert "APPROVE_FAIL_QC</b> is a separate human override" in html
    assert "automated QC/read-level status remains FAIL_QC" in html
    assert "Manual override recorded · automated QC remains FAIL_QC" in html
