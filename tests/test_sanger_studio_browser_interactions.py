from pathlib import Path

HTML = Path(__file__).parents[1] / "src" / "branchmanager" / "data" / "sanger_studio.html"


def test_candidate_visibility_is_enforced_in_dom_and_incompatible_buttons_disabled():
    html = HTML.read_text(encoding="utf-8")
    block = html[html.index("function applyPrimerReadChoices"):html.index("choiceBody?.addEventListener", html.index("function applyPrimerReadChoices"))]
    assert "card.hidden=!compatible" in block
    assert "card.setAttribute('aria-hidden',compatible?'false':'true')" in block
    assert ".candidate-card[hidden]{display:none!important}" in html
    assert "card.hidden=!compatible" in block  # compatibility is computed from selected reads and distinct primers


def test_candidate_cards_show_component_read_count_and_individual_read_chips():
    html = HTML.read_text(encoding="utf-8")
    assert "candidate-component-count" in html
    assert "component read(s)" in html
    assert "candidate-read-chip" in html
    assert "String(x.ReadIDs||'').split(';').map(id=>id.trim()).filter(Boolean).map(id=>`<span class=\"candidate-read-chip\">${esc(id)}</span>`).join('')" in html


def test_trace_selector_uses_candidate_reads_mapped_to_consensus_position():
    html = HTML.read_text(encoding="utf-8")
    start = html.index("function renderTraceReview")
    end = html.index("function showAlignmentPosition", start)
    block = html[start:end]
    assert "const rows=alignmentReads.map(read=>" in block
    assert "const rp=trimmedPositionAt(read,position)" in block
    assert "return rp>0?{read,rp,trace}:null" in block
    assert "const available=!!(trace&&trace.available&&(trace.channels||{}).A)" in block
    assert "available?'':'disabled'" in block


def test_trace_selector_checkbox_is_disabled_when_read_has_no_mapped_base_or_trace():
    html = HTML.read_text(encoding="utf-8")
    start = html.index("function renderTraceReview")
    end = html.index("function showAlignmentPosition", start)
    block = html[start:end]
    assert "const rows=alignmentReads.map(read=>" in block
    assert "const usable=rows.filter(x=>x.trace&&x.trace.available&&(x.trace.channels||{}).A)" in block
    assert "const available=!!(trace&&trace.available&&(trace.channels||{}).A)" in block
    assert "${available?'':'disabled'}" in block
