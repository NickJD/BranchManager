from pathlib import Path


HTML = Path(__file__).parents[1] / "src" / "branchmanager" / "data" / "sanger_studio.html"


def test_incompatible_candidate_cards_are_hidden_and_not_focusable():
    source = HTML.read_text(encoding="utf-8")
    assert ".candidate-card[hidden]{display:none!important}" in source
    assert "card.hidden=!compatible" in source
    assert "aria-hidden" in source


def test_chromatogram_hover_moves_a_visible_svg_cursor_and_reports_channel_support():
    source = HTML.read_text(encoding="utf-8")
    assert "function attachTraceHover" in source
    assert "pointermove" in source
    assert "mousemove" in source
    assert "data-trace-hover-cursor" in source
    assert "dominant channel" in source
    assert "pointerleave" in source


def test_clearing_read_selection_resets_and_disables_alignment_controls():
    source = HTML.read_text(encoding="utf-8")
    assert "emptyRange.max='1'" in source
    assert "emptyRange.disabled=true" in source
    assert "emptyStart.disabled=true" in source
    assert "alignmentWindowLabel" in source
    assert "No active candidate" in source


def test_trace_hover_cursor_and_readout_share_plot_margin_mapping():
    source = HTML.read_text(encoding='utf-8')
    block = source[source.index('function attachTraceHover'):source.index('function renderTraceReview')]
    assert 'const plotX=leftMargin+plotFraction*Math.max(1,vbWidth-leftMargin-rightMargin)' in block
    assert "cursor.setAttribute('x1',String(plotX))" in block
    assert 'raw samples ${rawSampleStart}–${rawSampleEnd}; max-pooled display' in block
