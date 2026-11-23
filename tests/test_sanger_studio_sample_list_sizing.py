from pathlib import Path

def test_sanger_studio_isolate_list_fills_sidebar_and_keeps_inspection_button_visible():
    html = Path(__file__).parents[1] / "src" / "branchmanager" / "data" / "sanger_studio.html"
    text = html.read_text()
    assert ".samples{flex:1 1 auto;min-height:0;max-height:none;overflow-y:auto;overflow-x:hidden;" in text
    assert ".sidebar-runbar{flex:0 0 auto;" in text
    assert ".sidebar>.cardbody{flex:1 1 auto;min-height:0;display:flex;flex-direction:column;overflow:hidden}" in text
    assert 'id="runBtnSidebar"' in text
    assert "$('#runBtnSidebar').onclick=inspect" in text
