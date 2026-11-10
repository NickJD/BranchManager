from pathlib import Path


def test_generic_input_folder_picker_is_primary_and_keeps_other_routes():
    html = (Path(__file__).resolve().parents[1] / 'src/branchmanager/data/sanger_studio.html').read_text()
    assert html.count('id="chooseInputFolderBtn"') == 1
    assert html.count('id="chooseMapBtn"') == 1  # standalone-map alternative
    assert html.count('id="chooseReadsBtn"') == 1  # direct AB1/ABI alternative
    assert 'id="chooseCohortFolderBtn"' not in html
    assert 'Choose input folder' in html
    assert 'Assistant cohort outputs, standalone sample maps, or raw AB1/ABI chromatograms' in html
    assert 'NOT_SCREENED' in html
    assert 'Assistant cohort mode' in html


def test_shutdown_control_is_present_and_calls_local_shutdown_api():
    html = (Path(__file__).resolve().parents[1] / 'src/branchmanager/data/sanger_studio.html').read_text()
    assert html.count('id="shutdownBtn"') == 1
    assert 'Shut down Studio' in html
    assert "api('/api/shutdown'" in html
    assert 'window.confirm(' in html
    assert 'saved outputs and review decisions remain on disk' in html
