from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_sanger_studio_uses_browser_file_selection_cross_platform():
    """Studio must not launch native GUI pickers from HTTP handlers."""
    implementations = [
        ROOT / "src" / "branchmanager" / "sanger_studio.py",
        ROOT / "src" / "branchmanager" / "pipeline" / "sanger_studio.py",
    ]
    for path in implementations:
        source = path.read_text(encoding="utf-8")
        assert "ThreadingHTTPServer" in source  # safe now that request handlers never create GUI windows
        assert "osascript" not in source
        assert "tk.Tk()" not in source
        assert "browser-native file selection" in source
        assert "/api/browser-upload" in source


def test_sanger_studio_html_has_cross_platform_browser_picker_inputs():
    html = (ROOT / "src" / "branchmanager" / "data" / "sanger_studio.html").read_text(encoding="utf-8")
    for picker in [
        "browserReadsInput", "browserFolderInput", "browserMapInput",
        "browserReviewInput", "browserTaxRefInput", "browserTaxTableInput",
    ]:
        assert picker in html
    assert "webkitdirectory" in html
    assert "FormData" in html
    assert "browserReadsInput').click()" in html


def test_browser_folder_upload_preserves_nested_paths_and_blocks_traversal(tmp_path):
    from branchmanager.sanger_studio import SangerStudio

    studio = SangerStudio(tmp_path, tmp_path / "runs")
    saved = studio.save_browser_uploads("folder", [
        ("partnerA/reads/sample_27F.ab1", b"chromatogram"),
        ("../../escape.txt", b"safe"),
    ])
    upload_root = (tmp_path / "runs" / "browser_uploads" / "folder").resolve()
    assert (upload_root / "partnerA" / "reads" / "sample_27F.ab1").read_bytes() == b"chromatogram"
    assert all(Path(p).resolve().is_relative_to(upload_root) for p in saved)
    assert not (tmp_path / "escape.txt").exists()
