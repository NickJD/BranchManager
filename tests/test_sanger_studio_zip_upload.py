import io
import zipfile

import pytest

from branchmanager.sanger_studio import SangerStudio


def _zip_bytes(entries):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, content in entries:
            archive.writestr(name, content)
    return buffer.getvalue()


def test_upload_zip_loads_nested_ab1_and_abi_files(tmp_path):
    studio = SangerStudio(tmp_path, tmp_path / "runs")
    payload = _zip_bytes([
        ("batch/QUB_0001_16S_27F.ab1", b"fake ab1"),
        ("batch/nested/QUB_0001_16S_1492R.abi", b"fake abi"),
        ("README.txt", b"ignored"),
    ])

    result = studio.load_uploaded_ab1_zip(payload, "reads.zip")

    assert result["input_kind"] == "raw_chromatograms"
    assert result["uploaded_archive"] == "reads.zip"
    assert result["read_count"] == 2
    assert len(result["samples"]) == 1
    assert result["samples"][0]["sequence_id"] == "QUB_0001_16S"
    assert all((tmp_path / "runs" / "input_uploads").rglob("*.ab1"))


def test_upload_zip_rejects_non_zip_filename(tmp_path):
    studio = SangerStudio(tmp_path, tmp_path / "runs")
    with pytest.raises(ValueError, match=r"\.zip"):
        studio.load_uploaded_ab1_zip(b"not a zip", "reads.tar")


def test_upload_zip_rejects_archive_without_chromatograms(tmp_path):
    studio = SangerStudio(tmp_path, tmp_path / "runs")
    payload = _zip_bytes([("notes.txt", b"nothing to process")])
    with pytest.raises(ValueError, match="No .ab1 or .abi"):
        studio.load_uploaded_ab1_zip(payload, "reads.zip")


def test_upload_zip_rejects_path_traversal(tmp_path):
    studio = SangerStudio(tmp_path, tmp_path / "runs")
    payload = _zip_bytes([("../outside.ab1", b"bad")])
    with pytest.raises(ValueError, match="Unsafe path"):
        studio.load_uploaded_ab1_zip(payload, "reads.zip")


def test_sanger_studio_page_exposes_zip_upload():
    from importlib.resources import files

    page = files("branchmanager").joinpath("data/sanger_studio.html").read_text(encoding="utf-8")
    assert 'id="uploadZipBtn"' in page
    assert 'id="ab1ZipInput"' in page
    assert "/api/upload-ab1-zip" in page


def test_upload_zip_rejects_windows_style_path_traversal(tmp_path):
    studio = SangerStudio(tmp_path, tmp_path / "runs")
    payload = _zip_bytes([("..\\\\outside.ab1", b"bad")])
    with pytest.raises(ValueError, match="Unsafe path"):
        studio.load_uploaded_ab1_zip(payload, "reads.zip")
