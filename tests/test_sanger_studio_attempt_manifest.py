import csv
from pathlib import Path

from branchmanager.sanger_studio import _read_normalised_read_map


def test_attempt_metadata_survives_normalised_manifest_loading(tmp_path):
    manifest = tmp_path / "manifest.csv"
    fields = [
        "sequence_id", "read_file", "direction", "primer", "attempt_id",
        "attempt_group", "source_path", "source_file", "source_sha256",
    ]
    rows = [
        ["SS172", "SS172_907R_attempt01.ab1", "reverse", "907R", "attempt01",
         "EC00170083", "/source/EC00170083/ss172_907R.ab1", "ss172_907R.ab1", "a" * 64],
        ["SS172", "SS172_907R_attempt02.ab1", "reverse", "907R", "attempt02",
         "EC00072402_ss_1st_round", "/source/EC00072402_ss_1st_round/ss172_907R.ab1",
         "ss172_907R.ab1", "b" * 64],
    ]
    with manifest.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(fields)
        writer.writerows(rows)

    by_read, by_sequence = _read_normalised_read_map(manifest)
    assert by_sequence == {
        "SS172": ["SS172_907R_attempt01", "SS172_907R_attempt02"]
    }
    assert by_read["SS172_907R_attempt01"]["attempt_id"] == "attempt01"
    assert by_read["SS172_907R_attempt02"]["attempt_group"] == "EC00072402_ss_1st_round"
    assert by_read["SS172_907R_attempt02"]["source_sha256"] == "b" * 64



def test_sanger_studio_read_quality_table_displays_attempt_provenance():
    html = (
        Path(__file__).resolve().parents[1]
        / "src" / "branchmanager" / "data" / "sanger_studio.html"
    ).read_text()
    assert "Attempt / submission" in html
    assert "x.AttemptID" in html
    assert "x.AttemptGroup" in html
