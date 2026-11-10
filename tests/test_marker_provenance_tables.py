import pytest

from branchmanager.marker_provenance import _read_table, _review_decisions


def test_short_tsv_uses_extension_without_sniffing(tmp_path):
    path = tmp_path / "review.tsv"
    path.write_text(
        "sequence_id\tdecision\treviewer\tnotes\n"
        "ISO_1\tAPPROVE\tReviewer A\tConsensus checked\n",
        encoding="utf-8",
    )
    rows = _read_table(path)
    assert rows == [{
        "sequence_id": "ISO_1",
        "decision": "APPROVE",
        "reviewer": "Reviewer A",
        "notes": "Consensus checked",
    }]
    assert _review_decisions(path)["ISO_1"]["decision"] == "APPROVE"


def test_csv_uses_extension(tmp_path):
    path = tmp_path / "review.csv"
    path.write_text(
        "sequence_id,decision,reviewer,notes\n"
        "ISO_1,REJECT,Reviewer A,No consensus\n",
        encoding="utf-8",
    )
    assert _review_decisions(path)["ISO_1"]["decision"] == "REJECT"


def test_unknown_extension_sniffer_failure_has_safe_fallback(tmp_path):
    path = tmp_path / "review.data"
    path.write_text("sequence_id\tdecision\nISO_1\tRESEQUENCE\n", encoding="utf-8")
    assert _review_decisions(path)["ISO_1"]["decision"] == "RESEQUENCE"


@pytest.mark.parametrize(
    "contents, message",
    [
        ("", "empty or contains no data rows"),
        ("foo\tbar\nvalue\tvalue\n", "no sequence identifier column"),
        ("sequence_id\treviewer\nISO_1\tReviewer A\n", "no decision column"),
    ],
)
def test_invalid_review_file_fails_with_actionable_error(tmp_path, contents, message):
    path = tmp_path / "review.tsv"
    path.write_text(contents, encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        _review_decisions(path)
