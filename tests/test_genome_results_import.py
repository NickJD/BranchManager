import pytest

from branchmanager.genome_evidence_import import _parse_qc_pass


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("1", True),
        ("true", True),
        ("passed", True),
        ("0", False),
        ("false", False),
        ("failed", False),
        ("", None),
        ("NA", None),
        ("not assessed", None),
    ],
)
def test_parse_qc_pass(value, expected):
    assert _parse_qc_pass(value, row_number=2) is expected


def test_parse_qc_pass_rejects_unrecognised_values():
    with pytest.raises(ValueError, match="row 7"):
        _parse_qc_pass("probably", row_number=7)
