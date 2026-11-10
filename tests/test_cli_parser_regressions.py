from branchmanager.cli.parser import build_parser


def test_parser_builds_without_duplicate_option_conflicts():
    parser = build_parser()
    assert parser is not None


def test_import_genome_results_parser_has_staging_output_argument():
    parser = build_parser()
    args = parser.parse_args([
        "import-genome-results",
        "--db", "project.db",
        "--results-tsv", "genomes.tsv",
        "--out", "genome_import",
        "--dry-run",
    ])
    assert args.db == "project.db"
    assert args.results_tsv == "genomes.tsv"
    assert args.out == "genome_import"
    assert args.dry_run is True


def test_import_genome_evidence_parser_builds_with_expected_arguments():
    parser = build_parser()
    args = parser.parse_args([
        "import-genome-evidence",
        "--db", "project.db",
        "--workbook", "partners.xlsm",
        "--genome-dir", "genomes",
        "--gtdbtk", "gtdbtk.tsv",
        "--checkm2", "checkm2.tsv",
        "--dataset", "AAFC_01",
        "--out", "genome_evidence",
    ])
    assert args.out == "genome_evidence"
    assert args.gtdbtk_version == "not_recorded"
    assert args.checkm2_version == "not_recorded"


import pytest
from branchmanager import __version__
from branchmanager.cli.parser import build_parser


@pytest.mark.parametrize("option", ["--version", "-v"])
def test_cli_version_option(option, capsys):
    parser = build_parser()

    with pytest.raises(SystemExit) as exc:
        parser.parse_args([option])

    assert exc.value.code == 0

    captured = capsys.readouterr()
    assert f"branchmanager {__version__}" in captured.out