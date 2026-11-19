# BranchManager tests

This directory contains the automated tests for BranchManager. The suite covers command-line behaviour, Sanger-Studio and AB1 processing, marker QC and assembly, taxonomy, novelty and baseline assessment, pangenome logic, genome evidence imports, reporting, project workflows, and the ground-truth validation toolkit.

## Requirements

- Python 3.12 or newer.
- BranchManager installed from the checked-out source.
- `pytest`.

Some tests exercise optional integrations or external bioinformatics tools. Depending on the test, these may require executables such as VSEARCH, MAFFT or FastTree, and relevant reference data. A test that depends on an unavailable external resource may be skipped or fail; check the test output rather than assuming the full suite is environment-independent.

## Install and run

From the repository root:

```bash
conda activate branchmanager
python -m pip install -e .
python -m pip install pytest
python -m pytest -q
```

Run a specific test module:

```bash
python -m pytest -q tests/test_primers_to_manifest.py
```

Run tests matching a name or keyword:

```bash
python -m pytest -q -k "attempt or manifest"
```

For more detail about failures:

```bash
python -m pytest -ra
```

## Relevant test areas

| Area | Example test modules |
|---|---|
| Sanger-Studio and chromatograms | `test_sanger_studio_*.py`, `test_v1_ab1_review.py` |
| AB1 manifest generation and resequencing | `test_primers_to_manifest.py`, `test_sanger_studio_attempt_manifest.py` |
| Taxonomy and reference handling | `test_taxonomy_screen.py`, `test_taxonomy_and_db.py`, `test_partner_taxonomy_validation.py` |
| Selection, novelty and pangenomes | `test_v1_selection.py`, `test_novelty_identity_ladder.py`, `test_pangenome_domain_fallback.py` |
| Reporting and dashboards | `test_reporting_dashboard.py`, `test_dashboard_cli.py`, `test_demo_dashboard_build.py` |
| Project workflows and imports | `test_operational_workflows.py`, `test_genome_results_import.py`, `test_v1_cohort.py` |
| Ground-truth validation | `test_validation_toolkit.py` |

The examples above are representative, not an exhaustive index. See the test filenames in this directory for the full coverage.

## Manual audit utility

`audit_manual_qc_overrides.py` is a standalone audit utility, not a pytest test module. Run it with its own `--help` option to see the arguments supported by the checked-out version:

```bash
python tests/audit_manual_qc_overrides.py --help
```

## Interpreting results

A passing test suite indicates that the behaviours exercised by those tests passed in the current environment. It does not, by itself, establish the biological correctness of a reference dataset or analysis result.

When reporting test results, include the command used, Python/BranchManager versions, the number of passed, failed and skipped tests, and any unavailable external tools or reference files. Investigate failures individually; do not treat failures caused by missing environment dependencies as equivalent to confirmed software regressions.
