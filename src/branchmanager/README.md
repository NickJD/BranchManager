# BranchManager source package

This directory contains the installable Python package. It is not intended as a user-facing workflow guide; start with the repository [README](../../README.md).

## Main components

| Path | Responsibility |
|---|---|
| `cli/` | Argument parser and CLI package entry points. |
| `cli_core.py` | Command orchestration and high-level workflow implementations. |
| `data/` | Packaged Sanger Studio interface and bundled reference-anchor resources. |
| `db/` | Database interface and schema. |
| `pipeline/` | Sequence QC, classification, novelty, phylogenetic and selection logic. |
| `tools/` | Standalone helper commands. |
| `utils/` | Shared sequence, text and file utilities. |
| `sanger_studio.py` | Local Sanger Studio server, API and upload handling. |
| `mailroom.py` | Sanger delivery inventory and sample-map preparation. |
| `onboarding.py` | Partner metadata and isolate validation. |
| `taxonomy.py` / `taxonomy_io.py` | Taxonomy handling and input/output. |
| `reporting.py` | Report and dashboard generation, including nearest-hit evidence from novelty metrics. |
| `run_manifest.py` | Run provenance and manifest support. |

## Dashboard regeneration and evidence fields

The `branchmanager dashboard` CLI command is a reporting-only wrapper around
`reporting.write_performance_review_dashboard()`. It resolves an existing run
directory (or its `assessment/` directory), reads the existing report files and
writes `performance_review_dashboard.html` into the run root. It must not invoke
the biological assessment pipeline or modify the source reports.

The dashboard's nearest-hit evidence is sourced preferentially from
`novelty_metrics.tsv`. The comparison prefixes currently used by the reporting
layer are `Hungate`, `SecondaryBaseline`, `CulturedProject`, `Project` and
`Reference`. Keep the identity, hit identifier, query coverage and alignment
length associated with the same comparison prefix. Legacy selection-summary
fields are fallbacks, not a replacement for the comparison-specific novelty
metrics.

When changing this reporting behaviour, add regression tests covering:
- precedence of novelty-metrics nearest-hit fields over stale summary aliases;
- each cultured and partner-project comparison pool;
- absent or incomplete novelty metrics;
- marker-QC warnings on project identity;
- the generated dashboard payload and visible detail panel.

## Development conventions

- Keep the package version in `branchmanager.__version__` synchronised with `pyproject.toml`.
- Keep CLI version reporting sourced from `__version__`.
- Preserve provenance and distinguish observed evidence from recommendations.
- Add regression tests for changes to scientific interpretation or selection eligibility.
- Avoid silent fallback between taxonomy evidence sources.
- Keep packaged data declarations in `pyproject.toml` and `MANIFEST.in` synchronised.

See [Development and testing](../../docs/development.md).
