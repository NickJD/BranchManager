# Tests

The test suite covers CLI behaviour, database handling, sequence QC, taxonomy, genome evidence, pangenome eligibility, reporting and Sanger Studio.

Run the complete suite from the repository root:

```bash
python -m pytest
```

Run one test module:

```bash
python -m pytest tests/test_selection_taxonomy_gating.py -q
```

## Test principles

Tests should protect both software behaviour and biological interpretation. In particular:

- Partner domain/pathway routing must remain separate from GTDB pangenome taxonomy.
- GTDB-unresolved records must not silently fall back to marker-classifier taxonomy for pangenome grouping.
- Supported genus assignments must not be forced to species.
- Hierarchy and sequencing recommendations must agree on eligibility.
- Identity, query coverage and reference coverage must remain separate.
- Manual review decisions and QC failures must not be overwritten by automatic processing.

Some integration tests require external executables or reference resources. Record these prerequisites clearly and distinguish passed, failed, skipped and unexecuted tests in release reports.
