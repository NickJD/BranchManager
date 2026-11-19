# BranchManager v1.2.0 test audit (2026-11-19)

## Changes made

- `pyproject.toml`: added pytest configuration so the source-layout package is importable from the repository root (`pythonpath = ["src"]`, `testpaths = ["tests"]`).
- `tests/test_genomic_16s_selection.py`: corrected the two fake VSEARCH executables. They were incorrectly trying to read/write the `--notmatched` argument, which the production command does not pass. This was a mock/test-harness bug, not a VSEARCH or production workflow issue.

## Environment

- Python 3.13.5
- pytest 9.0.2
- BranchManager 1.2.0 editable install, installed with `--no-build-isolation --no-deps`
- External bioinformatics tools are unavailable in this environment (VSEARCH, Barrnap, MAFFT, FastTree, CD-HIT, DIAMOND, MMseqs2, BLAST+)
- No project reference databases are available.

## Results

Initial full run before the fake-VSEARCH correction: **442 passed, 59 failed, 2 skipped, 44 subtests passed**.

After correcting the mock: `python -m pytest -q tests/test_genomic_16s_selection.py` → **4 passed**.

A subsequent full-suite run was attempted but did not complete within the execution window; it progressed to approximately 68% and was still running in `tests/test_taxonomy_and_db.py::OutputHelperTests::test_build_sequence_assessment_rows_reads_project_novelty_columns`. Therefore, no complete post-change full-suite count is claimed.

The initial failures included legacy tests that expect older APIs/output schemas, taxonomy/pangenome expectation mismatches, unavailable external executables/reference files, Sanger-Studio browser tests asserting old exact HTML/JS strings, and a cohort IT Desk validation failure. These must be individually triaged before a release verdict.

**Do not delete, skip, or weaken tests solely to obtain a green suite.** Update expected behaviour only when the changed scientific or workflow specification is confirmed and documented.

## Full-suite command in your provisioned environment

From the repository root, after activating the environment and installing BranchManager and external tools/reference data:

```bash
python -m pytest -q --junitxml=test-results-v1.2.0.xml
```

Verbose diagnostic run:

```bash
python -m pytest -ra --tb=long
```

Check tool availability:

```bash
for tool in vsearch barrnap mafft FastTree fasttree cd-hit diamond mmseqs blastn blastp; do
    printf '%-12s' "$tool"
    command -v "$tool" || echo MISSING
 done
```

Before release, record versions, a complete pass/fail/skip summary, and the disposition of every remaining failure.

## Follow-up audit and test-contract corrections (19 November 2026)

### Confirmed newer behaviour retained

- Pangenome grouping uses usable GTDB lineage only; partner catalogue values are domain/pathway expectations, not pangenome taxonomic authority.
- Sanger-Studio resequencing, per-primer selection, chromatogram channel mapping, attempt provenance, and candidate taxonomy on the selected sequence are retained.
- The minimum Sanger assembly overlap default is now **40 bp** in `PaperTrailQCPolicy`; the Paper Trail and Assistant CLI defaults inherit this value.

### Test corrections made

- Replaced stale partner-catalogue-authoritative pangenome assertions with tests for GTDB-authoritative grouping, unresolved/low-similarity GTDB rejection, tree conflict checks, and GTDB hierarchy provenance.
- Updated partner metadata tests for the current `selected_for_sequencing` / `selected_for_wgs` schema (the old `selection_provided` key is no longer emitted).
- Updated Mailroom CLI test arguments to require the current manifest interface.
- Corrected the rootstock atomic-publication test mock target to patch `branchmanager.cli_core`, where the implementation now resolves the Filing Cabinet helper.
- Replaced obsolete Sanger-Studio browser implementation-string tests with checks for current read-selection, trace mapping and hidden-card behaviour; added component-read counts and individual read chips to candidate cards, and CSS enforcing hidden incompatible cards.
- Made the isolated taxonomy-screen test fixture self-contained by stubbing its required utility modules; updated its evidence expectations to the current multi-reference support policy and VSEARCH `maxrejects=256` setting.
- Added tests asserting the 40 bp minimum-overlap default in the QC policy and both Paper Trail and Assistant CLI defaults.

### Verification performed

- Targeted updated/relevant test groups, including the new overlap boundary and isolated taxonomy-screen tests: **110 passed, 37 subtests passed** in 1.63 seconds.
- This targeted selection includes the updated Sanger-Studio browser regression groups, taxonomy-contract tests, partner-taxonomy validation, taxonomy-screen tests, and operational-workflow tests.
- The cohort-resume integration test still requires external tools and fixture completion; it is not claimed as passing.
- A full-suite run did not complete within the available execution window. Remaining failure clusters must be triaged; this audit does **not** certify v1.2.0 for release.
