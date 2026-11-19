# BranchManager v1.2.0 Test Repair Iteration — 19 November 2026

## Scope and baseline

Applied a test-maintenance pass to the 19 Nov 12:21 uploaded source archive. This pass is intentionally conservative: it does not change production decision logic or taxonomy/pangenome science.

## Changes made

- Updated stale CLI test patch targets from the former `branchmanager.cli` namespace to the current `branchmanager.cli_core` implementation boundary.
- Fixed the Interview delegation test so it mocks the actual `cmd_paper_trail` function instead of invoking the real QC workflow.
- Updated background-check mocks to the current classification-hit tuple and VSEARCH user-field layout.
- Updated the taxonomy merge expectation to `NA` where confidence is not source-provided; identity is not treated as a substitute confidence score.
- Removed an assertion for a deprecated partner-metadata key; the test still verifies current selected-for-sequencing and genome-availability states.
- Corrected the nearest-search test fixture to create the query FASTA, write to the per-round VSEARCH `--userout` file, assert the current bounded VSEARCH search parameters, and respect the current 80% partial-coverage tier rather than expecting a 58.82%-coverage hit to be accepted.
- Removed invalid assertions in the FASTA-only Paper Trail test that expected partner metadata without supplying any partner metadata; retained the consensus sequence assertion.
- Added dependency-aware skips for tree-building/cohort integration tests when external tools are absent.
- Corrected a test-only undefined-variable/reference issue.

## Verification

- Sanger and taxonomy-gating focused suite: **59 passed** (including 40 bp overlap policy tests).
- Taxonomy/database plus cohort subset after changes: **125 passed, 12 failed, 4 skipped** (2 subtests passed), 8.59 seconds.
- Full suite run was attempted but exceeded the 120-second execution limit before completing; no full-suite result is claimed.
- `python -m compileall -q src tests` completed successfully.

## Remaining failures

The remaining 12 failures are not resolved in this iteration. They include:

- Performance-review test requiring MAFFT/FastTree (environment/tool availability).
- Baseline-extension/pangenome selection-set expectations needing reconciliation with current selection logic and required fixture evidence.
- Placement-warning flags test using the legacy `LOW_CONFIDENCE` label.
- Three novelty metrics tests with mocks/output assumptions not aligned with the expanded per-reference TSV output and current command/output wiring.
- Tiered Paper Trail QC expectation (current result includes `IsoFail` in review-template output; needs policy review, not a blind expected-value change).
- Two taxonomy screening tests missing `gtdb.fna.gz` fixtures.
- Secondary-baseline extension fixture not providing the current required baseline evidence structure.
- Baseline genome-availability test expects `PRIMARY` when current algorithm assigns `DIVERSITY_CANDIDATE`; requires decision-contract review.

## Important limitation

This iteration contains **test changes only; no production-code fixes are claimed**. The unresolved behavioural cases need a separate evidence-led review before changing the production logic or weakening assertions.
