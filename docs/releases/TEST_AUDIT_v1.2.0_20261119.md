# BranchManager v1.2.0 Test Audit — 19 November 2026

## Baseline

Audited archive: `BranchManager(20261119-121036).zip`.

The supplied pytest log reports **471 passed, 27 failed, 1 skipped** (499 tests collected). Failures are concentrated in `tests/test_taxonomy_and_db.py` (26) and `tests/test_v1_cohort.py` (1). The focused regression subset rerun from the archive in the audit environment reports **59 passed**.

## Confirmed focused regression coverage

- Sanger minimum overlap policy (40 bp)
- Sanger attempt manifests and attempt-aware handling
- Sanger Studio browser interactions and regressions
- Resequencing/read selection
- Trace-channel mapping
- Selection taxonomy gating
- Pangenome domain fallback
- Partner taxonomy validation

No failures occurred in this focused subset.

## Remaining failures: assessment

### Likely stale test interface or expectations

- `test_background_check_reports_high_quality_taxonomy_disagreements`: mock hit records contain two values, but current background-check classification records require eight fields. Update the fixture to the current record contract; do not loosen production unpacking without evidence that a two-field record is a supported API.
- `test_interview_resolves_mailroom_output_and_delegates_to_shared_qc`: patches an obsolete CLI symbol/path and expects delegation that the current `cmd_interview` does not expose in that form. Update the mock to patch the current handler boundary and assert the intended delegation.
- Multiple `test_performance_review_*` tests: patches `branchmanager.cli.qc` / `branchmanager.cli.classify`, which no longer resolve to the current CLI module architecture. Patch the modules/functions actually imported by `cli_core` or use dependency injection at a stable boundary.
- `test_merge_combined_taxonomy_rows_overrides_existing_with_run_results`: expects numeric confidence from the classifier, but current merge policy gives partner catalogue taxonomy precedence and uses a provenance label rather than numeric classifier confidence for the catalogue identity. The test should assert taxonomy source/provenance, not force a numeric value unsupported by that source.
- `test_partner_metadata_separates_selection_commitment_from_available_genome`: expects legacy `selection_provided`; current normalised schema exposes `selected_for_sequencing`, `genome_available`, and `selected_for_wgs`. Align assertions with the current schema or explicitly add a documented compatibility alias if external consumers still require it.
- Novelty-metrics tests expecting the old short TSV layout should assert against current named columns and their values, rather than old positional substrings.
- Several Paper Trail tests have stale/incomplete fixtures (missing `gtdb.fna.gz`) or stale test code (`NameError: lines`). Repair fixture setup and test helpers before evaluating behaviour.
- `test_paper_trail_uses_tiered_ambiguity_policy`, baseline-extension role expectations, placement-warning flags, secondary-baseline labels, and sequencing-set roles appear to assert earlier policy/labels. Review each expected value against the current evidence-selection policy before changing code.

### Environment / integration dependencies

- Tree-building tests cannot execute without MAFFT and FastTree in the local audit environment; install required tools or isolate external-tool integration tests and skip with an explicit dependency check.
- `test_cohort_resume_joint_assessment_and_publish` fails IT Desk before the workflow can run because required external tools (`vsearch`, `mafft`, `FastTree`) are absent and the configured reference FASTA path does not exist. This is an environment/setup failure, not evidence of a workflow logic regression. Provide these dependencies and a valid test reference, or mock IT Desk only in a unit test while retaining a separate integration test.

## Recommended action

Do not change scientific/selection logic to satisfy these tests until their assertions and fixtures are reconciled with the current v1.2.0 contracts. Update stale mocks and assertions, preserve strict taxonomy evidence rules, and make external-tool integration requirements explicit. Then rerun the complete suite in the intended Conda environment with `vsearch`, `mafft`, and `FastTree` available and valid test references.

## Audit changes

This audit report is documentation-only. No production code or scientific decision logic was modified during this review.
