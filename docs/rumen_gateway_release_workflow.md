# Rumen Gateway release workflow (corrected 22 November 2026)

Use release-specific paths and the project profile. Each release owns its database, inputs, logs, and outputs. Do not write to the live project database during a release run.


## Fresh initial release

1. Initialise/validate the project once, if not already initialised.
2. Create a fresh release.
3. Stage only the initial marker datasets (QUB, UoG, BSI) and initial genome datasets (AAFC, IAS_CAAS) into that release's `inputs/`.
4. Build the baseline/Filing Cabinet once for this fresh release.
5. Register initial genome-only isolates and import their genome QC/GTDB results. Dry-run the results import first and inspect its report.
6. Run `assistant --cohort-manifest ... --cohort-prepare-only` using the release, project profile, catalogue and release-local manifest.
7. Review the generated marker-review template. Reconcile decisions by `sequence_id`; save the reviewed TSV in the release's `inputs/reviews/`.
8. Run `assistant --cohort-resume ... --marker-review <review.tsv>` against the same release and cohort preparation folder. This operation **does perform the joint marker assessment**. Its internal Performance Review outputs are initially written below `cohort/01_joint_assistant/03_performance_review_hiring_panel/`. The corrected BranchManager build also promotes assessment files to the release-level `assessment/` and the HTML dashboard to `dashboard/`, and marks the assessment/dashboard workflow states complete when the required promoted files are present.
9. Inspect the release-level assessment tables and dashboard. Do not infer scientific sign-off from workflow status alone.
10. Run `branchmanager release validate --project "$PROJECT_ROOT" --release "$RELEASE_ID"` and inspect the validation report.
11. Only mark the release `READY` after validation and review. BranchManager now refuses READY if assessment/dashboard workflow stages are not COMPLETE or if these files are absent/empty: `assessment/sequence_assessment.tsv`, `assessment/selection_summary.tsv` (complete pending SAB decision set; QC/resequencing cases included, completed sequencing commitments excluded), `dashboard/performance_review_dashboard.html`.


## Continuation release

1. Choose the most recent validated and reviewed release as `PARENT_ID`.
2. Create a child release with `--parent-release "$PARENT_ID" --copy-database`.
3. Stage only new marker/genome data, manifests, catalogues, and review files in the child release's `inputs/`.
4. Register/import new genome records/results. Dry-run result imports and inspect them before writing.
5. Prepare and review only the new marker cohort, then resume it in the child release. Do not rebuild baselines or re-import inherited data.
6. Confirm cumulative assessment and dashboard outputs are present in the child release's top-level `assessment/` and `dashboard/` folders. Review the differences from the parent.
7. Validate and scientifically review the child. Mark READY only when validation, assessment, dashboard and required output files are complete.


## Release readiness files

A READY release must have non-empty:


- `assessment/sequence_assessment.tsv`
- `assessment/selection_summary.tsv` (complete pending SAB decision set; QC/resequencing cases included, completed sequencing commitments excluded)
- `dashboard/performance_review_dashboard.html`

and workflow statuses `assessment=COMPLETE`, `dashboard=COMPLETE`, and `validation=COMPLETE`. SQLite integrity alone is not sufficient for release readiness.


## Failure handling

If a step fails, preserve its log and outputs. Do not blindly repeat database-writing steps. Inspect the run manifest and release manifest first. A cohort-preparation completion is not the same as a completed joint assessment; use the cohort-resume run record and the release-level assessment/dashboard outputs as confirmation.

## Project-wide current publication

After the release passes validation and is marked `READY`, BranchManager automatically updates:


```text
Releases/BranchManager_Project/current/
├── project.sqlite
├── assessment/                 # full assessment plus neighbourhood/support files
├── dashboard/                  # dashboard and any dashboard-local assets
└── current_release.json         # source release ID and database checksum
```

This is the replaceable latest publication. Historical releases in `versions/` are not changed. Do not use `current/` as an archival location; use the release directory for provenance and reproducibility.
