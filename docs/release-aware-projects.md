# Release-aware project management

BranchManager now provides project/release lifecycle commands while preserving the existing workflow commands and their `--db` / `--out` interfaces.


## Initialise a project once

```bash
branchmanager project init --root /path/to/Database_Live --name "Rumen Gateway" \
  --config-source /path/to/project_config.json
```
This writes `project_manifest.json` and records the project config and release-store location. By default the release store is `Releases/BranchManager_Project`; use `--releases-dir` to choose another relative or absolute store.

Project-local release stores can move with the project, including stores recorded
as absolute paths by older BranchManager versions. Pass the current project path
to `--project` or `--root`; existing manifests are resolved without rewriting
them. The `latest` alias and release database paths follow the relocated store.

An absolute store outside the project retains its configured location. If that
store moves too, specify its current location for each command, for example:

```bash
branchmanager release info --project /new/project --releases-dir /new/store --release latest
branchmanager project list-releases --root /new/project --releases-dir /new/store
branchmanager continue-project --project-root /new/project --releases-dir /new/store \
  --release latest --partner-catalogue partners.csv --genome-results genomes.tsv
```

Release-aware workflows also accept `--project-root` when the project cannot be
discovered from a release path or config. Relative paths in a project config use
the current project root in release context; absolute paths to external
scientific inputs remain unchanged and must be updated separately if they move.

## Create a fresh release

```bash
branchmanager release create --project /path/to/Database_Live
```
An optional `--release-id` can be supplied. Existing release IDs are never overwritten. The release receives its own database path, config snapshot, workflow directories, `release_manifest.json`, and the project store's `latest.json` pointer.

## Adopt an existing release directory

For a pre-existing release directory that has not yet been given a release manifest:

```bash
branchmanager release adopt --project /path/to/Database_Live --release-id 20261120_164056
```
Adoption preserves existing contents, creates missing standard workflow directories and records metadata. It refuses to overwrite a release manifest. Configure the project's release-store path during `project init` to match the existing release-store layout.

## Release-aware workflow commands

Database-backed workflows accept `--release <release-id-or-directory>`. This resolves the release database from the manifest. A conflicting explicit `--db` is rejected. Unless an output path is explicitly provided, supported commands choose a release-local output directory. Explicit output paths must remain within the release directory.

Examples:

```bash
branchmanager filing-cabinet --config /path/to/project_config.json \
  --release /path/to/release --build-tree --threads 10
branchmanager assistant --config /path/to/project_config.json \
  --release /path/to/release --cohort-manifest /path/to/manifest.csv \
  --cohort-prepare-only
```

The project config continues to define shared reference paths and scientific defaults. The release manifest owns the release ID, database path, lifecycle status and workflow status.

## Inspect and validate

```bash
branchmanager project info --root /path/to/Database_Live
branchmanager project validate --root /path/to/Database_Live
branchmanager project list-releases --root /path/to/Database_Live
branchmanager release info --project /path/to/Database_Live --release latest
branchmanager release validate --project /path/to/Database_Live --release 20261120_164056
branchmanager release status --project /path/to/Database_Live --release 20261120_164056 --set IN_PROGRESS
```

`release validate` checks the manifest, standard release directories, database path containment and SQLite integrity (when a database exists). A successful validation records `validation=COMPLETE` in the release workflow status. A release cannot be marked `READY` unless assessment and dashboard stages are recorded as `COMPLETE` and these non-empty release-local files exist:

- `assessment/sequence_assessment.tsv`
- `assessment/selection_summary.tsv`
- `dashboard/performance_review_dashboard.html`

For cohort-resume runs, the joint assessment is generated under `cohort/01_joint_assistant/03_performance_review_hiring_panel/`; BranchManager now promotes the assessment tables and dashboard to the release-level `assessment/` and `dashboard/` directories and records their workflow completion. The release manager never removes a release directory or overwrites a pre-existing release ID.

`continue-project` intentionally remains on its existing continuation/versioning path for this initial release-aware implementation. It will be unified with ReleaseManager in a follow-on refactor so its existing transactional publication and snapshot safeguards can be preserved and tested before changing its lifecycle.

## Automatic cohort review reconciliation

During cohort PREPARE, BranchManager first looks for an explicitly saved review-decision TSV matching `inputs/reviews/cohort_marker_review_decisions*.tsv` in the release directory or its ancestors. If none is found, it searches previous sibling version folders for the most recent populated `cohort_marker_review_template.tsv` and carries decisions forward by canonical `sequence_id`. You can override source discovery with `--marker-review`. If multiple explicit decision TSVs are found, PREPARE stops and asks you to select one rather than guessing.

This is an incremental review workflow: existing decisions are restored for matching IDs; new isolates remain undecided and are available for review; IDs present in the old cohort but absent from the new cohort are reported as informational rather than treated as a fatal mismatch. Decision mismatches, duplicate IDs, or malformed source/template rows still fail reconciliation. The audit writes:


- `cohort_review_reconciliation.tsv` — new/removed IDs and any discrepancies.
- `cohort_review_reconciliation.json` — status and detailed audit fields, including carried decisions.


If no prior decision source is found, the report status is `NO_SOURCE` and the template remains available for manual review.

## Project-wide current publication

When a release is successfully validated and marked `READY`, BranchManager publishes a replaceable project-wide bundle beside `versions/` at `Releases/BranchManager_Project/current/`. This contains `project.sqlite`, the complete `assessment/` directory (including neighbourhood/support files), the complete `dashboard/` directory, and `current_release.json` identifying the source release and database checksum. The publication is staged and swapped into place; a failed publication rolls the release status back rather than leaving it READY. The previous release directories remain immutable under `versions/`. Consumers should use `current/` for the most recently published release and `versions/<release_id>/` for historical/reproducible outputs.
