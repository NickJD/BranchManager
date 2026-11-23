# BranchManager 1.4.0 — Release-aware project management

## Added
- `branchmanager project init/info/validate/list-releases` project management commands.
- `branchmanager release create/adopt/info/validate/status` release lifecycle commands.
- Timestamped release creation and collision protection (never overwrites an existing release ID).
- Existing-directory adoption for releases that predate manifests; adoption preserves files and refuses to overwrite an existing manifest.
- Release manifests with project/release paths, status, timestamps, parent release, workflow status, and event history.
- Project `latest.json` pointer and `latest` release alias.
- SQLite backup API for creating a release from a parent release database; integrity check runs after backup.
- Standard release subdirectories and release-local database path validation.
- Release context `--release` for selected database-backed workflow commands; release DB is resolved from manifest, conflicting explicit DB path rejected, and default outputs are release-local.
- Release path validation on explicit output paths.

## Legacy compatibility
- Existing project config defaults remain supported.
- Existing standalone workflow commands continue to use explicit `--db`/`--out`.
- `continue-project` remains on its established versioning path pending a careful integration that preserves transaction and publication safeguards.

## Testing
- Added `tests/test_release_manager.py` covering project init, release isolation, ID collision, path containment, manifest status changes, database validation, SQLite database copy, adoption and latest alias.
- Targeted regression command used: `PYTHONPATH=src pytest -q tests/test_release_manager.py tests/test_cli_parser_regressions.py tests/test_continue_project_versioning.py tests/test_sanger_studio_attempt_manifest.py tests/test_dashboard_cli.py`.
