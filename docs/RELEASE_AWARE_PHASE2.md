# Release-aware workflows — Phase 2 (1.4.1)

Phase 2 deepens the release context integration while retaining legacy command-line paths.

## Release-aware continuation

`continue-project --release <release-id-or-path>` resolves the database to the selected release database. If `--out` is omitted, reports go under `<release>/continuation/`; if `--version-root` is omitted, continuation snapshots are placed under `<release>/continuation_versions/`. Both paths are validated to remain within the selected release. Explicit output/version-root paths remain supported but must be release-contained.

## Workflow provenance

Release manifests now receive workflow start/completion/failure events and run records, including the invoked command/arguments and SHA256 hashes for file arguments that can be identified as existing files. This is best-effort provenance; arbitrary inputs passed through config values or nested options are not necessarily captured.

## Status

Workflow status records are tracked separately from overall release status. A completed workflow is marked `COMPLETE`, while the release remains `IN_PROGRESS` until a separate validation/review step. Workflow exceptions record `FAILED`.

## CLI compatibility

For database-backed commands carrying `--release`, the parser allows database and output arguments to be omitted because release context supplies them. Legacy standalone commands still require the same explicit paths as before.
