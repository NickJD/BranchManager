# Sanger Studio and release validation fixes — 22 November 2026

## Sanger Studio

- The isolate list now grows to fill the available sidebar space and scrolls internally; the inspection action bar is non-shrinking so the **Inspect and assemble this isolate** button remains visible.
- Saving an unchanged retained candidate from an attached review TSV preserves its recorded candidate ID/type/length/source. Historical candidates are no longer rejected merely because they were not regenerated in the current Studio session.
- Newly changed candidate selections remain validated against the latest inspection and still require reviewer identity.
- Attaching a marker review with no reads loaded now prompts for AB1/ABI files and reports read-to-review reconciliation. Exact canonical IDs, explicitly supplied partner aliases, ambiguous aliases, and unmatched IDs are reported separately; aliases are not silently renamed.

## Release management

- Release validation now reports a missing SQLite project database as a validation issue.
- Validation results are recorded in the release manifest.
- A release cannot be marked `READY` until validation has been recorded; the release is revalidated when the status is changed to `READY`.
- Workflow event recording is protected by a manager-level lock to prevent concurrent threads from overwriting manifest updates.

## Verification

- Sanger Studio / release manager focused regression tests: passed.
- Candidate review, sidebar layout, review-first reconciliation, Sanger Studio browser tests, and release lifecycle/phase-two tests: **53 passed**.
- Cohort / evidence / selection / validation toolkit tests: **77 passed, 2 skipped**.
- Python compile check and JavaScript syntax check passed.
- Package wheel built successfully using the installed build dependencies (`--no-build-isolation`).
- The complete test suite was attempted but did not finish within the execution time available; no complete-suite pass is claimed.
