# Validation and ground truth

BranchManager's validation toolkit is an **additive, read-only facility**. It
is not called by Assistant, Filing Cabinet or Performance Review and does not
modify biological outputs, workflow state or the project SQLite database.

## 1. Create a provisional reference

Create a reference snapshot from a completed, reviewed run:

```bash
branchmanager ground-truth init \
  --run-dir /path/to/20261111_120106 \
  --out /path/to/ground_truth/rumen_gateway_reference \
  --label "Rumen Gateway reviewed cohort" \
  --db /path/to/project.sqlite
```

`--db` is optional when exactly one unambiguous project database is found in
the run directory. Supplying it explicitly is recommended, particularly when
an Assistant run contains staged or intermediate databases.

The command captures the required assessment and SAB summary outputs, plus
optional sequencing, novelty and marker-review outputs when present. It
records missing optional outputs, rejects blank or duplicate comparison keys,
checksums expected files, captures input provenance and records read-only
row-level fingerprints for selected project database tables.

A new reference is always `PROVISIONAL`. The source run is not automatically
considered scientifically correct. Review the expected outputs, database
state, input fingerprints, reference database releases and external tool
versions before approval. The input files recorded in the run manifest must
be accessible and their SHA256 fingerprints must verify before a reference
can be approved.

## 2. Approve a reviewed reference

After scientific and technical review:

```bash
branchmanager ground-truth approve \
  --reference /path/to/ground_truth/rumen_gateway_reference \
  --approved-by "Name or review role" \
  --reason "Reviewed against the accepted Rumen Gateway cohort"
```

Approval checks required outputs, row keys, input provenance and output
checksums. It writes a detached approval record containing the manifest and
expected-output digests. Subsequent edits to an approved manifest or expected
file are detected during validation. This is tamper-evident protection against
accidental or unauthorised changes when the approval record itself is kept
under the same access controls; it is not a digital signature against an
attacker able to rewrite every file in the reference directory.

Do not edit an approved reference. Create and review a new reference instead.

## 3. Validate a completed run

### Strict regression mode

```bash
branchmanager validate \
  --run-dir /path/to/new_run \
  --reference /path/to/ground_truth/rumen_gateway_reference \
  --db /path/to/project.sqlite \
  --mode STRICT_REGRESSION
```

This is the default. It is intended for a fixed cohort and reproducible
configuration. Input fingerprint changes, missing/added isolate keys and
changed decision-bearing values are reported as drift.

### Production compatibility mode

```bash
branchmanager validate \
  --run-dir /path/to/new_run \
  --reference /path/to/ground_truth/rumen_gateway_reference \
  --db /path/to/project.sqlite \
  --mode PRODUCTION_COMPATIBILITY
```

This mode is intended for evolving partner cohorts. It still verifies current
input files against the current run manifest. When the input fingerprint set
has legitimately changed, added or removed isolate keys and corresponding
database record additions/removals are reported as informational cohort
changes. Changes to values for shared isolates remain HIGH drift. If inputs
are unchanged, missing or added keys remain HIGH drift. This avoids treating
new partner submissions as software regressions without masking changes to
existing isolate decisions.

The validator distinguishes input provenance states:

- `VERIFIED`: required input files/directories exist and their recomputed
  checksums match the run manifest.
- `VERIFIED_FROM_RUN_MANIFEST`: the only missing required input is the
  temporary `project_database_before` state, and matching before/after
  database records are retained in the captured Performance Review manifest.
- `CHANGED`: the verified input fingerprints differ from the reference.
- `UNAVAILABLE`: required provenance or source files are unavailable.
- `INVALID`: a source file does not match the checksum recorded in its own
  run manifest.

Unavailable provenance cannot produce an authoritative PASS. Invalid
provenance is a critical failure. The captured workflow-manifest evidence is
also checked during approval and every subsequent validation.

## 4. Compared outputs and database state

The initial output set is:

- `sequence_assessment.tsv` (required)
- `selection_summary.tsv` (required)
- `sequencing_sets.tsv` (optional)
- `novelty_metrics.tsv` (optional)
- `cohort_review_application_audit.tsv` (optional)
- `cohort_joint_marker_qc.tsv` (optional)

Tables are compared by stable identifiers rather than row order. Numeric
values use the tolerance recorded in the reference manifest. Known run-specific
figure and source-manifest path columns are excluded from cell-wise comparison.
The database comparison checks SQLite integrity, foreign keys, schema version,
table columns and row-level fingerprints keyed by primary key for the selected
tables: `sequences`, `genome_records`, `isolate_status`, `sequencing_metadata`,
`marker_evidence` and `marker_review_decisions`. Tables without a primary key
are compared as order-independent row-hash collections.

## 5. Reports and exit codes

Reports are written to `RUN_DIR/validation/` by default; `--out` overrides this.

- `validation_summary.json` — machine-readable result and individual checks.
- `drift_report.tsv` — one row per comparison or invariant.
- `validation_report.html` — human-readable report.

Exit codes:

- `0`: approved reference; no HIGH or CRITICAL drift.
- `1`: HIGH drift detected.
- `2`: reference remains provisional.
- `3`: CRITICAL failure, including invalid input checksums or reference tampering.
- `4`: validation incomplete because required provenance is unavailable.

LOW/INFO differences are recorded and do not by themselves fail validation.

## 6. Safety and scope

Validation is observational. It opens SQLite databases read-only and writes
only to the validation report directory. It does not change project state,
repair results, rewrite TSVs or automatically update ground truth. An explicit
`--db` can point to a database outside the run directory; it is still opened
read-only.

The first reference created from a production run is a candidate baseline,
not certified biological truth. It should be reviewed by the project team
before approval. The toolkit currently provides strict regression and
production-cohort compatibility checks; rank-aware semantic taxonomy
comparison and dashboard history remain future extensions.
