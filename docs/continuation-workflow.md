# Versioned project continuation

`continue-project` now performs a cumulative **Continue → Reconcile → Reassess → Publish** cycle. It is intended to extend an existing project without replacing the historical database or assessment snapshots.


## What it does

1. Creates a timestamped version directory and snapshots the pre-continuation SQLite database.
2. Copies the previous Performance Review output tree when the database has a recorded, accessible Performance Review manifest.
3. Reconciles the cumulative partner catalogue (new registrations and metadata/status updates).
4. Imports optional partner metadata from Excel or CSV/TSV.
5. Imports completed genome results into the staged database. Genome records and their availability contribute to the reassessment and selection/pangenome metrics.
6. Optionally processes new/replacement marker sequences through the existing marker QC, provenance and replacement workflow.
7. Exports the cumulative active candidate marker collection and runs a full project Performance Review using the latest recorded assessment configuration.
8. Validates required assessment outputs and SQLite integrity.
9. Publishes the database only after the staged continuation has completed and validated. Failed runs retain their version directory for diagnostics but do not replace the current database.

Each continuation preserves a pre-run database snapshot and produces a new database and assessment snapshot. The timestamped version directories are retained; the `latest.json` pointer is updated only after publication.


## Inputs

Required:

- `--db`: current project database.
- `--partner-catalogue`: cumulative catalogue CSV/TSV.
- `--genome-results`: genome results CSV/TSV.

Optional:

- `--partner-workbook`: partner metadata workbook or CSV/TSV. CSV/TSV has no worksheet; the sheet option is ignored.
- `--marker-fasta`: new marker sequences, including sequences for previously markerless or failed isolates. Existing non-identical sequences enter the marker replacement decision workflow.
- `--marker-qc`, `--marker-review`: provenance/QC and manual review sidecars for the supplied marker FASTA.
- `--replacement-policy`: `prompt` (default), `approve_qc_pass`, `approve_all`, or `decline`. `approve_all` can accept lower-QC replacements and should be used deliberately. `approve_qc_pass` automatically accepts replacement markers but leaves lower-QC replacement decisions conservative.
- `--ref`, `--taxa`, `--main-ref`, `--mwl`, `--mwl-sheet`, `--threads`: overrides for inherited Performance Review configuration.
- `--baseline-fasta`, `--baseline-taxa-assignments`, `--chimera-ref`, `--anchors`, `--alt-ref`, `--alt-taxa`, `--alt-ref-name`, `--sequence-domain`: optional configuration overrides when a prior review's input paths are no longer available or configuration needs changing.
- `--version-root`: custom version storage location. Default is `<database-parent>/<database-stem>_versions`.

The continuation inherits assessment settings from the latest successful Performance Review manifest and its recorded command. If no prior Performance Review is found, `--ref` is required to establish the assessment configuration. Missing inherited input files are reported rather than silently ignored.


## Example

```bash
branchmanager continue-project \
  --db project.sqlite \
  --partner-catalogue Partner_Isolate_Catalogue.csv \
  --genome-results Genome_Results.tsv \
  --partner-workbook AarhusUniversity_DATA.csv \
  --marker-fasta new_or_resequenced_markers.fasta \
  --marker-qc marker_qc.tsv \
  --replacement-policy approve_qc_pass \
  --dry-run \
  -o continue_project_dry_run
```

Use `--dry-run` to stage imports and reassessment and generate the version/output files without publishing the current database. Remove it to publish the validated version.


## Version layout

A version folder includes the previous database snapshot, cumulative marker FASTA used for assessment, the new project database snapshot, the full assessment output, optional incoming-marker review output, continuation reports and `version_manifest.json`. The user-specified `-o` directory receives summary/import reports and principal assessment tables/dashboard for convenience.


## Marker replacements and genome status

New markers are checked against existing IDs. Identical sequences are retained with incoming provenance recorded. Non-identical sequences enter the established resequencing/replacement logic and failed or declined attempts remain auditable. The replacement policy controls automatic decisions; prompt mode is interactive.

Selection for sequencing and actually sequenced are distinct lifecycle states. Catalogue selection flags are reconciled, completed genome results are imported, and genomes are included in the reassessment's genome coverage, candidate-selection and pangenome calculations. Markerless genome records remain genome evidence only; the workflow does not fabricate marker sequences.


## Partner-owned IDs

The catalogue `partner_own_ids` field is retained as `partner_own_id`. Outputs and dashboard displays preserve both BranchManager/project IDs and partner IDs. Older records without partner-owned IDs remain valid.
