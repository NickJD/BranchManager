# Workflows and CLI commands

BranchManager can be run as a sequence of explicit stages or through the integrated `assistant` command. The appropriate route depends on whether you are preparing a new collection, processing one partner batch, reviewing an existing project or importing genome evidence.

Run `branchmanager --help` for the installed command list and `branchmanager <command> --help` for command-specific options.

## Recommended preflight

Before a production run:

```bash
branchmanager it-desk -o preflight_check
```

Resolve failures for required tools, modules and input files before proceeding.

## Sanger read preparation

A common Sanger route is:

1. **Mailroom** – inventory the delivery and reconcile reads with partner metadata.
2. **Interview** – perform standalone chromatogram QC and resequencing triage if needed.
3. **Onboarding** – validate isolate IDs and partner metadata.
4. **Paper Trail** – trim, assess and assemble reads into marker sequences.

Sanger Studio provides a graphical way to inspect and prepare reads. The command-line stages remain useful for reproducible batch processing and automated workflows.

Example Mailroom command:

```bash
branchmanager mailroom \
  --read-dir All_AB1 \
  --metadata supplier_metadata.csv \
  --dataset Partner_01 \
  -o mailroom/Partner_01
```

Review `ab1_map.tsv`, the inventory and the Mailroom report before using the map downstream.

## Establish a project baseline

Use Filing Cabinet to register a cultured baseline and, where configured, classify it and build the project backbone tree:

```bash
branchmanager filing-cabinet \
  --fasta baseline.fasta \
  --db project.db \
  --dataset Baseline \
  --ref gtdb_representatives.fna.gz \
  --classify \
  --build-tree \
  --threads 8 \
  -o filing_cabinet
```

The baseline should represent the collection against which redundancy and project-relative novelty are assessed. Keep the source FASTA, taxonomy assignments and provenance with the project.

## Run a partner batch with Assistant

The `assistant` command orchestrates the standard raw-trace-to-assessment workflow. A minimal illustrative command is:

```bash
branchmanager assistant \
  --sample-map ab1_map.tsv \
  --partner-metadata partner_metadata.tsv \
  --db project.db \
  --dataset Partner_01 \
  --ref gtdb_representatives.fna.gz \
  --threads 8 \
  -o runs/Partner_01
```

In practice, use a project profile to supply stable paths, baselines, taxonomy references, MWL and selection settings. Keep the batch-specific sample map, dataset label and output directory explicit.

## Process multiple partner batches together

For a multi-partner cohort, prepare a CSV/TSV manifest with one row per batch. A manifest can identify the partner, dataset, sample map, read directory and optional workbook. Paths are resolved relative to the manifest unless absolute.

Example:

```text
partner_id  dataset  sample_map  read_dir  partner_workbook
PartnerA    Batch_A  mailroom/A/ab1_map.tsv  reads/A  metadata/A.xlsm
PartnerB    Batch_B  mailroom/B/ab1_map.tsv  reads/B  metadata/B.xlsm
```

Run:

```bash
branchmanager assistant \
  --config configs/rumen_gateway.json \
  --cohort-manifest cohort.tsv \
  --dataset PROJECT_COHORT \
  -o runs/PROJECT_COHORT
```

Each batch is prepared and QC'd separately; passing marker evidence is then assessed together. Review the combined marker-review template and record decisions in a separate copy. Do not edit a generated template in place if the workflow will regenerate it.

## Later project stages

| Command | Use |
|---|---|
| `performance-review` | Assess new marker sequences against project baselines and references. |
| `quarterly-review` | Nominate a later sequencing tranche under a defined budget. |
| `status-meeting` | Import factual isolate lifecycle updates. |
| `records-update` | Import completed genome, QC, GTDB and ANI evidence. |
| `import-genome-16s` | Register genome-derived 16S marker evidence. |
| `import-genome-evidence` | Import genome FASTA and associated genome evidence using configured tools. |
| `import-genome-results` | Import precomputed genome result tables. |
| `annual-report` | Generate cumulative project reporting. |
| `org-chart` | Extract a focused tree for a taxon. |
| `label-maker` | Regenerate iTOL annotation files. |

## Regenerate an existing run dashboard

Use the `dashboard` command to rebuild the interactive HTML dashboard from reports
already present in a BranchManager run directory. This is a reporting-only operation:
it does not rerun marker QC, taxonomy assignment, novelty searches, phylogenetic
analysis, pangenome assessment or candidate selection, and it does not require the
project database.

Pass either the run output directory or its `assessment/` subdirectory:

```bash
# Supply the run root (normally contains assessment/ and run_manifest.json)
branchmanager dashboard /path/to/run_output

# Or supply the assessment directory itself
branchmanager dashboard /path/to/run_output/assessment
```

The command requires `selection_summary.tsv` in the supplied directory or in its
`assessment/` child. It writes the regenerated file to the run root:

```text
run_output/
├── assessment/
│   ├── selection_summary.tsv
│   ├── sequence_assessment.tsv
│   ├── novelty_metrics.tsv
│   ├── pangenome_hierarchy.tsv
│   └── genome_collection_assessment.tsv
├── run_manifest.json
└── performance_review_dashboard.html
```

`selection_summary.tsv` is the required input. The other report files and
`run_manifest.json` are optional to the command; when absent, BranchManager logs
warnings and generates the dashboard with the evidence available. Consequently,
a dashboard generated from incomplete inputs may omit corresponding sections or
details. For the nearest cultured and partner-project matches, retain
`assessment/novelty_metrics.tsv`: it provides the comparison-specific nearest-hit
identifiers, identity, query coverage and aligned length.

The command overwrites the dashboard HTML at the destination. It does not modify
the source TSV files or rerun the underlying analyses. To inspect command-specific
help, run:

```bash
branchmanager dashboard --help
```

## Profiles

Profiles are JSON files that keep project paths and workflow defaults together. Start from:

- `configs/project.template.json` for a broader project profile.
- `configs/assistant.template.json` for an Assistant-focused profile.
- `configs/rumen_gateway.json` as a Rumen Gateway-specific example.

Explicit command-line arguments override profile defaults. Do not place credentials or sensitive participant data in shared profiles. See [Configuration files](../configs/README.md) and [Collection profiles and validation](collection-profiles-and-validation.md).

## Reproducibility

For each run, retain:

- The exact BranchManager version.
- The project profile and its checksum/run manifest.
- Input sample maps and metadata tables.
- Reference and baseline FASTA files and taxonomy tables, including release/version and checksums.
- Relevant external-tool versions.
- QC settings and any manual review decisions.
- The complete output directory.

Do not overwrite original input data or generated review templates. Use separate reviewed copies and preserve the audit trail.
