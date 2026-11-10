# Examples and templates

This directory contains small example datasets, templates and helper scripts. They are starting points for preparing or testing inputs; they are not authoritative biological reference datasets.

| File | Purpose |
|---|---|
| `project_partner_metadata.example.tsv` | Example partner/isolate metadata structure. |
| `external_genomes.example.tsv` | Example external genome evidence table. |
| `benchmark_labels.template.tsv` | Template for independently adjudicated candidate benchmark labels. |
| `sample_diversity_metadata.tsv` | Example diversity metadata. |
| `rolling_partner_test.sh` | Example helper for a rolling partner workflow. |
| `Hungate_Plus_NewSequences_example.jpg` | Example visual output. |

Check the header names and expected values against the current command's `--help` and the relevant workflow documentation before adapting a template.

## Benchmarking

The benchmark label template is intended for independent outcome labels. Do not label a candidate as successful merely because BranchManager selected it. Define the outcome before examining rankings, include the full candidate pool and compare against a documented baseline under the same budget.

See the collection-profile and validation guide for a more complete scientific benchmarking plan.
