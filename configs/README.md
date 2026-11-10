# Configuration profiles

This directory contains JSON profiles for project-specific paths, reference datasets, baselines and command defaults.

## Files

| File | Purpose |
|---|---|
| `project.template.json` | General project profile template. |
| `assistant.template.json` | Assistant-focused template for recurring partner batches. |
| `assistant.json` | Example Assistant profile. |
| `rumen_gateway.json` | Rumen Gateway-specific project configuration. |
| `Rumen_Gateway_Assistant.json` | Additional Rumen Gateway Assistant configuration. |

Check the actual paths and settings in each file before use. A profile is not portable without the corresponding reference, baseline, metadata and database files.

## Recommended use

1. Copy the relevant template to a project-specific filename.
2. Set the project directory and data paths.
3. Configure the marker/domain, references, taxonomy sources and baseline datasets.
4. Set project-specific QC and selection defaults.
5. Validate the profile with the relevant command's `--help` and the `it-desk` preflight.
6. Keep the profile and referenced data versions with the analysis archive.

For example:

```bash
cp configs/assistant.template.json configs/my_collection.json
```

Use the profile:

```bash
branchmanager assistant \
  --config configs/my_collection.json \
  --sample-map batch_map.tsv \
  --dataset Collection_01 \
  -o runs/Collection_01
```

Explicit command-line arguments override profile defaults.

## Path handling and privacy

Relative paths are resolved using the profile's configured project directory or the profile file location, according to the profile settings. Environment substitutions can be used for local paths.

Do not commit credentials, personal data, confidential partner metadata or restricted-access information to a shared profile.

## Biological scope

A profile changes paths and defaults; it does not make every workflow biologically valid for every marker or domain. Before applying BranchManager to a new collection, validate the marker, reference databases, taxonomy authority, anchors, QC thresholds, chimera reference and selection unit.

See [Collection profiles and validation](../docs/collection-profiles-and-validation.md).
