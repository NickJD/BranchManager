# BranchManager documentation

This directory contains the maintained user and developer documentation. The root [README](../README.md) is the entry point; use the guides below for detail.

## User guides

- [Installation and dependencies](installation.md)
- [Sanger Studio](sanger-studio.md)
- [Workflows and CLI commands](workflows.md)
- [Taxonomy and pangenome interpretation](taxonomy-and-pangenome.md)
- [Outputs and interpretation](outputs-and-interpretation.md)
- [Collection profiles and validation](collection-profiles-and-validation.md)

## Developer guide

- [Development and testing](development.md)

## Release notes

- [BranchManager v1.0.0](releases/v1.0.0.md)

## Documentation policy

The dated patch notes supplied during development are useful historical records, but they are not all current user instructions. Where an older patch note conflicts with these maintained guides or with the current implementation, use the maintained guide and the current source code.

In particular, historical notes describing partner-catalogue genus/species taxonomy as the pangenome authority are superseded. The current pangenome grouping contract uses GTDB-derived taxonomy; the partner catalogue's `taxon` field controls domain/pathway routing.
