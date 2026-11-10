# Development and testing

This guide is for contributors working on the BranchManager source repository.

## Source layout

```text
src/branchmanager/
├── cli/                 # CLI parser and entry-point package
├── data/                # Bundled interface and reference-anchor data
├── db/                  # Database schema and persistence
├── examples/            # Packaged examples and templates
├── pipeline/             # QC, classification, novelty, trees and selection
├── tools/                # Standalone helper commands
├── utils/                # Shared sequence, file and text utilities
├── cli_core.py           # Command implementations and orchestration
├── sanger_studio.py      # Local Sanger Studio server and API
└── ...
```

The package is installed from `src/` using the project configuration in `pyproject.toml`.

## Development environment

```bash
conda create -n branchmanager-dev python=3.12 -y
conda activate branchmanager-dev
conda install -c conda-forge biopython vsearch mafft fasttree -y
python -m pip install -e .
```

Install any additional tools needed for the workflow being developed, such as IQ-TREE, Barrnap or GTDB-Tk.

## Run tests

Run the complete suite from the repository root:

```bash
python -m pytest
```

Run a focused test module:

```bash
python -m pytest tests/test_selection_taxonomy_gating.py -q
```

Run a group of related tests:

```bash
python -m pytest tests/test_sanger_studio_extensions.py tests/test_sanger_studio_zip_upload.py -q
```

Not every optional external tool or reference database is available in every development environment. Tests that require external executables or large reference resources should clearly identify those requirements; do not report skipped or unexecuted tests as passing.

## Static checks

Compile the package:

```bash
python -m compileall -q src/branchmanager
```

If editing the Sanger Studio browser interface, run the JavaScript syntax check where Node.js is available:

```bash
node --check src/check_sanger_studio.js
```

This checks the browser-side JavaScript syntax helper included in the repository.

## Versioning

The package version is declared in `pyproject.toml` and exposed as `branchmanager.__version__`. Keep these values synchronised when preparing a release. The CLI obtains its version from the package metadata:

```bash
branchmanager --version
```

Do not hardcode a separate CLI version string.

## Scientific regression tests

Changes to taxonomy, QC, novelty, pangenome grouping or selection should include tests that capture the biological contract, not only the internal function output.

Examples include:

- Domain routing must not become pangenome taxonomy.
- GTDB taxonomy is authoritative for pangenome grouping.
- Unresolved GTDB taxonomy must not silently fall back to another classifier.
- Genus-level evidence must not be promoted to species without support.
- Hierarchy eligibility and selection recommendations must agree.
- Alignment identity and query/reference coverage must remain distinct.
- Failed QC and manual review decisions must not be overwritten by automatic reruns.

## Release hygiene

Before creating a release:

1. Review `git status` and commit all intended source and documentation changes.
2. Confirm `pyproject.toml` and `__version__` agree.
3. Run focused tests and the full test suite in an environment with the required dependencies.
4. Build and inspect the source distribution and wheel.
5. Test installation in a clean Conda environment.
6. Run `branchmanager --version`, `branchmanager --help` and `branchmanager it-desk`.
7. Confirm packaged HTML, anchors and example resources are included.
8. Tag the exact commit that was tested.

A release should not claim a full-suite pass unless the full suite completed successfully.
