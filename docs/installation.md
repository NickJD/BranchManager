# Installation and dependencies

## Supported Python version

The current package metadata requires **Python 3.12 or newer**. Use a dedicated Conda environment to keep BranchManager and its dependencies separate from other projects.

## Recommended Conda installation

Install [Miniforge](https://github.com/conda-forge/miniforge), then run:

```bash
conda create -n branchmanager python=3.12 -y
conda activate branchmanager
```

Install the command-line tools used by the production workflow:

```bash
conda install -c conda-forge biopython vsearch mafft fasttree -y
```

Install BranchManager from PyPI once the v1.0.0 package has been published:

```bash
python -m pip install --upgrade BranchManager
```

For a development checkout:

```bash
git clone https://github.com/NickJD/BranchManager.git
cd BranchManager
python -m pip install -e .
```

Use the same active Conda environment for installation and execution.

## Required and optional tools

The `it-desk` preflight checks the Python environment and important external tools.

| Component | Role | Status |
|---|---|---|
| Python 3.12+ | Runtime | Required by the current package metadata |
| Biopython (`Bio`) | Sequence and chromatogram handling | Required by preflight |
| VSEARCH | Reference searches, sequence comparisons and related workflows | Required by preflight |
| MAFFT | Multiple sequence alignment | Required by preflight |
| FastTree | Default tree-building option | Required unless an IQ-TREE method is configured |
| IQ-TREE 2/3 | Alternative tree building | Optional; required when selected as the tree method |
| Barrnap | Genome-derived rRNA extraction | Optional for marker-only workflows |
| GTDB-Tk | Genome taxonomy assignment | Optional unless the genome-taxonomy workflow is used |
| CheckM2 | Genome quality assessment | Used by genome evidence workflows where configured |
| OpenPyXL | Excel workbook handling | Installed as a package dependency; needed for workbook workflows |

The package's Python dependencies are declared in `pyproject.toml`. Some external command-line programs and the Biopython module are checked separately by `it-desk`; installing the Python package alone should not be assumed to install every external executable.

## Check the installation

```bash
branchmanager --version
branchmanager --help
branchmanager it-desk -o preflight_check
```

The preflight report contains `it_desk_report.tsv` and `it_desk_summary.json`. A failed required check should be resolved before running a production workflow. Optional checks may be reported as warnings when the corresponding functionality is not being used.

For a project-specific check, supply the database, references and output path using the options shown by:

```bash
branchmanager it-desk --help
```

## Updating BranchManager

Activate the environment first:

```bash
conda activate branchmanager
python -m pip install --upgrade BranchManager
```

If you installed from a source checkout in editable mode, update that checkout with Git and reinstall only if package metadata or dependencies have changed.

## Troubleshooting

- **`branchmanager: command not found`**: activate the correct Conda environment and confirm that installation completed in that environment.
- **`No module named Bio`**: install Biopython in the active environment (`conda install -c conda-forge biopython`).
- **VSEARCH/MAFFT/FastTree not found**: install the missing executable and check that it is available on `PATH`.
- **Wrong version displayed**: check `which branchmanager` (Linux/macOS) or `where branchmanager` (Windows), then check `branchmanager --version`. Multiple Python environments can contain different installations.
- **Reference file errors**: confirm that the path exists, that the file is readable and that the reference is in the expected FASTA format.

The Rumen Gateway configuration and bundled rumen/gut anchors are project-specific examples. Other collections should configure and validate appropriate references and anchors before interpreting biological results.
