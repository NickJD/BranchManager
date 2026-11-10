# Sanger Studio

Sanger Studio is a local browser-based interface for inspecting Sanger chromatograms, organising reads and reviewing marker-sequence processing. It uses BranchManager's existing Paper Trail processing and assembly engine; it is not a separate base-caller or independent assembler.

## Start Sanger Studio

Activate the BranchManager environment and launch:

```bash
conda activate branchmanager
branchmanager sanger-studio
```

The application starts a local web server and normally opens the browser automatically. To choose a port or prevent automatic browser opening:

```bash
branchmanager sanger-studio --port 8765 --no-browser
```

Keep the terminal session running while using the interface. Stop the local server with `Ctrl+C`.

## Input options

Sanger Studio supports several ways to provide reads:

- Select individual `.ab1` or `.abi` chromatogram files.
- Select a directory containing chromatograms; nested directories are searched.
- Upload a ZIP archive containing chromatograms, including files inside nested folders.
- Load an existing sample map or cohort-preparation output where appropriate.

### Upload a ZIP archive

Choose the ZIP upload option in the input panel and select the archive. You do not need to extract it manually.

Example archive:

```text
Sanger_batch.zip
└── Sanger_batch/
    ├── Partner_A/
    │   ├── QUB_001_27F.ab1
    │   ├── QUB_001_1492R.ab1
    │   └── QUB_002_27F.ab1
    └── Partner_B/
        ├── QUB_002_1492R.ab1
        └── QUB_003_27F.ab1
```

Sanger Studio searches the extracted directory tree for `.ab1` and `.abi` files. Non-chromatogram files are not treated as reads. The upload implementation applies limits and path-safety checks, including a 500 MB compressed upload limit and limits on archive entries and chromatogram count.

The upload process creates a sample map using the existing filename-based inference rules. **Check the inferred isolate identifiers, primer names and read directions before processing.** Filename inference is a convenience, not a substitute for checking the sequencing manifest.

## File naming and sample maps

Clear, consistent filenames make read grouping easier. For example:

```text
QUB_001_27F.ab1
QUB_001_1492R.ab1
QUB_002_27F.ab1
QUB_002_1492R.ab1
```

The exact filename convention depends on the sample map and primer definitions used by the project. If names are ambiguous, edit or provide the sample map rather than relying on automatic inference.

Keep stable isolate identifiers distinct from supplier sequencing IDs and individual read IDs. One isolate may have several chromatograms, such as forward and reverse reads or multiple primer reads.

## Review and processing

A typical session is:

1. Load the chromatograms or ZIP archive.
2. Check the detected files and inferred sample/read assignments.
3. Confirm isolate IDs, primer names and forward/reverse directions.
4. Inspect chromatogram traces and sequence quality.
5. Review read overlaps, assembly placements and consensus evidence.
6. Inspect QC warnings and failures.
7. Correct metadata or review decisions where necessary.
8. Process the batch and retain the output reports.

Sanger Studio includes review and assembly visualisation features. The underlying assembly is overlap-based and uses the existing Paper Trail logic. Reads that do not meet the configured overlap criteria are reported as unmerged rather than being forced into the consensus.

## Sequence targets and primers

The interface supports target profiles such as 16S rRNA, ITS, 18S rRNA, 28S/LSU rRNA and custom targets. Profiles provide starting QC values; they are not universal biological acceptance thresholds. Review the settings for the assay and sequencing protocol in use.

Custom primer sequences can be supplied in the format described by the interface. For custom primers, assign the read direction explicitly when it cannot be inferred reliably from the primer name.

## Outputs

The exact output set depends on the selected workflow and options. It can include assembled FASTA sequences, per-read QC, assembly reports, review tables, chromatogram/assembly visualisations and a packaged batch result.

Retain the original AB1 files, sample map, run settings and generated reports together. These are needed to audit how a consensus was produced.

## Interpretation cautions

- A high-quality base call does not prove that the template is a pure culture or the expected organism.
- Mixed peaks, contamination, chimeras, poor coverage and primer-specific effects require biological review.
- A consensus sequence is only as reliable as the contributing reads and their assembly evidence.
- A sequence that fails QC should not be treated as a negative taxonomic result.
- Check all automatically inferred sample and primer metadata before processing a large batch.
