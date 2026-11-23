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

### Improved Sanger-Studio read review and resequencing support

Sanger-Studio now supports a more robust workflow for reviewing isolates with multiple sequencing attempts. Where resequencing has been performed, individual reads are retained and assessed separately, allowing alternative reads for the same primer to be compared before selecting the read to use for assembly. Only one read per primer is selected for assembly, while the original reads and their quality-control evidence remain available for review.

The multiple-read alignment interface provides a position-by-position comparison of selected reads against candidate sequences, with interactive chromatogram inspection to help assess disagreements. Chromatogram evidence is displayed in the appropriate nucleotide orientation, including reverse-direction reads, improving the interpretation of fluorescence signals when comparing raw sequencing data with aligned sequences.

These improvements help distinguish genuine sequence disagreements from orientation or coordinate-mapping issues and provide a more transparent basis for read selection and sequence assembly. Assembly candidates remain subject to inspection and user selection rather than being accepted automatically.

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

## Alternative and resequenced AB1 files

Partner primer-mapping workbooks may list more than one chromatogram for the same isolate and primer. These are treated as alternative sequencing reads for the same isolate, not as additional isolates.

The `primers_to_manifest` utility accepts the `Primer file mapping` worksheet format, including the standard primer columns and the generic `Other #1`–`Other #4` columns. For `Other` entries, it infers the primer label from the AB1 filename where possible. This includes common numeric bacterial 16S primers (for example, 27F, 63F, 515F, 785F, 806R, 907R and 1492R), as well as ITS, 18S, 28S and DS labels. Numeric primer names with separators such as `907-R` or `1492_R` are normalised to `907R` and `1492R`. The isolate ID in the workbook remains authoritative; it is not replaced by an identifier inferred from the filename.

Example:

```bash
branchmanager primers-to-manifest \\
  --input "primer_file_mapping.xlsx" \\
  --source ./all_AB1 \\
  --output ./sanger_input \\
  --manifest manifest.csv
```

A run log is written to `./sanger_input/primers_to_manifest.log` by default. It records the same progress messages, warnings and errors shown in the terminal. Use `--log /path/to/run.log` to choose a different log-file location.

Workbook file references may include submission subdirectories and may omit the `.ab1`/`.abi` extension. The resolver first checks the supplied relative path. For a qualified reference, it can accommodate an export-added prefix to the submission-folder name, but it will not fall back to a matching basename in an unrelated submission. Basename-only references must resolve uniquely. Ambiguous or missing references are reported rather than silently assigned.

The generated manifest contains one row per physical chromatogram. It retains:

- The stable isolate ID, read filename, primer and direction.
- An attempt identifier for repeated isolate/primer combinations.
- The source submission folder (`attempt_group`), original source path and filename.
- A SHA-256 checksum of the source trace.
- The workbook cell value, column and row used to map the file.

Repeated reads for the same isolate and primer receive unique filenames such as `SS172_907R_attempt01.ab1` and `SS172_907R_attempt02.ab1`. The suffix is a unique read label, not a claim about chronological order. The source folder is retained separately as provenance.

During assembly, reads with the same explicit primer are treated as alternatives: the best QC-passing read for that primer is retained for pairwise assembly, while reads from different primers remain eligible to provide complementary coverage. All original reads remain represented in the read-level QC and provenance outputs. This prevents two resequencing attempts of the same primer from being treated as independent overlapping primer reads.

Review the manifest and import warnings before processing. In particular, generic `Other` entries whose primer cannot be inferred retain an explicit `Other` label and should be checked manually. The workflow does not infer that a newer submission is necessarily better; read quality and assembly evidence determine selection.

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


### Resequencing read selection in Sanger Studio

When an isolate has multiple reads for the same primer, Sanger Studio treats them as alternative attempts rather than complementary reads. Review the per-read QC table and chromatogram evidence, then use **Select the read to retain for each primer** to choose one attempt per primer. The default is the highest-ranked QC-passing read (QC class, mean trimmed quality, then trimmed length). This choice does not alter the underlying QC classification.

Assembly candidates are generated only from distinct primers and are filtered to match the current per-primer selections. Reads from the same primer are never compared as assembly overlaps or combined in an assembly candidate. If the selected reads do not meet the configured overlap criteria, Studio reports that outcome rather than constructing a consensus from repeated reads of one primer.

## Primer selection, candidate taxonomy and chromatogram hover review

The primer-level read table allows one active read per primer. Unticking the active read excludes that primer from assembly; selecting another attempt replaces the previous attempt for that primer. A change to the read set invalidates any previous final-candidate choice. Candidate cards are filtered against the current read selections, and the displayed candidate must be explicitly selected before it can be downloaded or taxonomically screened.

Candidate taxonomy is a separate VSEARCH operation on the exact selected candidate sequence (either a multi-read assembly or an individual-read fallback). It uses the same conservative Paper Trail taxonomy-screen policy as the per-read screen. Component-read assignments are not copied onto the candidate. The candidate result records its candidate ID, sequence length, reference, assignment status, assigned rank, near-best support and reference hits in a candidate-specific output directory.

The chromatogram viewer supports cursor inspection. Hovering over an individual trace reports the nearest trace sample and the relative A/C/G/T channel signal at that point, alongside the called base and quality/secondary-peak evidence at the mapped base peak. Channel signal proportions are raw trace evidence, not calibrated posterior probabilities or a new base call. In peak-aligned overlay mode, the hover readout reports the corresponding signal support for each displayed read.

## Review-first cohort loading

When a cohort marker-review TSV is attached, Sanger Studio attempts to locate a cohort-level partner manifest in the associated cohort directory. It uses the review rows' `partner_id` and `dataset` values to select Sanger entries, then loads the listed sample maps. Relative AB1 paths are resolved first relative to the sample map and then against the manifest `read_dir`. If a manifest or matching map cannot be found, the review still loads and the status reports why reads were not automatically loaded. Standalone Sanger work continues to use manual folder, map, or AB1 selection.

### Retained candidate selection and marker-review persistence

When an inspected isolate has multiple sequence candidates, selecting a candidate in the inspection panel now synchronises that explicit retained choice into the attached cohort marker-review record. The saved TSV records `selected_candidate_id`, `selected_candidate_type` (`ASSEMBLY` or `INDIVIDUAL_READ`), `selected_candidate_length`, and `selected_candidate_source` (`MANUAL`). Candidate selection is validated against the most recent server-side inspection; candidate IDs supplied by the browser are not accepted unless they exist in that inspection. Manual retained-candidate choices require reviewer identity and update review provenance. The separate `selected_read_id` fallback field remains reserved for individual-read fallback decisions; when an explicit retained candidate is selected, the review table displays that selection instead of stale output/fallback text.
