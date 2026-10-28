**SYSTEM OVERVIEW & GUIDE**

**BranchManager**

Marker-Gene Quality Control, Novelty Scoring, & Genome-Candidate Prioritisation Toolkit

**Core Mission:** _Answers the critical question - "Has this microbial lineage already been isolated and characterised, or is it sufficiently novel and underrepresented to justify whole-genome sequencing (WGS) investment?"_

**1\. Executive Summary**

**BranchManager** is an end-to-end biological decision-support platform built for microbiology laboratories, research consortia, and biobanks managing large isolate collections. When processing thousands of physical isolates, determining which candidates warrant high-cost downstream sequencing is a major bottleneck.

BranchManager streamlines this process by evaluating marker-gene sequences (such as 16S rRNA), performing strict quality control on raw instrument chromatograms, classifying sequences across multiple global taxonomy databases, measuring multi-layered biological novelty, and proposing optimal primary and backup candidate sets for whole-genome sequencing.

| **🔬 Automated Quality Control**<br><br>Evaluates raw Sanger/AB1 chromatogram trace files, trims low-quality ends, detects mixed peaks, quality-weights overlapping reads, and isolates failed sequences for manual review. | **🌐 Multi-Lens Novelty**<br><br>Measures sequence divergence independently against internal physical lab baselines (The Hungate Collection), rolling multi-partner candidate collections, and global public references (GTDB, SILVA). | **🎯 Balanced Selection**<br><br>Maintains configurable pangenome targets per species, picking primary target candidates to fill genomic gaps while reserving backup strains for extraction failure resilience. |
| --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |

**2\. Core Operational Ethos & Principles**

**• Absence of a Hit ≠ Novelty:** A sequence lacking close reference hits is not automatically novel. It may represent poor trace quality, chimeric sequences, contamination, or incomplete coverage. BranchManager enforces strict QC filtering before novelty scoring.

**• Facts vs. Recommendations are Kept Distinct:** Factual states (e.g., physical DNA availability, confirmed genome presence, QC status) are maintained independently from recommendation states (e.g., Primary vs. Backup nomination). Recommendation flags never modify factual records automatically.

**• Fail-Closed Safety Principles:** Selection-critical processes fail closed. An incomplete tree build, unverified quality score, or broken novelty calculation will halt report publication rather than generating superficially valid recommendations.

**• Novelty is Relative to YOUR Project:** Novelty metrics specifically measure how unique an isolate is compared to what your project or lab has already characterised, helping to prevent redundant sequencing of strains already in your freezer.

**3\. The Workflow Stages ("Office Vocabulary")**

BranchManager uses an "office workflow" vocabulary to guide users through the multi-step candidate selection lifecycle from initial delivery to annual reporting:

| **Workflow Stage**      | **Scientific Purpose**                                                                                            | **Key Outputs & Deliverables**                                                           |
| ----------------------- | ----------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------- |
| **Mailroom**            | Inventories physical AB1 chromatogram deliveries, reconciles read IDs, and verifies primer metadata.              | Batch map (ab1_map.tsv), read inventory, missing/duplicate read reports.                 |
| **Interview**           | Standalone chromatogram QC and resequencing triage without project database modifications.                        | QC assembly report, resequencing recommendations, failed trace manifests.                |
| **Onboarding**          | Validates metadata ledgers, confirms isolate ownership, and enforces unique project ID mapping.                   | Normalised input datasets, partner ledger validation logs.                               |
| **Paper Trail & Merge** | Converts raw AB1/FASTQ traces, trims low-quality ends, masks low-confidence bases, and merges multi-primer reads. | Assembled FASTA, per-base error profiles, trace visual reports, manual review templates. |
| **Filing Cabinet**      | Registers cultured baseline datasets (e.g., Hungate) and constructs initial reference backbone trees.             | Baseline taxonomy tables, initial tree structures, baseline identity maps.               |
| **Performance Review**  | Classifies marker sequences, screens chimeras, scores multi-lens novelty, and updates phylogenetic trees.         | Sequence audit table, novelty summary, iTOL metadata strips, local clade figures.        |
| **Hiring Panel**        | Proposes group-level primary and backup isolate candidates based on current pangenome targets.                    | Selection summaries, primary/backup diversity plans, dashboard interface.                |
| **Quarterly Review**    | Reconsiders full collections across partner batches to select budgeted genome sequencing tranches.                | Tranche nomination lists, marginal distance audits, updated clade visual maps.           |
| **Status & Records**    | Imports factual lab progress (e.g., extraction success, completed WGS assemblies, ANI results).                   | Updated genome collections, factual isolate status logs, ANI quality gates.              |
| **Annual Report**       | Generates cumulative, audit-ready marker-to-genome close-out reports for institutional stakeholders.              | Executive HTML dashboard, complete historical audit trail, ledger close-outs.            |

**4\. Key System Tools & Operational Modules**

**Unified Golden-Path Runner (assistant)**

The primary standard execution route. It orchestrates Onboarding, Paper Trail quality control, Performance Review novelty calculations, and Hiring Panel candidate selection in an integrated sequence with zero manual intermediate handling.

**Chromatogram Quality Control Engine (mailroom, interview, paper-trail)**

Processes raw Sanger sequencing files. Performs Phred quality score evaluation, end-trimming, internal base masking, dye-channel peak analysis for mixed/heterozygous positions, and orientation-aware multi-read consensus assembly. Sequences failing QC are isolated in a dedicated triage folder with detailed diagnostic reason codes.

**Baseline & Reference Classifier (filing-cabinet)**

Loads physical baseline datasets (such as cultured strain biobanks) to anchor phylogenetic trees and establish the initial novelty threshold. Supports multi-database taxonomic reporting (GTDB, SILVA, GG2, NCBI).

**Novelty Scoring & Selection Engine (performance-review, quarterly-review)**

Evaluates incoming candidate datasets against the accumulating project database. Calculates neighbourhood crowding metrics, identifies baseline-redundant strains (≥99.8% identity across ≥95% alignment), and applies farthest-first phylogenetic branch distance calculations to pick optimal isolate subsets.

**Focused Taxonomy Slicing & Tree Styling (org-chart, label-maker)**

Allows researchers to instantly extract specific taxons (e.g., domain Archaea or family Lachnospiraceae) from existing alignments to generate high-resolution focused phylogenetic trees and customise interactive iTOL visual metadata strips without recomputing alignments.

**5\. Candidate Selection & Decision Categories**

BranchManager replaces subjective manual picking with a transparent, evidence-based selection model. Decisions are presented clearly for expert-human review:

| **Decision Category**            | **Status Type**        | **Scientific Rationale & Use Case**                                                                                                       |
| -------------------------------- | ---------------------- | ----------------------------------------------------------------------------------------------------------------------------------------- |
| **PRIORITISE - SET PRIMARY**     | **Primary Nomination** | High-value candidate filling an unrepresented species gap or missing pangenome slot. Highest priority for WGS.                            |
| **RESERVE - SET BACKUP**         | **Backup Candidate**   | Phylogenetically diverse strain retained as a secondary option in case the primary candidate fails DNA extraction or library preparation. |
| **EXCLUDE - BASELINE REDUNDANT** | **Excluded Strain**    | Strain sharing ≥99.8% 16S identity over ≥95% coverage with an already-cultured baseline strain. Omitted to prevent redundant sequencing.  |
| **REVIEW - PANGENOME BOUNDARY**  | **Manual Review**      | Sequence falls near species identity boundaries; requires manual review to confirm if it represents a distinct novel lineage.             |
| **LOWER PRIORITY - TARGET MET**  | **Lower Priority**     | The maximum target number of genomes for this exact GTDB species has already been satisfied by existing baseline or pending genomes.      |

**6\. Core Benefits & Strategic Value**

**• Sequencing Budget Optimisation:** Prevents spending WGS resources on isolates that duplicate strains already characterised in internal or public biobanks.

**• Multi-Partner Harmonisation:** Provides a single source of truth database across multiple contributing laboratories, ensuring consistent ID mapping and metadata tracking.

**• Maximum Diversity Capture:** Utilises branch-length algorithms to pick candidates that maximise phylogenetic spread across target clades.

**• Full Auditability & Reproducibility:** Every candidate recommendation is linked to raw chromatograms, specific quality metrics, component novelty scores, and paginated visual trace reports.

**• Scalable Tranche Management:** Supports phased sequencing rounds (e.g., quarterly budgets), automatically updating recommendations as newly sequenced genomes become available.


## Quick start

### Collection profiles

Use a JSON project profile to keep collection-specific names, reference databases, baseline datasets, anchor FASTA, domain, and selection policy together. A profile supplies defaults; an explicit command-line option always takes precedence. Copy [`configs/project.template.json`](configs/project.template.json) and edit it for a collection. Set the optional top-level `project_directory` to the working directory for that project; relative data and CLI-default paths are resolved from there. If omitted, paths remain relative to the profile file's directory. Environment substitutions use `${NAME}` or `${NAME:-fallback}`, allowing local paths to stay out of shared configuration files.

For an Assistant-only setup, copy [`configs/assistant.template.json`](configs/assistant.template.json). It puts the commonly reused database, partner metadata, marker reference(s), baselines, MWL, anchors, and workflow defaults under `defaults.assistant`. Then a batch run only needs its read map (or FASTA/PacBio map), dataset name, and output directory:

```bash
branchmanager assistant --config configs/my_collection_assistant.json \
  --sample-map batch_map.tsv --dataset Collection_01 -o runs/Collection_01
```

CLI flags can still override any configured default. Keep the per-batch input and dataset on the command line so each run remains explicit. The broader `project.template.json` remains available for collections that also want defaults for other BranchManager commands.

In Assistant profiles, list each cross-check reference as one object in `defaults.assistant.alt_databases`, with `name`, `ref`, and optional `taxa` (use `null` when taxonomy is embedded in the FASTA headers). This keeps each taxonomy table attached to its reference and display name. Do not combine `alt_databases` with the older parallel `alt_ref`, `alt_taxa`, and `alt_ref_name` lists in the same Assistant profile.

The checked-in [`configs/rumen_gateway.json`](configs/rumen_gateway.json) is the Rumen Gateway project profile. It holds the project database and partner catalogue paths, GTDB r232 plus GG2/NCBI/SILVA cross-check databases, the Hungate priority baseline and the Rumen Gateway `Secondary` baseline, the MWL, and selection defaults. Paths resolve relative to the configured Rumen Gateway `project_directory`. It omits a custom anchor path so BranchManager uses its bundled 26-sequence rumen/gut anchors by default. Use the profile once with Filing Cabinet to initialize both baselines, then with Assistant for each partner batch:

```bash
branchmanager assistant --config configs/rumen_gateway.json \
  --sample-map batch_map.tsv --dataset QUB_01 -o runs/QUB_01
```

The profile does not contain a batch's sample map, dataset label, or output directory. Supply those per run; command-line values override configured defaults. The `baselines` list is used by the Filing Cabinet command and is what registers the priority and secondary sequence pools in the project database.

When several partner mailrooms and metadata workbooks are ready, run one joint assessment with an Assistant cohort manifest. Mailroom remains a per-partner operation. The manifest is a CSV or TSV with one row per partner batch; paths are relative to the manifest file unless absolute. `partner_workbook` is optional when that workbook has already been imported. Each batch is validated and marker-QC'd separately, then all passing markers are assessed together in one selection run, so partner order does not affect the candidate pool. QC failures are retained as markerless isolates and reported for later replacement.

```tsv
partner_id	dataset	sample_map	read_dir	partner_workbook
QUB	QUB_01	mailroom/QUB_01/ab1_map.tsv	reads/QUB_01	metadata/QUB_01.xlsm
UoG	UoG_02	mailroom/UoG_02/ab1_map.tsv	reads/UoG_02	metadata/UoG_02.xlsm
```

```bash
branchmanager assistant --config configs/rumen_gateway.json \
  --cohort-manifest cohort.tsv --dataset RGW_COHORT_2026_10 \
  -o runs/RGW_COHORT_2026_10
```

The configured `--partner-metadata` remains the cumulative isolate catalogue used to validate ownership and prior sequencing commitments. Per-partner workbook rows in the manifest add workbook metadata to the project database before assessment. Give each partner batch a unique dataset label. The cohort run writes per-batch onboarding and Paper Trail reports under `00_cohort_preparation/`, a combined `cohort_marker_review_template.tsv` there, and the single combined selection report under `01_joint_assistant/`. To record manual decisions, copy the combined template to a separate file, fill `decision`, `reviewer`, and `notes` for reviewed isolates (`approve` or `reject`), then rerun Assistant with `--marker-review reviewed_cohort.tsv`. Do not edit the generated template in place because a rerun regenerates it.

The profile's `project` object is collection metadata and is copied into run manifests; `data` documents relevant files and, for Assistant, supplies run defaults for `database`, `dataset`, `sample_map`, and `output_directory`; `defaults` maps command names to CLI argument names. Explicit CLI values override profile defaults. Each run manifest records the profile file checksum. Additional columns in the partner ledger are retained per isolate in the database's `sequencing_metadata.project_metadata_json` field, while BranchManager interprets its standard ID, partner, selection, and genome-availability columns. The bundled rumen/gut anchor set is used unless a project config or command-line option supplies a custom set; projects outside that biological scope should configure vetted, domain-appropriate anchors. See [collection profiles and scope](docs/collection-profiles-and-validation.md) before interpreting results from a new marker/domain.

Profiles can declare multiple baselines under `baselines`, each with a unique `dataset`, `fasta`, `tier` (`priority` or `secondary`), and optional `taxonomy_assignments` or `classify`. A single Filing Cabinet invocation loads each baseline into its own dataset and tier, then builds one combined backbone from the available marker sequences:

```bash
branchmanager filing-cabinet --config configs/rumen_gateway.json \
  --db project.sqlite --build-tree -o filing_cabinet
```

The profile's `project.selection_budget` supplies the default number of primary genomes nominated in a Quarterly Review; override it with `--genome-budget`. Backups are controlled separately. In a validation table, precision@k is the fraction of the k selected candidates that meet the prespecified success outcome, while recall@k is the fraction of all successful candidates recovered among those k. These metrics should use the same k as the real sequencing budget.

For a complete partner submission, use the guarded golden path:

```bash
branchmanager assistant \
  --sample-map ab1_mapping.tsv \
  --partner-metadata partner_metadata.tsv \
  --db project.sqlite --dataset QUB_01 \
  --ref gtdb_ssu_reps.fna.gz --ref-name GTDB \
  --chimera-ref silva_gold_ssu.fasta.gz \
  --baseline-fasta hungate.fasta --baseline-dataset Hungate \
  --mwl MWL.csv --threads 10 \
  -o runs/01_QUB_01
```

`assistant` runs IT Desk checks, Onboarding, Paper Trail/Merge Meeting, Performance Review, and Hiring Panel in numbered directories. Each stage remains independently runnable:

```bash
# 1. Establish the cultured baseline Filing Cabinet and build a reference tree
branchmanager filing-cabinet \
  --fasta hungate.fasta \
  --db project.db \
  --dataset Hungate \
  --ref gtdb_r232_bac_arch_ssu_reps.fna \
  --classify --build-tree \
  -o filing_cabinet_out \
  --threads 10

# 2. Process new sequences — scores novelty against the baseline
branchmanager performance-review \
  --input new_sequences.fasta \
  --db project.db \
  --dataset Batch1 \
  --ref gtdb_r232_bac_arch_ssu_reps.fna \
  --previous-review filing_cabinet_out \
  -o batch1_out \
  --threads 10

# 3. Zoom into a specific taxon
branchmanager org-chart \
  --db project.db \
  --taxon archaea \
  --from-dir filing_cabinet_out \
  -o archaea_out

# 4. Regenerate iTOL files with new grouping options
branchmanager label-maker \
  --db project.db \
  --out filing_cabinet_out \
  --group-phyla archaea \
  --group-phyla "Bacillota,Bacillota_I,Bacillota_A"
```

---

## Core concepts

### Novelty is relative to YOUR submitted sequences

All novelty and neighbourhood-density calculations are made against the sequences submitted to BranchManager (the Filing Cabinet baseline plus all prior Performance Review datasets stored in the DB). They are **not** made against the full external reference (GTDB/SILVA).

This is intentional: BranchManager helps you judge whether a new sequence is worth investigating relative to what your lab or project has already characterised — not relative to all known biology.

Each successive Performance Review extends the reference pool, so novelty scores become increasingly precise as your project grows.

### Three-layer tree architecture

The phylogenetic tree is built from three independent layers:

| Layer | Source | Shown in iTOL? |
|---|---|---|
| 1 — Anchors | 26 bundled NCBI RefSeq rumen/gut anchor sequences, unless a custom `--anchors` FASTA is supplied | **No** — invisible topology scaffolding |
| 2 — Filing Cabinet | Your cultured baseline dataset (e.g. Hungate) | Yes |
| 3 — Performance Review | Each new partner marker batch | Yes |

Anchor sequences constrain phylum-level topology during MAFFT + FastTree inference, then are pruned from the stored newick. They are never shown in outputs, iTOL files, or novelty scoring.

---

## Subcommands

### `branchmanager filing-cabinet` 

Load a baseline FASTA dataset, classify against a reference, and build the backbone tree. This is always the first step.

```bash
branchmanager filing-cabinet \
  --fasta baseline.fasta \
  --db project.db \
  --dataset Hungate \
  --ref gtdb_r232_bac_arch_ssu_reps.fna \
  --classify --build-tree \
  -o filing_cabinet_out \
  --threads 10 \
  --collapse
```

| Parameter | Required | Default | Description |
|---|---|---|---|
| `--fasta` | ✓ | — | Input FASTA file containing the baseline sequences |
| `--db` | ✓ | — | Path to the BranchManager SQLite database (created if absent) |
| `-o / --out` | | `.` | Output directory for tree, iTOL files, and reports |
| `--dataset` | ✓ | — | Label stored in the DB (e.g. `Hungate`) |
| `--ref` | | — | Reference FASTA (GTDB/SILVA reps) for classification and tree orientation |
| `--classify` | | off | Classify sequences against `--ref` and store taxonomy |
| `--build-tree` | | off | Build the MAFFT + FastTree backbone tree |
| `--taxa` | | — | TSV mapping reference IDs to lineages (optional when `--ref` headers contain GTDB lineages) |
| `--taxa-assignments` | | — | Pre-computed taxonomy for the INPUT sequences (TSV or embedded-lineage FASTA) |
| `--shorten-ids / --no-shorten-ids` | | `--no-shorten-ids` | Preserve input headers exactly by default; use `--shorten-ids` only when compact generated IDs are desired |
| `--collapse` | | off | Collapse near-identical same-taxonomy sequences into one representative for the tree |
| `--collapse-threshold` | | `99.8` | Identity threshold (%) for collapsing |
| `--sequence-domain` | | `bacteria` | Domain/profile to process: `bacteria`, `archaea`, `fungi`, or `mixed`/`all`/`none` to disable filtering. Use matching references, previous reviews, and anchors for non-bacterial runs |
| `--anchors` | | bundled | Custom anchor FASTA for tree topology scaffolding |
| `--threads` | | `4` | CPU threads for MAFFT and VSEARCH |
| `--colours` | | — | CSV mapping sequence IDs to custom hex colours for iTOL (columns: `id`, `colour`) |
| `--group-phyla SPEC` | | — | Group phyla into a single colour in iTOL (repeatable; see *Phylum grouping*) |
| `--functional` | | — | TSV file mapping sequence IDs to functional attributes (first column = ID; subsequent columns = attributes). Generates one iTOL file per column (DATASET_BINARY / DATASET_SIMPLEBAR / DATASET_COLORSTRIP). |

---

### `branchmanager mailroom`

Prepare the immutable per-batch AB1 map before Onboarding. Mailroom scans every physical chromatogram, extracts ABIF sample/run/instrument and base-call evidence, reconciles supplier sequencing IDs, and reports missing, duplicated, ambiguous, or unmapped reads.

```bash
branchmanager mailroom \
  --read-dir UoG/All_AB1 \
  --metadata UoG/supplier_metadata.csv \
  --dataset UoG_01 \
  --forward-primer 63F \
  --reverse-primer 1492R \
  -o UoG/UoG_01
```

Only use `--forward-primer` and `--reverse-primer` after the oligonucleotides have been confirmed from the partner protocol or sequencing submission. Supplier metadata may instead contain a `primer` column. Mailroom uses embedded ABIF primer names when present, but never treats sequence-motif similarity as proof of the exact primer.

Preferred supplier-metadata columns:

```tsv
sequencing_id	isolate_number	read	primer	processing_mode
KKX994	SW_0016	Forward	63F	assemble
KKY011	SW_0016	Reverse	1492R	assemble
```

`primer` and `processing_mode` are optional. With `--processing-mode auto`, paired/multi-read isolates use `assemble` and single-read isolates use `best_read`.

| Output | Purpose |
|---|---|
| `ab1_map.tsv` | Validated one-row-per-read input for Onboarding |
| `ab1_inventory.tsv` | Physical file and extracted ABIF metadata inventory |
| `mailroom_report.tsv` | Missing, duplicate, ambiguous, unreadable, unmapped, and unresolved-primer findings |
| `mailroom_summary.json` | Machine-readable counts and `PASS`, `REVIEW_REQUIRED`, or `FAIL` status |
| `run_manifest.json` | Checksums, inputs, outputs, and stage status |

Mailroom returns `REVIEW_REQUIRED` when files reconcile but primer names remain unresolved, and `FAIL` for structural/file errors. Correct the supplier metadata or confirm primers before Onboarding.

For PacBio 16S, retain the same Mailroom entry point but select the technology explicitly. PacBio Mailroom validates one segmented HiFi FASTQ per isolate, reconciles it with the supplied `sequenceid`/`sequence_id` and `read_file`/`fastq_file` columns, writes `pacbio_map.tsv`, and automatically runs clustering and representative selection under `pacbio_16s_import/`. Sanger remains the default.

```bash
branchmanager mailroom \
  --technology pacbio \
  --read-dir UKNZ_01/All_Pacbio \
  --metadata UKNZ_01/UKNZ_01_pacbio_metadata.csv \
  --dataset UKNZ_01 \
  -o UKNZ_01/mailroom
```

---

### `branchmanager interview`

Run AB1 conversion and QC without Onboarding, taxonomy, tree building, candidate selection, or any project-database access. Interview consumes the complete Mailroom output directory and uses the same versioned QC policy as `assistant`.

```bash
branchmanager interview \
  --mailroom UoG/UoG_01 \
  --screen-ref gtdb_r232_bac_arch_ssu_reps.fna.gz \
  --threads 10 \
  -o UoG/UoG_01_interview
```

`--screen-ref` is optional and only checks whether separate primer reads from the same isolate have discordant taxonomic hits. A disagreement produces `PASS_WITH_WARNINGS / MANUAL_REVIEW` unless another QC criterion fails. The common comparison rank and each `read_id=taxon` assignment are written to the concordance, assembly, recommendation, and review tables. Without it, all chromatogram, trimming, assembly, quality, ambiguity, mixed-peak, and resequencing outputs are still produced.

Interview accepts either the Mailroom directory or its `ab1_map.tsv`. It requires the accompanying `mailroom_summary.json`, rejects a Mailroom `FAIL`, and stops on `REVIEW_REQUIRED` unless the batch has been reviewed and `--allow-mailroom-review` is supplied explicitly.

The principal hand-off files are `assembled.fasta`, `assembly_report.tsv`, `marker_review_template.tsv`, and `resequence_recommendations.tsv`. Once any manual-review decisions are completed, these can enter the full workflow without repeating chromatogram processing:

```bash
branchmanager assistant \
  --fasta UoG/UoG_01_interview/assembled.fasta \
  --marker-qc UoG/UoG_01_interview/assembly_report.tsv \
  --marker-review UoG/UoG_01_interview/marker_review_template.tsv \
  --partner-metadata Sequence_Metadata.csv \
  --db project.sqlite --dataset UoG_01 --ref gtdb_ssu_reps.fna.gz \
  -o runs/UoG_01
```

Do not supply the blank review template to Assistant. Add decisions only for the listed `PASS_WITH_WARNINGS` markers first; omit `--marker-review` when the template contains no rows.

---

### `branchmanager onboarding`

Every partner dataset is onboarded from exactly one marker source. AB1 submissions use a sample map; FASTA-only submissions use the FASTA directly. Both paths validate sequence IDs against the same cumulative project metadata ledger. The ledger can also record `addition_method` as `marker_sequence`, `maldi_tof`, or `both`; leave it blank when unknown. MALDI-only isolates retain a stable isolate ID and an empty marker slot, so a passing marker sequence can be attached to that same ID in a later batch. A failed marker submission remains as a markerless QC-failed record without classifier taxonomy when no accepted marker exists. If an accepted marker already exists, a failed replacement attempt is recorded while the accepted marker and assignment are retained.

```bash
# AB1 or multiple primer reads
branchmanager onboarding \
  --sample-map UoG/UoG_01/ab1_map.tsv \
  --read-dir UoG_01_AB1/ \
  --partner-metadata project_partner_metadata.tsv \
  --partner-id UoG --dataset UoG_01 \
  -o UoG_01/01_onboarding

# Partner-supplied marker FASTA; .gz is accepted
branchmanager onboarding \
  --fasta QUB_01_16S.fasta.gz \
  --partner-metadata project_partner_metadata.tsv \
  --partner-id QUB --dataset QUB_01 \
  -o QUB_01/01_onboarding
```

FASTA Onboarding writes `normalised_input.fasta`. It does not invent chromatogram quality evidence. Supply a reviewed `--marker-qc` table to Performance Review where one exists, or use `--accept-unverified-marker-qc` as an explicit audited acceptance; unverified marker evidence remains visible in the reports.

### PacBio 16S partner submissions

Use PacBio 16S only after SMRT Link has segmented/de-concatenated and demultiplexed the library molecules. One FASTQ is supplied per cultured isolate, containing many observations of the approximately 1.5 kb full-length 16S amplicon. BranchManager filters reads by length and mean quality, dereplicates them, clusters the unique sequences with vsearch at 99.5% identity, and selects the highest-abundance observed sequence from the dominant cluster. An isolate passes only when it has at least 20 passing reads and the dominant cluster contains at least 80% of them; adjust these defaults only with documented validation. The current default library provenance is `Kinnex 16S rRNA`; it is recorded in the marker QC reason field alongside the SMRT Link segmentation workflow.

```tsv
sequence_id	fastq_file
Iso001	fastq/Iso001.fastq.gz
Iso002	fastq/Iso002.fastq.gz
```

```bash
branchmanager pacbio-16s-import --pacbio-map QUB_PACBIO_01/pacbio_map.tsv -o QUB_PACBIO_01/00_pacbio_16s

# Or run it as the partner-sequenced input to the full workflow:
branchmanager assistant --pacbio-map QUB_PACBIO_01/pacbio_map.tsv \
  --partner-metadata project_partner_metadata.tsv --partner-id QUB --dataset QUB_PACBIO_01 \
  --db project.db --ref gtdb_16s.fasta -o QUB_PACBIO_01/run
```

PacBio 16S Import writes `pacbio_16s_representatives.fasta`, `marker_qc.tsv`, `read_qc.tsv`, and `pacbio_16s_qc.tsv`. Failed or mixed isolates are excluded from the marker FASTA and remain visible in the QC report; `read_qc.tsv` retains the source FASTQ link in the downstream provenance, while raw FASTQs and the per-isolate clustering work files remain available for audit. The PacBio QC sidecar is consumed automatically by `assistant`, so it is not accepted as unverified partner FASTA.

Partner acronyms and submission labels are separate. For example, `UoG_01` and `UoG_02` are distinct `--dataset` values, while both sets of ledger rows use `partner_id=UoG`. Sequence IDs must be unique and stable across the entire project.

Keep exactly one cumulative metadata ledger with one row per durable isolate/marker ID. It owns `partner_id`, `selected_for_genome_sequencing`, and `already_sequenced`; chromatogram run IDs do not belong in it. Keep a separate read map under each partner/batch, for example `UoG/UoG_01/ab1_mapping.csv`, which maps every physical trace filename to its durable isolate ID, direction, and `processing_mode`. Duplicate ledger IDs, partner mismatches, missing traces, and traces assigned to more than one isolate are Onboarding errors.

Onboarding reports only actionable errors and warnings. Normalised submission/read-map tables omit optional columns that are empty throughout the submitted batch.

---

### `branchmanager paper-trail` (AB1 conversion and Merge Meeting)

Convert Sanger AB1 chromatograms or already-basecalled primer reads into a trimmed FASTA for the Performance Review. When multiple reads belong to the same isolate, for example `27F` and `907R`, BranchManager orients reverse-primer reads and builds a quality-aware overlap consensus.

#### Interactive Sanger Studio

Use the local browser-based Sanger Studio to inspect raw reads and assembly decisions before submitting a batch. It calls the same Paper Trail base-caller, trimming, QC, and assembler used by the command line. It displays read-level QC, trimmed positions, primer orientation, overlap/consensus evidence, warnings, and the generated visual reports. Primer and direction fields can be corrected for an isolate and re-run; read mode can use the map setting, force assembly, or select the best individual read. A passing sequence with warnings can be recorded as approved or rejected for later marker review.

```bash
branchmanager sanger-studio \
  --root /path/to/partner_submission \
  --sample-map /path/to/partner_submission/mailroom/ab1_map.tsv \
  -o /path/to/partner_submission/sanger_studio_runs
```

The root must contain both the map and every read path referenced by it. Set `--root` to a common parent if reads are stored in a separate directory. If no map is supplied, use **Browse** to load one, or select raw AB1/ABI files to create a quick editable map. The GUI binds only to localhost; use `--no-browser` to start it without opening a browser automatically, and `--port` to choose another local port. **Process and download all** runs the loaded map through Paper Trail once and downloads a ZIP containing the assembled FASTA, `assembly_report.tsv`, read-level reports, visualizations, failed-QC records, resolved read map, run manifest, and any saved `marker_review_decisions.tsv`.

Sanger Studio is an optional review/preparation step; it does not read or change the project database and it never modifies or replaces an Assistant run. It shares the same Paper Trail engine and QC defaults. There are two supported ways to continue into the main workflow:

- For Assistant to process the raw reads itself as part of the normal end-to-end run, give it the Mailroom map with `branchmanager assistant --sample-map ...`. Assistant runs Paper Trail/Merge Meeting again under its own run directory and manifest.
- To use the exact sequence set already reviewed in Studio, extract the downloaded ZIP and give Assistant the included `assembled.fasta` with `--fasta`, its adjacent `assembly_report.tsv` with `--marker-qc`, and the Studio `marker_review_decisions.tsv` with `--marker-review` (if review decisions were recorded). Assistant then starts at marker provenance and taxonomy assessment, rather than reprocessing the traces.

These are separate runs. Choose the first route when you want the standard fully auditable Assistant run from raw traces; choose the second when Studio's inspected assembly is the intended hand-off. QC defaults are shared, but changing Studio-only settings does not silently change a later Assistant run.

For a joint `assistant --cohort-manifest` run, each cohort entry currently expects a partner sample map and Assistant performs Paper Trail itself. Studio-produced FASTAs are accepted by the single-partner `assistant --fasta` route, but are not yet a direct input type for cohort-manifest entries. To preserve the cohort's joint selection behavior, use Studio as a preflight review and then provide the original partner maps to the cohort Assistant.

The default policy allows an individual usable primer read of at least 300 bp while requiring the final marker sequence or assembled consensus to be at least 800 bp. Therefore, a good 700 bp read can contribute to a passing consensus if it overlaps its partner read; the 800 bp threshold applies to the resulting marker, not each primer read. A low-quality read excluded from assembly cannot contribute length. Passing with warnings means the evidence needs review; it does not automatically mean resequencing is necessary. Failed QC results include the read-level or assembly reasons to help distinguish likely poor culture/purity from trimming, mapping, or coverage issues.

```bash
branchmanager paper-trail \
  --input ab1_reads/ \
  -o paper_trail_out \
  --min-quality 20 \
  --min-mean-quality 20 \
  --min-length 800 \
  --min-overlap 80 \
  --min-overlap-identity 0.98

branchmanager performance-review \
  --input paper_trail_out/assembled.fasta \
  --partner-metadata new_sequences_metadata.tsv \
  --db project.db \
  --dataset Batch1 \
  --ref gtdb_r232_bac_arch_ssu_reps.fna \
  -o batch1_out
```

If filenames are structured like `Iso001_27F.ab1` and `Iso001_907R.ab1`, BranchManager can infer the read information. Production partner submissions should nevertheless use one explicit `--sample-map` per batch, preferably with one row per physical read:

```tsv
sequence_id	dataset	read_file	primer	direction	processing_mode
Iso001	UoG_01	well_A01.ab1	27F	forward	assemble
Iso001	UoG_01	well_A02.ab1	907R	reverse	assemble
Iso002	UoG_01	well_B01.ab1	27F	forward	best_read
```

Segmented PacBio HiFi 16S imports use `vsearch --cluster_fast` inside kinnex. `branchmanager pacbio-16s-import` exposes `--threads`, `--vsearch-strand`, `--vsearch-query-cov`, `--vsearch-maxaccepts`, and `--vsearch-maxrejects`, and `branchmanager mailroom --technology pacbio` forwards the same controls with `--pacbio-...` prefixes when it automatically runs kinnex after validating the FASTQ map.

Multiple rows may share the same `sequence_id` when an isolate has multiple primer reads. Relative paths are resolved from the batch-map location. A wide one-row-per-isolate representation is also accepted when laboratories naturally supply separate primer columns:

```tsv
isolate_id	27F	907R	processing_mode
Iso001	well_A01.ab1	well_A02.ab1	assemble
Iso002	well_B01.ab1	well_B02.ab1	best_read
```

Run it with:

```bash
branchmanager paper-trail --sample-map sample_reads.tsv -o paper_trail_out
```

Onboarding writes `normalised_read_map.tsv` as the validated hand-off to Paper Trail. It is a generated run artefact, not another source ledger and not a file users should maintain.

`assemble`/`merge`/`consensus` orients the primer reads and tries to build one longer sequence. `best_read`/`highest_quality`/`select_best`/`independent` converts each read independently and writes only the highest-quality passing read for that isolate.

BranchManager applies the same versioned Sanger QC policy in `interview`, `paper-trail`, and `assistant`:

- `--min-quality` is a Phred cutoff used for Mott-style end trimming.
- `--mask-quality` masks internal bases below Q20 to `N` before assembly; masking is reported but does not by itself force manual review.
- Primer sequences are removed only after a confident IUPAC-aware leading match; customise them with repeatable `--primer-sequence NAME=SEQUENCE`.
- Mixed chromatogram peaks are measured from `PLOC` and four `DATA9`-`DATA12` dye channels. `--secondary-peak-ratio` and `--max-mixed-peak-percent` control the retained-region gate.
- Overlap consensus bases use posterior base probabilities from both Phred observations; qualities are not added as if they were independent certainty scores.
- `--screen-ref` optionally classifies each primer read independently. Family/genus discordance triggers explicit manual review and reports the differing read-level assignments; it is not a QC failure by itself.
- Final sequences must be at least 800 bp. Individual primer reads may be at least 300 bp so complementary reads can assemble into a valid final marker.
- Ambiguous bases use a tiered rule: more than 3% requires manual review and more than 5% fails QC by default.
- `--min-mean-quality`, expected-error limits, ambiguity, internal low-quality runs, mixed peaks, and overlap conflict density determine `PASS_HIGH_CONFIDENCE`, `PASS_WITH_WARNINGS`, or `FAIL_QC`.
- Final `FAIL_QC` sequences are withheld from `assembled.fasta` and listed as `RESEQUENCE`.
- Final `FAIL_QC` sequences and failed reads are retained under `failed_qc_sequences/` for manual review.
- `PASS_WITH_WARNINGS` sequences are included but listed as `MANUAL_REVIEW`.
- Visual reports are automatically split into numbered PNG pages before they exceed 2400 pixels high. Adjust the ceiling with `--max-report-image-height` (minimum 600 pixels).

The defaults follow Phred error probabilities and published Sanger guidance that Q20 is a reliable base call and that, as a general rule, trimmed sequence should contain less than 5% ambiguous bases. The 3% review boundary is deliberately more conservative for phylogenetic candidate selection. See [Ewing and Green 1998](https://doi.org/10.1101/gr.8.3.186) and [Crossley et al. 2020](https://doi.org/10.1177/1040638720905833).

Outputs:

| File | Description |
|---|---|
| `assembled.fasta` | One trimmed/assembled sequence per isolate, suitable for `branchmanager performance-review --input` |
| `failed_qc_sequences/failed_final_sequences.fasta` | Final isolate consensus/best-read sequences that failed final QC |
| `failed_qc_sequences/failed_read_sequences.fasta` | Individual read sequences that failed read-level QC, oriented when direction is known |
| `failed_qc_sequences/failed_qc_manifest.tsv` | One row per failed isolate, combining the final outcome with read IDs, source files, read-failure evidence, and the best recoverable read |
| `failed_qc_sequences/failed_read_manifest.tsv` | One row per failed physical read with trimming and QC evidence; kept separate to avoid duplicate isolate rows |
| `failed_qc_sequences/failed_qc_guide.txt` | Plain-language explanation of failure reason codes and the distinction between isolate- and read-level records |
| `trimmed_oriented_reads.fasta` | Individual reads after quality trimming and primer-direction orientation |
| `raw_reads.fasta` | Raw base calls extracted from AB1 or input sequence files |
| `read_qc.tsv` | Per-read trimming, quality, expected error, length, and filter status |
| `per_base_error.tsv` | Per-base quality/error probability table with left/right trim and retained-base status |
| `visual_reports/read_error_profiles/read_error_profiles_page_*.png` | Paginated per-read quality/error profiles showing retained trim windows |
| `visual_reports/trace_chromatograms/trace_chromatograms_page_*.png` | Paginated dye-channel chromatograms with retained windows and mixed-peak evidence |
| `read_taxonomy_concordance.tsv` | Independent primer-read classifications, common comparison rank, `read_id=taxon` assignments, and within-isolate agreement when `--screen-ref` is used |
| `assembly_report.tsv` | Per-isolate assembly status, overlap identity, conflicts, unmerged reads, contributing read IDs, and any read-level taxonomy disagreement with the compared rank and assignments |
| `assembly_read_placements.tsv` | Per-read consensus coordinates, aligned blocks, and whether the read contributed to the selected consensus |
| `resequence_recommendations.tsv` | Per-isolate `ACCEPT`, `MANUAL_REVIEW`, or `RESEQUENCE` decision with reason codes, taxonomy-concordance evidence, and suggested action |
| `marker_review_template.tsv` | Pre-filled rows for manual-review markers, including conflicting read-level taxonomy assignments; complete `decision`, `reviewer`, and `notes`, then provide it as `--marker-review` |
| `paper_trail_qc_policy.tsv` | Versioned thresholds actually used for the run, for reproducibility |
| `visual_reports/assembly_overviews/assembly_overview_page_*.png` | Paginated per-isolate overviews with reads placed on consensus coordinates, contribution status, and assembly diagnostics |
| `visual_reports/visual_report_manifest.tsv` | Page index with report type, record range, dimensions, and continuation notes |
| `paper_trail_summary.txt` | Short run summary |
| `run_manifest.json` / `.tsv` | Input checksums, software/tool versions, thresholds, stages, outputs, and final status |

Supported inputs are `.ab1`, `.abi`, `.fasta`, `.fa`, `.fna`, `.fastq`, and `.fq`, with `.gz` accepted for AB1/ABI, FASTA, and FASTQ files.

---

### `branchmanager performance-review`

Process new sequences against the baseline; score novelty and update the phylogenetic tree.

```bash
branchmanager performance-review \
  --input new_sequences.fasta \
  --db project.db \
  --dataset Batch1 \
  --ref gtdb_r232_bac_arch_ssu_reps.fna \
  --previous-review filing_cabinet_out \
  -o batch1_out \
  --threads 10 \
  --collapse
```

| Parameter | Required | Default | Description |
|---|---|---|---|
| `--input` | ✓ | — | FASTA file of new sequences to analyse |
| `--db` | ✓ | — | Path to the BranchManager SQLite database |
| `-o / --out` | ✓ | — | Output directory (assessment TSV, novelty metrics, tree, iTOL files) |
| `--dataset` | ✓ | — | Label for this batch (used in iTOL dataset-membership strip) |
| `--ref` | | — | Reference FASTA for classification and tree orientation |
| `--previous-review` | | — | Previous Filing Cabinet or Performance Review output directory used to seed an incremental MSA. When the same `-o` directory is reused, `tree/current_alignment.fasta` is detected automatically |
| `--taxa` | | — | TSV mapping reference IDs to lineages |
| `--taxa-assignments` | | — | Pre-computed taxonomy for the input sequences |
| `--partner-metadata / --sequencing-metadata` | Performance Review | — | Cumulative CSV/TSV ledger with sequence IDs, partner acronyms, optional `selected_for_genome_sequencing`, and required `already_sequenced` status |
| `--shorten-ids / --no-shorten-ids` | | `--no-shorten-ids` | Preserve input headers exactly by default; use `--shorten-ids` only when compact generated IDs are desired |
| `--min-len` | | `800` | Minimum sequence length to retain (bp) |
| `--max-n-percent` | | `5.0` | Maximum percentage of ambiguous (N) bases allowed; matches Paper Trail/Interview final-sequence QC |
| `--marker-qc` | | auto | Paper Trail/Merge Meeting `assembly_report.tsv`; auto-discovered beside the FASTA |
| `--marker-review` | | — | Reviewed accept/reject decisions for warning or unverified marker evidence |
| `--accept-unverified-marker-qc` | | off | Explicit, audited acceptance of independently validated FASTA without Paper Trail provenance |
| `--chimera-ref` | | primary ref | Curated chimera-free marker reference used by reference UCHIME |
| `--skip-chimera-check` | | off | Explicitly omit UCHIME; affected isolates remain review-required |
| `--collapse` | | off | Collapse near-identical same-taxonomy sequences for the tree |
| `--collapse-threshold` | | `99.8` | Identity threshold (%) for collapsing |
| `--sequence-domain` | | `bacteria` | Domain/profile to process: `bacteria`, `archaea`, `fungi`, or `mixed`/`all`/`none` to disable filtering. Run archaea/fungi separately with matching `--ref`, `--alt-ref`, `--baseline-fasta`, `--previous-review`, and `--anchors` |
| `--phylum` | | — | Filter iTOL output to a specific phylum (does not affect novelty scoring) |
| `--target` | | — | FASTA to measure novelty against instead of the DB |
| `--baseline-fasta` | | — | Baseline/context FASTA to load before evaluating, e.g. Hungate |
| `--baseline-dataset` | | `Baseline` | Dataset label for `--baseline-fasta` and default cultured-baseline novelty pool |
| `--novelty-baseline-dataset` | | — | Existing DB dataset to include in the baseline novelty pool; repeatable |
| `--baseline-shorten-ids / --no-baseline-shorten-ids` | | `--no-baseline-shorten-ids` | Preserve baseline IDs exactly by default |
| `--force-rebuild` | | off | Force a full tree rebuild from scratch |
| `--anchors` | | bundled | Custom anchor FASTA for tree scaffolding |
| `--threads` | | `4` | CPU threads |
| `--neighbourhood-format` | | `png` | Local-clade figure format; PNG is the only supported output |
| `--user-colours` | | — | CSV mapping sequence IDs to custom hex colours for iTOL |
| `--group-phyla SPEC` | | — | Group phyla into a single colour in iTOL (repeatable) |
| `--functional` | | — | TSV file mapping sequence IDs to functional attributes (first column = ID; subsequent columns = attributes). Generates one iTOL file per column (DATASET_BINARY / DATASET_SIMPLEBAR / DATASET_COLORSTRIP). |

### Rolling project runs

Use one persistent SQLite database for the project and give each partner submission a unique `--dataset` label. The input FASTA may contain only the genuinely new isolates; `--partner-metadata` may be a cumulative table containing both new and previously stored partner IDs, so changed genome-selection flags are applied before the collection is reassessed.

For an auditable history, use a new output directory for every Performance Review and pass the preceding review with `--previous-review`. Reusing the same output directory is also supported and automatically resumes from `tree/current_alignment.fasta`, but replaces that directory's previous reports. Keep reference-database versions, domain, QC settings, collapse settings, anchors, and tree method constant when comparing reviews. Do not use `--force-rebuild` unless a deliberate full reanalysis is required.

Sequence IDs are persistent keys. Resubmitting an existing ID with the same sequence is allowed; resubmitting it with a different sequence is rejected rather than silently changing the project record. External-reference and cultured-baseline evidence should remain stable when inputs and settings are unchanged. Project crowding, nearest-project hit, genome-collection coverage, pangenome gap, selection-set role, and recommendation are deliberately recalculated because new isolates and updated WGS commitments change those questions.

Incremental runs retain the old MSA columns and add new sequences with MAFFT. The tree builder still estimates a new tree from the expanded alignment, so small changes in old topology or branch lengths remain scientifically possible. Use the preceding Performance Review directory as `--previous-review`, keep parameters/reference versions fixed, inspect `assessment/decision_changes.tsv`, and reserve `--force-rebuild` for deliberate recalibration. Strictly fixed old placements require a validated fixed-backbone phylogenetic-placement workflow and are not claimed by BranchManager.

Domain/profile handling:

BranchManager defaults to `--sequence-domain bacteria` for Filing Cabinet and Performance Review. This keeps bacterial 16S runs clean by filtering non-bacterial assignments before DB insertion and tree building where taxonomy is available. Process archaea and fungi as separate reviews with their own DB, outputs, baselines, references, and anchors:

```bash
branchmanager performance-review ... --sequence-domain archaea --ref archaeal_16s_refs.fasta --previous-review archaea_filing_cabinet -o archaea_review
branchmanager performance-review ... --sequence-domain fungi --ref fungal_its_or_18s_refs.fasta --anchors fungal_anchors.fasta -o fungi_review
```

Use `--sequence-domain mixed` when you intentionally want to keep all domains in one run.

---

### `branchmanager org-chart`

Extract all sequences matching a given taxon from the DB and build a focused phylogenetic tree for that group only.

**Fast path** (recommended): if `--from-dir` points to a directory containing `current_alignment.fasta`, sequences are sliced from the pre-built MSA and only FastTree is run — seconds to minutes for any sized group.

**Slow path**: if no existing alignment is found, a full MAFFT + FastTree build is performed.

```bash
# Domain-level (auto-detected keyword)
branchmanager org-chart --db project.db --taxon archaea      --from-dir filing_cabinet_out -o archaea_out
branchmanager org-chart --db project.db --taxon bacteria     --from-dir filing_cabinet_out -o bacteria_out

# Phylum — plain name (rank auto-detected) or GTDB-prefixed
branchmanager org-chart --db project.db --taxon Bacteroidota       --from-dir filing_cabinet_out -o bact_out
branchmanager org-chart --db project.db --taxon p__Bacillota       --from-dir filing_cabinet_out -o firm_out

# Family
branchmanager org-chart --db project.db --taxon f__Lachnospiraceae --from-dir filing_cabinet_out -o lachno_out

# Genus
branchmanager org-chart --db project.db --taxon g__Ruminococcus    --from-dir filing_cabinet_out -o rumino_out
```

Accepted taxon formats:

| Input | Detected rank | Notes |
|---|---|---|
| `archaea` / `bacteria` | domain | Auto-matched case-insensitively |
| `Bacillota` | phylum (fallback) | No prefix → defaults to phylum |
| `p__Bacillota` | phylum | GTDB prefix auto-detected |
| `f__Lachnospiraceae` | family | |
| `g__Ruminococcus` | genus | |
| `d__Archaea` | domain | Explicit GTDB domain prefix |

| Parameter | Required | Default | Description |
|---|---|---|---|
| `--db` | ✓ | — | Path to the BranchManager SQLite database |
| `-o / --out` | ✓ | — | Output directory |
| `--taxon` | ✓ | — | Taxon to extract (see table above) |
| `--rank` | | `auto` | Override rank detection: `domain` / `phylum` / `family` / `genus` / `species` |
| `--from-dir` | | — | Existing Filing Cabinet or Performance Review output directory with `current_alignment.fasta` or `tree/current_alignment.fasta` (fast path) |
| `--ref` | | — | Reference FASTA for orientation correction (slow path only) |
| `--anchors` | | bundled | Custom anchor FASTA |
| `--threads` | | `4` | CPU threads |
| `--min-seqs` | | `3` | Minimum matching sequences required to build a tree |
| `--no-tree` | | off | Skip tree building; only write taxonomy TSV and iTOL colour files |
| `--group-phyla SPEC` | | — | Group phyla into a single colour in iTOL (repeatable) |
| `--functional` | | — | TSV file mapping sequence IDs to functional attributes (first column = ID; subsequent columns = attributes). Generates one iTOL file per column (DATASET_BINARY / DATASET_SIMPLEBAR / DATASET_COLORSTRIP). |

Org Chart outputs:

| File | Contents |
|---|---|
| `org_chart_tree.nwk` | Focused newick tree (anchor-free, nodes labelled `NODE####`) |
| `org_chart_alignment.fasta` | Filtered alignment slice used for the tree |
| `org_chart_combined_taxonomy.tsv` | ID → taxonomy → confidence for matched sequences |
| `org_chart_sequence_list.tsv` | ID, taxonomy, confidence, dataset for matched sequences |
| `itol_phylum_colours.itol` | iTOL colour strip by phylum |
| `itol_dataset_membership.itol` | iTOL strip showing which dataset each sequence came from |

---

### `branchmanager label-maker`

Regenerate all iTOL colour files from taxonomy already stored in the DB — without re-classifying or rebuilding the tree. Useful when changing phylum groupings.

```bash
branchmanager label-maker \
  --db project.db \
  --out filing_cabinet_out \
  --group-phyla archaea \
  --group-phyla "Bacillota,Bacillota_I,Bacillota_A"
```

| Parameter | Required | Default | Description |
|---|---|---|---|
| `--db` | ✓ | — | Path to the BranchManager SQLite database |
| `-o / --out` | ✓ | — | Output directory (use your Filing Cabinet or Performance Review output directory) |
| `--include-datasets` | | all | Comma-separated dataset names to include |
| `--sequence-domain` | | — | Include only sequences from this domain/profile |
| `--group-phyla SPEC` | | — | Group phyla into a single colour (repeatable) |
| `--functional` | | — | TSV file mapping sequence IDs to functional attributes (first column = ID; subsequent columns = attributes). Generates one iTOL file per column (DATASET_BINARY / DATASET_SIMPLEBAR / DATASET_COLORSTRIP). |

---

## Phylum grouping for iTOL (`--group-phyla`)

Available on all four subcommands. Groups multiple phyla into a single colour entry in the iTOL legend.

| Spec format | Effect |
|---|---|
| `archaea` | All archaeal phyla (detected by `d__Archaea` domain) → one colour labelled **Archaea** |
| `bacteria` | All bacterial phyla → one colour labelled **Bacteria** |
| `"Bacillota,Bacillota_I,Bacillota_A"` | These three phyla → one colour; label = first name (`Bacillota`) |
| `"Firmicutes:Bacillota,Bacillota_I"` | These two phyla → one colour labelled **Firmicutes** |

Multiple `--group-phyla` arguments are independent and additive:

```bash
branchmanager performance-review ... \
  --group-phyla archaea \
  --group-phyla "Bacillota,Bacillota_I,Bacillota_A" \
  --group-phyla "Bacteroidota,Bacteroidota_A"
```

---

## Multiple dataset stacking

Datasets accumulate in the DB. Each Performance Review sees the growing baseline and project collection:

```bash
# Step 1 — baseline
branchmanager filing-cabinet --fasta hungate.fasta --db project.db --dataset Hungate ...

# Step 2 — first batch: novelty scored against Hungate
branchmanager performance-review --input batch1.fasta --db project.db --dataset Batch1 ...

# Step 3 — second batch: novelty scored against Hungate + Batch1
branchmanager performance-review --input batch2.fasta --db project.db --dataset Batch2 ...
```

---

## Clustering and tree redundancy reduction (`--collapse`)

When `--collapse` is enabled, BranchManager groups sequences sharing ≥ `--collapse-threshold` identity and the same taxonomy, keeping only one **cluster representative** per group for tree building.

#### Column Explanation

| Column | Meaning |
|---|---|
| ID | User-supplied sequence ID (preserved unless ID shortening is explicitly requested) |
| Taxonomy | Full GTDB lineage assigned by the classifier |
| BestHit | Closest reference genome accession (VSEARCH best hit) |
| ClassificationIdentity | % identity to the BestHit reference (VSEARCH alignment) |
| ClassificationConfidence | Confidence of the taxonomy assignment (from the taxa TSV); NA when taxonomy is parsed directly from reference FASTA headers (no confidence column available) |
| NearestHit | Closest sequence among YOUR previously submitted sequences (Filing Cabinet + prior Performance Reviews) |
| NearestIdentity | % identity to NearestHit |
| MatchesGE99 / GE97 / GE95 | Count of YOUR submitted sequences within 99% / 97% / 95% identity |
| NoveltyScore | 100 − NearestIdentity (higher = more novel vs. your collection) |
| Crowding | Neighbourhood density: crowded / moderate / sparse based on how many of your sequences are nearby |
| SequencingPriority | HIGH / MEDIUM / LOW — suggested priority for follow-up based on novelty + crowding |
| Reference* | Parallel nearest hit, identity, score, crowding, and priority against the selected external reference FASTA, usually GTDB |
| InTree | Yes = entered the phylogenetic tree; No = excluded (see ClusterRepresentative) |
| ClusterRepresentative | `self` = this sequence is in the tree; an ID = represented by that leaf after optional collapse; `duplicate` = retained isolate omitted from the final tree as an exact marker duplicate |
| ClusterSize | Total sequences in this cluster (1 = singleton) |
| ClusteredMembers | Semicolon-separated IDs of OTHER sequences collapsed under this representative |
| PlacementFlags | Warnings: LOW_CLASSIFICATION_IDENTITY, LOW_NEAREST_IDENTITY, NOVEL_BUT_ASSIGNED, etc. |

### `novelty_metrics.tsv`

This table keeps the novelty lenses separate:

- `Hungate*`: nearest hit, identity, query coverage, novelty, and density against priority Hungate baselines.
- `SecondaryBaseline*`: the same evidence against secondary cultured-rumen baseline datasets.
- `CulturedRumen*`: closest evidence across Hungate and secondary baseline tiers combined.
- `Project*`: nearest hit and density against all partner candidates accumulated across runs, including the current collection but excluding each query's self-hit.
- `Reference*`: nearest hit and density against the selected external reference FASTA supplied with `--ref`, usually GTDB.
- `GenomeCollection*` / `Pangenome*`: rolling genome coverage. Every baseline isolate counts because its genome is available; partner isolates count once selected. The default target is nine committed genomes per exact GTDB species or unresolved local clade.
- `SelectionGroupType`: `HUNGATE_BASELINE_EXTENSION` means a Hungate genome anchors the exact GTDB species, `SECONDARY_BASELINE_EXTENSION` means only secondary cultured-rumen genomes anchor it, and `CANDIDATE_ONLY_GROUP` has no baseline-genome anchor.
- `BaselineExtension*`: membership gate for a baseline-pangenome extension. The candidate and nearest cultured-rumen baseline must have the same exact GTDB species, with >=98.65% 16S identity across >=95% of the query by default. A failed gate is `PANGENOME_BOUNDARY_REVIEW`, not an automatic rejection: it may represent a separate candidate lineage.
- `BaselineRedundancy*`: hard eligibility gate for uncommitted candidates. By default, a nearest cultured-rumen hit at >=99.8% identity across >=95% of the query is reported as `BASELINE_REDUNDANT` and receives no panel rank.
- `SelectionDiversityDistance`: marginal patristic distance from baseline and already committed genome markers when tree distances are available. This drives farthest-first diversity capture; cultured-baseline, project, and GTDB-reference divergence provide the fallback, with MWL evidence used as an additional tie-break signal.

A species represented by one baseline genome therefore starts at `1/9`, leaving eight primary genome slots. Eligible partner isolates are ranked to extend that baseline pangenome; additional ranked isolates can be retained as extraction-failure/diversity backups. MWL evidence can order eligible isolates but cannot override the species/coverage boundary or the near-identical baseline exclusion.

---

## All output files

### Filing Cabinet outputs

| File | Description |
|---|---|
| `current_tree.nwk` | Backbone phylogenetic tree (newick, nodes labelled `NODE####`) |
| `current_alignment.fasta` | Full multiple sequence alignment (includes anchor sequences) |
| `filing_cabinet_combined_taxonomy.tsv` | ID → taxonomy → confidence for all registered baseline sequences |
| `taxonomy.tsv` | Raw classifier output; hit, identity, taxonomy, confidence, and alignment metric headers include the reference database name (for example, `BestHit_GTDB_r232`). |
| `itol_phylum_colours.itol` | iTOL colour strip by phylum |
| `itol_family_colours.itol` | iTOL colour strip by family |
| `itol_genus_colours.itol` | iTOL colour strip by genus |
| `filing_cabinet_dataset.itol` | iTOL strip marking all sequences with the dataset colour |
| `filing_cabinet_id_map.tsv` | Short ID → original FASTA header mapping |
| `filing_cabinet_collapsed_map.tsv` | Cluster representative → taxonomy → count |
| `filing_cabinet_collapsed_members.tsv` | Member → representative mapping |
| `tree_build_warnings.tsv` | Warnings about sequence quality or alignment issues |
| `tree_orientation_summary.tsv` | Per-sequence orientation audit (forward / RC / unknown) |
| `OUTPUT_EXPLANATIONS.tsv` | Manifest describing each output file |

### Performance Review and Hiring Panel outputs

| File | Description |
|---|---|
| `performance_review_dashboard.html` | Compact linked view of the current Hiring Panel recommendations |
| `assessment/sequence_assessment.tsv` | **Current-dataset full audit table.** Per-sequence novelty, taxonomy, crowding, priority, clustering, MWL matches, and placement flags for only the dataset supplied to this Assistant/Performance Review run |
| `assessment/selection_summary.tsv` | **Current-dataset SAB decision table.** Recommendation, evidence quality, component novelty evidence, same-species genome coverage, candidate-set role, rationale, and local-tree link |
| `assessment/sequencing_sets.tsv` | **Current-dataset nine-member diversity plan.** Prior datasets and baselines inform coverage and diversity, but only current input candidates receive rows in this run |
| `assessment/novelty_metrics.tsv` | Current-dataset novelty and crowding summary; searches use the cumulative baseline and project collection as comparison context |
| `assessment/neighbourhoods/clade_*.png` | Labelled local phylogenetic neighbourhoods (default). Nearby assessed isolates are grouped into one figure rather than duplicated across figures |
| `assessment/neighbourhoods/clade_*_pairwise_pident.tsv` | Complete long-form MSA percent-identity table for every pair of displayed tree leaves, including compared-column counts and the identity definition |
| `assessment/neighbourhoods/neighbourhood_manifest.tsv` | Maps every assessed sequence to its image and pident table, including the P1 identity anchor, forced nearest-baseline hits, assessed isolates, baseline leaves, and already-sequenced genomes shown |
| `assessment/mwl_matches.tsv` | Most Wanted List hits when `--mwl` is supplied |
| `assessment/cluster_summary.tsv` | Cluster-level prioritisation report when clustering/tree reports are available |
| `assessment/clusters.csv` | One consolidated membership table containing all clusters and all member isolates; replaces the former one-file-per-cluster directory |
| `assessment/backup_candidates.tsv` | Ranked alternative isolates within clusters when a primary candidate cannot be sequenced |
| `baseline/baseline_hits.tsv` | Nearest-hit report against Hungate or other provided baseline datasets |
| `baseline/nearest_baseline_hits_raw.tsv` | Raw nearest-hit output used for baseline novelty scoring |
| `taxonomy/<DB>.tsv` | Taxonomic assignment report for each reference database, e.g. `taxonomy/GTDB.tsv`, `taxonomy/GG2.tsv` |
| `taxonomy/all_databases.tsv` | Multi-database assignment summary; every database, including GTDB, has database-suffixed hit, identity, taxonomy, confidence, and alignment metric columns. |
| `taxonomy/tree_taxonomy.tsv` | ID → taxonomy → confidence table used for tree/iTOL metadata |
| `tree/current_tree.nwk` | Cumulative project-context tree incorporating baseline, prior partner datasets, and the current dataset |
| `tree/current_alignment.fasta` | Cumulative MSA used to build the project-context tree |
| `tree/tree_build_warnings.tsv` | Tree-quality warnings |
| `itol/phylum.itol` | iTOL colour strip by phylum |
| `itol/family.itol` | iTOL colour strip by family |
| `itol/genus.itol` | iTOL colour strip by genus |
| `itol/dataset_membership.itol` | iTOL strip showing which dataset each sequence belongs to |
| `itol/novelty.itol` | iTOL strip showing nearest-hit novelty |
| `ids/user_id_map.tsv` | Short ID → original FASTA header mapping for this run |
| `intermediate/` | QC, optional collapse, and classifier scratch outputs retained for debugging; exact marker duplicates retain separate isolate IDs in assessment reports |
| `logs/branchmanager.log` | Run log |

BranchManager keeps iTOL `DATASET_COLORSTRIP` files, one per metadata type. Older `TREE_COLORS` branch/range files and symbol-strip variants are removed because they encoded the same metadata in additional visual styles.

### `selection_summary.tsv`

This is the board-facing table for scientific advisory board discussions. It avoids an opaque adjusted score and keeps the evidence components visible: cultured-baseline novelty, rolling project coverage, GTDB-reference divergence, rank-aware MWL evidence, exact same-species genome counts, marker evidence quality, candidate-set role, and local-tree figure.

`SelectedForGenomeSequencing` and `GenomeAlreadySequenced` keep pending commitments distinct from completed genomes. `SelectedPendingGenomesSameSpecies` contributes to `CommittedGenomesSameSpecies` but not to nearest-available-genome identity.

The decision labels are:

| Decision | Meaning |
|---|---|
| `PRIORITISE - SET PRIMARY` | Fills one of the missing genomes needed to reach the per-species pangenome target |
| `RESERVE - SET BACKUP` | Additional phylogenetically spread isolate retained for extraction failure and strain-diversity capture |
| `STRONG CANDIDATE` / `SECONDARY CANDIDATE` | Useful evidence, but not currently assigned a primary/backup place in the working set |
| `SECONDARY - STRAIN DIVERSITY` | The numerical target is met, but the isolate adds a ranked diversity direction not represented by baseline/committed markers |
| `REVIEW - PANGENOME BOUNDARY` | The species, identity, or coverage evidence does not justify placing the isolate in the baseline-pangenome extension; review it as a possible separate lineage |
| `EXCLUDE - BASELINE REDUNDANT` | >=99.8% identity to a cultured baseline marker across >=95% query coverage; retained in reports but excluded from ranking |
| `LOWER PRIORITY - TARGET MET` | Available genomes plus committed pending selections meet the requested exact-GTDB-species target |
| `REVIEW BEFORE SELECTION` | Marker classification, QC, disagreement, or required comparison evidence needs review |
| `ALREADY SELECTED - GENOME PENDING` | The isolate is already committed to sequencing but no usable genome is available yet |
| `ALREADY SEQUENCED` | The current isolate already has an available genome and is not a new sequencing recommendation |

### Project-wide expansion rounds

The default Performance Review targets nine committed genomes per exact GTDB species, or per local phylogenetic clade when species assignment is unresolved. It writes up to nine ranked eligible partner candidates per group. Baseline and completed/pending partner genomes determine the remaining target gap, while candidates outside that gap remain ranked as backups or post-target diversity options.

Near-identical cultured-baseline matches are excluded only when both identity and alignment extent support the claim: >=99.8% marker identity and >=95% query coverage by default. This prevents a short high-identity fragment from being treated as redundant. Diversity is ranked first using marginal tree distance from baseline and committed genomes where available, then baseline/project/reference divergence and phylogenetic isolation. MWL rank and score add priority among otherwise comparable diversity choices but never override the baseline-redundancy gate.

After partner submissions have accumulated and sequencing statuses have been updated, run a Quarterly Review to choose the next tranche across the whole project:

```bash
branchmanager quarterly-review \
  --db project.sqlite \
  --genome-budget 24 \
  --backups-per-primary 1 \
  --tree runs/final/tree/current_tree.nwk \
  --alignment runs/final/tree/current_alignment.fasta \
  --partner-metadata partner_status.tsv \
  -o quarterly_review_01
```

`--genome-budget` is the number of new primary nominations, not including backups. Quarterly Review fills remaining nine-genome coverage gaps first, filters baseline-redundant candidates, then balances diversity across species. Within those tiers it uses marginal tree distance where available, nearest-genome marker identity, cultured-baseline novelty, GTDB-reference context, phylogenetic isolation, and MWL evidence. When `--alignment` is supplied, nearest available-genome identity is recomputed against the current genome collection, so newly completed genomes affect the next round immediately. It does not change `already_sequenced`; recommendations become available genomes only after metadata is updated in a later round.

Every current Performance Review stores its normalised assessment rows in the project database. For an older review made before assessment snapshots were introduced, import its full audit table once:

```bash
branchmanager quarterly-review \
  --db project.sqlite \
  --assessment runs/01_Batch1/assessment/sequence_assessment.tsv \
  --genome-budget 24 \
  --tree runs/01_Batch1/tree/current_tree.nwk \
  -o quarterly_review_01
```

Quarterly Review outputs are:

| File | Description |
|---|---|
| `next_genome_set.tsv` | Current `PRIMARY` nominations and their `BACKUP` isolates |
| `quarterly_review_summary.tsv` | Full audit containing all candidates, already-sequenced isolates, review exclusions, tiers, marginal distances, and reasons |
| `quarterly_review_manifest.tsv` | Round parameters, recommendation counts, and scientific scope |
| `neighbourhoods/clade_*.png` | Round-specific local-clade figures whose stars show this Quarterly Review's primary and backup nominations |

Quarterly Review is marker-gene decision support. Once genomes are available, strain-level expansion should be checked with ANI, phylogenomics, and pangenome structure rather than treating 16S distance as a strain boundary.

### `sequence_assessment.tsv` columns

| Column | Meaning |
|---|---|
| `ID` | Sequence ID supplied by the user; IDs are preserved unless shortening is explicitly requested |
| `GTDBTaxonomy` | Full authoritative GTDB lineage assigned by the classifier |
| `GTDBClassificationHit` | Closest GTDB reference accession used for classification |
| `GTDBClassificationIdentity` | % identity to the classification hit |
| `GTDBClassificationConfidence` | Confidence of the GTDB assignment |
| `BaselineNearestHit` | Closest sequence in the cultured baseline, e.g. Hungate |
| `BaselineNearestHitDataset` | Dataset label for that baseline hit |
| `BaselineNearestIdentity` | % identity to the cultured-baseline hit |
| `BaselineMatchesGE99 / GE97 / GE95` | Count of baseline sequences within 99% / 97% / 95% identity |
| `BaselineNoveltyScore` | 0–100 score; higher = more divergent and less crowded versus cultured isolates |
| `BaselineCrowding` | Baseline neighbourhood density: `crowded` / `moderate` / `sparse` / `isolated` |
| `BaselinePriority` | `HIGH` / `MEDIUM` / `LOW` heuristic from baseline novelty and density |
| `Project*` | Parallel novelty hit, identity, score, crowding, and priority against the rolling partner-candidate collection, excluding self |
| `GTDBReference*` | Parallel nearest hit, identity, score, and crowding against the authoritative GTDB reference |
| `PartnerID` | Partner acronym loaded from `--partner-metadata`, e.g. `QUB` or `UoG` |
| `SelectedForGenomeSequencing` | Whether the project has committed this isolate to sequencing; the genome may still be pending |
| `GenomeAlreadySequenced` | Whether the current partner isolate already has an available genome; a recommendation alone does not set this field |
| `GenomeCollection*` | Nearest hit and density across all baseline genomes plus already-sequenced partner genomes |
| `BaselineGenomesSameAssessmentSpecies` / `SequencedPartnerGenomesSameAssessmentSpecies` | Existing exact-species genome coverage by source; assessment species means GTDB for bacterial/archaeal Performance Reviews |
| `CommittedGenomesSameAssessmentSpecies` / `SelectedPendingGenomesSameAssessmentSpecies` / `PangenomeTarget` / `PangenomeGap` | Rolling completed-genome coverage, pending commitments, and the number of additional selections required |
| `SequencingSet*` | Stable group ID and proposed primary/backup/alternate action for the isolate |
| `LocalNeighbourhoodFigure` | Relative path to the local-clade PNG containing this assessed sequence |
| `TreeContextLeafCount` | Number of tree leaves displayed in that figure; this is display context, not a crowding score |
| `AssessedSequencesInTreeContext` | Number of current assessed isolates sharing the figure |
| `InTree` | `Yes` = in tree; `No` = excluded (see `ClusterRepresentative`) |
| `ClusterRepresentative` | `self`, a representative ID, or `duplicate` |
| `ClusterSize` | Sequences in this cluster |
| `ClusteredMembers` | Other sequences collapsed under this representative |
| `PlacementFlags` | Warning codes (e.g. `LOW_CLASSIFICATION_IDENTITY`) |

---

## Reference inputs — which flag to use

| Scenario | Recommended flags |
|---|---|
| Classify against GTDB reps FASTA | `--ref gtdb.fna --classify` |
| Input FASTA has embedded GTDB lineages in headers | `--taxa-assignments input.fasta` |
| Separate taxonomy TSV for the input sequences | `--taxa-assignments assignments.tsv` |
| GTDB FASTA + separate taxonomy TSV for the reference | `--ref gtdb.fna --taxa gtdb_taxonomy.tsv --classify` |

`--ref` = the reference database used for classification and tree orientation.  
`--taxa-assignments` = pre-computed taxonomy for the *input* sequences themselves.

### Partner sequencing metadata

`branchmanager performance-review` requires `--partner-metadata` / `--sequencing-metadata`. Provide this as a simple sidecar `.csv`, `.tsv`, or gzipped CSV/TSV alongside the FASTA. It must contain:

- A sequence ID column such as `sequence_id`, `isolate_id`, `sample_id`, or `ID`. Values may match IDs in the current FASTA or partner isolates already stored in the same project SQLite database.
- A partner acronym column such as `partner_id`, `partner`, or `partner_acronym`, with values like `QUB` or `UoG`.
- An optional selection-commitment column such as `selected_for_genome_sequencing`. Use `yes` after the project has actually committed that isolate to sequencing.
- An optional `addition_method` column: `marker_sequence`, `maldi_tof`, or `both`. Blank means unknown. This records how the isolate entered the catalogue; it does not itself grant a taxonomic assignment. MALDI-only isolates can receive a marker sequence later under the same durable `sequence_id`.
- An already-sequenced/genome-available column such as `already_sequenced`, `genome_available`, or `genome_sequenced`. This column is required and must be explicit `yes` or `no`.

Example:

```tsv
sequence_id	partner_id	selected_for_genome_sequencing	already_sequenced
Iso001	QUB	yes	no
Iso002	UoG	yes	yes
Iso003	UoG	no	no
```

For example, set `addition_method` to `maldi_tof` for an isolate first added from MALDI-TOF. When a marker is added later, update it to `marker_sequence` or `both` as appropriate. If marker QC fails and no previously accepted marker exists, BranchManager keeps the isolate without classifier taxonomy; a later passing marker submission can fill its marker slot.

Set `selected_for_genome_sequencing=yes` only after the advisory/project decision has been made. It reserves that isolate in rolling pangenome commitments but does not put it in the nearest-genome comparison pool. Set `already_sequenced=yes` only when a usable genome is genuinely available; that state does enter WGS coverage and nearest-genome calculations. `PRIMARY` and `BACKUP` are BranchManager recommendations and do not update either factual state automatically.

Keep one cumulative ledger for the project. It may contain rows from earlier submissions even when the current FASTA contains only new sequences. Each Performance Review refreshes matched stored isolates before recalculating commitments, available-genome coverage, and recommendations.

---

## ID shortening

By default, BranchManager preserves the FASTA IDs exactly as supplied. This is especially important for Hungate/baseline datasets and partner-provided isolate IDs.

- `--no-shorten-ids` (default) — keep supplied IDs.
- `--shorten-ids` — generate compact IDs when explicitly requested.
- `--no-baseline-shorten-ids` (default for Performance Review baselines) — keep baseline IDs exactly.
- `--baseline-shorten-ids` — generate compact baseline IDs when explicitly requested.

When IDs are shortened, BranchManager writes an ID map under `ids/`. When IDs are preserved, the map is still useful as an audit trail but should normally be identity-to-identity.

---

## Anchor sequences

Anchor sequences are **invisible scaffolding** constraining phylum-level topology during tree construction. They are pruned from the stored newick after each build and never appear in outputs, iTOL files, or novelty scoring.

By default, BranchManager uses the bundled anchor set at `src/branchmanager/data/rumen-gut_reference_anchors.fasta`: 26 NCBI RefSeq 16S sequences forming a rumen/gut topology scaffold. A companion metadata table (`src/branchmanager/data/rumen-gut_reference_anchors.tsv`) explains what each anchor represents, its rumen relevance, and whether it is a core rumen anchor or broader topology scaffold. The bundled set is used when `anchors` is `null` or omitted. Users can supply a custom FASTA with `--anchors /path/to/custom_anchors.fasta` (or configure that path in an Assistant profile) to replace the bundled set for a collection.

| Phylum | Type strain | Accession |
|---|---|---|
| Bacillota | *Ruminococcus albus* 7 | NR_025930 |
| Bacillota | *Clostridium butyricum* ATCC 19398 | NR_074545 |
| Bacillota | *Bacillus subtilis* 168 | NR_027552 |
| Bacillota | *Lactobacillus acidophilus* ATCC 4356 | NR_075051 |
| Bacillota | *Streptococcus equinus* NBRC 12553 | NR_113594 |
| Bacillota | *Butyrivibrio fibrisolvens* D1 | NR_044858 |
| Bacteroidota | *Bacteroides fragilis* ATCC 25285 | NR_041386 |
| Bacteroidota | *Prevotella bryantii* B14 | NR_044825 |
| Bacteroidota | *Porphyromonas gingivalis* ATCC 33277 | NR_040847 |
| Pseudomonadota | *Escherichia coli* K-12 | NR_102804 |
| Pseudomonadota | *Helicobacter pylori* 26695 | NR_073694 |
| Pseudomonadota | *Wolinella succinogenes* DSM 1740 | NR_043184 |
| Actinomycetota | *Bifidobacterium longum* NCC2705 | NR_040783 |
| Actinomycetota | *Streptomyces griseus* ATCC 23345 | NR_043823 |
| Spirochaetota | *Treponema pallidum* str. Nichols | NR_027243 |
| Spirochaetota | *Borrelia burgdorferi* B31 | NR_025890 |
| Fibrobacterota | *Fibrobacter succinogenes* S85 | NR_041558 |
| Fusobacteriota | *Fusobacterium nucleatum* ATCC 25586 | NR_026043 |
| Planctomycetota | *Planctomyces maris* DSM 8797 | NR_043399 |
| Verrucomicrobiota | *Akkermansia muciniphila* ATCC BAA-835 | NR_042817 |
| Thermotogota | *Thermotoga maritima* MSB8 | NR_043084 |
| Deinococcota | *Deinococcus radiodurans* R1 | NR_036779 |
| Cyanobacteriota | *Synechococcus elongatus* PCC 6301 | NR_043317 |
| Archaea (outgroup) | *Methanobrevibacter ruminantium* M1 | NR_044812 |
| Archaea (outgroup) | *Methanobacterium thermoautotrophicum* Marburg | NR_028243 |
| Archaea (outgroup) | *Sulfolobus acidocaldarius* DSM 639 | NR_043325 |

This is a pragmatic rumen/gut scaffold rather than a final rumen-only reference panel. The TSV flags which anchors are core rumen representatives and which are broader topology anchors that can be replaced as better vetted rumen-specific sequences become available.

To use a custom anchor file:
```bash
branchmanager assistant --config configs/my_collection_assistant.json \
  --anchors /path/to/my_anchors.fasta
```
Custom anchor headers must begin with `BRANCHMANAGER_REF_`:
```
>BRANCHMANAGER_REF_MyPhylum accession=NR_XXXXXX source=SILVA138
```

To refresh the bundled anchors from NCBI:
```bash
branchmanager build-anchors --email your@institution.ac.uk
```

---

## Interpreting outputs

### Selection evidence

| Metric | Meaning |
|---|---|
| `CulturedGap` | `LARGE` below 97%, `MODERATE` from 97% to <98.65%, and `SMALL` at >=98.65% against the cultured baseline |
| `ProjectCoverage` | Number of other rolling partner candidates at >=97% identity: uncovered, sparse (1-3), moderate (4-10), or dense (>10) |
| `ReferenceContext` | `HIGH DIVERGENCE` below 94.5%, `MODERATE DIVERGENCE` from 94.5% to <98.65%, or `CLOSE REFERENCE` at >=98.65% |
| `GenomeCoverage` | Exact same-GTDB-species available-genome count and remaining target gap; nearest-genome identity is shown separately as context |
| `EvidenceQuality` | Whether classification identity and placement warnings support using the marker result for selection |

The 94.5% and 98.65% values are commonly used full-length 16S heuristics for genus- and species-range divergence. BranchManager uses them as decision-support categories, not as declarations of a new genus or species. Sequence quality, alignment coverage, marker copy variation, phylogenetic placement, and ultimately genome-based analysis still matter.

Threshold basis: [Yarza et al. (2014), genus and higher-rank 16S boundaries](https://doi.org/10.1038/nrmicro3330) and [Kim et al. (2014), 98.65% 16S/ANI species-demarcation correspondence](https://pubmed.ncbi.nlm.nih.gov/24505072/). Their use in rumen-focused taxonomy is discussed by [Henderson et al. (2019)](https://doi.org/10.7717/peerj.6496), including limitations for unusually diverse genera.

The local-clade images expand around each assessed tree leaf to show a compact context (normally 8-30 leaves). Overlapping contexts are merged when their common clade remains compact. The nearest cultured-baseline hit for each grouped context is forced into the image even when reaching it crosses the normal leaf limit; less informative non-target leaves are pruned first. An orange diamond and `[SEQUENCED]` mark a confirmed available genome. A filled purple star marks the current sequence-next recommendation (`[P1]`, or another primary if more than one genome is missing); outlined purple stars mark recommended backups (`[B2]` and later ranks). `[ALT]` isolates are contextual alternatives, not current recommendations.

Each leaf label reports MSA pident relative to the figure's P1 sequence, and `clade_*_pairwise_pident.tsv` contains all displayed leaf-to-leaf comparisons. MSA pident is defined as identical A/C/G/T bases divided by alignment columns where both sequences contain an unambiguous A/C/G/T base; `ComparableACGTColumns` reports the overlap supporting the value. Terminal overhangs, gaps, and ambiguous bases are excluded rather than being misreported as substitutions. Internal nodes are inferred ancestors without observed marker sequences, so BranchManager does not assign pident to them; branch geometry and the scale bar remain modelled substitutions per site. These plots support visual inspection, while quantitative crowding remains the vsearch `MatchesGE*` count against explicitly named pools.

### Placement flags

| Flag | Meaning |
|---|---|
| `LOW_CLASSIFICATION_IDENTITY` | < 95% identity to best reference hit |
| `LOW_CONFIDENCE` | Classification confidence below threshold |
| `LOW_NEAREST_IDENTITY` | < 95% identity to nearest sequence in your collection |
| `NOVEL_BUT_ASSIGNED` | High novelty score despite a confident taxonomy assignment |

### Tree quality warnings

The tree is a context plot, not absolute proof of novelty. Always pair tree inspection with `sequence_assessment.tsv`. Warnings in `tree_build_warnings.tsv`:

| Warning | Meaning |
|---|---|
| `PARTIAL_16S_SEQUENCES` | Sequences shorter than 1200 bp reduce placement accuracy |
| `HIGH_N_CONTENT` | High ambiguous-base content can distort branch lengths |
| `VERY_SHORT_ALIGNED_FRAGMENTS` | Some aligned rows are mostly gaps after MAFFT |
| `MISSING_REFERENCE_ANCHORS` | Tree was built without topology scaffolding |

If a large `unknown` sector appears in iTOL, check:
- `combined_taxonomy.tsv` or `filing_cabinet_combined_taxonomy.tsv`
- `taxonomy.tsv` (raw classifier output)
- `taxonomy_input_warnings.tsv`

---

## Operational updates and project close-out

Recommendations are not genome facts. Use `status-meeting` for laboratory progress and `records-update` only when genome evidence exists.

```bash
# Laboratory/SAB lifecycle changes
branchmanager status-meeting --db project.sqlite --input isolate_status.tsv -o updates/status_01

# Completed assemblies, with automatic 90% completeness / 5% contamination defaults
branchmanager records-update --db project.sqlite --input genome_results.tsv -o updates/genomes_01

# Re-open the complete collection after the initial target-nine selection phase
branchmanager quarterly-review --db project.sqlite --genome-budget 24 \
  --from-dir runs/latest/03_performance_review_hiring_panel -o quarterly_reviews/round_02

# Final cumulative report after sequencing and genome QC
branchmanager annual-report --db project.sqlite -o final_annual_report
```

`isolate_status.tsv` requires `sequence_id` and `status`; supported states include `RECEIVED`, `TRACE_REVIEW`, `MARKER_QC_PASSED`, `PROPOSED`, `SAB_APPROVED`, `DNA_EXTRACTION_FAILED`, `SEQUENCED`, `GENOME_QC_FAILED`, and `GENOME_QC_PASSED`.

`genome_results.tsv` requires `sequence_id` and `genome_id`/`accession`. Recommended fields are `genome_status`, `completeness`, `contamination`, `gtdb_taxonomy`, `ani_cluster`, and `genome_path`. A `SEQUENCED` row with both quality metrics is automatically assigned genome-QC pass/fail using the configured thresholds. Explicit QC pass values that conflict with supplied metrics are rejected. Only QC-passed genomes update the available genome collection.

### Add genomes sequenced outside the selection workflow

Use `add-external-genomes` to register a project-level baseline genome collection. Supply one row per isolate/genome. Genome-level taxonomy is retained even when the isolate has no 16S marker. If marker sequences are available, pass a FASTA with headers matching `sequence_id`:

```bash
branchmanager add-external-genomes --config configs/rumen_gateway.json --db project.sqlite \
  --input external_genomes.tsv --marker-fasta external_16s.fasta \
  --dataset ExternalArchive --baseline-tier secondary -o imports/external_archive
```

Required columns are `sequence_id` and `genome_id` or `accession`. Optional fields include `gtdb_taxonomy`, `genome_qc_pass`, completeness/contamination, `ani_cluster`, and per-row `dataset`/`baseline_tier`. Genome availability is recorded as confirmed. Markerless isolates contribute genome/taxonomy coverage and are omitted from 16S similarity searches and trees.

### Add partner isolates identified by MALDI-TOF

Use `add-maldi-identifications` when a partner isolate has a MALDI-TOF identification but no marker sequence. Required columns are `sequence_id`, `partner_id`, and `taxon_name` or `taxon_id`. Include `taxonomic_rank` (`species`, `genus`, or `family`) with `taxon_name` to enable exact-rank soft context in later marker-based reviews. Optional columns are `identification_score`, `instrument`, `library_name`, `identification_date`, and `notes`.

```bash
branchmanager add-maldi-identifications --db project.sqlite \
  --input maldi_identifications.tsv --dataset UoG_MALDI -o imports/UoG_MALDI
```

The importer adds partner records with an empty marker slot and stores MALDI results, taxon IDs/names, and source provenance in a separate table. Matching MALDI taxon names can appear as a soft note in `SpeciesContext`; they do not become GTDB/classifier taxonomy and do not affect marker similarity, novelty, genome-coverage counts, or sequencing scores. The import report lists accepted and rejected rows.

For the Rumen Gateway `.xlsm` partner submission template, use `import-partner-workbook`. It reads the `Bulk Upload Template` sheet, uses `record_id` as the isolate ID, checks the MALDI method checkbox, and imports the MALDI result, instrument/database provenance, and a whitelist of non-identifying host and sample metadata. It can also register new markerless isolate records so a later passing marker submission can use the same ID. Existing genome availability and sequencing-selection decisions are preserved. Workbook taxon values in the field labelled `Taxon ID` are treated cautiously: numeric values are retained as taxon IDs, simple genus/species names may supply soft same-taxon context, and partner culture IDs or ambiguous strings are retained as unresolved observations for review. No workbook result is inserted into classifier taxonomy or selection scores.

```bash
branchmanager import-partner-workbook --db project.sqlite \
  --input JP_October_09.xlsm --dataset QUB_Workbook -o imports/QUB_workbook
```

The command writes `partner_workbook_import_report.tsv`. Its `metadata_details` column contains a JSON object of the workbook metadata saved for each isolate, including host/sample/culture fields and workbook provenance. Review rows marked `MALDI_RESULT_NEEDS_REVIEW`, `MALDI_SELECTED_NO_RESULT`, or `MALDI_FIELDS_WITHOUT_SELECTION` before relying on those records.

## Benchmarking recommendations

Before using selection output as a scientific conclusion, benchmark it against independently adjudicated outcomes on held-out cohorts. The repository includes a schema and evaluator in [`examples/benchmark_labels.template.tsv`](examples/benchmark_labels.template.tsv) and `branchmanager benchmark-selection`. Include all eligible candidates for each case, label the predeclared outcome independently of BranchManager, and include competing methods under the same `case_id` values. Evaluate multiple selection budgets (`precision@k`, `recall@k`) and overall rank quality (`mean average precision`); report cohort-level results and uncertainty, not just pooled values.

```bash
branchmanager benchmark-selection --input adjudicated_rankings.tsv \
  --out benchmark_summary.tsv --k 1 5 10
```

The output includes case-level bootstrap 95% percentile intervals (2,000 resamples by default; set `--bootstrap 0` to disable). This evaluator summarizes rankings; it does not supply validation data or establish biological validity. A journal submission still needs real, documented, held-out collection data; a comparator; sensitivity analyses for thresholds/reference versions; and error analysis for QC failures, unresolved taxonomy, and missing or biased references. The validation protocol and scientific scope checklist are in [collection profiles and scope](docs/collection-profiles-and-validation.md).

Every Performance Review, Quarterly Review, and Records Update works on a locked staged database and atomically publishes the result after validation. Annual Report writes an HTML dashboard, isolate and genome ledgers, workflow history, and the latest decision-change audit.

Office-stage reports are named `onboarding_report.tsv`, `onboarding_summary.json`, `status_meeting_report.tsv`, `records_update_report.tsv`, `annual_report.html`, `it_desk_report.tsv`, and `it_desk_summary.json`.

Before a production run:

```bash
branchmanager it-desk --db project.sqlite --ref gtdb_ssu_reps.fna.gz \
  --tree-method fasttree --strict -o it_desk
```

## What BranchManager does

BranchManager is aimed at ranking marker-gene sequence evidence for isolate follow-up. It helps identify lineages that may have been missed because of primer bias, sparse reference coverage, taxonomy lag, or conservative filtering. Many targets are not new in nature — they are often **new to the reference record** or **underrepresented in existing collections**.

That is why BranchManager combines reference-aware classification with novelty scoring, Sanger QC, baseline/reference comparisons, partner metadata, and warning layers, rather than treating every long branch as a discovery.
