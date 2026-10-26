# Collection profiles, supported scope, and validation

## Applying BranchManager to another isolate collection

BranchManager's collection-neutral data model uses a stable isolate/marker ID, a batch or dataset label, collection/partner metadata, a marker sequence (or supported raw-read input), and optional genome-availability evidence. Keep the durable isolate ID distinct from sequencing-run/read IDs. `partner_metadata` is a cumulative project ledger with one row per isolate ID; preserve source fields that are meaningful to the collection (host, environment, geography, collection date, culture conditions, accession, and consent/accession restrictions) in that ledger or a linked source table.

Start by copying `configs/project.template.json`. The profile has three parts:

* `project`: free-form collection tags and declared scope (organism/domain, marker, taxonomy authority, selection unit).
* `data`: named reference, baseline, partner-ledger, and anchor file locations for people maintaining the profile.
* `defaults`: command-to-CLI defaults for `assistant`, `filing-cabinet`, `performance-review`, and `quarterly-review`.
* `baselines`: one or more named baseline FASTAs, each assigned `priority` or `secondary`, with optional precomputed taxonomy or a request to classify.

The supported profile settings are validated against each command's CLI argument names. Set top-level `project_directory` when paths should be relative to a project working directory rather than the profile's location. Otherwise relative paths are resolved from the profile directory. Absolute paths and environment substitutions such as `${PROJECT_REF}` and `${PROJECT_REF:-/data/ref.fasta}` are also supported. Do not put credentials or private participant data in a profile. Use CLI options to override a profile for an individual run. Run manifests preserve the profile checksum and its `project` tags; keep the profile and reference files with the study's analysis archive.

Before using a new profile, decide and record:

1. What exactly counts as an isolate, a sequence, a collection baseline, and an eligible candidate; how replicate markers and repeated isolation events are linked.
2. Which marker and sequence technology are in scope, which taxonomic database/version is authoritative, and how taxonomy conflicts are handled.
3. Which baseline datasets represent physical isolates versus public references, and which are eligible for novelty and genome-coverage calculations.
4. The target unit and budget (species, lineage, functional group, or another defined unit), redundancy criteria, review states, and who can approve exceptions.
5. Which metadata are required, how missing values are represented, permitted data sharing, and identifiers/accessions needed for provenance.

The current candidate grouping and classification evidence are largely designed around bacterial/archaeal marker-gene taxonomies and GTDB-style rank strings. A profile changes input paths and defaults; it does not by itself make every workflow biologically valid for another marker or domain. Fungal/eukaryotic and non-16S collections require independent evaluation of classifiers, rank parsing, marker QC/chimera references, anchors, thresholds, and selection interpretation before use.

## Rumen Gateway profile

`configs/rumen_gateway.json` captures the current Rumen Gateway working setup rather than asserting that its choices are universal. It sets Hungate as the priority baseline and the 159-sequence Rumen Gateway `Secondary` collection as a secondary baseline with its supplied taxonomy assignments. It uses GTDB r232 as primary taxonomy, Greengenes2 2024.09, NCBI 16S, and SILVA 144 as cross-check references, the current partner catalogue and Most Wanted List, a default 24-primary quarterly budget, a nine-genome target/candidate panel, and the current marker redundancy and extension cutoffs. Relative paths resolve from the profile's `project_directory`. The profile omits a custom anchor path so BranchManager uses its bundled 26-sequence rumen/gut anchors by default; collections can provide their own versioned FASTA on the command line.

Use `branchmanager filing-cabinet --config configs/rumen_gateway.json --build-tree -o filing_cabinet` once to initialize/update both configured baseline datasets in the configured project database. Then use the profile with `branchmanager assistant --config configs/rumen_gateway.json ...` for partner submissions. The Assistant defaults are reusable reference, metadata, database, and assessment settings; the varying sample map, dataset name, and output directory remain per-run inputs. `branchmanager import-partner-workbook` separately imports partner workbook metadata and MALDI observations; it does not assign workbook-reported MALDI names as marker taxonomy or alter existing genome/selection decisions.

The bundled anchor file is a pragmatic 16S topology scaffold whose metadata explicitly flags rumen relevance. It contains both rumen/gut representatives and broad topology anchors, including representatives that are not appropriate ecological references for many collections. It is not a substitute for a collection-specific, versioned reference set. A publication should report the actual reference checksums/accessions, taxonomy release, anchor composition, and any replacements used.

## Scientific validation plan

The benchmark evaluator expects one row per candidate per method and case. `relevant=1` must reflect an outcome defined before reviewing BranchManager's ranking (for example, an independently approved candidate that yielded a QC-passing genome, if that is the intended endpoint). Do not label an isolate relevant merely because the toolkit nominated it. Include every candidate considered by each method, identify missing method/candidate combinations, and document exclusions.

For journal-level evidence, prepare a benchmark with:

* **Independent cohorts:** use held-out batches, sites, or time periods. Prevent isolate duplicates and near-identical marker sequences from leaking between development and evaluation partitions.
* **A comparator:** compare to a clearly defined current laboratory process and at least one simple baseline (such as random selection under the same eligibility and budget constraints). Apply identical input data and candidate budgets.
* **Predeclared endpoints:** report precision and recall at realistic budgets, average precision/ranking performance, and operational outcomes such as QC-passing genomes recovered per sequencing slot. State how unresolved or unavailable outcomes are handled.
* **Uncertainty:** report cohort-level estimates and confidence intervals (for example, bootstrap resampling by cohort or collection batch), not only a single pooled point estimate.
* **Sensitivity and ablation:** rerun across plausible identity/coverage cutoffs, target sizes, reference releases, and with/without each evidence layer. Show whether recommendations materially change.
* **Error analysis:** inspect false nominations and missed positives, poor-quality/chimeric sequences, taxa absent from references, taxonomy disagreements, metadata errors, and any collection/subgroup performance differences.
* **Reproducibility:** archive profile and input checksums, software commit/version, external-tool versions, reference releases, exact commands, exclusion rules, and the adjudication protocol. Keep the outcome labels separate from features used to generate recommendations.

The evaluator computes per-case `precision@k`, `recall@k`, and average precision, then averages each metric across cases for each method and cutoff. Here `k` is the number of candidates the lab can actually take forward under its sequencing budget. Precision@k is the proportion of those k choices that meet the predeclared outcome; recall@k is the proportion of all outcome-positive candidates in that case recovered within the k choices. It also reports case-level bootstrap 95% percentile intervals (2,000 replicates by default). Cases with no positive labels are omitted from recall and average-precision means; inspect the cohort data and report the number of cases contributing to those metrics. This script does not correct for repeated isolates or establish the validity of the labels; use a collection-aware resampling unit when cases are clustered by site, batch, or lineage.

`project.selection_budget` sets the default Quarterly Review primary-genome budget (`--genome-budget` overrides it); backup nominations are additional. Choose a benchmark cutoff that matches the actual primary sequencing budget, and report any backup outcome separately.

## External genomes and marker availability

Use `add-external-genomes` for genomes already available from another effort. This importer can add a genome-only isolate with an empty marker sequence, retain its genome accession and GTDB taxonomy, and mark the genome as available. If the isolate already exists as a candidate, it remains in its existing dataset and the genome evidence is linked to it. An optional marker FASTA attaches 16S records by isolate ID; only records with nonempty markers enter marker searches and trees. Genome-level taxonomy remains available for exact-species genome coverage even when no marker gene was recovered.

Use `add-maldi-identifications` for partner isolates that have a MALDI-TOF identification but no marker sequence. The taxon and its MALDI provenance are stored separately from classifier taxonomy. Exact matches at a supplied species, genus, or family rank may appear as soft context in later marker-based reviews; these records are not genome coverage and do not affect marker novelty or selection scores.

## Interpretation limits to state in a manuscript

BranchManager prioritizes follow-up based on the supplied sequences, metadata, baseline, taxonomy references, and configuration. Lack of a close reference hit is not evidence that an organism is new. Marker identity is not a strain boundary and cannot establish a novel species or pangenome structure. Selection thresholds are decision rules whose performance depends on the collection and marker. Human review, genome-quality checks, and genome-based comparisons remain necessary for biological claims. Publish both the intended decision use and the cases where recommendations should not be used.
