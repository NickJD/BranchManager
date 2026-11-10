# Outputs and interpretation

BranchManager produces machine-readable tables, run manifests, visualisations and decision-support reports. Exact filenames vary by command and workflow options.

## Evidence dimensions

Do not collapse the following measurements into one informal notion of “confidence”.

| Field/concept | What it means | What it does not mean |
|---|---|---|
| Identity | Percentage of identical bases within the aligned segment. | It does not say how much of the query or reference aligned. |
| Query coverage | Fraction of the query sequence covered by the alignment. | It is not reference coverage. |
| Reference/target coverage | Fraction of the reference sequence covered by the alignment. | It is not query coverage. |
| Alignment length | Number of aligned bases. | It does not by itself establish taxonomic identity. |
| Taxonomy support | Support for a lineage among the candidate reference hits under the configured consensus policy. | It is not automatically a calibrated probability. |
| Source classifier confidence | Confidence supplied by the classifier, when available. | It should not be reconstructed from percentage identity. |
| Assignment evidence level | A reporting category that combines configured evidence rules. | It is not a species-boundary definition or calibrated probability. |

Always read identity alongside query coverage, reference coverage, alignment length, assignment status and supporting-hit count.

## Taxonomic evidence

Classification outputs retain evidence about the reference hits and the consensus assignment. Depending on the workflow, useful files include:

- `taxonomy_<DATABASE>.tsv`
- `taxonomy_all_dbs.tsv`
- `classification_diagnostics_<DATABASE>.tsv`
- `classification_supporting_hits_<DATABASE>.tsv`
- `taxonomy_conflicts_<DATABASE>.tsv`
- `vsearch_search_audit.tsv`

A provisional or unresolved assignment is not equivalent to a confident species assignment. Diagnostic-only hits should remain labelled as diagnostic evidence and should not be treated as ordinary production-ladder matches.

## Novelty is project-relative

BranchManager's novelty and neighbourhood-density metrics compare new isolates against the configured project sequence pools, including the registered cultured baseline and previously assessed project sequences. External databases provide taxonomic/reference context and nearest-hit evidence, but a lack of a close external hit is not, by itself, a declaration of novelty.

Novelty can be affected by:

- Baseline composition and sampling depth.
- Reference database release and taxonomic coverage.
- Marker choice, primer bias and amplicon length.
- Sequence quality and query coverage.
- Thresholds and search policy.
- Whether the isolate has already been represented in a baseline or genome collection.

A novelty category should be interpreted with the underlying hit tables and evidence provenance.

## Nearest cultured and project sequence matches

The Performance Review dashboard reports nearest-hit evidence from the
comparison pools recorded in `assessment/novelty_metrics.tsv`. The comparison
names distinguish the sequence collections being searched; they should not be
treated as interchangeable baselines.

| Dashboard comparison | Novelty metrics prefix | Interpretation |
|---|---|---|
| Hungate cultured baseline | `Hungate` | Nearest matching sequence in the configured Hungate cultured baseline. |
| Secondary cultured baseline | `SecondaryBaseline` | Nearest matching sequence in the configured secondary cultured collection. |
| Combined cultured collection | `CulturedProject` | Nearest match across the combined cultured collection represented by the run's configuration. |
| Partner project marker collection | `Project` | Nearest matching marker sequence in the partner/project sequence collection. |
| GTDB reference genomes | `Reference` | Nearest reference sequence in the configured reference collection; shown in the novelty comparison evidence. |

For each comparison, the novelty metrics table can record:

- **Nearest hit** (`<PREFIX>NearestHit`): the identifier of the sequence reported as the closest match in that comparison pool.
- **Nearest identity** (`<PREFIX>NearestIdentity`): percentage identity over the aligned region.
- **Nearest query coverage** (`<PREFIX>NearestQueryCoverage`): percentage of the query sequence covered by that alignment.
- **Nearest alignment length** (`<PREFIX>NearestAlignmentLength`): number of aligned bases.

The dashboard's individual-isolate evidence panel displays the nearest-hit
identifier alongside identity, query coverage and aligned length. This is
intended to make the evidence actionable: a percentage alone does not tell the
reader which cultured isolate or project marker produced the match.

The dashboard's headline **Closest cultured identity** and **Closest project
identity** fields are concise summaries. The detailed comparison panel retains
the separate Hungate, secondary-baseline, combined-cultured and partner-project
results. In particular, a hit in the combined cultured collection should not
be assumed to originate from Hungate; inspect the hit identifier and the
individual source comparisons.

### Interpreting a nearest hit

A nearest-hit identifier identifies the sequence record returned by the
configured search and comparison workflow. It does not, by itself, establish
that two isolates are the same strain, species or biological entity. Interpret
it alongside:

1. Identity, query coverage and aligned length.
2. Marker QC status and any warning attached to the sequence.
3. The provenance and composition of the relevant baseline/project collection.
4. Taxonomic evidence from the configured reference databases.
5. The original FASTA and novelty metrics table when a result is important to a decision.

A project identity associated with `PASS_WITH_WARNINGS` or `FAIL_QC` is displayed
with a QC warning. The match is retained for inspection, but should not be
interpreted as equally reliable to a match from a `PASS_HIGH_CONFIDENCE` marker.
The warning concerns the quality of the query marker evidence; it does not
rewrite the search result or independently validate the identity of the hit.

If the nearest-hit identifier or identity is absent, check whether the relevant
comparison was run and whether its fields are present in `novelty_metrics.tsv`.
Do not infer that no biological match exists solely from a missing dashboard
value.

## Pangenome reports

`pangenome_hierarchy.tsv` describes GTDB-derived groups and their eligibility/representation. `sequencing_sets.tsv` contains candidate-oriented recommendations.

A group that lacks sufficient GTDB or tree/reference context should be treated as unresolved or requiring review. Do not read a missing target or an ineligible group as zero genome coverage.

For the source-of-taxonomy rules, see [Taxonomy and pangenome interpretation](taxonomy-and-pangenome.md).

## Selection recommendations

Candidate labels such as primary, backup, lower priority or manual review are recommendations under a configured policy and budget. They are not factual records of sequencing completion or biological novelty.

Keep factual lifecycle evidence separate from recommendations:

- Physical isolate availability.
- DNA extraction status.
- Library preparation and sequencing status.
- Genome assembly and QC status.
- GTDB taxonomy and ANI evidence.
- Manual review decisions.

Use the relevant status/records workflow to update facts. Do not infer that a nominated isolate has been sequenced merely because it appears in a selection output.

## Dashboards

Interactive HTML dashboards provide a navigable summary of the run, including collection overview, taxonomy, pangenome coverage, candidate evidence and relevant visual outputs where available. The Performance Review dashboard's individual-isolate panel also reports the nearest cultured and partner-project sequence hits, with identity, query coverage and aligned length where those fields are available.

Use `branchmanager dashboard /path/to/run_output` to regenerate the dashboard from an existing run. This rebuilds the HTML from the existing reports; it does not rerun the underlying biological analyses.

A dashboard is a presentation layer over the generated reports. When a dashboard and a source TSV appear to disagree, inspect the source TSVs and run manifest and report the discrepancy rather than relying on the visual summary alone.

## Reproducibility checklist

Retain the following with each analysis:

- BranchManager version (`branchmanager --version`).
- Run manifest and configuration profile.
- Original sample map and partner metadata.
- Original chromatograms or source FASTA.
- Reference/baseline FASTA and taxonomy files, versions and checksums.
- External tool versions.
- QC and selection settings.
- Manual review decisions.
- All TSV, JSON, FASTA, tree and HTML outputs.

Do not overwrite raw reads or generated review templates. Keep reviewed copies and audit decisions separately.
