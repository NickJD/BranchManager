# Taxonomy and pangenome interpretation

This guide describes the evidence-source contract used by the current BranchManager pangenome workflow.

## Two separate questions

BranchManager must keep domain/pathway routing separate from taxonomic identity.

### 1. Partner catalogue: which analytical pathway?

The partner catalogue's `taxon` field describes the expected biological domain or pathway for an isolate. For example:

| Partner `taxon` | Intended routing |
|---|---|
| `bacteria` | Bacterial workflow |
| `archaea` | Archaeal workflow |
| `eukaryota` | Eukaryotic workflow where supported |
| Missing/unknown | Apply configured defaults or request review |

A value such as `bacteria` is a domain-level routing instruction. It is not a genus or species assignment and must not be used to group all bacterial isolates into one pangenome taxon.

### 2. GTDB: where does the isolate belong taxonomically?

GTDB/reference genome taxonomy supplies the lineage used for pangenome grouping. The current pangenome logic reads the GTDB taxonomy fields and does not silently substitute partner taxonomy or a marker classifier when GTDB evidence is missing or unusable.

The source taxonomy and its provenance should remain visible in the assessment outputs.

## Grouping rank

The pangenome grouping rank follows the deepest usable rank supported by the GTDB assignment and the configured hierarchy. The system should not force a species assignment when the evidence only supports genus or a higher rank.

For example:

| Isolate | Partner taxon | GTDB assignment | Pangenome grouping |
|---|---|---|---|
| Isolate_001 | bacteria | `g__Winkia; s__Winkia sp...` | GTDB-supported species where usable; otherwise genus |
| Isolate_002 | bacteria | `g__Prevotella; s__...` | GTDB-supported species where usable; otherwise genus |
| Isolate_003 | bacteria | `g__Winkia; s__...` | GTDB-supported species where usable; otherwise genus |

The partner designation is identical in all three rows, but it does not determine their pangenome groups.

## Use lineage context, not only display names

Taxonomic names can be ambiguous or reused. Where the GTDB lineage is available, BranchManager uses the lineage through the grouping rank as part of the grouping key. This avoids merging unrelated taxa solely because they share a genus or species label.

The current database stores GTDB taxonomy strings rather than a separate normalised taxon-ID field for every rank. Consequently, lineage-derived keys are used to preserve context. Do not interpret these keys as official GTDB taxon IDs.

## GTDB assignment unavailable

If GTDB taxonomy is missing, unresolved or fails the configured assignment-status checks:

- The pangenome taxonomy is reported as unresolved.
- Marker classifier taxonomy may remain available as supporting/diagnostic evidence.
- Partner taxonomy remains available for pathway routing and expected-identity checks.
- The missing GTDB assignment must not be silently replaced with another source for pangenome grouping.
- The isolate must not receive a pangenome-gap-based sequencing priority until eligibility is established.

Unresolved is not equivalent to zero pangenome coverage and is not evidence that the isolate is novel.

## Pangenome eligibility is not pathway eligibility

An isolate may be correctly routed through the bacterial workflow but still be ineligible for pangenome assessment. Pangenome eligibility depends on usable GTDB taxonomy and the relevant reference-genome/tree context.

Keep these statuses separate in interpretation:

1. **Pathway status** – which domain-specific workflow the isolate follows.
2. **Taxonomic assignment status** – what GTDB supports for the isolate.
3. **Pangenome eligibility** – whether the relevant group can be assessed against reference-genome context.
4. **Selection status** – whether the isolate is recommended, held for review or not prioritised.

## Genome representation and gap calculations

Genome representation should use genome records with appropriate QC status and GTDB taxonomy. Marker-only classifier assignments and partner-catalogue text should not be counted as GTDB genome representation.

Genome-only records can contribute to taxon-level genome counts when their genome QC and GTDB evidence are suitable. However, a genome-only record without a marker/tree leaf does not create local marker-tree context by itself. This distinction matters when the pangenome target requires a local clade or tree neighbourhood.

A reported gap is meaningful only when the target group is eligible and its reference context is available. An unresolved group should be reported as requiring evidence/review, not as an empty taxon with a confirmed gap.

## What to inspect

When reviewing a pangenome recommendation, compare:

- `sequence_assessment.tsv` – isolate-level taxonomy, evidence and eligibility fields.
- `pangenome_hierarchy.tsv` – GTDB-supported groups, ranks, reference counts and eligibility.
- `sequencing_sets.tsv` – decision-facing candidate panels and priorities.
- The corresponding GTDB taxonomy and genome QC records.
- Tree/neighbourhood outputs where local tree context is required.

An isolate receiving a pangenome-gap-based priority should have a corresponding eligible GTDB group. If the hierarchy cannot establish that group, the recommendation should not rely on a pangenome gap.

## Scientific caution

GTDB is the configured authority for this project's pangenome grouping; this does not mean every GTDB assignment is necessarily correct for every isolate or marker. Review assignment status, identity, query coverage, reference coverage, supporting hits and taxonomy conflicts. Record the GTDB release and reference checksums used in the analysis.
