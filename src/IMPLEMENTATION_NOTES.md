# Expected-taxonomy validation update

This update separates partner-declared expected taxonomy from sequence-derived taxonomy in the Performance Review workflow.

## Partner metadata

The cumulative partner catalogue can now include an optional `taxon` column. The following aliases are accepted:

- `taxon`
- `expected_taxon`
- `organism` / `organism_name`
- `scientific_name`
- `taxon_name`

The value can be a species, genus, family, or a rank-prefixed lineage. The partner value is retained in sequencing metadata as `expected_taxon`. The isolate-ledger workflow records `expected_domain=Bacteria` by default.

## Validation behaviour

The pre-insertion classification step now applies these rules:

- No usable VSEARCH taxonomy: retain the isolate as `UNCLASSIFIED_EXPECTED_BACTERIAL`. This is an expected-domain assumption, not a fabricated sequence-derived assignment.
- Explicit non-bacterial domain: reject as `DOMAIN_DISAGREEMENT`.
- Partner taxon and VSEARCH lineage compatible: retain as `CONFIRMED_COMPATIBLE`.
- VSEARCH lineage compatible but less specific than the partner taxon: retain as `COMPATIBLE_LESS_SPECIFIC`.
- Explicit conflicting taxon: reject as `TAXON_DISAGREEMENT`.
- No partner taxon supplied, but a bacterial lineage is assigned: retain as `ASSIGNED_REQUESTED_DOMAIN`.
- Taxonomy present without an explicit domain or a usable taxon comparison: retain for review as `UNRESOLVED_DOMAIN_EXPECTED_BACTERIAL`.

Taxonomic validation remains separate from marker-sequence QC. A sequence can pass QC while its taxonomy remains unverified, or fail QC independently of its taxonomic assignment.

## Audit output

`assessment/expected_taxonomy_validation.tsv` records the expected domain, partner-declared taxon, observed lineage, validation status, retention decision, and explanation for every sequence processed by the early classification step.

## Comparison limitations

The current comparison is conservative and name/rank based. It handles common rank-prefixed GTDB/SILVA-style lineages, genus/species comparisons, and broader observed assignments. It does not resolve taxonomic synonyms or historical name changes. A synonym can therefore be reported as a disagreement and should be reviewed rather than treated as a definitive biological contradiction.

## Validation performed on this source update

- All Python modules compile successfully.
- Smoke tests cover exact agreement, species-level disagreement within a shared genus, genus-only less-specific assignments, and unclassified sequences.
- A synthetic classification-table test confirms that an unclassified sequence is retained, while an explicit archaeal assignment is rejected under the bacterial isolate profile.
