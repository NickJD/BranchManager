"""BranchManager CLI — argument parser for all subcommands.

Extracted from cli.py to reduce file size.
"""
import argparse
from branchmanager import __version__

from branchmanager.pipeline import (
    paper_trail as _paper_trail_module,
)
from branchmanager.pipeline import selection_sets as _selection_sets

# Re-exported from cli.py so the parser module is self-contained.
SEQUENCE_DOMAIN_CHOICES = (
    'bacteria', 'bacterial',
    'archaea', 'archaeal',
    'fungi', 'fungal',
    'eukaryota', 'eukarya',
    'mixed', 'all', 'none',
)


def build_parser():
    _Fmt = type('_Fmt', (
        argparse.ArgumentDefaultsHelpFormatter,
        argparse.RawDescriptionHelpFormatter,
    ), {'max_help_position': 40, 'width': 120})
    paper_trail_policy = _paper_trail_module.DEFAULT_QC_POLICY
    selection_target = _selection_sets.DEFAULT_PANGENOME_TARGET
    selection_panel_size = _selection_sets.DEFAULT_CANDIDATE_SET_SIZE
    baseline_redundancy_identity = _selection_sets.DEFAULT_BASELINE_REDUNDANCY_IDENTITY
    baseline_redundancy_coverage = _selection_sets.DEFAULT_BASELINE_REDUNDANCY_QUERY_COVERAGE
    baseline_extension_identity = _selection_sets.DEFAULT_BASELINE_EXTENSION_MIN_IDENTITY
    baseline_extension_coverage = _selection_sets.DEFAULT_BASELINE_EXTENSION_MIN_QUERY_COVERAGE
    parser = argparse.ArgumentParser(
        prog='branchmanager',
        description=(
            'BranchManager — marker-gene QC, taxonomy, novelty scoring, and isolate prioritisation toolkit.\n\n'
            'Subcommands:\n'
            '  mailroom           Inventory AB1 deliveries and build a batch map.\n'
            '  interview          Run standalone AB1 QC from a Mailroom batch.\n'
            '  background-check   Pre-classify reference collections once and reuse the evidence.\n'
            '  genomic-16s-selection Extract genome-identifiable 16S representatives with Barrnap + VSEARCH + optional GTDB-Tk comparison.\n'
            '  filing-cabinet     Register the cultured baseline and backbone tree.\n'
            '  onboarding         Validate a partner submission before project state changes.\n'
            '  paper-trail        Interpret chromatograms and assemble primer reads.\n'
            '  sanger-studio      Inspect and assemble Sanger reads in a local browser GUI.\n'
            '  performance-review Assess taxonomy, novelty, phylogeny, and sequencing candidates.\n'            '  dashboard          Regenerate the interactive dashboard from an existing run output directory.\n'
            '  status-meeting     Import factual isolate lifecycle changes.\n'
            '  records-update     Import completed-genome, QC, GTDB, and ANI evidence.\n'
            '  add-maldi-identifications Register partner isolates with MALDI-TOF taxonomy and no marker.\n'
            '  import-genome-16s  Register genomic 16S (GENOME_BARRNAP) evidence and existing genomes.\n'
            '  import-genome-evidence Import partner genome FASTAs with GTDB-Tk and CheckM2 evidence.\n'
            '  marker-review-apply Persist manual marker-review decisions in the database.\n'
            '  import-partner-workbook Import partner metadata and MALDI status from an .xlsm workbook.\n'
            '  quarterly-review   Nominate a later genome tranche from the complete project.\n'
            '  exit-interview     Remove isolates from active project state with a retained audit.\n'
            '  annual-report      Produce the cumulative end-of-project report.\n'
            '  it-desk            Check dependencies and project inputs before production.\n'
            '  assistant          Run the complete raw-trace-to-Hiring-Panel workflow.\n'
            '  org-chart          Extract a focused tree and iTOL files for a specific taxon.\n'
            '  label-maker        Regenerate iTOL annotation files from stored taxonomy.\n\n'
            'Typical workflow:\n'
            '  P. branchmanager mailroom --read-dir All_AB1 --manifest manifest.csv --dataset QUB_01 --forward-primer 63F -o QUB_01\n'
            '  Q. branchmanager interview --mailroom QUB_01 -o QUB_01_interview\n'
            '  0. branchmanager background-check --dataset hungate16s=hungate.fasta --ref gtdb.fna --taxa gtdb_tax.tsv -o background_check_out\n'
            '  1. branchmanager filing-cabinet --fasta baseline.fasta --db project.db --dataset Hungate --taxa-assignments background_check_out/pipeline_taxonomy.tsv --build-tree -o filing_cabinet_out\n'
            '  2. branchmanager performance-review --input new_seqs.fasta --partner-metadata new_seqs_metadata.tsv --db project.db --dataset Batch1 --ref gtdb.fna --mwl MWL.xlsx -o review_out\n'
            '  3. branchmanager quarterly-review --db project.db --genome-budget 24 --tree review_out/tree/current_tree.nwk -o quarterly_review_01\n'
            '  4. branchmanager org-chart --db project.db --taxon archaea --from-dir filing_cabinet_out -o archaea_out\n'
            '  5. branchmanager genomic-16s-selection --genome-dir assemblies --reference gtdb_r232_bac_arch_ssu_reps.fna.gz --gtdb-summary gtdbtk.bac120.summary.tsv -o genomic_16s\n'

        ),
        formatter_class=_Fmt,
    )
    parser.add_argument('-v', '--version',
        action='version', version=f'%(prog)s {__version__}', help='Show the installed BranchManager version and exit.')
    parser.add_argument('--config', help='Project profile JSON. Profile values provide defaults; command-line options override them.')
    sub = parser.add_subparsers(dest='command')

    def add_project_config_option(command_parser):
        command_parser.add_argument('--config', help='Project profile JSON. Profile values provide defaults; command-line options override them.')

    filing_cabinet_parser = sub.add_parser(
        'filing-cabinet',
        help='Filing Cabinet: register a cultured baseline and build the reference tree.',
        description=(
            'Filing Cabinet registers a baseline FASTA dataset (e.g. Hungate 16S) in the DB, optionally classifies\n'
            'sequences against a reference (GTDB/SILVA), collapse near-identical sequences,\n'
            'and build the backbone phylogenetic tree.\n\n'
            'This is always the first step. All subsequent Performance Reviews measure novelty against\n'
            'sequences stored by this command.'
        ),
        formatter_class=_Fmt,
    )
    add_project_config_option(filing_cabinet_parser)
    filing_cabinet_parser.add_argument('--fasta', required=False,
        help='Input FASTA file containing the baseline sequences to load.')
    filing_cabinet_parser.add_argument('--db', required=True,
        help='Path to the BranchManager SQLite database (created if it does not exist).')
    filing_cabinet_parser.add_argument('-o', '--out', required=False, default='.',
        help='Output directory for tree, iTOL files, and reports (default: current directory).')
    filing_cabinet_parser.add_argument('--dataset', required=False,
        help='Label for this dataset stored in the DB (e.g. Hungate). Used to colour iTOL strips.')
    filing_cabinet_parser.add_argument('--baseline-tier', choices=['priority', 'secondary'], default=None,
        help='Cultured baseline tier for this dataset: priority for Hungate, secondary for other project baseline genomes. Inferred from --dataset when omitted.')
    filing_cabinet_parser.add_argument('--shorten-ids', dest='shorten_ids',
        action=argparse.BooleanOptionalAction, default=False,
        help='Replace input headers with compact IDs (e.g. HUN001). Default is to preserve the IDs exactly as supplied.')
    filing_cabinet_parser.add_argument('--classify', action='store_true',
        help='Classify sequences against --ref and store taxonomy in the DB. Requires --ref.')
    filing_cabinet_parser.add_argument('--build-tree', action='store_true',
        help='Build the backbone MAFFT + FastTree phylogenetic tree after loading.')
    filing_cabinet_parser.add_argument('--ref', required=False,
        help='Reference FASTA (GTDB/SILVA reps) for classification and tree orientation. Preferred over --taxa-assignments for externally classified inputs.')
    filing_cabinet_parser.add_argument('--taxa', required=False,
        help='Taxonomy table matching IDs in --ref (TSV/CSV, optionally .gz: id<TAB>lineage or id,lineage). Optional when --ref FASTA headers already contain lineages.')
    filing_cabinet_parser.add_argument('--ref-name', dest='ref_name', required=False, default=None,
        help='Display name for the primary reference database (default: derived from --ref filename). Used to label taxonomy columns.')
    filing_cabinet_parser.add_argument('--alt-ref', dest='alt_ref', action='append', default=None, metavar='FASTA',
        help='Additional reference FASTA to classify against (repeatable). Produces extra taxonomy columns in output files.')
    filing_cabinet_parser.add_argument('--alt-taxa', dest='alt_taxa', action='append', default=None, metavar='TABLE',
        help='Taxonomy TSV/CSV for the corresponding --alt-ref (positionally paired; repeatable; .gz accepted).')
    filing_cabinet_parser.add_argument('--alt-ref-name', dest='alt_ref_name', action='append', default=None, metavar='NAME',
        help='Display name for the corresponding --alt-ref (positionally paired; repeatable). Default: derived from filename.')
    filing_cabinet_parser.add_argument('--main-ref', dest='main_ref', required=False, default=None,
        help='Name of the reference database to use as the primary taxonomy source (default: primary --ref). Must match one of the --ref-name / --alt-ref-name values.')
    filing_cabinet_parser.add_argument('--taxa-assignments',
        dest='taxa_assignments', required=False,
        help='Pre-computed taxonomy assignments for the INPUT sequences (TSV/CSV, optionally .gz: query_id + lineage, or a FASTA with embedded lineages). Use this instead of --classify when you already have taxonomy.')
    filing_cabinet_parser.add_argument('--collapse', action='store_true',
        help='Collapse sequences that share ≥ --collapse-threshold identity AND the same taxonomy into a single representative for the tree. Saves time and reduces visual clutter.')
    filing_cabinet_parser.add_argument('--collapse-threshold', type=float, default=99.8,
        help='Identity threshold (percent) for collapsing duplicate-like sequences (default: 99.8).')
    filing_cabinet_parser.add_argument('--sequence-domain', '--organism-domain', dest='sequence_domain',
        choices=SEQUENCE_DOMAIN_CHOICES, default=None,
        help=(
            'Sequence/domain profile for this Filing Cabinet. Default behaviour is bacteria. '
            'Use archaea for archaeal 16S, fungi for fungal/eukaryotic runs with suitable refs/anchors, '
            'or mixed/all/none to disable domain filtering.'
        ))
    filing_cabinet_parser.add_argument('--anchors', required=False, default=None,
        help='Custom reference anchor FASTA for tree topology scaffolding. By default uses the bundled 26-sequence rumen/gut anchor set (src/branchmanager/data/rumen-gut_reference_anchors.fasta).')
    filing_cabinet_parser.add_argument('--threads', type=int, required=False, default=4,
        help='Number of CPU threads for MAFFT and VSEARCH (default: 4).')
    filing_cabinet_parser.add_argument('--tree-method', dest='tree_method',
        choices=['fasttree', 'iqtree', 'iqtree-fast'], default='fasttree',
        help=(
            'Phylogenetic tree-building backend (default: fasttree). '
            'fasttree: approximate ML, GTR+CAT — fast. '
            'iqtree: full ML, GTR+G+I — more accurate and stable topology (slower; '
            'recommended for publication-quality trees or when FastTree produces unstable clades). '
            'iqtree-fast: IQ-TREE 2 with -fast flag — good compromise for exploratory runs. '
            'Requires iqtree2 in PATH when using iqtree/iqtree-fast.'
        ),
    )
    filing_cabinet_parser.add_argument('--colours', required=False,
        help='CSV file mapping sequence IDs to custom hex colours for iTOL (columns: id, colour).')
    filing_cabinet_parser.add_argument(
        '--group-phyla', dest='group_phyla', action='append', default=None, metavar='SPEC',
        help=(
            'Collapse multiple phyla into a single colour in iTOL legends. Repeatable. '
            'Formats: "archaea" (all archaeal phyla), "bacteria" (all bacterial phyla), '
            '"Bacillota,Bacillota_I" (explicit list; label = first name), '
            '"Firmicutes:Bacillota,Bacillota_I" (named group).'
        ),
    )
    filing_cabinet_parser.add_argument('--functional', dest='functional', required=False, default=None,
        help=(
            'TSV file mapping sequence IDs to functional attributes (pathways, functions, '
            'traits, scores, etc.). Header row required; first column = sequence ID; '
            'subsequent columns = one functional attribute each. '
            'One iTOL file is generated per column: binary (0/1/yes/no) → DATASET_BINARY, '
            'numeric → DATASET_SIMPLEBAR, categorical → DATASET_COLORSTRIP.'
        ),
    )
    filing_cabinet_parser.add_argument('--draft-rumen-functions', dest='draft_rumen_functions',
        action='store_true', default=False,
        help=(
            'Auto-generate a draft rumen functional-group annotation from the output taxonomy. '
            'Maps each sequence to a broad ruminant microbiome functional category '
            '(e.g. Cellulolytic/Fibrolytic, Methanogenic Archaea, Butyrate Producers) '
            'and writes rumen_functions_draft.tsv + itol_func_Rumen_Functional_Group.itol. '
            'The draft TSV can be edited and re-supplied via --functional in future runs.'
        ),
    )

    performance_review_parser = sub.add_parser(
        'performance-review',
        help='Performance Review: assess taxonomy, novelty, phylogeny, and sequencing candidates.',
        description=(
            'Assess new partner 16S isolate sequences against the project baseline.\n\n'
            'The workflow classifies against GTDB (primary), optionally cross-checks NCBI/GG2/SILVA\n'
            'as --alt-ref databases, scores novelty against separate cultured-baseline and rolling partner\n'
            'collections, updates the tree, and optionally matches GTDB taxonomy against\n'
            'the Most Wanted List via --mwl.\n\n'
            'Provide --baseline-fasta for context datasets such as Hungate when they have not already\n'
            'been registered with `branchmanager filing-cabinet`. Novelty is always relative to YOUR submitted data,\n'
            'not the full external reference. Each successive Performance Review extends the partner collection and\n'
            'updates same-species genome coverage and candidate sets.'
        ),
        formatter_class=_Fmt,
    )
    add_project_config_option(performance_review_parser)
    performance_review_parser.add_argument('--input', required=True,
        help='FASTA file of new sequences to analyse.')
    performance_review_parser.add_argument('--marker-source-type',
        choices=['auto', 'SANGER_AB1', 'PARTNER_FASTA', 'PACBIO_HIFI', 'MIXED'],
        default='auto',
        help='Provenance of the marker sequence. Auto treats direct --input FASTA as PARTNER_FASTA; Assistant resolves Sanger/PacBio automatically.')
    performance_review_parser.add_argument('--db', required=True,
        help='Path to the BranchManager SQLite database (must have been initialised with `branchmanager filing-cabinet`).')
    performance_review_parser.add_argument('-o', '--out', required=True,
        help='Output directory for this run (sequence_assessment.tsv, novelty_metrics.tsv, tree, iTOL files, etc.).')
    performance_review_parser.add_argument('--dataset', required=True,
        help='Label for this batch of sequences stored in the DB (e.g. Batch1). Used in iTOL dataset-membership strip.')
    performance_review_parser.add_argument('--ref', required=False,
        help='Reference FASTA (GTDB/SILVA reps) used for classification and tree orientation. Use the same file as the Filing Cabinet.')
    performance_review_parser.add_argument('--taxa', required=False,
        help='Taxonomy table matching IDs in --ref (TSV/CSV, optionally .gz: id<TAB>lineage or id,lineage).')
    performance_review_parser.add_argument('--ref-name', dest='ref_name', required=False, default=None,
        help='Display name for the primary reference database (default: derived from --ref filename).')
    performance_review_parser.add_argument('--alt-ref', dest='alt_ref', action='append', default=None, metavar='FASTA',
        help='Additional reference FASTA to classify against (repeatable). Adds extra taxonomy columns to sequence_assessment.tsv and taxonomy_all_dbs.tsv.')
    performance_review_parser.add_argument('--alt-taxa', dest='alt_taxa', action='append', default=None, metavar='TABLE',
        help='Taxonomy TSV/CSV for the corresponding --alt-ref (positionally paired; repeatable; .gz accepted).')
    performance_review_parser.add_argument('--alt-ref-name', dest='alt_ref_name', action='append', default=None, metavar='NAME',
        help='Display name for the corresponding --alt-ref (positionally paired; repeatable).')
    performance_review_parser.add_argument('--main-ref', dest='main_ref', required=False, default=None,
        help='Name of the authoritative assessment database. Bacterial/archaeal Performance Review requires GTDB; other databases are reporting cross-checks.')
    performance_review_parser.add_argument('--mwl', dest='mwl', required=False, default=None,
        help='Most Wanted List workbook/TSV/CSV. GTDB taxonomy is matched against this list and MWL columns are added to sequence_assessment.tsv.')
    performance_review_parser.add_argument('--mwl-sheet', dest='mwl_sheet', required=False, default='MWL_V1',
        help='Sheet name to read from an MWL .xlsx workbook (default: MWL_V1).')
    performance_review_parser.add_argument('--mwl-min-rank', dest='mwl_min_rank', required=False, default='p',
        choices=['domain', 'd', 'phylum', 'p', 'class', 'c', 'order', 'o', 'family', 'f', 'genus', 'g', 'species', 's'],
        help='Minimum matched rank required for an MWL hit (default: phylum). Domain-only MWL entries still match at domain.')
    performance_review_parser.add_argument('--partner-metadata', '--sequencing-metadata', dest='partner_metadata', required=False, default=None,
        help='Cumulative CSV/TSV ledger with sequence IDs, partner acronyms, optional taxon/expected_taxon and taxon_rank, optional selected_for_genome_sequencing, and required already_sequenced status. For partner isolates, the catalogue taxon is authoritative; classifier disagreements are diagnostic only. selection is a commitment; already_sequenced means a genome is available. .gz is accepted.')
    performance_review_parser.add_argument('--marker-qc', dest='marker_qc', default=None,
        help='Paper Trail/Merge Meeting assembly_report.tsv sidecar. Auto-discovered beside --input when omitted.')
    performance_review_parser.add_argument('--marker-review', dest='marker_review', default=None,
        help='CSV/TSV manual-review decisions for PASS_WITH_WARNINGS or unverified marker sequences.')
    performance_review_parser.add_argument('--accept-unverified-marker-qc', action='store_true', default=False,
        help='Permit processing FASTA without QC provenance; this does not approve it for selection and requires the explicit authorizer/reason below.')
    performance_review_parser.add_argument('--unverified-override-authorizer', default='',
        help='Person authorizing --accept-unverified-marker-qc (required with the flag).')
    performance_review_parser.add_argument('--unverified-override-reason', default='',
        help='Auditable reason for --accept-unverified-marker-qc (required with the flag).')
    performance_review_parser.add_argument('--baseline-fasta', dest='baseline_fasta', required=False, default=None,
        help='Optional baseline/context FASTA to load before evaluating the new sequences (e.g. Hungate 16S).')
    performance_review_parser.add_argument('--baseline-dataset', dest='baseline_dataset', required=False, default='Baseline',
        help='Dataset label for --baseline-fasta in the DB and default cultured-baseline novelty pool (default: Baseline). Must differ from --dataset.')
    performance_review_parser.add_argument('--baseline-tier', choices=['priority', 'secondary'], default=None,
        help='Cultured baseline tier for --baseline-fasta: priority for Hungate, secondary for other project baseline genomes. Inferred from --baseline-dataset when omitted.')
    performance_review_parser.add_argument('--novelty-baseline-dataset', dest='novelty_baseline_datasets', action='append', default=[],
        help='Existing DB dataset label to include in the baseline/cultured novelty pool (repeatable; useful for Hungate plus other cultured isolate sets).')
    performance_review_parser.add_argument('--baseline-taxa-assignments', dest='baseline_taxa_assignments', required=False, default=None,
        help='Pre-computed taxonomy for --baseline-fasta (TSV/CSV, optionally .gz: sequence_id + lineage + optional confidence, or embedded-lineage FASTA). Skips baseline classification.')
    performance_review_parser.add_argument('--baseline-skip-classify', dest='baseline_skip_classify', action='store_true', default=False,
        help='Load --baseline-fasta into the DB without classifying it. Novelty still uses the baseline sequences, but taxonomy/iTOL context may be sparse.')
    performance_review_parser.add_argument('--baseline-colours', dest='baseline_colours', required=False, default=None,
        help='Optional CSV with baseline sequence colours, in the same format as filing-cabinet --colours.')
    performance_review_parser.add_argument('--baseline-shorten-ids', dest='baseline_shorten_ids',
        action=argparse.BooleanOptionalAction, default=False,
        help='Replace baseline FASTA headers with compact IDs. Default is to preserve the IDs exactly as supplied.')
    performance_review_parser.add_argument('--taxa-assignments',
        dest='taxa_assignments', required=False,
        help='Pre-computed taxonomy for the INPUT sequences (TSV/CSV, optionally .gz: query_id + lineage, or embedded-lineage FASTA).')
    performance_review_parser.add_argument('--previous-review', dest='previous_review', required=False,
        help='Path to the previous Filing Cabinet or Performance Review output directory. Used to seed the tree alignment so only new sequences need aligning.')
    performance_review_parser.add_argument('--shorten-ids', dest='shorten_ids',
        action=argparse.BooleanOptionalAction, default=False,
        help='Replace input headers with compact IDs. Default is to preserve the IDs exactly as supplied.')
    performance_review_parser.add_argument('--min-len', dest='min_len', type=int, default=800,
        help='Minimum sequence length to retain (bp, default: 800). Shorter sequences are filtered out.')
    performance_review_parser.add_argument('--max-n-percent', dest='max_n_percent', type=float, default=5.0,
        help='Maximum percentage of ambiguous (N) bases allowed (default: 5.0).')
    performance_review_parser.add_argument('--chimera-ref', dest='chimera_ref', default=None,
        help='Curated chimera-free marker reference for UCHIME. Defaults to the primary classification reference.')
    performance_review_parser.add_argument('--skip-chimera-check', dest='skip_chimera_check', action='store_true', default=False,
        help='Explicitly skip reference-based chimera screening. The omission is recorded and marker evidence is downgraded.')
    performance_review_parser.add_argument('--collapse', action='store_true',
        help='Collapse near-identical same-taxonomy sequences into representatives for the tree.')
    performance_review_parser.add_argument('--collapse-threshold', type=float, default=99.8,
        help='Identity threshold (percent) for collapsing (default: 99.8).')
    performance_review_parser.add_argument('--sequence-domain', '--organism-domain', dest='sequence_domain',
        choices=SEQUENCE_DOMAIN_CHOICES, default=None,
        help=(
            'Sequence/domain profile for this Performance Review. Omitted means bacteria. '
            'Use archaea for archaeal runs, fungi for fungal/eukaryotic runs with suitable references, '
            'or mixed/all/none to disable domain filtering. Provide domain-specific --ref/--alt-ref, '
            '--baseline-fasta, --previous-review, and --anchors as needed.'
        ))
    performance_review_parser.add_argument('--phylum', required=False,
        help='Filter iTOL output to sequences assigned to this phylum (e.g. Bacillota). Does not affect novelty scoring.')
    performance_review_parser.add_argument('--target', required=False, default=None,
        help=(
            'FASTA of sequences to measure novelty against instead of the DB. '
            'Leave unset to use the full rolling partner-candidate collection, including this batch with self-hits removed (recommended).'
        ),
    )
    performance_review_parser.add_argument('--force-rebuild', '--rebuild-tree', dest='force_rebuild', action='store_true', default=False,
        help='Rebuild the entire tree from scratch even when an existing alignment is present. '
             'When combined with --previous-review, ignores the previous alignment and jointly estimates '
             'tree topology across all datasets. (--rebuild-tree is an alias for this flag.)')
    performance_review_parser.add_argument('--anchors', required=False, default=None,
        help='Custom reference anchor FASTA for tree scaffolding. By default uses the bundled 26-sequence rumen/gut anchor set.')
    performance_review_parser.add_argument('--threads', dest='threads', type=int, default=4,
        help='CPU threads for MAFFT and VSEARCH (default: 4).')
    performance_review_parser.add_argument('--tree-method', dest='tree_method',
        choices=['fasttree', 'iqtree', 'iqtree-fast'], default='fasttree',
        help=(
            'Phylogenetic tree-building backend (default: fasttree). '
            'fasttree: approximate ML, GTR+CAT. '
            'iqtree: full ML, GTR+G+I (recommended for production/publication runs). '
            'iqtree-fast: IQ-TREE 2 with -fast flag (good for exploratory incremental runs).'
        ),
    )
    performance_review_parser.add_argument(
        '--neighbourhood-format', dest='neighbourhood_format',
        choices=['png'], default='png',
        help='Image format for local phylogenetic-neighbourhood figures (PNG only).',
    )
    performance_review_parser.add_argument(
        '--pangenome-target', dest='pangenome_target', type=int, default=selection_target,
        help='Target number of committed genomes per GTDB species/local pangenome group (default: 9). Baseline isolates count because their genomes are available.',
    )
    performance_review_parser.add_argument(
        '--candidate-set-size', dest='candidate_set_size', type=int, default=selection_panel_size,
        help='Maximum ranked diversity-panel candidates per GTDB species/local clade (default: 9).',
    )
    performance_review_parser.add_argument(
        '--baseline-redundancy-identity', type=float, default=baseline_redundancy_identity,
        help='Exclude uncommitted candidates at or above this nearest cultured-baseline marker identity when coverage also passes its threshold (default: 99.8).',
    )
    performance_review_parser.add_argument(
        '--baseline-redundancy-min-query-coverage', type=float, default=baseline_redundancy_coverage,
        help='Minimum query coverage required before baseline near-identity exclusion is applied (default: 95).',
    )
    performance_review_parser.add_argument(
        '--baseline-extension-min-identity', type=float, default=baseline_extension_identity,
        help='Minimum cultured-baseline 16S identity for a same-species candidate to enter a baseline-pangenome extension (default: 98.65).',
    )
    performance_review_parser.add_argument(
        '--baseline-extension-min-query-coverage', type=float, default=baseline_extension_coverage,
        help='Minimum query coverage for admission to a baseline-pangenome extension (default: 95).',
    )
    performance_review_parser.add_argument('--user-colours', dest='user_colours', required=False,
        help='CSV file mapping sequence IDs to custom hex colours for iTOL (columns: id, colour).')
    performance_review_parser.add_argument(
        '--group-phyla', dest='group_phyla', action='append', default=None, metavar='SPEC',
        help=(
            'Collapse multiple phyla into one colour in iTOL legends. Repeatable. '
            'Formats: "archaea", "bacteria", "Bacillota,Bacillota_I", "Firmicutes:Bacillota,Bacillota_I".'
        ),
    )
    performance_review_parser.add_argument('--functional', dest='functional', required=False, default=None,
        help=(
            'TSV file mapping sequence IDs to functional attributes. '
            'Header row required; first column = sequence ID; subsequent columns = functional attributes. '
            'Generates one iTOL file per column (DATASET_BINARY / DATASET_SIMPLEBAR / DATASET_COLORSTRIP).'
        ),
    )
    performance_review_parser.add_argument('--draft-rumen-functions', dest='draft_rumen_functions',
        action='store_true', default=False,
        help=(
            'Auto-generate a draft rumen functional-group annotation from the output taxonomy. '
            'Writes rumen_functions_draft.tsv and itol_func_Rumen_Functional_Group.itol. '
            'The draft TSV can be edited and re-supplied via --functional in later runs.'
        ),
    )

    quarterly_review_parser = sub.add_parser(
        'quarterly-review',
        help='Select the next project-wide genome tranche after rolling Performance Reviews.',
        description=(
            'Reconsider every accumulated partner isolate in one auditable selection round.\n\n'
            'Residual nine-genome coverage gaps are filled first. Remaining budget is allocated\n'
            'across species using marginal tree/marker diversity, existing genome representation,\n'
            'cultured-baseline novelty, GTDB-reference context, and MWL evidence. Recommendations\n'
            'never change already_sequenced status. Use genome ANI/phylogenomics to validate\n'
            'strain-level diversity once genome data are available.'
        ),
        formatter_class=_Fmt,
    )
    add_project_config_option(quarterly_review_parser)
    quarterly_review_parser.add_argument('--db', required=True,
        help='Rolling BranchManager SQLite database containing all partner and baseline sequences.')
    quarterly_review_parser.add_argument('-o', '--out', required=True,
        help='Output directory for quarterly_review_summary.tsv, next_genome_set.tsv, and quarterly_review_manifest.tsv.')
    quarterly_review_parser.add_argument('--genome-budget', type=int, required=False,
        help='Number of new PRIMARY genomes to nominate in this round. This is intentionally explicit.')
    quarterly_review_parser.add_argument('--backups-per-primary', type=int, default=1,
        help='Number of extraction-failure/diversity backups to nominate per primary (default: 1).')
    quarterly_review_parser.add_argument('--pangenome-target', type=int, default=selection_target,
        help='Initial exact-GTDB-species/local-group coverage target before expansion (default: 9).')
    quarterly_review_parser.add_argument('--baseline-redundancy-identity', type=float, default=baseline_redundancy_identity,
        help='Exclude uncommitted candidates at or above this cultured-baseline marker identity (default: 99.8).')
    quarterly_review_parser.add_argument('--baseline-redundancy-min-query-coverage', type=float, default=baseline_redundancy_coverage,
        help='Minimum query coverage required for baseline redundancy exclusion (default: 95).')
    quarterly_review_parser.add_argument('--baseline-extension-min-identity', type=float, default=baseline_extension_identity,
        help='Minimum same-species cultured-baseline marker identity for baseline-pangenome extension membership (default: 98.65).')
    quarterly_review_parser.add_argument('--baseline-extension-min-query-coverage', type=float, default=baseline_extension_coverage,
        help='Minimum query coverage for baseline-pangenome extension membership (default: 95).')
    quarterly_review_parser.add_argument('--assessment', action='append', default=None, metavar='TSV',
        help='Import a full sequence_assessment.tsv before selection (repeatable; later files supersede earlier rows). Future Performance Reviews store snapshots automatically.')
    quarterly_review_parser.add_argument('--partner-metadata', '--sequencing-metadata', dest='partner_metadata', default=None,
        help='Optional cumulative metadata TSV/CSV used to refresh confirmed already_sequenced status before this round.')
    quarterly_review_parser.add_argument('--tree', default=None,
        help='Latest cumulative Newick tree. Enables marginal patristic-diversity ranking.')
    quarterly_review_parser.add_argument('--alignment', default=None,
        help='Latest cumulative MSA. Recomputes nearest currently available genome context after metadata updates.')
    quarterly_review_parser.add_argument('--from-dir', default=None,
        help='Latest Performance Review output directory; used to locate tree/current_tree.nwk when --tree is omitted.')
    quarterly_review_parser.add_argument('--round-id', default=None,
        help='Stable round label. Default: UTC timestamp such as quarterly_review_20260713T120000Z.')
    quarterly_review_parser.add_argument('--neighbourhood-format', choices=['png'], default='png',
        help='Image format for Quarterly Review local-neighbourhood figures (PNG only).')
    quarterly_review_parser.add_argument('--include-moderate-evidence', action='store_true', default=False,
        help='Allow MODERATE marker-evidence rows into the Quarterly Review. By default they are REVIEW only.')

    label_maker_parser = sub.add_parser(
        'label-maker',
        help='Label Maker: regenerate iTOL annotation files from stored taxonomy.',
        description=(
            'Re-generate all iTOL colour strips (phylum, family, genus, dataset membership)\n'
            'from the taxonomy already stored in the DB. Useful after changing --group-phyla\n'
            'options or after manually editing the database.'
        ),
        formatter_class=_Fmt,
    )
    label_maker_parser.add_argument('--db', required=True,
        help='Path to the BranchManager SQLite database.')
    label_maker_parser.add_argument('-o', '--out', required=True,
        help='Output directory where iTOL files will be written (normally a Filing Cabinet or Performance Review output directory).')
    label_maker_parser.add_argument('--include-datasets', required=False,
        help='Comma-separated list of dataset names to include (default: all datasets in the DB).')
    label_maker_parser.add_argument('--sequence-domain', '--organism-domain', dest='sequence_domain',
        choices=SEQUENCE_DOMAIN_CHOICES, default=None,
        help='Optional domain profile filter for regenerated outputs: bacteria, archaea, fungi, or mixed/all/none.')
    label_maker_parser.add_argument(
        '--group-phyla', dest='group_phyla', action='append', default=None, metavar='SPEC',
        help=(
            'Collapse multiple phyla into one colour in iTOL legends. Repeatable. '
            'Formats: "archaea", "bacteria", "Bacillota,Bacillota_I", "Firmicutes:Bacillota,Bacillota_I".'
        ),
    )
    label_maker_parser.add_argument('--functional', dest='functional', required=False, default=None,
        help=(
            'TSV file mapping sequence IDs to functional attributes. '
            'Generates one iTOL file per column (DATASET_BINARY / DATASET_SIMPLEBAR / DATASET_COLORSTRIP).'
        ),
    )
    label_maker_parser.add_argument('--draft-rumen-functions', dest='draft_rumen_functions',
        action='store_true', default=False,
        help='Auto-generate rumen functional-group iTOL annotation from stored taxonomy.',
    )

    org_chart_parser = sub.add_parser(
        'org-chart',
        help='Org Chart: build a focused tree and iTOL files for a specific taxon.',
        description=(
            'Extract all sequences matching a given taxon from the DB and build a focused\n'
            'phylogenetic tree for that group only.\n\n'
            'Fast path: if --from-dir points to a directory containing current_alignment.fasta\n'
            '(a Filing Cabinet or Performance Review output), sequences are sliced from the pre-built alignment and\n'
            'FastTree is run directly — no MAFFT re-alignment needed (~seconds for hundreds of seqs).\n\n'
            'Slow path: if no existing alignment is found, a full MAFFT + FastTree build is run.\n\n'
            'Taxon formats accepted:\n'
            '  archaea, bacteria          → all sequences at domain level\n'
            '  Bacillota, Bacteroidota    → phylum name (auto-detected)\n'
            '  p__Bacillota               → GTDB-prefixed phylum\n'
            '  f__Lachnospiraceae         → GTDB-prefixed family\n'
            '  g__Ruminococcus            → GTDB-prefixed genus\n'
        ),
        formatter_class=_Fmt,
    )
    org_chart_parser.add_argument('--db', required=True,
        help='Path to the BranchManager SQLite database.')
    org_chart_parser.add_argument('-o', '--out', required=True,
        help='Output directory for the focused-tree results.')
    org_chart_parser.add_argument('--taxon', required=True,
        help='Taxon to extract (see description above for accepted formats).')
    org_chart_parser.add_argument(
        '--rank', required=False, default='auto',
        choices=['auto', 'domain', 'd', 'phylum', 'p', 'class', 'c',
                 'order', 'o', 'family', 'f', 'genus', 'g', 'species', 's'],
        help='Taxonomic rank to filter on. Default "auto" detects the rank from the taxon name or prefix.',
    )
    org_chart_parser.add_argument('--from-dir', dest='from_dir', required=False, default=None,
        help=(
            'Existing Filing Cabinet or Performance Review output directory containing current_alignment.fasta. '
            'Enables the fast path (sequences extracted from the existing MSA; no re-alignment).'
        ),
    )
    org_chart_parser.add_argument('--ref', required=False,
        help='Reference FASTA for orientation correction (slow-path full build only).')
    org_chart_parser.add_argument('--anchors', required=False, default=None,
        help='Custom reference anchor FASTA. By default uses the bundled 26-sequence rumen/gut anchor set (NCBI RefSeq).')
    org_chart_parser.add_argument('--threads', type=int, default=4,
        help='CPU threads for FastTree / MAFFT (default: 4).')
    org_chart_parser.add_argument('--min-seqs', dest='min_seqs', type=int, default=3,
        help='Minimum sequences required to proceed with tree building (default: 3).')
    org_chart_parser.add_argument('--no-tree', dest='no_tree', action='store_true', default=False,
        help='Skip tree building; only write taxonomy TSV, sequence list, and iTOL colour files.')
    org_chart_parser.add_argument(
        '--group-phyla', dest='group_phyla', action='append', default=None, metavar='SPEC',
        help=(
            'Collapse multiple phyla into one colour in iTOL legends. Repeatable. '
            'Formats: "archaea", "bacteria", "Bacillota,Bacillota_I", "Firmicutes:Bacillota,Bacillota_I".'
        ),
    )
    org_chart_parser.add_argument('--functional', dest='functional', required=False, default=None,
        help=(
            'TSV file mapping sequence IDs to functional attributes. Header row required; first column = sequence ID; subsequent columns = functional attributes. '
            'Generates one iTOL file per column (DATASET_BINARY / DATASET_SIMPLEBAR / DATASET_COLORSTRIP).'
        ),
    )
    org_chart_parser.add_argument('--draft-rumen-functions', dest='draft_rumen_functions',
        action='store_true', default=False,
        help='Auto-generate rumen functional-group iTOL annotation from stored taxonomy.',
    )

    paper_trail_parser = sub.add_parser(
        'paper-trail',
        help='Paper Trail and Merge Meeting: convert AB1 reads and assemble primer sequences.',
        description=(
            'Process Sanger chromatogram reads before Performance Review.\n\n'
            'Inputs may be AB1/ABI files, FASTA, or FASTQ. AB1 files are base-called from '
            'PBAS/PCON tags, quality-trimmed, oriented by primer direction, and optionally '
            'assembled into one consensus 16S sequence per isolate. For example, 27F reads '
            'are kept forward and 907R reads are reverse-complemented before overlap assembly.\n\n'
            'If --sample-map is omitted, sequence IDs and primer names are '
            'inferred from filenames such as Iso001_27F.ab1 and Iso001_907R.ab1.'
        ),
        formatter_class=_Fmt,
    )
    sanger_studio_parser = sub.add_parser(
        'sanger-studio',
        help='Open the local browser-based Sanger chromatogram inspector and assembler.',
        description=(
            'Launch Sanger Studio. Select sample maps or chromatograms, choose an output '
            'directory, inspect reads, process batches, and complete cohort marker reviews '
            'from the browser interface. No input or output paths are required at startup.'
        ),
        formatter_class=_Fmt,
    )
    sanger_studio_parser.add_argument('--port', type=int, default=8765,
        help='Local web server port (default: 8765).')
    sanger_studio_parser.add_argument('--no-browser', action='store_true',
        help='Start the local server without opening a browser window.')


    paper_trail_parser.add_argument(
        '--input', nargs='+', required=False, default=[],
        help='AB1/ABI/FASTA/FASTQ files or directories to process. Directories are searched recursively by default. Optional when --sample-map lists the read files.',
    )
    paper_trail_parser.add_argument('-o', '--out', required=True,
        help='Output directory for assembled.fasta, read_qc.tsv, visual reports, and assembly_report.tsv.')
    paper_trail_parser.add_argument(
        '--sample-map', required=False, default=None,
        help=(
            'One CSV/TSV per partner batch. Prefer one row per read with sequence_id, read_file, '
            'primer, direction, and processing_mode. Multiple rows may share a sequence_id when '
            'different primer reads belong to one isolate. Relative paths are resolved next to the map.'
        ),
    )
    paper_trail_parser.add_argument(
        '--primer', dest='primers', action='append', default=None,
        help='Primer name to recognise in filenames (repeatable). Defaults include common 16S primers such as 27F, 907R, and 1492R.',
    )
    paper_trail_parser.add_argument(
        '--primer-sequence', dest='primer_sequences', action='append', default=None, metavar='NAME=SEQUENCE',
        help='Primer oligonucleotide used for confident leading-primer removal (repeatable). Built-in common 16S primers are used by default.',
    )
    paper_trail_parser.add_argument('--trim-primers', dest='trim_primers',
        action=argparse.BooleanOptionalAction, default=True,
        help='Remove confidently matched primer sequence from the retained read (default: enabled).')
    paper_trail_parser.add_argument('--secondary-peak-ratio', dest='secondary_peak_ratio', type=float, default=paper_trail_policy.secondary_peak_ratio,
        help='Secondary/called chromatogram peak ratio considered mixed (default: 0.33).')
    paper_trail_parser.add_argument('--max-mixed-peak-percent', dest='max_mixed_peak_percent', type=float, default=paper_trail_policy.max_mixed_peak_percent,
        help='Maximum percent retained high-quality positions with mixed peaks before read QC failure (default: 15).')
    paper_trail_parser.add_argument('--mixed-peak-min-quality', dest='mixed_peak_min_quality', type=int, default=paper_trail_policy.mixed_peak_min_quality,
        help='Minimum Phred score for a mixed-peak position to count (default: 20).')
    paper_trail_parser.add_argument('--screen-ref', dest='screen_ref', default=None,
        help='Optional marker reference FASTA for independent per-primer taxonomy concordance screening; disagreements require manual review and retain the read-level assignments.')
    paper_trail_parser.add_argument('--screen-taxa', dest='screen_taxa', default=None,
        help='Optional taxonomy table corresponding to --screen-ref; FASTA header taxonomy is used when omitted.')
    paper_trail_parser.add_argument('--threads', type=int, default=4,
        help='CPU threads for optional primer-read taxonomy screening (default: 4).')
    paper_trail_parser.add_argument('--min-quality', dest='min_quality', type=int, default=paper_trail_policy.min_quality,
        help='Phred cutoff for Mott-style end trimming (default: 20).')
    paper_trail_parser.add_argument('--min-length', dest='min_length', type=int, default=paper_trail_policy.min_final_length,
        help='Minimum final sequence length to write to assembled.fasta (default: 800 bp).')
    paper_trail_parser.add_argument('--min-read-length', dest='min_read_length', type=int, default=paper_trail_policy.min_read_length,
        help='Minimum trimmed read length allowed to contribute to an assembly or best-read decision (default: 300 bp).')
    paper_trail_parser.add_argument('--min-mean-quality', dest='min_mean_quality', type=float, default=paper_trail_policy.min_mean_quality,
        help='Minimum mean Phred score after trimming/masking for read and final QC (default: 20).')
    paper_trail_parser.add_argument('--mask-quality', dest='mask_quality', type=int, default=paper_trail_policy.mask_quality,
        help='Mask internal bases below this Phred score to N before assembly (default: 20).')
    paper_trail_parser.add_argument('--max-read-expected-errors', dest='max_read_expected_errors', type=float, default=paper_trail_policy.max_read_expected_errors,
        help='Maximum expected base-call errors allowed per retained read (default: 8).')
    paper_trail_parser.add_argument('--max-output-expected-errors', dest='max_output_expected_errors', type=float, default=paper_trail_policy.max_output_expected_errors,
        help='Maximum expected base-call errors allowed in the final isolate sequence (default: 5).')
    paper_trail_parser.add_argument('--warn-n-percent', dest='warn_n_percent', type=float, default=paper_trail_policy.warn_n_percent,
        help='Percent N above which a passing sequence requires manual review (default: 3.0).')
    paper_trail_parser.add_argument('--max-n-percent', dest='max_n_percent', type=float, default=paper_trail_policy.max_n_percent,
        help='Maximum percent N allowed after masking before QC failure (default: 5.0).')
    paper_trail_parser.add_argument('--warn-internal-lowq-run', dest='warn_internal_low_quality_run', type=int, default=paper_trail_policy.warn_internal_low_quality_run,
        help='Internal low-quality run length above which manual review is required (default: 5 bp).')
    paper_trail_parser.add_argument('--max-internal-lowq-run', dest='max_internal_low_quality_run', type=int, default=paper_trail_policy.max_internal_low_quality_run,
        help='Maximum internal low-quality/ambiguous run length before read failure (default: 20 bp).')
    paper_trail_parser.add_argument('--max-conflict-density', dest='max_conflict_density', type=float, default=paper_trail_policy.max_conflict_density,
        help='Maximum overlap conflicts per 100 final bases before final QC failure (default: 1.0).')
    paper_trail_parser.add_argument('--max-ambiguous-overlap-conflicts-without-review', dest='max_ambiguous_overlap_conflicts_without_review', type=int, default=paper_trail_policy.max_ambiguous_overlap_conflicts_without_review,
        help='Unresolved overlap conflicts tolerated as N before manual review is required (default: 2).')
    paper_trail_parser.add_argument('--quality-difference', dest='quality_difference', type=int, default=paper_trail_policy.quality_difference,
        help='Minimum Phred-scaled posterior odds required to resolve a conflicting overlap base (default: 10).')
    paper_trail_parser.add_argument('--allow-missing-quality', dest='allow_missing_quality',
        action='store_true', default=False,
        help='Allow AB1 reads missing PCON quality scores to pass with warnings. By default they fail QC.')
    paper_trail_parser.add_argument('--min-overlap', dest='min_overlap', type=int, default=paper_trail_policy.min_overlap,
        help='Minimum overlap length for assembling multiple primer reads (default: 80 bp).')
    paper_trail_parser.add_argument('--min-overlap-identity', dest='min_overlap_identity', type=float, default=paper_trail_policy.min_overlap_identity,
        help='Minimum overlap identity for assembly, 0-1 (default: 0.98).')
    paper_trail_parser.add_argument('--assemble', dest='assemble',
        action=argparse.BooleanOptionalAction, default=True,
        help='Assemble multiple reads per sequence_id when possible (default). Use --no-assemble to keep the best read.')
    paper_trail_parser.add_argument('--recursive', dest='recursive',
        action=argparse.BooleanOptionalAction, default=True,
        help='Search input directories recursively (default).')
    paper_trail_parser.add_argument('--max-report-image-height', dest='max_report_image_height',
        type=int, default=2400,
        help='Maximum height in pixels for each visual-report PNG; larger reports are split into numbered pages (minimum: 600, default: 2400).')

    replacement_parser = sub.add_parser(
        'replace-selected-markers',
        help='Replace selected future-genome markers with newly assembled Sanger 16S consensus sequences.',
        formatter_class=_Fmt,
    )
    replacement_parser.add_argument('--sample-map', required=True,
        help='CSV/TSV mapping selected isolate IDs to new forward/reverse Sanger AB1 reads.')
    replacement_parser.add_argument('--db', required=True,
        help='Existing BranchManager SQLite project database.')
    replacement_parser.add_argument('--primer', dest='primers', action='append', default=None)
    replacement_parser.add_argument('--primer-sequence', dest='primer_sequences', action='append', default=None,
        metavar='NAME=SEQUENCE')
    replacement_parser.add_argument('--min-overlap', type=int, default=paper_trail_policy.min_overlap)
    replacement_parser.add_argument('--min-overlap-identity', type=float, default=paper_trail_policy.min_overlap_identity)
    replacement_parser.add_argument('--threads', type=int, default=4)
    replacement_parser.add_argument('--allow-missing-quality', action='store_true', default=False)
    replacement_parser.add_argument('-o', '--out', required=True,
        help='Output directory for Sanger QC, replacement audit, and replacement FASTA.')

    from branchmanager.tools import (
        benchmark_selection as _tool_benchmark, build_anchor_fasta as _tool_anchors,
        id_sanitation as _tool_ids, metadata_to_diversity as _tool_diversity,
        primers_to_manifest as _tool_manifest, genomic_16s_selection as _tool_genomic_16s,
    )
    for _name, _module, _help in (
        ('primers-to-manifest', _tool_manifest,
         "Tool: convert a partner 'Primer file mapping' Excel sheet into an AB1 manifest.csv (with ID sanitation)."),
        ('metadata-to-diversity', _tool_diversity,
         'Tool: extract diversity-metadata columns from a partner Excel template into a TSV.'),
        ('benchmark-selection', _tool_benchmark,
         'Tool: benchmark isolate rankings against adjudicated labels (precision@k, recall@k, MAP).'),
        ('build-anchors', _tool_anchors,
         'Tool: download NCBI RefSeq 16S tree-anchor sequences into rumen-gut_reference_anchors.fasta.'),
        ('check-ids', _tool_ids,
         'Tool: check (and optionally sanitise) sequence/isolate IDs in a CSV/TSV table or FASTA.'),
        ('genomic-16s-selection', _tool_genomic_16s,
         'Tool: run Barrnap, VSEARCH against GTDB SSU representatives, optional GTDB-Tk concordance, and select one 16S representative per genome.'),

    ):
        _tool_parser = sub.add_parser(_name, help=_help, description=_help, formatter_class=_Fmt)
        _module.add_arguments(_tool_parser)
        _tool_parser.set_defaults(_tool=_module, _tool_parser=_tool_parser)

    mailroom_parser = sub.add_parser(
        'mailroom',
        help='Mailroom: validate a Sanger or PacBio partner delivery and build its batch map.',
        description=(
            'Sanger (the default) reconciles AB1/ABI chromatograms and writes ab1_map.tsv. '
            'PacBio reconciles one segmented HiFi FASTQ per isolate, writes pacbio_map.tsv, '
            'then runs PacBio 16S clustering and representative selection automatically.'
        ),
        formatter_class=_Fmt,
    )
    mailroom_parser.add_argument('--read-dir', required=True,
        help='Directory containing Sanger AB1/ABI or segmented PacBio HiFi FASTQ files for one partner batch.')
    mailroom_parser.add_argument('--manifest', required=True,
        help='Manifest file. Sanger requires read ID, isolate ID, and direction; PacBio requires isolate ID and FASTQ filename.')
    mailroom_parser.add_argument('--technology', choices=['sanger', 'pacbio'], default='sanger',
        help='Delivery technology (default: sanger). PacBio triggers FASTQ clustering after validation.')
    mailroom_parser.add_argument('--pacbio-library-prep', default='Kinnex 16S rRNA Kit',
        help='PacBio library-preparation protocol retained in provenance.')
    mailroom_parser.add_argument('--pacbio-threads', type=int, default=8,
        help='CPU threads for PacBio vsearch clustering during automatic kinnex import (default: 8).')
    mailroom_parser.add_argument('--pacbio-vsearch-strand', choices=['plus', 'both'], default='plus',
        help='vsearch strand mode for PacBio clustering (default: plus).')
    mailroom_parser.add_argument('--pacbio-vsearch-query-cov', type=float, default=0.0,
        help='Optional PacBio vsearch --query_cov threshold from 0 to 1; 0 disables the flag (default: 0).')
    mailroom_parser.add_argument('--pacbio-vsearch-maxaccepts', type=int, default=0,
        help='PacBio vsearch --maxaccepts value for clustering; 0 keeps all accepted hits (default: 0).')
    mailroom_parser.add_argument('--pacbio-vsearch-maxrejects', type=int, default=0,
        help='PacBio vsearch --maxrejects value for clustering; 0 keeps the vsearch default behaviour for no reject cap (default: 0).')
    mailroom_parser.add_argument('--dataset', required=True,
        help='Stable batch label written to every map row, for example UoG_01.')
    mailroom_parser.add_argument('--forward-primer', default=None,
        help='Confirmed primer name applied to forward rows lacking a supplied/embedded primer.')
    mailroom_parser.add_argument('--reverse-primer', default=None,
        help='Confirmed primer name applied to reverse rows lacking a supplied/embedded primer.')
    mailroom_parser.add_argument('--processing-mode', choices=['auto', 'assemble', 'best_read'], default='auto',
        help='Default read handling. Auto selects assemble for multi-read isolates and best_read for single-read isolates.')
    mailroom_parser.add_argument('--recursive', action=argparse.BooleanOptionalAction, default=True,
        help='Search the AB1 directory recursively (default: enabled).')
    mailroom_parser.add_argument('-o', '--out', required=True,
        help='Output directory for ab1_map.tsv, AB1 inventory, discrepancy report, and summary.')

    interview_parser = sub.add_parser(
        'interview',
        help='Interview: run standalone AB1 conversion, assembly, and QC from Mailroom output.',
        description=(
            'Run the same versioned AB1 QC policy used by Assistant, starting from a completed '
            'Mailroom output directory or its ab1_map.tsv. Interview writes converted and assembled '
            'FASTA files, QC tables, chromatogram and assembly PNGs, manual-review records, and '
            'resequencing recommendations. It does not read or modify a BranchManager project database.'
        ),
        formatter_class=_Fmt,
    )
    interview_parser.add_argument(
        '--mailroom', required=True,
        help='Mailroom output directory, or its ab1_map.tsv file.',
    )
    interview_parser.add_argument(
        '--primer-sequence', dest='primer_sequences', action='append', default=None,
        metavar='NAME=SEQUENCE',
        help='Confirmed primer oligonucleotide for leading-primer removal (repeatable).',
    )
    interview_parser.add_argument(
        '--trim-primers', dest='trim_primers', action=argparse.BooleanOptionalAction,
        default=True, help='Remove confidently matched leading primer sequences (default: enabled).',
    )
    interview_parser.add_argument(
        '--screen-ref', dest='screen_ref', default=None,
        help='Optional marker reference FASTA for independent per-read taxonomy concordance screening; disagreements trigger manual review and report each read assignment.',
    )
    interview_parser.add_argument(
        '--screen-taxa', dest='screen_taxa', default=None,
        help='Optional taxonomy table corresponding to --screen-ref.',
    )
    interview_parser.add_argument(
        '--threads', type=int, default=4,
        help='CPU threads for optional taxonomy screening (default: 4).',
    )
    interview_parser.add_argument(
        '--allow-missing-quality', action='store_true', default=False,
        help='Allow traces missing PCON quality scores to pass with warnings.',
    )
    interview_parser.add_argument(
        '--allow-mailroom-review', action='store_true', default=False,
        help='Process a REVIEW_REQUIRED Mailroom batch after an explicit review; FAIL batches are always rejected.',
    )
    interview_parser.add_argument(
        '--max-report-image-height', type=int, default=2400,
        help='Maximum PNG height before automatic pagination (minimum: 600, default: 2400).',
    )
    interview_parser.add_argument(
        '-o', '--out', required=True,
        help='Output directory for standalone AB1 QC, converted FASTA files, and visual reports.',
    )

    onboarding_parser = sub.add_parser(
        'onboarding',
        help='Onboarding: validate partner IDs, metadata, and either raw-read ownership or a supplied FASTA.',
        formatter_class=_Fmt,
    )
    onboarding_inputs = onboarding_parser.add_mutually_exclusive_group(required=True)
    onboarding_inputs.add_argument('--sample-map',
        help='CSV/TSV mapping isolate IDs to raw read files for an AB1/primer-read submission.')
    onboarding_inputs.add_argument('--fasta',
        help='Partner-supplied marker FASTA to validate without running Paper Trail/Merge Meeting.')
    onboarding_parser.add_argument('--partner-metadata', required=True,
        help='Cumulative project metadata CSV/TSV. Required for every submission and updated as selection/genome statuses change.')
    onboarding_parser.add_argument('--partner-id', default=None,
        help='Expected partner acronym for this submission, used to validate cumulative-ledger ownership.')
    onboarding_parser.add_argument('--dataset', default=None,
        help='Unique partner batch label, for example QUB_01 or UoG_02.')
    onboarding_parser.add_argument('--read-dir', default=None,
        help='Optional base directory for relative read paths in the sample map.')
    onboarding_parser.add_argument('--primer', dest='primers', action='append', default=None,
        help='Primer column/name to recognise (repeatable).')
    onboarding_parser.add_argument('-o', '--out', required=True,
        help='Output directory for the normalised submission and Onboarding report.')

    kinnex_parser = sub.add_parser(
        'pacbio-16s-import',
        help='PacBio 16S Import: cluster segmented HiFi FASTQs into one marker per isolate.',
        formatter_class=_Fmt,
    )
    kinnex_parser.add_argument('--pacbio-map', required=True,
        help='CSV/TSV with one row per isolate: sequence_id (or isolate_id) and file/fastq_file.')
    kinnex_parser.add_argument('-o', '--out', required=True,
        help='Output directory for representatives FASTA, marker-QC sidecar, and PacBio QC report.')
    kinnex_parser.add_argument('--min-read-length', type=int, default=1000)
    kinnex_parser.add_argument('--max-read-length', type=int, default=2000)
    kinnex_parser.add_argument('--min-mean-quality', type=float, default=20.0)
    kinnex_parser.add_argument('--cluster-identity', type=float, default=0.995)
    kinnex_parser.add_argument('--min-reads', type=int, default=20)
    kinnex_parser.add_argument('--min-dominant-fraction', type=float, default=0.80)
    kinnex_parser.add_argument('--threads', type=int, default=1,
        help='CPU threads for vsearch clustering (default: 1).')
    kinnex_parser.add_argument('--vsearch-strand', choices=['plus', 'both'], default='plus',
        help='vsearch strand mode for cluster_fast (default: plus).')
    kinnex_parser.add_argument('--vsearch-query-cov', type=float, default=0.0,
        help='Optional vsearch --query_cov threshold from 0 to 1; 0 disables the flag (default: 0).')
    kinnex_parser.add_argument('--vsearch-maxaccepts', type=int, default=0,
        help='vsearch --maxaccepts value for clustering; 0 keeps all accepted hits (default: 0).')
    kinnex_parser.add_argument('--vsearch-maxrejects', type=int, default=0,
        help='vsearch --maxrejects value for clustering; 0 keeps the vsearch default behaviour for no reject cap (default: 0).')
    kinnex_parser.add_argument('--library-prep', default='Kinnex 16S rRNA Kit',
        help='Library-preparation protocol to retain in PacBio QC provenance.')

    status_meeting_parser = sub.add_parser(
        'status-meeting',
        help='Status Meeting: import factual isolate lifecycle changes into the project ledger.',
        formatter_class=_Fmt,
    )
    status_meeting_parser.add_argument('--db', required=True)
    status_meeting_parser.add_argument('--input', required=True,
        help='CSV/TSV with sequence_id/status/detail, or the returned selection_summary.tsv with SABDecision filled as APPROVE, REJECT or DEFER.')
    status_meeting_parser.add_argument('-o', '--out', required=True)

    records_update_parser = sub.add_parser(
        'records-update',
        help='Records Update: import completed genome/QC/GTDB/ANI evidence.',
        formatter_class=_Fmt,
    )
    records_update_parser.add_argument('--db', required=True)
    records_update_parser.add_argument('--input', required=True,
        help='CSV/TSV with sequence_id, genome_id/accession, genome_status, QC, and optional GTDB/ANI fields.')
    records_update_parser.add_argument('--min-completeness', type=float, default=90.0,
        help='Minimum estimated genome completeness for automatic QC pass (default: 90).')
    records_update_parser.add_argument('--max-contamination', type=float, default=5.0,
        help='Maximum estimated contamination for automatic QC pass (default: 5).')
    records_update_parser.add_argument('-o', '--out', required=True)

    external_genomes_parser = sub.add_parser(
        'add-external-genomes',
        help='Register previously sequenced collection genomes, with optional 16S markers.',
        formatter_class=_Fmt,
    )
    add_project_config_option(external_genomes_parser)
    external_genomes_parser.add_argument('--db', required=True)
    external_genomes_parser.add_argument('--input', required=True,
        help='CSV/TSV with sequence_id, genome_id/accession, optional dataset/tier and genome-level taxonomy/QC fields.')
    external_genomes_parser.add_argument('--marker-fasta', default=None,
        help='Optional marker FASTA keyed by sequence_id; records without a marker remain genome-only.')
    external_genomes_parser.add_argument('--dataset', default='ExternalGenomeBaseline')
    external_genomes_parser.add_argument('--baseline-tier', choices=['priority', 'secondary'], default='secondary')
    external_genomes_parser.add_argument('-o', '--out', required=True)

    maldi_parser = sub.add_parser(
        'add-maldi-identifications',
        help='Register partner isolates identified by MALDI-TOF without marker sequences.',
        formatter_class=_Fmt,
    )
    maldi_parser.add_argument('--db', required=True)
    maldi_parser.add_argument('--input', required=True,
        help='CSV/TSV with sequence_id, partner_id, and taxon_name or taxon_id; optional rank, score, instrument, library, date, and notes.')
    maldi_parser.add_argument('--dataset', default='MALDI_TOF_Partner_Isolates')
    maldi_parser.add_argument('-o', '--out', required=True)

    partner_workbook_parser = sub.add_parser(
        'import-partner-workbook',
        help='Import isolate metadata and MALDI status from a partner .xlsm workbook.',
        formatter_class=_Fmt,
    )
    partner_workbook_parser.add_argument('--db', required=True)
    partner_workbook_parser.add_argument('--input', required=True,
        help='Partner .xlsm workbook with a Bulk Upload Template sheet.')
    partner_workbook_parser.add_argument('--sheet', default='Bulk Upload Template',
        help='Worksheet containing the partner isolate records.')
    partner_workbook_parser.add_argument('--dataset', default='PartnerWorkbook',
        help='Dataset label for newly created markerless isolate records.')
    partner_workbook_parser.add_argument('--genome-only', action='store_true',
        help='Register isolate IDs and partner metadata only; ignore all marker/16S fields. Genome assemblies and QC are imported separately.')
    partner_workbook_parser.add_argument('-o', '--out', required=True)

    exit_interview_parser = sub.add_parser(
        'exit-interview',
        help='Exit Interview: withdraw sequences from active project analyses with an audit trail.',
        description=(
            'Preview or apply removal of one or more sequences from the project database. '
            'Applying an Exit Interview removes active sequence, taxonomy, distance, assessment, '
            'selection, provenance, status, and metadata records while retaining an immutable '
            'tombstone. Baseline and genome-backed records are protected by default.'
        ),
        formatter_class=_Fmt,
    )
    exit_source = exit_interview_parser.add_mutually_exclusive_group(required=True)
    exit_source.add_argument('--sequence-id', action='append',
        help='Sequence/isolate ID to remove; repeat for multiple IDs.')
    exit_source.add_argument('--input',
        help='CSV/TSV, optionally gzipped, with sequence_id and reason columns.')
    exit_interview_parser.add_argument('--db', required=True)
    exit_interview_parser.add_argument('--reason', default='',
        help='Required reason for direct IDs, or fallback reason for input rows.')
    exit_interview_parser.add_argument('--apply', action='store_true', default=False,
        help='Apply the reviewed plan. Without this flag, Exit Interview is preview-only.')
    exit_interview_parser.add_argument('--allow-baseline', action='store_true', default=False,
        help='Explicitly allow removal from a registered cultured baseline.')
    exit_interview_parser.add_argument('--allow-genome-records', action='store_true', default=False,
        help='Explicitly allow removal of an isolate that has genome evidence.')
    exit_interview_parser.add_argument('--backup', action=argparse.BooleanOptionalAction, default=True,
        help='Write project_before_exit_interview.sqlite before applying (default: enabled).')
    exit_interview_parser.add_argument('-o', '--out', required=True)

    marker_evidence_migrate_parser = sub.add_parser(
        'marker-evidence-migrate',
        help='Create/backfill the source-aware marker evidence registry for selection database v1.',
        formatter_class=_Fmt,
    )
    marker_evidence_migrate_parser.add_argument('--db', required=True)
    marker_evidence_migrate_parser.add_argument('--force', action='store_true', default=False,
        help='Re-run the legacy backfill. Normally only run once after reviewing the existing records.')
    marker_evidence_migrate_parser.add_argument('-o', '--out', required=True)

    register_marker_evidence_parser = sub.add_parser(
        'register-marker-evidence',
        help='Register additional marker observations (e.g. Barrnap 16S) without replacing the active marker.',
        formatter_class=_Fmt,
    )
    register_marker_evidence_parser.add_argument('--db', required=True)
    register_marker_evidence_parser.add_argument('--fasta', required=True,
        help='FASTA of supplemental genomic 16S sequences; first header token maps through metadata fasta_id/genome_id.')
    register_marker_evidence_parser.add_argument('--metadata', required=True,
        help='TSV mapping fasta_id/genome_id to an existing BranchManager sequence_id, source/tier, and genome/16S taxonomy evidence.')
    register_marker_evidence_parser.add_argument('--source-type', default='GENOME_BARRNAP',
        choices=['GENOME_BARRNAP'])
    register_marker_evidence_parser.add_argument('--evidence-tier', default='MODERATE', choices=['MODERATE', 'LOW'])
    register_marker_evidence_parser.add_argument('-o', '--out', required=True)

    genome_16s_import_parser = sub.add_parser(
        'import-genome-16s',
        help='Register genomic-16s-selection outputs: genomes, GENOME_BARRNAP supplementary markers, GTDB-Tk concordance.',
        formatter_class=_Fmt,
    )
    genome_16s_import_parser.add_argument('--db', required=True)
    genome_16s_import_parser.add_argument('--genomic-dir', required=True,
        help='Output directory of genomic-16s-selection.')
    genome_16s_import_parser.add_argument('--genome-metadata', default=None,
        help='TSV mapping each genome_id to its stable isolate sequence_id; may also include partner_id, dataset, accession, completeness, contamination, genome_qc_pass, and genome_path.')
    genome_16s_import_parser.add_argument('--genome-id-isolate-id', action='store_true',
        help='Use genome IDs as isolate IDs only when this is an explicit, stable naming convention.')
    genome_16s_import_parser.add_argument('--dataset', default='Genomes')
    genome_16s_import_parser.add_argument('--partner-id', default='')
    genome_16s_import_parser.add_argument('-o', '--out', required=True)

    genome_evidence_parser = sub.add_parser(
        'import-genome-evidence',
        help='Import partner genome FASTAs with GTDB-Tk and CheckM2 evidence using partner culture IDs.',
        formatter_class=_Fmt,
    )
    genome_evidence_parser.add_argument('--db', required=True)
    genome_evidence_parser.add_argument('--workbook', required=True,
        help='Partner Bulk Upload Template workbook used to map culture IDs to BranchManager isolate IDs.')
    genome_evidence_parser.add_argument('--sheet', default='Bulk Upload Template')
    genome_evidence_parser.add_argument('--genome-dir', required=True,
        help='Directory containing partner-named .fa/.fna/.fasta files (optionally gzipped).')
    genome_evidence_parser.add_argument('--gtdbtk', required=True,
        help='GTDB-Tk bac120/ar53 summary TSV (delimiter is detected automatically).')
    genome_evidence_parser.add_argument('--checkm2', required=True,
        help='CheckM2 quality_report.tsv.')
    genome_evidence_parser.add_argument('--dataset', required=True,
        help='Candidate dataset used when the genome-only workbook was registered.')
    genome_evidence_parser.add_argument('--min-completeness', type=float, default=90.0)
    genome_evidence_parser.add_argument('--max-contamination', type=float, default=5.0)

    genome_results_parser = sub.add_parser(
        'import-genome-results',
        help='Import precomputed per-isolate genome QC and GTDB taxonomy results (no marker evidence required).',
        formatter_class=_Fmt,
    )

    genome_results_parser.add_argument('--db', required=True)
    genome_results_parser.add_argument(
        '--results-tsv',
        required=True,
        help='TSV with sequence_id, genome_id, genome_status, genome_qc_pass, completeness, contamination and gtdb_taxonomy.',
    )
    genome_results_parser.add_argument('--dry-run', action='store_true')
    genome_results_parser.add_argument(
        '--out',
        required=True,
        help='Output directory for staged genome-results import and audit reports.',
    )
    genome_evidence_parser.add_argument('--gtdbtk-version', default='not_recorded')
    genome_evidence_parser.add_argument('--checkm2-version', default='not_recorded')
    genome_evidence_parser.add_argument('--dry-run', action='store_true',
        help='Validate mappings and generate a report without changing the project database.')
    genome_evidence_parser.add_argument('-o', '--out', required=True)

    marker_review_apply_parser = sub.add_parser(
        'marker-review-apply',
        help='Persist manual marker-review decisions (APPROVE, APPROVE_FAIL_QC, APPROVE_UNVERIFIED, REJECT, RESEQUENCE) in the database.',
        formatter_class=_Fmt,
    )
    marker_review_apply_parser.add_argument('--db', required=True)
    marker_review_apply_parser.add_argument('--decisions', required=True,
        help='Completed marker_review_template.tsv (sequence_id, decision, reviewer, notes).')
    marker_review_apply_parser.add_argument('--dry-run', action='store_true', default=False)
    marker_review_apply_parser.add_argument('-o', '--out', required=True)

    annual_report_parser = sub.add_parser(
        'annual-report',
        help='Annual Report: create a point-in-time project overview, ledgers, and decision-change report.',
        formatter_class=_Fmt,
    )
    annual_report_parser.add_argument('--db', required=True)
    annual_report_parser.add_argument('-o', '--out', required=True)

    it_desk_parser = sub.add_parser(
        'it-desk',
        help='IT Desk: check dependencies, references, database integrity, and output location.',
        formatter_class=_Fmt,
    )
    it_desk_parser.add_argument('--db', default=None)
    it_desk_parser.add_argument('--ref', dest='references', action='append', default=[])
    it_desk_parser.add_argument('--tree-method', choices=['fasttree', 'iqtree', 'iqtree-fast'], default='fasttree')
    it_desk_parser.add_argument('-o', '--out', default='branchmanager_it_desk')
    it_desk_parser.add_argument('--strict', action='store_true', default=False,
        help='Return a non-zero exit status when any required check fails.')

    assistant_parser = sub.add_parser(
        'assistant',
        help='Assistant to the Branch Manager: run Onboarding through the Hiring Panel.',
        formatter_class=_Fmt,
    )
    add_project_config_option(assistant_parser)
    # Config profiles may provide a per-run input path through data.sample_map
    # (or assistant defaults), so enforce presence after config defaults load.
    assistant_input = assistant_parser.add_mutually_exclusive_group(required=False)
    assistant_input.add_argument('--sample-map', help='AB1/primer-read sample map; runs Paper Trail and Merge Meeting.')
    assistant_input.add_argument('--fasta', help='Partner-supplied marker FASTA; bypasses Paper Trail and Merge Meeting.')
    assistant_input.add_argument('--pacbio-map', help='One segmented PacBio HiFi FASTQ per isolate; clusters reads then continues as partner-sequenced markers.')
    assistant_input.add_argument('--cohort-manifest',
        help='Typed CSV/TSV (partner_id, dataset, input_type[sanger|partner_fasta|pacbio_hifi|genome_16s], modality-specific paths; optional genome_id_isolate_id=true) for a joint mixed-partner assessment.')
    assistant_input.add_argument('--cohort-resume', metavar='COHORT_OUT_DIR',
        help='Resume a prepared cohort: validate --marker-review decisions, assess all eligible markers jointly, then publish the database.')
    assistant_parser.add_argument('--cohort-prepare-only', action='store_true', default=False,
        help='With --cohort-manifest: run only the PREPARE stage (QC, review template, evidence manifest); no Hiring Panel, no database publish.')
    assistant_parser.add_argument(
        '--partner-metadata', required=True,
        help='Cumulative partner isolate catalogue. Its taxon/expected_taxon value is the authoritative partner-isolate identity; sequence-classifier disagreements are reported but do not replace the catalogue value or exclude the isolate.',
    )
    assistant_parser.add_argument('--partner-id', default=None,
        help='Expected partner acronym for validation, for example QUB or UoG.')
    assistant_parser.add_argument('--marker-source-type',
        choices=['auto', 'SANGER_AB1', 'PARTNER_FASTA', 'PACBIO_HIFI', 'MIXED'], default='auto',
        help='Evidence provenance for the marker imported by Assistant (default: infer from input mode).')
    assistant_parser.add_argument('--read-dir', default=None)
    assistant_parser.add_argument('--db', required=True)
    assistant_parser.add_argument('--dataset', required=True)
    assistant_parser.add_argument('--ref', required=True)
    assistant_parser.add_argument('--taxa', default=None)
    assistant_parser.add_argument('--ref-name', default='GTDB')
    assistant_parser.add_argument('--alt-ref', action='append', default=[])
    assistant_parser.add_argument('--alt-taxa', action='append', default=[])
    assistant_parser.add_argument('--alt-ref-name', action='append', default=[])
    assistant_parser.add_argument('--baseline-fasta', default=None)
    assistant_parser.add_argument('--baseline-dataset', default='Baseline')
    assistant_parser.add_argument('--baseline-tier', choices=['priority', 'secondary'], default=None)
    assistant_parser.add_argument('--baseline-taxa-assignments', default=None)
    assistant_parser.add_argument('--anchors', default=None,
        help='Custom reference anchor FASTA for tree scaffolding. By default uses the bundled 26-sequence rumen/gut anchor set; configure this option to use a collection-specific set.')
    assistant_parser.add_argument('--mwl', default=None)
    assistant_parser.add_argument('--sequence-domain', choices=SEQUENCE_DOMAIN_CHOICES, default='bacteria')
    assistant_parser.add_argument('--threads', type=int, default=4,
        help='CPU threads for vsearch and MAFFT (default: 4).')
    assistant_parser.add_argument('--tree-method', choices=['fasttree', 'iqtree', 'iqtree-fast'], default='fasttree')
    assistant_parser.add_argument('--pangenome-target', type=int, default=selection_target)
    assistant_parser.add_argument('--candidate-set-size', type=int, default=selection_panel_size)
    assistant_parser.add_argument('--baseline-redundancy-identity', type=float, default=baseline_redundancy_identity)
    assistant_parser.add_argument('--baseline-redundancy-min-query-coverage', type=float, default=baseline_redundancy_coverage)
    assistant_parser.add_argument('--baseline-extension-min-identity', type=float, default=baseline_extension_identity)
    assistant_parser.add_argument('--baseline-extension-min-query-coverage', type=float, default=baseline_extension_coverage)
    assistant_parser.add_argument('--primer', dest='primers', action='append', default=None)
    assistant_parser.add_argument('--primer-sequence', dest='primer_sequences', action='append', default=None,
        metavar='NAME=SEQUENCE', help='Primer sequence used for IUPAC-aware trimming; repeatable.')
    assistant_parser.add_argument('--trim-primers', dest='trim_primers',
        action=argparse.BooleanOptionalAction, default=True)
    assistant_parser.add_argument('--min-quality', type=int, default=paper_trail_policy.min_quality,
        help='Phred cutoff for Mott-style end trimming (default: 20).')
    assistant_parser.add_argument('--min-marker-length', dest='min_length', type=int, default=paper_trail_policy.min_final_length,
        help='Minimum final assembled sequence length in bp (default: 800).')
    assistant_parser.add_argument('--min-read-length', type=int, default=paper_trail_policy.min_read_length,
        help='Minimum trimmed primer-read length allowed to contribute to an assembly (default: 300 bp).')
    assistant_parser.add_argument('--min-mean-quality', type=float, default=paper_trail_policy.min_mean_quality,
        help='Minimum mean Phred quality after trimming (default: 20.0).')
    assistant_parser.add_argument('--mask-quality', type=int, default=paper_trail_policy.mask_quality,
        help='Mask internal bases below this Phred score to N (default: 20).')
    assistant_parser.add_argument('--max-read-expected-errors', type=float, default=paper_trail_policy.max_read_expected_errors,
        help='Maximum expected errors per retained read (default: 8.0).')
    assistant_parser.add_argument('--max-output-expected-errors', type=float, default=paper_trail_policy.max_output_expected_errors,
        help='Maximum expected errors in the final assembled sequence (default: 5.0).')
    assistant_parser.add_argument('--warn-n-percent', type=float, default=paper_trail_policy.warn_n_percent,
        help='Percent N above which a passing sequence requires manual review (default: 3.0).')
    assistant_parser.add_argument('--max-n-percent', type=float, default=paper_trail_policy.max_n_percent,
        help='Maximum percent N bases allowed after masking before QC failure (default: 5.0).')
    assistant_parser.add_argument('--warn-internal-lowq-run', dest='warn_internal_low_quality_run', type=int, default=paper_trail_policy.warn_internal_low_quality_run,
        help='Internal low-quality run length above which manual review is required (default: 5 bp).')
    assistant_parser.add_argument('--max-internal-lowq-run', dest='max_internal_low_quality_run', type=int, default=paper_trail_policy.max_internal_low_quality_run,
        help='Maximum internal low-quality/ambiguous run before QC failure (default: 20 bp).')
    assistant_parser.add_argument('--secondary-peak-ratio', type=float, default=paper_trail_policy.secondary_peak_ratio,
        help='Secondary/called peak ratio threshold for mixed-base detection (default: 0.33).')
    assistant_parser.add_argument('--max-mixed-peak-percent', type=float, default=paper_trail_policy.max_mixed_peak_percent,
        help='Maximum percent mixed-peak positions before read QC failure (default: 15.0).')
    assistant_parser.add_argument('--mixed-peak-min-quality', dest='mixed_peak_min_quality', type=int, default=paper_trail_policy.mixed_peak_min_quality,
        help='Minimum Phred score for a mixed-peak position to count (default: 20).')
    assistant_parser.add_argument('--max-conflict-density', type=float, default=paper_trail_policy.max_conflict_density,
        help='Maximum overlap conflicts per 100 final bases before QC failure (default: 1.0).')
    assistant_parser.add_argument('--max-ambiguous-overlap-conflicts-without-review', type=int, default=paper_trail_policy.max_ambiguous_overlap_conflicts_without_review,
        help='Unresolved overlap conflicts tolerated as N before manual review is required (default: 2).')
    assistant_parser.add_argument('--quality-difference', type=int, default=paper_trail_policy.quality_difference,
        help='Phred-scaled posterior odds required to resolve an overlap conflict (default: 10).')
    assistant_parser.add_argument('--allow-missing-quality', action='store_true', default=False,
        help='Allow traces missing PCON quality scores to pass with warnings.')
    assistant_parser.add_argument('--max-report-image-height', type=int, default=2400,
        help='Maximum Paper Trail visual-report PNG height before automatic pagination (minimum: 600, default: 2400).')
    assistant_parser.add_argument('--min-overlap', type=int, default=paper_trail_policy.min_overlap)
    assistant_parser.add_argument('--min-overlap-identity', type=float, default=paper_trail_policy.min_overlap_identity)
    assistant_parser.add_argument('--chimera-ref', default=None,
        help='Curated reference FASTA for UCHIME; defaults to the primary reference.')
    assistant_parser.add_argument('--skip-chimera-check', action='store_true', default=False,
        help='Audited override: continue with all markers marked for review.')
    assistant_parser.add_argument('--marker-qc', default=None,
        help='Marker-QC sidecar for --fasta submissions, such as a reviewed assembly_report.tsv.')
    assistant_parser.add_argument('--marker-review', default=None,
        help='CSV/TSV cohort review decisions (including APPROVE_FAIL_QC); during PREPARE these are carried into the regenerated template, and during RESUME they are validated and applied.')
    assistant_parser.add_argument('--accept-unverified-marker-qc', action='store_true', default=False,
        help='Permit processing a partner FASTA without QC provenance; this does not approve it for selection.')
    assistant_parser.add_argument('--unverified-override-authorizer', default='',
        help='Person authorizing --accept-unverified-marker-qc (required with the flag).')
    assistant_parser.add_argument('--unverified-override-reason', default='',
        help='Auditable reason for --accept-unverified-marker-qc (required with the flag).')
    assistant_parser.add_argument('--pacbio-min-read-length', type=int, default=1000)
    assistant_parser.add_argument('--pacbio-max-read-length', type=int, default=2000)
    assistant_parser.add_argument('--pacbio-min-mean-quality', type=float, default=20.0)
    assistant_parser.add_argument('--pacbio-cluster-identity', type=float, default=0.995)
    assistant_parser.add_argument('--pacbio-min-reads', type=int, default=20)
    assistant_parser.add_argument('--pacbio-min-dominant-fraction', type=float, default=0.80)
    assistant_parser.add_argument('--pacbio-threads', type=int, default=1)
    assistant_parser.add_argument('--pacbio-vsearch-strand', choices=['plus', 'both'], default='plus')
    assistant_parser.add_argument('--pacbio-vsearch-query-cov', type=float, default=0.0)
    assistant_parser.add_argument('--pacbio-vsearch-maxaccepts', type=int, default=0)
    assistant_parser.add_argument('--pacbio-vsearch-maxrejects', type=int, default=0)
    assistant_parser.add_argument('--pacbio-library-prep', default='Kinnex 16S rRNA Kit')
    assistant_parser.add_argument('-o', '--out', required=True)

    background_check_parser = sub.add_parser(
        'background-check',
        help='Background Check: pre-classify reference collections for reuse.',
        description=(
            'Classify one or more reference FASTA collections against a reference database\n'
            'using vsearch and save the results so the main pipeline can reuse them without\n'
            're-classifying on every run.\n\n'
            'Known dataset names (use with --dataset NAME=FASTA):\n'
            '  hungate16s / hungate   Hungate1000 curated rumen 16S sequences\n'
            '  silva                  SILVA ribosomal RNA database\n'
            '  rdp                    Ribosomal Database Project\n'
            '  homd                   Human Oral Microbiome Database\n'
            '  greengenes2 / gg2      GreenGenes2 16S reference\n'
            '  ncbi16s                NCBI 16S rRNA RefSeq collection\n'
            '  gtdb                   GTDB 16S representative sequences\n'
            '  (any other name)       Treated as a custom dataset\n\n'
            'Outputs written to --out:\n'
            '  {name}_classification.tsv   Full classification (human-readable)\n'
            '  {name}_taxonomy.tsv         Condensed per-dataset taxonomy\n'
            '  {name}_taxonomic_disagreement.tsv  High-quality hits with conflicting taxa\n'
            '  combined_taxonomy.tsv       All datasets merged (with Dataset column)\n'
            '  pipeline_taxonomy.tsv       All datasets merged; pass to --taxa-assignments\n'
            '  background_check_summary.txt  Plain-text summary with usage examples\n\n'
            'Example:\n'
            '  branchmanager background-check \\\n'
            '    --dataset hungate16s=/data/hungate.fasta \\\n'
            '    --dataset silva=/data/silva_16s.fasta \\\n'
            '    --ref /data/gtdb_ssu_reps.fna \\\n'
            '    --taxa /data/gtdb_taxonomy.tsv.gz \\\n'
            '    --threads 8 -o background_check_out/\n\n'
            'Then use the output in the Filing Cabinet:\n'
            '  branchmanager filing-cabinet --fasta hungate.fasta \\\n'
            '    --taxa-assignments background_check_out/pipeline_taxonomy.tsv \\\n'
            '    --db project.db --dataset Hungate -o filing_cabinet_out/'
        ),
        formatter_class=_Fmt,
    )
    background_check_parser.add_argument(
        '--dataset', dest='datasets', nargs='+', metavar='NAME FASTA',
        required=True,
        help=(
            'One or more NAME FASTA pairs, e.g.: '
            '--dataset Hungate /data/hun.fasta QUB /data/qub.fasta. '
            'Alternatively use NAME=FASTA pairs (repeatable): '
            '--dataset Hungate=/data/hun.fasta --dataset QUB=/data/qub.fasta. '
            'Known names: hungate16s, silva, rdp, homd, greengenes2, ncbi16s, gtdb. '
            'Any other name is treated as a custom dataset.'
        ),
    )
    background_check_parser.add_argument(
        '--ref', required=True,
        help='Reference FASTA (e.g. GTDB/SILVA reps) used by vsearch for classification.',
    )
    background_check_parser.add_argument(
        '--taxa', required=False, default=None,
        help='Taxonomy table matching IDs in --ref (TSV/CSV, optionally .gz: id<TAB>lineage or id,lineage). '
             'If omitted, taxonomy is parsed directly from reference FASTA headers.',
    )
    background_check_parser.add_argument(
        '-o', '--out', required=True,
        help='Output directory where all classification files will be written.',
    )
    background_check_parser.add_argument(
        '--threads', type=int, default=4,
        help='CPU threads for vsearch (default: 4).',
    )
    background_check_parser.add_argument(
        '--min-identity', dest='min_identity', type=float, default=0.80,
        help='Minimum vsearch alignment identity threshold 0–1 (default: 0.80 = 80%%).',
    )
    background_check_parser.add_argument(
        '--max-hits', dest='max_hits', type=int, default=10,
        help=(
            'Number of candidate hits vsearch collects per query '
            '(--maxaccepts / --maxhits).  The best hit by %% identity is '
            'selected after collection.  Higher values are more thorough but '
            'slower (default: 10).'
        ),
    )
    background_check_parser.add_argument(
        '--max-rejects', dest='max_rejects', type=int, default=256,
        help=(
            'vsearch --maxrejects: maximum number of non-matching candidate '
            'sequences examined before giving up on a query.  Raising this '
            'above vsearch\'s default of 32 helps classify ambiguous sequences '
            '(e.g. those with many N\'s) (default: 256).'
        ),
    )
    background_check_parser.add_argument(
        '--low-confidence-threshold', dest='low_confidence_threshold',
        type=float, default=0.97,
        help=(
            'Identity threshold below which a classified hit is flagged as low-confidence '
            'and written to *_low_confidence.tsv for manual review. '
            '0–1 (default: 0.97 = 97%%, the traditional species-level cutoff).'
        ),
    )

    dashboard_parser = sub.add_parser(
        'dashboard',
        help='Regenerate the interactive HTML dashboard from an existing run output directory.',
        description=(
            'Rebuild the interactive performance-review dashboard from the reports already present in a '
            'BranchManager run directory. This command does not rerun taxonomy, novelty, QC, phylogeny, '
            'pangenome or selection analyses. Optionally supply --db to overlay current SAB commitments and genome records without rerunning the assessment.\n\n'
            'Provide either the run output directory or its assessment/ subdirectory. The dashboard is '
            'written to the run output directory as performance_review_dashboard.html.'
        ),
        formatter_class=_Fmt,
    )
    dashboard_parser.add_argument(
        'output_dir',
        help='Existing BranchManager run output directory (or its assessment/ subdirectory).',
    )
    dashboard_parser.add_argument('--db', default=None, help='Optional project SQLite database; read-only overlay of current SAB decisions and registered genomes.')

    regenerate_parser = sub.add_parser(
        'regenerate-assessment',
        help='Regenerate selection, pangenome and dashboard reports from an existing assessment TSV.',
        formatter_class=_Fmt,
    )
    regenerate_parser.add_argument('--assessment-dir', required=True,
        help='Directory containing sequence_assessment.tsv (normally the assessment/ directory).')
    regenerate_parser.add_argument('--db', required=True,
        help='Existing BranchManager SQLite database; opened read-only for selection evidence.')
    regenerate_parser.add_argument('--tree', default=None,
        help='Existing Newick tree. If omitted, current_tree.nwk is searched in the assessment directory and its parent.')
    regenerate_parser.add_argument('--pangenome-target', type=int, default=selection_target)
    regenerate_parser.add_argument('--candidate-set-size', type=int, default=selection_panel_size)
    regenerate_parser.add_argument('--baseline-redundancy-identity', type=float, default=baseline_redundancy_identity)
    regenerate_parser.add_argument('--baseline-redundancy-min-query-coverage', type=float, default=baseline_redundancy_coverage)
    regenerate_parser.add_argument('--baseline-extension-min-identity', type=float, default=baseline_extension_identity)
    regenerate_parser.add_argument('--baseline-extension-min-query-coverage', type=float, default=baseline_extension_coverage)
    regenerate_parser.add_argument('--no-dashboard', action='store_true',
        help='Regenerate TSV reports only; do not rebuild the HTML dashboard.')

    # Additive validation toolkit. Existing workflow parsers and execution
    # paths remain unchanged.
    from branchmanager.validation import add_cli_arguments as _add_validation_arguments
    _add_validation_arguments(sub)

    return parser
