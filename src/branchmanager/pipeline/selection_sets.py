"""Rolling, clade-aware candidate sets for genome-sequencing decisions."""

from __future__ import annotations

import csv
import hashlib
import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from branchmanager.pipeline.evidence_policy import evidence_quality as _shared_evidence_quality
from branchmanager.pipeline.neighbourhood import pairwise_leaf_distances, _parse_tree, _mrca, _leaf_nodes
from branchmanager.taxonomy import domain_designation, infer_taxon_rank, normalise_taxon_name, parse_taxon_string


logger = logging.getLogger(__name__)

RANKS = [('s', 'species'), ('g', 'genus'), ('f', 'family'), ('o', 'order'),
         ('c', 'class'), ('p', 'phylum'), ('d', 'domain')]
MWL_STRENGTH = {
    'species': 5, 's': 5, 'genus': 4, 'g': 4, 'family': 3, 'f': 3,
    'order': 2, 'o': 2, 'class': 1, 'c': 1, 'phylum': 0, 'p': 0,
    'domain': 0, 'd': 0,
}
DIVERSITY_STRENGTH = {'high': 2, 'moderate': 1, 'low': 0, 'na': 0}
EVIDENCE_STRENGTH = {'HIGH': 2, 'MODERATE': 1, 'LOW': 0}
DEFAULT_PANGENOME_TARGET = 9
DEFAULT_CANDIDATE_SET_SIZE = 9
DEFAULT_BASELINE_REDUNDANCY_IDENTITY = 99.8
DEFAULT_BASELINE_REDUNDANCY_QUERY_COVERAGE = 95.0
DEFAULT_BASELINE_EXTENSION_MIN_IDENTITY = 98.65
DEFAULT_BASELINE_EXTENSION_MIN_QUERY_COVERAGE = 95.0
SELECTION_BLOCKING_PLACEMENT_FLAGS = frozenset({
    'NO_REFERENCE_HIT', 'NO_CLASSIFICATION',
    'CHIMERA', 'CHIMERA_CONFIRMED', 'CHIMERA_DETECTED',
    'CHIMERA_INDETERMINATE', 'CHIMERA_CHECK_SKIPPED', 'CHIMERA_NOT_RUN',
    'VERY_SHORT', 'HIGH_N_CONTENT', 'MARKER_QC_FAILED',
    'MARKER_QC_REVIEW_REQUIRED', 'MARKER_QC_UNVERIFIED',
    'MARKER_REPLACEMENT_QC_FAILED', 'MANUAL_APPROVAL_SEVERE_QC_OVERRIDE',
})


def _float(value, default: Optional[float] = None) -> Optional[float]:
    try:
        if value in (None, '', 'NA', 'None'):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _int(value, default: int = 0) -> int:
    try:
        if value in (None, '', 'NA', 'None'):
            return default
        return int(float(value))
    except (TypeError, ValueError):
        return default


def _selected(row: dict) -> bool:
    return str(row.get('selected_for_genome_sequencing') or '').strip().lower() == 'true'


def _available(row: dict) -> bool:
    return str(row.get('already_sequenced') or '').strip().lower() == 'true'


def _committed(row: dict) -> bool:
    return _selected(row) or _available(row)


def _evidence_quality(row: dict) -> str:
    return _shared_evidence_quality(row)


def _species(taxonomy: object) -> str:
    value = str(parse_taxon_string(str(taxonomy or '')).get('s') or '').strip()
    return value if value and value.lower() not in ('na', 'none', 'unclassified') else ''


def _assignment_status(row: dict) -> str:
    """Return the primary classifier status across supported assessment schemas.

    Regeneration reads the presentation-formatted sequence_assessment.tsv and
    converts headers to snake_case. Depending on the producer/version, the
    primary status is named classification_assignment_status,
    gtdb_assignment_status, or taxonomy_assignment_status.  Prefer an
    informative value over a placeholder so a default UNRESOLVED field cannot
    mask a populated GTDB status.
    """
    keys = (
        'classification_assignment_status', 'gtdb_assignment_status',
        'assignment_status_gtdb_r232', 'assignment_status_gtdb',
        'taxonomy_assignment_status', 'assignment_status',
    )
    values = [str(row.get(key) or '').strip().upper() for key in keys]
    missing = {'', 'NA', 'NONE', 'NAN', 'NOT_AVAILABLE'}
    # A generic unresolved field must not mask an informative database-specific
    # result (the common case in regenerated assessments).
    for value in values:
        if value and value not in missing | {'UNRESOLVED', 'NOT_CONFIDENT', 'NOT_ASSESSED'}:
            return value
    for value in values:
        if value and value not in missing:
            return value
    return 'UNRESOLVED'


def _taxonomic_concordance(row: dict) -> tuple[bool, str]:
    """Diagnostic classifier concordance; it is not an authority over the catalogue."""
    if _assignment_status(row) != 'CONFIDENT':
        return False, f'primary classifier assignment is {_assignment_status(row)}'
    primary_species = _species(row.get('taxonomy'))
    if not primary_species:
        return False, 'primary classifier assignment has no species rank'
    alt_tax_keys = sorted(key for key in row if key.startswith('alt_tax_'))
    for tax_key in alt_tax_keys:
        suffix = tax_key[len('alt_tax_'):]
        status = str(row.get(f'alt_assignment_status_{suffix}') or 'UNRESOLVED').upper()
        if status != 'CONFIDENT':
            return False, f'{suffix} classifier assignment is {status}'
        alt_species = _species(row.get(tax_key))
        if not alt_species:
            return False, f'{suffix} confident classifier assignment has no species rank'
        if normalise_taxon_name(alt_species) != normalise_taxon_name(primary_species):
            return False, f'{suffix} species {alt_species} conflicts with primary {primary_species}'
    return True, 'configured classifier assignments are confident and species-concordant'


def _classifier_database_concordance_status(row: dict) -> str:
    """Summarise cross-database classifier evidence without involving the catalogue."""
    if _assignment_status(row) != 'CONFIDENT':
        return 'NOT_CONFIDENT'
    concordant, reason = _taxonomic_concordance(row)
    if concordant:
        return 'CONCORDANT'
    if 'no species rank' in reason:
        return 'NOT_COMPARABLE_AT_SPECIES_RANK'
    if 'assignment is' in reason:
        return 'NOT_CONFIDENT_ALT_DATABASE'
    return 'DISCORDANT'


def _partner_taxon_lineage(row: dict) -> dict[str, str]:
    """Return rank/name pairs derived only from the partner catalogue value.

    A ranked lineage is preserved as supplied. For free-text catalogue values,
    conservative inference supports binomial species names and single-token
    genus names, plus unambiguous common higher-rank suffixes. No classifier
    result is used to fill missing catalogue ranks.
    """
    raw = str(
        row.get('partner_taxon') or row.get('partner_taxon_raw')
        or row.get('expected_taxon') or ''
    ).strip()
    if raw in {'', 'NA', 'None', 'unclassified'}:
        return {}
    parsed = parse_taxon_string(raw)
    if parsed:
        return parsed
    domain = domain_designation(raw)
    if domain:
        return {'d': domain}
    name = ' '.join(raw.split())
    rank = infer_taxon_rank(name, row.get('partner_taxon_rank'))
    if rank == 's' and len(name.split()) >= 2:
        parts = name.split()
        # The first two tokens are the binomial; strain/subspecies text is
        # retained in the original partner_taxon field, not promoted to rank.
        name = ' '.join(parts[:2])
    return {rank: name}




def _gtdb_lineage_key(lineage: dict, target_rank: str | None = None) -> str:
    """Stable composite GTDB taxon key at *target_rank* and broader ranks.

    The database stores GTDB taxonomy strings rather than separate GTDB taxon
    IDs. The lineage prefix through the requested rank disambiguates homonyms
    without splitting a genus into separate groups for each assigned species.
    """
    ranks = RANKS
    if target_rank:
        rank_codes = [rank for rank, _rank_name in RANKS]
        if target_rank in rank_codes:
            ranks = RANKS[rank_codes.index(target_rank):]
    return '|'.join(
        f'{rank}:{normalise_taxon_name(lineage.get(rank, ""))}'
        for rank, _rank_name in ranks if lineage.get(rank)
    )


def _pangenome_lineage(row: dict) -> dict[str, str]:
    """Return GTDB-derived lineage only; never substitute partner or marker taxonomy."""
    raw = str(row.get('gtdb_taxonomy') or row.get('gtdbtk_taxonomy') or '').strip()
    if not raw or raw.casefold() in {'na', 'none', 'unresolved', 'unclassified'}:
        return {}
    status = str(row.get('gtdb_assignment_status') or row.get('gtdb_taxonomy_assignment_status') or '').strip().upper()
    if status in {'FAIL', 'FAILED', 'UNRESOLVED', 'LOW_SIMILARITY', 'NOT_ASSESSED', 'NOT_CONFIDENT'}:
        return {}
    return parse_taxon_string(raw)


def _pangenome_lineage_source(row: dict) -> str:
    """Report whether a usable authoritative GTDB lineage is present."""
    return 'GTDB' if _pangenome_lineage(row) else 'UNRESOLVED'


def _partner_rank_taxon(row: dict) -> tuple[str, str, str]:
    lineage = _pangenome_lineage(row)
    for rank, rank_name in RANKS:
        value = str(lineage.get(rank) or '').strip()
        if value and value.lower() not in {'na', 'none', 'unclassified', 'uncultured'}:
            return rank, rank_name, value
    return '', '', ''


def _rank_value(row: dict, rank: str) -> str:
    """Read a rank for diagnostic matching; GTDB is preferred, marker taxonomy is secondary evidence."""
    lineage = _pangenome_lineage(row)
    if not lineage and row.get('taxonomy'):
        lineage = parse_taxon_string(str(row.get('taxonomy') or ''))
    if rank in lineage:
        return str(lineage[rank]).strip()
    if rank == 'g' and lineage.get('s'):
        return str(lineage['s']).strip().split()[0]
    return ''


def _group_for_row(row: dict) -> tuple[str, str, str]:
    """Group candidates by their most specific available GTDB rank and tree context."""
    figure = str(row.get('local_neighbourhood_figure') or '').strip()
    tree_context = figure if figure not in {'', 'NA', 'None', 'no_tree_context'} else 'no_tree_context'
    rank, rank_name, taxon = _partner_rank_taxon(row)
    if rank and taxon:
        lineage = _pangenome_lineage(row)
        # Include the lineage through this rank to avoid merging homonymous taxa.
        lineage_key = _gtdb_lineage_key(lineage, rank)
        return (
            f'gtdb:{rank}:{lineage_key}|clade:{tree_context}',
            f'GTDB {rank_name} + local tree clade',
            taxon,
        )
    status = _assignment_status(row)
    return (
        f'unresolved:{tree_context}:{status}',
        f'GTDB taxonomy unresolved; classifier status {status}',
        'UNASSIGNED [GTDB_TAXONOMY_UNRESOLVED]',
    )


def _stable_set_id(group_key: str, prefix: str = 'BMSET') -> str:
    digest = hashlib.sha1(group_key.encode('utf-8')).hexdigest()[:8].upper()
    return f'{prefix}_{digest}'


def _novel_looking(rows: Iterable[dict]) -> bool:
    for row in rows:
        evidence_class = str(row.get('baseline_evidence_class') or '').upper()
        if evidence_class in {
            'NOVEL TO BOTH BASELINES',
            'NOVEL TO HUNGATE / SECONDARY NOT ASSESSED',
            'NOVEL TO SECONDARY BASELINE / HUNGATE NOT ASSESSED',
            'HUNGATE GAP / SECONDARY COVERED',
        }:
            return True
        baseline = _float(row.get('nearest_identity'))
        baseline_coverage = _float(row.get('nearest_query_coverage'))
        reference = _float(row.get('reference_nearest_identity'))
        reference_coverage = _float(row.get('reference_nearest_query_coverage'))
        if (
            baseline is not None and baseline < 98.65
            and baseline_coverage is not None
            and baseline_coverage >= DEFAULT_BASELINE_EXTENSION_MIN_QUERY_COVERAGE
        ):
            return True
        if (
            reference is not None and reference < 98.65
            and reference_coverage is not None
            and reference_coverage >= DEFAULT_BASELINE_EXTENSION_MIN_QUERY_COVERAGE
        ):
            return True
    return False


def _first_present(*values):
    for value in values:
        if value not in (None, '', 'NA', 'None'):
            return value
    return None


def _baseline_hit_records(row: dict, *, include_combined: bool = True) -> list[dict]:
    records = [
        {
            'tier': 'priority',
            'label': 'Hungate',
            'hit': row.get('hungate_nearest_hit'),
            'identity': _float(row.get('hungate_nearest_identity')),
            'coverage': _float(row.get('hungate_nearest_query_coverage')),
            'taxonomy': row.get('hungate_nearest_hit_taxonomy'),
            'dataset': row.get('hungate_nearest_hit_dataset'),
        },
        {
            'tier': 'secondary',
            'label': 'SecondaryBaseline',
            'hit': row.get('secondary_baseline_nearest_hit'),
            'identity': _float(row.get('secondary_baseline_nearest_identity')),
            'coverage': _float(row.get('secondary_baseline_nearest_query_coverage')),
            'taxonomy': row.get('secondary_baseline_nearest_hit_taxonomy'),
            'dataset': row.get('secondary_baseline_nearest_hit_dataset'),
        },
    ]
    if include_combined:
        records.append({
            'tier': 'cultured_project',
            'label': 'CulturedProject',
            'hit': _first_present(row.get('cultured_project_nearest_hit'), row.get('cultured_rumen_nearest_hit'), row.get('nearest_hit')),
            'identity': _float(_first_present(row.get('cultured_project_nearest_identity'), row.get('cultured_rumen_nearest_identity'), row.get('nearest_identity'))),
            'coverage': _float(_first_present(row.get('cultured_project_nearest_query_coverage'), row.get('cultured_rumen_nearest_query_coverage'), row.get('nearest_query_coverage'))),
            'taxonomy': _first_present(row.get('cultured_project_nearest_hit_taxonomy'), row.get('cultured_rumen_nearest_hit_taxonomy'), row.get('nearest_hit_taxonomy')),
            'dataset': _first_present(row.get('cultured_project_nearest_hit_dataset'), row.get('cultured_rumen_nearest_hit_dataset'), row.get('nearest_hit_dataset')),
        })
    deduped = []
    seen = set()
    for record in records:
        hit = str(record.get('hit') or '').strip()
        if hit in {'', 'NA', 'None'}:
            continue
        key = (record['label'], hit)
        if key in seen:
            continue
        seen.add(key)
        deduped.append(record)
    return deduped


def _baseline_evidence_rank(row: dict) -> int:
    evidence_class = str(row.get('baseline_evidence_class') or '').strip().upper()
    ranks = {
        'NOVEL TO BOTH BASELINES': 4,
        'NOVEL TO HUNGATE / SECONDARY NOT ASSESSED': 3,
        'NOVEL TO SECONDARY BASELINE / HUNGATE NOT ASSESSED': 3,
        'HUNGATE GAP / SECONDARY COVERED': 3,
        'HUNGATE COVERED': 2,
        'SECONDARY COVERED / HUNGATE NOT ASSESSED': 2,
        'NO CULTURED BASELINE AVAILABLE': 1,
        'CULTURED BASELINE REDUNDANT': 0,
    }
    if evidence_class in ranks:
        return ranks[evidence_class]
    identity = _float(row.get('nearest_identity'))
    coverage = _float(row.get('nearest_query_coverage'))
    if identity is not None and coverage is not None and identity >= DEFAULT_BASELINE_REDUNDANCY_IDENTITY and coverage >= DEFAULT_BASELINE_REDUNDANCY_QUERY_COVERAGE:
        return 0
    if (
        identity is not None and coverage is not None
        and identity < DEFAULT_BASELINE_EXTENSION_MIN_IDENTITY
        and coverage >= DEFAULT_BASELINE_EXTENSION_MIN_QUERY_COVERAGE
    ):
        return 4
    return 1


def baseline_redundancy_evidence(
    row: dict,
    identity_threshold: float = DEFAULT_BASELINE_REDUNDANCY_IDENTITY,
    min_query_coverage: float = DEFAULT_BASELINE_REDUNDANCY_QUERY_COVERAGE,
) -> tuple[bool, Optional[float], Optional[float]]:
    """Return whether a candidate is near-identical to either cultured-project baseline tier."""
    hits = _baseline_hit_records(row, include_combined=False)
    if not hits:
        hits = _baseline_hit_records(row, include_combined=True)
    if not hits:
        identity = _float(_first_present(row.get('cultured_project_nearest_identity'), row.get('cultured_rumen_nearest_identity'), row.get('nearest_identity')))
        coverage = _float(_first_present(row.get('cultured_project_nearest_query_coverage'), row.get('cultured_rumen_nearest_query_coverage'), row.get('nearest_query_coverage')))
        redundant = bool(
            identity is not None and coverage is not None
            and identity >= float(identity_threshold)
            and coverage >= float(min_query_coverage)
        )
        return redundant, identity, coverage
    observed = [record for record in hits if record.get('identity') is not None]
    best = max(observed, key=lambda record: record['identity'], default=None)
    redundant_hits = [
        record for record in hits
        if record.get('identity') is not None
        and record.get('coverage') is not None
        and record['identity'] >= float(identity_threshold)
        and record['coverage'] >= float(min_query_coverage)
    ]
    if redundant_hits:
        winner = max(redundant_hits, key=lambda record: (record['identity'], record['coverage'] or 0.0))
        return True, winner.get('identity'), winner.get('coverage')
    return False, (best or {}).get('identity'), (best or {}).get('coverage')


def baseline_extension_evidence(
    row: dict,
    min_identity: float = DEFAULT_BASELINE_EXTENSION_MIN_IDENTITY,
    min_query_coverage: float = DEFAULT_BASELINE_EXTENSION_MIN_QUERY_COVERAGE,
) -> tuple[bool, str, str]:
    """Test whether a candidate can extend either cultured-project baseline tier."""
    candidate_species = _species(row.get('taxonomy'))
    if not candidate_species:
        return False, 'CANDIDATE_SPECIES_UNRESOLVED', 'candidate GTDB species is unresolved'
    hits = _baseline_hit_records(row, include_combined=False)
    if not hits:
        hits = _baseline_hit_records(row, include_combined=True)
    if not hits:
        return False, 'BASELINE_HIT_UNAVAILABLE', 'no cultured-baseline hit is available'

    failures: list[tuple[str, str]] = []
    for record in hits:
        hit = str(record.get('hit') or '').strip()
        baseline_species = _species(record.get('taxonomy'))
        label = str(record.get('label') or 'Baseline')
        status_prefix = {
            'Hungate': 'HUNGATE_BASELINE',
            'SecondaryBaseline': 'SECONDARY_BASELINE',
            'CulturedProject': 'BASELINE',
        }.get(label, 'BASELINE')
        eligible_status = {
            'Hungate': 'ELIGIBLE_HUNGATE_BASELINE_EXTENSION',
            'SecondaryBaseline': 'ELIGIBLE_SECONDARY_BASELINE_EXTENSION',
            'CulturedProject': 'ELIGIBLE_BASELINE_PANGENOME_EXTENSION',
        }.get(label, 'ELIGIBLE_BASELINE_PANGENOME_EXTENSION')
        if not baseline_species:
            failures.append((
                f'{status_prefix}_SPECIES_UNRESOLVED',
                f'{label} hit {hit} has unresolved species',
            ))
            continue
        if normalise_taxon_name(candidate_species) != normalise_taxon_name(baseline_species):
            failures.append((
                f'{status_prefix}_SPECIES_MISMATCH',
                f'candidate species {candidate_species} differs from {label} species {baseline_species}',
            ))
            continue
        identity = record.get('identity')
        if identity is None or identity < float(min_identity):
            observed = 'unavailable' if identity is None else f'{identity:.2f}%'
            failures.append((
                f'{status_prefix}_IDENTITY_BELOW_EXTENSION_THRESHOLD',
                f'{label} identity {observed} is below {float(min_identity):.2f}%',
            ))
            continue
        coverage = record.get('coverage')
        if coverage is None or coverage < float(min_query_coverage):
            observed = 'unavailable' if coverage is None else f'{coverage:.2f}%'
            failures.append((
                f'{status_prefix}_COVERAGE_BELOW_EXTENSION_THRESHOLD',
                f'{label} query coverage {observed} is below {float(min_query_coverage):.2f}%',
            ))
            continue
        return (
            True,
            eligible_status,
            f'exact species agreement with {hit}; {identity:.2f}% identity across {coverage:.2f}% of the query',
        )

    if failures:
        return False, failures[0][0], failures[0][1]
    return False, 'BASELINE_EXTENSION_CRITERIA_NOT_MET', 'no cultured-baseline tier met extension criteria'


def _evidence_tuple(row: dict) -> tuple:
    baseline = _float(row.get('nearest_identity'), 100.0)
    project = _float(row.get('project_nearest_identity'), 100.0)
    reference = _float(row.get('reference_nearest_identity'), 100.0)
    isolation = _float(row.get('phylo_isolation'), 0.0)
    diversity_rank = DIVERSITY_STRENGTH.get(str(row.get('diversity_metadata_strength') or '').strip().lower(), 0)
    diversity_score = _float(row.get('diversity_metadata_score'), 0.0)
    mwl = MWL_STRENGTH.get(str(row.get('mwl_matched_rank') or '').lower(), 0)
    mwl_score = _float(row.get('mwl_score'), 0.0)
    return (
        _baseline_evidence_rank(row),
        EVIDENCE_STRENGTH.get(_evidence_quality(row), 0),
        100.0 - float(baseline),
        100.0 - float(project),
        100.0 - float(reference),
        float(isolation),
        diversity_rank,
        float(diversity_score),
        mwl,
        float(mwl_score),
        str(row.get('id') or ''),
    )


def _ranking_reason(row: dict, opening: str, spread: Optional[float]) -> str:
    parts = [opening]
    if spread is not None:
        parts.append(f'marginal patristic distance {spread:.5f}')
    baseline = _float(row.get('nearest_identity'))
    if baseline is not None and baseline > 0:
        parts.append(f'nearest cultured baseline {baseline:.2f}%')
    project = _float(row.get('project_nearest_identity'))
    if project is not None and project > 0:
        parts.append(f'nearest project marker {project:.2f}%')
    diversity_strength = str(row.get('diversity_metadata_strength') or '').strip().upper()
    if diversity_strength in {'HIGH', 'MODERATE'}:
        reason = str(row.get('diversity_metadata_reason') or '').strip()
        parts.append(f'diversity metadata {diversity_strength.lower()} confidence' + (f' ({reason})' if reason else ''))
    mwl_rank = str(row.get('mwl_matched_rank') or '').strip()
    if mwl_rank not in {'', 'NA', 'None'}:
        parts.append(f'MWL {mwl_rank}-level context')
    return '; '.join(parts)


def _tree_leaf(row: dict) -> str:
    sid = str(row.get('id') or '')
    representative = str(row.get('cluster_representative') or '')
    if str(row.get('in_tree') or '').lower() == 'yes':
        return sid
    return representative if representative not in ('', 'self', 'duplicate', 'N/A', 'NA') else sid


def _distance(left: str, right: str, distances: Dict[tuple[str, str], float]) -> Optional[float]:
    return distances.get(tuple(sorted((left, right))))


def _tree_reference_rows(db, tree_path, assessment_rows: List[dict]) -> list[dict]:
    """Load GTDB taxonomy for non-candidate reference leaves in the tree."""
    if db is None or not tree_path or not Path(tree_path).exists():
        return []
    try:
        _root, _parents, leaves = _parse_tree(tree_path)
        with db.connect() as conn:
            records = conn.execute(
                'SELECT sequence_id, gtdb_taxonomy, genome_qc_pass FROM genome_records '
                'WHERE sequence_id IS NOT NULL'
            ).fetchall()
        assessed_ids = {str(row.get('id')) for row in assessment_rows}
        result = []
        for sid, taxonomy, qc_pass in records:
            sid = str(sid or '')
            if (not sid or sid in assessed_ids or sid not in leaves or not taxonomy
                    or str(qc_pass).strip().casefold() not in {'1', 'true', 'yes', 'pass', 'passed'}):
                continue
            result.append({
                'id': sid, 'gtdb_taxonomy': str(taxonomy),
                'gtdb_assignment_status': 'CONFIDENT',
                '_tree_taxonomy_source': 'GTDB_REFERENCE_GENOME',
            })
        return result
    except Exception:
        return []


def _tree_taxon_clade_agrees(tree_path, members: List[dict], all_rows: List[dict], rank: str, taxon: str) -> tuple[bool, str]:
    """Check MRCA exclusivity against other catalogue taxa and confident baseline assignments at the same rank."""
    if not tree_path or not Path(tree_path).exists():
        return False, 'current phylogenetic tree is unavailable'
    try:
        _root, parents, leaves = _parse_tree(tree_path)
        member_leaf_ids = [_tree_leaf(row) for row in members]
        missing_leaf_ids = sorted({leaf_id for leaf_id in member_leaf_ids if not leaf_id or leaf_id not in leaves})
        if missing_leaf_ids:
            preview = ';'.join(missing_leaf_ids[:10])
            suffix = '…' if len(missing_leaf_ids) > 10 else ''
            return False, f'{len(missing_leaf_ids)} candidate member(s) are absent from the tree: {preview}{suffix}'
        target_leaf_ids = set(member_leaf_ids)
        if not target_leaf_ids:
            return False, 'candidate GTDB taxon has no resolvable tree leaves'
        mrca = _mrca([leaves[leaf_id] for leaf_id in sorted(target_leaf_ids)], parents)
        descendant_ids = {node.name.strip() for node in _leaf_nodes(mrca) if node.name}
        conflicting = set()
        for row in all_rows:
            leaf_id = _tree_leaf(row)
            if leaf_id not in descendant_ids or leaf_id in target_leaf_ids:
                continue
            other_lineage = _pangenome_lineage(row)
            other = str(other_lineage.get(rank) or '').strip()
            if not other and rank == 'g' and other_lineage.get('s'):
                other = str(other_lineage['s']).split()[0]
            if other and normalise_taxon_name(other) != normalise_taxon_name(taxon):
                conflicting.add(other)
        if conflicting:
            return False, f'tree MRCA also contains other assigned GTDB {rank} taxa: ' + ';'.join(sorted(conflicting))
        return True, f'tree MRCA contains no conflicting assigned {rank} taxa at this rank'
    except Exception as exc:
        return False, f'tree clade concordance could not be assessed: {exc}'


def _farthest_first(rows: List[dict], anchors: List[str], distances: Dict[tuple[str, str], float]) -> List[tuple[dict, Optional[float]]]:
    remaining = list(rows)
    chosen: List[tuple[dict, Optional[float]]] = []
    selected_leaves = [anchor for anchor in anchors if anchor]
    while remaining:
        ranked = []
        for row in remaining:
            leaf = _tree_leaf(row)
            observed = [
                value for anchor in selected_leaves
                if (value := _distance(leaf, anchor, distances)) is not None
            ]
            spread = min(observed) if observed else None
            ranked.append((spread is not None, spread or 0.0, _evidence_tuple(row), row))
        has_spread, spread, _, winner = max(ranked, key=lambda item: (item[0], item[1], item[2]))
        chosen.append((winner, spread if has_spread else None))
        selected_leaves.append(_tree_leaf(winner))
        remaining.remove(winner)
    return chosen


def _committed_ids_by_rank(db) -> dict[tuple[str, str], list[str]]:
    """Index committed genomes exclusively by their stored GTDB taxonomy."""
    if db is None:
        return {}
    try:
        with db.connect() as conn:
            tables = {str(row[0]) for row in conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            ).fetchall()}
            if 'genome_records' not in tables:
                return {}
            records = conn.execute(
                "SELECT sequence_id, genome_qc_pass, gtdb_taxonomy FROM genome_records "
                "WHERE sequence_id IS NOT NULL"
            ).fetchall()
    except Exception:
        return {}
    indexed = defaultdict(set)
    for sid, qc_pass, taxonomy in records:
        sid = str(sid or '').strip()
        if not sid or str(qc_pass).strip().casefold() not in {'1', 'true', 'yes', 'pass', 'passed'}:
            continue
        lineage = parse_taxon_string(str(taxonomy or ''))
        for rank, _rank_name in RANKS:
            value = str(lineage.get(rank) or '').strip()
            if not value and rank == 'g' and lineage.get('s'):
                value = str(lineage['s']).split()[0]
            if value:
                indexed[(rank, normalise_taxon_name(value), _gtdb_lineage_key(lineage, rank))].add(sid)
    return {key: sorted(values) for key, values in indexed.items()}


def write_pangenome_hierarchy_tsv(
    assessment_rows: List[dict], output_path: str | Path, *, tree_path=None,
    db=None, pangenome_target: int = DEFAULT_PANGENOME_TARGET,
) -> str:
    """Write independent pangenome assessments at every GTDB-supported rank."""
    target = max(1, int(pangenome_target or DEFAULT_PANGENOME_TARGET))
    grouped = defaultdict(list)
    rank_names = dict(RANKS)
    pangenome_rows = list(assessment_rows)
    # Genome-only isolates have no marker sequence/tree leaf, but their genome
    # taxonomy and QC are valid inputs to genome representation counts.
    if db is not None:
        try:
            with db.connect() as conn:
                genome_rows = conn.execute(
                    "SELECT gr.sequence_id, s.dataset, gr.genome_id, gr.genome_status, "
                    "gr.genome_qc_pass, gr.completeness, gr.contamination, gr.gtdb_taxonomy "
                    "FROM genome_records gr LEFT JOIN sequences s ON s.id=gr.sequence_id"
                ).fetchall()
            assessment_ids = {str(row.get('id') or row.get('sequence_id') or '') for row in pangenome_rows}
            for sid, dataset, genome_id, status, qc_pass, completeness, contamination, taxonomy in genome_rows:
                sid = str(sid or '')
                qc_passed = str(qc_pass).strip().casefold() in {'1', 'true', 'yes', 'pass', 'passed'}
                if not sid or sid in assessment_ids or not taxonomy or not qc_passed:
                    continue
                pangenome_rows.append({
                    'id': sid, 'dataset': str(dataset or ''),
                    'taxonomy': str(taxonomy or ''),
                    'gtdb_taxonomy': str(taxonomy or ''),
                    'gtdb_assignment_status': 'CONFIDENT',
                    'classification_assignment_status': 'PROVISIONAL',
                    'already_sequenced': 'true', 'selected_for_genome_sequencing': 'false',
                    'local_neighbourhood_figure': 'no_tree_context',
                    'genome_id': str(genome_id or ''), 'genome_qc_pass': qc_pass,
                    'genome_completeness': completeness, 'genome_contamination': contamination,
                    'evidence_source': 'genome_only',
                })
        except Exception as exc:
            logger.warning('[PANGENOME] Could not load genome-only records: %s', exc)
    for row in pangenome_rows:
        lineage = _pangenome_lineage(row)
        lineage_source = _pangenome_lineage_source(row)
        figure = str(row.get('local_neighbourhood_figure') or '').strip()
        context = figure if figure not in {'', 'NA', 'None', 'no_tree_context'} else 'no_tree_context'
        for rank, rank_name in RANKS:
            value = str(lineage.get(rank) or '').strip()
            if not value and rank == 'g' and lineage.get('s'):
                value = str(lineage['s']).split()[0]
            if value and value.lower() not in {'na','none','unclassified','uncultured'}:
                lineage_key = _gtdb_lineage_key(lineage, rank)
                grouped[(rank, normalise_taxon_name(value), context, lineage_source, lineage_key)].append(row)
    committed = _committed_ids_by_rank(db)
    tree_context_rows = _tree_reference_rows(db, tree_path, pangenome_rows)
    all_tree_rows = list(pangenome_rows) + tree_context_rows
    output = []
    for (rank, _normalised, context, lineage_source, _lineage_key), members in sorted(grouped.items()):
        taxon = str(_pangenome_lineage(members[0]).get(rank) or '').strip()
        member_leaf_ids = [_tree_leaf(row) for row in members]
        try:
            tree_leaves = _parse_tree(tree_path)[2] if tree_path and Path(tree_path).exists() else {}
        except Exception:
            tree_leaves = {}
        missing_tree_leaves = sorted({leaf for leaf in member_leaf_ids if leaf in {'', 'NA', 'None'} or leaf not in tree_leaves})
        leaves_ok = context != 'no_tree_context' and bool(tree_leaves) and not missing_tree_leaves
        tree_ok = False
        tree_reason = 'tree/local-clade context unavailable'
        if leaves_ok:
            # The generic checker also verifies all members are represented
            # and excludes competing GTDB assignments at this rank.
            tree_ok, tree_reason = _tree_taxon_clade_agrees(
                tree_path, members, all_tree_rows, rank, taxon,
            )
        committed_ids = committed.get((rank, _normalised, _lineage_key), [])
        # Pangenome targets are defined at GTDB rank/taxon, not at the
        # individual local-tree-clade level.  Candidate-set reporting already
        # counts committed isolates across every clade for this taxon; use the
        # same universe here so hierarchy rows cannot report a different gap
        # merely because a committed isolate belongs to another neighbourhood.
        current_committed = set()
        for candidate in pangenome_rows:
            candidate_lineage = _pangenome_lineage(candidate)
            candidate_value = str(candidate_lineage.get(rank) or '').strip()
            if not candidate_value and rank == 'g' and candidate_lineage.get('s'):
                candidate_value = str(candidate_lineage['s']).split()[0]
            if (
                candidate_value
                and normalise_taxon_name(candidate_value) == _normalised
                and _gtdb_lineage_key(candidate_lineage, rank) == _lineage_key
                and _committed(candidate)
            ):
                candidate_id = str(candidate.get('id') or candidate.get('sequence_id') or '').strip()
                if candidate_id:
                    current_committed.add(candidate_id)
        all_committed = sorted(set(committed_ids) | current_committed)
        # Formal targets/gaps require a usable GTDB lineage and validated tree context.
        eligible = bool(taxon and lineage_source == 'GTDB' and leaves_ok and tree_ok)
        if lineage_source != 'GTDB':
            tree_reason = f'{lineage_source.lower().replace("_", " ")} lineage is exploratory only; usable GTDB lineage is not supplied'
        elif missing_tree_leaves:
            tree_reason = f'{len(missing_tree_leaves)} member(s) are absent from the phylogenetic tree'
        count = len(all_committed)
        disagreements = []
        unresolved = []
        for row in members:
            classifier_value = str(parse_taxon_string(str(row.get('taxonomy') or '')).get(rank) or '').strip()
            if not classifier_value:
                unresolved.append(str(row.get('id')))
            elif normalise_taxon_name(classifier_value) != _normalised:
                disagreements.append(str(row.get('id')))
        output.append({
            'PangenomeRank': rank_names.get(rank, rank).upper(),
            'PangenomeTaxon': taxon or 'NA',
            'TaxonomyProvenance': lineage_source,
            'LocalTreeContext': context,
            'CandidateCount': len(members),
            'CommittedGenomeCount': count,
            'CommittedGenomeIDs': ';'.join(all_committed) or 'None',
            'PangenomeTarget': target if eligible else 'NA',
            'PangenomeGap': max(0, target-count) if eligible else 'NA',
            'PangenomeEligible': 'Yes' if eligible else 'No',
            'EligibilityReason': tree_reason if eligible else f'GTDB rank/taxon or reference tree context unresolved; {tree_reason}',
            'ClassifierMarkerGTDBDisagreementCount': len(disagreements),
            'ClassifierMarkerGTDBDisagreementIDs': ';'.join(sorted(disagreements)) or 'None',
            'ClassifierUnresolvedCount': len(unresolved),
            'ClassifierUnresolvedIDs': ';'.join(sorted(unresolved)) or 'None',
            'CandidateIDs': ';'.join(sorted(str(row.get('id')) for row in members)),
        })
    fields = [
        'PangenomeRank','PangenomeTaxon','TaxonomyProvenance','LocalTreeContext','CandidateCount',
        'CommittedGenomeCount','CommittedGenomeIDs','PangenomeTarget','PangenomeGap',
        'PangenomeEligible','EligibilityReason','ClassifierMarkerGTDBDisagreementCount',
        'ClassifierMarkerGTDBDisagreementIDs','ClassifierUnresolvedCount',
        'ClassifierUnresolvedIDs','CandidateIDs',
    ]
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(output)
    return str(path)


def build_sequencing_sets(
    assessment_rows: List[dict],
    output_path: str | Path,
    *,
    tree_path: Optional[str | Path] = None,
    db=None,
    pangenome_target: int = DEFAULT_PANGENOME_TARGET,
    candidate_set_size: int = DEFAULT_CANDIDATE_SET_SIZE,
    baseline_redundancy_identity: float = DEFAULT_BASELINE_REDUNDANCY_IDENTITY,
    baseline_redundancy_min_query_coverage: float = DEFAULT_BASELINE_REDUNDANCY_QUERY_COVERAGE,
    baseline_extension_min_identity: float = DEFAULT_BASELINE_EXTENSION_MIN_IDENTITY,
    baseline_extension_min_query_coverage: float = DEFAULT_BASELINE_EXTENSION_MIN_QUERY_COVERAGE,
) -> str:
    """Assign rolling candidate panels using GTDB rank and local tree context."""
    target = max(1, int(pangenome_target or DEFAULT_PANGENOME_TARGET))
    set_size = max(target, int(candidate_set_size or DEFAULT_CANDIDATE_SET_SIZE))
    redundancy_identity = float(baseline_redundancy_identity)
    redundancy_coverage = float(baseline_redundancy_min_query_coverage)
    extension_identity = float(baseline_extension_min_identity)
    extension_coverage = float(baseline_extension_min_query_coverage)
    if not 0.0 <= redundancy_identity <= 100.0:
        raise ValueError('baseline redundancy identity must be between 0 and 100')
    if not 0.0 <= redundancy_coverage <= 100.0:
        raise ValueError('baseline redundancy minimum query coverage must be between 0 and 100')
    if not 0.0 <= extension_identity <= redundancy_identity:
        raise ValueError('baseline extension minimum identity must be between 0 and the redundancy identity threshold')
    if not 0.0 <= extension_coverage <= 100.0:
        raise ValueError('baseline extension minimum query coverage must be between 0 and 100')
    grouped: Dict[str, List[dict]] = defaultdict(list)
    group_metadata = {}
    for row in assessment_rows:
        key, basis, taxon = _group_for_row(row)
        grouped[key].append(row)
        group_metadata[key] = (basis, taxon)

    committed_by_rank = _committed_ids_by_rank(db)
    tree_context_rows = _tree_reference_rows(db, tree_path, assessment_rows)
    all_tree_rows = list(assessment_rows) + tree_context_rows
    leaf_ids = [_tree_leaf(row) for row in assessment_rows]
    leaf_ids.extend(sid for committed_ids in committed_by_rank.values() for sid in committed_ids)
    distances = {}
    tree_leaf_ids = set()
    if tree_path and Path(tree_path).exists():
        try:
            distances = pairwise_leaf_distances(tree_path, leaf_ids)
            tree_leaf_ids = {
                left for (left, right) in distances if left == right
            }
        except Exception:
            distances = {}
            tree_leaf_ids = set()
    output_rows = []
    for group_key in sorted(grouped):
        members = grouped[group_key]
        basis, taxon = group_metadata[group_key]
        rank, rank_name, partner_taxon = _partner_rank_taxon(members[0])
        gtdb_taxon_present = bool(rank and partner_taxon)
        member_taxa_agree = all(
            _partner_rank_taxon(row)[0] == rank
            and normalise_taxon_name(_partner_rank_taxon(row)[2]) == normalise_taxon_name(partner_taxon)
            for row in members
        ) if gtdb_taxon_present else False
        classifier_disagreements = [
            str(row.get('id')) for row in members
            if partner_taxon
            and parse_taxon_string(str(row.get('taxonomy') or '')).get(rank)
            and normalise_taxon_name(parse_taxon_string(str(row.get('taxonomy') or '')).get(rank))
                != normalise_taxon_name(partner_taxon)
        ]
        concordance_reason = (
            f'{len(classifier_disagreements)} marker-classifier/GTDB disagreement(s); GTDB retained as pangenome authority'
            if classifier_disagreements else 'usable GTDB taxon is available'
        )
        tree_context_consistent = bool(tree_path and Path(tree_path).exists()) and all(
            str(row.get('local_neighbourhood_figure') or '').strip()
            not in {'', 'NA', 'None', 'no_tree_context'}
            and _tree_leaf(row) in tree_leaf_ids
            for row in members
        )
        tree_clade_agrees, tree_clade_reason = (
            _tree_taxon_clade_agrees(tree_path, members, all_tree_rows, rank, partner_taxon)
            if tree_context_consistent and gtdb_taxon_present else (False, 'GTDB taxon or resolved local tree context is unavailable')
        )
        pangenome_eligible = gtdb_taxon_present and member_taxa_agree and tree_context_consistent and tree_clade_agrees
        group_lineage_key = _gtdb_lineage_key(_pangenome_lineage(members[0]), rank)
        taxon_key = (rank, normalise_taxon_name(partner_taxon), group_lineage_key) if pangenome_eligible else ('', '', '')
        # Count all current committed candidates at this GTDB rank before
        # treating the remaining committed IDs as reference genomes.
        current_committed_for_taxon = set()
        if pangenome_eligible:
            for assessment_row in assessment_rows:
                candidate_rank, _candidate_rank_name, candidate_taxon = _partner_rank_taxon(assessment_row)
                if (
                    _committed(assessment_row) and candidate_rank == rank
                    and normalise_taxon_name(candidate_taxon) == normalise_taxon_name(partner_taxon)
                    and _gtdb_lineage_key(_pangenome_lineage(assessment_row), rank) == group_lineage_key
                ): 
                    current_committed_for_taxon.add(str(assessment_row.get('id')))
        committed_ids = (
            sorted(set(committed_by_rank.get(taxon_key, [])) | current_committed_for_taxon)
            if pangenome_eligible else sorted(str(row.get('id')) for row in members if _committed(row))
        )
        # Counts are drawn from the rank-matched committed index. This avoids
        # reusing species-only novelty metrics for genus/family groups.
        committed_count = len(committed_ids) if pangenome_eligible else len([row for row in members if _committed(row)])
        baseline_count = max(0, len(set(committed_ids) - current_committed_for_taxon)) if pangenome_eligible else 0
        # Keep source-specific baseline counts separate from the generic
        # rank-matched total; the rank index combines all registered baselines.
        hungate_count = 0
        secondary_count = 0
        selected_count = len([row for row in members if _selected(row)]) if pangenome_eligible else 0
        pending_count = selected_count
        gap = max(0, target - committed_count) if pangenome_eligible else target
        reported_target = target if pangenome_eligible else 'NA'
        reported_gap = gap if pangenome_eligible else 'NA'
        for member in members:
            member['pangenome_eligible'] = 'Yes' if pangenome_eligible else 'No'
            member['pangenome_target'] = reported_target
            member['pangenome_gap'] = reported_gap
            member['pangenome_rank'] = rank_name.upper() if rank_name else 'UNASSIGNED'
            member['pangenome_committed_count'] = committed_count
            member['pangenome_committed_ids'] = ';'.join(committed_ids) or 'None'
            member['pangenome_baseline_count'] = baseline_count
            member['pangenome_taxon'] = partner_taxon or 'NA'
            classifier_rank_value = str(parse_taxon_string(str(member.get('taxonomy') or '')).get(rank) or '') if rank else ''
            member['classifier_gtdb_disagreement'] = (
                'UNRESOLVED' if not classifier_rank_value else
                'Yes' if normalise_taxon_name(classifier_rank_value) != normalise_taxon_name(partner_taxon) else 'No'
            )
            member['classifier_database_concordance_status'] = _classifier_database_concordance_status(member)
            member['pangenome_taxonomy_provenance'] = 'GTDB' if gtdb_taxon_present else 'UNRESOLVED'
            if not gtdb_taxon_present:
                member['pangenome_eligibility_status'] = 'NOT_ELIGIBLE_NO_USABLE_GTDB_TAXON'
            elif not member_taxa_agree:
                member['pangenome_eligibility_status'] = 'NOT_ELIGIBLE_MIXED_GTDB_TAXA'
            elif not tree_context_consistent:
                member['pangenome_eligibility_status'] = 'NOT_ELIGIBLE_TREE_LEAF_OR_CONTEXT_MISSING'
            elif not tree_clade_agrees:
                member['pangenome_eligibility_status'] = 'NOT_ELIGIBLE_TREE_CLADE_CONFLICT_OR_UNRESOLVED'
            else:
                member['pangenome_eligibility_status'] = 'ELIGIBLE_GTDB_TAXON_TREE_CONCORDANT'
            member['pangenome_eligibility_reason'] = (
                f'GTDB {rank_name} assignment and rank-matched local tree context available; {concordance_reason}'
                if pangenome_eligible else
                f'{concordance_reason}; {tree_clade_reason}'
            )
        novel_clade = _novel_looking(members)
        if pangenome_eligible and rank == 's' and baseline_count > 0:
            group_type = 'BASELINE_PANGENOME_EXTENSION'
            basis = 'GTDB species group with rank-matched cultured-baseline context'
            set_prefix = 'BMEXT'
        else:
            group_type = 'CANDIDATE_ONLY_GROUP'
            set_prefix = 'BMSET'
        is_baseline_extension = group_type != 'CANDIDATE_ONLY_GROUP'
        set_id = _stable_set_id(group_key, set_prefix)

        def same_species_hits(hit_key: str, tax_key: str) -> set[str]:
            if not rank or not partner_taxon:
                return set()
            target_value = _rank_value({'gtdb_taxonomy': members[0].get('gtdb_taxonomy')}, rank)
            return {
                str(row.get(hit_key)) for row in members
                if str(row.get(hit_key) or '') not in {'', 'NA', 'None'}
                and _rank_value({'taxonomy': row.get(tax_key)}, rank)
                and normalise_taxon_name(_rank_value({'taxonomy': row.get(tax_key)}, rank)) == normalise_taxon_name(target_value)
            }

        hungate_anchor_ids = sorted(same_species_hits('hungate_nearest_hit', 'hungate_nearest_hit_taxonomy'))
        secondary_anchor_ids = sorted(same_species_hits('secondary_baseline_nearest_hit', 'secondary_baseline_nearest_hit_taxonomy'))
        cultured_anchor_ids = sorted(same_species_hits('nearest_hit', 'nearest_hit_taxonomy'))
        baseline_anchor_ids = sorted(set(hungate_anchor_ids + secondary_anchor_ids + cultured_anchor_ids))

        for row in members:
            row.update({
                'sequencing_set_id': set_id,
                'sequencing_set_role': 'ALTERNATE',
                'sequencing_set_rank': 'NA',
                'sequencing_set_reason': 'not selected into the proposed working set',
                'selection_group_basis': basis,
                'selection_group_type': group_type,
                'selection_group_taxon': taxon,
                'novel_looking_clade': 'Yes' if novel_clade else 'No',
                'baseline_redundancy_status': 'ELIGIBLE',
                'baseline_redundancy_identity_threshold': f'{redundancy_identity:.2f}',
                'baseline_redundancy_min_query_coverage': f'{redundancy_coverage:.2f}',
                'baseline_extension_status': 'NOT_APPLICABLE_CANDIDATE_GROUP',
                'baseline_extension_min_identity': f'{extension_identity:.2f}',
                'baseline_extension_min_query_coverage': f'{extension_coverage:.2f}',
                'selection_diversity_distance': 'NA',
            })

        eligible = []
        redundant_count = 0
        for row in members:
            if _available(row):
                row['sequencing_set_role'] = 'SEQUENCED'
                row['sequencing_set_reason'] = 'genome already sequenced and available'
                row['baseline_redundancy_status'] = 'NOT_APPLICABLE_COMMITTED'
                row['baseline_extension_status'] = 'NOT_APPLICABLE_COMMITTED'
                continue
            if _selected(row):
                row['sequencing_set_role'] = 'COMMITTED'
                row['sequencing_set_reason'] = 'already selected for genome sequencing; genome not yet available'
                row['baseline_redundancy_status'] = 'NOT_APPLICABLE_COMMITTED'
                row['baseline_extension_status'] = 'NOT_APPLICABLE_COMMITTED'
                continue

            evidence_quality = _evidence_quality(row)
            extension_eligible = False
            extension_status = row.get('baseline_extension_status', 'NOT_APPLICABLE_CANDIDATE_GROUP')
            extension_reason = 'not a baseline-pangenome extension group'
            if is_baseline_extension:
                classifier_species = _species(row.get('taxonomy'))
                if rank == 's' and (
                    not classifier_species
                    or normalise_taxon_name(classifier_species) != normalise_taxon_name(partner_taxon)
                ):
                    # A marker-classifier disagreement is diagnostic only. It
                    # must not move the isolate out of its catalogue-defined
                    # species pangenome group or force boundary review.
                    extension_eligible = True
                    extension_status = 'CATALOGUE_TAXON_RETAINED_CLASSIFIER_DISCORDANT'
                    extension_reason = 'GTDB species grouping retained despite marker-classifier discordance/unresolved species'
                else:
                    extension_eligible, extension_status, extension_reason = baseline_extension_evidence(
                        row, extension_identity, extension_coverage,
                    )
                row['baseline_extension_status'] = extension_status

            redundant, identity, coverage = baseline_redundancy_evidence(
                row, redundancy_identity, redundancy_coverage,
            )
            retained_for_gap = bool(
                redundant and is_baseline_extension and extension_eligible and gap > 0
            )
            if retained_for_gap:
                row['baseline_redundancy_status'] = 'RETAINED_FOR_PANGENOME_GAP'

            if evidence_quality == 'LOW':
                row['sequencing_set_role'] = 'REVIEW_EVIDENCE'
                row['sequencing_set_reason'] = 'marker evidence requires review before selection'
                if redundant and not retained_for_gap:
                    row['baseline_redundancy_status'] = 'PENDING_REVIEW_NEAR_IDENTICAL_BASELINE'
                    if is_baseline_extension:
                        row['baseline_extension_status'] = 'PENDING_REVIEW_NEAR_IDENTICAL_BASELINE'
                    row['sequencing_set_reason'] += (
                        f'; provisional near-identical cultured-baseline match ({identity:.2f}% identity, '
                        f'{coverage:.2f}% query coverage) must not drive exclusion before review'
                    )
                elif retained_for_gap:
                    row['baseline_redundancy_status'] = 'RETAINED_FOR_PANGENOME_GAP_PENDING_REVIEW'
                    row['sequencing_set_reason'] += '; GTDB-rank pangenome gap remains'
                else:
                    row['baseline_redundancy_status'] = 'NOT_EVALUATED_EVIDENCE'
                    row['baseline_extension_status'] = 'NOT_EVALUATED_EVIDENCE'
            elif redundant and not retained_for_gap:
                redundant_count += 1
                row['sequencing_set_role'] = 'BASELINE_REDUNDANT'
                row['baseline_redundancy_status'] = 'EXCLUDED_NEAR_IDENTICAL_BASELINE'
                if is_baseline_extension:
                    row['baseline_extension_status'] = 'EXCLUDED_NEAR_IDENTICAL_BASELINE'
                row['sequencing_set_reason'] = (
                    f'nearest cultured baseline {identity:.2f}% across {coverage:.2f}% of the query; '
                    f'excluded at >={redundancy_identity:.2f}% identity and '
                    f'>={redundancy_coverage:.2f}% query coverage; MWL evidence does not override redundancy'
                )
            elif is_baseline_extension and not extension_eligible:
                row['sequencing_set_role'] = 'PANGENOME_BOUNDARY_REVIEW'
                row['baseline_redundancy_status'] = 'NOT_EVALUATED_PANGENOME_BOUNDARY'
                row['sequencing_set_reason'] = (
                    f'not admitted to baseline-pangenome extension: {extension_reason}; '
                    'review as a possible separate candidate lineage'
                )
            else:
                eligible.append(row)

        ranked = _farthest_first(eligible, committed_ids, distances)
        recommended = ranked[:min(set_size, len(ranked))]
        for rank, (row, spread) in enumerate(recommended, start=1):
            role = 'DIVERSITY_CANDIDATE' if gap == 0 else ('PRIMARY' if rank <= gap else 'BACKUP')
            row['sequencing_set_role'] = role
            row['sequencing_set_rank'] = str(rank)
            row['selection_diversity_distance'] = f'{spread:.6f}' if spread is not None else 'NA'
            if role == 'PRIMARY':
                opening = f'fills genome {committed_count + rank} of target {target}'
            elif role == 'DIVERSITY_CANDIDATE':
                opening = f'ranked diversity expansion after the {target}-genome count target was met'
            else:
                opening = 'ranked backup for DNA-extraction failure and additional within-group diversity'
            row['sequencing_set_reason'] = _ranking_reason(row, opening, spread)

        ranked_count = len(recommended)
        unfilled_ranks = max(0, set_size - ranked_count)

        for row in members:
            output_rows.append({
                'SetID': set_id,
                'GroupType': group_type,
                'GroupBasis': basis,
                'AssessmentTaxon': taxon,
                'BaselineAnchorIDs': ';'.join(baseline_anchor_ids) or 'None',
                'HungateAnchorIDs': ';'.join(hungate_anchor_ids) or 'None',
                'SecondaryBaselineAnchorIDs': ';'.join(secondary_anchor_ids) or 'None',
                'NovelLookingClade': 'Yes' if novel_clade else 'No',
                'HungateBaselineGenomes': hungate_count,
                'SecondaryBaselineGenomes': secondary_count,
                'CulturedProjectBaselineGenomes': baseline_count,
                'BaselineGenomes': baseline_count,
                'SequencedPartnerGenomes': selected_count,
                'SelectedPendingGenomes': pending_count,
                'CommittedGenomeCount': committed_count,
                'CommittedGenomeIDs': ';'.join(committed_ids) or 'None',
                'PangenomeTarget': reported_target,
                'PangenomeGap': reported_gap,
                'PangenomeEligible': 'Yes' if pangenome_eligible else 'No',
                'PangenomeEligibilityReason': (
                    f'GTDB {rank_name} assignment and rank-matched tree context available; {concordance_reason}'
                    if pangenome_eligible else f'{concordance_reason}; {tree_clade_reason}'
                ),
                'PangenomeRank': rank_name.upper() if rank_name else 'UNASSIGNED',
                'PangenomeTaxon': partner_taxon or 'NA',
                'PartnerTaxon': row.get('partner_taxon', 'NA'),
                'ClassifierGTDBDisagreement': row.get('classifier_gtdb_disagreement', 'NA'),
                'TaxonomyAssignmentStatus': _assignment_status(row),
                'ClassifierDatabaseConcordanceStatus': row.get('classifier_database_concordance_status', 'NOT_ASSESSED'),
                'DiversityPanelSize': set_size,
                'EligibleCandidateCount': len(eligible),
                'RankedCandidateCount': ranked_count,
                'BaselineRedundantCount': redundant_count,
                'UnfilledPanelRanks': unfilled_ranks,
                'CandidateID': row.get('id', 'NA'),
                'PartnerID': row.get('partner_id', 'NA'),
                'SetRole': row.get('sequencing_set_role', 'NA'),
                'SetRank': row.get('sequencing_set_rank', 'NA'),
                'EvidenceQuality': _evidence_quality(row),
                'BaselineEvidenceClass': row.get('baseline_evidence_class', 'NA'),
                'HungateNearestIdentity': row.get('hungate_nearest_identity', 'NA'),
                'HungateNearestQueryCoverage': row.get('hungate_nearest_query_coverage', 'NA'),
                'SecondaryBaselineNearestIdentity': row.get('secondary_baseline_nearest_identity', 'NA'),
                'SecondaryBaselineNearestQueryCoverage': row.get('secondary_baseline_nearest_query_coverage', 'NA'),
                'CulturedProjectNearestIdentity': row.get('cultured_project_nearest_identity', row.get('nearest_identity', 'NA')),
                'CulturedProjectNearestQueryCoverage': row.get('cultured_project_nearest_query_coverage', row.get('nearest_query_coverage', 'NA')),
                'BaselineNearestIdentity': row.get('nearest_identity', 'NA'),
                'BaselineNearestQueryCoverage': row.get('nearest_query_coverage', 'NA'),
                'BaselineRedundancyStatus': row.get('baseline_redundancy_status', 'NA'),
                'BaselineRedundancyIdentityThreshold': f'{redundancy_identity:.2f}',
                'BaselineRedundancyMinQueryCoverage': f'{redundancy_coverage:.2f}',
                'BaselineExtensionStatus': row.get('baseline_extension_status', 'NA'),
                'BaselineExtensionMinIdentity': f'{extension_identity:.2f}',
                'BaselineExtensionMinQueryCoverage': f'{extension_coverage:.2f}',
                'ProjectNearestIdentity': row.get('project_nearest_identity', 'NA'),
                'GTDBReferenceIdentity': row.get('reference_nearest_identity', 'NA'),
                'MWLMatchedRank': row.get('mwl_matched_rank', 'NA'),
                'MWLScore': row.get('mwl_score', 'NA'),
                'PhyloIsolation': row.get('phylo_isolation', 'NA'),
                'SelectionDiversityDistance': row.get('selection_diversity_distance', 'NA'),
                'LocalTreeFigure': row.get('local_neighbourhood_figure', 'NA'),
                'SelectionReason': row.get('sequencing_set_reason', 'NA'),
            })

    role_order = {'PRIMARY': 0, 'BACKUP': 1, 'DIVERSITY_CANDIDATE': 2,
                  'COMMITTED': 3, 'SEQUENCED': 4, 'ALTERNATE': 5,
                  'REVIEW_EVIDENCE': 6, 'PANGENOME_BOUNDARY_REVIEW': 7,
                  'BASELINE_REDUNDANT': 8, 'TARGET_MET': 9}
    output_rows.sort(key=lambda row: (
        row['SetID'], role_order.get(str(row['SetRole']), 99),
        _int(row['SetRank'], 999), str(row['CandidateID']),
    ))
    fields = list(output_rows[0]) if output_rows else [
        'SetID', 'GroupType', 'GroupBasis', 'AssessmentTaxon', 'BaselineAnchorIDs',
        'HungateAnchorIDs', 'SecondaryBaselineAnchorIDs',
        'NovelLookingClade', 'HungateBaselineGenomes', 'SecondaryBaselineGenomes',
        'CulturedProjectBaselineGenomes', 'BaselineGenomes',
        'SequencedPartnerGenomes', 'SelectedPendingGenomes', 'CommittedGenomeCount', 'CommittedGenomeIDs',
        'PangenomeTarget', 'PangenomeGap', 'PangenomeEligible', 'PangenomeEligibilityReason',
        'PangenomeRank', 'PangenomeTaxon', 'PartnerTaxon', 'ClassifierGTDBDisagreement',
        'TaxonomyAssignmentStatus', 'ClassifierDatabaseConcordanceStatus',
        'DiversityPanelSize', 'EligibleCandidateCount',
        'RankedCandidateCount', 'BaselineRedundantCount', 'UnfilledPanelRanks',
        'CandidateID', 'PartnerID', 'SetRole', 'SetRank', 'EvidenceQuality',
        'BaselineEvidenceClass', 'HungateNearestIdentity', 'HungateNearestQueryCoverage',
        'SecondaryBaselineNearestIdentity', 'SecondaryBaselineNearestQueryCoverage',
        'CulturedProjectNearestIdentity', 'CulturedProjectNearestQueryCoverage',
        'BaselineNearestIdentity', 'BaselineNearestQueryCoverage', 'BaselineRedundancyStatus',
        'BaselineRedundancyIdentityThreshold', 'BaselineRedundancyMinQueryCoverage',
        'BaselineExtensionStatus', 'BaselineExtensionMinIdentity',
        'BaselineExtensionMinQueryCoverage',
        'ProjectNearestIdentity', 'GTDBReferenceIdentity', 'MWLMatchedRank', 'MWLScore', 'PhyloIsolation',
        'SelectionDiversityDistance', 'LocalTreeFigure', 'SelectionReason',
    ]
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(output_rows)
    return str(path)
