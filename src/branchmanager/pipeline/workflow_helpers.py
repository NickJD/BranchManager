from __future__ import annotations

from branchmanager.pipeline.evidence_policy import evidence_quality as _shared_evidence_quality
import gzip
import json
import os
import csv
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

from branchmanager.taxonomy_io import iter_taxonomy_assignment_rows
from branchmanager.taxonomy import (
    normalise_taxon_name,
    parse_reference_header_taxonomy,
    parse_taxon_string,
    domain_designation,
    infer_taxon_rank,
    taxonomy_matches_kingdom,
)
from branchmanager.utils.fasta import read_fasta


ClassificationRow = Dict[str, object]


def _optional_numeric(parts: Sequence[str], index: int, *, integer: bool = False):
    if index >= len(parts) or parts[index] in ('', 'NA', 'None'):
        return None
    try:
        value = float(parts[index])
        return int(value) if integer else value
    except (TypeError, ValueError):
        return None


def _assignment_source_is_fasta(assignments_path: str, source_fasta_path: Optional[str] = None) -> bool:
    try:
        if source_fasta_path and os.path.samefile(assignments_path, source_fasta_path):
            return True
    except Exception:
        pass  # samefile may fail on different mounts or symlinks

    try:
        if str(assignments_path).endswith('.gz'):
            handle_ctx = gzip.open(assignments_path, 'rt')
        else:
            handle_ctx = open(assignments_path, 'rt')
        with handle_ctx as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped:
                    continue
                return stripped.startswith('>')
    except Exception:
        logger.debug('Could not read file %s to check if it is a FASTA', assignments_path)
        return False
    return False


def build_orig_to_short_map(alias_entries, db) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    if not alias_entries:
        return mapping
    for short, orig in alias_entries:
        mapping[orig] = short
        mapping[short] = short
        try:
            cid = db._canonical_from_header(orig)
            if cid:
                mapping[cid] = short
        except Exception:
            logger.debug('Could not resolve canonical ID for alias %r', orig)
    return mapping


def resolve_short_id(qid: str, orig_to_short: Dict[str, str], db) -> Optional[str]:
    if qid in orig_to_short:
        return orig_to_short[qid]
    try:
        cid = db._canonical_from_header(qid)
    except Exception:
        logger.debug('Could not resolve canonical ID for %r', qid)
        cid = None
    if isinstance(cid, str) and cid in orig_to_short:
        return orig_to_short[cid]
    try:
        lk = qid.lower()
        for key, value in orig_to_short.items():
            if isinstance(key, str) and key.lower() == lk:
                return value
    except Exception:
        logger.debug('Case-insensitive ID lookup failed for %r', qid)
    return None


def iter_classification_rows(classification_tsv: str) -> Iterator[ClassificationRow]:
    with open(classification_tsv) as handle:
        header = next(handle, None)
        header_fields = header.rstrip('\n').split('\t') if header else []
        source_confidence_index = next(
            (i for i, name in enumerate(header_fields)
             if 'sourceclassifierconfidence' in name.lower().replace('_', '')),
            None,
        )
        for line in handle:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 1 or not parts[0]:
                continue
            try:
                identity = float(parts[2]) if len(parts) > 2 and parts[2] not in ('', 'NA') else None
            except Exception:
                identity = None
            try:
                confidence = (
                    float(parts[source_confidence_index])
                    if source_confidence_index is not None
                    and len(parts) > source_confidence_index
                    and parts[source_confidence_index] not in ('', 'NA')
                    else None
                )
            except Exception:
                confidence = None
            yield {
                'qid': parts[0],
                'best': parts[1] if len(parts) > 1 and parts[1] != '' else None,
                'identity': identity,
                'tax': parts[3] if len(parts) > 3 and parts[3] != '' else None,
                'confidence': confidence,
                'query_coverage': _optional_numeric(parts, 5),
                'target_coverage': _optional_numeric(parts, 6),
                'alignment_length': _optional_numeric(parts, 7, integer=True),
                'query_length': _optional_numeric(parts, 8, integer=True),
                'target_length': _optional_numeric(parts, 9, integer=True),
                'mismatches': _optional_numeric(parts, 10, integer=True),
                'gaps': _optional_numeric(parts, 11, integer=True),
                'assignment_status': parts[12] if len(parts) > 12 else None,
                'taxonomy_support_pct': _optional_numeric(parts, 13),
                'taxonomy_lead_pct': _optional_numeric(parts, 14),
                'assigned_rank': parts[15] if len(parts) > 15 else None,
                'evidence_tier': parts[16] if len(parts) > 16 else None,
                'taxonomy_candidate_count': _optional_numeric(parts, 17, integer=True),
                'alignment_evidence_status': parts[18] if len(parts) > 18 else None,
                'assignment_evidence_level': parts[19] if len(parts) > 19 else None,
                'evidence_interpretation': parts[20] if len(parts) > 20 else None,
            }


def iter_assignment_rows(assignments_path: str, source_fasta_path: Optional[str] = None) -> Iterator[dict[str, object]]:
    """Yield normalised taxonomy-assignment rows from table files or FASTA headers.

    If `assignments_path` is the same file as `source_fasta_path`, taxonomy is
    parsed from GTDB-style FASTA headers. Otherwise TSV/CSV and .gz variants are
    parsed as ID-to-taxonomy assignment tables.
    """
    use_fasta_headers = _assignment_source_is_fasta(assignments_path, source_fasta_path=source_fasta_path)

    if use_fasta_headers:
        for header, _ in read_fasta(assignments_path):
            _, tax = parse_reference_header_taxonomy(header)
            yield {
                'qid': header,
                'tax': tax,
                'confidence': None,
            }
        return

    for row in iter_taxonomy_assignment_rows(assignments_path):
        yield {
            'qid': row.get('id'),
            'tax': row.get('taxonomy'),
            'confidence': row.get('confidence'),
        }


def load_taxonomy_entries_from_assignments(taxa_file: str, orig_to_short: Dict[str, str], db, dataset: str, source_fasta_path: Optional[str] = None):
    entries = []
    for row in iter_assignment_rows(taxa_file, source_fasta_path=source_fasta_path):
        mapped = resolve_short_id(str(row['qid']), orig_to_short, db)
        if not mapped:
            continue
        entries.append((mapped, row.get('tax'), row.get('confidence'), dataset))
    return entries


def load_classification_results_for_dataset(classification_tsv: str, orig_to_short: Dict[str, str], db, dataset: str):
    tax_entries = []
    dist_entries = []
    for row in iter_classification_rows(classification_tsv):
        mapped = resolve_short_id(str(row['qid']), orig_to_short, db)
        if not mapped:
            continue
        tax_entries.append((mapped, row['tax'], row['confidence'], dataset))
        dist_entries.append((mapped, dataset, row['best'], row['identity']))
    return tax_entries, dist_entries


def classification_ids_matching_kingdom(classification_tsv: str, kingdom: str):
    matched = set()
    if not classification_tsv or not kingdom:
        return matched
    for row in iter_classification_rows(classification_tsv):
        raw_tax = row.get('tax')
        tax = str(raw_tax) if raw_tax is not None else None
        if taxonomy_matches_kingdom(tax, kingdom):
            matched.add(str(row['qid']))
    return matched


def prune_dataset_by_kingdom(db, dataset: str, kingdom: str | None):
    if not kingdom:
        return 0
    to_delete = []
    with db.connect() as conn:
        cur = conn.cursor()
        cur.execute(
            "SELECT s.id, t.taxonomy FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id WHERE s.dataset = ?",
            (dataset,),
        )
        for iid, tax in cur.fetchall():
            if tax is None:
                continue
            tax_text = str(tax).strip()
            if not tax_text or tax_text.lower() in ('na', 'none'):
                continue
            if not taxonomy_matches_kingdom(tax_text, kingdom):
                to_delete.append(iid)
        for iid in to_delete:
            cur.execute('DELETE FROM distances WHERE id = ?', (iid,))
            cur.execute('DELETE FROM taxonomy WHERE id = ?', (iid,))
            cur.execute('DELETE FROM colours WHERE id = ?', (iid,))
            cur.execute('DELETE FROM sequences WHERE id = ?', (iid,))
        conn.commit()
    return len(to_delete)


def read_combined_taxonomy_ids(path: str | Path):
    ids = []
    p = Path(path)
    if not p.exists():
        return ids
    with open(p) as handle:
        next(handle, None)
        for line in handle:
            iid = line.strip().split('\t')[0]
            if iid:
                ids.append(iid)
    return ids


def collect_db_taxonomy_rows(db, ids: Optional[Sequence[str]] = None):
    with db.connect() as conn:
        cur = conn.cursor()
        if ids:
            placeholders = ','.join('?' for _ in ids)
            cur.execute(
                f"SELECT s.id, t.taxonomy, t.confidence FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id WHERE s.id IN ({placeholders})",
                tuple(ids),
            )
        else:
            cur.execute("SELECT s.id, t.taxonomy, t.confidence FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id")
        return cur.fetchall()


def merge_combined_taxonomy_rows(base_rows, classification_tsv: str, orig_to_short: Dict[str, str], db):
    merged = {}
    order: List[str] = []
    for rid, tax, conf in base_rows:
        merged[rid] = (tax if tax is not None else 'NA', conf if conf is not None else 'NA')
        order.append(rid)
    if classification_tsv:
        for row in iter_classification_rows(classification_tsv):
            mapped = resolve_short_id(str(row['qid']), orig_to_short, db)
            if not mapped:
                continue
            if mapped not in merged:
                order.append(mapped)
            merged[mapped] = (
                row['tax'] if row['tax'] is not None else 'NA',
                row['confidence'] if row['confidence'] is not None else 'NA',
            )
    # The combined tree/iTOL taxonomy is partner-catalogue-led for partner
    # isolates. Classifier taxonomy remains available in the assessment tables
    # and database evidence; it must never overwrite the catalogue identity.
    try:
        with db.connect() as conn:
            metadata_rows = conn.execute(
                'SELECT id, project_metadata_json FROM sequencing_metadata'
            ).fetchall()
        for sid, metadata_json in metadata_rows:
            metadata = json.loads(metadata_json) if metadata_json else {}
            raw = str(metadata.get('partner_taxon_raw') or metadata.get('expected_taxon') or '').strip()
            if not raw or str(sid) not in merged:
                continue
            parsed = parse_taxon_string(raw)
            if not parsed and domain_designation(raw):
                parsed = {'d': domain_designation(raw)}
            if not parsed:
                clean = ' '.join(raw.split())
                rank = infer_taxon_rank(clean, metadata.get('partner_taxon_rank'))
                if rank == 's' and len(clean.split()) >= 2:
                    genus, species = clean.split()[:2]
                    parsed = {'g': genus, 's': f'{genus} {species}'}
                else:
                    parsed = {rank: clean}
            lineage = ';'.join(f'{rank}__{value}' for rank, value in parsed.items())
            if lineage:
                merged[str(sid)] = (lineage, 'PARTNER_CATALOGUE_GOLD_STANDARD')
    except Exception:
        pass
    return [(iid, merged[iid][0], merged[iid][1]) for iid in order]


def write_combined_taxonomy_tsv(path: str | Path, rows: Iterable[Tuple[str, object, object]]):
    p = Path(path)
    with open(p, 'w') as outf:
        outf.write('ID\tTaxon\tConfidence\n')
        for iid, tax, conf in rows:
            outf.write(f"{iid}\t{tax if tax is not None else 'NA'}\t{conf if conf is not None else 'NA'}\n")
    return str(p)


def iter_novelty_rows(novelty_tsv: str):
    with open(novelty_tsv) as handle:
        next(handle, None)
        for line in handle:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            try:
                nearest_identity = float(parts[1]) if parts[1] not in ('', 'NA') else None
            except Exception:
                nearest_identity = None
            yield {
                'qid': parts[0],
                'nearest_identity': nearest_identity,
                'nearest_hit': parts[2] if len(parts) > 2 else None,
                'novel': parts[3] if len(parts) > 3 else None,
            }


def build_placement_warning_rows(classification_tsv: str, novelty_tsv: Optional[str], orig_to_short: Dict[str, str], db):
    novelty_by_qid = {}
    if novelty_tsv:
        try:
            for row in iter_novelty_rows(novelty_tsv):
                novelty_by_qid[str(row['qid'])] = row
        except Exception:
            novelty_by_qid = {}

    warning_rows = []
    for row in iter_classification_rows(classification_tsv):
        qid = str(row['qid'])
        mapped = resolve_short_id(qid, orig_to_short, db) or qid
        flags = []
        best = row.get('best')
        identity = row.get('identity')
        query_coverage = row.get('query_coverage')
        target_coverage = row.get('target_coverage')
        tax = row.get('tax')
        nov = novelty_by_qid.get(qid, {})
        nearest_identity = nov.get('nearest_identity')

        if not best or best == 'NA' or identity in (None, 0, 0.0):
            flags.append('NO_REFERENCE_HIT')
        if isinstance(identity, (int, float)) and float(identity) < 95.0:
            flags.append('LOW_CLASSIFICATION_IDENTITY')
        if isinstance(query_coverage, (int, float)) and float(query_coverage) < 90.0:
            flags.append('LOW_QUERY_COVERAGE')
        if isinstance(target_coverage, (int, float)) and float(target_coverage) < 90.0:
            flags.append('LOW_REFERENCE_COVERAGE')
        if isinstance(nearest_identity, (int, float)) and float(nearest_identity) < 97.0:
            flags.append('LOW_NEAREST_IDENTITY')
            if tax not in (None, '', 'NA'):
                flags.append('NOVEL_BUT_ASSIGNED')
        if nov and nearest_identity in (None, 0, 0.0) and tax not in (None, '', 'NA'):
            flags.append('ASSIGNED_WITHOUT_PLACEMENT_SUPPORT')

        if flags:
            warning_rows.append({
                'id': mapped,
                'query_id': qid,
                'taxonomy': tax if tax is not None else 'NA',
                'classification_identity': identity if identity is not None else 'NA',
                'classification_confidence': row.get('confidence') if row.get('confidence') is not None else 'NA',
                'nearest_identity': nearest_identity if nearest_identity is not None else 'NA',
                'nearest_hit': nov.get('nearest_hit', 'NA') if nov else 'NA',
                'flags': ';'.join(flags),
            })
    return warning_rows


def write_placement_warning_tsv(path: str | Path, warning_rows):
    p = Path(path)
    with open(p, 'w') as fh:
        fh.write('ID\tQueryID\tTaxonomy\tClassificationIdentity\tSourceClassifierConfidence\tNearestIdentity\tNearestHit\tFlags\n')
        for row in warning_rows:
            fh.write(
                f"{row['id']}\t{row['query_id']}\t{row['taxonomy']}\t{row['classification_identity']}\t"
                f"{row['classification_confidence']}\t{row['nearest_identity']}\t{row['nearest_hit']}\t{row['flags']}\n"
            )
    return str(p)


def _taxonomy_agreement(query_taxonomy, comparison_taxonomy) -> tuple[str, str]:
    """Return deepest shared rank and first explicit conflict for two lineages."""
    ranks = [('d', 'domain'), ('p', 'phylum'), ('c', 'class'), ('o', 'order'),
             ('f', 'family'), ('g', 'genus'), ('s', 'species')]
    query = parse_taxon_string(str(query_taxonomy or ''))
    comparison = parse_taxon_string(str(comparison_taxonomy or ''))
    if not query or not comparison:
        return 'NA', 'NA'
    deepest = 'none'
    for key, name in ranks:
        left = normalise_taxon_name(query.get(key, ''))
        right = normalise_taxon_name(comparison.get(key, ''))
        if not left or not right:
            continue
        if left != right:
            return deepest, f'{name}:{query.get(key)}!={comparison.get(key)}'
        deepest = name
    return deepest, 'none'


def build_sequence_assessment_rows(
    run_ids,
    class_out: str,
    novelty_metrics_tsv: str,
    placement_warning_rows,
    orig_to_short: Dict[str, str],
    db,
    member_to_rep: Optional[Dict[str, str]] = None,
    rep_to_members: Optional[Dict[str, List[str]]] = None,
    tree_ids: Optional[set] = None,
    alt_taxonomies: Optional[Dict[str, Dict[str, tuple]]] = None,
    alt_ref_dbs: Optional[List[str]] = None,
    alt_classification_evidence: Optional[Dict[str, Dict[str, dict]]] = None,
):
    """Build per-sequence assessment rows.

    Optional clustering info (``member_to_rep``, ``rep_to_members``) is
    included when ``--collapse`` was used.  Tree membership (``tree_ids``)
    records whether a sequence entered the phylogenetic tree or was filtered
    out because a cluster representative was used in its place.

    ``alt_taxonomies`` is an optional dict ``{seq_id: {ref_db: (tax, conf, hit, identity)}}``
    produced when multiple reference databases are used for classification.
    ``alt_ref_dbs`` is the ordered list of alternative ref-db names to include.

    Cluster columns added to each row
    -----------------------------------
    InTree              : 'Yes' / 'No' / 'Unknown'
    ClusterRepresentative : seq-id of the representative for this cluster,
                            or 'self' if this sequence IS in the tree as its own entry,
                            or 'duplicate' if omitted from the final tree as an exact duplicate,
                            or 'N/A' if no clustering was done.
    ClusterSize         : number of sequences in this cluster (1 if singleton).
    ClusteredMembers    : semi-colon-separated list of member IDs (empty for singletons).
    """
    run_id_set = {str(iid) for iid in run_ids}
    class_by_id = {}
    for row in iter_classification_rows(class_out):
        qid = str(row['qid'])
        mapped = resolve_short_id(qid, orig_to_short, db)
        if not mapped and qid in run_id_set:
            mapped = qid
        if mapped:
            class_by_id[mapped] = row

    # Reload stored evidence when a submitted ID was already present in the
    # project database or was not represented in the current classifier output.
    if run_id_set:
        try:
            placeholders = ','.join('?' for _ in run_id_set)
            with db.connect() as conn:
                taxonomy_rows = conn.execute(
                    f"SELECT id, taxonomy, confidence FROM taxonomy "
                    f"WHERE id IN ({placeholders}) ORDER BY rowid",
                    tuple(sorted(run_id_set)),
                ).fetchall()
                distance_rows = conn.execute(
                    f"SELECT id, nearest, identity FROM distances "
                    f"WHERE id IN ({placeholders}) ORDER BY rowid",
                    tuple(sorted(run_id_set)),
                ).fetchall()
            db_classification = {}
            for sid, taxonomy, confidence in taxonomy_rows:
                entry = db_classification.setdefault(str(sid), {'qid': str(sid)})
                if taxonomy not in (None, ''):
                    entry['tax'] = taxonomy
            for sid, nearest, identity in distance_rows:
                entry = db_classification.setdefault(str(sid), {'qid': str(sid)})
                if nearest not in (None, ''):
                    entry['best'] = nearest
                if identity is not None:
                    entry['identity'] = identity
            for sid, entry in db_classification.items():
                class_by_id.setdefault(sid, entry)
            evidence_rows = db.get_classification_evidence(run_id_set)
            evidence_by_id = {}
            for (sid, _ref_db), evidence in evidence_rows.items():
                evidence_by_id.setdefault(sid, evidence)
            for sid, evidence in evidence_by_id.items():
                class_by_id.setdefault(sid, {'qid': sid}).update({
                    key: evidence.get(key) for key in (
                        'query_coverage', 'target_coverage', 'alignment_length', 'query_length',
                        'target_length', 'mismatches', 'gaps',
                    )
                })
                for key in (
                    'source_classifier_confidence', 'assignment_status',
                    'taxonomy_support_pct', 'taxonomy_candidate_count', 'evidence_tier',
                    'alignment_evidence_status', 'assignment_evidence_level',
                    'evidence_interpretation',
                ):
                    target_key = 'confidence' if key == 'source_classifier_confidence' else key
                    if class_by_id[sid].get(target_key) in (None, ''):
                        class_by_id[sid][target_key] = evidence.get(key)
        except Exception:
            logger.debug('Could not read evidence TSV %s', evidence_tsv)

    novelty_by_id = {}
    try:
        with open(novelty_metrics_tsv) as fh:
            reader = csv.DictReader(fh, delimiter='	')
            for row in reader:
                qid = row.get('ID') or row.get('id')
                if not qid:
                    continue

                def val(*names, default='NA'):
                    for name in names:
                        got = row.get(name)
                        if got not in (None, ''):
                            return got
                    return default

                def source_available(*names) -> bool:
                    source = val(*names, default='')
                    return source not in ('', 'none', 'NA')

                cultured_available = source_available(
                    'CulturedProjectDensitySource', 'CulturedRumenDensitySource', 'BaselineDensitySource',
                )
                hungate_available = source_available('HungateDensitySource')
                secondary_available = source_available('SecondaryBaselineDensitySource')

                def pool_val(prefix: str, field: str, *, available: bool, fallback_prefix=None, default='NA'):
                    if not available:
                        return default
                    names = [prefix + field]
                    if fallback_prefix:
                        prefixes = (
                            [fallback_prefix] if isinstance(fallback_prefix, str)
                            else list(fallback_prefix)
                        )
                        names.extend(item + field for item in prefixes)
                    return val(*names, default=default)

                novelty_by_id[qid] = {
                    'baseline_evidence_class': val('BaselineEvidenceClass'),
                    'hungate_nearest_identity': pool_val('Hungate', 'NearestIdentity', available=hungate_available),
                    'hungate_nearest_hit': pool_val('Hungate', 'NearestHit', available=hungate_available),
                    'hungate_nearest_query_coverage': pool_val('Hungate', 'NearestQueryCoverage', available=hungate_available),
                    'hungate_nearest_alignment_length': pool_val('Hungate', 'NearestAlignmentLength', available=hungate_available),
                    'hungate_novelty_score': pool_val('Hungate', 'NoveltyScore', available=hungate_available),
                    'hungate_density_source': pool_val('Hungate', 'DensitySource', available=hungate_available),
                    'hungate_density_count_method': pool_val('Hungate', 'DensityCountMethod', available=hungate_available),
                    'secondary_baseline_nearest_identity': pool_val('SecondaryBaseline', 'NearestIdentity', available=secondary_available),
                    'secondary_baseline_nearest_hit': pool_val('SecondaryBaseline', 'NearestHit', available=secondary_available),
                    'secondary_baseline_nearest_query_coverage': pool_val('SecondaryBaseline', 'NearestQueryCoverage', available=secondary_available),
                    'secondary_baseline_nearest_alignment_length': pool_val('SecondaryBaseline', 'NearestAlignmentLength', available=secondary_available),
                    'secondary_baseline_novelty_score': pool_val('SecondaryBaseline', 'NoveltyScore', available=secondary_available),
                    'secondary_baseline_density_source': pool_val('SecondaryBaseline', 'DensitySource', available=secondary_available),
                    'secondary_baseline_density_count_method': pool_val('SecondaryBaseline', 'DensityCountMethod', available=secondary_available),
                    'cultured_project_nearest_identity': pool_val(
                        'CulturedProject', 'NearestIdentity', available=cultured_available,
                        fallback_prefix=('CulturedRumen', 'Baseline'),
                    ),
                    'cultured_project_nearest_hit': pool_val(
                        'CulturedProject', 'NearestHit', available=cultured_available,
                        fallback_prefix=('CulturedRumen', 'Baseline'),
                    ),
                    'cultured_project_nearest_query_coverage': pool_val(
                        'CulturedProject', 'NearestQueryCoverage', available=cultured_available,
                        fallback_prefix=('CulturedRumen', 'Baseline'),
                    ),
                    'cultured_project_nearest_alignment_length': pool_val(
                        'CulturedProject', 'NearestAlignmentLength', available=cultured_available,
                        fallback_prefix=('CulturedRumen', 'Baseline'),
                    ),
                    'cultured_project_novelty_score': pool_val(
                        'CulturedProject', 'NoveltyScore', available=cultured_available,
                        fallback_prefix=('CulturedRumen', 'Baseline'),
                    ),
                    'cultured_project_density_source': pool_val(
                        'CulturedProject', 'DensitySource', available=cultured_available,
                        fallback_prefix=('CulturedRumen', 'Baseline'),
                    ),
                    'nearest_identity': pool_val(
                        'CulturedProject', 'NearestIdentity', available=cultured_available,
                        fallback_prefix=('CulturedRumen', 'Baseline'),
                    ),
                    'nearest_hit': pool_val(
                        'CulturedProject', 'NearestHit', available=cultured_available,
                        fallback_prefix=('CulturedRumen', 'Baseline'),
                    ),
                    'nearest_query_coverage': pool_val(
                        'CulturedProject', 'NearestQueryCoverage', available=cultured_available,
                        fallback_prefix=('CulturedRumen', 'Baseline'),
                    ),
                    'nearest_alignment_length': pool_val(
                        'CulturedProject', 'NearestAlignmentLength', available=cultured_available,
                        fallback_prefix=('CulturedRumen', 'Baseline'),
                    ),
                    'novel': pool_val('CulturedProject', 'Novel', available=cultured_available, fallback_prefix=('CulturedRumen', 'Baseline')),
                    'matches_ge_99': pool_val('CulturedProject', 'MatchesGE99', available=cultured_available, fallback_prefix=('CulturedRumen', 'Baseline')),
                    'matches_ge_97': pool_val('CulturedProject', 'MatchesGE97', available=cultured_available, fallback_prefix=('CulturedRumen', 'Baseline')),
                    'matches_ge_95': pool_val('CulturedProject', 'MatchesGE95', available=cultured_available, fallback_prefix=('CulturedRumen', 'Baseline')),
                    'novelty_score': pool_val('CulturedProject', 'NoveltyScore', available=cultured_available, fallback_prefix=('CulturedRumen', 'Baseline')),
                    'crowding': pool_val('CulturedProject', 'Crowding', available=cultured_available, fallback_prefix=('CulturedRumen', 'Baseline')),
                    'sequencing_priority': pool_val('CulturedProject', 'SequencingPriority', available=cultured_available, fallback_prefix=('CulturedRumen', 'Baseline')),
                    'density_source': pool_val('CulturedProject', 'DensitySource', available=cultured_available, fallback_prefix=('CulturedRumen', 'Baseline')),
                    'density_count_method': pool_val('CulturedProject', 'DensityCountMethod', available=cultured_available, fallback_prefix=('CulturedRumen', 'Baseline')),
                    'project_nearest_identity': val('ProjectNearestIdentity'),
                    'project_nearest_hit': val('ProjectNearestHit'),
                    'project_nearest_query_coverage': val('ProjectNearestQueryCoverage'),
                    'project_nearest_alignment_length': val('ProjectNearestAlignmentLength'),
                    'project_novel': val('ProjectNovel'),
                    'project_matches_ge_99': val('ProjectMatchesGE99'),
                    'project_matches_ge_97': val('ProjectMatchesGE97'),
                    'project_matches_ge_95': val('ProjectMatchesGE95'),
                    'project_novelty_score': val('ProjectNoveltyScore'),
                    'project_crowding': val('ProjectCrowding'),
                    'project_sequencing_priority': val('ProjectSequencingPriority'),
                    'project_density_source': val('ProjectDensitySource'),
                    'project_density_count_method': val('ProjectDensityCountMethod'),
                    'reference_nearest_identity': val('ReferenceNearestIdentity'),
                    'reference_nearest_hit': val('ReferenceNearestHit'),
                    'reference_nearest_query_coverage': val('ReferenceNearestQueryCoverage'),
                    'reference_nearest_alignment_length': val('ReferenceNearestAlignmentLength'),
                    'reference_novel': val('ReferenceNovel'),
                    'reference_matches_ge_99': val('ReferenceMatchesGE99'),
                    'reference_matches_ge_97': val('ReferenceMatchesGE97'),
                    'reference_matches_ge_95': val('ReferenceMatchesGE95'),
                    'reference_novelty_score': val('ReferenceNoveltyScore'),
                    'reference_crowding': val('ReferenceCrowding'),
                    'reference_sequencing_priority': val('ReferenceSequencingPriority'),
                    'reference_density_source': val('ReferenceDensitySource'),
                    'reference_density_count_method': val('ReferenceDensityCountMethod'),
                    'partner_id': val('PartnerID'),
                    'selected_for_genome_sequencing': val('AlreadySelectedForGenomeSequencing', 'SelectedForGenomeSequencing'),
                    'already_sequenced': val('GenomeAlreadySequenced'),
                    'nearest_genome_hit': val('NearestGenomeHit'),
                    'nearest_genome_identity': val('NearestGenomeIdentity'),
                    'genome_collection_matches_ge_99': val('GenomeCollectionMatchesGE99'),
                    'genome_collection_matches_ge_97': val('GenomeCollectionMatchesGE97'),
                    'genome_collection_matches_ge_95': val('GenomeCollectionMatchesGE95'),
                    'related_genome_clade_ge_97': val('RelatedGenomeCladeGE97'),
                    'hungate_genome_count_same_species': val('HungateGenomeCountSameAssessmentSpecies', default='0'),
                    'secondary_baseline_genome_count_same_species': val('SecondaryBaselineGenomeCountSameAssessmentSpecies', default='0'),
                    'cultured_project_genome_count_same_species': val(
                        'CulturedProjectGenomeCountSameAssessmentSpecies',
                        default=val('BaselineGenomeCountSameAssessmentSpecies', default='0'),
                    ),
                    'genome_committed_count_same_species': val(
                        'CommittedGenomeCountSameAssessmentSpecies',
                        default=val('AvailableGenomeCountSameAssessmentSpecies',
                                    default=val('GenomeCommittedCountSameAssessmentSpecies')),
                    ),
                    'genome_available_count_same_species': val(
                        'CulturedProjectGenomeCountSameAssessmentSpecies',
                        default=val('BaselineGenomeCountSameAssessmentSpecies',
                                    default=val('GenomeAvailableCountSameAssessmentSpecies')),
                    ),
                    'genome_selected_count_same_species': val(
                        'SequencedPartnerGenomeCountSameAssessmentSpecies',
                        default=val('GenomeSelectedCountSameAssessmentSpecies'),
                    ),
                    'genome_pending_count_same_species': val(
                        'SelectedPendingGenomeCountSameAssessmentSpecies',
                    ),
                    'pangenome_target': val('PangenomeTarget', default='NA'),
                    'pangenome_gap': val('PangenomeGap'),
                    'pangenome_eligibility_status': val('PangenomeEligibilityStatus', default='DEFERRED_TO_TREE_VALIDATED_SELECTION_SETS'),
                    'genome_sequencing_metadata_source': val('GenomeSequencingMetadataSource'),
                    'closest_project_isolate': val('ClosestProjectIsolate', default='None'),
                    'closest_project_isolate_identity': val('ClosestProjectIsolateIdentity', default='0.00'),
                    'project_isolate_only_match': val('ProjectIsolateOnlyMatch', default='False'),
                }
    except FileNotFoundError:
        import logging as _logging
        _logging.getLogger(__name__).warning(
            "[ASSESSMENT] Novelty metrics file not found: %s — novelty columns will be NA",
            novelty_metrics_tsv,
        )

    nearest_meta = {}
    nearest_meta = {}
    nearest_ids_set = set()
    for row in novelty_by_id.values():
        for key in (
            'nearest_hit', 'hungate_nearest_hit', 'secondary_baseline_nearest_hit',
            'cultured_project_nearest_hit', 'project_nearest_hit', 'reference_nearest_hit',
            'closest_project_isolate',
        ):
            hit = row.get(key)
            if hit not in (None, '', 'NA', 'None'):
                nearest_ids_set.add(str(hit))
    nearest_ids = sorted(nearest_ids_set)
    if nearest_ids:
        try:
            with db.connect() as conn:
                cur = conn.cursor()
                placeholders = ','.join('?' for _ in nearest_ids)
                cur.execute(
                    "SELECT s.id, s.dataset, t.taxonomy "
                    "FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id "
                    f"WHERE s.id IN ({placeholders})",
                    tuple(nearest_ids),
                )
                for hit_id, dataset, taxonomy in cur.fetchall():
                    nearest_meta[str(hit_id)] = {
                        'dataset': dataset or 'NA',
                        'taxonomy': taxonomy or 'NA',
                    }
        except Exception:
            nearest_meta = {}

    warnings_by_id = {row['id']: row for row in placement_warning_rows}

    # Partner-catalogue taxonomy is the project-level gold standard. Keep it
    # alongside (never in place of) the sequence-derived classifier lineage.
    partner_metadata_by_id = {}
    try:
        partner_metadata_by_id = db.get_sequencing_metadata_for_ids(run_id_set)
    except Exception:
        logger.warning('Could not load partner-catalogue taxonomy metadata for assessment rows')

    # GTDB genome taxonomy is a separate, authoritative evidence channel for
    # pangenome grouping. Do not infer it from the marker classifier taxonomy.
    gtdb_genome_by_id = {}
    try:
        if run_id_set:
            placeholders = ','.join('?' for _ in run_id_set)
            with db.connect() as conn:
                gtdb_rows = conn.execute(
                    f"SELECT sequence_id, gtdb_taxonomy, genome_qc_pass FROM genome_records "
                    f"WHERE sequence_id IN ({placeholders})",
                    tuple(sorted(run_id_set)),
                ).fetchall()
            for sid, taxonomy, qc_pass in gtdb_rows:
                if taxonomy:
                    gtdb_genome_by_id[str(sid)] = {
                        'taxonomy': str(taxonomy),
                        'qc_pass': qc_pass,
                    }
    except Exception:
        logger.warning('Could not load GTDB genome taxonomy for assessment rows')

    # If the configured primary classifier itself is GTDB, its stored taxonomy
    # is valid GTDB evidence. Never promote SILVA/other classifier assignments.
    try:
        classification_evidence = db.get_classification_evidence(run_id_set) if run_id_set else {}
        taxonomy_by_id = {}
        if run_id_set:
            placeholders = ','.join('?' for _ in run_id_set)
            with db.connect() as conn:
                taxonomy_by_id = {
                    str(sid): str(taxonomy or '')
                    for sid, taxonomy in conn.execute(
                        f"SELECT id, taxonomy FROM taxonomy WHERE id IN ({placeholders})",
                        tuple(sorted(run_id_set)),
                    ).fetchall()
                }
        for (sid, ref_db), evidence in classification_evidence.items():
            ref_label = ''.join(char for char in str(ref_db or '').casefold() if char.isalnum())
            if not (ref_label == 'gtdb' or ref_label.startswith('gtdbtk') or ref_label.startswith('gtdbr')):
                continue
            if str(sid) in gtdb_genome_by_id:
                continue
            taxonomy = taxonomy_by_id.get(str(sid), '')
            if taxonomy:
                gtdb_genome_by_id[str(sid)] = {
                    'taxonomy': taxonomy,
                    'qc_pass': '',
                    'assignment_status': str(evidence.get('assignment_status') or 'UNRESOLVED'),
                }
    except Exception:
        logger.debug('Could not load GTDB classifier evidence for assessment rows')

    # Normalise cluster dicts
    _m2r: Dict[str, str] = member_to_rep or {}
    _r2m: Dict[str, List[str]] = rep_to_members or {}
    _tree: Optional[set] = tree_ids  # may be None (collapse not used)

    assessment_rows = []
    for iid in run_ids:
        c = class_by_id.get(iid, {})
        n = novelty_by_id.get(iid, {})
        w = warnings_by_id.get(iid, {})
        partner_meta = partner_metadata_by_id.get(str(iid), {})
        partner_project_meta = partner_meta.get('project_metadata') or {}
        partner_taxon = str(
            partner_project_meta.get('partner_taxon_raw')
            or partner_project_meta.get('expected_taxon')
            or ''
        ).strip()
        partner_taxon_rank = str(
            partner_project_meta.get('partner_taxon_rank')
            or partner_project_meta.get('taxonomic_rank')
            or ''
        ).strip()

        if _m2r or _r2m:
            if iid in _r2m:
                # This sequence IS the cluster representative
                members = _r2m[iid]
                cluster_rep = 'self'
                cluster_size = str(len(members) + 1)
                clustered_members = ';'.join(sorted(members))
            elif iid in _m2r:
                # This sequence was collapsed into a representative
                cluster_rep = _m2r[iid]
                # find total size from rep's member list
                members = _r2m.get(cluster_rep, [])
                cluster_size = str(len(members) + 1)
                clustered_members = ''
            else:
                # Singleton — no cluster
                cluster_rep = 'self'
                cluster_size = '1'
                clustered_members = ''
        else:
            cluster_rep = 'N/A'
            cluster_size = '1'
            clustered_members = ''

        if _tree is not None:
            in_tree = 'Yes' if iid in _tree else 'No'
        elif cluster_rep not in ('N/A', 'self'):
            in_tree = 'No'  # was collapsed; rep is in tree
        else:
            in_tree = 'Yes'  # if no collapse, all sequences enter the tree

        # Resolve None values returned by iter_classification_rows to 'NA'.
        _conf = c.get('confidence')
        _ident = c.get('identity')
        _best = c.get('best')
        _tax = c.get('tax')
        reference_hit = n.get('reference_nearest_hit', 'NA')
        reference_hit_taxonomy = nearest_meta.get(str(reference_hit), {}).get('taxonomy', 'NA')
        if reference_hit not in (None, '', 'NA', 'None') and _best is not None and str(reference_hit) == str(_best):
            reference_hit_taxonomy = _tax if _tax is not None else 'NA'

        # Some tree builders can omit an exact duplicate even though every
        # isolate remains in the assessment and project database.
        if in_tree == 'No' and cluster_rep == 'self':
            cluster_rep = 'duplicate'

        hungate_hit = n.get('hungate_nearest_hit', 'NA')
        secondary_hit = n.get('secondary_baseline_nearest_hit', 'NA')
        cultured_hit = n.get('cultured_project_nearest_hit', n.get('nearest_hit', 'NA'))
        hungate_hit_taxonomy = nearest_meta.get(str(hungate_hit), {}).get('taxonomy', 'NA')
        secondary_hit_taxonomy = nearest_meta.get(str(secondary_hit), {}).get('taxonomy', 'NA')
        cultured_hit_taxonomy = nearest_meta.get(str(cultured_hit), {}).get('taxonomy', 'NA')
        baseline_hit_taxonomy = cultured_hit_taxonomy
        baseline_agreement_rank, baseline_taxonomy_conflict = _taxonomy_agreement(
            _tax,
            baseline_hit_taxonomy,
        )

        # Make it unmistakable when a candidate's only close match is another
        # submitted project isolate rather than genuine baseline/reference evidence —
        # otherwise an apparently "not novel" candidate could just be resembling an
        # unconfirmed sibling submission.
        closest_project_isolate = n.get('closest_project_isolate', 'NA')
        project_isolate_only_match = str(n.get('project_isolate_only_match', 'False')).strip().lower() == 'true'
        placement_flags = w.get('flags', '')
        if project_isolate_only_match and closest_project_isolate not in (None, '', 'NA', 'None'):
            extra_flag = f'ONLY_MATCHES_PROJECT_ISOLATE:{closest_project_isolate}'
            placement_flags = ';'.join(f for f in (placement_flags, extra_flag) if f)

        assessment_rows.append({
            'id': iid,
            # Immutable partner-catalogue identity; classifier evidence is
            # deliberately kept in separate fields below.
            'partner_taxon': partner_taxon or 'NA',
            'partner_taxon_rank': partner_taxon_rank or 'NA',
            'partner_taxonomy_source': partner_project_meta.get('partner_taxonomy_source_file', 'PARTNER_METADATA') if partner_taxon else 'NOT_PROVIDED',
            'partner_pathway_domain': domain_designation(partner_taxon) or 'NOT_SPECIFIED',
            'partner_taxon_raw': partner_taxon or 'NA',
            'partner_id': partner_meta.get('partner_id', n.get('partner_id', 'NA')),
            'gtdb_taxonomy': gtdb_genome_by_id.get(str(iid), {}).get('taxonomy', ''),
            'gtdb_assignment_status': (
                'CONFIDENT' if gtdb_genome_by_id.get(str(iid), {}).get('taxonomy')
                and str(gtdb_genome_by_id.get(str(iid), {}).get('qc_pass', '')).strip().casefold()
                    in {'1', 'true', 'yes', 'pass', 'passed'}
                else str(gtdb_genome_by_id.get(str(iid), {}).get('assignment_status') or 'PROVISIONAL')
                if gtdb_genome_by_id.get(str(iid), {}).get('taxonomy')
                else 'UNRESOLVED'
            ),
            # ClassificationHit = best vsearch hit in the REFERENCE DB (GTDB/SILVA) → drives Taxonomy
            'taxonomy': _tax if _tax is not None else 'NA',
            'taxonomy_assignment_label': (
                f"{_tax} [{str(c.get('assignment_status') or 'UNRESOLVED').upper()}]"
                if _tax not in (None, '', 'NA') else 'UNRESOLVED'
            ),
            'classification_hit': _best if _best is not None else 'NA',
            'classification_identity': _ident if _ident is not None else 'NA',
            'classification_confidence': _conf if _conf is not None else 'NA',
            'classification_assignment_status': c.get('assignment_status', 'NA'),
            'classification_taxonomy_support_pct': c.get('taxonomy_support_pct', 'NA'),
            'classification_taxonomy_lead_pct': c.get('taxonomy_lead_pct', 'NA'),
            'classification_assigned_rank': c.get('assigned_rank', 'NA'),
            'classification_evidence_tier': c.get('evidence_tier', 'NA'),
            'classification_taxonomy_candidate_count': c.get('taxonomy_candidate_count', 'NA'),
            'classification_alignment_evidence_status': c.get('alignment_evidence_status', 'NA'),
            'classification_assignment_evidence_level': c.get('assignment_evidence_level', 'NA'),
            'classification_evidence_interpretation': c.get('evidence_interpretation', 'NA'),
            'classification_query_coverage': c.get('query_coverage', 'NA'),
            'classification_target_coverage': c.get('target_coverage', 'NA'),
            'classification_alignment_length': c.get('alignment_length', 'NA'),
            'classification_query_length': c.get('query_length', 'NA'),
            'classification_target_length': c.get('target_length', 'NA'),
            'classification_mismatches': c.get('mismatches', 'NA'),
            'classification_gaps': c.get('gaps', 'NA'),
            'baseline_evidence_class': n.get('baseline_evidence_class', 'NA'),
            'hungate_nearest_hit': hungate_hit,
            'hungate_nearest_hit_dataset': nearest_meta.get(str(hungate_hit), {}).get('dataset', 'NA'),
            'hungate_nearest_hit_taxonomy': hungate_hit_taxonomy,
            'hungate_nearest_identity': n.get('hungate_nearest_identity', 'NA'),
            'hungate_nearest_query_coverage': n.get('hungate_nearest_query_coverage', 'NA'),
            'hungate_nearest_alignment_length': n.get('hungate_nearest_alignment_length', 'NA'),
            'hungate_novelty_score': n.get('hungate_novelty_score', 'NA'),
            'hungate_density_source': n.get('hungate_density_source', 'NA'),
            'secondary_baseline_nearest_hit': secondary_hit,
            'secondary_baseline_nearest_hit_dataset': nearest_meta.get(str(secondary_hit), {}).get('dataset', 'NA'),
            'secondary_baseline_nearest_hit_taxonomy': secondary_hit_taxonomy,
            'secondary_baseline_nearest_identity': n.get('secondary_baseline_nearest_identity', 'NA'),
            'secondary_baseline_nearest_query_coverage': n.get('secondary_baseline_nearest_query_coverage', 'NA'),
            'secondary_baseline_nearest_alignment_length': n.get('secondary_baseline_nearest_alignment_length', 'NA'),
            'secondary_baseline_novelty_score': n.get('secondary_baseline_novelty_score', 'NA'),
            'secondary_baseline_density_source': n.get('secondary_baseline_density_source', 'NA'),
            'cultured_project_nearest_hit': cultured_hit,
            'cultured_project_nearest_hit_dataset': nearest_meta.get(str(cultured_hit), {}).get('dataset', 'NA'),
            'cultured_project_nearest_hit_taxonomy': cultured_hit_taxonomy,
            'cultured_project_nearest_identity': n.get('cultured_project_nearest_identity', 'NA'),
            'cultured_project_nearest_query_coverage': n.get('cultured_project_nearest_query_coverage', 'NA'),
            'cultured_project_nearest_alignment_length': n.get('cultured_project_nearest_alignment_length', 'NA'),
            'cultured_project_novelty_score': n.get('cultured_project_novelty_score', 'NA'),
            'cultured_project_density_source': n.get('cultured_project_density_source', 'NA'),
            # Older Baseline* output aliases mirror the combined cultured-project pool.
            'nearest_hit': n.get('nearest_hit', 'NA'),
            'nearest_hit_dataset': nearest_meta.get(str(n.get('nearest_hit')), {}).get('dataset', 'NA'),
            'nearest_hit_taxonomy': baseline_hit_taxonomy,
            'baseline_taxonomy_agreement_rank': baseline_agreement_rank,
            'baseline_taxonomy_conflict': baseline_taxonomy_conflict,
            'nearest_identity': n.get('nearest_identity', 'NA'),
            'nearest_query_coverage': n.get('nearest_query_coverage', 'NA'),
            'nearest_alignment_length': n.get('nearest_alignment_length', 'NA'),
            'novel': n.get('novel', 'NA'),
            'matches_ge_99': n.get('matches_ge_99', 'NA'),
            'matches_ge_97': n.get('matches_ge_97', 'NA'),
            'matches_ge_95': n.get('matches_ge_95', 'NA'),
            'novelty_score': n.get('novelty_score', 'NA'),
            'crowding': n.get('crowding', 'NA'),
            'sequencing_priority': n.get('sequencing_priority', 'NA'),
            'density_source': n.get('density_source', 'NA'),
            'project_nearest_hit': n.get('project_nearest_hit', 'NA'),
            'project_nearest_hit_dataset': nearest_meta.get(str(n.get('project_nearest_hit')), {}).get('dataset', 'NA'),
            'project_nearest_hit_taxonomy': nearest_meta.get(str(n.get('project_nearest_hit')), {}).get('taxonomy', 'NA'),
            'project_nearest_identity': n.get('project_nearest_identity', 'NA'),
            'project_nearest_query_coverage': n.get('project_nearest_query_coverage', 'NA'),
            'project_nearest_alignment_length': n.get('project_nearest_alignment_length', 'NA'),
            'project_novel': n.get('project_novel', 'NA'),
            'project_matches_ge_99': n.get('project_matches_ge_99', 'NA'),
            'project_matches_ge_97': n.get('project_matches_ge_97', 'NA'),
            'project_matches_ge_95': n.get('project_matches_ge_95', 'NA'),
            'project_novelty_score': n.get('project_novelty_score', 'NA'),
            'project_crowding': n.get('project_crowding', 'NA'),
            'project_sequencing_priority': n.get('project_sequencing_priority', 'NA'),
            'project_density_source': n.get('project_density_source', 'NA'),
            'reference_nearest_hit': reference_hit,
            'reference_nearest_hit_taxonomy': reference_hit_taxonomy,
            'reference_nearest_identity': n.get('reference_nearest_identity', 'NA'),
            'reference_nearest_query_coverage': n.get('reference_nearest_query_coverage', 'NA'),
            'reference_nearest_alignment_length': n.get('reference_nearest_alignment_length', 'NA'),
            'reference_novel': n.get('reference_novel', 'NA'),
            'reference_matches_ge_99': n.get('reference_matches_ge_99', 'NA'),
            'reference_matches_ge_97': n.get('reference_matches_ge_97', 'NA'),
            'reference_matches_ge_95': n.get('reference_matches_ge_95', 'NA'),
            'reference_novelty_score': n.get('reference_novelty_score', 'NA'),
            'reference_crowding': n.get('reference_crowding', 'NA'),
            'reference_sequencing_priority': n.get('reference_sequencing_priority', 'NA'),
            'reference_density_source': n.get('reference_density_source', 'NA'),
            'partner_id': n.get('partner_id', 'NA'),
            'selected_for_genome_sequencing': n.get('selected_for_genome_sequencing', 'NA'),
            'already_sequenced': n.get('already_sequenced', 'NA'),
            'nearest_genome_hit': n.get('nearest_genome_hit', 'NA'),
            'nearest_genome_identity': n.get('nearest_genome_identity', 'NA'),
            'closest_project_isolate': n.get('closest_project_isolate', 'NA'),
            'closest_project_isolate_dataset': nearest_meta.get(str(n.get('closest_project_isolate')), {}).get('dataset', 'NA'),
            'closest_project_isolate_identity': n.get('closest_project_isolate_identity', 'NA'),
            'project_isolate_only_match': n.get('project_isolate_only_match', 'NA'),
            'genome_collection_matches_ge_99': n.get('genome_collection_matches_ge_99', 'NA'),
            'genome_collection_matches_ge_97': n.get('genome_collection_matches_ge_97', 'NA'),
            'genome_collection_matches_ge_95': n.get('genome_collection_matches_ge_95', 'NA'),
            'related_genome_clade_ge_97': n.get('related_genome_clade_ge_97', 'NA'),
            'genome_committed_count_same_species': n.get('genome_committed_count_same_species', 'NA'),
            'hungate_genome_count_same_species': n.get('hungate_genome_count_same_species', 'NA'),
            'secondary_baseline_genome_count_same_species': n.get('secondary_baseline_genome_count_same_species', 'NA'),
            'cultured_project_genome_count_same_species': n.get('cultured_project_genome_count_same_species', 'NA'),
            'genome_available_count_same_species': n.get('genome_available_count_same_species', 'NA'),
            'genome_selected_count_same_species': n.get('genome_selected_count_same_species', 'NA'),
            'genome_pending_count_same_species': n.get('genome_pending_count_same_species', 'NA'),
            'pangenome_target': n.get('pangenome_target', '9'),
            'pangenome_gap': n.get('pangenome_gap', 'NA'),
            'genome_sequencing_metadata_source': n.get('genome_sequencing_metadata_source', 'NA'),
            'placement_flags': placement_flags,
            'in_tree': in_tree,
            'cluster_representative': cluster_rep,
            'cluster_size': cluster_size,
            'clustered_members': clustered_members,
            # Phylogenetic isolation + composite investigation score
            # (populated by cluster_report.generate_cluster_reports; default 'NA' here)
            'phylo_isolation': 'NA',
            'investigation_score': 'NA',
            # Alt-db taxonomy (keys: 'alt_tax_<refdb>', 'alt_class_hit_<refdb>', etc.)
            **_build_alt_tax_cols(
                iid, alt_taxonomies or {}, alt_ref_dbs or [],
                alt_classification_evidence or {},
            ),
        })
    try:
        maldi_rows = db.get_maldi_identifications()
    except Exception:
        maldi_rows = []
    by_rank_and_name = {}
    for maldi in maldi_rows:
        rank = str(maldi.get('taxonomic_rank') or '').strip().lower()
        taxon_name = str(maldi.get('taxon_name') or '').strip()
        sequence_id = str(maldi.get('sequence_id') or '').strip()
        key = normalise_taxon_name(taxon_name)
        if rank in {'s', 'g', 'f'} and key and sequence_id:
            by_rank_and_name.setdefault((rank, key), set()).add(sequence_id)
    rank_names = {'s': 'species', 'g': 'genus', 'f': 'family'}
    for row in assessment_rows:
        taxonomy = parse_taxon_string(str(row.get('taxonomy') or ''))
        context = ''
        for rank in ('s', 'g', 'f'):
            rank_value = taxonomy.get(rank)
            if not rank_value:
                continue
            ids = by_rank_and_name.get((rank, normalise_taxon_name(rank_value)), set())
            ids = ids - {str(row.get('id') or '').strip()}
            if ids:
                visible_ids = sorted(ids)[:5]
                suffix = f'; +{len(ids) - 5} more' if len(ids) > 5 else ''
                context = (
                    f'MALDI-TOF soft context: {len(ids)} partner isolate(s) reported as '
                    f'{rank_value} ({rank_names[rank]}; IDs {", ".join(visible_ids)}{suffix}); '
                    'not marker-confirmed or genome evidence'
                )
                break
        row['maldi_taxon_context'] = context
    try:
        sequencing_metadata = db.get_sequencing_metadata_for_ids(run_ids)
    except Exception:
        sequencing_metadata = {}
    selection_context_fields = (
        'addition_method', 'host_species', 'host_breed', 'host_age_weeks',
        'production_type', 'sample_location', 'ruminant_diet', 'maldi_instrument',
        'maldi_library', 'maldi_taxon_raw',
    )
    for row in assessment_rows:
        metadata = sequencing_metadata.get(str(row.get('id') or ''), {}).get('project_metadata', {})
        for field in selection_context_fields:
            row[field] = metadata.get(field) or 'NA'
    return assessment_rows


def _build_alt_tax_cols(
    seq_id: str,
    alt_taxonomies: Dict[str, Dict[str, tuple]],
    alt_ref_dbs: List[str],
    alt_classification_evidence: Optional[Dict[str, Dict[str, dict]]] = None,
) -> dict:
    """Return per-database taxonomy and alignment evidence columns."""
    import re as _re
    cols: dict = {}
    per_seq = alt_taxonomies.get(seq_id, {})
    evidence_by_db = (alt_classification_evidence or {}).get(seq_id, {})
    for ref_db in alt_ref_dbs:
        safe = _re.sub(r'[^A-Za-z0-9]+', '_', ref_db).strip('_')
        entry = per_seq.get(ref_db, (None, None, None, None))
        tax, conf, hit, ident = entry if len(entry) == 4 else (None, None, None, None)
        evidence = evidence_by_db.get(ref_db, {})
        cols[f'alt_tax_{safe}'] = tax if tax is not None else 'NA'
        cols[f'alt_class_hit_{safe}'] = hit if hit is not None else 'NA'
        cols[f'alt_ident_{safe}'] = f"{float(ident):.2f}" if ident is not None else 'NA'
        source_confidence = evidence.get('source_classifier_confidence')
        cols[f'alt_conf_{safe}'] = f"{float(source_confidence):.4f}" if source_confidence is not None else 'NA'
        for field in (
            'query_coverage', 'target_coverage', 'alignment_length',
            'assignment_status', 'taxonomy_support_pct', 'taxonomy_candidate_count',
            'evidence_tier', 'alignment_evidence_status',
            'assignment_evidence_level', 'evidence_interpretation',
        ):
            value = evidence.get(field)
            cols[f'alt_{field}_{safe}'] = value if value not in (None, '') else 'NA'
    return cols


def write_sequence_assessment_tsv(path: str | Path, rows, assessment_db_name: str = 'GTDB'):
    """Write the complete, component-based per-isolate assessment audit."""
    p = Path(path)
    rows = list(rows)
    label = 'GTDB' if 'gtdb' in str(assessment_db_name).lower() else ''.join(
        char if char.isalnum() else '_' for char in str(assessment_db_name or 'Assessment')
    ).strip('_') or 'Assessment'
    alt_ref_safes: List[str] = []
    include_mwl = False
    for row in rows:
        alt_ref_safes = sorted(
            k[len('alt_tax_'):]
            for k in row
            if k.startswith('alt_tax_')
        )
        include_mwl = any(k.startswith('mwl_') or k == 'evaluation_score' for k in row)
        break
    fields = [
        ('ID', 'id'), ('PartnerID', 'partner_id'),
        ('AdditionMethod', 'addition_method'), ('HostSpecies', 'host_species'),
        ('HostBreed', 'host_breed'), ('HostAgeWeeks', 'host_age_weeks'),
        ('ProductionType', 'production_type'), ('SampleLocation', 'sample_location'),
        ('RuminantDiet', 'ruminant_diet'), ('MALDIInstrument', 'maldi_instrument'),
        ('MALDIDatabase', 'maldi_library'), ('MALDITaxonRaw', 'maldi_taxon_raw'),
        (f'{label}Taxonomy', 'taxonomy'),
        ('PartnerTaxon', 'partner_taxon'), ('PartnerTaxonRank', 'partner_taxon_rank'),
        ('PartnerTaxonomySource', 'partner_taxonomy_source'),
        ('PartnerPathwayDomain', 'partner_pathway_domain'),
        ('GTDBGenomeTaxonomy', 'gtdb_taxonomy'), ('GTDBGenomeAssignmentStatus', 'gtdb_assignment_status'),
        ('PangenomeRank', 'pangenome_rank'), ('PangenomeTaxon', 'pangenome_taxon'),
        ('ClassifierGTDBDisagreement', 'classifier_gtdb_disagreement'),
        ('ClassifierDatabaseConcordanceStatus', 'classifier_database_concordance_status'),
        (f'{label}TaxonomyAssignmentLabel', 'taxonomy_assignment_label'),
        (f'{label}ClassificationHit', 'classification_hit'),
        (f'{label}ClassificationIdentity', 'classification_identity'),
        (f'{label}SourceClassifierConfidence', 'classification_confidence'),
        (f'{label}AssignmentStatus', 'classification_assignment_status'),
        (f'{label}TaxonomySupportPct', 'classification_taxonomy_support_pct'),
        (f'{label}TaxonomyLeadPct', 'classification_taxonomy_lead_pct'),
        (f'{label}AssignedRank', 'classification_assigned_rank'),
        (f'{label}EvidenceTier', 'classification_evidence_tier'),
        (f'{label}TaxonomyCandidateCount', 'classification_taxonomy_candidate_count'),
        (f'{label}AlignmentEvidenceStatus', 'classification_alignment_evidence_status'),
        (f'{label}AssignmentEvidenceLevel', 'classification_assignment_evidence_level'),
        (f'{label}EvidenceInterpretation', 'classification_evidence_interpretation'),
        (f'{label}QueryCoverage', 'classification_query_coverage'),
        (f'{label}TargetCoverage', 'classification_target_coverage'),
        (f'{label}AlignmentLength', 'classification_alignment_length'),
        (f'{label}QueryLength', 'classification_query_length'),
        (f'{label}TargetLength', 'classification_target_length'),
        (f'{label}Mismatches', 'classification_mismatches'),
        (f'{label}Gaps', 'classification_gaps'),
        ('MALDITaxonContext', 'maldi_taxon_context'),
        ('BaselineEvidenceClass', 'baseline_evidence_class'),
        ('HungateNearestHit', 'hungate_nearest_hit'), ('HungateNearestHitDataset', 'hungate_nearest_hit_dataset'),
        ('HungateNearestHitTaxonomy', 'hungate_nearest_hit_taxonomy'),
        ('HungateNearestIdentity', 'hungate_nearest_identity'),
        ('HungateNearestQueryCoverage', 'hungate_nearest_query_coverage'),
        ('HungateNearestAlignmentLength', 'hungate_nearest_alignment_length'),
        ('HungateNoveltyScore', 'hungate_novelty_score'), ('HungateSource', 'hungate_density_source'),
        ('HungateDensityCountMethod', 'hungate_density_count_method'),
        ('SecondaryBaselineNearestHit', 'secondary_baseline_nearest_hit'),
        ('SecondaryBaselineNearestHitDataset', 'secondary_baseline_nearest_hit_dataset'),
        ('SecondaryBaselineNearestHitTaxonomy', 'secondary_baseline_nearest_hit_taxonomy'),
        ('SecondaryBaselineNearestIdentity', 'secondary_baseline_nearest_identity'),
        ('SecondaryBaselineNearestQueryCoverage', 'secondary_baseline_nearest_query_coverage'),
        ('SecondaryBaselineNearestAlignmentLength', 'secondary_baseline_nearest_alignment_length'),
        ('SecondaryBaselineNoveltyScore', 'secondary_baseline_novelty_score'),
        ('SecondaryBaselineSource', 'secondary_baseline_density_source'),
        ('SecondaryBaselineDensityCountMethod', 'secondary_baseline_density_count_method'),
        ('CulturedProjectNearestHit', 'cultured_project_nearest_hit'),
        ('CulturedProjectNearestHitDataset', 'cultured_project_nearest_hit_dataset'),
        ('CulturedProjectNearestHitTaxonomy', 'cultured_project_nearest_hit_taxonomy'),
        ('CulturedProjectNearestIdentity', 'cultured_project_nearest_identity'),
        ('CulturedProjectNearestQueryCoverage', 'cultured_project_nearest_query_coverage'),
        ('CulturedProjectNearestAlignmentLength', 'cultured_project_nearest_alignment_length'),
        ('CulturedProjectNoveltyScore', 'cultured_project_novelty_score'),
        ('CulturedProjectSource', 'cultured_project_density_source'),
        ('CulturedProjectDensityCountMethod', 'density_count_method'),
        ('BaselineNearestHit', 'nearest_hit'), ('BaselineNearestHitDataset', 'nearest_hit_dataset'),
        ('BaselineNearestHitTaxonomy', 'nearest_hit_taxonomy'),
        ('BaselineTaxonomyAgreementRank', 'baseline_taxonomy_agreement_rank'),
        ('BaselineTaxonomyConflict', 'baseline_taxonomy_conflict'),
        ('BaselineNearestIdentity', 'nearest_identity'), ('BaselineNovel', 'novel'),
        ('BaselineNearestQueryCoverage', 'nearest_query_coverage'),
        ('BaselineNearestAlignmentLength', 'nearest_alignment_length'),
        ('BaselineMatchesGE99', 'matches_ge_99'), ('BaselineMatchesGE97', 'matches_ge_97'),
        ('BaselineMatchesGE95', 'matches_ge_95'), ('BaselineNoveltyScore', 'novelty_score'),
        ('BaselineCrowding', 'crowding'), ('BaselinePriority', 'sequencing_priority'),
        ('BaselineSource', 'density_source'),
        ('BaselineDensityCountMethod', 'density_count_method'),
        ('ProjectNearestHit', 'project_nearest_hit'),
        ('ProjectNearestHitDataset', 'project_nearest_hit_dataset'),
        ('ProjectNearestHitTaxonomy', 'project_nearest_hit_taxonomy'),
        ('ProjectNearestIdentity', 'project_nearest_identity'), ('ProjectNovel', 'project_novel'),
        ('ProjectNearestQueryCoverage', 'project_nearest_query_coverage'),
        ('ProjectNearestAlignmentLength', 'project_nearest_alignment_length'),
        ('ProjectMatchesGE99', 'project_matches_ge_99'), ('ProjectMatchesGE97', 'project_matches_ge_97'),
        ('ProjectMatchesGE95', 'project_matches_ge_95'), ('ProjectNoveltyScore', 'project_novelty_score'),
        ('ProjectCrowding', 'project_crowding'), ('ProjectPriority', 'project_sequencing_priority'),
        ('ProjectSource', 'project_density_source'),
        ('ProjectDensityCountMethod', 'project_density_count_method'),
        (f'{label}ReferenceNearestHit', 'reference_nearest_hit'),
        (f'{label}ReferenceNearestHitTaxonomy', 'reference_nearest_hit_taxonomy'),
        (f'{label}ReferenceNearestIdentity', 'reference_nearest_identity'),
        (f'{label}ReferenceNearestQueryCoverage', 'reference_nearest_query_coverage'),
        (f'{label}ReferenceNearestAlignmentLength', 'reference_nearest_alignment_length'),
        (f'{label}ReferenceNovel', 'reference_novel'),
        (f'{label}ReferenceMatchesGE99', 'reference_matches_ge_99'),
        (f'{label}ReferenceMatchesGE97', 'reference_matches_ge_97'),
        (f'{label}ReferenceMatchesGE95', 'reference_matches_ge_95'),
        (f'{label}ReferenceNoveltyScore', 'reference_novelty_score'),
        (f'{label}ReferenceCrowding', 'reference_crowding'),
        (f'{label}ReferenceSource', 'reference_density_source'),
        (f'{label}ReferenceDensityCountMethod', 'reference_density_count_method'),
        ('AlreadySelectedForGenomeSequencing', 'selected_for_genome_sequencing'),
        ('GenomeAlreadySequenced', 'already_sequenced'),
        ('NearestGenomeHit', 'nearest_genome_hit'), ('NearestGenomeIdentity', 'nearest_genome_identity'),
        ('ClosestProjectIsolate', 'closest_project_isolate'),
        ('ClosestProjectIsolateDataset', 'closest_project_isolate_dataset'),
        ('ClosestProjectIsolateIdentity', 'closest_project_isolate_identity'),
        ('ProjectIsolateOnlyMatch', 'project_isolate_only_match'),
        ('GenomeCollectionMatchesGE99', 'genome_collection_matches_ge_99'),
        ('GenomeCollectionMatchesGE97', 'genome_collection_matches_ge_97'),
        ('GenomeCollectionMatchesGE95', 'genome_collection_matches_ge_95'),
        ('RelatedGenomeCladeGE97', 'related_genome_clade_ge_97'),
        ('HungateGenomesSameAssessmentSpecies', 'hungate_genome_count_same_species'),
        ('SecondaryBaselineGenomesSameAssessmentSpecies', 'secondary_baseline_genome_count_same_species'),
        ('CulturedProjectGenomesSameAssessmentSpecies', 'cultured_project_genome_count_same_species'),
        ('BaselineGenomesSameAssessmentSpecies', 'genome_available_count_same_species'),
        ('SequencedPartnerGenomesSameAssessmentSpecies', 'genome_selected_count_same_species'),
        ('SelectedPendingGenomesSameAssessmentSpecies', 'genome_pending_count_same_species'),
        ('CommittedGenomesSameAssessmentSpecies', 'genome_committed_count_same_species'),
        ('PangenomeTarget', 'pangenome_target'), ('PangenomeGap', 'pangenome_gap'),
        ('PangenomeEligibilityStatus', 'pangenome_eligibility_status'),
        ('PangenomeEligibilityReason', 'pangenome_eligibility_reason'),
        ('PangenomeCommittedGenomeCount', 'pangenome_committed_count'),
        ('PangenomeBaselineGenomeCount', 'pangenome_baseline_count'),
        ('PangenomeCommittedGenomeIDs', 'pangenome_committed_ids'),
        ('GenomeCollectionSource', 'genome_sequencing_metadata_source'),
        ('MarkerQCClass', 'marker_qc_class'),
        ('MarkerQCRecommendation', 'marker_qc_recommendation'),
        ('MarkerManualReviewStatus', 'marker_manual_review_status'),
        ('MarkerQCFlag', 'marker_qc_flag'),
        ('MarkerQCReasons', 'marker_qc_reasons'),
        ('MarkerSourceManifest', 'marker_source_manifest'),
        ('ChimeraCall', 'chimera_call'), ('UCHIMEScore', 'chimera_score'),
    ]
    if include_mwl:
        fields.extend([
            ('MWLMatch', 'mwl_match'), ('MWLID', 'mwl_id'),
            ('MWLMatchedRank', 'mwl_matched_rank'), ('MWLMatchedTaxon', 'mwl_matched_taxon'),
            ('MWLTaxonomicScore', 'mwl_taxonomic_score'), ('MWLIdentity', 'mwl_identity'),
            ('MWLScore', 'mwl_score'), ('MWLRole', 'mwl_role'),
        ])
    fields.extend([
        ('SequencingSetID', 'sequencing_set_id'), ('SequencingSetRole', 'sequencing_set_role'),
        ('SequencingSetRank', 'sequencing_set_rank'),
        ('SelectionDiversityDistance', 'selection_diversity_distance'),
        ('SelectionGroupType', 'selection_group_type'),
        ('BaselineRedundancyStatus', 'baseline_redundancy_status'),
        ('BaselineRedundancyIdentityThreshold', 'baseline_redundancy_identity_threshold'),
        ('BaselineRedundancyMinQueryCoverage', 'baseline_redundancy_min_query_coverage'),
        ('BaselineExtensionStatus', 'baseline_extension_status'),
        ('BaselineExtensionMinIdentity', 'baseline_extension_min_identity'),
        ('BaselineExtensionMinQueryCoverage', 'baseline_extension_min_query_coverage'),
        ('SequencingSetReason', 'sequencing_set_reason'),
        ('SelectionGroupBasis', 'selection_group_basis'), ('SelectionGroupTaxon', 'selection_group_taxon'),
        ('NovelLookingClade', 'novel_looking_clade'),
        ('PhyloIsolation', 'phylo_isolation'), ('LocalNeighbourhoodFigure', 'local_neighbourhood_figure'),
        ('LocalPairwisePidentTable', 'local_pairwise_pident_table'),
        ('TreeContextLeafCount', 'tree_context_leaf_count'),
        ('AssessedSequencesInTreeContext', 'assessed_sequences_in_tree_context'),
        ('InTree', 'in_tree'), ('ClusterRepresentative', 'cluster_representative'),
        ('ClusterSize', 'cluster_size'), ('ClusteredMembers', 'clustered_members'),
        ('PlacementFlags', 'placement_flags'),
    ])
    for safe in alt_ref_safes:
        fields.extend([
            (f'Taxonomy_{safe}', f'alt_tax_{safe}'),
            (f'ClassificationHit_{safe}', f'alt_class_hit_{safe}'),
            (f'Identity_{safe}', f'alt_ident_{safe}'),
            (f'SourceClassifierConfidence_{safe}', f'alt_conf_{safe}'),
            (f'QueryCoverage_{safe}', f'alt_query_coverage_{safe}'),
            (f'TargetCoverage_{safe}', f'alt_target_coverage_{safe}'),
            (f'AlignmentLength_{safe}', f'alt_alignment_length_{safe}'),
            (f'AssignmentStatus_{safe}', f'alt_assignment_status_{safe}'),
            (f'TaxonomySupportPct_{safe}', f'alt_taxonomy_support_pct_{safe}'),
            (f'TaxonomyCandidateCount_{safe}', f'alt_taxonomy_candidate_count_{safe}'),
            (f'EvidenceTier_{safe}', f'alt_evidence_tier_{safe}'),
            (f'AlignmentEvidenceStatus_{safe}', f'alt_alignment_evidence_status_{safe}'),
            (f'AssignmentEvidenceLevel_{safe}', f'alt_assignment_evidence_level_{safe}'),
            (f'EvidenceInterpretation_{safe}', f'alt_evidence_interpretation_{safe}'),
        ])

    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t', lineterminator='\n')
        writer.writerow([header for header, _ in fields])
        for row in rows:
            writer.writerow([row.get(key, 'NA') for _, key in fields])
    return str(p)


def _as_float(value, default=None):
    try:
        if value in (None, '', 'NA', 'None'):
            return default
        return float(value)
    except Exception:
        return default


def _as_int(value, default=None):
    try:
        if value in (None, '', 'NA', 'None'):
            return default
        return int(float(value))
    except Exception:
        return default


def _pool_available(row: dict, source_key: str) -> bool:
    return str(row.get(source_key) or '').lower() not in ('', 'na', 'none')


def _cultured_gap(row: dict) -> str:
    identity = _as_float(row.get('nearest_identity'))
    if not _pool_available(row, 'density_source') or identity is None:
        return 'UNKNOWN - no cultured baseline'
    if identity < 97.0:
        return 'LARGE - no match >=97%'
    if identity < 98.65:
        return 'MODERATE - nearest 97-98.65%'
    return 'SMALL - close match >=98.65%'


def _project_coverage(row: dict) -> str:
    if not _pool_available(row, 'project_density_source'):
        return 'UNKNOWN - no project collection'
    count = _as_int(row.get('project_matches_ge_97'), 0)
    if count == 0:
        return 'UNCOVERED - 0 neighbours >=97%'
    if count <= 3:
        return f'SPARSE - {count} neighbour(s) >=97%'
    if count <= 10:
        return f'MODERATE - {count} neighbours >=97%'
    return f'DENSE - {count} neighbours >=97%'


def _reference_context(row: dict) -> str:
    identity = _as_float(row.get('reference_nearest_identity'))
    if not _pool_available(row, 'reference_density_source') or identity is None:
        return 'UNKNOWN - no external reference comparison'
    if identity < 94.5:
        return 'HIGH DIVERGENCE - <94.5%'
    if identity < 98.65:
        return 'MODERATE DIVERGENCE - 94.5-98.65%'
    return 'CLOSE REFERENCE - >=98.65%'


def _genome_coverage(row: dict) -> str:
    if str(row.get('already_sequenced') or '').lower() == 'true':
        return 'CURRENT ISOLATE ALREADY SEQUENCED'
    if str(row.get('selected_for_genome_sequencing') or '').lower() == 'true':
        return 'CURRENT ISOLATE SELECTED - GENOME PENDING'
    committed = _as_int(row.get('pangenome_committed_count'))
    target = _as_int(row.get('pangenome_target'), 9)
    gap = _as_int(row.get('pangenome_gap'))
    rank = str(row.get('pangenome_rank') or 'catalogue-rank').lower()
    taxon = str(row.get('pangenome_taxon') or '').strip()
    group_label = f'{rank} {taxon}'.strip()
    if committed is not None and gap is not None:
        if gap <= 0:
            return f'TARGET MET - {committed}/{target} {group_label} genomes'
        return f'PANGENOME GAP - {committed}/{target} {group_label} genomes committed; {gap} required'
    identity = _as_float(row.get('nearest_genome_identity'))
    if identity is None or identity <= 0:
        return 'NOT REPRESENTED'
    if identity >= 98.65:
        return f'NEAR-DUPLICATE AVAILABLE GENOME - {identity:.2f}%'
    if identity >= 97.0:
        return f'RELATED AVAILABLE GENOME - {identity:.2f}%'
    return f'NOT CLOSELY REPRESENTED - nearest {identity:.2f}%'


def _species_context(row: dict) -> str:
    pangenome_gap = _as_int(row.get('pangenome_gap'))
    if pangenome_gap is None or pangenome_gap <= 0:
        return ''
    nearest = _as_float(row.get('nearest_genome_identity'))
    if nearest is None or nearest <= 0:
        return ''
    cluster_size = _as_int(row.get('cluster_size'))
    project_identity = _as_float(row.get('project_nearest_identity'))
    project_hits = _as_int(row.get('project_matches_ge_97'))
    if nearest >= 98.65:
        if cluster_size is not None and cluster_size > 1:
            return f'POSSIBLE TAXONOMY SPLIT - {nearest:.2f}% nearest genome within a {cluster_size}-member cluster; labels may split near-identical genomes'
        return f'POSSIBLE TAXONOMY SPLIT - {nearest:.2f}% nearest genome; labels may split near-identical genomes'
    if project_identity is not None and project_identity >= 97.0 and project_hits is not None and project_hits > 0:
        return f'POSSIBLE TAXONOMY SPLIT - {project_hits} project neighbour(s) at >=97%; labels may split near-identical genomes'
    if cluster_size is not None and cluster_size > 1:
        return f'CLUSTER CONTEXT - {cluster_size}-member local cluster'
    return ''


def _evidence_quality(row: dict) -> str:
    return _shared_evidence_quality(row)


def build_selection_decision(row: dict) -> dict:
    """Build transparent, categorical genome-sequencing recommendation evidence."""
    if row.get('no_marker_sequence'):
        reason = str(row.get('marker_qc_reasons') or 'every read failed QC before assembly').strip()
        return {
            'decision': 'RESEQUENCE REQUIRED - NO SEQUENCE',
            'evidence_quality': 'NONE',
            'cultured_gap': 'UNKNOWN - no marker sequence',
            'project_coverage': 'UNKNOWN - no marker sequence',
            'reference_context': 'UNKNOWN - no marker sequence',
            'genome_coverage': 'NOT REPRESENTED - no marker sequence',
            'species_context': '',
            'decision_reason': f'no marker sequence was assembled for this isolate; {reason}',
        }
    cultured_gap = _cultured_gap(row)
    project_coverage = _project_coverage(row)
    reference_context = _reference_context(row)
    genome_coverage = _genome_coverage(row)
    species_context_parts = [
        value for value in (_species_context(row), str(row.get('maldi_taxon_context') or '').strip())
        if value
    ]
    species_context = '; '.join(species_context_parts)
    evidence_quality = _evidence_quality(row)
    selected = genome_coverage == 'CURRENT ISOLATE ALREADY SEQUENCED'
    pending = genome_coverage == 'CURRENT ISOLATE SELECTED - GENOME PENDING'
    pangenome_gap = _as_int(row.get('pangenome_gap'))
    target_met = pangenome_gap == 0
    set_role = str(row.get('sequencing_set_role') or '').upper()
    baseline_redundant = set_role == 'BASELINE_REDUNDANT'
    boundary_review = set_role == 'PANGENOME_BOUNDARY_REVIEW'
    mwl_match = str(row.get('mwl_match') or '').lower() == 'yes'
    mwl_rank = str(row.get('mwl_matched_rank') or '').lower()
    strong_mwl = mwl_match and mwl_rank in ('species', 's', 'genus', 'g', 'family', 'f')
    moderate_mwl = mwl_match and mwl_rank in ('order', 'o')
    diversity_strength = str(row.get('diversity_metadata_strength') or '').strip().upper()
    strong_diversity = diversity_strength in {'HIGH', 'MODERATE'}
    baseline_identity = _as_float(row.get('nearest_identity'))
    project_identity = _as_float(row.get('project_nearest_identity'))
    reference_identity = _as_float(row.get('reference_nearest_identity'))

    positive = []
    caution = []
    if str(row.get('project_isolate_only_match') or '').strip().lower() == 'true':
        isolate = row.get('closest_project_isolate', 'NA')
        caution.insert(0, f'only close match is another project isolate ({isolate}), not a confirmed baseline/reference hit')
    if selected:
        caution.append('current isolate already has a genome available')
    elif pending:
        caution.append('current isolate is already selected; genome pending')
    if strong_mwl:
        matched = row.get('mwl_matched_taxon', 'NA')
        positive.append(f'MWL target match ({matched})')
    elif moderate_mwl:
        positive.append(f'MWL order-level context ({row.get("mwl_matched_taxon", "NA")})')
    if strong_diversity:
        reason = str(row.get('diversity_metadata_reason') or '').strip()
        positive.append('diversity metadata supports an underrepresented isolate' + (f' ({reason})' if reason else ''))
    if baseline_identity is not None and _pool_available(row, 'density_source') and baseline_identity < 97.0:
        positive.append(f'no close cultured isolate ({baseline_identity:.2f}%)')
    if project_identity is not None and _pool_available(row, 'project_density_source') and project_identity < 97.0:
        positive.append(f'no close project isolate ({project_identity:.2f}%)')
    if reference_identity is not None and _pool_available(row, 'reference_density_source') and reference_identity < 98.65:
        positive.append(f'divergent from external reference ({reference_identity:.2f}%)')
    if evidence_quality != 'HIGH':
        caution.append(f'{evidence_quality.lower()} marker evidence')
    if baseline_redundant:
        caution.append('near-identical to a cultured baseline marker at high query coverage')
    if boundary_review:
        status = str(row.get('baseline_extension_status') or 'extension criteria not met')
        caution.append(f'not admitted to the baseline-pangenome extension ({status})')
    if target_met:
        caution.append('GTDB-rank pangenome target already met')
    elif pangenome_gap is not None:
        positive.append(f'{pangenome_gap} genome(s) still required for pangenome target')
    if project_coverage.startswith('DENSE'):
        caution.append('dense prior project neighbourhood')
    if cultured_gap.startswith('SMALL'):
        caution.append('close cultured-baseline match')
    if cultured_gap.startswith('UNKNOWN') and project_coverage.startswith('UNKNOWN'):
        caution.append('baseline and prior-project comparisons unavailable')

    if selected:
        decision = 'ALREADY SEQUENCED'
    elif pending:
        decision = 'ALREADY SELECTED - GENOME PENDING'
    elif evidence_quality == 'LOW':
        decision = 'REVIEW BEFORE SELECTION'
    elif baseline_redundant:
        decision = 'EXCLUDE - BASELINE REDUNDANT'
    elif boundary_review:
        decision = 'REVIEW - PANGENOME BOUNDARY'
    elif set_role == 'PRIMARY':
        decision = 'PRIORITISE - SET PRIMARY'
    elif set_role == 'BACKUP':
        decision = 'RESERVE - SET BACKUP'
    elif set_role == 'DIVERSITY_CANDIDATE' and evidence_quality == 'HIGH' and (strong_mwl or len(positive) >= 2):
        # Meeting the numerical pangenome target should not suppress a strong
        # strain-diversity candidate. The set role describes why it was ranked;
        # the recommendation must still reflect independent biological evidence.
        decision = 'STRONG CANDIDATE - STRAIN DIVERSITY'
    elif set_role == 'DIVERSITY_CANDIDATE':
        decision = 'SECONDARY - STRAIN DIVERSITY'
    elif target_met:
        decision = 'LOWER PRIORITY - TARGET MET'
    elif set_role == 'ALTERNATE':
        decision = 'SECONDARY CANDIDATE'
    elif strong_mwl or len(positive) >= 2:
        decision = 'STRONG CANDIDATE'
    elif positive:
        decision = 'SECONDARY CANDIDATE'
    elif cultured_gap.startswith('SMALL') and project_coverage.startswith(('MODERATE', 'DENSE')):
        decision = 'LOWER PRIORITY - LIMITED ADDED VALUE'
    elif cultured_gap.startswith('UNKNOWN') and project_coverage.startswith('UNKNOWN'):
        decision = 'REVIEW BEFORE SELECTION'
    else:
        decision = 'SECONDARY CANDIDATE'

    reasons = positive[:3] + caution[:2]
    if not reasons:
        reasons.append('no strong novelty or redundancy signal')
    return {
        'decision': decision,
        'evidence_quality': evidence_quality,
        'cultured_gap': cultured_gap,
        'project_coverage': project_coverage,
        'reference_context': reference_context,
        'genome_coverage': genome_coverage,
        'species_context': species_context,
        'decision_reason': '; '.join(reasons),
    }


def _board_recommendation(row: dict) -> tuple[str, str]:
    """Compatibility helper returning the transparent decision and rationale."""
    support = build_selection_decision(row)
    return support['decision'], support['decision_reason']


def _marker_review_display(row: dict) -> str:
    status = str(row.get('marker_manual_review_status') or 'NOT_REVIEWED').strip().upper()
    qc_class = str(row.get('marker_qc_class') or '').strip().upper()
    recommendation = str(row.get('marker_qc_recommendation') or '').strip().upper()
    if status in {'', 'NA', 'NONE'}:
        status = 'NOT_REVIEWED'
    if (
        status == 'NOT_REVIEWED'
        and qc_class == 'PASS_HIGH_CONFIDENCE'
        and recommendation in {'', 'ACCEPT', 'ACCEPTED'}
    ):
        return 'NOT_REQUIRED'
    return status


def write_selection_summary_tsv(path: str | Path, rows, assessment_db_name: str = 'GTDB'):
    """Write the SAB-facing shortlist of proposed genome-sequencing candidates.

    Only positive or provisional recommendations are included. Isolates already
    sequenced/committed, explicitly excluded as baseline-redundant, below the
    recommendation threshold, or lacking a marker sequence remain in the full
    sequence_assessment.tsv audit but are not presented as SAB candidates.

    The summary deliberately exposes the key biological and taxonomic context
    needed for a decision; detailed classifier and novelty evidence remains in
    sequence_assessment.tsv.
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    headers = [
        'SequenceID',
        'PartnerID',
        'HostSpecies',
        'MWLID',
        'MWLMatchedTaxon',
        'BaselineNearestHit',
        'BaselineNearestHitGTDBTaxonomy',
        'BaselineNearestIdentity',
        'Taxonomy_NCBI',
        'Taxonomy_SILVA',
        'Taxonomy_GG2',
        'MarkerQC',
        'EvidenceQuality',
        'Recommendation',
        'RecommendationReason',
        'SABDecision',
    ]
    # These categories represent an actual proposal or a provisional candidate
    # requiring discussion. Explicit exclusions, already completed commitments,
    # no-marker resequencing cases and lower-priority/no-added-value cases are
    # retained in the comprehensive assessment, not the SAB decision shortlist.
    sab_candidate_decisions = {
        'PRIORITISE - SET PRIMARY',
        'RESERVE - SET BACKUP',
        'STRONG CANDIDATE',
        'STRONG CANDIDATE - STRAIN DIVERSITY',
        'SECONDARY - STRAIN DIVERSITY',
        'SECONDARY CANDIDATE',
        'REVIEW BEFORE SELECTION',
        'REVIEW - PANGENOME BOUNDARY',
    }
    decision_order = {
        'PRIORITISE - SET PRIMARY': 0,
        'RESERVE - SET BACKUP': 1,
        'STRONG CANDIDATE': 2,
        'STRONG CANDIDATE - STRAIN DIVERSITY': 2,
        'SECONDARY - STRAIN DIVERSITY': 3,
        'SECONDARY CANDIDATE': 4,
        'REVIEW BEFORE SELECTION': 5,
        'REVIEW - PANGENOME BOUNDARY': 6,
    }

    def display_value(row, key):
        value = row.get(key)
        if value is None or str(value).strip() in {'', 'None', 'nan'}:
            return 'NA'
        return value

    prepared = []
    for row in rows:
        support = build_selection_decision(row)
        if support['decision'] not in sab_candidate_decisions:
            continue
        prepared.append((support, row))
    prepared.sort(key=lambda item: (
        decision_order.get(item[0]['decision'], 99),
        str(item[1].get('partner_id', '')),
        str(item[1].get('id', '')),
    ))
    with open(p, 'w', encoding='utf-8', newline='') as fh:
        fh.write('\t'.join(headers) + '\n')
        for support, row in prepared:
            values = [
                display_value(row, 'id'),
                display_value(row, 'partner_id'),
                display_value(row, 'host_species'),
                display_value(row, 'mwl_id'),
                display_value(row, 'mwl_matched_taxon'),
                display_value(row, 'nearest_hit'),
                display_value(row, 'nearest_hit_taxonomy'),
                display_value(row, 'nearest_identity'),
                display_value(row, 'alt_tax_NCBI_16S_2026_08'),
                display_value(row, 'alt_tax_SILVA_144'),
                display_value(row, 'alt_tax_GG2_2024_09'),
                display_value(row, 'marker_qc_class') or 'QUALITY_UNVERIFIED',
                support['evidence_quality'],
                support['decision'],
                support['decision_reason'],
                '',  # Deliberately blank: the SAB records its decision in this column.
            ]
            safe = [str(value if value is not None else 'NA').replace('\t', ' ').replace('\n', ' ') for value in values]
            fh.write('\t'.join(safe) + '\n')
    return str(p)


def write_baseline_hits_tsv(path: str | Path, rows):
    """Write a concise tier-aware nearest cultured-baseline hit report."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        'ID', 'BaselineEvidenceClass',
        'HungateNearestHit', 'HungateDataset', 'HungateNearestIdentity',
        'HungateNearestQueryCoverage', 'HungateNearestTaxonomy', 'HungateNoveltyScore',
        'SecondaryBaselineNearestHit', 'SecondaryBaselineDataset', 'SecondaryBaselineNearestIdentity',
        'SecondaryBaselineNearestQueryCoverage', 'SecondaryBaselineNearestTaxonomy', 'SecondaryBaselineNoveltyScore',
        'CulturedProjectNearestHit', 'CulturedProjectDataset', 'CulturedProjectNearestIdentity',
        'CulturedProjectNearestQueryCoverage', 'CulturedProjectNearestTaxonomy', 'CulturedProjectNoveltyScore',
    ]
    with open(p, 'w') as fh:
        fh.write('\t'.join(fields) + '\n')
        for row in rows:
            cultured_hit = row.get('cultured_project_nearest_hit') or row.get('nearest_hit', 'NA')
            if cultured_hit in (None, '', 'NA', 'None'):
                continue
            values = [
                row.get('id', 'NA'), row.get('baseline_evidence_class', 'NA'),
                row.get('hungate_nearest_hit', 'NA'), row.get('hungate_nearest_hit_dataset', 'NA'),
                row.get('hungate_nearest_identity', 'NA'), row.get('hungate_nearest_query_coverage', 'NA'),
                row.get('hungate_nearest_hit_taxonomy', 'NA'), row.get('hungate_novelty_score', 'NA'),
                row.get('secondary_baseline_nearest_hit', 'NA'), row.get('secondary_baseline_nearest_hit_dataset', 'NA'),
                row.get('secondary_baseline_nearest_identity', 'NA'), row.get('secondary_baseline_nearest_query_coverage', 'NA'),
                row.get('secondary_baseline_nearest_hit_taxonomy', 'NA'), row.get('secondary_baseline_novelty_score', 'NA'),
                cultured_hit, row.get('cultured_project_nearest_hit_dataset', row.get('nearest_hit_dataset', 'NA')),
                row.get('cultured_project_nearest_identity', row.get('nearest_identity', 'NA')),
                row.get('cultured_project_nearest_query_coverage', row.get('nearest_query_coverage', 'NA')),
                row.get('cultured_project_nearest_hit_taxonomy', row.get('nearest_hit_taxonomy', 'NA')),
                row.get('cultured_project_novelty_score', row.get('novelty_score', 'NA')),
            ]
            safe = [str(value if value is not None else 'NA').replace('\t', ' ').replace('\n', ' ') for value in values]
            fh.write('\t'.join(safe) + '\n')
    return str(p)
