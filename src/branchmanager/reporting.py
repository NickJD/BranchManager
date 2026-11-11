"""Cumulative project overviews and decision deltas."""

from __future__ import annotations

import csv
import hashlib
import html
import json
import os
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path

from branchmanager.pipeline.workflow_helpers import build_selection_decision


def _dashboard_number(value):
    try:
        if value is None or str(value).strip() in {'', 'NA', 'None', 'nan'}:
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def _dashboard_first_number(*values):
    """Return the first parseable numeric value from alternative TSV fields."""
    for value in values:
        number = _dashboard_number(value)
        if number is not None:
            return number
    return None


def _dashboard_first_text(*values):
    """Return the first meaningful text value, ignoring common missing tokens."""
    for value in values:
        if value is None:
            continue
        text = str(value).strip()
        if text and text.casefold() not in {'na', 'none', 'nan', 'null'}:
            return text
    return ''


def _dashboard_taxon(value: str) -> str:
    """Return the deepest *assigned* taxon from a GTDB lineage.

    GTDB lineages commonly retain empty rank placeholders (for example,
    ``g__Blautia;s__``).  Taking the final semicolon-delimited component
    therefore loses a perfectly valid genus assignment.  Ignore empty and
    placeholder values, and include the rank in the display label so that
    users can distinguish a genus-level assignment from a species-level one.
    """
    if value is None:
        return 'Not assigned'
    text = str(value).strip()
    if not text or text.casefold() in {'na', 'none', 'nan', 'unclassified'}:
        return 'Not assigned'
    rank_labels = {
        'd': 'Domain', 'p': 'Phylum', 'c': 'Class', 'o': 'Order',
        'f': 'Family', 'g': 'Genus', 's': 'Species',
    }
    assigned = []
    for component in text.split(';'):
        component = component.strip()
        if not component:
            continue
        if '__' in component:
            prefix, name = component.split('__', 1)
            rank = prefix.strip().casefold()
            name = name.strip()
        else:
            # Retain compatibility with plain-text taxon values.
            rank, name = '', component
        if not name or name.casefold() in {
            'na', 'none', 'nan', 'unclassified', 'uncultured', 'unknown',
        }:
            continue
        if name in {'-', '___'}:
            continue
        assigned.append((rank_labels.get(rank, rank.upper()), name))
    if not assigned:
        return 'Not assigned'
    rank, name = assigned[-1]
    return f'{rank}: {name}' if rank else name



def write_genome_collection_assessment(db, output_path: str | Path) -> str:
    """Export registered genome records for the Performance Review dashboard.

    This is a reporting export only: it does not alter genome QC or selection
    decisions. PartnerID is taken from the sequence dataset identifier.
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with db.connect() as conn:
        rows = conn.execute(
            "SELECT g.sequence_id, COALESCE(s.dataset,''), g.genome_id, "
            "g.genome_status, g.genome_qc_pass, g.completeness, g.contamination, "
            "g.gtdb_taxonomy, g.genome_path, g.notes "
            "FROM genome_records g LEFT JOIN sequences s ON s.id=g.sequence_id "
            "ORDER BY g.sequence_id, g.genome_id"
        ).fetchall()
    fields = [
        'sequence_id', 'partner_id', 'genome_id', 'genome_status',
        'genome_qc_pass', 'completeness', 'contamination',
        'gtdb_taxonomy', 'genome_path', 'notes',
    ]
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(fields)
        writer.writerows(rows)
    return str(path)


def _read_dashboard_database(database_path: str | Path) -> tuple[dict, list[dict], dict]:
    """Read current sequencing decisions and genome records without mutating SQLite."""
    path = Path(database_path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f'Project database does not exist: {path}')
    states: dict[str, dict] = {}
    genomes: list[dict] = []
    identity = {'schema_version': 'unknown', 'sha256': '', 'path': str(path)}
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    identity['sha256'] = digest.hexdigest()
    uri = f'file:{path}?mode=ro'
    with sqlite3.connect(uri, uri=True, timeout=10) as conn:
        conn.row_factory = sqlite3.Row
        tables = {row[0] for row in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        if 'schema_migrations' in tables:
            version = conn.execute('SELECT MAX(version) FROM schema_migrations').fetchone()[0]
            if version is not None:
                identity['schema_version'] = str(version)
        if 'sequencing_metadata' in tables:
            for row in conn.execute(
                'SELECT id, selected_for_sequencing, selected_for_wgs, operational_status '
                'FROM sequencing_metadata'
            ):
                states[str(row['id'])] = {
                    'selected': bool(row['selected_for_sequencing']),
                    'sequenced': bool(row['selected_for_wgs']),
                    'status': str(row['operational_status'] or ''),
                }
        if 'isolate_status' in tables:
            for row in conn.execute('SELECT sequence_id, status FROM isolate_status'):
                state = states.setdefault(str(row['sequence_id']), {'selected': False, 'sequenced': False, 'status': ''})
                state['status'] = str(row['status'] or state.get('status', ''))
        if 'genome_records' in tables:
            query = (
                "SELECT g.sequence_id, COALESCE(s.dataset,''), g.genome_id, g.genome_status, "
                "g.genome_qc_pass, g.completeness, g.contamination, g.gtdb_taxonomy, g.genome_path, g.notes "
                "FROM genome_records g LEFT JOIN sequences s ON s.id=g.sequence_id "
                "ORDER BY g.sequence_id, g.genome_id"
            )
            for row in conn.execute(query):
                sequence_id = str(row[0] or '')
                genomes.append({
                    'sequence_id': sequence_id, 'partner_id': str(row[1] or ''),
                    'genome_id': str(row[2] or ''), 'genome_status': str(row[3] or ''),
                    'genome_qc_pass': row[4], 'completeness': row[5], 'contamination': row[6],
                    'gtdb_taxonomy': str(row[7] or ''), 'genome_path': str(row[8] or ''),
                    'notes': str(row[9] or ''),
                })
                state = states.setdefault(sequence_id, {'selected': False, 'sequenced': False, 'status': ''})
                state['sequenced'] = True
    for state in states.values():
        if state.get('status') == 'SAB_APPROVED':
            state['selected'] = True
        elif state.get('status') in {'SAB_DECLINED', 'SAB_DEFERRED'}:
            state['selected'] = False
        if state.get('sequenced'):
            state['selected'] = bool(state.get('selected'))
    return states, genomes, identity

def write_performance_review_dashboard(
    selection_summary: str | Path,
    outdir: str | Path,
    *,
    sequence_assessment: str | Path | None = None,
    pangenome_hierarchy: str | Path | None = None,
    novelty_metrics: str | Path | None = None,
    run_manifest: str | Path | None = None,
    genome_assessment: str | Path | None = None,
    database: str | Path | None = None,
) -> str:
    # Read existing TSV/JSON outputs; do not recalculate assessment decisions.
    selection_path = Path(selection_summary)
    output_dir = Path(outdir)
    output_dir.mkdir(parents=True, exist_ok=True)
    assessment_path = Path(sequence_assessment) if sequence_assessment else selection_path.with_name('sequence_assessment.tsv')
    pangenome_path = Path(pangenome_hierarchy) if pangenome_hierarchy else selection_path.with_name('pangenome_hierarchy.tsv')
    novelty_path = Path(novelty_metrics) if novelty_metrics else selection_path.with_name('novelty_metrics.tsv')
    manifest_path = Path(run_manifest) if run_manifest else output_dir / 'run_manifest.json'

    selection_rows = _read_tsv_rows(selection_path)
    database_states, database_genome_rows, database_identity = (
        _read_dashboard_database(database) if database else ({}, [], {})
    )
    assessment_rows = _read_tsv_rows(assessment_path)
    pangenome_rows = _read_tsv_rows(pangenome_path)
    novelty_rows = _read_tsv_rows(novelty_path)
    # Novelty exports have used both ID and SequenceID across releases.
    # Accept either (and tolerate whitespace/BOM-normalised headers) so a
    # valid novelty_metrics.tsv cannot silently render as an empty table.
    novelty_by_id = {}
    for row in novelty_rows:
        novelty_id = str(row.get('ID') or row.get('SequenceID') or row.get('sequence_id') or '').strip()
        if novelty_id:
            novelty_by_id[novelty_id] = row
    # Resolve each candidate to its most specific usable hierarchy row. The
    # hierarchy contains parallel rows at domain through species; taking the
    # first/last arbitrary row can incorrectly attach DOMAIN/Bacteria to every
    # isolate and make its pangenome gap appear to be zero.
    rank_depth = {'DOMAIN': 0, 'PHYLUM': 1, 'CLASS': 2, 'ORDER': 3,
                  'FAMILY': 4, 'GENUS': 5, 'SPECIES': 6}
    pangenome_by_candidate = {}
    for pg in pangenome_rows:
        rank = str(pg.get('PangenomeRank', '') or '').strip().upper()
        taxon = str(pg.get('PangenomeTaxon', '') or '').strip()
        if not rank or not taxon:
            continue
        eligible = str(pg.get('PangenomeEligible', '') or '').strip().lower() == 'yes'
        candidate_ids = str(pg.get('CandidateIDs', '') or pg.get('CandidateID', '') or '')
        for candidate_id in candidate_ids.split(';'):
            candidate_id = candidate_id.strip()
            if not candidate_id:
                continue
            candidate = {
                'rank': rank, 'taxon': taxon,
                'gap': _dashboard_number(pg.get('PangenomeGap', '')),
                'target': _dashboard_number(pg.get('PangenomeTarget', '')),
                'eligible': eligible,
                'committed': _dashboard_number(pg.get('CommittedGenomeCount', '')),
            }
            old = pangenome_by_candidate.get(candidate_id)
            score = (eligible, rank_depth.get(rank, -1))
            old_score = (old['eligible'], rank_depth.get(old['rank'], -1)) if old else None
            if old is None or score > old_score:
                pangenome_by_candidate[candidate_id] = candidate
    # Candidate-level selection decisions take precedence over hierarchy
    # context. The hierarchy may contain exploratory rows at a more specific
    # classifier-derived rank; those must not overwrite a validated
    # partner-catalogue target carried by the selection summary.
    for source_row in selection_rows:
        candidate_id = str(source_row.get('SequenceID', '') or '').strip()
        rank = str(source_row.get('PangenomeRank', '') or '').strip().upper()
        taxon = str(source_row.get('PangenomeTaxon', '') or '').strip()
        if not candidate_id or not rank or not taxon:
            continue
        eligible = str(source_row.get('PangenomeEligible', '') or '').strip().lower() in {'yes', 'true', '1'}
        if eligible:
            pangenome_by_candidate[candidate_id] = {
                'rank': rank, 'taxon': taxon,
                'gap': _dashboard_number(source_row.get('PangenomeGap', '')),
                'target': _dashboard_number(source_row.get('PangenomeTarget', '')),
                'eligible': True,
                'committed': _dashboard_number(source_row.get('PangenomeCommittedGenomeCount', '')),
            }

    assessment_by_id = {
        str(row.get('ID', row.get('SequenceID', ''))): row
        for row in assessment_rows
        if row.get('ID', row.get('SequenceID', ''))
    }

    databases = [
        ('GTDB r232', 'GTDB_r232'),
        ('Genome Taxonomy Database (GG2 2024-09)', 'GG2_2024_09'),
        ('NCBI 16S 2026-08', 'NCBI_16S_2026_08'),
        ('SILVA 144', 'SILVA_144'),
    ]
    # The classifier exports are the authoritative per-database evidence.
    # sequence_assessment.tsv uses a different schema and historically used
    # Taxon_<database>, not Taxonomy_<database>. Load the original tables by
    # isolate ID so the dashboard cannot silently turn valid assignments into
    # "Not available" merely because of a column-name mismatch.
    taxonomy_tables = {}
    taxonomy_dir = output_dir / 'taxonomy'
    for _label, database_key in databases:
        taxonomy_file = taxonomy_dir / f'taxonomy_{database_key}.tsv'
        rows = _read_tsv_rows(taxonomy_file) if taxonomy_file.is_file() else []
        taxonomy_tables[database_key] = {
            str(row.get('ID', '')).strip(): row
            for row in rows if str(row.get('ID', '')).strip()
        }
    records = []
    for source in selection_rows:
        sequence_id = str(source.get('SequenceID', '')).strip()
        if not sequence_id:
            continue
        evidence = assessment_by_id.get(sequence_id, {})
        partner = str(source.get('PartnerID', 'Unknown') or 'Unknown')
        recommendation = str(source.get('Recommendation', 'Not classified') or 'Not classified')
        qc = str(source.get('MarkerQC', evidence.get('MarkerQCClass', 'Not available')) or 'Not available')
        taxonomy = []
        for label, prefix in databases:
            classifier = taxonomy_tables.get(prefix, {}).get(sequence_id, {})
            if prefix == 'GTDB_r232':
                # Prefer the standalone GTDB classifier export, then fall back
                # to fields carried through the assessment summary.
                taxon = classifier.get('Taxon_GTDB_r232', classifier.get('Taxon_GTDB', '')) or evidence.get('GTDBTaxonomy', source.get('GTDBTaxonomy', source.get('Taxonomy', '')))
                status = classifier.get('AssignmentStatus_GTDB_r232', classifier.get('AssignmentStatus_GTDB', '')) or evidence.get('GTDBAssignmentStatus', source.get('GTDBAssignmentStatus', ''))
                identity = classifier.get('Identity_GTDB_r232', classifier.get('Identity_GTDB', '')) or evidence.get('GTDBClassificationIdentity', source.get('GTDBClassificationIdentity', ''))
                qcov = classifier.get('QueryCoverage_GTDB_r232', classifier.get('QueryCoverage_GTDB', '')) or evidence.get('GTDBQueryCoverage', source.get('GTDBQueryCoverage', ''))
                tcov = classifier.get('TargetCoverage_GTDB_r232', classifier.get('TargetCoverage_GTDB', '')) or evidence.get('GTDBTargetCoverage', '')
                explanation = classifier.get('EvidenceInterpretation_GTDB_r232', classifier.get('EvidenceInterpretation_GTDB', '')) or evidence.get('GTDBEvidenceInterpretation', source.get('GTDBEvidenceInterpretation', ''))
                rank = classifier.get('AssignedRank_GTDB_r232', classifier.get('AssignedRank_GTDB', '')) or evidence.get('GTDBAssignedRank', '')
            else:
                taxon = classifier.get(f'Taxon_{prefix}', '') or evidence.get(f'Taxon_{prefix}', evidence.get(f'Taxonomy_{prefix}', ''))
                status = classifier.get(f'AssignmentStatus_{prefix}', '') or evidence.get(f'AssignmentStatus_{prefix}', '')
                identity = classifier.get(f'Identity_{prefix}', '') or evidence.get(f'Identity_{prefix}', '')
                qcov = classifier.get(f'QueryCoverage_{prefix}', '') or evidence.get(f'QueryCoverage_{prefix}', '')
                tcov = classifier.get(f'TargetCoverage_{prefix}', '') or evidence.get(f'TargetCoverage_{prefix}', '')
                explanation = classifier.get(f'EvidenceInterpretation_{prefix}', '') or evidence.get(f'EvidenceInterpretation_{prefix}', '')
                rank = classifier.get(f'AssignedRank_{prefix}', '') or evidence.get(f'AssignedRank_{prefix}', '')
            taxonomy.append({
                'database': label,
                'status': str(status or 'Not available'),
                'taxon': _dashboard_taxon(str(taxon or '')),
                'lineage': str(taxon or ''),
                'identity': _dashboard_number(identity),
                'query_coverage': _dashboard_number(qcov),
                'reference_coverage': _dashboard_number(tcov),
                'rank': str(rank or ''),
                'explanation': str(explanation or ''),
            })

        novelty_source = novelty_by_id.get(sequence_id, {})
        novelty_comparisons = []
        for key, label in [
            ('Hungate', 'Hungate cultured baseline'),
            ('SecondaryBaseline', 'Secondary cultured baseline'),
            ('CulturedProject', 'Combined cultured collection'),
            ('Project', 'Partner project marker collection'),
            ('Reference', 'GTDB reference genomes'),
        ]:
            novelty_comparisons.append({
                'key': key, 'label': label,
                'nearest_identity': _dashboard_number(novelty_source.get(f'{key}NearestIdentity', '')),
                'nearest_hit': str(novelty_source.get(f'{key}NearestHit', '') or ''),
                'matches_99': _dashboard_number(novelty_source.get(f'{key}MatchesGE99', '')),
                'matches_97': _dashboard_number(novelty_source.get(f'{key}MatchesGE97', '')),
                'matches_95': _dashboard_number(novelty_source.get(f'{key}MatchesGE95', '')),
                'novelty_score': _dashboard_number(novelty_source.get(f'{key}NoveltyScore', '')),
                'crowding': str(novelty_source.get(f'{key}Crowding', '') or ''),
                'priority': str(novelty_source.get(f'{key}SequencingPriority', '') or ''),
                'query_coverage': _dashboard_number(novelty_source.get(f'{key}NearestQueryCoverage', '')),
                'alignment_length': _dashboard_number(novelty_source.get(f'{key}NearestAlignmentLength', '')),
                'density_source': str(novelty_source.get(f'{key}DensitySource', '') or ''),
                'density_method': str(novelty_source.get(f'{key}DensityCountMethod', '') or ''),
                'novel_flag': str(novelty_source.get(f'{key}Novel', '') or ''),
            })
        novelty_by_key = {item['key']: item for item in novelty_comparisons}
        cultured_project_match = novelty_by_key.get('CulturedProject', {})
        hungate_match = novelty_by_key.get('Hungate', {})
        secondary_baseline_match = novelty_by_key.get('SecondaryBaseline', {})
        project_match = novelty_by_key.get('Project', {})

        records.append({
            'id': sequence_id,
            'partner': partner,
            'recommendation': recommendation,
            'sequencing_role': str(source.get('SequencingSetRole', evidence.get('SequencingSetRole', '')) or ''),
            'marker_qc': qc,
            'marker_review': str(source.get('MarkerReview', evidence.get('MarkerManualReviewStatus', '')) or ''),
            'evidence_quality': str(source.get('EvidenceQuality', '') or ''),
            'partner_taxon': str(source.get('PartnerTaxon', source.get('PartnerExpectedTaxon', evidence.get('PartnerTaxon', ''))) or ''),
            'already_selected_for_genome_sequencing': str(source.get('AlreadySelectedForGenomeSequencing', '') or '').strip().lower() in {'true', '1', 'yes'},
            'genome_already_sequenced': str(source.get('GenomeAlreadySequenced', '') or '').strip().lower() in {'true', '1', 'yes'},
            'sab_decision': str(source.get('SABDecision', '') or ''),
            'gtdb_taxonomy': str(source.get('GTDBTaxonomy', source.get('Taxonomy', evidence.get('GTDBTaxonomy', ''))) or ''),
            'cultured_gap': str(source.get('CulturedGap', '') or ''),
            'baseline_class': str(source.get('BaselineEvidenceClass', evidence.get('BaselineEvidenceClass', '')) or ''),
            # Prefer the explicit combined cultured-pool field. Older and
            # newer exports may populate this while the legacy Baseline*
            # alias is empty or marked NA.
            'baseline_nearest_identity': _dashboard_first_number(
                cultured_project_match.get('nearest_identity'),
                source.get('CulturedProjectNearestIdentity'),
                source.get('BaselineNearestIdentity'),
                evidence.get('cultured_project_nearest_identity'),
                evidence.get('CulturedProjectNearestIdentity'),
                evidence.get('nearest_identity'),
                evidence.get('BaselineNearestIdentity'),
            ),
            'baseline_nearest_hit': _dashboard_first_text(
                cultured_project_match.get('nearest_hit'),
                source.get('CulturedProjectNearestHit'),
                source.get('BaselineNearestHit'),
                evidence.get('cultured_project_nearest_hit'),
                evidence.get('CulturedProjectNearestHit'),
                evidence.get('nearest_hit'),
                evidence.get('BaselineNearestHit'),
            ),
            'project_nearest_identity': _dashboard_first_number(
                project_match.get('nearest_identity'),
                source.get('ProjectNearestIdentity'),
                evidence.get('project_nearest_identity'),
                evidence.get('ProjectNearestIdentity'),
            ),
            'project_nearest_hit': _dashboard_first_text(
                project_match.get('nearest_hit'),
                source.get('ProjectNearestHit'),
                evidence.get('project_nearest_hit'),
                evidence.get('ProjectNearestHit'),
            ),
            'cultured_match_details': [
                {'key': 'Hungate', 'label': 'Hungate cultured baseline',
                 'identity': hungate_match.get('nearest_identity'),
                 'hit': hungate_match.get('nearest_hit'),
                 'query_coverage': hungate_match.get('query_coverage'),
                 'alignment_length': hungate_match.get('alignment_length')},
                {'key': 'SecondaryBaseline', 'label': 'Secondary cultured baseline',
                 'identity': secondary_baseline_match.get('nearest_identity'),
                 'hit': secondary_baseline_match.get('nearest_hit'),
                 'query_coverage': secondary_baseline_match.get('query_coverage'),
                 'alignment_length': secondary_baseline_match.get('alignment_length')},
                {'key': 'CulturedProject', 'label': 'Combined cultured collection',
                 'identity': cultured_project_match.get('nearest_identity'),
                 'hit': cultured_project_match.get('nearest_hit'),
                 'query_coverage': cultured_project_match.get('query_coverage'),
                 'alignment_length': cultured_project_match.get('alignment_length')},
            ],
            'project_match_details': {
                'label': 'Partner project marker collection',
                'identity': project_match.get('nearest_identity'),
                'hit': project_match.get('nearest_hit'),
                'query_coverage': project_match.get('query_coverage'),
                'alignment_length': project_match.get('alignment_length'),
            },
            'project_identity_warning': qc.strip().upper() != 'PASS_HIGH_CONFIDENCE',
            'project_identity_warning_reason': (
                f'Marker QC is {qc}; interpret this project-marker identity cautiously.'
                if qc.strip().upper() != 'PASS_HIGH_CONFIDENCE' else ''
            ),
            'project_close_neighbours': _dashboard_number(source.get('ProjectCloseNeighboursGE97', evidence.get('ProjectMatchesGE97', ''))),
            'pangenome_rank': pangenome_by_candidate.get(sequence_id, {}).get('rank') or str(source.get('PangenomeRank', evidence.get('PangenomeRank', '')) or ''),
            'pangenome_taxon': pangenome_by_candidate.get(sequence_id, {}).get('taxon') or str(source.get('PangenomeTaxon', evidence.get('PangenomeTaxon', '')) or ''),
            'pangenome_gap': pangenome_by_candidate.get(sequence_id, {}).get('gap') if sequence_id in pangenome_by_candidate else _dashboard_number(source.get('PangenomeGap', evidence.get('PangenomeGap', ''))),
            'pangenome_target': pangenome_by_candidate.get(sequence_id, {}).get('target') if sequence_id in pangenome_by_candidate else _dashboard_number(source.get('PangenomeTarget', evidence.get('PangenomeTarget', ''))),
            'recommendation_reason': str(source.get('RecommendationReason', '') or ''),
            'species_context': str(source.get('SpeciesContext', '') or ''),
            'genome_coverage': str(source.get('GenomeCoverage', '') or (
                'Genome already sequenced' if str(source.get('GenomeAlreadySequenced', '')).strip().lower() in {'true', '1', 'yes'}
                else 'Selected for sequencing; genome pending' if str(source.get('AlreadySelectedForGenomeSequencing', '')).strip().lower() in {'true', '1', 'yes'}
                else 'Not currently selected for sequencing'
            )),
            'local_tree_figure': str(source.get('LocalTreeFigure', evidence.get('LocalNeighbourhoodFigure', '')) or ''),
            'taxonomy': taxonomy,
            'novelty': novelty_comparisons,
            'genome_nearest_hit': str(novelty_source.get('NearestGenomeHit', '') or ''),
            'genome_nearest_identity': _dashboard_number(novelty_source.get('NearestGenomeIdentity', '')),
            'genome_matches_97': _dashboard_number(novelty_source.get('GenomeCollectionMatchesGE97', '')),
            'genome_matches_95': _dashboard_number(novelty_source.get('GenomeCollectionMatchesGE95', '')),
            'genome_representation_interpretation': str(novelty_source.get('GenomeRepresentationInterpretation', '') or ''),
        })

    # An optional live database overlay lets the user update a historical
    # dashboard after the SAB meeting without rerunning biological analyses.
    if database_states:
        decision_labels = {
            'SAB_APPROVED': 'APPROVED', 'SAB_DECLINED': 'REJECTED',
            'SAB_DEFERRED': 'DEFERRED',
        }
        for record in records:
            state = database_states.get(record['id'])
            if state is None:
                continue
            selected = bool(state.get('selected'))
            sequenced = bool(state.get('sequenced'))
            record['already_selected_for_genome_sequencing'] = selected
            record['genome_already_sequenced'] = sequenced
            record['sab_decision'] = decision_labels.get(state.get('status', ''), record.get('sab_decision', ''))
            if sequenced:
                record['recommendation'] = 'ALREADY SEQUENCED'
                record['sequencing_role'] = 'SEQUENCED'
                record['genome_coverage'] = 'Genome already sequenced'
            elif selected:
                record['recommendation'] = 'ALREADY SELECTED - GENOME PENDING'
                record['sequencing_role'] = 'COMMITTED'
                record['genome_coverage'] = 'Selected for sequencing; genome pending'

    counts = Counter(row['recommendation'] for row in records)
    qc_counts = Counter(row['marker_qc'] for row in records)
    partner_counts = Counter(row['partner'] for row in records)
    evidence_counts = Counter(row['evidence_quality'] or 'Not available' for row in records)
    taxonomy_counts = {
        db_label: Counter(
            (next((item['status'] for item in row['taxonomy'] if item['database'] == db_label), '') or 'Not available').upper()
            for row in records
        ) for db_label, _ in databases
    }

    manifest = {}
    if manifest_path.is_file():
        try:
            manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
        except (OSError, json.JSONDecodeError):
            manifest = {}
    run_info = {
        'workflow': manifest.get('workflow', 'performance_review'),
        'status': manifest.get('status', 'Not recorded'),
        'started_at': manifest.get('started_at', ''),
        'completed_at': manifest.get('completed_at', ''),
        'version': manifest.get('branchmanager_version', ''),
        'workflow_version': manifest.get('workflow_version', 'Not recorded'),
        'database_schema_version': (manifest.get('database_identity') or {}).get('schema_version', 'Not recorded'),
        'database_sha256': (manifest.get('database_identity') or {}).get('sha256', ''),
        'git_revision': manifest.get('git_revision', ''),
        'database_overlay': bool(database),
        'database_path': database_identity.get('path', ''),
    }
    if database:
        run_info['database_schema_version'] = database_identity.get('schema_version', 'unknown')
        run_info['database_sha256'] = database_identity.get('sha256', '')

    # Aggregate the hierarchy at rank + taxon level. Candidate IDs are
    # deduplicated across repeated local-tree/context rows; when an older
    # hierarchy export has no IDs, retain the largest row count rather than
    # summing repeated group-level counts. Committed counts are group-level
    # values and are therefore also reconciled by maximum, never summed.
    pangenome_groups = {}
    for row in pangenome_rows:
        rank = str(row.get('PangenomeRank', '') or '').strip().upper()
        taxon = str(row.get('PangenomeTaxon', '') or '').strip()
        if not rank or not taxon:
            continue
        key = (rank, taxon.casefold())
        item = pangenome_groups.setdefault(key, {
            'rank': rank, 'taxon': taxon,
            'provenance_values': [],
            'candidate_ids': set(), 'committed_ids': set(),
            'candidate_count_fallback': 0,
            'committed': None, 'target': None,
            'eligible_values': [], 'reasons': [], 'source_rows': 0,
        })
        item['source_rows'] += 1
        provenance = str(row.get('TaxonomyProvenance', '') or '').strip()
        if provenance and provenance not in item['provenance_values']:
            item['provenance_values'].append(provenance)

        raw_candidates = str(row.get('CandidateIDs', '') or '').strip()
        if raw_candidates and raw_candidates not in {'None', 'NA'}:
            item['candidate_ids'].update(
                value.strip() for value in raw_candidates.split(';') if value.strip()
            )
        else:
            count = _dashboard_number(row.get('CandidateCount', ''))
            if count is not None:
                item['candidate_count_fallback'] += int(count)

        raw_committed = str(row.get('CommittedGenomeIDs', '') or '').strip()
        if raw_committed and raw_committed not in {'None', 'NA'}:
            item['committed_ids'].update(
                value.strip() for value in raw_committed.split(';') if value.strip()
            )

        committed = _dashboard_number(row.get('CommittedGenomeCount', ''))
        if committed is not None:
            item['committed'] = (
                int(committed) if item['committed'] is None
                else max(item['committed'], int(committed))
            )
        target = _dashboard_number(row.get('PangenomeTarget', ''))
        if target is not None:
            item['target'] = (
                int(target) if item['target'] is None
                else max(item['target'], int(target))
            )
        eligible = str(row.get('PangenomeEligible', '') or '').strip()
        if eligible:
            item['eligible_values'].append(eligible)
        reason = str(
            row.get('PangenomeEligibilityReason', '')
            or row.get('EligibilityReason', '') or ''
        ).strip()
        if reason and reason not in item['reasons']:
            item['reasons'].append(reason)

    # The selection summary is the candidate-level source for validated
    # partner-catalogue/tree eligibility. It can also supply member IDs and
    # the group-level committed count when the hierarchy export lacks them.
    selection_pangenome = {}
    for row in selection_rows:
        rank = str(row.get('PangenomeRank', '') or '').strip().upper()
        taxon = str(row.get('PangenomeTaxon', '') or '').strip()
        if not rank or not taxon:
            continue
        key = (rank, taxon.casefold())
        group = selection_pangenome.setdefault(key, {
            'candidate_ids': set(), 'eligible_ids': set(),
            'target': None, 'committed': None, 'reasons': [],
        })
        sid = str(row.get('SequenceID', '') or '').strip()
        if sid:
            group['candidate_ids'].add(sid)
        if str(row.get('PangenomeEligible', '') or '').strip().lower() in {'yes', 'true', '1'}:
            if sid:
                group['eligible_ids'].add(sid)
            target = _dashboard_number(row.get('PangenomeTarget', ''))
            committed = _dashboard_number(row.get('PangenomeCommittedGenomeCount', ''))
            if target is not None:
                group['target'] = target if group['target'] is None else max(group['target'], target)
            if committed is not None:
                group['committed'] = committed if group['committed'] is None else max(group['committed'], committed)
            reason = str(row.get('PangenomeEligibilityReason', '') or '').strip()
            if reason and reason not in group['reasons']:
                group['reasons'].append(reason)

    for key, selected in selection_pangenome.items():
        if not selected['eligible_ids']:
            continue
        item = pangenome_groups.setdefault(key, {
            'rank': key[0],
            'taxon': next((str(r.get('PangenomeTaxon', '')).strip() for r in selection_rows
                           if str(r.get('PangenomeRank', '')).strip().upper() == key[0]
                           and str(r.get('PangenomeTaxon', '')).strip().casefold() == key[1]), key[1]),
            'provenance_values': [], 'candidate_ids': set(), 'committed_ids': set(),
            'candidate_count_fallback': 0, 'committed': None, 'target': None,
            'eligible_values': [], 'reasons': [], 'source_rows': 0,
        })
        item['candidate_ids'].update(selected['candidate_ids'])
        item['eligible_values'].append('Yes')
        if selected['target'] is not None:
            item['target'] = selected['target'] if item['target'] is None else max(item['target'], selected['target'])
        if selected['committed'] is not None:
            item['committed'] = selected['committed'] if item['committed'] is None else max(item['committed'], selected['committed'])
        for reason in selected['reasons']:
            if reason not in item['reasons']:
                item['reasons'].append(reason)

    pangenome = []
    for item in pangenome_groups.values():
        eligible_values = {str(value).strip().lower() for value in item.pop('eligible_values')}
        item['eligible'] = (
            'Yes' if 'yes' in eligible_values
            else 'No' if 'no' in eligible_values
            else 'Not recorded'
        )
        item['reason'] = '; '.join(item.pop('reasons'))
        provenance_values = item.pop('provenance_values')
        item['provenance'] = (
            provenance_values[0] if len(provenance_values) == 1
            else 'MULTIPLE_SOURCE_ROWS' if provenance_values
            else 'Taxonomy source not recorded'
        )
        item['candidate_ids'] = sorted(item['candidate_ids'])
        item['committed_ids'] = sorted(item['committed_ids'])
        item['candidates'] = max(
            len(item['candidate_ids']), int(item.pop('candidate_count_fallback', 0))
        )
        item['committed'] = int(item['committed'] or 0)
        # A supplied gap can be stale or inconsistent with the committed
        # count. Derive it from the reconciled target and committed count.
        item['gap'] = (
            max(0, int(item['target']) - item['committed'])
            if item['eligible'] == 'Yes' and item['target'] is not None
            else None
        )
        pangenome.append(item)
    pangenome.sort(key=lambda x: (
        ['DOMAIN','PHYLUM','CLASS','ORDER','FAMILY','GENUS','SPECIES'].index(x['rank'])
        if x['rank'] in ['DOMAIN','PHYLUM','CLASS','ORDER','FAMILY','GENUS','SPECIES'] else 99,
        x['taxon'].casefold()
    ))
    pangenome_source_rows = len(pangenome_rows)

    genome_path = Path(genome_assessment) if genome_assessment else selection_path.with_name('genome_collection_assessment.tsv')
    genome_rows = database_genome_rows if database else _read_tsv_rows(genome_path)
    genome_records = []
    for row in genome_rows:
        genome_records.append({
            'id': str(row.get('sequence_id', '') or ''),
            'partner': str(row.get('partner_id', '') or ''),
            'genome_id': str(row.get('genome_id', '') or ''),
            'status': str(row.get('genome_status', '') or ''),
            'qc': str(row.get('genome_qc_pass', '') or ''),
            'completeness': _dashboard_number(row.get('completeness')),
            'contamination': _dashboard_number(row.get('contamination')),
            'taxonomy': _dashboard_taxon(str(row.get('gtdb_taxonomy', '') or '')),
            'lineage': str(row.get('gtdb_taxonomy', '') or ''),
            'path': str(row.get('genome_path', '') or ''),
            'notes': str(row.get('notes', '') or ''),
        })
    record_by_id = {row['id']: row for row in records}
    genome_by_id = {row['id']: row for row in genome_records if row.get('id')}
    for group in pangenome:
        member_ids = sorted(set(group.get('candidate_ids', [])) | set(group.get('committed_ids', [])))
        members = []
        for member_id in member_ids:
            isolate = record_by_id.get(member_id)
            genome = genome_by_id.get(member_id)
            if isolate:
                kind = 'Assessed marker isolate'
                partner = isolate.get('partner') or 'Not recorded'
                recommendation = isolate.get('recommendation') or ''
                genome_id = genome.get('genome_id', '') if genome else ''
                qc = isolate.get('marker_qc', '')
                taxonomy = isolate.get('gtdb_taxonomy', '')
            elif genome:
                kind = 'Genome collection record'
                partner = genome.get('partner') or 'Not recorded'
                recommendation = ''
                genome_id = genome.get('genome_id', '')
                qc = 'Genome QC passed' if str(genome.get('qc', '')).strip().lower() in {'1','true','yes','pass'} else 'Genome QC not passed / not recorded'
                taxonomy = genome.get('lineage', '')
            else:
                kind = 'Registered sequence / baseline (dataset not available in dashboard inputs)'
                partner = ''
                recommendation = ''
                genome_id = ''
                qc = ''
                taxonomy = ''
            members.append({
                'id': member_id, 'kind': kind, 'partner': partner,
                'recommendation': recommendation, 'genome_id': genome_id,
                'qc': qc, 'taxonomy': taxonomy,
                'candidate': member_id in set(group.get('candidate_ids', [])),
                'committed': member_id in set(group.get('committed_ids', [])),
            })
        group['members'] = members
        group['member_count'] = len(members)

    # The dashboard lives in the report directory; resolve the manifest from
    # the assessment directory beside selection_summary, not from the process
    # working directory.  This also supports callers passing relative paths.
    tree_manifest = selection_path.parent / 'neighbourhoods' / 'neighbourhood_manifest.tsv'
    if not tree_manifest.is_file():
        tree_manifest = output_dir / 'assessment' / 'neighbourhoods' / 'neighbourhood_manifest.tsv'
    tree_rows = []
    seen_tree_figures = set()
    for row in _read_tsv_rows(tree_manifest):
        figure = str(row.get('Figure') or '').strip()
        if not figure or figure in seen_tree_figures or str(row.get('Status') or '').lower() != 'rendered':
            continue
        image_path = _resolve_visual_file(tree_manifest, figure)
        if image_path and image_path.is_file():
            seen_tree_figures.add(figure)
            tree_rows.append({
                'clade': str(row.get('TreeLeafID') or row.get('SequenceID') or ''),
                'label': str(row.get('SequenceID') or 'Clade neighbourhood'),
                'figure': _html_path_target(image_path, output_dir / 'performance_review_dashboard.html'),
                'href': _html_path_target(image_path, output_dir / 'performance_review_dashboard.html'),
                'members': str(row.get('AssessedSequencesInFigure') or ''),
            })

    data = {
        'records': records,
        'genome_records': genome_records,
        'tree_rows': tree_rows,
        'pangenome': pangenome,
        'pangenome_source_rows': pangenome_source_rows,
        'novelty_source_rows': len(novelty_rows),
        'counts': {
            'total': len(records),
            'recommendations': dict(counts),
            'qc': dict(qc_counts),
            'partners': dict(partner_counts),
            'evidence': dict(evidence_counts),
            'taxonomy': {key: dict(value) for key, value in taxonomy_counts.items()},
            'already_sequenced': int(counts.get('ALREADY SEQUENCED', 0)),
            'selected_pending': int(counts.get('ALREADY SELECTED - GENOME PENDING', 0)),
            'genome_only_records': len(genome_records),
            'genome_qc_pass': sum(
                str(row.get('genome_qc_pass', '')).strip().casefold() in {'1', 'true', 'yes', 'pass'}
                for row in genome_rows
            ),
        },
        'run': run_info,
        'sources': {
            'selection': 'assessment/selection_summary.tsv',
            'assessment': 'assessment/sequence_assessment.tsv',
            'pangenome': 'assessment/pangenome_hierarchy.tsv',
            'novelty': 'assessment/novelty_metrics.tsv',
            'taxonomy': 'taxonomy/all_databases.tsv',
        },
    }
    payload = json.dumps(data, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    output = output_dir / 'performance_review_dashboard.html'
    html_doc = r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>BranchManager | Genome assessment</title>
<style>
:root{--ink:#17251f;--muted:#617168;--line:#dce5df;--paper:#fff;--wash:#f4f7f5;--green:#176b4a;--green2:#e7f2ec;--amber:#8b5b0a;--amber2:#fff4dc;--red:#a43b35;--red2:#fbeceb;--blue:#315e85;--blue2:#eaf1f7;--shadow:0 8px 28px #20392c0c}
*{box-sizing:border-box}html{scroll-behavior:smooth;scroll-padding-top:20px}body{margin:0;background:var(--wash);color:var(--ink);font:14px/1.5 system-ui,-apple-system,Segoe UI,Roboto,sans-serif}a{color:var(--green)}.shell{max-width:none;margin:0;padding:28px 34px 56px 286px}.top{display:flex;justify-content:space-between;gap:20px;align-items:flex-start;margin-bottom:24px}.brand{font-size:12px;letter-spacing:.12em;text-transform:uppercase;font-weight:750;color:var(--green)}h1{font-size:clamp(26px,4vw,38px);line-height:1.12;margin:7px 0}.sub{color:var(--muted)}.status{border:1px solid var(--line);background:white;padding:8px 12px;border-radius:999px;font-size:12px;white-space:nowrap}.nav{position:fixed;left:0;top:0;bottom:0;width:250px;padding:26px 16px;background:#fff;border-right:1px solid var(--line);display:flex;flex-direction:column;gap:6px;margin:0;z-index:4;overflow:auto}.nav .navbrand{font-weight:800;color:var(--green);font-size:12px;letter-spacing:.08em;text-transform:uppercase;padding:0 10px 18px}.nav a{padding:10px 12px;border:1px solid transparent;border-radius:9px;background:transparent;text-decoration:none;font-size:13px;font-weight:650;color:var(--muted)}.nav a:hover,.nav a.active{background:var(--green2);border-color:#d0e5d8;color:var(--green)}.nav a{padding:8px 12px;border:1px solid var(--line);border-radius:999px;background:white;text-decoration:none;font-size:12px;font-weight:650}.grid{display:grid;gap:14px}.metrics{grid-template-columns:repeat(4,minmax(0,1fr))}.card,.panel{background:var(--paper);border:1px solid var(--line);border-radius:14px;box-shadow:var(--shadow)}.metric{padding:18px}.metric strong{font-size:30px;display:block;line-height:1.1}.metric span{color:var(--muted);font-size:12px}.section{margin-top:26px;scroll-margin-top:16px}.section-head{display:flex;justify-content:space-between;align-items:end;gap:12px;margin-bottom:12px}h2{font-size:20px;margin:0}h3{font-size:15px;margin:0 0 10px}.hint{color:var(--muted);font-size:12px}.two{grid-template-columns:1fr 1fr}.panel{padding:18px}.barlist{display:grid;gap:11px}.barrow{display:grid;grid-template-columns:minmax(130px,220px) 1fr 42px;gap:10px;align-items:center;font-size:12px}.track{height:10px;border-radius:99px;background:#edf1ee;overflow:hidden}.fill{height:100%;border-radius:99px;background:var(--green)}.fill.amber{background:#c18a2e}.fill.red{background:#b85b54}.fill.blue{background:#6486a5}.legend{display:flex;flex-wrap:wrap;gap:12px;color:var(--muted);font-size:11px;margin-top:12px}.dot{display:inline-block;width:8px;height:8px;border-radius:50%;margin-right:5px;background:var(--green)}.callout{background:var(--green2);border:1px solid #d0e5d8;padding:14px 16px;border-radius:12px;color:#254d38}.callout p{margin:4px 0 0}.tablewrap{overflow:auto;border:1px solid var(--line);border-radius:12px;background:white}table{border-collapse:collapse;width:100%;min-width:850px}th,td{text-align:left;padding:10px 11px;border-bottom:1px solid #edf1ee;vertical-align:top;font-size:12px}th{position:sticky;top:0;background:#f2f6f3;color:#46564d;z-index:1}tr:hover td{background:#fafcfb}.filters{display:flex;flex-wrap:wrap;gap:8px;margin:12px 0}.filters input,.filters select{padding:9px 11px;border:1px solid var(--line);border-radius:9px;background:white;color:var(--ink);font:inherit;min-width:150px}.filters input{flex:1;min-width:220px}.btn{border:1px solid var(--line);background:white;border-radius:9px;padding:8px 11px;font:inherit;cursor:pointer}.btn.primary{background:var(--green);border-color:var(--green);color:white}.pill{display:inline-block;padding:3px 7px;border-radius:999px;background:#eef2ef;font-size:10px;font-weight:700;white-space:nowrap}.pill.good{background:var(--green2);color:var(--green)}.pill.warn{background:var(--amber2);color:var(--amber)}.pill.bad{background:var(--red2);color:var(--red)}.pill.info{background:var(--blue2);color:var(--blue)}.detail{position:fixed;inset:0;background:#10221966;display:none;z-index:5;padding:20px;overflow:auto}.detail.open{display:block}.drawer{max-width:950px;margin:20px auto;background:white;border-radius:16px;padding:24px;box-shadow:0 20px 80px #0003}.drawerhead{display:flex;justify-content:space-between;gap:16px;align-items:start}.detailgrid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:12px;margin:16px 0}.qc-warning{display:inline-block;margin-left:5px;padding:2px 5px;border-radius:5px;background:#fff1d6;color:#8a4b00;font-size:10px;font-weight:700;white-space:nowrap}.qc-warning-text{margin:8px 0 0;color:#8a4b00;font-size:12px;line-height:1.45}.detailbox{border:1px solid var(--line);border-radius:10px;padding:12px}.detailbox small{color:var(--muted);display:block}.taxgrid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:10px}.taxcard{border:1px solid var(--line);border-radius:10px;padding:12px}.taxcard p{margin:5px 0;font-size:12px}.taxcard .lineage{color:var(--muted);overflow-wrap:anywhere}.empty{padding:24px;text-align:center;color:var(--muted)}.to-top{position:fixed;right:22px;bottom:22px;z-index:3;border:0;border-radius:999px;background:var(--green);color:white;padding:11px 15px;font-weight:700;box-shadow:0 5px 18px #0002;cursor:pointer;display:none}.to-top.visible{display:block}.foot{margin-top:28px;color:var(--muted);font-size:11px;border-top:1px solid var(--line);padding-top:16px}.smalllinks{display:flex;flex-wrap:wrap;gap:12px}.smalllinks a{font-size:12px}@media(max-width:850px){.shell{padding:82px 18px 40px}.nav{right:0;bottom:auto;width:auto;height:64px;top:0;padding:8px 12px;flex-direction:row;align-items:center;overflow-x:auto;border-right:0;border-bottom:1px solid var(--line)}.nav .navbrand{display:none}.nav a{white-space:nowrap;padding:8px 10px}.metrics{grid-template-columns:repeat(2,minmax(0,1fr))}.two{grid-template-columns:1fr}.top{display:block}.status{display:inline-block;margin-top:12px}}@media(max-width:560px){.shell{padding:78px 12px 40px}.metrics{grid-template-columns:1fr 1fr}.metric{padding:13px}.metric strong{font-size:25px}.detailgrid,.taxgrid{grid-template-columns:1fr}.barrow{grid-template-columns:105px 1fr 30px}}
.section:not(#overview){height:min(76vh,820px);min-height:420px;overflow:auto;overscroll-behavior:contain;padding:16px;background:var(--paper);border:1px solid var(--line);border-radius:14px;box-shadow:var(--shadow)}.section:not(#overview)>.section-head{position:sticky;top:-16px;z-index:2;background:var(--paper);padding:10px 0 12px;margin:-6px 0 12px;border-bottom:1px solid var(--line)}.section:not(#overview)>.section-head h2{line-height:1.25}.section:not(#overview)>.callout{margin-top:14px}.section:not(#overview) .tablewrap{max-height:calc(76vh - 190px);overflow:auto}.section#phylogeny .gallery{display:grid;grid-template-columns:repeat(auto-fill,minmax(280px,1fr));gap:12px}.section#phylogeny .tree-card img{height:210px}.section#genomes .tablewrap{max-height:calc(76vh - 180px)}@media(max-width:850px){.section:not(#overview){height:72vh;min-height:360px;padding:12px}.section:not(#overview)>.section-head{top:-12px}.section:not(#overview) .tablewrap{max-height:calc(72vh - 180px)}}
.pangenome-row{cursor:pointer}.pangenome-row:hover td{background:#edf6f0!important}.member-table{min-width:720px}.tree-card{display:flex;flex-direction:column;gap:8px;padding:10px;text-decoration:none;min-width:0}.tree-card img{width:100%;height:190px;object-fit:contain;background:#f8faf8;border:1px solid var(--line);border-radius:8px}.tree-card strong{color:var(--ink)}.tree-card .hint{color:var(--muted)}/* Dashboard layout refinements: compact narrative panels; scroll only data-heavy regions. */
.section{min-height:0!important;height:auto!important;max-height:none!important;overflow:visible!important;padding:18px 20px!important;margin-bottom:16px!important}
.callout{padding:10px 14px!important;margin:10px 0!important}.callout p{margin:5px 0!important}#overviewText{min-height:0!important}.tablewrap{max-height:440px;overflow:auto}#treeGallery{max-height:560px;overflow:auto;align-content:start}#coverage .tablewrap{max-height:520px}#explorer .tablewrap,#novelty .tablewrap{max-height:520px}#genomes .tablewrap{max-height:560px}#phylogeny .gallery{max-height:560px;overflow:auto}
</style></head><body><main class="shell">
<header class="top"><div><div class="brand">BranchManager · Rumen Gateway</div><h1>Isolate assessment overview</h1><div class="sub" id="runSubtitle">Performance review · Interactive report</div></div><div class="status" id="runStatus">Assessment status</div></header>
<nav class="nav" aria-label="Report sections"><div class="navbrand">BranchManager<br><span style="font-weight:500;color:var(--muted);letter-spacing:0">Assessment report</span></div><a href="#overview" class="active">Overview</a><a href="#quality">Marker quality</a><a href="#taxonomy">Taxonomic evidence</a><a href="#novelty">Novelty explorer</a><a href="#coverage">Pangenome coverage</a><a href="#genomes">Genome collection</a><a href="#phylogeny">Phylogeny</a><a href="#explorer">Isolate explorer</a><a href="#sources">Methods & files</a></nav>
<section id="overview"><div class="grid metrics" id="metrics"></div><div class="grid two section"><div class="panel"><h3>Sequencing recommendations</h3><div id="recommendationChart" class="barlist"></div><p class="hint">These are selection recommendations, not marker quality classifications.</p></div><div class="panel"><h3>Sequences by partner</h3><div id="partnerChart" class="barlist"></div><p class="hint">Counts are based on records in the selection summary.</p></div></div><div class="callout section" id="overviewText"></div></section>
<section class="section" id="quality"><div class="section-head"><div><h2>Marker quality</h2><div class="hint">Quality-control status is separate from the genome sequencing recommendation.</div></div></div><div class="grid two"><div class="panel"><h3>Marker QC classifications</h3><div id="qcChart" class="barlist"></div></div><div class="panel"><h3>Evidence quality</h3><div id="evidenceChart" class="barlist"></div><p class="hint">Evidence quality is the assessment's own category; it is not a probability of correct identification.</p></div></div><div class="callout section"><strong>How to read this:</strong><p>“Pass with warnings” means the marker was retained with qualifications. “Fail QC” means it did not meet the normal marker quality requirements. Any manual review or override should be examined in the underlying evidence before a sequencing decision is made.</p></div></section>
<section class="section" id="taxonomy"><div class="section-head"><div><h2>Taxonomic evidence</h2><div class="hint">Assignments from each reference database are shown independently.</div></div></div><div class="panel"><div id="taxonomyChart" class="barlist"></div><div class="legend"><span><i class="dot"></i>Confident</span><span><i class="dot" style="background:#c18a2e"></i>Contested / provisional</span><span><i class="dot" style="background:#b85b54"></i>Low similarity / unresolved</span></div></div><div class="callout section"><strong>Important interpretation:</strong><p>A high sequence identity does not automatically guarantee an unambiguous taxonomic assignment. Query coverage describes how much of the submitted sequence aligns; reference coverage describes how much of the reference sequence is covered. Conflicting near-best reference assignments can result in a contested classification even when the alignment is long and highly similar.</p></div></section>

<section class="section" id="novelty"><div class="section-head"><div><h2>Novelty explorer</h2><div class="hint">Compare each marker against the cultured baselines, wider partner collection and external reference set. These are sequence-similarity comparisons, not species novelty claims.</div></div><button class="btn" id="noveltyExport">Export novelty results</button></div><div class="grid metrics" id="noveltyMetrics"></div><div class="callout"><strong>How to interpret these comparisons:</strong><p>Nearest identity is the highest sequence identity found in that comparison set. The ≥99%, ≥97% and ≥95% counts describe the number of reported neighbours meeting those operational similarity cut-offs. They are not formal species boundaries. A low identity to a cultured baseline can coexist with a high identity to an external reference; that means the isolate may be poorly represented in the cultured collection while having a close external reference match.</p></div><div class="panel"><div class="filters"><input id="noveltySearch" placeholder="Search isolate, partner, closest hit or taxonomy"><select id="noveltyPartner"><option value="">All partners</option></select><select id="noveltyComparison"><option value="CulturedProject">Combined cultured collection</option><option value="Hungate">Hungate baseline</option><option value="SecondaryBaseline">Secondary cultured baseline</option><option value="Project">Partner project markers</option><option value="Reference">GTDB reference genomes</option></select><select id="noveltySort"><option value="identity_asc">Lowest nearest identity first</option><option value="identity_desc">Highest nearest identity first</option><option value="matches_desc">Most ≥97% neighbours first</option><option value="score_desc">Highest novelty score first</option><option value="id_asc">Isolate ID</option></select><select id="noveltyPageSize"><option>25</option><option selected>50</option><option>100</option><option>250</option></select></div><div class="hint" id="noveltyCount"></div><div class="tablewrap"><table><thead><tr><th>Isolate</th><th>Partner</th><th>Nearest hit</th><th>Identity</th><th>≥99%</th><th>≥97%</th><th>≥95%</th><th>Crowding</th><th>Priority</th><th>Details</th></tr></thead><tbody id="noveltyBody"></tbody></table></div><div class="pager" id="noveltyPager"></div></div><div class="panel" id="noveltyDetail" style="display:none;margin-top:14px"><div class="section-head"><div><h3 id="noveltyDetailTitle">Isolate comparison</h3><div class="hint" id="noveltyDetailSubtitle"></div></div><button class="btn" id="closeNoveltyDetail">Close</button></div><div class="tablewrap"><table><thead><tr><th>Comparison set</th><th>Nearest hit</th><th>Identity</th><th>Query coverage</th><th>≥99%</th><th>≥97%</th><th>≥95%</th><th>Crowding</th><th>Priority</th></tr></thead><tbody id="noveltyDetailBody"></tbody></table></div><div class="hint" id="noveltyGenomeContext"></div></div></section>
<section class="section" id="genomes"><div class="section-head"><div><h2>Genome collection assessment</h2><div class="hint">Genome-only and marker-associated assemblies are assessed independently of marker QC.</div></div></div><div class="grid metrics" id="genomeMetrics"></div><div class="panel"><div class="tablewrap"><table><thead><tr><th>Isolate</th><th>Partner</th><th>Genome</th><th>Genome quality status</th><th>Completeness</th><th>Contamination</th><th>GTDB taxonomy</th></tr></thead><tbody id="genomeBody"></tbody></table></div></div></section>
<section class="section" id="coverage"><div class="section-head"><div><h2>Pangenome coverage</h2><div class="hint">Assessed isolates are unique marker records in the group. Genome commitments are taken from the assessment hierarchy; open a row to distinguish existing genome records from selected isolates awaiting sequencing. A blank target means the group is not eligible for a formal target.</div></div><a class="btn" href="assessment/pangenome_hierarchy.tsv">Full hierarchy TSV</a></div><div class="grid metrics" id="pangenomeMetrics"></div><div class="panel"><div class="filters"><input id="pangenomeSearch" placeholder="Search rank, taxon or eligibility reason"><select id="rankFilter"><option value="">All ranks</option></select></div><div class="hint" id="pangenomeCount"></div><div class="tablewrap"><table><thead><tr><th>Rank</th><th>Taxon</th><th>Assessed isolates</th><th>Genome commitments</th><th>Target</th><th>Remaining gap</th><th>Target assessable?</th><th>Why a target is or is not available</th></tr></thead><tbody id="pangenomeBody"></tbody></table><div class="pager" id="pangenomePager"></div></div></div></section>
<section class="section" id="phylogeny"><div class="section-head"><div><h2>Phylogeny and clade neighbourhoods</h2><div class="hint">Local clade figures generated by the assessment; select a figure to inspect it at full size.</div></div></div><div class="gallery" id="treeGallery"></div></section>
<section class="section" id="explorer"><div class="section-head"><div><h2>Isolate explorer</h2><div class="hint">Search and filter the assessed sequences, then select a row for its evidence summary.</div></div><button class="btn" id="exportBtn">Export filtered CSV</button></div><div class="filters"><input id="search" placeholder="Search isolate, taxonomy or recommendation"><select id="partnerFilter"><option value="">All partners</option></select><select id="recommendationFilter"><option value="">All recommendations</option></select><select id="qcFilter"><option value="">All marker QC</option></select><select id="pageSize"><option>25</option><option selected>50</option><option>100</option><option>250</option></select><button class="btn" id="clearFilters">Clear filters</button></div><div class="hint" id="resultCount"></div><div class="tablewrap"><table><thead><tr><th>Isolate</th><th>Partner</th><th>Marker QC</th><th>Recommendation</th><th>GTDB assignment</th><th>Closest cultured identity</th><th>Closest project identity</th></tr></thead><tbody id="isolateBody"></tbody></table></div><div class="pager" id="isolatePager"></div></section>
<section class="section" id="sources"><div class="section-head"><div><h2>Methods, provenance and source files</h2><div class="hint">The dashboard summarises existing assessment outputs; it does not rerun analyses.</div></div></div><div class="panel"><div id="runDetails" class="detailgrid"></div><div class="smalllinks"><a href="assessment/selection_summary.tsv">Selection summary TSV</a><a href="assessment/sequence_assessment.tsv">Full sequence assessment TSV</a><a href="assessment/novelty_metrics.tsv">Novelty metrics TSV</a><a href="assessment/pangenome_hierarchy.tsv">Pangenome hierarchy TSV</a><a href="assessment/genome_collection_assessment.tsv">Genome collection TSV</a><a href="taxonomy/all_databases.tsv">All-database taxonomy TSV</a><a href="run_manifest.json">Run manifest</a></div></div></section>
<footer class="foot">Generated by BranchManager. This report is decision support, not a substitute for reviewing the underlying evidence. Taxonomic and novelty statements are limited to the configured reference collections, thresholds and analyses used in this run.</footer>
</main><button class="to-top" id="toTop" aria-label="Back to top" title="Back to top">↑ <span>Top</span></button><div class="detail" id="detail" role="dialog" aria-modal="true"><div class="drawer"><div class="drawerhead"><div><div class="brand">Individual isolate evidence</div><h2 id="detailTitle"></h2><div id="detailSubtitle" class="hint"></div></div><button class="btn" id="closeDetail">Close ✕</button></div><div id="detailContent"></div></div></div><div class="detail" id="pangenomeDetail" role="dialog" aria-modal="true"><div class="drawer"><div class="drawerhead"><div><div class="brand">Pangenome group membership</div><h2 id="pangenomeDetailTitle"></h2><div id="pangenomeDetailSubtitle" class="hint"></div></div><button class="btn" id="closePangenomeDetail">Close ✕</button></div><div id="pangenomeDetailContent"></div></div></div>
<script id="dashboard-data" type="application/json">__DATA__</script>
<script>
(()=>{'use strict';const D=JSON.parse(document.getElementById('dashboard-data').textContent),$=id=>document.getElementById(id),fmt=n=>n==null?'—':Number.isInteger(n)?String(n):n.toFixed(1),pct=n=>n==null?'—':fmt(n)+'%',esc=s=>String(s??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c])),pill=(s)=>{let v=String(s||'').toUpperCase(),c=v.includes('FAIL')||v.includes('LOW_SIMILARITY')?'bad':v.includes('CONTEST')||v.includes('PROVISION')||v.includes('WARN')||v.includes('REVIEW')?'warn':v.includes('CONFIDENT')||v.includes('PASS')?'good':'info';return `<span class="pill ${c}">${esc(s||'Not available')}</span>`};
function pager(id,page,pages,total,onPage){
 const el=$(id); if(!el)return;
 if(!total){el.innerHTML='';return;}
 const button=(label,target,disabled=false,current=false)=>`<button type="button" class="btn pager-btn${current?' active':''}" data-page="${target}" ${disabled?'disabled':''}>${label}</button>`;
 let html=`<span class="hint">${total.toLocaleString()} records · Page ${page} of ${pages}</span>`;
 html+=button('Previous',page-1,page<=1);
 const start=Math.max(1,page-2),end=Math.min(pages,page+2);
 for(let p=start;p<=end;p++)html+=button(String(p),p,false,p===page);
 html+=button('Next',page+1,page>=pages);
 el.innerHTML=html;
 el.querySelectorAll('button[data-page]').forEach(b=>b.addEventListener('click',()=>onPage(Number(b.dataset.page))));
}

const genomeRows=D.genome_records||[];const qcPass=x=>['yes','1','true','pass'].includes(String(x.qc??'').trim().toLowerCase());const genomeQCLabel=x=>{const v=String(x.qc??'').trim().toLowerCase();if(['yes','1','true','pass'].includes(v))return 'Pass';if(['no','0','false','fail'].includes(v))return 'Does not pass';return x.status||'Not assessed'};const gqc=genomeRows.filter(qcPass).length;const gcomp=genomeRows.filter(x=>x.completeness!=null);$('genomeMetrics').innerHTML=[['Genome records',genomeRows.length],['QC passed',gqc],['Mean completeness',gcomp.length?fmt(gcomp.reduce((a,x)=>a+x.completeness,0)/gcomp.length)+'%':'—'],['Genome-only partners',[...new Set(genomeRows.map(x=>x.partner))].join(', ')||'—']].map(x=>`<div class="metric"><div class="metric-value">${esc(x[1])}</div><div class="metric-label">${esc(x[0])}</div></div>`).join('');$('genomeBody').innerHTML=genomeRows.map(x=>`<tr><td>${esc(x.id)}</td><td>${esc(x.partner)}</td><td>${esc(x.genome_id)}</td><td>${pill(genomeQCLabel(x))}</td><td>${pct(x.completeness)}</td><td>${pct(x.contamination)}</td><td title="${esc(x.lineage)}">${esc(x.taxonomy)}</td></tr>`).join('');const trees=D.tree_rows||[];$('treeGallery').innerHTML=trees.length?trees.map(x=>`<a class="panel tree-card" href="${esc(x.href)}" target="_blank"><img loading="lazy" src="${esc(x.figure)}" alt="${esc(x.label)}"><strong>${esc(x.label)}</strong><span class="hint">${esc(x.clade)} ${esc(x.members)}</span></a>`).join(''):'<p class="muted">No clade/neighbourhood figures were found.</p>';const total=D.counts.total; $('runSubtitle').textContent=`Performance review · ${total.toLocaleString()} assessed marker sequences${D.run.completed_at?' · completed '+D.run.completed_at.replace('T',' ').replace('Z',' UTC'):''}`;$('runStatus').textContent=`${D.run.status||'Status not recorded'}${D.run.version?' · v'+D.run.version:''}`;
const metrics=[['Marker isolates assessed',total,'Partner marker sequences in this assessment'],['Already sequenced',D.counts.already_sequenced||0,'Marker isolates with an existing genome'],['Selected · genome pending',D.counts.selected_pending||0,'Selected for sequencing; genome not yet recorded'],['Genome-only partner isolates',D.counts.genome_only_records||0,'Registered genome records without a marker-cohort entry'],['High-confidence marker QC',D.counts.qc.PASS_HIGH_CONFIDENCE||0,'Passed without significant QC warnings'],['Marker QC needs attention',(D.counts.qc.PASS_WITH_WARNINGS||0)+(D.counts.qc.FAIL_QC||0),'Warnings or QC failures'],['Partner collections',Object.keys(D.counts.partners).length,'Contributing marker collections'],['Genome QC passed',D.counts.genome_qc_pass||0,'Registered genomes passing the recorded QC flag']];$('metrics').innerHTML=metrics.map(x=>`<div class="card metric"><strong>${Number(x[1]).toLocaleString()}</strong><span>${esc(x[0])}</span><div class="hint">${esc(x[2])}</div></div>`).join('');
function bars(target,obj,colour=''){let entries=Object.entries(obj||{}).sort((a,b)=>b[1]-a[1]),max=Math.max(1,...entries.map(x=>x[1]));$(target).innerHTML=entries.map(([k,v])=>`<div class="barrow"><span title="${esc(k)}">${esc(k.replaceAll('_',' '))}</span><div class="track"><div class="fill ${colour}" style="width:${100*v/max}%"></div></div><b>${v}</b></div>`).join('')||'<div class="empty">No data available</div>'}
bars('recommendationChart',D.counts.recommendations);bars('partnerChart',D.counts.partners,'blue');bars('qcChart',D.counts.qc);bars('evidenceChart',D.counts.evidence,'blue');
$('overviewText').innerHTML=`<strong>What this assessment says</strong><p>${total.toLocaleString()} marker isolates from ${Object.keys(D.counts.partners).length} partner collections were assessed. ${Number(D.counts.already_sequenced||0).toLocaleString()} are already represented by genome data, and ${Number(D.counts.selected_pending||0).toLocaleString()} have been selected but are still awaiting a genome. A further ${Number(D.counts.genome_only_records||0).toLocaleString()} genome-only partner isolates are registered separately. ${Number(D.counts.qc.PASS_HIGH_CONFIDENCE||0).toLocaleString()} marker sequences passed QC without significant warnings. Marker quality and sequencing priority are different questions; use the isolate and pangenome sections together when considering new sequencing.</p>`;
let taxHtml='';for(const [db,counts] of Object.entries(D.counts.taxonomy)){let entries=Object.entries(counts).sort((a,b)=>b[1]-a[1]),sum=entries.reduce((a,x)=>a+x[1],0)||1;taxHtml+=`<div class="barrow"><span>${esc(db)}</span><div class="track" style="display:flex">${entries.map(([k,v])=>`<div title="${esc(k)}: ${v}" style="width:${100*v/sum}%;background:${k.includes('CONFIDENT')?'#27805a':k.includes('CONTEST')||k.includes('PROVISION')?'#c18a2e':k.includes('LOW')?'#b85b54':'#8b9a91'}"></div>`).join('')}</div><b>${sum}</b></div><div class="hint" style="margin:-5px 0 8px 0">${entries.map(([k,v])=>`${esc(k)}: ${v}`).join(' · ')}</div>`}$('taxonomyChart').innerHTML=taxHtml||'<div class="empty">No taxonomy evidence available</div>';
const ranks=[...new Set(D.pangenome.map(x=>x.rank).filter(Boolean))].sort((a,b)=>['DOMAIN','PHYLUM','CLASS','ORDER','FAMILY','GENUS','SPECIES'].indexOf(a.toUpperCase())-['DOMAIN','PHYLUM','CLASS','ORDER','FAMILY','GENUS','SPECIES'].indexOf(b.toUpperCase()));
$('rankFilter').innerHTML='<option value="">All ranks</option>'+ranks.map(x=>`<option>${esc(x)}</option>`).join('');
$('pangenomeMetrics').innerHTML=[['Taxon groups',D.pangenome.length],['Targets assessable',D.pangenome.filter(x=>x.eligible==='Yes').length],['Groups needing taxon/tree review',D.pangenome.filter(x=>x.eligible==='No').length],['Groups with a known gap',D.pangenome.filter(x=>x.gap!=null&&x.gap>0).length]].map(x=>`<div class="metric"><div class="metric-value">${Number(x[1]).toLocaleString()}</div><div class="metric-label">${esc(x[0])}</div></div>`).join('');
let pangenomeFiltered=[],pangenomePage=1;
function renderPangenome(){const q=$('pangenomeSearch').value.trim().toLowerCase(),r=$('rankFilter').value; pangenomeFiltered=D.pangenome.filter(x=>(!r||x.rank===r)&&(!q||[x.rank,x.taxon,x.reason].join(' ').toLowerCase().includes(q)));pangenomeFiltered.sort((a,b)=>(b.gap??-1)-(a.gap??-1)||a.taxon.localeCompare(b.taxon));const size=50,pages=Math.max(1,Math.ceil(pangenomeFiltered.length/size));pangenomePage=Math.min(pangenomePage,pages);const shown=pangenomeFiltered.slice((pangenomePage-1)*size,pangenomePage*size);$('pangenomeCount').textContent=`Showing ${shown.length} of ${pangenomeFiltered.length} taxon groups`; $('pangenomeBody').innerHTML=shown.map(x=>`<tr class="pangenome-row" data-pangenome-open="${D.pangenome.indexOf(x)}" tabindex="0" aria-label="View members of ${esc(x.rank)} ${esc(x.taxon)}"><td>${esc(x.rank)}</td><td><button class="btn" data-pangenome-open="${D.pangenome.indexOf(x)}">${esc(x.taxon)}</button><div class="hint">${esc(x.provenance||'Taxonomy source not recorded')}</div></td><td>${fmt(x.candidates)}</td><td>${fmt(x.committed)}</td><td>${fmt(x.target)}</td><td>${fmt(x.gap)}</td><td>${pill(x.eligible)}</td><td title="${esc(x.reason)}">${esc(x.reason)}</td></tr>`).join('')||'<tr><td colspan="8" class="empty">No pangenome groups match these filters</td></tr>';pager('pangenomePager',pangenomePage,pages,pangenomeFiltered.length,p=>{pangenomePage=p;renderPangenome()})}
['pangenomeSearch','rankFilter'].forEach(id=>$(id).addEventListener(id==='pangenomeSearch'?'input':'change',()=>{pangenomePage=1;renderPangenome()}));renderPangenome();
function openPangenomeGroup(g){if(!g)return;$('pangenomeDetailTitle').textContent=`${g.rank} · ${g.taxon}`;$('pangenomeDetailSubtitle').textContent=`${g.member_count||0} unique member(s) · ${g.eligible==='Yes'?'formal target assessable':'exploratory / not currently target-assessable'}`;const members=g.members||[];const rows=members.map(m=>`<tr><td><b>${esc(m.id)}</b><div class="hint">${esc(m.kind)}</div></td><td>${esc(m.partner||'—')}</td><td>${m.candidate?'Candidate':''}${m.candidate&&m.committed?' · ':''}${m.committed?'Genome commitment':''}</td><td>${esc(m.genome_id||'—')}</td><td>${esc(m.qc||'—')}</td><td>${esc(m.recommendation||'—')}</td><td>${esc(m.taxonomy||'—')}</td></tr>`).join('');$('pangenomeDetailContent').innerHTML=`<div class="detailgrid">${box('Target assessable',g.eligible)}${box('Assessed isolates',fmt(g.candidates))}${box('Genome commitments',fmt(g.committed))}${box('Target / remaining',fmt(g.target)+' / '+fmt(g.gap))}${box('Taxonomy provenance',g.provenance)}${box('Source hierarchy rows',fmt(g.source_rows))}</div><div class="detailbox"><small>Eligibility rationale</small><p>${esc(g.reason||'No rationale recorded.')}</p></div><div class="detailbox" style="margin-top:12px"><small>Isolates and genome records in this group</small><div class="tablewrap"><table class="member-table"><thead><tr><th>Record</th><th>Partner</th><th>Group role</th><th>Genome ID</th><th>QC</th><th>Recommendation</th><th>GTDB taxonomy</th></tr></thead><tbody>${rows||'<tr><td colspan="7" class="empty">No member identifiers were recorded in the hierarchy export.</td></tr>'}</tbody></table></div></div>`;$('pangenomeDetail').classList.add('open');document.body.style.overflow='hidden'}
$('pangenomeBody').addEventListener('click',e=>{const b=e.target.closest('[data-pangenome-open]');if(!b)return;const group=D.pangenome[Number(b.dataset.pangenomeOpen)];if(group)openPangenomeGroup(group)});
$('pangenomeBody').addEventListener('keydown',e=>{if((e.key==='Enter'||e.key===' ')&&e.target.matches('tr[data-pangenome-open]')){e.preventDefault();const group=D.pangenome[Number(e.target.dataset.pangenomeOpen)];if(group)openPangenomeGroup(group)}});
$('closePangenomeDetail').addEventListener('click',()=>{$('pangenomeDetail').classList.remove('open');document.body.style.overflow=''});
$('pangenomeDetail').addEventListener('click',e=>{if(e.target===$('pangenomeDetail'))$('closePangenomeDetail').click()});

$('rankFilter').innerHTML+ = '';
</script></body></html>'''
    # Complete the interactive script separately to keep the HTML template readable.
    html_doc = html_doc.replace("$('rankFilter').innerHTML+ = '';", r'''const noveltyPartners=[...new Set(D.records.map(x=>x.partner))].sort();$('noveltyPartner').innerHTML+=noveltyPartners.map(x=>`<option>${esc(x)}</option>`).join('');
const noveltyLabels={Hungate:'Hungate cultured baseline',SecondaryBaseline:'Secondary cultured baseline',CulturedProject:'Combined cultured collection',Project:'Partner project marker collection',Reference:'GTDB reference genomes'};
const getNovelty=(x,key)=>x.novelty?.find(n=>n.key===key)||{};
const baseNovelty=D.records.map(x=>getNovelty(x,'CulturedProject'));
const novelCount=baseNovelty.filter(n=>n.nearest_identity!=null).length;
const highIdentity=baseNovelty.filter(n=>n.nearest_identity!=null&&n.nearest_identity>=97).length;
const no97=baseNovelty.filter(n=>n.nearest_identity!=null&&(n.matches_97||0)===0).length;
$('noveltyMetrics').innerHTML=[['Compared to cultured collection',novelCount,'Isolates with a cultured comparison'],['Nearest cultured identity ≥97%',highIdentity,'At least one close cultured match'],['No cultured neighbours ≥97%',no97,'No reported neighbour at this cut-off'],['Genome collection comparisons',D.records.filter(x=>x.genome_nearest_hit||x.genome_matches_97!=null).length,'Separate from marker-only comparisons']].map(x=>`<div class="card metric"><strong>${Number(x[1]).toLocaleString()}</strong><span>${esc(x[0])}</span><div class="hint">${esc(x[2])}</div></div>`).join('');
let noveltyFiltered=[],noveltyPage=1;
function getNoveltyFiltered(){const q=$('noveltySearch').value.trim().toLowerCase(),p=$('noveltyPartner').value,k=$('noveltyComparison').value,sort=$('noveltySort').value;let a=D.records.filter(x=>!p||x.partner===p).filter(x=>{const n=getNovelty(x,k);return !q||[x.id,x.partner,x.partner_taxon,n.nearest_hit,n.density_source,x.recommendation].join(' ').toLowerCase().includes(q)});a.sort((x,y)=>{const n=getNovelty(x,k),m=getNovelty(y,k);if(sort==='identity_asc')return (n.nearest_identity??Infinity)-(m.nearest_identity??Infinity);if(sort==='identity_desc')return (m.nearest_identity??-Infinity)-(n.nearest_identity??-Infinity);if(sort==='matches_desc')return (m.matches_97??-1)-(n.matches_97??-1);if(sort==='score_desc')return (m.novelty_score??-Infinity)-(n.novelty_score??-Infinity);return x.id.localeCompare(y.id)});return a}
function renderNovelty(){noveltyFiltered=getNoveltyFiltered();const key=$('noveltyComparison').value,size=Number($('noveltyPageSize').value||50),pages=Math.max(1,Math.ceil(noveltyFiltered.length/size));noveltyPage=Math.min(noveltyPage,pages);const shown=noveltyFiltered.slice((noveltyPage-1)*size,noveltyPage*size);$('noveltyCount').textContent=`Showing ${shown.length.toLocaleString()} of ${noveltyFiltered.length.toLocaleString()} isolates`;$('noveltyBody').innerHTML=shown.map(x=>{const n=getNovelty(x,key);return `<tr><td><b>${esc(x.id)}</b></td><td>${esc(x.partner)}</td><td>${esc(n.nearest_hit||'No hit recorded')}</td><td>${pct(n.nearest_identity)}</td><td>${fmt(n.matches_99)}</td><td>${fmt(n.matches_97)}</td><td>${fmt(n.matches_95)}</td><td>${esc(n.crowding||'—')}</td><td>${pill(n.priority||'Not recorded')}</td><td><button class="btn" data-novelty-open="${esc(x.id)}">Compare</button></td></tr>`}).join('')||'<tr><td colspan="10" class="empty">No novelty records match these filters</td></tr>';pager('noveltyPager',noveltyPage,pages,noveltyFiltered.length,p=>{noveltyPage=p;renderNovelty()})}
['noveltySearch','noveltyPartner','noveltyComparison','noveltySort','noveltyPageSize'].forEach(id=>$(id).addEventListener(id==='noveltySearch'?'input':'change',()=>{noveltyPage=1;renderNovelty()}));
function openNovelty(x){if(!x)return;const rows=x.novelty||[];$('noveltyDetailTitle').textContent=`${x.id} · comparison profile`;$('noveltyDetailSubtitle').textContent=`${x.partner} · ${x.recommendation}`;$('noveltyDetailBody').innerHTML=rows.map(n=>`<tr><td>${esc(n.label)}</td><td>${esc(n.nearest_hit||'No hit recorded')}</td><td>${pct(n.nearest_identity)}</td><td>${pct(n.query_coverage)}</td><td>${fmt(n.matches_99)}</td><td>${fmt(n.matches_97)}</td><td>${fmt(n.matches_95)}</td><td>${esc(n.crowding||'—')}</td><td>${pill(n.priority||'Not recorded')}</td></tr>`).join('');$('noveltyGenomeContext').textContent=`Genome collection: nearest genome hit ${x.genome_nearest_hit||'not recorded'}; identity ${pct(x.genome_nearest_identity)}; ≥97% matches ${fmt(x.genome_matches_97)}. ${x.genome_representation_interpretation||''}`;$('noveltyDetail').style.display='block';$('noveltyDetail').scrollIntoView({behavior:'smooth',block:'nearest'})}
$('noveltyBody').addEventListener('click',e=>{const b=e.target.closest('[data-novelty-open]');if(b)openNovelty(D.records.find(x=>x.id===b.dataset.noveltyOpen))});$('closeNoveltyDetail').addEventListener('click',()=>{$('noveltyDetail').style.display='none'});$('noveltyExport').addEventListener('click',()=>{const key=$('noveltyComparison').value,cols=['id','partner','recommendation'];const rows=[['Isolate','Partner','Recommendation','Comparison','Nearest hit','Identity','Matches >=99%','Matches >=97%','Matches >=95%','Novelty score','Crowding','Priority','Query coverage','Density source','Density count method'].join(',')];for(const x of noveltyFiltered){const n=getNovelty(x,key);const vals=[x.id,x.partner,x.recommendation,noveltyLabels[key],n.nearest_hit,n.nearest_identity,n.matches_99,n.matches_97,n.matches_95,n.novelty_score,n.crowding,n.priority,n.query_coverage,n.density_source,n.density_method];rows.push(vals.map(v=>'"'+String(v??'').replaceAll('"','""')+'"').join(','))}const a=document.createElement('a'),u=URL.createObjectURL(new Blob([rows.join('\r\n')],{type:'text/csv;charset=utf-8'}));a.href=u;a.download='branchmanager_novelty_comparison.csv';a.click();URL.revokeObjectURL(u)});renderNovelty();
const partners=[...new Set(D.records.map(x=>x.partner))].sort(),recs=[...new Set(D.records.map(x=>x.recommendation))].sort(),qcs=[...new Set(D.records.map(x=>x.marker_qc))].sort();$('partnerFilter').innerHTML+=partners.map(x=>`<option>${esc(x)}</option>`).join('');$('recommendationFilter').innerHTML+=recs.map(x=>`<option>${esc(x)}</option>`).join('');$('qcFilter').innerHTML+=qcs.map(x=>`<option>${esc(x)}</option>`).join('');
let filtered=[],isolatePage=1;
function getFiltered(){
 let q=$('search').value.trim().toLowerCase(),p=$('partnerFilter').value,r=$('recommendationFilter').value,c=$('qcFilter').value;
 return D.records.filter(x=>(!p||x.partner===p)&&(!r||x.recommendation===r)&&(!c||x.marker_qc===c)&&(!q||[
 x.id,x.partner,x.recommendation,x.partner_taxon,x.gtdb_taxonomy,x.baseline_class,x.recommendation_reason,x.species_context,x.genome_coverage,
 ...(x.taxonomy||[]).flatMap(t=>[t.database,t.status,t.taxon,t.lineage,t.rank,t.explanation])
 ].join(' ').toLowerCase().includes(q)));
}
function renderTable(){
 filtered=getFiltered();const size=Number($('pageSize').value||50),pages=Math.max(1,Math.ceil(filtered.length/size));isolatePage=Math.min(isolatePage,pages);
 const shown=filtered.slice((isolatePage-1)*size,isolatePage*size);
 $('resultCount').textContent=`Showing ${shown.length.toLocaleString()} isolates (${((isolatePage-1)*size+1).toLocaleString()}–${Math.min(isolatePage*size,filtered.length).toLocaleString()}) of ${filtered.length.toLocaleString()} matches · ${total.toLocaleString()} total`;
 $('isolateBody').innerHTML=shown.map(x=>`<tr><td><button class="btn" data-open="${esc(x.id)}">${esc(x.id)}</button></td><td>${esc(x.partner)}</td><td>${pill(x.marker_qc)}</td><td>${esc(x.recommendation)}</td><td>${pill((x.taxonomy.find(t=>t.database==='GTDB r232')||{}).status)}</td><td>${pct(x.baseline_nearest_identity)}</td><td>${pct(x.project_nearest_identity)}${x.project_identity_warning?` <span class="qc-warning" title="${esc(x.project_identity_warning_reason)}" aria-label="Marker QC warning">⚠ QC</span>`:''}</td></tr>`).join('')||'<tr><td colspan="7" class="empty">No isolates match these filters</td></tr>';
 pager('isolatePager',isolatePage,pages,filtered.length,p=>{isolatePage=p;renderTable()});
}
['search','partnerFilter','recommendationFilter','qcFilter'].forEach(id=>$(id).addEventListener(id==='search'?'input':'change',()=>{isolatePage=1;renderTable()}));
$('pageSize').addEventListener('change',()=>{isolatePage=1;renderTable()});
$('clearFilters').addEventListener('click',()=>{['search','partnerFilter','recommendationFilter','qcFilter'].forEach(id=>$(id).value='');isolatePage=1;renderTable()});
$('isolateBody').addEventListener('click',e=>{let b=e.target.closest('[data-open]');if(b)openDetail(D.records.find(x=>x.id===b.dataset.open))});renderTable();
function box(label,value){return `<div class="detailbox"><small>${esc(label)}</small><strong>${esc(value==null||value===''?'Not available':value)}</strong></div>`}
function openDetail(x){if(!x)return;$('detailTitle').textContent=x.id;$('detailSubtitle').textContent=`${x.partner} · ${x.recommendation}`;let tax=x.taxonomy.map(t=>`<div class="taxcard"><strong>${esc(t.database)}</strong><div>${pill(t.status)}</div><p><b>Best taxon:</b> ${esc(t.taxon)}</p><p><b>Identity:</b> ${pct(t.identity)} · <b>Query coverage:</b> ${pct(t.query_coverage)} · <b>Reference coverage:</b> ${pct(t.reference_coverage)}</p><p><b>Assigned rank:</b> ${esc(t.rank||'Not available')}</p>${t.explanation?`<p>${esc(t.explanation)}</p>`:''}<p class="lineage">${esc(t.lineage)}</p></div>`).join('');let fig='';if(x.local_tree_figure&&!x.local_tree_figure.match(/^(\.\.|\/)/))fig=`<p><a href="${esc(x.local_tree_figure.startsWith('assessment/')?x.local_tree_figure:'assessment/'+x.local_tree_figure)}" target="_blank" rel="noopener">Open local phylogenetic figure</a></p>`;$('detailContent').innerHTML=`<div class="detailgrid">${box('Marker QC',x.marker_qc)}${box('Manual review',x.marker_review)}${box('Evidence quality',x.evidence_quality)}${box('SAB decision',x.sab_decision||'Not recorded')}${box('Sequencing role',x.sequencing_role)}${box('Partner-supplied taxon',x.partner_taxon)}${box('Cultured baseline category',x.baseline_class)}<div class="detailbox"><small>Closest cultured identity</small><strong>${esc(pct(x.baseline_nearest_identity))}</strong><p><b>Nearest hit:</b> ${esc(x.baseline_nearest_hit||'No cultured hit recorded')}</p></div><div class="detailbox"><small>Closest project identity</small><strong>${esc(pct(x.project_nearest_identity))}</strong><p><b>Nearest hit:</b> ${esc(x.project_nearest_hit||'No project hit recorded')}</p>${x.project_identity_warning?`<p class="qc-warning-text">⚠ ${esc(x.project_identity_warning_reason)}</p>`:''}</div>${box('Pangenome group',x.pangenome_rank+' · '+x.pangenome_taxon)}${box('Pangenome gap / target',fmt(x.pangenome_gap)+' / '+fmt(x.pangenome_target))}</div><div class="detailbox" style="margin-top:12px"><small>Closest cultured and project sequence matches</small><div class="taxgrid">${[...(x.cultured_match_details||[]),{key:'Project',...x.project_match_details}].map(m=>`<div class="taxcard"><strong>${esc(m.label)}</strong><p><b>Identity:</b> ${pct(m.identity)}</p><p><b>Nearest hit:</b> ${esc(m.hit||'No hit recorded')}</p><p><b>Query coverage:</b> ${pct(m.query_coverage)} · <b>Aligned length:</b> ${fmt(m.alignment_length)} bp</p></div>`).join('')}</div></div><div class="detailbox"><small>Recommendation explanation</small><p>${esc(({ 'SECONDARY - STRAIN DIVERSITY':'This isolate is identified as a candidate that may add strain diversity to the current collection.','REVIEW BEFORE SELECTION':'The available evidence indicates that this isolate should be reviewed before a sequencing decision is made.','EXCLUDE - BASELINE REDUNDANT':'The assessment identifies a close match in the cultured baseline, so this isolate is not prioritised as an additional genome at present.','LOWER PRIORITY - TARGET MET':'The configured pangenome target is reported as met for the relevant group, so this isolate has a lower sequencing priority.','ALREADY SELECTED - GENOME PENDING':'This isolate is already selected for genome sequencing, but the genome is recorded as pending.','ALREADY SEQUENCED':'A genome is already recorded for this isolate.'})[x.recommendation]||'No plain-language definition is available for this recommendation.')}</p><p><b>Assessment rationale:</b> ${esc(x.recommendation_reason||'No detailed rationale recorded.')}</p></div><div class="detailbox" style="margin-top:12px"><small>Taxonomic evidence across reference databases</small><div class="taxgrid">${tax}</div></div><div class="detailbox" style="margin-top:12px"><small>Additional biological context</small><p>${esc(x.species_context||'No additional species context recorded.')}</p><p>${esc(x.genome_coverage||'')}</p>${fig}</div>`;$('detail').classList.add('open');document.body.style.overflow='hidden'}
$('closeDetail').addEventListener('click',()=>{$('detail').classList.remove('open');document.body.style.overflow='' });$('detail').addEventListener('click',e=>{if(e.target===$('detail'))$('closeDetail').click()});document.addEventListener('keydown',e=>{if(e.key==='Escape'){if($('pangenomeDetail').classList.contains('open'))$('closePangenomeDetail').click();else $('closeDetail').click()}});
$('exportBtn').addEventListener('click',()=>{let rows=filtered,cols=['id','partner','marker_qc','recommendation','sequencing_role','evidence_quality','partner_taxon','baseline_class','baseline_nearest_identity','baseline_nearest_hit','project_nearest_identity','project_nearest_hit','pangenome_rank','pangenome_taxon','pangenome_gap','recommendation_reason'];let csv=[cols.join(','),...rows.map(x=>cols.map(k=>'"'+String(x[k]??'').replaceAll('"','""')+'"').join(','))].join('\r\n');let a=document.createElement('a'),u=URL.createObjectURL(new Blob([csv],{type:'text/csv;charset=utf-8'}));a.href=u;a.download='branchmanager_filtered_isolates.csv';a.click();URL.revokeObjectURL(u)});
const navLinks=[...document.querySelectorAll('.nav a')],sections=navLinks.map(a=>document.querySelector(a.getAttribute('href'))).filter(Boolean);
if('IntersectionObserver' in window){const io=new IntersectionObserver(entries=>entries.forEach(entry=>{if(entry.isIntersecting){navLinks.forEach(a=>a.classList.toggle('active',a.getAttribute('href')==='#'+entry.target.id))}}),{rootMargin:'-18% 0px -72% 0px'});sections.forEach(s=>io.observe(s))}
const topButton=$('toTop');window.addEventListener('scroll',()=>topButton.classList.toggle('visible',window.scrollY>450),{passive:true});topButton.addEventListener('click',()=>window.scrollTo({top:0,behavior:'smooth'}));
document.addEventListener('keydown',e=>{if(e.key==='/'&&!['INPUT','SELECT','TEXTAREA'].includes(document.activeElement.tagName)){e.preventDefault();$('search').focus()}});
$('runDetails').innerHTML=[['Workflow',D.run.workflow],['Status',D.run.status],['Started',D.run.started_at],['Completed',D.run.completed_at],['BranchManager version',D.run.version],['Workflow version',D.run.workflow_version],['Database schema version',D.run.database_schema_version],['Database file SHA256',D.run.database_sha256],['Live database overlay',D.run.database_overlay?'Yes':'No'],['Git revision',D.run.git_revision]].map(x=>`<div class="detailbox"><small>${esc(x[0])}</small><strong>${esc(x[1]||'Not recorded')}</strong></div>`).join('');
})();
</script></body></html>''')
    html_doc = html_doc.replace('__DATA__', payload)
    output.write_text(html_doc, encoding='utf-8')
    return str(output)


def _latest_two_snapshots(db):
    with db.connect() as conn:
        ids = [row[0] for row in conn.execute(
            'SELECT snapshot_id FROM assessment_snapshots GROUP BY snapshot_id '
            'ORDER BY MAX(created_at) DESC, MAX(rowid) DESC LIMIT 2'
        ).fetchall()]
        snapshots = []
        for snapshot_id in ids:
            rows = conn.execute(
                'SELECT sequence_id, assessment_json FROM assessment_snapshots WHERE snapshot_id = ?',
                (snapshot_id,),
            ).fetchall()
            snapshots.append((snapshot_id, {
                str(sequence_id): json.loads(payload) for sequence_id, payload in rows
            }))
    return snapshots


def write_decision_changes(db, path: str | Path) -> str:
    snapshots = _latest_two_snapshots(db)
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    headers = [
        'SequenceID', 'PreviousSnapshot', 'CurrentSnapshot', 'PreviousDecision',
        'CurrentDecision', 'DecisionChanged', 'PreviousSetRole', 'CurrentSetRole',
        'PreviousPangenomeGap', 'CurrentPangenomeGap', 'ChangeReason',
    ]
    with open(output, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter='\t')
        writer.writeheader()
        if len(snapshots) < 2:
            return str(output)
        (current_id, current), (previous_id, previous) = snapshots
        for sequence_id in sorted(set(previous) | set(current)):
            old = previous.get(sequence_id, {})
            new = current.get(sequence_id, {})
            old_decision = build_selection_decision(old)['decision'] if old else 'NOT_ASSESSED'
            new_decision = build_selection_decision(new)['decision'] if new else 'NOT_ASSESSED'
            changes = []
            for label, key in (
                ('set role', 'sequencing_set_role'), ('pangenome gap', 'pangenome_gap'),
                ('nearest genome', 'nearest_genome_hit'), ('marker QC', 'marker_qc_flag'),
            ):
                if str(old.get(key, 'NA')) != str(new.get(key, 'NA')):
                    changes.append(f'{label}: {old.get(key, "NA")} -> {new.get(key, "NA")}')
            writer.writerow({
                'SequenceID': sequence_id,
                'PreviousSnapshot': previous_id,
                'CurrentSnapshot': current_id,
                'PreviousDecision': old_decision,
                'CurrentDecision': new_decision,
                'DecisionChanged': 'yes' if old_decision != new_decision else 'no',
                'PreviousSetRole': old.get('sequencing_set_role', 'NA'),
                'CurrentSetRole': new.get('sequencing_set_role', 'NA'),
                'PreviousPangenomeGap': old.get('pangenome_gap', 'NA'),
                'CurrentPangenomeGap': new.get('pangenome_gap', 'NA'),
                'ChangeReason': '; '.join(changes) or 'recalculation did not change tracked evidence',
            })
    return str(output)


def _project_root_from_db(db) -> Path | None:
    db_path = str(getattr(db, 'path', '') or '')
    if not db_path or db_path == ':memory:':
        return None
    try:
        return Path(db_path).expanduser().resolve().parent
    except Exception:
        return None


def _query_dataset_rows(conn, datasets=None) -> list[tuple[str, str, str, int, int, int]]:
    params = ()
    where = ''
    if datasets is not None:
        names = [str(name) for name in datasets if str(name or '').strip()]
        if not names:
            return []
        where = 'WHERE s.dataset IN (' + ','.join('?' for _ in names) + ') '
        params = tuple(names)
    role_columns = {row[1] for row in conn.execute('PRAGMA table_info(dataset_roles)').fetchall()}
    tier_expr = 'r.baseline_tier' if 'baseline_tier' in role_columns else "''"
    rows = conn.execute(
        'SELECT s.dataset, COALESCE(r.role, "unregistered"), COALESCE(' + tier_expr + ', ""), '
        'CASE WHEN COALESCE(r.genomes_available, 0) = 1 THEN COUNT(s.id) '
        '     ELSE COALESCE(SUM(COALESCE(m.selected_for_wgs, 0)), 0) END, '
        "COALESCE(SUM(CASE WHEN UPPER(COALESCE(i.status, '')) = 'PROPOSED' "
        'THEN 1 ELSE 0 END), 0), '
        'COUNT(s.id) FROM sequences s '
        'LEFT JOIN dataset_roles r ON r.dataset=s.dataset '
        'LEFT JOIN sequencing_metadata m ON m.id=s.id '
        'LEFT JOIN isolate_status i ON i.sequence_id=s.id '
        f'{where}'
        'GROUP BY s.dataset, r.role, r.genomes_available, ' + tier_expr + ' ORDER BY s.dataset',
        params,
    ).fetchall()
    return [
        (
            str(dataset or ''),
            str(role or 'unregistered'),
            str(baseline_tier or ''),
            int(genomes or 0),
            int(proposed or 0),
            int(sequences or 0),
        )
        for dataset, role, baseline_tier, genomes, proposed, sequences in rows
        if str(dataset or '').strip()
    ]


def _query_status_rows_for_datasets(conn, datasets) -> list[tuple[str, int]]:
    names = [str(name) for name in datasets if str(name or '').strip()]
    if not names:
        return []
    rows = conn.execute(
        "SELECT COALESCE(i.status, 'RECEIVED'), COUNT(*) FROM sequences s "
        'LEFT JOIN isolate_status i ON i.sequence_id=s.id '
        'WHERE s.dataset IN (' + ','.join('?' for _ in names) + ') '
        "GROUP BY COALESCE(i.status, 'RECEIVED') "
        "ORDER BY COALESCE(i.status, 'RECEIVED')",
        tuple(names),
    ).fetchall()
    return [(str(status or 'RECEIVED'), int(count or 0)) for status, count in rows]


def _query_isolate_rows(conn, datasets=None):
    params = ()
    where = ''
    if datasets is not None:
        names = [str(name) for name in datasets if str(name or '').strip()]
        if not names:
            return []
        where = 'WHERE s.dataset IN (' + ','.join('?' for _ in names) + ') '
        params = tuple(names)
    return conn.execute(
        'SELECT s.id, s.dataset, COALESCE(m.partner_id, ""), '
        'COALESCE(i.status, "RECEIVED"), COALESCE(p.marker_qc_class, "UNVERIFIED"), '
        'COALESCE(m.selected_for_sequencing, 0), COALESCE(m.selected_for_wgs, 0) FROM sequences s '
        'LEFT JOIN sequencing_metadata m ON m.id=s.id '
        'LEFT JOIN isolate_status i ON i.sequence_id=s.id '
        'LEFT JOIN sequence_provenance p ON p.sequence_id=s.id '
        f'{where}'
        'ORDER BY s.id',
        params,
    ).fetchall()


def _command_value(command, flag: str) -> str:
    if not isinstance(command, list):
        return ''
    prefix = f'{flag}='
    for index, token in enumerate(command):
        value = str(token)
        if value == flag and index + 1 < len(command):
            return str(command[index + 1])
        if value.startswith(prefix):
            return value[len(prefix):]
    return ''


def _dataset_from_path_hint(value: object) -> str:
    raw = str(value or '').strip()
    if not raw:
        return ''
    try:
        parts = [part for part in Path(raw).parts if part not in {'', '/'}]
    except TypeError:
        return ''
    for index, part in enumerate(parts):
        if part == 'mailroom' and index > 0:
            return parts[index - 1]
    for part in reversed(parts):
        if part.startswith('01_interview_') and len(part) > len('01_interview_'):
            return part[len('01_interview_'):]
    stage_names = {
        '00_it_desk', '01_onboarding', '02_paper_trail_merge_meeting',
        '03_performance_review_hiring_panel',
    }
    for index, part in enumerate(parts):
        if part in stage_names and index > 0:
            return parts[index - 1]
    return ''


def _manifest_dataset(data: dict, path: Path | None = None) -> str:
    command_dataset = str(_command_value(data.get('command'), '--dataset') or '').strip()
    if command_dataset:
        return command_dataset
    direct = str(data.get('dataset') or '').strip()
    if direct:
        return direct
    mailroom = _command_value(data.get('command'), '--mailroom')
    dataset = _dataset_from_path_hint(mailroom)
    if dataset:
        return dataset
    for record in data.get('inputs', []):
        if str(record.get('role') or '') in {'sample_map', 'mailroom_summary'}:
            dataset = _dataset_from_path_hint(record.get('path'))
            if dataset:
                return dataset
    if path is not None:
        dataset = _dataset_from_path_hint(path)
        if dataset:
            return dataset
    return ''


def _manifest_run_id(path: Path, data: dict, project_root: Path) -> str:
    run_id = str(data.get('run_id') or '').strip()
    if run_id:
        return run_id
    workflow = str(data.get('workflow') or 'workflow').strip()
    dataset = _manifest_dataset(data, path)
    started = str(data.get('started_at') or '').strip()
    try:
        relative = path.relative_to(project_root).as_posix()
    except ValueError:
        relative = path.as_posix()
    return ':'.join(part for part in (workflow, dataset, started, relative) if part)


def _discover_run_manifests(project_root: Path | None, outdir: Path):
    if project_root is None or not project_root.exists():
        return []
    output = outdir.expanduser().resolve()
    rows = []
    for path in sorted(project_root.rglob('run_manifest.json')):
        resolved = path.resolve()
        try:
            resolved.relative_to(output)
            continue
        except ValueError:
            pass
        try:
            data = json.loads(path.read_text())
        except (OSError, ValueError, TypeError):
            continue
        workflow = str(data.get('workflow') or '').strip()
        if not workflow or workflow == 'annual_report':
            continue
        rows.append((
            _manifest_run_id(path, data, project_root),
            workflow,
            _manifest_dataset(data, path),
            str(data.get('status') or '').strip(),
            str(data.get('started_at') or '').strip(),
            str(data.get('completed_at') or '').strip(),
            str(path),
        ))
    return rows


def _resolve_manifest_path(value: object, manifest_path: Path, project_root: Path | None) -> Path:
    raw = str(value or '').strip()
    candidate = Path(raw).expanduser()
    if candidate.exists():
        return candidate
    if project_root is not None:
        parts = list(candidate.parts)
        if project_root.name in parts:
            index = parts.index(project_root.name)
            remapped = project_root.joinpath(*parts[index + 1:])
            if remapped.exists():
                return remapped
    local = manifest_path.parent / candidate.name
    if local.exists():
        return local
    return candidate


def _output_paths(data: dict, role: str, manifest_path: Path, project_root: Path | None) -> list[Path]:
    return [
        _resolve_manifest_path(record.get('path'), manifest_path, project_root)
        for record in data.get('outputs', [])
        if str(record.get('role') or '') == role and str(record.get('path') or '').strip()
    ]


def _first_output_path(data: dict, role: str, manifest_path: Path, project_root: Path | None) -> Path | None:
    paths = _output_paths(data, role, manifest_path, project_root)
    return paths[0] if paths else None


def _read_tsv_rows(path: Path | None) -> list[dict]:
    if path is None or not path.exists() or not path.is_file():
        return []
    try:
        with open(path, newline='', encoding='utf-8-sig') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            if not reader.fieldnames:
                return []
            reader.fieldnames = [str(name or '').strip().lstrip('\ufeff') for name in reader.fieldnames]
            return [
                {str(key or '').strip().lstrip('\ufeff'): (value.strip() if isinstance(value, str) else value)
                 for key, value in row.items()}
                for row in reader
            ]
    except (OSError, csv.Error):
        return []


def _manifest_input_signature(data: dict) -> str:
    for role in ('sample_map', 'mailroom_summary'):
        for record in data.get('inputs', []):
            if str(record.get('role') or '') == role and record.get('sha256'):
                return str(record.get('sha256'))
    return ''


def _sequence_signature(rows: list[dict]) -> str:
    ids = sorted({str(row.get('SequenceID') or row.get('sequence_id') or '').strip() for row in rows})
    ids = [sid for sid in ids if sid]
    if not ids:
        return ''
    digest = hashlib.sha1('\n'.join(ids).encode('utf-8')).hexdigest()
    return f'sequence_ids:{digest}'


def _qc_summary_from_manifest(path: Path, data: dict, project_root: Path | None) -> dict | None:
    workflow = str(data.get('workflow') or '').strip()
    if workflow not in {'interview', 'paper_trail'}:
        return None
    dataset = _manifest_dataset(data, path)
    assembly_path = _first_output_path(data, 'assembly_tsv', path, project_root)
    read_qc_path = _first_output_path(data, 'read_qc_tsv', path, project_root)
    failed_path = _first_output_path(data, 'failed_manifest_tsv', path, project_root)
    failed_read_path = _first_output_path(data, 'failed_read_manifest_tsv', path, project_root)
    visual_manifest_path = _first_output_path(data, 'visual_manifest_tsv', path, project_root)
    assembly_rows = _read_tsv_rows(assembly_path)
    read_rows = _read_tsv_rows(read_qc_path)
    failed_rows = _read_tsv_rows(failed_path)
    failed_read_rows = _read_tsv_rows(failed_read_path)
    visual_rows = _read_tsv_rows(visual_manifest_path)
    if not any((assembly_rows, read_rows, failed_rows, failed_read_rows, visual_rows)):
        return None
    passed_sequences = sum(
        1 for row in assembly_rows
        if str(row.get('QCClass') or '').strip().upper().startswith('PASS')
    )
    accepted_sequences = sum(
        1 for row in assembly_rows
        if str(row.get('Recommendation') or '').strip().upper() == 'ACCEPT'
    )
    review_sequences = sum(
        1 for row in assembly_rows
        if str(row.get('Recommendation') or '').strip().upper() == 'MANUAL_REVIEW'
        or str(row.get('QCClass') or '').strip().upper() == 'PASS_WITH_WARNINGS'
    )
    failed_sequences = len(failed_rows) or sum(
        1 for row in assembly_rows
        if str(row.get('Recommendation') or '').strip().upper() not in {'', 'ACCEPT'}
    )
    failed_read_files = len(failed_read_rows) or sum(
        1 for row in read_rows
        if str(row.get('Status') or '').strip().lower() not in {'', 'kept'}
    )
    report_counts = Counter(str(row.get('Report') or 'visual_pages') for row in visual_rows)
    sequence_signature = _sequence_signature(assembly_rows)
    signature = sequence_signature or _manifest_input_signature(data) or str(path)
    preference = (
        1 if workflow != 'interview' else 0,
        1 if dataset else 0,
        passed_sequences + failed_sequences,
        len(read_rows),
    )
    return {
        'dataset': dataset,
        'workflow': workflow,
        'started_at': str(data.get('started_at') or ''),
        'manifest_path': str(path),
        'signature': signature,
        'preference': preference,
        'passed_sequence_files': passed_sequences,
        'accepted_sequence_files': accepted_sequences,
        'review_sequence_files': review_sequences,
        'failed_sequence_files': failed_sequences,
        'read_files': len(read_rows),
        'failed_read_files': failed_read_files,
        'visual_pages': len(visual_rows),
        'read_error_pages': report_counts.get('read_error_profiles', len(_output_paths(data, 'read_error_pngs', path, project_root))),
        'chromatogram_pages': report_counts.get('trace_chromatograms', len(_output_paths(data, 'chromatogram_pngs', path, project_root))),
        'assembly_pages': report_counts.get('assembly_overviews', len(_output_paths(data, 'assembly_pngs', path, project_root))),
        'visual_manifest_path': str(visual_manifest_path or ''),
    }


def _discover_qc_summaries(project_root: Path | None, outdir: Path) -> list[dict]:
    if project_root is None or not project_root.exists():
        return []
    output = outdir.expanduser().resolve()
    summaries: dict[tuple[str, str], dict] = {}
    for path in sorted(project_root.rglob('run_manifest.json')):
        resolved = path.resolve()
        try:
            resolved.relative_to(output)
            continue
        except ValueError:
            pass
        try:
            data = json.loads(path.read_text())
        except (OSError, ValueError, TypeError):
            continue
        summary = _qc_summary_from_manifest(path, data, project_root)
        if not summary:
            continue
        key = (summary['dataset'], summary['signature'])
        previous = summaries.get(key)
        if previous is None or summary['preference'] > previous['preference']:
            summaries[key] = summary
    return sorted(
        summaries.values(),
        key=lambda row: (str(row.get('dataset') or ''), str(row.get('started_at') or ''), str(row.get('manifest_path') or '')),
    )


def _query_dataset_marker_qc_rows(conn) -> dict[str, dict[str, int]]:
    rows = conn.execute(
        "SELECT s.dataset, "
        "SUM(CASE WHEN UPPER(COALESCE(p.marker_qc_class, '')) LIKE 'PASS%' THEN 1 ELSE 0 END), "
        "SUM(CASE WHEN UPPER(COALESCE(p.marker_qc_class, '')) IN "
        "('PASS_WITH_WARNINGS', 'MARKER_QC_REVIEW_REQUIRED') "
        "OR UPPER(COALESCE(i.status, '')) = 'TRACE_REVIEW' THEN 1 ELSE 0 END), "
        "SUM(CASE WHEN UPPER(COALESCE(p.marker_qc_class, '')) LIKE 'FAIL%' "
        "OR UPPER(COALESCE(i.status, '')) LIKE '%FAILED%' THEN 1 ELSE 0 END) "
        "FROM sequences s "
        "LEFT JOIN sequence_provenance p ON p.sequence_id=s.id "
        "LEFT JOIN isolate_status i ON i.sequence_id=s.id "
        "GROUP BY s.dataset"
    ).fetchall()
    return {
        str(dataset or ''): {
            'active_marker_qc_passed': int(passed or 0),
            'active_marker_qc_review': int(review or 0),
            'active_marker_qc_failed': int(failed or 0),
        }
        for dataset, passed, review, failed in rows
        if str(dataset or '').strip()
    }


def _aggregate_qc_summaries(qc_summaries: list[dict]) -> dict[str, dict[str, int]]:
    totals = defaultdict(lambda: {
        'passed_sequence_files': 0,
        'accepted_sequence_files': 0,
        'review_sequence_files': 0,
        'failed_sequence_files': 0,
        'read_files': 0,
        'failed_read_files': 0,
        'visual_pages': 0,
    })
    for row in qc_summaries:
        dataset = str(row.get('dataset') or '').strip()
        if not dataset:
            continue
        for key in totals[dataset]:
            totals[dataset][key] += int(row.get(key) or 0)
    return dict(totals)


def _html_path_target(path: object, output: Path) -> str:
    raw = str(path or '').strip()
    if not raw:
        return ''
    p = Path(raw).expanduser()
    if p.exists():
        try:
            # ``output`` is the HTML file, not its containing directory.
            # Relative browser URLs are resolved from the HTML file's parent.
            base = output.resolve().parent if output.suffix.lower() in {'.html', '.htm'} else output.resolve()
            return Path(os.path.relpath(p.resolve(), base)).as_posix()
        except (OSError, ValueError):
            return p.resolve().as_posix()
    return raw


def _html_path_link(path: object, output: Path, label: str | None = None) -> str:
    raw = str(path or '').strip()
    if not raw:
        return ''
    target = _html_path_target(raw, output)
    return f'<a href="{html.escape(target)}">{html.escape(label or Path(raw).name or raw)}</a>'


def _resolve_visual_file(manifest_path: Path, file_value: object) -> Path | None:
    raw = str(file_value or '').strip()
    if not raw:
        return None
    candidate = Path(raw).expanduser()
    candidates = [candidate] if candidate.is_absolute() else []
    if not candidate.is_absolute():
        candidates.extend([
            manifest_path.parent.parent / candidate,
            manifest_path.parent / candidate,
            manifest_path.parent / candidate.name,
        ])
    for item in candidates:
        if item.exists() and item.is_file() and item.suffix.lower() == '.png':
            return item
    return None


def _paper_trail_gallery_rows(
    dataset: str, workflow: str, manifest_path: Path, output: Path, *, max_per_report: int = 1,
) -> list[dict]:
    rows = []
    counts: Counter[str] = Counter()
    for row in _read_tsv_rows(manifest_path):
        report = str(row.get('Report') or 'visual_report').strip()
        if counts[report] >= max_per_report:
            continue
        image = _resolve_visual_file(manifest_path, row.get('File'))
        if image is None:
            continue
        page = str(row.get('Page') or '').strip()
        page_count = str(row.get('PageCount') or '').strip()
        first_record = str(row.get('FirstRecord') or '').strip()
        last_record = str(row.get('LastRecord') or '').strip()
        range_text = ''
        if first_record and last_record:
            range_text = f' ({first_record} to {last_record})'
        caption = report.replace('_', ' ').title()
        if page:
            caption += f' page {page}'
            if page_count:
                caption += f' of {page_count}'
        caption += range_text
        rows.append({
            'dataset': dataset,
            'workflow': workflow,
            'report': report.replace('_', ' ').title(),
            'caption': caption,
            'image_path': image.as_posix(),
            'link_path': image.as_posix(),
            'src': _html_path_target(image, output),
            'href': _html_path_target(image, output),
        })
        counts[report] += 1
    return rows


def _neighbourhood_gallery_rows(
    dataset: str, workflow: str, manifest_path: Path, output: Path, *, max_images: int = 2,
) -> list[dict]:
    rows = []
    seen = set()
    for row in _read_tsv_rows(manifest_path):
        if str(row.get('Status') or '').strip().lower() != 'rendered':
            continue
        figure = str(row.get('Figure') or '').strip()
        if not figure or figure in seen:
            continue
        image = _resolve_visual_file(manifest_path, figure)
        if image is None:
            continue
        assessed = str(row.get('AssessedSequencesInFigure') or '').strip()
        leaves = str(row.get('TreeLeavesShown') or '').strip()
        caption = Path(figure).stem.replace('_', ' ').title()
        if assessed:
            caption += f' ({assessed})'
        if leaves:
            caption += f'; {leaves} leaves shown'
        rows.append({
            'dataset': dataset,
            'workflow': workflow,
            'report': 'Local Neighbourhood',
            'caption': caption,
            'image_path': image.as_posix(),
            'link_path': image.as_posix(),
            'src': _html_path_target(image, output),
            'href': _html_path_target(image, output),
        })
        seen.add(figure)
        if len(rows) >= max_images:
            break
    return rows


def _build_visual_gallery_rows(visual_report_rows: list[tuple], output: Path, *, max_images: int = 12) -> list[dict]:
    gallery_rows = []
    for dataset, workflow, report_name, _items, path_value in visual_report_rows:
        manifest_path = Path(str(path_value or '')).expanduser()
        if not manifest_path.exists():
            continue
        if str(report_name) == 'Paper Trail visual pages':
            gallery_rows.extend(_paper_trail_gallery_rows(str(dataset), str(workflow), manifest_path, output))
        elif str(report_name) == 'Local neighbourhood figures':
            gallery_rows.extend(_neighbourhood_gallery_rows(str(dataset), str(workflow), manifest_path, output))
        if len(gallery_rows) >= max_images:
            return gallery_rows[:max_images]
    return gallery_rows[:max_images]


def _visual_gallery_html(rows: list[dict]) -> str:
    if not rows:
        return '<p class="muted">No embeddable PNG previews were discovered for this report.</p>'
    cards = []
    for row in rows:
        title = ' / '.join(part for part in (row.get('dataset'), row.get('report')) if part)
        cards.append(
            '<figure class="visual-card">'
            f'<a href="{html.escape(str(row.get("href") or ""))}">'
            f'<img src="{html.escape(str(row.get("src") or ""))}" alt="{html.escape(str(row.get("caption") or title))}" loading="lazy"></a>'
            f'<figcaption><b>{html.escape(title)}</b><span>{html.escape(str(row.get("caption") or ""))}</span></figcaption>'
            '</figure>'
        )
    return '<div class="gallery">' + ''.join(cards) + '</div>'


def _annual_report_readme_markdown() -> str:
    return """# BranchManager Annual Report Guide

This folder is a point-in-time project overview built from the BranchManager database, workflow manifests, and run outputs.

## How to Read the HTML

- Summary cards give project-wide counts for baseline markers, active candidate markers, proposed sequencing commitments, marker-QC outcomes, and visual pages.
- The Datasets table combines database state with Paper Trail / Interview file-level QC summaries. `QC-passed sequence files` includes accepted files and files that passed with review warnings; `Accepted sequence files` are ready for downstream review without manual marker review; `Review sequence files` passed with warnings; `QC-failed sequence files` should be resequenced or excluded.
- Current recommendations come from the latest Performance Review assessment snapshot. Baseline-redundant candidates are near-identical to cultured baseline markers at the configured identity and coverage thresholds.
- Embedded visual previews are a curated subset of the available PNG outputs. Click any preview to open the full-size PNG, and use `project_visual_reports.tsv` for the full manifest of visual artifacts.
- Workflow runs show reproducibility context: workflow name, inferred dataset, status, timestamps, and manifest path.

## Key TSV Ledgers

- `project_dataset_summary.tsv`: dataset-level counts, marker-QC summaries, and visual page counts.
- `project_candidate_overview.tsv`: current candidate-level recommendations and nearest-hit context.
- `project_visual_reports.tsv`: visual report manifests, dashboards, and local-neighbourhood figure manifests.
- `project_visual_gallery.tsv`: the exact PNG previews embedded in the annual HTML.
- `project_isolate_ledger.tsv`: active isolate/project state from the database.
- `project_genome_ledger.tsv`: genome records and genome QC state.
- `decision_changes.tsv`: latest tracked decision deltas between assessment snapshots.

## Interpretation Notes

BranchManager treats marker quality as decision evidence. Low marker evidence is held for review before sequencing selection; high versus moderate marker evidence is also used when otherwise comparable candidates are ranked. Visuals are supporting evidence, not replacement evidence: the TSV ledgers remain the auditable source of the numerical decisions.
"""


def _annual_report_readme_html() -> str:
    return (
        '<!doctype html><html><head><meta charset="utf-8"><title>BranchManager Annual Report Guide</title>'
        '<style>body{font:15px system-ui;margin:32px;color:#17202a;max-width:920px;line-height:1.45}'
        'h1,h2{color:#123b2d}li{margin:6px 0}code{background:#eef4f1;padding:1px 4px;border-radius:4px}</style></head><body>'
        '<h1>BranchManager Annual Report Guide</h1>'
        '<p>This folder is a point-in-time project overview built from the BranchManager database, workflow manifests, and run outputs.</p>'
        '<h2>How To Read The HTML</h2>'
        '<ul>'
        '<li>Summary cards give project-wide counts for baseline markers, active candidate markers, proposed sequencing commitments, marker-QC outcomes, and visual pages.</li>'
        '<li>The Datasets table combines database state with Paper Trail / Interview file-level QC summaries. QC-passed sequence files include accepted files and files that passed with review warnings.</li>'
        '<li>Current recommendations come from the latest Performance Review assessment snapshot. Baseline-redundant candidates are near-identical to cultured baseline markers at the configured identity and coverage thresholds.</li>'
        '<li>Embedded visual previews are a curated subset of the available PNG outputs. Click any preview to open the full-size PNG.</li>'
        '<li>Workflow runs show reproducibility context: workflow name, inferred dataset, status, timestamps, and manifest path.</li>'
        '</ul>'
        '<h2>Key TSV Ledgers</h2>'
        '<ul>'
        '<li><code>project_dataset_summary.tsv</code>: dataset-level counts, marker-QC summaries, and visual page counts.</li>'
        '<li><code>project_candidate_overview.tsv</code>: current candidate-level recommendations and nearest-hit context.</li>'
        '<li><code>project_visual_reports.tsv</code>: visual report manifests, dashboards, and local-neighbourhood figure manifests.</li>'
        '<li><code>project_visual_gallery.tsv</code>: the exact PNG previews embedded in the annual HTML.</li>'
        '<li><code>project_isolate_ledger.tsv</code>, <code>project_genome_ledger.tsv</code>, and <code>decision_changes.tsv</code>: audit ledgers for active isolates, genome records, and decision changes.</li>'
        '</ul>'
        '<h2>Interpretation Notes</h2>'
        '<p>BranchManager treats marker quality as decision evidence. Low marker evidence is held for review before sequencing selection; high versus moderate marker evidence is also used when otherwise comparable candidates are ranked. Visuals are supporting evidence, not replacement evidence: the TSV ledgers remain the auditable source of the numerical decisions.</p>'
        '</body></html>\n'
    )


def _discover_visual_report_rows(project_root: Path | None, outdir: Path, qc_summaries: list[dict]) -> list[tuple[str, str, str, int, str]]:
    rows: list[tuple[str, str, str, int, str]] = []
    seen: set[tuple[str, str]] = set()
    for summary in qc_summaries:
        visual_manifest = str(summary.get('visual_manifest_path') or '').strip()
        if visual_manifest:
            rows.append((
                str(summary.get('dataset') or ''),
                str(summary.get('workflow') or ''),
                'Paper Trail visual pages',
                int(summary.get('visual_pages') or 0),
                visual_manifest,
            ))
            seen.add(('paper_trail_visuals', visual_manifest))
    if project_root is None or not project_root.exists():
        return rows
    output = outdir.expanduser().resolve()
    for path in sorted(project_root.rglob('run_manifest.json')):
        resolved = path.resolve()
        try:
            resolved.relative_to(output)
            continue
        except ValueError:
            pass
        try:
            data = json.loads(path.read_text())
        except (OSError, ValueError, TypeError):
            continue
        workflow = str(data.get('workflow') or '').strip()
        dataset = _manifest_dataset(data, path)
        for dashboard in _output_paths(data, 'performance_review_dashboard', path, project_root):
            key = ('performance_review_dashboard', dashboard.as_posix())
            if key in seen:
                continue
            rows.append((dataset, workflow, 'Performance Review dashboard', 1, dashboard.as_posix()))
            seen.add(key)
        neighbourhood_manifest = path.parent / 'assessment' / 'neighbourhoods' / 'neighbourhood_manifest.tsv'
        if neighbourhood_manifest.exists():
            key = ('local_neighbourhoods', neighbourhood_manifest.as_posix())
            if key not in seen:
                manifest_rows = _read_tsv_rows(neighbourhood_manifest)
                rendered = sum(1 for row in manifest_rows if str(row.get('Status') or '').lower() == 'rendered')
                rows.append((dataset, workflow, 'Local neighbourhood figures', rendered, neighbourhood_manifest.as_posix()))
                seen.add(key)
    return sorted(rows, key=lambda row: (str(row[0]), str(row[1]), str(row[2]), str(row[4])))


def _infer_dataset_role(dataset: str, workflow: str = '') -> str:
    name = str(dataset or '').lower()
    if workflow == 'filing_cabinet' or 'hungate' in name:
        return 'baseline'
    return 'candidate'


def _infer_baseline_tier(dataset: str, role: str = '') -> str:
    if str(role or '').lower() != 'baseline':
        return ''
    return 'priority' if 'hungate' in str(dataset or '').lower() else 'secondary'


def _merge_run_rows(primary_rows, manifest_rows):
    rows = [tuple(row) for row in primary_rows]
    seen = {
        (str(row[1] or ''), str(row[2] or ''), str(row[4] or ''), str(row[5] or ''))
        for row in rows
    }
    for row in manifest_rows:
        key = (str(row[1] or ''), str(row[2] or ''), str(row[4] or ''), str(row[5] or ''))
        if key not in seen:
            rows.append(tuple(row))
            seen.add(key)
    return sorted(rows, key=lambda row: (str(row[4] or ''), str(row[0] or '')))


def _merge_dataset_rows(primary_rows, staged_rows, manifest_rows):
    rows = {str(row[0]): tuple(row) for row in primary_rows if str(row[0] or '').strip()}
    for row in staged_rows:
        dataset = str(row[0] or '')
        if dataset and dataset not in rows:
            rows[dataset] = tuple(row)
    for run in manifest_rows:
        dataset = str(run[2] or '')
        if dataset and dataset not in rows:
            role = _infer_dataset_role(dataset, str(run[1] or ''))
            rows[dataset] = (dataset, role, _infer_baseline_tier(dataset, role), 0, 0, 0)
    return sorted(rows.values(), key=lambda row: str(row[0]))


def _merge_status_rows(primary_rows, extra_rows):
    counts = {str(status): int(count or 0) for status, count in primary_rows}
    for status, count in extra_rows:
        key = str(status or 'RECEIVED')
        counts[key] = counts.get(key, 0) + int(count or 0)
    return sorted(counts.items())


def _discover_staged_state(project_root: Path | None, known_datasets: set[str]):
    if project_root is None or not project_root.exists():
        return [], [], []
    dataset_rows = []
    isolate_rows = []
    status_counts = {}
    added = set()
    for path in sorted(project_root.rglob('.*_project.sqlite')):
        try:
            if not path.is_file() or path.stat().st_size == 0:
                continue
            with sqlite3.connect(str(path)) as conn:
                tables = {
                    str(row[0]) for row in conn.execute(
                        "SELECT name FROM sqlite_master WHERE type='table'"
                    ).fetchall()
                }
                if 'sequences' not in tables:
                    continue
                rows = _query_dataset_rows(conn)
                missing = [
                    str(row[0]) for row in rows
                    if str(row[0]) not in known_datasets and str(row[0]) not in added
                ]
                if not missing:
                    continue
                dataset_rows.extend(row for row in rows if str(row[0]) in missing)
                isolate_rows.extend(_query_isolate_rows(conn, missing))
                for status, count in _query_status_rows_for_datasets(conn, missing):
                    status_counts[status] = status_counts.get(status, 0) + int(count or 0)
                added.update(missing)
        except (OSError, sqlite3.Error):
            continue
    return dataset_rows, sorted(status_counts.items()), isolate_rows


def write_annual_report(db, outdir: str | Path) -> dict:
    output = Path(outdir)
    output.mkdir(parents=True, exist_ok=True)
    with db.connect() as conn:
        dataset_rows = _query_dataset_rows(conn)
        marker_qc_rows = _query_dataset_marker_qc_rows(conn)
        status_rows = conn.execute(
            'SELECT status, COUNT(*) FROM isolate_status GROUP BY status ORDER BY status'
        ).fetchall()
        genome_rows = conn.execute(
            'SELECT genome_id, sequence_id, accession, genome_status, genome_qc_pass, '
            'completeness, contamination, gtdb_taxonomy, ani_cluster FROM genome_records '
            'ORDER BY sequence_id, genome_id'
        ).fetchall()
        run_rows = conn.execute(
            'SELECT run_id, workflow, dataset, status, started_at, completed_at, manifest_path '
            'FROM project_runs ORDER BY started_at, run_id'
        ).fetchall()
        round_rows = conn.execute(
            'SELECT round_id, mode, created_at FROM selection_rounds ORDER BY created_at, round_id'
        ).fetchall()
        isolate_rows = _query_isolate_rows(conn)
        removal_rows = conn.execute(
            'SELECT sequence_id, original_dataset, partner_id, sequence_length, '
            'taxonomy, reason, removed_at FROM sequence_removals '
            'ORDER BY removed_at, removal_id'
        ).fetchall()

    project_root = _project_root_from_db(db)
    manifest_rows = _discover_run_manifests(project_root, output)
    qc_summaries = _discover_qc_summaries(project_root, output)
    qc_rows_by_dataset = _aggregate_qc_summaries(qc_summaries)
    visual_report_rows = _discover_visual_report_rows(project_root, output, qc_summaries)
    visual_gallery_rows = _build_visual_gallery_rows(visual_report_rows, output)
    staged_dataset_rows, staged_status_rows, staged_isolate_rows = _discover_staged_state(
        project_root,
        {str(row[0]) for row in dataset_rows},
    )
    run_rows = _merge_run_rows(run_rows, manifest_rows)
    dataset_rows = _merge_dataset_rows(dataset_rows, staged_dataset_rows, manifest_rows)
    status_rows = _merge_status_rows(status_rows, staged_status_rows)
    if staged_isolate_rows:
        isolate_rows = sorted(
            [*isolate_rows, *staged_isolate_rows], key=lambda row: str(row[0] or '')
        )

    active = {
        str(row[0]): {
            'dataset': str(row[1] or ''), 'partner_id': str(row[2] or ''),
            'status': str(row[3] or ''), 'marker_qc': str(row[4] or ''),
        }
        for row in isolate_rows
    }
    dataset_roles = db.get_dataset_roles()
    for dataset, role, baseline_tier, genomes_available, _proposed, _sequences in dataset_rows:
        dataset = str(dataset or '')
        if dataset and dataset not in dataset_roles:
            inferred_role = role if role in {'baseline', 'candidate'} else _infer_dataset_role(dataset)
            dataset_roles[dataset] = {
                'role': inferred_role,
                'genomes_available': bool(genomes_available),
                'baseline_tier': baseline_tier or _infer_baseline_tier(dataset, inferred_role),
            }
    latest_assessments = db.get_latest_assessment_rows()
    candidate_rows = []
    recommendation_counts = {}
    for sequence_id, metadata in active.items():
        if dataset_roles.get(metadata['dataset'], {}).get('role') != 'candidate':
            continue
        assessment = latest_assessments.get(sequence_id, {})
        decision = build_selection_decision(assessment) if assessment else {
            'decision': 'NOT YET ASSESSED', 'evidence_quality': 'NA',
            'decision_reason': 'no stored Performance Review assessment',
        }
        recommendation = decision['decision']
        recommendation_counts[recommendation] = recommendation_counts.get(recommendation, 0) + 1
        candidate_rows.append((
            sequence_id, metadata['dataset'], metadata['partner_id'], metadata['status'],
            metadata['marker_qc'], recommendation, decision.get('evidence_quality', 'NA'),
            assessment.get('sequencing_set_role', 'NA'), assessment.get('taxonomy', 'NA'),
            assessment.get('nearest_hit', 'NA'), assessment.get('nearest_identity', 'NA'),
            assessment.get('project_nearest_hit', 'NA'),
            assessment.get('project_nearest_identity', 'NA'),
            assessment.get('investigation_score', 'NA'),
            decision.get('decision_reason', ''),
        ))

    isolate_path = output / 'project_isolate_ledger.tsv'
    with open(isolate_path, 'w') as handle:
        handle.write('SequenceID\tDataset\tPartnerID\tOperationalStatus\tMarkerQC\tSelectedForGenomeSequencing\tGenomeAlreadySequenced\n')
        for row in isolate_rows:
            handle.write('\t'.join(str(value) for value in row) + '\n')
    genome_path = output / 'project_genome_ledger.tsv'
    with open(genome_path, 'w') as handle:
        handle.write('GenomeID\tSequenceID\tAccession\tGenomeStatus\tGenomeQCPass\tCompleteness\tContamination\tGTDBTaxonomy\tANICluster\n')
        for row in genome_rows:
            handle.write('\t'.join('NA' if value is None else str(value) for value in row) + '\n')
    candidate_path = output / 'project_candidate_overview.tsv'
    with open(candidate_path, 'w') as handle:
        handle.write(
            'SequenceID\tDataset\tPartnerID\tOperationalStatus\tMarkerQC\tRecommendation\t'
            'EvidenceQuality\tSequencingSetRole\tGTDBTaxonomy\tBaselineNearestHit\t'
            'BaselineNearestIdentity\tProjectNearestHit\tProjectNearestIdentity\t'
            'InvestigationScore\tRecommendationReason\n'
        )
        for row in candidate_rows:
            handle.write('\t'.join(str(value) for value in row) + '\n')
    removals_path = output / 'sequence_removal_ledger.tsv'
    with open(removals_path, 'w') as handle:
        handle.write('SequenceID\tOriginalDataset\tPartnerID\tSequenceLength\tTaxonomy\tReason\tRemovedAt\n')
        for row in removal_rows:
            handle.write('\t'.join('NA' if value is None else str(value) for value in row) + '\n')
    changes_path = write_decision_changes(db, output / 'decision_changes.tsv')

    # Explicit marker evidence ledger: multiple marker observations may exist
    # for one isolate, while project_isolate_ledger remains one row per isolate.
    marker_evidence_path = output / 'project_marker_evidence_ledger.tsv'
    try:
        from branchmanager.marker_evidence import write_marker_evidence_report
        write_marker_evidence_report(db, marker_evidence_path)
    except Exception as exc:
        marker_evidence_path.write_text(
            'sequence_id\tdataset\tsource_type\tevidence_tier\tqc_class\tmanual_review_status\tError\n'
            f'\t\t\t\t\t\t{str(exc).replace(chr(9), " ").replace(chr(10), " ")}\n'
        )

    # MALDI observations remain a separate soft-identification ledger; they are
    # never inserted into the marker sequence/evidence registry or VSEARCH pool.
    maldi_evidence_path = output / 'project_maldi_ledger.tsv'
    try:
        maldi_rows = db.get_maldi_identifications() if hasattr(db, 'get_maldi_identifications') else []
        maldi_rows = [dict(row) for row in (maldi_rows or []) if isinstance(row, dict)]
        maldi_fields = sorted({key for row in maldi_rows for key in row})
        if not maldi_fields:
            maldi_fields = ['sequence_id']
        with open(maldi_evidence_path, 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=maldi_fields, delimiter='\t', extrasaction='ignore')
            writer.writeheader()
            writer.writerows(maldi_rows)
    except Exception as exc:
        maldi_evidence_path.write_text(
            'sequence_id\tError\n\t' + str(exc).replace('\t', ' ').replace('\n', ' ') + '\n'
        )

    v1_paths = {}
    v1_summary = []
    try:
        from branchmanager.v1_ledgers import write_v1_ledgers, v1_summary_counts
        v1_paths = write_v1_ledgers(db, output)
        v1_summary = v1_summary_counts(db, v1_paths.pop('concordance_rows', None))
        v1_summary_path = output / 'project_v1_summary.tsv'
        with open(v1_summary_path, 'w') as handle:
            handle.write('Section\tLabel\tCount\n')
            for section, label, count in v1_summary:
                handle.write(f'{section}\t{label}\t{count}\n')
        v1_paths['v1_summary'] = str(v1_summary_path)
    except Exception as exc:
        (output / 'project_v1_ledgers_error.txt').write_text(str(exc))

    dataset_summary_headers = (
        'Dataset', 'Role', 'Baseline tier', 'Genomes available', 'Proposed sequences', 'Sequences',
        'Active marker-QC passed', 'Active marker-QC review', 'Active marker-QC failed',
        'QC-passed sequence files', 'Accepted sequence files', 'Review sequence files',
        'QC-failed sequence files', 'Read files', 'Failed read files', 'Visual report pages',
    )
    dataset_summary_rows = []
    for dataset, role, baseline_tier, genomes_available, proposed, sequences in dataset_rows:
        marker = marker_qc_rows.get(str(dataset), {})
        qc = qc_rows_by_dataset.get(str(dataset), {})
        dataset_summary_rows.append((
            dataset, role, baseline_tier, genomes_available, proposed, sequences,
            marker.get('active_marker_qc_passed', 0),
            marker.get('active_marker_qc_review', 0),
            marker.get('active_marker_qc_failed', 0),
            qc.get('passed_sequence_files', 0),
            qc.get('accepted_sequence_files', 0),
            qc.get('review_sequence_files', 0),
            qc.get('failed_sequence_files', 0),
            qc.get('read_files', 0),
            qc.get('failed_read_files', 0),
            qc.get('visual_pages', 0),
        ))
    dataset_summary_path = output / 'project_dataset_summary.tsv'
    with open(dataset_summary_path, 'w') as handle:
        handle.write('\t'.join(dataset_summary_headers) + '\n')
        for row in dataset_summary_rows:
            handle.write('\t'.join(str(value) for value in row) + '\n')

    visual_report_path = output / 'project_visual_reports.tsv'
    with open(visual_report_path, 'w') as handle:
        handle.write('Dataset\tWorkflow\tReport\tItems\tPath\n')
        for row in visual_report_rows:
            handle.write('\t'.join(str(value) for value in row) + '\n')
    visual_gallery_path = output / 'project_visual_gallery.tsv'
    with open(visual_gallery_path, 'w') as handle:
        handle.write('Dataset\tWorkflow\tReport\tCaption\tImagePath\tLinkPath\n')
        for row in visual_gallery_rows:
            handle.write('\t'.join(str(row.get(key, '')) for key in (
                'dataset', 'workflow', 'report', 'caption', 'image_path', 'link_path',
            )) + '\n')
    readme_path = output / 'README.md'
    readme_path.write_text(_annual_report_readme_markdown())
    readme_html_path = output / 'annual_report_guide.html'
    readme_html_path.write_text(_annual_report_readme_html())

    baseline_count = sum(row[5] for row in dataset_rows if row[1] == 'baseline')
    hungate_count = sum(row[5] for row in dataset_rows if row[1] == 'baseline' and row[2] == 'priority')
    secondary_count = sum(row[5] for row in dataset_rows if row[1] == 'baseline' and row[2] == 'secondary')
    candidate_count = sum(row[5] for row in dataset_rows if row[1] == 'candidate')
    proposed_count = sum(row[4] for row in dataset_rows if row[1] == 'candidate')
    qc_passed_sequence_files = sum(row.get('passed_sequence_files', 0) for row in qc_summaries)
    qc_review_sequence_files = sum(row.get('review_sequence_files', 0) for row in qc_summaries)
    qc_failed_sequence_files = sum(row.get('failed_sequence_files', 0) for row in qc_summaries)
    visual_page_count = sum(row.get('visual_pages', 0) for row in qc_summaries)
    cards = [
        ('Cultured baseline markers', baseline_count),
        ('Hungate baseline markers', hungate_count),
        ('Secondary baseline markers', secondary_count),
        ('Active partner candidates', candidate_count),
        ('Proposed sequences', proposed_count),
        ('QC-passed sequence files', qc_passed_sequence_files),
        ('QC-review sequence files', qc_review_sequence_files),
        ('QC-failed sequence files', qc_failed_sequence_files),
        ('Visual report pages', visual_page_count),
        ('Withdrawn sequences', len(removal_rows)),
        ('QC-passed genomes', sum(bool(row[4]) for row in genome_rows)),
        ('Genome QC failures', sum(row[3] == 'GENOME_QC_FAILED' for row in genome_rows)),
        ('Selection rounds', len(round_rows)),
    ]
    def table(headers, rows, *, raw_columns=None):
        raw_columns = set(raw_columns or ())
        head = ''.join(f'<th>{html.escape(str(item))}</th>' for item in headers)

        def cell(value):
            return '' if value is None else str(value)

        body_rows = []
        for row in rows:
            cells = []
            for header, value in zip(headers, row):
                rendered = cell(value)
                if header in raw_columns:
                    cells.append(f'<td>{rendered}</td>')
                else:
                    cells.append(f'<td>{html.escape(rendered)}</td>')
            body_rows.append('<tr>' + ''.join(cells) + '</tr>')
        return f'<table><thead><tr>{head}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>'

    report = output / 'annual_report.html'
    report.write_text(
        '<!doctype html><html><head><meta charset="utf-8"><title>BranchManager Project Overview</title>'
        '<style>body{font:15px system-ui;margin:32px;color:#17202a;max-width:1200px;line-height:1.35}'
        'h1,h2{color:#123b2d}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:12px}'
        '.card{border:1px solid #cad5cf;padding:14px;border-radius:6px}.card b{font-size:26px;display:block}'
        '.guide{border:1px solid #cad5cf;background:#f6faf8;border-radius:6px;padding:14px 16px;margin:18px 0 24px}'
        '.guide p{margin:6px 0}.muted{color:#596861}.gallery{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:14px;margin:10px 0 28px}'
        '.visual-card{border:1px solid #cad5cf;border-radius:6px;margin:0;background:white;overflow:hidden}.visual-card img{display:block;width:100%;height:230px;object-fit:contain;background:#f8faf9}'
        '.visual-card figcaption{border-top:1px solid #d9e1dd;padding:9px 10px}.visual-card figcaption b{display:block;color:#123b2d}.visual-card figcaption span{display:block;color:#44534d;font-size:13px;margin-top:3px}'
        'table{border-collapse:collapse;width:100%;margin:10px 0 28px}th,td{border-bottom:1px solid #d9e1dd;text-align:left;padding:7px;vertical-align:top}th{background:#eef4f1}</style>'
        '</head><body><h1>Cumulative Project Overview</h1><p>BranchManager point-in-time marker-to-genome project report.</p>'
        '<div class="cards">' + ''.join(f'<div class="card"><b>{value}</b>{html.escape(label)}</div>' for label, value in cards) + '</div>'
        '<section class="guide"><h2>How To Read This Report</h2>'
        '<p>Use the summary cards for project-level counts, the Datasets table for dataset/QC file outcomes, and Current recommendations for the latest selection decisions.</p>'
        '<p>Embedded visual previews are sampled from the PNG outputs; click any preview to open the full-size image. The TSV ledgers remain the auditable source for numerical decisions.</p>'
        '<p><a href="annual_report_guide.html">Open the full guide</a> or read <a href="README.md">README.md</a>.</p></section>'
        '<h2>Datasets</h2>' + table(dataset_summary_headers, dataset_summary_rows) +
        '<h2>Isolate lifecycle</h2>' + table(('Status', 'Count'), status_rows) +
        '<h2>Current recommendations</h2>' + table(('Recommendation', 'Candidates'), sorted(recommendation_counts.items())) +
        '<h2>Visual reports</h2>' + table(
            ('Dataset', 'Workflow', 'Report', 'Items', 'Link'),
            [
                (dataset, workflow, report_name, items, _html_path_link(path_value, output))
                for dataset, workflow, report_name, items, path_value in visual_report_rows
            ],
            raw_columns={'Link'},
        ) +
        '<h2>v1 evidence summary</h2>' + table(('Section', 'Label', 'Count'), v1_summary) +
        '<h2>Embedded visual previews</h2>' + _visual_gallery_html(visual_gallery_rows) +
        '<h2>Workflow runs</h2>' + table(('Run', 'Workflow', 'Dataset', 'Status', 'Started', 'Completed', 'Manifest'), run_rows) +
        '<h2>Selection rounds</h2>' + table(('Round', 'Mode', 'Created'), round_rows) +
        '<h2>Exit Interviews</h2>' + table(('Sequence', 'Dataset', 'Partner', 'Length', 'Taxonomy', 'Reason', 'Removed'), removal_rows) +
        '<p>Detailed ledgers: <a href="project_dataset_summary.tsv">datasets</a>, '
        '<a href="project_isolate_ledger.tsv">isolates</a>, '
        '<a href="project_candidate_overview.tsv">current candidate assessments</a>, '
        '<a href="project_genome_ledger.tsv">genomes</a>, '
        '<a href="project_marker_evidence_ledger.tsv">marker evidence</a>, '
        '<a href="project_maldi_ledger.tsv">MALDI observations</a>, '
        '<a href="project_taxonomy_concordance.tsv">taxonomy concordance</a>, '
        '<a href="project_unreviewed_markers.tsv">unreviewed markers</a>, '
        '<a href="project_visual_reports.tsv">visual reports</a>, '
        '<a href="project_visual_gallery.tsv">visual gallery</a>, '
        '<a href="annual_report_guide.html">report guide</a>, '
        '<a href="sequence_removal_ledger.tsv">Exit Interviews</a>, '
        '<a href="decision_changes.tsv">latest decision changes</a>.</p></body></html>\n'
    )
    return {
        'html': str(report), 'readme': str(readme_path), 'readme_html': str(readme_html_path),
        'datasets': str(dataset_summary_path), 'isolates': str(isolate_path),
        'genomes': str(genome_path), 'candidates': str(candidate_path),
        'marker_evidence': str(marker_evidence_path), 'maldi': str(maldi_evidence_path),
        'visual_reports': str(visual_report_path), 'visual_gallery': str(visual_gallery_path),
        'removals': str(removals_path), 'decision_changes': changes_path,
        **v1_paths,
    }
