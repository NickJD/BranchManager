import json
import re
import sqlite3
from pathlib import Path

from branchmanager.reporting import write_performance_review_dashboard


def _write_tsv(path: Path, rows):
    headers = list(rows[0]) if rows else []
    with path.open('w', encoding='utf-8', newline='') as handle:
        handle.write('\t'.join(headers) + '\n')
        for row in rows:
            handle.write('\t'.join(str(row.get(key, '')) for key in headers) + '\n')


def test_dashboard_embeds_assessment_and_is_interactive(tmp_path):
    assessment = tmp_path / 'assessment'
    assessment.mkdir()
    _write_tsv(assessment / 'selection_summary.tsv', [{
        'SequenceID': 'QUB_0001', 'PartnerID': 'QUB',
        'Recommendation': 'SECONDARY - STRAIN DIVERSITY',
        'MarkerQC': 'PASS_HIGH_CONFIDENCE', 'EvidenceQuality': 'HIGH',
        'GTDBAssignmentStatus': 'CONTESTED', 'GTDBTaxonomy': 'd__Bacteria;g__Example;s__Example species',
        'GTDBClassificationIdentity': '99.5', 'GTDBQueryCoverage': '100',
        'RecommendationReason': 'No close project isolate',
        'LocalTreeFigure': 'neighbourhoods/clade_1.png',
    }])
    _write_tsv(assessment / 'sequence_assessment.tsv', [{
        'ID': 'QUB_0001', 'GTDBTargetCoverage': '92.0',
        'Taxonomy_GG2_2024_09': 'd__Bacteria;g__Example;s__Example species',
        'AssignmentStatus_GG2_2024_09': 'CONFIDENT',
        'Identity_GG2_2024_09': '99.4', 'QueryCoverage_GG2_2024_09': '100',
        'TargetCoverage_GG2_2024_09': '91',
    }])
    _write_tsv(assessment / 'pangenome_hierarchy.tsv', [{
        'PangenomeRank': 'SPECIES', 'PangenomeTaxon': 'Example species',
        'CandidateCount': '1', 'CommittedGenomeCount': '2',
        'PangenomeTarget': '9', 'PangenomeGap': '7', 'PangenomeEligible': 'Yes',
        'EligibilityReason': 'No conflicting domain assignment', 'CandidateIDs': 'QUB_0001',
    }])
    (tmp_path / 'run_manifest.json').write_text(json.dumps({
        'workflow': 'performance_review', 'status': 'COMPLETE',
        'branchmanager_version': '0.1.7',
        'workflow_version': '1.1.0',
        'database_identity': {'schema_version': '6', 'sha256': 'abc123'},
    }))

    output = Path(write_performance_review_dashboard(
        assessment / 'selection_summary.tsv', tmp_path,
    ))
    page = output.read_text(encoding='utf-8')
    embedded = re.search(
        r'<script id="dashboard-data" type="application/json">(.*?)</script>',
        page, re.S,
    )
    data = json.loads(embedded.group(1))

    assert data['counts']['total'] == 1
    assert data['run']['workflow_version'] == '1.1.0'
    assert data['run']['database_schema_version'] == '6'
    assert data['run']['database_sha256'] == 'abc123'
    assert data['records'][0]['taxonomy'][0]['reference_coverage'] == 92.0
    assert data['pangenome'][0]['gap'] == 7.0
    assert data['records'][0]['pangenome_rank'] == 'SPECIES'
    assert data['records'][0]['pangenome_taxon'] == 'Example species'
    assert 'isolateBody' in page
    assert 'exportBtn' in page
    assert 'QUB_0001' in embedded.group(1)




def test_dashboard_can_overlay_current_sab_and_genome_status_from_database(tmp_path):
    assessment = tmp_path / 'assessment'
    assessment.mkdir()
    _write_tsv(assessment / 'selection_summary.tsv', [
        {'SequenceID': 'QUB_0001', 'PartnerID': 'QUB', 'Recommendation': 'PRIORITISE - SET PRIMARY',
         'AlreadySelectedForGenomeSequencing': 'False', 'GenomeAlreadySequenced': 'False'},
        {'SequenceID': 'QUB_0002', 'PartnerID': 'QUB', 'Recommendation': 'PRIORITISE - SET PRIMARY',
         'AlreadySelectedForGenomeSequencing': 'False', 'GenomeAlreadySequenced': 'False'},
        {'SequenceID': 'QUB_0003', 'PartnerID': 'QUB', 'Recommendation': 'PRIORITISE - SET PRIMARY',
         'AlreadySelectedForGenomeSequencing': 'False', 'GenomeAlreadySequenced': 'False'},
    ])
    _write_tsv(assessment / 'sequence_assessment.tsv', [
        {'ID': 'QUB_0001'}, {'ID': 'QUB_0002'}, {'ID': 'QUB_0003'},
    ])
    db_path = tmp_path / 'project.sqlite'
    with sqlite3.connect(db_path) as conn:
        conn.executescript("""
            CREATE TABLE schema_migrations (version INTEGER);
            INSERT INTO schema_migrations VALUES (6);
            CREATE TABLE sequencing_metadata (
                id TEXT PRIMARY KEY, selected_for_sequencing INTEGER,
                selected_for_wgs INTEGER, operational_status TEXT
            );
            CREATE TABLE isolate_status (sequence_id TEXT PRIMARY KEY, status TEXT);
            CREATE TABLE sequences (id TEXT PRIMARY KEY, dataset TEXT);
            CREATE TABLE genome_records (
                sequence_id TEXT, genome_id TEXT, genome_status TEXT,
                genome_qc_pass INTEGER, completeness REAL, contamination REAL,
                gtdb_taxonomy TEXT, genome_path TEXT, notes TEXT
            );
            INSERT INTO sequencing_metadata VALUES ('QUB_0001', 1, 0, 'SAB_APPROVED');
            INSERT INTO sequencing_metadata VALUES ('QUB_0002', 0, 0, 'SAB_DEFERRED');
            INSERT INTO sequencing_metadata VALUES ('QUB_0003', 0, 0, 'GENOME_QC_FAILED');
            INSERT INTO isolate_status VALUES ('QUB_0001', 'SAB_APPROVED');
            INSERT INTO isolate_status VALUES ('QUB_0002', 'SAB_DEFERRED');
            INSERT INTO isolate_status VALUES ('QUB_0003', 'GENOME_QC_FAILED');
            INSERT INTO sequences VALUES ('QUB_0003', 'QUB');
            INSERT INTO genome_records VALUES ('QUB_0003', 'GCA_000003.1', 'GENOME_QC_FAILED', 0, 70, 8, '', '', 'QC failed');
        """)

    output = Path(write_performance_review_dashboard(
        assessment / 'selection_summary.tsv', tmp_path, database=db_path,
    ))
    page = output.read_text(encoding='utf-8')
    embedded = re.search(r'<script id="dashboard-data" type="application/json">(.*?)</script>', page, re.S)
    data = json.loads(embedded.group(1))
    records = {row['id']: row for row in data['records']}
    assert records['QUB_0001']['recommendation'] == 'ALREADY SELECTED - GENOME PENDING'
    assert records['QUB_0001']['sab_decision'] == 'APPROVED'
    assert records['QUB_0002']['sab_decision'] == 'DEFERRED'
    assert records['QUB_0003']['recommendation'] == 'ALREADY SEQUENCED'
    assert records['QUB_0003']['genome_already_sequenced'] is True
    assert data['counts']['selected_pending'] == 1
    assert data['counts']['already_sequenced'] == 1
    assert len(data['genome_records']) == 1
    assert data['run']['database_overlay'] is True
    assert data['run']['database_schema_version'] == '6'


def test_dashboard_uses_classifier_taxonomy_exports_with_taxon_column(tmp_path):
    assessment = tmp_path / 'assessment'
    taxonomy = tmp_path / 'taxonomy'
    assessment.mkdir()
    taxonomy.mkdir()
    _write_tsv(assessment / 'selection_summary.tsv', [{
        'SequenceID': 'QUB_0007', 'PartnerID': 'QUB',
        'Recommendation': 'REVIEW BEFORE SELECTION',
    }])
    _write_tsv(assessment / 'sequence_assessment.tsv', [{
        'ID': 'QUB_0007', 'PartnerTaxon': 'bacteria',
    }])
    _write_tsv(taxonomy / 'taxonomy_NCBI_16S_2026_08.tsv', [{
        'ID': 'QUB_0007',
        'Taxon_NCBI_16S_2026_08': 'd__Bacteria;p__Bacillota;g__Streptococcus;s__Streptococcus equinus',
        'AssignmentStatus_NCBI_16S_2026_08': 'PROVISIONAL',
        'Identity_NCBI_16S_2026_08': '99.8',
        'QueryCoverage_NCBI_16S_2026_08': '100.0',
        'TargetCoverage_NCBI_16S_2026_08': '94.0',
        'AssignedRank_NCBI_16S_2026_08': 'species',
        'EvidenceInterpretation_NCBI_16S_2026_08': 'Provisional species-level consensus',
    }])
    output = Path(write_performance_review_dashboard(
        assessment / 'selection_summary.tsv', tmp_path,
    ))
    page = output.read_text(encoding='utf-8')
    embedded = re.search(
        r'<script id="dashboard-data" type="application/json">(.*?)</script>',
        page, re.S,
    )
    data = json.loads(embedded.group(1))
    ncbi = next(x for x in data['records'][0]['taxonomy'] if x['database'] == 'NCBI 16S 2026-08')
    assert ncbi['status'] == 'PROVISIONAL'
    assert ncbi['taxon'] == 'Species: Streptococcus equinus'
    assert ncbi['identity'] == 99.8
    assert ncbi['query_coverage'] == 100.0


def test_dashboard_handles_missing_optional_evidence(tmp_path):
    assessment = tmp_path / 'assessment'
    assessment.mkdir()
    _write_tsv(assessment / 'selection_summary.tsv', [{
        'SequenceID': 'BSI_0001', 'PartnerID': 'BSI',
        'Recommendation': 'REVIEW BEFORE SELECTION',
    }])
    output = Path(write_performance_review_dashboard(
        assessment / 'selection_summary.tsv', tmp_path,
    ))
    assert output.is_file()
    assert 'Not available' in output.read_text(encoding='utf-8')


def test_dashboard_aggregates_candidate_level_hierarchy_rows(tmp_path):
    assessment = tmp_path / 'assessment'
    assessment.mkdir()
    _write_tsv(assessment / 'selection_summary.tsv', [{
        'SequenceID': 'QUB_0001', 'PartnerID': 'QUB',
        'Recommendation': 'REVIEW BEFORE SELECTION',
    }])
    _write_tsv(assessment / 'pangenome_hierarchy.tsv', [
        {'PangenomeRank': 'DOMAIN', 'PangenomeTaxon': 'Bacteria',
         'CandidateCount': '1', 'CommittedGenomeCount': '1040',
         'PangenomeTarget': '9', 'PangenomeGap': '0', 'PangenomeEligible': 'Yes'},
        {'PangenomeRank': 'DOMAIN', 'PangenomeTaxon': 'Bacteria',
         'CandidateCount': '2', 'CommittedGenomeCount': '1040',
         'PangenomeTarget': '9', 'PangenomeGap': '0', 'PangenomeEligible': 'Yes'},
    ])
    output = Path(write_performance_review_dashboard(
        assessment / 'selection_summary.tsv', tmp_path,
    ))
    page = output.read_text(encoding='utf-8')
    embedded = re.search(
        r'<script id="dashboard-data" type="application/json">(.*?)</script>',
        page, re.S,
    )
    data = json.loads(embedded.group(1))
    assert len(data['pangenome']) == 1
    assert data['pangenome'][0]['candidates'] == 3
    assert data['pangenome'][0]['committed'] == 1040
    assert data['pangenome'][0]['source_rows'] == 2
    assert data['pangenome_source_rows'] == 2
    assert 'class="nav"' in page
    assert 'id="toTop"' in page
    assert 'id="pageSize"' in page
    assert 'isolatePager' in page


def test_dashboard_includes_novelty_comparisons_and_keeps_pangenome_panel(tmp_path):
    assessment = tmp_path / 'assessment'
    assessment.mkdir()
    _write_tsv(assessment / 'selection_summary.tsv', [{
        'SequenceID': 'QUB_0007', 'PartnerID': 'QUB',
        'Recommendation': 'SECONDARY - STRAIN DIVERSITY',
    }])
    _write_tsv(assessment / 'novelty_metrics.tsv', [{
        'ID': 'QUB_0007',
        'HungateNearestIdentity': '86.0', 'HungateNearestHit': 'HUN426',
        'HungateMatchesGE97': '0', 'HungateNearestQueryCoverage': '98.9',
        'SecondaryBaselineNearestIdentity': '86.3',
        'SecondaryBaselineNearestHit': 'H1000-1',
        'CulturedProjectNearestIdentity': '86.3',
        'CulturedProjectNearestHit': 'H1000-1',
        'ProjectNearestIdentity': '83.3', 'ProjectNearestHit': 'QUB_0168',
        'ReferenceNearestIdentity': '99.9', 'ReferenceNearestHit': 'RS_GCF_1',
        'ReferenceMatchesGE97': '2', 'NearestGenomeHit': 'None',
        'GenomeCollectionMatchesGE97': '0',
    }])
    _write_tsv(assessment / 'pangenome_hierarchy.tsv', [{
        'PangenomeRank': 'DOMAIN', 'PangenomeTaxon': 'Bacteria',
        'CandidateCount': '1', 'CommittedGenomeCount': '1040',
        'PangenomeTarget': '9', 'PangenomeGap': '0', 'PangenomeEligible': 'Yes',
    }])
    output = Path(write_performance_review_dashboard(
        assessment / 'selection_summary.tsv', tmp_path,
    ))
    page = output.read_text(encoding='utf-8')
    embedded = re.search(
        r'<script id="dashboard-data" type="application/json">(.*?)</script>',
        page, re.S,
    )
    data = json.loads(embedded.group(1))
    comparisons = data['records'][0]['novelty']
    assert len(comparisons) == 5
    assert comparisons[0]['nearest_identity'] == 86.0
    assert comparisons[-1]['nearest_identity'] == 99.9
    assert 'id="novelty"' in page
    assert 'noveltyPartners' in page
    assert "ranks[0]?.toUpperCase()==='DOMAIN'" not in page
    assert 'id="coverage"' in page
    assert 'id="pangenomeBody"' in page


def test_genome_collection_export_includes_database_records(tmp_path):
    from branchmanager.db.interface import Database
    from branchmanager.reporting import write_genome_collection_assessment

    db = Database(str(tmp_path / 'project.sqlite'))
    db.initialise()
    with db.connect() as conn:
        conn.execute("INSERT INTO sequences (id, dataset) VALUES ('QUB_0001', 'QUB')")
    db.upsert_genome_records([{
        'genome_id': 'GCA_000001.1', 'sequence_id': 'QUB_0001',
        'genome_status': 'QC_PASSED', 'genome_qc_pass': True,
        'completeness': 99.1, 'contamination': 0.4,
        'gtdb_taxonomy': 'd__Bacteria;p__Bacillota;g__Example',
        'genome_path': '/genomes/GCA_000001.1.fna', 'notes': 'test record',
    }])
    output = tmp_path / 'assessment' / 'genome_collection_assessment.tsv'
    write_genome_collection_assessment(db, output)
    import csv
    with output.open(newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert len(rows) == 1
    assert rows[0]['sequence_id'] == 'QUB_0001'
    assert rows[0]['partner_id'] == 'QUB'
    assert rows[0]['genome_qc_pass'] == '1'
    assert rows[0]['completeness'] == '99.1'


def test_dashboard_joins_novelty_metrics_using_sequence_id_header(tmp_path):
    assessment = tmp_path / 'assessment'
    assessment.mkdir()
    _write_tsv(assessment / 'selection_summary.tsv', [{
        'SequenceID': 'QUB_0099', 'PartnerID': 'QUB',
        'Recommendation': 'REVIEW BEFORE SELECTION',
    }])
    _write_tsv(assessment / 'novelty_metrics.tsv', [{
        'SequenceID': 'QUB_0099', 'CulturedProjectNearestHit': 'HUN_001',
        'CulturedProjectNearestIdentity': '96.7',
        'CulturedProjectMatchesGE97': '0',
        'CulturedProjectMatchesGE95': '2',
    }])
    output = Path(write_performance_review_dashboard(
        assessment / 'selection_summary.tsv', tmp_path,
    ))
    page = output.read_text(encoding='utf-8')
    embedded = re.search(
        r'<script id="dashboard-data" type="application/json">(.*?)</script>',
        page, re.S,
    )
    data = json.loads(embedded.group(1))
    novelty = next(x for x in data['records'][0]['novelty'] if x['key'] == 'CulturedProject')
    assert novelty['nearest_hit'] == 'HUN_001'
    assert novelty['nearest_identity'] == 96.7
    assert novelty['matches_97'] == 0



def test_dashboard_uses_explicit_cultured_project_identity_when_legacy_alias_is_missing(tmp_path):
    assessment = tmp_path / 'assessment'
    assessment.mkdir()
    _write_tsv(assessment / 'selection_summary.tsv', [{
        'SequenceID': 'BSI_0031',
        'PartnerID': 'BSI',
        'Recommendation': 'PRIORITISE - SET PRIMARY',
        'MarkerQC': 'PASS_HIGH_CONFIDENCE',
        'BaselineEvidenceClass': '100.00',
        'BaselineNearestIdentity': 'NA',
        'BaselineNearestHit': 'ELIGIBLE',
        'CulturedProjectNearestIdentity': '98.72',
        'CulturedProjectNearestHit': 'HUN106',
        'ProjectNearestIdentity': '44.8',
        'ProjectNearestHit': 'BSI_0073',
    }])
    _write_tsv(assessment / 'sequence_assessment.tsv', [{
        'ID': 'BSI_0031',
        'project_nearest_identity': '44.8',
    }])

    output = Path(write_performance_review_dashboard(
        assessment / 'selection_summary.tsv', tmp_path,
    ))
    page = output.read_text(encoding='utf-8')
    embedded = re.search(
        r'<script id="dashboard-data" type="application/json">(.*?)</script>',
        page, re.S,
    )
    record = json.loads(embedded.group(1))['records'][0]

    assert record['baseline_nearest_identity'] == 98.72
    assert record['baseline_nearest_hit'] == 'HUN106'
    assert record['project_nearest_identity'] == 44.8
    assert record['project_nearest_hit'] == 'BSI_0073'
    assert record['project_identity_warning'] is False


def test_dashboard_flags_project_identity_when_marker_qc_has_warnings_or_fails(tmp_path):
    assessment = tmp_path / 'assessment'
    assessment.mkdir()
    _write_tsv(assessment / 'selection_summary.tsv', [
        {
            'SequenceID': 'BSI_0015',
            'PartnerID': 'BSI',
            'Recommendation': 'REVIEW BEFORE SELECTION',
            'MarkerQC': 'PASS_WITH_WARNINGS',
            'CulturedProjectNearestIdentity': '97.4',
            'ProjectNearestIdentity': '0.0',
        },
        {
            'SequenceID': 'BSI_0016',
            'PartnerID': 'BSI',
            'Recommendation': 'EXCLUDE',
            'MarkerQC': 'FAIL_QC',
            'CulturedProjectNearestIdentity': '96.2',
            'ProjectNearestIdentity': '88.1',
        },
    ])

    output = Path(write_performance_review_dashboard(
        assessment / 'selection_summary.tsv', tmp_path,
    ))
    page = output.read_text(encoding='utf-8')
    embedded = re.search(
        r'<script id="dashboard-data" type="application/json">(.*?)</script>',
        page, re.S,
    )
    records = {row['id']: row for row in json.loads(embedded.group(1))['records']}

    assert records['BSI_0015']['project_nearest_identity'] == 0.0
    assert records['BSI_0015']['project_identity_warning'] is True
    assert 'PASS_WITH_WARNINGS' in records['BSI_0015']['project_identity_warning_reason']
    assert records['BSI_0016']['project_identity_warning'] is True
    assert 'FAIL_QC' in records['BSI_0016']['project_identity_warning_reason']
    assert 'qc-warning' in page
    assert 'qc-warning-text' in page



def test_dashboard_displays_nearest_cultured_and_project_hits_from_novelty_metrics(tmp_path):
    assessment = tmp_path / 'assessment'
    assessment.mkdir()
    _write_tsv(assessment / 'selection_summary.tsv', [{
        'SequenceID': 'BSI_0072',
        'PartnerID': 'BSI',
        'Recommendation': 'PRIORITISE - SET PRIMARY',
        'MarkerQC': 'PASS_HIGH_CONFIDENCE',
        # Deliberately stale legacy values: novelty_metrics.tsv is the
        # source containing the actual nearest-hit identifiers.
        'BaselineNearestIdentity': 'NA',
        'BaselineNearestHit': 'ELIGIBLE',
        'ProjectNearestIdentity': '0',
        'ProjectNearestHit': 'BSI_0065',
    }])
    _write_tsv(assessment / 'novelty_metrics.tsv', [{
        'SequenceID': 'BSI_0072',
        'HungateNearestIdentity': '99.4',
        'HungateNearestHit': 'HUN040',
        'HungateNearestQueryCoverage': '100',
        'HungateNearestAlignmentLength': '825',
        'SecondaryBaselineNearestIdentity': '98.3',
        'SecondaryBaselineNearestHit': 'H1000-P-2558860977',
        'SecondaryBaselineNearestQueryCoverage': '100',
        'SecondaryBaselineNearestAlignmentLength': '825',
        'CulturedProjectNearestIdentity': '99.4',
        'CulturedProjectNearestHit': 'HUN040',
        'CulturedProjectNearestQueryCoverage': '100',
        'CulturedProjectNearestAlignmentLength': '825',
        'ProjectNearestIdentity': '100',
        'ProjectNearestHit': 'BSI_0065',
        'ProjectNearestQueryCoverage': '100',
        'ProjectNearestAlignmentLength': '825',
        'ReferenceNearestIdentity': '99.3',
        'ReferenceNearestHit': 'RS_GCF_900142515.1',
    }])

    output = Path(write_performance_review_dashboard(
        assessment / 'selection_summary.tsv', tmp_path,
        novelty_metrics=assessment / 'novelty_metrics.tsv',
    ))
    html = output.read_text(encoding='utf-8')
    embedded = re.search(
        r'<script id="dashboard-data" type="application/json">(.*?)</script>',
        html, re.S,
    )
    record = json.loads(embedded.group(1))['records'][0]

    assert record['baseline_nearest_identity'] == 99.4
    assert record['baseline_nearest_hit'] == 'HUN040'
    assert record['project_nearest_identity'] == 100.0
    assert record['project_nearest_hit'] == 'BSI_0065'
    assert record['cultured_match_details'][0]['hit'] == 'HUN040'
    assert record['cultured_match_details'][1]['hit'] == 'H1000-P-2558860977'
    assert record['project_match_details']['hit'] == 'BSI_0065'
    assert 'Nearest hit:' in html
    assert 'Closest cultured and project sequence matches' in html
    assert 'H1000-P-2558860977' in html
    assert 'BSI_0065' in html
