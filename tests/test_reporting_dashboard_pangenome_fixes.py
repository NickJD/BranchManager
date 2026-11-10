import json
import re
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
    assert data['records'][0]['taxonomy'][0]['reference_coverage'] == 92.0
    assert data['pangenome'][0]['gap'] == 7.0
    assert data['records'][0]['pangenome_rank'] == 'SPECIES'
    assert data['records'][0]['pangenome_taxon'] == 'Example species'
    assert 'isolateBody' in page
    assert 'exportBtn' in page
    assert 'QUB_0001' in embedded.group(1)


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


def test_dashboard_uses_plain_language_genome_qc_and_collection_overview(tmp_path):
    assessment = tmp_path / 'assessment'
    assessment.mkdir()
    _write_tsv(assessment / 'selection_summary.tsv', [
        {'SequenceID': 'QUB_0001', 'PartnerID': 'QUB',
         'Recommendation': 'ALREADY SEQUENCED', 'MarkerQC': 'PASS_HIGH_CONFIDENCE'},
        {'SequenceID': 'BSI_0001', 'PartnerID': 'BSI',
         'Recommendation': 'ALREADY SELECTED - GENOME PENDING', 'MarkerQC': 'PASS_WITH_WARNINGS'},
    ])
    _write_tsv(assessment / 'genome_collection_assessment.tsv', [{
        'sequence_id': 'AAFC_0001', 'partner_id': 'AAFC_01', 'genome_id': 'GCF_0001',
        'genome_status': 'registered', 'genome_qc_pass': '1',
        'completeness': '98.5', 'contamination': '1.2',
        'gtdb_taxonomy': 'd__Bacteria;p__Bacillota;g__Blautia',
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
    assert data['counts']['already_sequenced'] == 1
    assert data['counts']['selected_pending'] == 1
    assert data['counts']['genome_only_records'] == 1
    assert data['counts']['genome_qc_pass'] == 1
    assert 'Genome-only partner isolates' in page
    assert 'Selected · genome pending' in page
    assert "genomeQCLabel(x)" in page


def test_dashboard_reads_pangenome_eligibility_reason_and_provenance(tmp_path):
    assessment = tmp_path / 'assessment'
    assessment.mkdir()
    _write_tsv(assessment / 'selection_summary.tsv', [{
        'SequenceID': 'QUB_0001', 'PartnerID': 'QUB',
        'Recommendation': 'REVIEW BEFORE SELECTION',
    }])
    _write_tsv(assessment / 'pangenome_hierarchy.tsv', [{
        'PangenomeRank': 'GENUS', 'PangenomeTaxon': 'Blautia',
        'TaxonomyProvenance': 'CLASSIFIER_LINEAGE',
        'CandidateCount': '1', 'CommittedGenomeCount': '3',
        'PangenomeTarget': 'NA', 'PangenomeGap': 'NA',
        'PangenomeEligible': 'No',
        'EligibilityReason': 'partner catalogue rank not supplied',
        'CandidateIDs': 'QUB_0001',
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
    assert data['pangenome'][0]['eligible'] == 'No'
    assert data['pangenome'][0]['reason'] == 'partner catalogue rank not supplied'
    assert data['pangenome'][0]['provenance'] == 'CLASSIFIER_LINEAGE'
    assert data['pangenome'][0]['committed'] == 3
    assert data['pangenome'][0]['target'] is None
    assert data['pangenome'][0]['gap'] is None


def test_dashboard_defines_pager_before_renderers_and_reconciles_selection_eligibility(tmp_path):
    assessment = tmp_path / 'assessment'
    assessment.mkdir()
    _write_tsv(assessment / 'selection_summary.tsv', [{
        'SequenceID': 'BSI_0012', 'PartnerID': 'BSI',
        'Recommendation': 'PRIORITISE - SET PRIMARY',
        'PangenomeRank': 'GENUS', 'PangenomeTaxon': 'Pseudobutyrivibrio',
        'PangenomeEligible': 'Yes', 'PangenomeTarget': '9',
        'PangenomeGap': '8', 'PangenomeCommittedGenomeCount': '1',
        'PangenomeEligibilityReason': 'partner-catalogue genus assignment and rank-matched local tree context available',
    }])
    _write_tsv(assessment / 'pangenome_hierarchy.tsv', [{
        'PangenomeRank': 'GENUS', 'PangenomeTaxon': 'Pseudobutyrivibrio',
        'TaxonomyProvenance': 'CLASSIFIER_INFERRED',
        'CandidateCount': '1', 'CommittedGenomeCount': '1',
        'PangenomeTarget': 'NA', 'PangenomeGap': 'NA',
        'PangenomeEligible': 'No',
        'EligibilityReason': 'classifier inferred lineage is exploratory only',
        'CandidateIDs': 'BSI_0012',
    }])
    output = Path(write_performance_review_dashboard(
        assessment / 'selection_summary.tsv', tmp_path,
        pangenome_hierarchy=assessment / 'pangenome_hierarchy.tsv',
    ))
    page = output.read_text(encoding='utf-8')
    embedded = re.search(
        r'<script id="dashboard-data" type="application/json">(.*?)</script>',
        page, re.S,
    )
    data = json.loads(embedded.group(1))
    group = next(x for x in data['pangenome']
                 if x['rank'] == 'GENUS' and x['taxon'] == 'Pseudobutyrivibrio')
    assert group['eligible'] == 'Yes'
    assert group['target'] == 9
    assert group['gap'] == 8
    assert 'partner-catalogue genus assignment' in group['reason']
    assert 'function pager(id,page,pages,total,onPage)' in page
    assert page.index('function pager(id,page,pages,total,onPage)') < page.index('function renderPangenome()')


def test_dashboard_pangenome_group_opens_member_inventory_and_exposes_gap(tmp_path):
    assessment = tmp_path / 'assessment'
    assessment.mkdir()
    _write_tsv(assessment / 'selection_summary.tsv', [{
        'SequenceID': 'BSI_0097', 'PartnerID': 'BSI',
        'Recommendation': 'PRIORITISE - SET PRIMARY',
        'PangenomeRank': 'FAMILY', 'PangenomeTaxon': 'Selenomonadaceae',
        'PangenomeEligible': 'Yes', 'PangenomeTarget': '9',
        'PangenomeGap': '9', 'PangenomeCommittedGenomeCount': '0',
        'PangenomeEligibilityReason': 'Partner-catalogue family assignment and rank-matched tree context available',
    }])
    _write_tsv(assessment / 'pangenome_hierarchy.tsv', [{
        'PangenomeRank': 'FAMILY', 'PangenomeTaxon': 'Selenomonadaceae',
        'TaxonomyProvenance': 'PARTNER_CATALOGUE',
        'CandidateCount': '1', 'CommittedGenomeCount': '0',
        'PangenomeTarget': '9', 'PangenomeGap': '9', 'PangenomeEligible': 'Yes',
        'EligibilityReason': 'partner-catalogue family assignment and tree context available',
        'CandidateIDs': 'BSI_0097', 'CommittedGenomeIDs': 'HUN459;RUM001',
    }])
    _write_tsv(assessment / 'genome_collection_assessment.tsv', [{
        'sequence_id': 'RUM001', 'partner_id': 'Rumen_Isolates',
        'genome_id': 'GCF_RUM001', 'genome_status': 'registered',
        'genome_qc_pass': '1', 'completeness': '99.0',
        'contamination': '0.5', 'gtdb_taxonomy': 'd__Bacteria;f__Selenomonadaceae',
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
    group = next(x for x in data['pangenome']
                 if x['rank'] == 'FAMILY' and x['taxon'] == 'Selenomonadaceae')
    assert group['eligible'] == 'Yes'
    assert group['target'] == 9
    assert group['gap'] == 9
    assert {m['id'] for m in group['members']} == {'BSI_0097', 'HUN459', 'RUM001'}
    assert next(m for m in group['members'] if m['id'] == 'RUM001')['genome_id'] == 'GCF_RUM001'
    isolate = next(x for x in data['records'] if x['id'] == 'BSI_0097')
    assert isolate['pangenome_rank'] == 'FAMILY'
    assert isolate['pangenome_taxon'] == 'Selenomonadaceae'
    assert isolate['pangenome_gap'] == 9
    assert isolate['pangenome_target'] == 9
    assert 'data-pangenome-open' in page
    assert 'Pangenome group membership' in page


def test_dashboard_pangenome_renderer_uses_embedded_dataset_for_group_indices(tmp_path):
    assessment = tmp_path / 'assessment'
    assessment.mkdir()
    _write_tsv(assessment / 'selection_summary.tsv', [{
        'SequenceID': 'QUB_0001', 'PartnerID': 'QUB',
        'Recommendation': 'REVIEW BEFORE SELECTION',
    }])
    _write_tsv(assessment / 'pangenome_hierarchy.tsv', [{
        'PangenomeRank': 'GENUS', 'PangenomeTaxon': 'Blautia',
        'CandidateCount': '1', 'CommittedGenomeCount': '0',
        'PangenomeTarget': '9', 'PangenomeGap': '9',
        'PangenomeEligible': 'Yes', 'CandidateIDs': 'QUB_0001',
    }])
    output = Path(write_performance_review_dashboard(
        assessment / 'selection_summary.tsv', tmp_path,
        pangenome_hierarchy=assessment / 'pangenome_hierarchy.tsv',
    ))
    page = output.read_text(encoding='utf-8')
    assert 'D.pangenome.indexOf(x)' in page
    assert 'pangenome.indexOf(x)' not in page.replace('D.pangenome.indexOf(x)', '')


def test_dashboard_deduplicates_pangenome_members_and_recalculates_gap(tmp_path):
    assessment = tmp_path / 'assessment'
    assessment.mkdir()
    _write_tsv(assessment / 'selection_summary.tsv', [{
        'SequenceID': 'QUB_0001', 'PartnerID': 'QUB',
        'Recommendation': 'REVIEW BEFORE SELECTION',
    }])
    _write_tsv(assessment / 'pangenome_hierarchy.tsv', [
        {
            'PangenomeRank': 'FAMILY', 'PangenomeTaxon': 'Anaerovoracaceae',
            'TaxonomyProvenance': 'PARTNER_CATALOGUE',
            'CandidateCount': '3', 'CandidateIDs': 'QUB_0001;QUB_0002;QUB_0003',
            'CommittedGenomeCount': '2', 'CommittedGenomeIDs': 'QUB_0901;QUB_0902',
            'PangenomeTarget': '9', 'PangenomeGap': '9', 'PangenomeEligible': 'Yes',
            'EligibilityReason': 'Partner catalogue and rank-matched tree support target',
        },
        {
            'PangenomeRank': 'FAMILY', 'PangenomeTaxon': 'Anaerovoracaceae',
            'TaxonomyProvenance': 'PARTNER_CATALOGUE',
            'CandidateCount': '3', 'CandidateIDs': 'QUB_0002;QUB_0003;QUB_0004',
            'CommittedGenomeCount': '2', 'CommittedGenomeIDs': 'QUB_0901;QUB_0902',
            'PangenomeTarget': '9', 'PangenomeGap': '9', 'PangenomeEligible': 'Yes',
            'EligibilityReason': 'Partner catalogue and rank-matched tree support target',
        },
    ])
    output = Path(write_performance_review_dashboard(
        assessment / 'selection_summary.tsv', tmp_path,
        pangenome_hierarchy=assessment / 'pangenome_hierarchy.tsv',
    ))
    page = output.read_text(encoding='utf-8')
    embedded = re.search(
        r'<script id="dashboard-data" type="application/json">(.*?)</script>',
        page, re.S,
    )
    data = json.loads(embedded.group(1))
    group = next(x for x in data['pangenome'] if x['taxon'] == 'Anaerovoracaceae')

    assert group['candidates'] == 4
    assert group['committed'] == 2
    assert group['target'] == 9
    assert group['gap'] == 7
    assert group['member_count'] == 6
    assert {m['id'] for m in group['members']} == {
        'QUB_0001', 'QUB_0002', 'QUB_0003', 'QUB_0004', 'QUB_0901', 'QUB_0902'
    }


def test_dashboard_pangenome_click_and_keyboard_handlers_use_embedded_data():
    import inspect
    from branchmanager import reporting

    source = inspect.getsource(reporting.write_performance_review_dashboard)
    assert "openPangenomeGroup(group)" in source
    assert "const group=D.pangenome[Number(b.dataset.pangenomeOpen)]" in source
    assert "const group=D.pangenome[Number(e.target.dataset.pangenomeOpen)]" in source
    assert "openPangenomeGroup(pangenome[" not in source


def test_dashboard_recalculates_gap_from_committed_count_and_target(tmp_path):
    assessment = tmp_path / 'assessment'
    assessment.mkdir()
    _write_tsv(assessment / 'selection_summary.tsv', [{
        'SequenceID': 'QUB_0001', 'PartnerID': 'QUB',
        'Recommendation': 'REVIEW BEFORE SELECTION',
    }])
    _write_tsv(assessment / 'pangenome_hierarchy.tsv', [{
        'PangenomeRank': 'FAMILY', 'PangenomeTaxon': 'Anaerovoracaceae',
        'CandidateCount': '1', 'CandidateIDs': 'QUB_0001',
        'CommittedGenomeCount': '2', 'PangenomeTarget': '9',
        'PangenomeGap': '9', 'PangenomeEligible': 'Yes',
    }])
    output = Path(write_performance_review_dashboard(
        assessment / 'selection_summary.tsv', tmp_path,
    ))
    page = output.read_text(encoding='utf-8')
    embedded = re.search(
        r'<script id="dashboard-data" type="application/json">(.*?)</script>',
        page, re.S,
    )
    group = json.loads(embedded.group(1))['pangenome'][0]
    assert group['committed'] == 2
    assert group['target'] == 9
    assert group['gap'] == 7
