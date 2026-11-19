"""Regression tests for the current GTDB-authoritative pangenome contract."""
from branchmanager.pipeline.selection_sets import (
    _group_for_row, _pangenome_lineage, _assignment_status,
    build_sequencing_sets, write_pangenome_hierarchy_tsv,
)



def _row(sid, partner_taxon='bacteria', gtdb_taxonomy=None, taxonomy=None, **extra):
    gtdb_taxonomy = gtdb_taxonomy or 'd__Bacteria;p__Bacillota;g__Roseburia;s__Roseburia intestinalis'
    return {
        'id': sid,
        'partner_taxon': partner_taxon,
        'partner_taxonomy_source': 'Partner_Isolate_Catalogue',
        'taxonomy': taxonomy or 'd__Bacteria;p__Firmicutes;g__Lactobacillus;s__Lactobacillus acidophilus',
        'gtdb_taxonomy': gtdb_taxonomy,
        'gtdb_assignment_status': 'CONFIDENT',
        'classification_assignment_status': 'CONTESTED',
        'local_neighbourhood_figure': 'clade_001',
        'in_tree': 'Yes',
        'evidence_quality': 'HIGH',
        'nearest_identity': '90',
        'project_nearest_identity': '90',
        'reference_nearest_identity': '90',
        **extra,
    }


def _read_tsv(path):
    lines = path.read_text().splitlines()
    return [dict(zip(lines[0].split('\t'), line.split('\t'))) for line in lines[1:]]


def test_pangenome_grouping_uses_gtdb_not_partner_catalogue_taxon():
    row = _row('QUB_1', partner_taxon='Bacteria',
               gtdb_taxonomy='d__Bacteria;p__Bacillota;g__Roseburia;s__Roseburia intestinalis')
    key, basis, label = _group_for_row(row)
    assert key.startswith('gtdb:s:')
    assert 'GTDB species' in basis
    assert label == 'Roseburia intestinalis'
    assert 'partner' not in key.casefold()


def test_classifier_taxonomy_does_not_override_gtdb_grouping():
    row = _row('QUB_1', partner_taxon='Escherichia coli',
               gtdb_taxonomy='d__Bacteria;p__Bacillota;g__Roseburia;s__Roseburia intestinalis',
               taxonomy='d__Bacteria;p__Proteobacteria;g__Escherichia;s__Escherichia coli')
    key, _basis, label = _group_for_row(row)
    assert key.startswith('gtdb:s:')
    assert label == 'Roseburia intestinalis'
    assert _pangenome_lineage(row)['s'] == 'Roseburia intestinalis'


def test_missing_or_unresolved_gtdb_is_not_filled_from_marker_or_partner_taxonomy():
    row = _row('QUB_1', partner_taxon='Escherichia coli',
               gtdb_taxonomy='', gtdb_assignment_status='UNRESOLVED')
    key, basis, label = _group_for_row(row)
    assert key.startswith('unresolved:')
    assert 'unresolved' in basis.casefold()
    assert label == 'UNASSIGNED [GTDB_TAXONOMY_UNRESOLVED]'
    assert _pangenome_lineage(row) == {}


def test_unusable_gtdb_status_rejects_lineage_even_if_taxonomy_is_present():
    for status in ('UNRESOLVED', 'LOW_SIMILARITY', 'NOT_CONFIDENT', 'FAILED'):
        row = _row('QUB_1', gtdb_assignment_status=status)
        assert _pangenome_lineage(row) == {}
        assert _group_for_row(row)[0].startswith('unresolved:')


def test_pangenome_gap_is_not_reported_without_tree_context(tmp_path):
    row = _row('QUB_1')
    out = tmp_path / 'sets.tsv'
    build_sequencing_sets([row], out, tree_path=None, db=None)
    data = _read_tsv(out)
    assert data[0]['PangenomeEligible'] == 'No'
    assert data[0]['PangenomeTarget'] == 'NA'
    assert data[0]['PangenomeGap'] == 'NA'


def test_tree_conflict_check_uses_gtdb_lineage_of_other_leaf(tmp_path):
    from branchmanager.pipeline.selection_sets import _tree_taxon_clade_agrees
    members = [_row('QUB_1'), _row('QUB_2')]
    other = _row('HUN_1', gtdb_taxonomy='d__Bacteria;p__Bacillota;g__Lactobacillus;s__Lactobacillus acidophilus')
    tree = tmp_path / 'tree.nwk'
    tree.write_text('((QUB_1:0.1,HUN_1:0.1):0.1,QUB_2:0.1);')
    ok, reason = _tree_taxon_clade_agrees(tree, members, members + [other], 's', 'Roseburia intestinalis')
    assert ok is False
    assert 'Lactobacillus acidophilus' in reason
    assert 'GTDB' in reason


def test_hierarchy_reports_ranks_from_gtdb_not_catalogue_taxon(tmp_path):
    rows = [
        _row('QUB_1', partner_taxon='bacteria', gtdb_taxonomy='d__Bacteria;p__Bacillota;g__Roseburia;s__Roseburia intestinalis'),
        _row('QUB_2', partner_taxon='archaea', gtdb_taxonomy='d__Bacteria;p__Bacillota;g__Roseburia;s__Roseburia intestinalis'),
    ]
    tree = tmp_path / 'tree.nwk'
    tree.write_text('(QUB_1:0.1,QUB_2:0.2);')
    out = tmp_path / 'hierarchy.tsv'
    write_pangenome_hierarchy_tsv(rows, out, tree_path=tree)
    data = _read_tsv(out)
    ranks = {row['PangenomeRank'] for row in data}
    assert {'SPECIES', 'GENUS', 'PHYLUM', 'DOMAIN'}.issubset(ranks)
    assert all(row['TaxonomyProvenance'] == 'GTDB' for row in data)
    species = [row for row in data if row['PangenomeRank'] == 'SPECIES']
    assert len(species) == 1 and species[0]['PangenomeTaxon'] == 'Roseburia intestinalis'


def test_classifier_status_prefers_database_specific_gtdb_status():
    row = {'taxonomy_assignment_status': 'UNRESOLVED', 'gtdb_assignment_status': 'PROVISIONAL'}
    assert _assignment_status(row) == 'PROVISIONAL'
