from branchmanager.pipeline.selection_sets import (
    _group_for_row, build_sequencing_sets, write_pangenome_hierarchy_tsv,
)


def _row(sid, partner_taxon='Streptococcus equinus', taxonomy=None, **extra):
    return {
        'id': sid,
        'partner_taxon': partner_taxon,
        'partner_taxonomy_source': 'Partner_Isolate_Catalogue',
        'taxonomy': taxonomy or 'd__Bacteria;p__Firmicutes;g__Streptococcus;s__Streptococcus equinus',
        'classification_assignment_status': 'CONTESTED',
        'local_neighbourhood_figure': 'clade_001.png',
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


def test_partner_catalogue_taxon_drives_group_even_when_classifier_is_contested():
    row = _row('QUB_1')
    key, basis, label = _group_for_row(row)
    assert key.startswith('partner:s:')
    assert 'partner-catalogue species' in basis
    assert label == 'Streptococcus equinus'


def test_partner_genus_assignment_is_grouped_at_genus_level():
    row = _row('QUB_1', partner_taxon='Streptococcus')
    key, basis, label = _group_for_row(row)
    assert key.startswith('partner:g:')
    assert 'genus' in basis
    assert label == 'Streptococcus'


def test_classifier_disagreement_does_not_change_catalogue_taxon_or_block_tree_valid_group(tmp_path):
    rows = [
        _row('QUB_1', taxonomy='d__Bacteria;p__Firmicutes;g__Lactobacillus;s__Lactobacillus acidophilus'),
        _row('QUB_2', taxonomy='d__Bacteria;p__Firmicutes;g__Streptococcus;s__Streptococcus equinus'),
    ]
    tree = tmp_path / 'tree.nwk'
    tree.write_text('(QUB_1:0.1,QUB_2:0.2);')
    out = tmp_path / 'sets.tsv'
    build_sequencing_sets(rows, out, tree_path=tree, db=None)
    data = _read_tsv(out)
    assert all(row['PangenomeEligible'] == 'Yes' for row in data)
    assert all(row['PangenomeTaxon'] == 'Streptococcus equinus' for row in data)
    assert data[0]['PartnerTaxon'] == 'Streptococcus equinus'
    assert any(row['ClassifierCatalogueDisagreement'] == 'Yes' for row in data)


def test_pangenome_gap_is_not_reported_without_tree_context(tmp_path):
    row = _row('QUB_1')
    out = tmp_path / 'sets.tsv'
    build_sequencing_sets([row], out, tree_path=None, db=None)
    data = _read_tsv(out)
    assert data[0]['PangenomeEligible'] == 'No'
    assert data[0]['PangenomeTarget'] == 'NA'
    assert data[0]['PangenomeGap'] == 'NA'


def test_tree_mrca_with_another_catalogue_species_blocks_pangenome(tmp_path):
    rows = [_row('QUB_1'), _row('QUB_2')]
    rows.append(_row(
        'QUB_3', partner_taxon='Lactobacillus acidophilus',
        taxonomy='d__Bacteria;p__Firmicutes;g__Lactobacillus;s__Lactobacillus acidophilus',
        local_neighbourhood_figure='clade_002.png',
    ))
    tree = tmp_path / 'tree.nwk'
    tree.write_text('((QUB_1:0.1,QUB_3:0.1):0.1,QUB_2:0.1);')
    out = tmp_path / 'sets.tsv'
    build_sequencing_sets(rows, out, tree_path=tree, db=None)
    data = _read_tsv(out)
    candidates = [row for row in data if row['CandidateID'] in {'QUB_1', 'QUB_2'}]
    assert candidates
    assert all(row['PangenomeEligible'] == 'No' for row in candidates)
    assert all('other assigned s taxa' in row['PangenomeEligibilityReason'] for row in candidates)


def test_hierarchy_reports_species_and_genus_levels_from_catalogue_lineage(tmp_path):
    rows = [
        _row('QUB_1', partner_taxon='d__Bacteria;p__Firmicutes;g__Streptococcus;s__Streptococcus equinus'),
        _row('QUB_2', partner_taxon='d__Bacteria;p__Firmicutes;g__Streptococcus;s__Streptococcus equinus'),
    ]
    tree = tmp_path / 'tree.nwk'
    tree.write_text('(QUB_1:0.1,QUB_2:0.2);')
    out = tmp_path / 'hierarchy.tsv'
    write_pangenome_hierarchy_tsv(rows, out, tree_path=tree)
    data = _read_tsv(out)
    ranks = {row['PangenomeRank'] for row in data}
    assert {'SPECIES', 'GENUS', 'PHYLUM', 'DOMAIN'}.issubset(ranks)
    assert all(row['PangenomeTaxon'] != 'NA' for row in data)


def test_genus_level_pangenome_uses_catalogue_even_when_16s_assigns_another_genus(tmp_path):
    rows = [
        _row('QUB_1', partner_taxon='Streptococcus', taxonomy='d__Bacteria;p__Firmicutes;g__Lactobacillus;s__Lactobacillus acidophilus'),
        _row('QUB_2', partner_taxon='Streptococcus', taxonomy='d__Bacteria;p__Firmicutes;g__Enterococcus;s__Enterococcus faecalis'),
    ]
    tree = tmp_path / 'tree.nwk'
    tree.write_text('(QUB_1:0.1,QUB_2:0.2);')
    out = tmp_path / 'sets.tsv'
    build_sequencing_sets(rows, out, tree_path=tree, db=None)
    data = _read_tsv(out)
    assert all(row['PangenomeEligible'] == 'Yes' for row in data)
    assert all(row['PangenomeRank'] == 'GENUS' for row in data)
    assert all(row['PangenomeTaxon'] == 'Streptococcus' for row in data)
    assert all(row['ClassifierCatalogueDisagreement'] == 'Yes' for row in data)


def test_confident_baseline_tree_leaf_with_conflicting_species_blocks_clade(tmp_path):
    from branchmanager.pipeline.selection_sets import _tree_taxon_clade_agrees

    members = [_row('QUB_1'), _row('QUB_2')]
    baseline_leaf = {
        'id': 'HUN_1',
        'taxonomy': 'd__Bacteria;p__Firmicutes;g__Lactobacillus;s__Lactobacillus acidophilus',
        '_tree_taxonomy_source': 'CONFIDENT_BASELINE_CLASSIFIER',
    }
    tree = tmp_path / 'tree.nwk'
    tree.write_text('((QUB_1:0.1,HUN_1:0.1):0.1,QUB_2:0.1);')
    ok, reason = _tree_taxon_clade_agrees(
        tree, members, members + [baseline_leaf], 's', 'Streptococcus equinus',
    )
    assert ok is False
    assert 'Lactobacillus acidophilus' in reason


def test_committed_rank_index_uses_catalogue_for_partners_and_confident_classifier_for_baselines():
    import sqlite3
    from branchmanager.pipeline.selection_sets import _committed_ids_by_rank

    connection = sqlite3.connect(':memory:')
    connection.executescript('''
        CREATE TABLE sequences (id TEXT, dataset TEXT);
        CREATE TABLE taxonomy (id TEXT, taxonomy TEXT);
        CREATE TABLE sequencing_metadata (id TEXT, selected_for_wgs INTEGER, selected_for_sequencing INTEGER);
        INSERT INTO sequences VALUES ('HUN_BAD','Hungate'),('HUN_GOOD','Hungate'),('QUB_1','Partner');
        INSERT INTO taxonomy VALUES
            ('HUN_BAD','d__Bacteria;p__Firmicutes;g__Streptococcus;s__Streptococcus equinus'),
            ('HUN_GOOD','d__Bacteria;p__Firmicutes;g__Streptococcus;s__Streptococcus equinus'),
            ('QUB_1','d__Bacteria;p__Firmicutes;g__Lactobacillus;s__Lactobacillus acidophilus');
        INSERT INTO sequencing_metadata VALUES ('HUN_BAD',0,0),('HUN_GOOD',0,0),('QUB_1',1,0);
    ''')

    class FakeDB:
        def connect(self):
            return connection
        def get_dataset_roles(self):
            return {'Hungate': {'role': 'baseline', 'genomes_available': True}, 'Partner': {'role': 'candidate'}}
        def get_sequencing_metadata_for_ids(self, ids):
            return {'QUB_1': {'project_metadata': {'partner_taxon_raw': 'Streptococcus equinus'}}}
        def get_classification_evidence(self, ids):
            return {
                ('HUN_BAD', 'GTDB'): {'assignment_status': 'PROVISIONAL'},
                ('HUN_GOOD', 'GTDB'): {'assignment_status': 'CONFIDENT'},
            }
        def get_alt_ref_dbs(self):
            return []

    indexed = _committed_ids_by_rank(FakeDB())
    species_key = ('s', 'streptococcusequinus')
    genus_key = ('g', 'streptococcus')
    assert indexed[species_key] == ['HUN_GOOD', 'QUB_1']
    assert indexed[genus_key] == ['HUN_GOOD', 'QUB_1']
    assert 'HUN_BAD' not in indexed[species_key]


def test_confident_genus_only_classifier_is_not_misreported_as_discordant():
    from branchmanager.pipeline.selection_sets import _classifier_database_concordance_status

    row = {
        'classification_assignment_status': 'CONFIDENT',
        'taxonomy': 'd__Bacteria;p__Firmicutes;g__Streptococcus',
        'alt_tax_GG2': 'd__Bacteria;p__Firmicutes;g__Streptococcus',
        'alt_assignment_status_GG2': 'CONFIDENT',
    }
    assert _classifier_database_concordance_status(row) == 'NOT_COMPARABLE_AT_SPECIES_RANK'


def test_pangenome_lineage_uses_classifier_when_partner_only_declares_domain():
    from branchmanager.pipeline.selection_sets import _pangenome_lineage

    row = {
        "partner_taxon": "bacteria",
        "taxonomy": "d__Bacteria;p__Bacillota;c__Clostridia;o__Lachnospirales;"
                    "f__Lachnospiraceae;g__Roseburia;s__Roseburia bovis",
        "classification_assignment_status": "PROVISIONAL",
    }
    lineage = _pangenome_lineage(row)
    assert lineage.get("g") == "Roseburia"
    assert lineage.get("s") == "Roseburia bovis"


def test_pangenome_lineage_preserves_informative_partner_taxonomy():
    from branchmanager.pipeline.selection_sets import _pangenome_lineage

    row = {
        "partner_taxon": "Roseburia",
        "partner_taxon_rank": "genus",
        "taxonomy": "d__Bacteria;p__Bacillota;c__Clostridia;o__Lachnospirales;"
                    "f__Lachnospiraceae;g__OtherGenus;s__OtherGenus species",
        "classification_assignment_status": "CONFIDENT",
    }
    lineage = _pangenome_lineage(row)
    assert lineage.get("g") == "Roseburia"


def test_hierarchy_uses_taxon_wide_committed_count_across_local_clades(tmp_path):
    """A taxon-level target must not change with the candidate's tree partition."""
    rows = [
        _row('QUB_1', local_neighbourhood_figure='clade_001.png', already_sequenced='true'),
        _row('QUB_2', local_neighbourhood_figure='clade_002.png', already_sequenced='false'),
    ]
    tree = tmp_path / 'tree.nwk'
    tree.write_text('(QUB_1:0.1,QUB_2:0.2);')
    out = tmp_path / 'hierarchy.tsv'
    write_pangenome_hierarchy_tsv(rows, out, tree_path=tree)
    data = _read_tsv(out)
    species_rows = [row for row in data if row['PangenomeRank'] == 'SPECIES']
    assert len(species_rows) == 2
    assert {row['CommittedGenomeCount'] for row in species_rows} == {'1'}
    assert {row['PangenomeGap'] for row in species_rows} == {'8'}


def test_regenerated_assessment_uses_database_specific_taxonomy_status():
    from branchmanager.pipeline.selection_sets import _assignment_status, _pangenome_lineage

    row = {
        'taxonomy_assignment_status': 'UNRESOLVED',
        'gtdb_assignment_status': 'PROVISIONAL',
        'taxonomy': 'd__Bacteria;p__Bacillota;g__Streptococcus;s__Streptococcus equinus',
    }
    assert _assignment_status(row) == 'PROVISIONAL'
    assert _pangenome_lineage(row)['s'] == 'Streptococcus equinus'
