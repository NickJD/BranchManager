"""The partner catalogue controls pathway routing; GTDB alone controls pangenome lineage."""
from branchmanager.pipeline.selection_sets import _pangenome_lineage, _partner_rank_taxon


def test_domain_only_partner_taxon_does_not_supply_pangenome_lineage():
    row = {
        'partner_taxon': 'Bacteria',
        'partner_taxon_rank': 'domain',
        'gtdb_taxonomy': 'd__Bacteria;p__Bacillota;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia;s__Roseburia intestinalis',
        'gtdb_assignment_status': 'CONFIDENT',
        'taxonomy': 'd__Bacteria;p__Firmicutes;g__WrongGenus;s__WrongGenus species',
    }
    assert _pangenome_lineage(row)['g'] == 'Roseburia'
    assert _partner_rank_taxon(row) == ('s', 'species', 'Roseburia intestinalis')


def test_informative_partner_taxonomy_does_not_override_gtdb_lineage():
    row = {
        'partner_taxon': 'g__Blautia',
        'partner_taxon_rank': 'genus',
        'gtdb_taxonomy': 'd__Bacteria;p__Bacillota;g__Roseburia;s__Roseburia intestinalis',
        'gtdb_assignment_status': 'CONFIDENT',
        'taxonomy': 'd__Bacteria;p__Bacillota;g__OtherGenus;s__OtherGenus species',
    }
    assert _partner_rank_taxon(row) == ('s', 'species', 'Roseburia intestinalis')


def test_unresolved_gtdb_does_not_fall_back_to_partner_or_marker_taxonomy():
    row = {
        'partner_taxon': 'g__Blautia',
        'partner_taxon_rank': 'genus',
        'gtdb_taxonomy': 'd__Bacteria;p__Bacillota;g__Roseburia',
        'gtdb_assignment_status': 'UNRESOLVED',
        'taxonomy': 'd__Bacteria;p__Bacillota;g__MarkerGenus;s__MarkerGenus species',
    }
    assert _partner_rank_taxon(row) == ('', '', '')
    assert _pangenome_lineage(row) == {}
