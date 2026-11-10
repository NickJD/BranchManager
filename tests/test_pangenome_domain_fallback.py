from branchmanager.pipeline.selection_sets import _partner_rank_taxon, _pangenome_lineage


def test_domain_only_partner_taxon_uses_retained_classifier_lineage():
    row = {
        'partner_taxon': 'Bacteria',
        'partner_taxon_rank': 'domain',
        'taxonomy': 'd__Bacteria;p__Bacillota;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia;s__Roseburia intestinalis',
        'classification_assignment_status': 'CONFIDENT',
    }
    assert _pangenome_lineage(row)['g'] == 'Roseburia'
    assert _partner_rank_taxon(row) == ('s', 'species', 'Roseburia intestinalis')


def test_informative_partner_taxonomy_remains_authoritative():
    row = {
        'partner_taxon': 'g__Blautia',
        'partner_taxon_rank': 'genus',
        'taxonomy': 'd__Bacteria;p__Bacillota;g__Roseburia;s__Roseburia intestinalis',
        'classification_assignment_status': 'CONFIDENT',
    }
    assert _partner_rank_taxon(row) == ('g', 'genus', 'Blautia')


def test_domain_only_partner_does_not_fall_back_on_unresolved_classifier():
    row = {
        'partner_taxon': 'Bacteria',
        'partner_taxon_rank': 'domain',
        'taxonomy': 'd__Bacteria;p__Bacillota;g__Roseburia',
        'classification_assignment_status': 'UNRESOLVED',
    }
    assert _partner_rank_taxon(row) == ('d', 'domain', 'Bacteria')
