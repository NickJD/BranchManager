import unittest

from branchmanager.taxonomy import (
    get_domain_or_kingdom,
    parse_reference_header_taxonomy,
    parse_taxon_string,
)


class SilvaReferenceHeaderTaxonomyTests(unittest.TestCase):
    def test_silva_bacterial_header_is_converted_to_ranked_lineage(self):
        header = (
            'HG374866.9.1343 Bacteria;Pseudomonadati;Pseudomonadota;'
            'Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;'
            'Pseudomonas;Rhizoctonia solani'
        )
        ref_id, taxonomy = parse_reference_header_taxonomy(header)

        self.assertEqual(ref_id, 'HG374866.9.1343')
        parsed = parse_taxon_string(taxonomy)
        self.assertEqual(parsed['d'], 'Bacteria')
        self.assertEqual(parsed['k'], 'Pseudomonadati')
        self.assertEqual(parsed['p'], 'Pseudomonadota')
        self.assertEqual(parsed['c'], 'Gammaproteobacteria')
        self.assertEqual(parsed['o'], 'Pseudomonadales')
        self.assertEqual(parsed['f'], 'Pseudomonadaceae')
        self.assertEqual(parsed['g'], 'Pseudomonas')
        self.assertEqual(parsed['s'], 'Rhizoctonia solani')
        self.assertEqual(get_domain_or_kingdom(taxonomy), 'Bacteria')

    def test_silva_eukaryotic_header_is_recognised(self):
        header = (
            'CU874468.1.1374 Eukaryota;Amorphea;Obazoa;Opisthokonta;'
            'Nucletmycea;Fungi;Dikarya;Ascomycota;Pezizomycotina;'
            'Sordariomycetes;Sordariales;Lasiosphaeriaceae;Podospora;'
            'Podospora anserina'
        )
        ref_id, taxonomy = parse_reference_header_taxonomy(header)
        self.assertEqual(ref_id, 'CU874468.1.1374')
        self.assertTrue(taxonomy.startswith('d__Eukaryota;k__Amorphea;'))

    def test_gtdb_ranked_header_remains_unchanged(self):
        header = 'RS_GCF_000001405.40 d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria'
        ref_id, taxonomy = parse_reference_header_taxonomy(header)
        self.assertEqual(ref_id, 'RS_GCF_000001405.40')
        self.assertEqual(
            taxonomy,
            'd__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria',
        )

    def test_unrecognised_description_is_not_mistaken_for_taxonomy(self):
        ref_id, taxonomy = parse_reference_header_taxonomy('ABC123 hypothetical protein')
        self.assertEqual(ref_id, 'ABC123')
        self.assertIsNone(taxonomy)


if __name__ == '__main__':
    unittest.main()
