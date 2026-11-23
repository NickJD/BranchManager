import csv
from pathlib import Path

from branchmanager.sanger_studio import (
    _discover_cohort_manifest_for_review,
    _read_cohort_manifest_for_review,
    _read_map,
)


def test_review_tsv_locates_cohort_manifest_and_matching_partner_dataset(tmp_path):
    cohort = tmp_path / 'cohort'
    prep = cohort / '00_cohort_preparation'
    prep.mkdir(parents=True)
    review = prep / 'cohort_marker_review_template.tsv'
    review.write_text('sequence_id\tpartner_id\tdataset\nAU_0001\tAU\tAU_01\n', encoding='utf-8')
    reads = tmp_path / 'partner_reads'
    reads.mkdir()
    (reads / 'ss172_907R.ab1').write_bytes(b'fake ab1')
    sample_map = cohort / 'AU' / 'manifest.csv'
    sample_map.parent.mkdir()
    sample_map.write_text('sequence_id,read_file,primer,direction\nAU_0001,ss172_907R.ab1,907R,reverse\n', encoding='utf-8')
    manifest = cohort / 'partner_marker_cohort_manifest.csv'
    manifest.write_text(
        'partner_id,dataset,input_type,sample_map,read_dir\n'
        f'AU,AU_01,sanger,{sample_map.relative_to(cohort)},{reads}\n'
        'QUB,QUB_01,sanger,missing.csv,missing_reads\n', encoding='utf-8'
    )

    found, root = _discover_cohort_manifest_for_review(review)
    assert found == manifest
    assert root == cohort
    entries = _read_cohort_manifest_for_review(found, [{'partner_id': 'AU', 'dataset': 'AU_01'}])
    assert len(entries) == 1
    assert entries[0]['sample_map'] == sample_map.resolve()
    assert entries[0]['read_dir'] == reads.resolve()
    rows = _read_map(entries[0]['sample_map'], entries[0]['read_dir'])
    assert rows[0]['sequence_id'] == 'AU_0001'
    assert rows[0]['file'] == str((reads / 'ss172_907R.ab1').resolve())
