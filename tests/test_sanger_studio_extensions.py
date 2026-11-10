import random

from branchmanager.pipeline.paper_trail import (
    SangerRead,
    assemble_read_group,
    diagnose_overlap,
    infer_read_info,
)


def _read(read_id, sequence, direction, primer):
    return SangerRead(
        read_id=read_id,
        sequence_id='isolate_1',
        primer=primer,
        direction=direction,
        source_file=f'{read_id}.ab1',
        raw_sequence=sequence,
        raw_qualities=[35] * len(sequence),
        status='kept',
        oriented_sequence=sequence,
        oriented_qualities=[35] * len(sequence),
    )


def test_its_and_non_suffix_primer_directions_are_inferred():
    assert infer_read_info('sample_ITS1.ab1')[1:] == ('ITS1', 'forward')
    assert infer_read_info('sample_ITS2.ab1')[1:] == ('ITS2', 'reverse')
    assert infer_read_info('sample_NS2.ab1')[1:] == ('NS2', 'reverse')
    assert infer_read_info('sample_LR0R.ab1')[1:] == ('LR0R', 'reverse')


def test_two_forward_two_reverse_reads_can_form_one_consensus():
    random.seed(20261106)
    target = ''.join(random.choice('ACGT') for _ in range(2000))
    fragments = [target[0:800], target[500:1300], target[1000:1800], target[1500:2000]]
    reads = [
        _read('read_f1', fragments[0], 'forward', 'ITS1'),
        _read('read_f2', fragments[1], 'forward', 'ITS3'),
        _read('read_r1', fragments[2], 'reverse', 'ITS2'),
        _read('read_r2', fragments[3], 'reverse', 'ITS4'),
    ]
    consensus, stats = assemble_read_group(reads, min_overlap=80, min_identity=0.98)
    assert stats['status'] == 'assembled'
    assert stats['used_reads'] == 4
    assert set(stats['used_read_ids']) == {r.read_id for r in reads}
    assert consensus == target


def test_ab1_directory_input_groups_multiple_its_reads(tmp_path):
    from branchmanager.sanger_studio import SangerStudio

    trace_dir = tmp_path / 'chromatograms'
    trace_dir.mkdir()
    for name in (
        'isolateA_ITS1.ab1', 'isolateA_ITS2.ab1',
        'isolateA_ITS3.ab1', 'isolateA_ITS4.ab1',
    ):
        (trace_dir / name).write_bytes(b'')

    studio = SangerStudio(root=tmp_path, out=tmp_path / 'studio_output')
    loaded = studio.load_input_directory(trace_dir)

    assert loaded['input_kind'] == 'raw_chromatograms'
    assert loaded['read_count'] == 4
    assert len(loaded['samples']) == 1
    sample = loaded['samples'][0]
    assert sample['sequence_id'] == 'isolateA'
    assert sample['read_count'] == 4
    assert {(row['primer'], row['direction']) for row in sample['reads']} == {
        ('ITS1', 'forward'), ('ITS2', 'reverse'),
        ('ITS3', 'forward'), ('ITS4', 'reverse'),
    }


def test_custom_primer_sequence_parser_accepts_iupac_and_rejects_invalid_bases():
    import pytest
    from branchmanager.sanger_studio import _parse_custom_primer_sequences

    assert _parse_custom_primer_sequences('ITS1F2\tGAAGGTGAAGTCGTAACAAGG\nCUSTOM_R NNNRY') == {
        'ITS1F2': 'GAAGGTGAAGTCGTAACAAGG',
        'CUSTOM_R': 'NNNRY',
    }
    with pytest.raises(ValueError, match='non-IUPAC'):
        _parse_custom_primer_sequences('BAD_PRIMER ACCTZ')


def test_custom_primer_sequence_is_used_for_leading_primer_trimming():
    from branchmanager.pipeline.paper_trail import trim_leading_primer

    primer = 'ACGTTGCAACGTTGCAACGT'
    target = 'GATTACAGATTACAGATTACAGATTACA'
    trimmed, qualities, removed = trim_leading_primer(
        primer + target,
        [35] * (len(primer) + len(target)),
        'CUSTOM_F',
        primer_sequences={'CUSTOM_F': primer},
    )
    assert trimmed == target
    assert len(qualities) == len(target)
    assert removed == len(primer)


def test_overlap_diagnostics_report_rejected_candidate_reason():
    result = diagnose_overlap('A' * 120, 'A' * 120, min_overlap=150, min_identity=0.98)
    assert result['accepted'] is False
    assert result['compared_bases'] == 120
    assert result['identity'] == 1.0
    assert result['reason'] == 'OVERLAP_BELOW_MINIMUM'


def test_overlap_diagnostics_report_identity_failure():
    left = 'A' * 100 + 'C' * 100
    right = 'A' * 100 + 'G' * 90 + 'C' * 10
    result = diagnose_overlap(left, right, min_overlap=80, min_identity=0.98)
    assert result['accepted'] is False
    assert 'IDENTITY_BELOW_MINIMUM' in result['reason']


def test_consensus_evidence_does_not_count_ambiguous_n_as_agreement():
    from branchmanager.pipeline.paper_trail import _classify_consensus_observations

    disagreement, evidence = _classify_consensus_observations([('read1', 'G'), ('read2', 'N')])
    assert disagreement is False
    assert evidence == 'PARTIAL_AMBIGUOUS_SUPPORT'

    disagreement, evidence = _classify_consensus_observations([('read1', 'G'), ('read2', 'G')])
    assert disagreement is False
    assert evidence == 'MULTI_READ_AGREEMENT'

    disagreement, evidence = _classify_consensus_observations([('read1', 'G'), ('read2', 'A')])
    assert disagreement is True
    assert evidence == 'READ_DISAGREEMENT'


def test_trimmed_position_maps_back_to_raw_ab1_coordinates_for_both_directions():
    from branchmanager.sanger_studio import _raw_ab1_index_for_trimmed_position

    # Forward: trimmed base 1 follows the 12-base left trim.
    assert _raw_ab1_index_for_trimmed_position(1, 'forward', 12, 92, 100) == 12
    assert _raw_ab1_index_for_trimmed_position(80, 'forward', 12, 92, 100) == 91

    # Reverse: oriented base 1 is the last retained raw base.
    assert _raw_ab1_index_for_trimmed_position(1, 'reverse', 8, 88, 100) == 87
    assert _raw_ab1_index_for_trimmed_position(80, 'reverse', 8, 88, 100) == 8

    assert _raw_ab1_index_for_trimmed_position(81, 'reverse', 8, 88, 100) is None



def test_studio_and_assembler_share_global_read_to_consensus_alignment():
    from branchmanager.pipeline.paper_trail import (
        _align_read_to_consensus,
        _pairwise_aligner,
    )

    aligner = _pairwise_aligner()
    assert aligner.mode == 'global'
    assert aligner.match_score == 2.0
    assert aligner.mismatch_score == -3.0
    assert aligner.open_internal_insertion_score == -6.0
    assert aligner.open_internal_deletion_score == -6.0
    assert aligner.extend_internal_insertion_score == -1.0
    assert aligner.extend_internal_deletion_score == -1.0
    assert aligner.open_end_insertion_score == 0.0
    assert aligner.open_end_deletion_score == 0.0
    assert aligner.wildcard == 'N'

    # A read covering an internal section of a longer candidate must retain
    # its true candidate coordinates; terminal overhangs are free, not clipped
    # by a separate local-alignment policy.
    candidate = 'AAAACCCCGGGGTTTT'
    read = 'CCCCGGGG'
    assert _align_read_to_consensus(candidate, read) == [(4, 12, 0, 8)]


def test_sanger_studio_does_not_define_a_second_local_aligner():
    from pathlib import Path
    import branchmanager.sanger_studio as studio

    source = Path(studio.__file__).read_text()
    assert "aligner.mode = 'local'" not in source
    assert 'paper_trail._align_read_to_consensus(candidate_seq, read_seq)' in source
