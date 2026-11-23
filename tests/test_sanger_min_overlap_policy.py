from branchmanager.pipeline.paper_trail import DEFAULT_QC_POLICY
from branchmanager.cli.parser import build_parser


def test_default_minimum_sanger_overlap_is_40_bp_in_qc_policy():
    assert DEFAULT_QC_POLICY.min_overlap == 40


def test_paper_trail_cli_default_minimum_overlap_is_40_bp():
    args = build_parser().parse_args(['paper-trail', '--input', 'reads', '-o', 'out'])
    assert args.min_overlap == 40


def test_assistant_cli_default_minimum_overlap_is_40_bp():
    args = build_parser().parse_args([
        'assistant', '--fasta', 'markers.fasta', '--partner-metadata', 'metadata.tsv',
        '--db', 'project.sqlite', '--dataset', 'Partner', '--ref', 'gtdb.fa', '-o', 'out',
    ])
    assert args.min_overlap == 40


def test_default_overlap_assembler_accepts_exactly_40_matching_bases_but_not_39():
    from branchmanager.pipeline.paper_trail import find_best_overlap

    assert find_best_overlap('A' * 100, 'A' * 40) is not None
    assert find_best_overlap('A' * 100, 'A' * 39) is None
