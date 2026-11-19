from types import SimpleNamespace

import pytest

from branchmanager.pipeline import paper_trail


class _FakeRecord:
    seq = 'ACGT'
    letter_annotations = {'phred_quality': [40, 40, 40, 40]}

    def __init__(self, raw):
        self.annotations = {'abif_raw': raw}


def _patch_abi_reader(monkeypatch, raw):
    from Bio import SeqIO

    monkeypatch.setattr(SeqIO, 'read', lambda handle, fmt: _FakeRecord(raw))


def test_ab1_channel_order_is_taken_from_each_records_fwo1(tmp_path, monkeypatch):
    # DATA9..DATA12 are G,A,T,C for this record, not A,C,G,T.
    raw = {
        'FWO_1': b'GATC',
        'PLOC2': [0, 1, 2, 3],
        'DATA9': [90, 0, 0, 0],
        'DATA10': [0, 80, 0, 0],
        'DATA11': [0, 0, 70, 0],
        'DATA12': [0, 0, 0, 60],
    }
    _patch_abi_reader(monkeypatch, raw)
    path = tmp_path / 'test.ab1'
    path.write_bytes(b'fake')

    details = paper_trail._read_ab1_details(path)

    assert details['trace_available'] is True
    assert details['trace_order'] == 'GATC'
    assert details['trace_order_valid'] is True
    assert details['trace_channels'] == {
        'G': [90, 0, 0, 0],
        'A': [0, 80, 0, 0],
        'T': [0, 0, 70, 0],
        'C': [0, 0, 0, 60],
    }


def test_invalid_or_missing_fwo_does_not_guess_channel_order(tmp_path, monkeypatch):
    raw = {
        'FWO_1': b'NNNN',
        'PLOC2': [1, 2, 3, 4],
        'DATA9': [90, 0, 0, 0],
        'DATA10': [0, 80, 0, 0],
        'DATA11': [0, 0, 70, 0],
        'DATA12': [0, 0, 0, 60],
    }
    _patch_abi_reader(monkeypatch, raw)
    path = tmp_path / 'invalid.ab1'
    path.write_bytes(b'fake')

    details = paper_trail._read_ab1_details(path)

    assert details['trace_available'] is False
    assert details['trace_order'] == ''
    assert details['trace_order_valid'] is False
    assert details['trace_channels'] == {}


def test_trace_peak_evidence_uses_the_called_base_channel_as_denominator():
    ratios, secondary = paper_trail._trace_peak_evidence(
        'T', [0], {'A': [54], 'C': [1], 'G': [0], 'T': [45]}
    )
    assert secondary == ['A']
    assert ratios[0] == pytest.approx(54 / 45)


def test_hover_readout_distinguishes_raw_ab1_call_from_oriented_alignment_call():
    from pathlib import Path

    html = (Path(__file__).parents[1] / 'src' / 'branchmanager' / 'data' / 'sanger_studio.html').read_text()
    assert 'raw AB1 call' in html
    assert 'aligned call' in html
    assert 'peak_channel_signals' in html
    assert 'secondary oriented channel' in html
    assert "const traceComplement={A:'T',C:'G',G:'C',T:'A'}" in html
    assert 'rawChannels[isReverse?(traceComplement[base]||base):base]' in html
    assert 'orientedTraceBase' in html


def test_trace_is_unavailable_when_peak_locations_are_shorter_than_calls(tmp_path, monkeypatch):
    raw = {
        'FWO_1': b'ACGT',
        'PLOC2': [1, 2],
        'DATA9': [10, 10, 10, 10],
        'DATA10': [10, 10, 10, 10],
        'DATA11': [10, 10, 10, 10],
        'DATA12': [10, 10, 10, 10],
    }
    _patch_abi_reader(monkeypatch, raw)
    path = tmp_path / 'short_ploc.ab1'
    path.write_bytes(b'fake')

    details = paper_trail._read_ab1_details(path)

    assert details['trace_available'] is False
    assert details['secondary_peak_ratios'] == [0.0] * len(details['sequence'])


def test_trace_is_unavailable_when_peak_location_is_outside_channel_arrays(tmp_path, monkeypatch):
    raw = {
        'FWO_1': b'ACGT',
        'PLOC2': [1, 99, 2, 3],
        'DATA9': [10, 10, 10, 10],
        'DATA10': [10, 10, 10, 10],
        'DATA11': [10, 10, 10, 10],
        'DATA12': [10, 10, 10, 10],
    }
    _patch_abi_reader(monkeypatch, raw)
    path = tmp_path / 'bad_ploc.ab1'
    path.write_bytes(b'fake')

    details = paper_trail._read_ab1_details(path)

    assert details['trace_available'] is False
    assert details['secondary_peak_bases'] == [''] * len(details['sequence'])
