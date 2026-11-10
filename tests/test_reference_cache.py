import gzip
import hashlib
import json
from pathlib import Path

import pytest

from branchmanager.utils.reference import ensure_uncompressed_reference, validate_fasta


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def test_atomic_decompression_is_byte_identical_and_cached(tmp_path):
    source = tmp_path / "reference.fna.gz"
    payload = b">seq1 description\nACGTNRY\n>seq2\nTTTTACGT\n"
    with gzip.open(source, "wb") as handle:
        handle.write(payload)

    result = Path(ensure_uncompressed_reference(source, tmp_path / "out"))
    assert result.read_bytes() == payload
    assert validate_fasta(result) == (2, 15)

    manifest = result.with_suffix(result.suffix + ".manifest.json")
    metadata = json.loads(manifest.read_text())
    assert metadata["output_sha256"] == _sha(result)

    # A second call reuses the validated cache.
    assert ensure_uncompressed_reference(source, tmp_path / "out") == str(result)


def test_corrupt_cache_is_rebuilt(tmp_path):
    source = tmp_path / "reference.fna.gz"
    payload = b">seq1\nACGTACGT\n"
    with gzip.open(source, "wb") as handle:
        handle.write(payload)

    result = Path(ensure_uncompressed_reference(source, tmp_path / "out"))
    result.write_bytes(b">seq1\nACGT\x01CGT\n")

    rebuilt = Path(ensure_uncompressed_reference(source, tmp_path / "out"))
    assert rebuilt.read_bytes() == payload
    assert validate_fasta(rebuilt) == (1, 8)


def test_invalid_iupac_character_reports_record_and_position(tmp_path):
    fasta = tmp_path / "bad.fasta"
    fasta.write_bytes(b">bad_record\nACGT\x15TT\n")
    with pytest.raises(ValueError, match="bad_record.*sequence position 5"):
        validate_fasta(fasta)


def test_uncompressed_reference_is_validated(tmp_path):
    fasta = tmp_path / "bad.fasta"
    fasta.write_bytes(b">bad\nACGTuP\n")
    with pytest.raises(ValueError, match="invalid nucleotide byte"):
        ensure_uncompressed_reference(fasta, tmp_path / "out")
