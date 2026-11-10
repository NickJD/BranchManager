"""Safe, validated and atomic handling of compressed nucleotide reference FASTA files."""
from __future__ import annotations

import gzip
import hashlib
import json
import logging
import os
import tempfile
import time
from pathlib import Path

logger = logging.getLogger(__name__)

# IUPAC nucleotide alphabet, including U and alignment gap characters.
_VALID = frozenset(b"ACGTURYSWKMBDHVN-.")
_CHUNK = 8 * 1024 * 1024


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(_CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_fasta(path: str | Path) -> tuple[int, int]:
    """Validate FASTA structure and IUPAC nucleotide symbols.

    Returns (record_count, total_sequence_length). Reports the first invalid
    byte with its record identifier and 1-based sequence position.
    """
    path = Path(path)
    records = 0
    total = 0
    current_id = None
    position = 0
    with path.open("rb") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(b">"):
                if len(line) == 1:
                    raise ValueError(f"{path}: empty FASTA header at line {line_number}")
                records += 1
                current_id = line[1:].split(None, 1)[0].decode("utf-8", "replace")
                position = 0
                continue
            if current_id is None:
                raise ValueError(f"{path}: sequence data before first FASTA header at line {line_number}")
            for byte in line.upper():
                position += 1
                if byte not in _VALID:
                    raise ValueError(
                        f"{path}: invalid nucleotide byte 0x{byte:02x} in record "
                        f"{current_id!r}, sequence position {position} (line {line_number})"
                    )
            total += len(line)
    if records == 0 or total == 0:
        raise ValueError(f"{path}: FASTA contains no records or no sequence data")
    return records, total


def ensure_uncompressed_reference(
    ref_fasta: str | Path,
    outdir: str | Path,
    out_name: str = "ref_uncompressed.fasta",
) -> str:
    """Return a validated uncompressed reference, using an atomic cache.

    A JSON manifest binds the cache to the exact compressed source (SHA256,
    size and mtime). Cached FASTA content is independently hashed and validated.
    New output is written to a unique temporary file in the destination
    directory and atomically installed only after validation.
    """
    source = Path(ref_fasta)
    if source.suffix.lower() != ".gz":
        validate_fasta(source)
        return str(source)

    if not source.is_file():
        raise FileNotFoundError(f"Reference FASTA does not exist: {source}")

    destination_dir = Path(outdir)
    destination_dir.mkdir(parents=True, exist_ok=True)
    destination = destination_dir / out_name
    manifest = destination.with_suffix(destination.suffix + ".manifest.json")

    source_stat = source.stat()
    source_sha = _sha256(source)
    source_stat_after_hash = source.stat()
    if (source_stat.st_size, source_stat.st_mtime_ns) != (
        source_stat_after_hash.st_size, source_stat_after_hash.st_mtime_ns
    ):
        raise RuntimeError(f"Compressed reference changed while being checksummed: {source}")

    expected = {
        "source": str(source.resolve()),
        "source_size": source_stat.st_size,
        "source_mtime_ns": source_stat.st_mtime_ns,
        "source_sha256": source_sha,
    }

    if destination.is_file() and manifest.is_file():
        try:
            metadata = json.loads(manifest.read_text(encoding="utf-8"))
            if all(metadata.get(k) == v for k, v in expected.items()):
                records, bases = validate_fasta(destination)
                actual_sha = _sha256(destination)
                if metadata.get("output_sha256") == actual_sha:
                    logger.info(
                        "[REFERENCE] Reusing validated cache %s (%d records, %d bases; SHA256 %s)",
                        destination, records, bases, actual_sha,
                    )
                    return str(destination)
                logger.warning("[REFERENCE] Cached reference checksum mismatch; rebuilding %s", destination)
            else:
                logger.info("[REFERENCE] Cache provenance mismatch; rebuilding %s", destination)
        except Exception as exc:
            logger.warning("[REFERENCE] Existing cache failed validation (%s); rebuilding %s", exc, destination)

    fd, temporary_name = tempfile.mkstemp(
        prefix=destination.name + ".", suffix=".tmp", dir=str(destination_dir)
    )
    os.close(fd)
    temporary = Path(temporary_name)
    started = time.time()
    try:
        with gzip.open(source, "rb") as src, temporary.open("wb") as dst:
            while True:
                chunk = src.read(1 << 20)
                if not chunk:
                    break
                dst.write(chunk)
            dst.flush()
            os.fsync(dst.fileno())

        records, bases = validate_fasta(temporary)
        output_sha = _sha256(temporary)

        # Ensure the source was not replaced/modified during decompression.
        final_stat = source.stat()
        if (final_stat.st_size, final_stat.st_mtime_ns) != (
            source_stat.st_size, source_stat.st_mtime_ns
        ) or _sha256(source) != source_sha:
            raise RuntimeError(f"Compressed reference changed during decompression: {source}")

        os.replace(temporary, destination)
        metadata = {
            **expected,
            "output_sha256": output_sha,
            "records": records,
            "bases": bases,
            "validated_at_unix": time.time(),
        }
        manifest_fd, manifest_tmp_name = tempfile.mkstemp(
            prefix=manifest.name + ".", suffix=".tmp", dir=str(destination_dir)
        )
        try:
            with os.fdopen(manifest_fd, "w", encoding="utf-8") as handle:
                json.dump(metadata, handle, indent=2, sort_keys=True)
                handle.write("\n")
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(manifest_tmp_name, manifest)
        finally:
            if os.path.exists(manifest_tmp_name):
                os.unlink(manifest_tmp_name)

        logger.info(
            "[REFERENCE] Installed validated reference %s (%d records, %d bases, %.1f MB, %.1fs; SHA256 %s)",
            destination, records, bases, destination.stat().st_size / (1024 * 1024),
            time.time() - started, output_sha,
        )
        return str(destination)
    except Exception as exc:
        temporary.unlink(missing_ok=True)
        raise RuntimeError(f"Reference decompression/validation failed for {source}: {exc}") from exc
