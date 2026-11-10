"""Shared file I/O utilities used across the BranchManager codebase.

Consolidates copy-pasted helpers that were duplicated across 6+ files:
_open_text, _write_tsv, and _sha256.
"""

import csv
import gzip
import hashlib
from pathlib import Path


def open_text(path, mode="rt", newline="", encoding="utf-8"):
    """Open a file path, transparently handling .gz compressed files.

    Replaces the many local ``_open_text`` copies scattered across the
    codebase.  The signature matches the most common variant
    (newline-aware text mode).
    """
    if str(path).lower().endswith('.gz'):
        return gzip.open(path, mode, newline=newline, encoding=encoding)
    return open(path, mode, newline=newline, encoding=encoding)


def write_tsv(path: Path, fieldnames, rows) -> None:
    """Write *rows* as a tab-separated file with the given *fieldnames*.

    Creates parent directories if needed and uses Unix line endings.
    Replaces the many local ``_write_tsv`` copies.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter='\t',
            extrasaction='ignore', lineterminator='\n',
        )
        writer.writeheader()
        writer.writerows(rows)


def sha256sum(path) -> str:
    """Return the hex SHA-256 digest of a file.

    Replaces the many local ``_sha256`` copies.
    """
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def sha256sum_gzip(path) -> str:
    """Return the hex SHA-256 digest of a file, transparently handling .gz.

    Used where the input may be a compressed trace file (e.g. AB1).
    """
    digest = hashlib.sha256()
    opener = gzip.open if str(path).lower().endswith('.gz') else open
    with opener(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()
