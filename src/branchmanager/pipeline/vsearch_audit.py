"""Auditable execution wrapper for BranchManager VSEARCH searches."""
from __future__ import annotations

import csv
import hashlib
import shlex
import subprocess
import time
from pathlib import Path

from branchmanager.utils.subprocess import run_cmd

_REFERENCE_CACHE: dict[tuple[str, int, int], tuple[int, str]] = {}
_VERSION_CACHE: str | None = None


def _fasta_count(path: str | Path) -> int:
    p = Path(path)
    if not p.exists():
        return 0
    try:
        with p.open(encoding='utf-8', errors='replace') as handle:
            return sum(1 for line in handle if line.startswith('>'))
    except OSError:
        return 0


def _reference_facts(path: str | Path) -> tuple[int, str]:
    p = Path(path)
    try:
        stat = p.stat()
        key = (str(p.resolve()), stat.st_size, stat.st_mtime_ns)
    except OSError:
        return 0, 'NA'
    if key in _REFERENCE_CACHE:
        return _REFERENCE_CACHE[key]
    count = _fasta_count(p)
    digest = hashlib.sha256()
    try:
        with p.open('rb') as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b''):
                digest.update(block)
        checksum = digest.hexdigest()
    except OSError:
        checksum = 'NA'
    _REFERENCE_CACHE[key] = (count, checksum)
    return count, checksum


def _vsearch_version() -> str:
    global _VERSION_CACHE
    if _VERSION_CACHE is not None:
        return _VERSION_CACHE
    try:
        result = subprocess.run(
            ['vsearch', '--version'], capture_output=True, text=True, check=False,
        )
        text = (result.stdout or result.stderr or '').strip().splitlines()
        _VERSION_CACHE = text[0] if text else f'exit_{result.returncode}'
    except Exception as exc:  # version metadata must not mask the actual search
        _VERSION_CACHE = f'unavailable: {exc}'
    return _VERSION_CACHE


def run_vsearch_audited(
    command: str,
    *,
    audit_path: str | Path,
    stage: str,
    purpose: str,
    query_fasta: str | Path,
    reference_fasta: str | Path,
    runner=None,
):
    """Run VSEARCH and append reproducibility/performance metadata after success."""
    tokens = shlex.split(command)
    options = {}
    for flag in ('--id', '--query_cov', '--maxaccepts', '--maxrejects', '--maxhits', '--threads'):
        try:
            options[flag.lstrip('-')] = tokens[tokens.index(flag) + 1]
        except (ValueError, IndexError):
            options[flag.lstrip('-')] = 'default'
    try:
        userout = Path(tokens[tokens.index('--userout') + 1])
    except (ValueError, IndexError):
        userout = None

    ref_count, ref_sha256 = _reference_facts(reference_fasta)
    started = time.monotonic()
    result = (runner or run_cmd)(command)
    elapsed = time.monotonic() - started
    raw_hits = 0
    if userout and userout.exists():
        with userout.open(encoding='utf-8', errors='replace') as handle:
            raw_hits = sum(1 for line in handle if line.strip())

    fields = [
        'stage', 'purpose', 'query_fasta', 'reference_fasta', 'query_sequences',
        'reference_sequences', 'reference_sha256', 'identity_threshold',
        'query_coverage_threshold', 'maxaccepts', 'maxrejects', 'maxhits',
        'threads', 'runtime_seconds', 'raw_hit_count', 'vsearch_version',
        'userout',
    ]
    row = {
        'stage': stage, 'purpose': purpose, 'query_fasta': str(query_fasta),
        'reference_fasta': str(reference_fasta), 'query_sequences': _fasta_count(query_fasta),
        'reference_sequences': ref_count, 'reference_sha256': ref_sha256,
        'identity_threshold': options['id'], 'query_coverage_threshold': options['query_cov'],
        'maxaccepts': options['maxaccepts'], 'maxrejects': options['maxrejects'],
        'maxhits': options['maxhits'], 'threads': options['threads'],
        'runtime_seconds': f'{elapsed:.3f}', 'raw_hit_count': raw_hits,
        'vsearch_version': _vsearch_version(), 'userout': str(userout or 'NA'),
    }
    destination = Path(audit_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    write_header = not destination.exists() or destination.stat().st_size == 0
    with destination.open('a', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        if write_header:
            writer.writeheader()
        writer.writerow(row)
    return result
