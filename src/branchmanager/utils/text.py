"""Shared text normalisation utilities used across the BranchManager codebase.

Consolidates copy-pasted helpers that were duplicated across 4+ files:
_norm_id, _normalise, _normalise_column, _normalise_for_matching.
"""

import re
from typing import Optional


def norm_id(x: str) -> str:
    """Light normalisation for FASTA/sequence IDs.

    Take the first whitespace token, prefer the last pipe field, strip
    coordinate suffixes, and collapse non-alphanumeric characters to
    underscores.  This is the canonical version shared by
    ``pipeline/classify.py`` and ``pipeline/background_check.py``.
    """
    if not x:
        return ""
    x = str(x).split()[0]
    if "|" in x:
        x = x.split("|")[-1]
    if "#" in x:
        x = x.split("#")[0]
    if ":" in x:
        x = x.split(":")[0]
    x = x.strip("()[]{}").strip()
    x = re.sub(r"[^0-9A-Za-z]+", "_", x).strip("_")
    return x


def canon_id(x: Optional[str]) -> Optional[str]:
    """Aggressive canonicalisation for FASTA/sequence IDs.

    Remove node fragments and coordinates, prefer the last pipe field,
    and collapse to underscores.  Companion to :func:`norm_id`.
    """
    if x is None:
        return x
    y = str(x).split()[0]
    if "#" in y:
        y = y.split("#")[0]
    y = re.sub(r":\d+-\d+\(.*\)$", "", y)
    y = re.sub(r":\d+-\d+$", "", y)
    if "|" in y:
        y = y.split("|")[-1]
    y = re.sub(r"[^0-9A-Za-z]+", "_", y).strip("_")
    return y


def normalise_column(name: object) -> str:
    """Normalise a TSV/CSV column name for consistent lookup.

    Identical to the copies in ``partner_metadata.py``,
    ``pipeline/diversity_metadata.py`` (simpler variant), and used as
    the basis for ``_normalise_column`` in ``pipeline/paper_trail.py``.
    """
    return str(name or '').strip().lower().replace('-', '_')


def normalise_for_matching(name: object) -> str:
    """Normalise a name for fuzzy matching (no separators).

    Identical to the copies in ``partner_metadata.py`` and
    ``pipeline/diversity_metadata.py``.
    """
    return normalise_column(name).replace('_', '').replace(' ', '')


def normalise_header(value: object) -> str:
    """Normalise a header value with BOM stripping and space replacement.

    Combines the identical ``_normalise`` functions in
    ``marker_provenance.py`` and ``project_state.py``, and the
    ``_normalise_header`` in ``personnel.py`` and ``pipeline/paper_trail.py``.
    """
    return str(value or '').strip().lstrip('\ufeff').lower().replace('-', '_').replace(' ', '_')
