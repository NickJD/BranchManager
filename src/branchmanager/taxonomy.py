from __future__ import annotations

import re
from typing import Dict, Optional

RANK_ALIASES = {
    'domain': 'd',
    'superkingdom': 'd',
    'kingdom': 'k',
    'phylum': 'p',
    'class': 'c',
    'order': 'o',
    'family': 'f',
    'genus': 'g',
    'species': 's',
}

DOMAIN_QUERY_ALIASES = {
    'bacterial': 'bacteria',
    'bacterium': 'bacteria',
    'bacteria': 'bacteria',
    'archaeal': 'archaea',
    'archaea': 'archaea',
    'fungal': 'fungi',
    'fungus': 'fungi',
    'fungi': 'fungi',
    'eukaryotic': 'eukaryota',
    'eukarya': 'eukaryota',
    'eukaryota': 'eukaryota',
}


def _clean_taxon_value(value: str) -> str:
    if value is None:
        return ''
    value = str(value).strip().strip('"\'')
    value = re.sub(r'\s+', ' ', value)
    return value


def normalise_taxon_name(value: str) -> str:
    value = _clean_taxon_value(value)
    if not value:
        return ''
    value = value.lower()
    value = re.sub(r'[^a-z0-9]+', '', value)
    return value


def normalise_domain_query(value: str) -> str:
    norm = normalise_taxon_name(value)
    return DOMAIN_QUERY_ALIASES.get(norm, norm)


_DOMAIN_DISPLAY_NAMES = {
    'bacteria': 'Bacteria',
    'archaea': 'Archaea',
    'fungi': 'Fungi',
    'eukaryota': 'Eukaryota',
}


def domain_designation(value: Optional[str]) -> Optional[str]:
    """Return a canonical domain name when *value* is domain-level only.

    Partner catalogues sometimes use values such as ``bacterial`` or
    ``Bacteria`` in a free-text ``taxon`` column to indicate broad domain
    membership rather than a genus/species identification. Such values must
    not be compared as specific taxa.
    """
    if value is None:
        return None
    normalised = normalise_domain_query(str(value))
    return _DOMAIN_DISPLAY_NAMES.get(normalised)



_KNOWN_BACTERIAL_PHYLA = {
    'firmicutes', 'bacteroidetes', 'proteobacteria', 'actinobacteria',
    'spirochaetes', 'fusobacteria', 'cyanobacteria', 'verrucomicrobia',
    'chloroflexi', 'planctomycetes', 'deinococcusthermus', 'tenericutes',
    'acidobacteria', 'chlamydiae', 'thermotogae', 'synergistetes',
    'deferribacteres', 'fibrobacteres', 'elusimicrobia', 'armatimonadetes',
    'bacillota', 'bacteroidota', 'pseudomonadota', 'actinomycetota',
    'verrucomicrobiota', 'spirochaetota', 'fusobacteriota', 'cyanobacteriota',
    'chloroflexota', 'planctomycetota', 'thermotogota', 'acidobacteriota',
    'chlamydiota', 'synergistota', 'fibrobacterota', 'elusimicrobiota',
}


def infer_taxon_rank(value: str, rank_hint: str | None = None) -> str:
    """Infer a short taxonomic rank for a partner-supplied free-text taxon.

    Explicit rank hints take precedence. Name-based inference is deliberately
    conservative and intended for common catalogue values, not as a substitute
    for a fully ranked lineage.
    """
    aliases = {
        'domain': 'd', 'superkingdom': 'd', 'kingdom': 'k',
        'phylum': 'p', 'class': 'c', 'order': 'o', 'family': 'f',
        'genus': 'g', 'species': 's',
        'd': 'd', 'k': 'k', 'p': 'p', 'c': 'c', 'o': 'o',
        'f': 'f', 'g': 'g', 's': 's',
    }
    hinted = aliases.get(_clean_taxon_value(rank_hint or '').lower())
    if hinted:
        return hinted
    if domain_designation(value):
        return 'd'
    clean = _clean_taxon_value(value)
    norm = normalise_taxon_name(clean)
    if norm in _KNOWN_BACTERIAL_PHYLA:
        return 'p'
    low = clean.lower()
    if low.endswith('aceae'):
        return 'f'
    if low.endswith('ales'):
        return 'o'
    if low.endswith(('ota', 'phyta')):
        return 'p'
    if low.endswith(('ia', 'etes')):
        return 'c'
    if len(clean.split()) >= 2:
        return 's'
    return 'g'

def parse_taxon_string(taxon: str) -> Dict[str, str]:
    """Parse lineage strings like d__Bacteria; p__Firmicutes into rank->name.

    Supports both short prefixes (d__/p__) and long prefixes (domain__/phylum__).
    Returns short one-letter rank keys where possible and also preserves unknown
    prefixes verbatim.
    """
    result: Dict[str, str] = {}
    if taxon is None:
        return result
    for raw_part in str(taxon).split(';'):
        part = raw_part.strip()
        if not part or '__' not in part:
            continue
        raw_rank, raw_value = part.split('__', 1)
        rank = raw_rank.strip().lower()
        value = _clean_taxon_value(raw_value)
        if not value:
            continue
        rank = RANK_ALIASES.get(rank, rank[:1] if rank else rank)
        result[rank] = value
    return result


def get_rank_value(taxon: str, rank: str) -> Optional[str]:
    parsed = parse_taxon_string(taxon)
    key = RANK_ALIASES.get(rank.lower(), rank.lower()[:1])
    return parsed.get(key)


def get_domain_or_kingdom(taxon: str) -> Optional[str]:
    parsed = parse_taxon_string(taxon)
    return parsed.get('d') or parsed.get('k')


def taxonomy_matches_kingdom(taxon: Optional[str], kingdom: str) -> bool:
    """Return True if lineage belongs to the requested kingdom/domain.

    Matching is rank-aware and prefers the domain/kingdom rank rather than a
    loose substring search across the whole lineage.
    """
    if not kingdom:
        return True
    wanted = normalise_domain_query(kingdom)
    if not wanted:
        return True
    parsed = parse_taxon_string(taxon)
    domain = normalise_domain_query(parsed.get('d', ''))
    kingdom_value = normalise_domain_query(parsed.get('k', ''))
    candidates = [value for value in (domain, kingdom_value) if value]
    if candidates:
        return wanted in candidates
    # fallback for malformed lineage strings with no explicit rank prefixes
    raw = normalise_domain_query(taxon)
    if not raw:
        return False
    return raw == wanted



def compare_expected_taxon(expected_taxon: Optional[str], observed_taxonomy: Optional[str]) -> str:
    """Compare a partner-declared taxon with a sequence-derived lineage.

    Returns one of:
      MATCH                 same taxon or a more-specific compatible assignment
      LESS_SPECIFIC         observed lineage supports only a broader rank
      DISAGREEMENT          observed lineage contradicts the declared taxon
      UNCLASSIFIED          no usable sequence-derived taxonomy

    This is deliberately conservative and name-based. It does not claim to
    resolve nomenclatural synonyms; those require a versioned taxonomy mapping.
    """
    expected = _clean_taxon_value(expected_taxon)
    observed = _clean_taxon_value(observed_taxonomy)
    if not observed or observed.lower() in {'na', 'none', 'unresolved', 'unclassified'}:
        return 'UNCLASSIFIED'
    if not expected:
        return 'MATCH'

    expected_lineage = parse_taxon_string(expected)
    observed_lineage = parse_taxon_string(observed)

    # If the partner supplied a ranked lineage, compare every shared rank.
    if expected_lineage:
        shared = [rank for rank in expected_lineage if rank in observed_lineage]
        if any(
            normalise_taxon_name(expected_lineage[rank]) != normalise_taxon_name(observed_lineage[rank])
            for rank in shared
        ):
            return 'DISAGREEMENT'
        if shared:
            return 'MATCH' if len(observed_lineage) >= len(expected_lineage) else 'LESS_SPECIFIC'

    expected_norm = normalise_taxon_name(expected)
    if not expected_norm:
        return 'MATCH'
    if not observed_lineage:
        observed_norm = normalise_taxon_name(observed)
        if expected_norm == observed_norm or observed_norm.startswith(expected_norm):
            return 'MATCH'
        expected_tokens = _clean_taxon_value(expected).split()
        if len(expected_tokens) >= 2 and observed_norm == normalise_taxon_name(expected_tokens[0]):
            return 'LESS_SPECIFIC'
        return 'DISAGREEMENT'

    observed_names = [normalise_taxon_name(value) for value in observed_lineage.values()]
    observed_names = [value for value in observed_names if value]
    # A free-text partner taxon can be a genus, species binomial, or higher taxon.
    # Check rank values and adjacent genus/species components without requiring
    # the partner and reference database to use identical lineage formatting.
    genus = observed_lineage.get('g', '')
    species = observed_lineage.get('s', '')
    genus_species = normalise_taxon_name(f'{genus} {species}')
    if expected_norm in observed_names or (genus_species and expected_norm == genus_species):
        return 'MATCH'
    if genus and expected_norm == normalise_taxon_name(genus):
        return 'MATCH'
    if genus_species and genus_species.startswith(expected_norm):
        return 'MATCH'

    # A species-level expectation can be compatible with a genus-only result.
    expected_tokens = _clean_taxon_value(expected).split()
    if len(expected_tokens) >= 2 and genus and normalise_taxon_name(expected_tokens[0]) == normalise_taxon_name(genus):
        if not species:
            return 'LESS_SPECIFIC'
        return 'DISAGREEMENT'

    # A broader expected name may be present as a lineage component even if
    # the reference has a more specific classification.
    if expected_norm in observed_names:
        return 'MATCH'
    return 'DISAGREEMENT'

def canonicalise_sequence_id(header: str) -> Optional[str]:
    """Conservative sequence-id canonicalisation.

    Keeps biologically meaningful suffixes such as `_1` / `_2` (ASV/OTU labels),
    while stripping transport/coordinate artefacts from FASTA headers.
    """
    if header is None:
        return None
    token = str(header).strip().split()[0]
    if not token:
        return None
    if '|' in token:
        token = token.split('|')[-1]
    if '#' in token:
        token = token.split('#', 1)[0]
    token = re.sub(r':\d+-\d+(?:\([^)]*\))?$', '', token)
    token = re.sub(r':\d+-\d+$', '', token)
    token = token.strip('()[]{}')
    return token or None


def reference_lookup_keys(header: Optional[str]) -> list[str]:
    """Return resilient lookup keys for reference/taxonomy ID matching."""
    if header is None:
        return []
    raw = str(header).strip()
    if not raw:
        return []
    keys = []
    first = raw.split()[0]
    for candidate in (raw, first, canonicalise_sequence_id(raw), canonicalise_sequence_id(first)):
        if candidate and candidate not in keys:
            keys.append(candidate)
    if '|' in first:
        last = first.split('|')[-1]
        if last and last not in keys:
            keys.append(last)
        canon_last = canonicalise_sequence_id(last)
        if canon_last and canon_last not in keys:
            keys.append(canon_last)
    return keys


def parse_reference_header_taxonomy(header: Optional[str]) -> tuple[Optional[str], Optional[str]]:
    """Extract ``(reference_id, taxonomy)`` from GTDB or SILVA FASTA headers.

    GTDB headers carry rank-prefixed lineages, for example::

        RS_GCF_... d__Bacteria;p__Pseudomonadota;...

    SILVA SSURef headers carry an accession followed by an unprefixed,
    semicolon-delimited lineage, for example::

        HG374866.9.1343 Bacteria;Pseudomonadati;Pseudomonadota;...

    SILVA's positional lineage is converted to the rank-prefixed representation
    used internally by BranchManager. The second SILVA field is retained as
    kingdom (``k__``), followed by phylum, class, order, family, genus and
    species. This preserves the SILVA lineage while allowing existing
    rank-aware consumers to read domain, phylum, genus and species consistently.
    """
    if header is None:
        return None, None

    parts = str(header).strip().split(None, 1)
    ref_id = parts[0] if parts else None
    remainder = parts[1] if len(parts) > 1 else ''
    if not remainder:
        return ref_id, None

    # GTDB and other explicitly rank-prefixed lineages remain unchanged.
    lineage = remainder.split('[', 1)[0].strip()
    lineage = re.sub(r'\s*;\s*', ';', lineage)
    if re.search(r'(?:^|;)\s*(?:d|domain)__', lineage, flags=re.IGNORECASE):
        return ref_id, lineage or None

    # SILVA SSURef: accession followed by a plain semicolon-delimited lineage.
    # Do not mistake arbitrary descriptions for taxonomy: the first component
    # must be a recognised domain-level name.
    components = [_clean_taxon_value(value) for value in lineage.split(';')]
    components = [value for value in components if value]
    if not components:
        return ref_id, None

    domain = normalise_taxon_name(components[0])
    if domain not in {'bacteria', 'archaea', 'eukaryota'}:
        return ref_id, None

    # SILVA's lineage has an intermediate kingdom/supergroup component between
    # domain and phylum (e.g. Bacteria;Pseudomonadati;Pseudomonadota;...).
    # Positional mapping is intentional: incomplete lineages simply omit later
    # ranks, rather than shifting a genus into the species rank.
    ranks = ('d', 'k', 'p', 'c', 'o', 'f', 'g', 's')
    ranked = [f'{rank}__{value}' for rank, value in zip(ranks, components)]
    return ref_id, ';'.join(ranked) if ranked else None