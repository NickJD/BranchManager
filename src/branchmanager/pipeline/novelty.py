"""
novelty.py — Novelty detection for BranchManager.

Novelty is expressed as the identity percentage of the best hit against
the reference database. This gives a gradient from 0 (no hit at all, or
very distant) to 100 (identical to a known sequence).

The distinction from classify.py is important:
  - classify.py performs taxonomic assignment using a high-stringency adaptive
    identity/coverage ladder with bounded candidate collection
  - novelty.py uses a descending 90→50% identity ladder, preferring >=95%
    query coverage and retaining >=80% partial-coverage evidence separately;
    nearest-hit searches are bounded, while density searches remain exhaustive
    because they calculate exact qualifying-neighbour counts

The two vsearch runs serve different purposes and MUST NOT share matches.tsv.
classify.py  → outdir/matches.tsv         (adaptive high-stringency taxonomy)
novelty.py   → outdir/novelty_matches.tsv (nearest-neighbour evidence)

Output columns: ID\tNearestIdentity\tNearestHit\tNovel
  NearestIdentity : float 0–100, identity to best reference hit
                    0.0 means no hit found at any identity level
  NearestHit      : accession of best reference sequence
  Novel           : True if NearestIdentity < id_threshold * 100
"""

from __future__ import annotations

import logging
import os
import shlex
import time
from pathlib import Path
from typing import Optional

from branchmanager.pipeline import tree as tree_pipeline
from branchmanager.pipeline.selection_sets import (
    DEFAULT_BASELINE_EXTENSION_MIN_IDENTITY,
    DEFAULT_BASELINE_EXTENSION_MIN_QUERY_COVERAGE,
    DEFAULT_BASELINE_REDUNDANCY_IDENTITY,
    DEFAULT_BASELINE_REDUNDANCY_QUERY_COVERAGE,
    DEFAULT_PANGENOME_TARGET,
)
from branchmanager.taxonomy import normalise_taxon_name, parse_taxon_string
from branchmanager.utils.fasta import read_fasta, write_fasta
from branchmanager.utils.subprocess import run_cmd
from branchmanager.pipeline.vsearch_audit import run_vsearch_audited
from branchmanager.utils.reference import ensure_uncompressed_reference

logger = logging.getLogger(__name__)
VSEARCH_USERFIELDS = 'query+target+id+alnlen+mism+gaps+qlo+qhi+tlo+thi+ql+tl+qcov+tcov'
MIN_NEAREST_QUERY_COVERAGE = DEFAULT_BASELINE_EXTENSION_MIN_QUERY_COVERAGE
PARTIAL_NEAREST_QUERY_COVERAGE = 80.0
# Density is a prioritisation signal, not a census. Bound accepted hits to
# reduce repeated all-vs-all work as partner/reference collections grow.
DENSITY_HIT_CAP = 100
NEAREST_IDENTITY_LADDER = (90, 80, 70, 60, 50)


def _alignment_evidence(parts: list[str]) -> tuple[Optional[float], Optional[int]]:
    """Return VSEARCH query coverage and alignment length from a userout row."""
    if len(parts) < 12:
        return None, None
    try:
        alignment_length = int(float(parts[3]))
        query_length = int(float(parts[10]))
    except (TypeError, ValueError):
        return None, None
    try:
        query_coverage = float(parts[12]) if len(parts) >= 14 else (
            100.0 * alignment_length / query_length if query_length > 0 else None
        )
    except (TypeError, ValueError, ZeroDivisionError):
        query_coverage = None
    if query_coverage is not None:
        query_coverage = max(0.0, min(100.0, query_coverage))
    return query_coverage, alignment_length


def _write_genome_derived_context(input_fasta: str, outdir: str, db, qids, threads=None) -> Optional[str]:
    """Compare queries with Barrnap-derived 16S evidence without changing primary scores."""
    if db is None:
        return None
    try:
        pool_path, _source, datasets, count = _build_db_pool_fasta(
            db=db, outdir=outdir, name='genome_derived_16s', evidence_pool='genome_derived',
        )
    except Exception as exc:
        logger.warning('[NOVELTY] Could not build the genome-derived 16S context pool: %s', exc)
        return None
    if not pool_path:
        return None
    matches_path = Path(outdir) / 'genome_derived_16s_matches.tsv'
    best = _run_nearest_search(input_fasta, pool_path, matches_path, threads=threads, exclude_self=True)
    summary_path = Path(outdir) / 'genome_derived_16s_context.tsv'
    with open(summary_path, 'w', newline='') as handle:
        writer = __import__('csv').DictWriter(
            handle,
            fieldnames=['sequence_id', 'genome_derived_nearest_hit', 'identity', 'query_coverage',
                        'alignment_length', 'context_datasets', 'context_sequence_count'],
            delimiter='\t',
        )
        writer.writeheader()
        for sequence_id in qids:
            hit = best.get(sequence_id)
            writer.writerow({
                'sequence_id': sequence_id,
                'genome_derived_nearest_hit': hit[1] if hit else 'NA',
                'identity': f'{hit[0]:.3f}' if hit else 'NA',
                'query_coverage': f'{hit[2]:.3f}' if hit and hit[2] is not None else 'NA',
                'alignment_length': hit[3] if hit and hit[3] is not None else 'NA',
                'context_datasets': ';'.join(datasets),
                'context_sequence_count': count,
            })
    logger.info('[NOVELTY] Wrote supplementary Barrnap-derived 16S context for %d query sequence(s) to %s', len(best), summary_path)
    return str(summary_path)


def run_novelty(
    input_fasta: str,
    ref_fasta: str,
    outdir: str,
    db=None,
    run_dataset: Optional[str] = None,
    id_threshold: float = 0.97,
    threads: Optional[int] = None,
    target_fasta: Optional[str] = None,
) -> str:
    """
    Compute novelty scores for all sequences in input_fasta.

    Novelty is measured against previously submitted sequences (Filing Cabinet + other
    runs already in the DB).  This means a sequence is "novel" if it is more
    than ``1 - id_threshold`` different from everything the user has already
    submitted to BranchManager, not from a global reference.

    If ``target_fasta`` is given, novelty is measured against those sequences
    instead of (or as a supplement to) the DB-derived sequences.

    Parameters
    ----------
    input_fasta    : Query sequences (dereplicated, QC-passed)
    ref_fasta      : Reference database FASTA (accepted for older callers;
                     novelty is measured against submitted sequences)
    outdir         : Output directory
    db             : Database instance to pull other submitted sequences from
    run_dataset    : Current run dataset name (used to exclude self from comparison)
    id_threshold   : Sequences with identity below this are marked Novel=True
                     (default 0.97 = 97% identity, classic species boundary)
    threads        : VSEARCH thread count
    target_fasta   : Explicit FASTA of target sequences; overrides DB lookup

    Returns
    -------
    Path to novelty.tsv
    """
    output = f"{outdir}/novelty.tsv"
    matches_path = f"{outdir}/novelty_matches.tsv"
    query_ids = [
        header for header, _sequence in read_fasta(input_fasta)
        if not tree_pipeline.is_ref_anchor(header)
    ]
    if not query_ids:
        with open(output, 'w') as handle:
            handle.write('ID\tNearestIdentity\tNearestHit\tNovel\n')
        Path(matches_path).write_text('')
        return output

    # Preferred path: reuse the authoritative Stage B project nearest-hit
    # results produced by build_reference_novelty_metrics. This prevents the
    # legacy report from launching a second VSEARCH against a different pool.
    project_matches = Path(outdir) / 'novelty_project_nearest_matches.tsv'
    if project_matches.exists():
        import shutil
        shutil.copyfile(project_matches, matches_path)
        best_hits = _parse_best_hits(str(project_matches))
        novel_threshold_pct = id_threshold * 100.0
        with open(output, "w") as fh:
            fh.write("ID\tNearestIdentity\tNearestHit\tNovel\n")
            for header, _sequence in read_fasta(input_fasta):
                if tree_pipeline.is_ref_anchor(header):
                    continue
                hit = best_hits.get(header, (0.0, "None", None, None))
                fh.write(
                    f"{header}\t{hit[0]:.2f}\t{hit[1]}\t"
                    f"{hit[0] < novel_threshold_pct}\n"
                )
        logger.info(
            "[NOVELTY] Derived legacy novelty outputs from Stage B project matches; "
            "no additional VSEARCH search was run"
        )
        return output

    # Compatibility fallback for direct callers that have not run the
    # consolidated reference/project novelty workflow first.
    # Priority: explicit target_fasta > DB-derived submitted sequences > nothing
    db_fasta = Path(outdir) / "submitted_sequences.fasta"
    submission_db_exists = False

    if target_fasta and Path(target_fasta).exists():
        # Explicit target sequences override DB lookup
        db_fasta = Path(target_fasta)
        submission_db_exists = True
        logger.info("[NOVELTY] Using explicit --target FASTA for novelty comparison: %s", target_fasta)

    elif db is not None:
        try:
            # Use the v1 provenance-aware pool once the registry has been
            # migrated. Before migration, _query_db_pool intentionally retains
            # the legacy behaviour so existing installations can transition.
            records, pool_datasets = _query_db_pool(
                db, run_dataset=run_dataset, evidence_pool='primary',
            )
            if records:
                write_fasta(records, str(db_fasta))
                submission_db_exists = True
            logger.info(
                "[NOVELTY] Primary comparison pool contains %d eligible marker(s) "
                "from %d dataset(s): %s",
                len(records), len(pool_datasets), pool_datasets,
            )
        except Exception as e:
            raise RuntimeError(f'failed to build the submitted-sequence novelty pool: {e}') from e

    # Run vsearch against submission database if available
    if submission_db_exists and db_fasta.exists():
        thread_flag = f" --threads {int(threads)}" if threads and int(threads) > 0 else ""
        cmd = (
            f"vsearch --usearch_global {shlex.quote(str(input_fasta))} "
            f"--db {shlex.quote(str(db_fasta))} "
            f"--id 0.5 "                          # low floor — we want real distances
            f"--strand both "
            f"--userout {shlex.quote(str(matches_path))} "
            f"--userfields {VSEARCH_USERFIELDS} "
            f"--maxaccepts 32 --maxrejects 256 --maxhits 32 "
            f"--query_cov {MIN_NEAREST_QUERY_COVERAGE / 100.0:.2f}"
            f"{thread_flag}"
        )
        logger.info("[NOVELTY] Running vsearch against previously submitted sequences")
        try:
            run_vsearch_audited(
                cmd, audit_path=Path(outdir) / 'vsearch_search_audit.tsv',
                stage='B', purpose='legacy_project_novelty', query_fasta=input_fasta,
                reference_fasta=db_fasta, runner=run_cmd,
            )
            logger.info("[NOVELTY] vsearch done → %s", matches_path)
        except Exception as e:
            raise RuntimeError(f'vsearch against the submitted-sequence pool failed: {e}') from e
    else:
        logger.info(
            "[NOVELTY] No previously submitted sequences available and no --target FASTA provided; "
            "all input sequences will be treated as novel. Consider running `branchmanager filing-cabinet` first "
            "to build a baseline database."
        )

    # Parse best hits: query_id -> (pct_identity, hit_id)
    best_hits = _parse_best_hits(matches_path)

    novel_threshold_pct = id_threshold * 100.0

    with open(output, "w") as fh:
        fh.write("ID\tNearestIdentity\tNearestHit\tNovel\n")
        for header, _ in read_fasta(input_fasta):
            # Skip reference anchors — they are scaffolding, not results
            if tree_pipeline.is_ref_anchor(header):
                continue
            hit = best_hits.get(header, (0.0, "None", None, None))
            pct_id, hit_id = hit[0], hit[1]
            is_novel = pct_id < novel_threshold_pct
            fh.write(f"{header}\t{pct_id:.2f}\t{hit_id}\t{is_novel}\n")

    n_total = sum(1 for _ in read_fasta(input_fasta) if not tree_pipeline.is_ref_anchor(_[0]))
    n_novel = sum(1 for v in best_hits.values() if v[0] < novel_threshold_pct)
    n_unmatched = n_total - len(best_hits)

    logger.info(
        "[NOVELTY] %d sequences: %d novel (<%d%%), %d unmatched (treated as novel), %d known (by submission history)",
        n_total,
        n_novel + n_unmatched,
        int(novel_threshold_pct),
        n_unmatched,
        n_total - n_novel - n_unmatched,
    )
    return output


def _parse_best_hits(matches_path: str) -> dict[str, tuple[float, str, Optional[float], Optional[int]]]:
    """
    Parse VSEARCH blast6 output.
    Returns {query_id: (pct_identity_float, hit_id_str)}.
    Only the best (highest identity) hit per query is kept.
    """
    best = {}
    try:
        with open(matches_path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                qid = parts[0]
                hit_id = parts[1]
                try:
                    pct_id = float(parts[2])
                except ValueError:
                    continue
                query_coverage, alignment_length = _alignment_evidence(parts)
                if qid not in best or pct_id > best[qid][0]:
                    best[qid] = (pct_id, hit_id, query_coverage, alignment_length)
    except FileNotFoundError:
        logger.warning("[NOVELTY] matches file not found: %s", matches_path)
    return best


def _ensure_uncompressed(ref_fasta: str, outdir: str) -> str:
    """Return a validated, source-specific uncompressed reference cache."""
    source = Path(ref_fasta)
    safe_stem = ''.join(ch if ch.isalnum() else '_' for ch in source.name.removesuffix('.gz'))[:48]
    return ensure_uncompressed_reference(
        ref_fasta, outdir, out_name=f"ref_uncompressed_{safe_stem}.fasta"
    )


def build_novelty_itol(
    novelty_tsv: str,
    outdir: str,
) -> str:
    """
    Generate an iTOL DATASET_COLORSTRIP file showing novelty as a colour
    gradient from red (<90% identity) through orange/yellow to green (100%).

    This replaces the binary True/False novelty flag with a proper gradient
    that matches the Hungate-style figure.
    """
    output = os.path.join(outdir, "itol_novelty.itol")

    # Identity bins → colour (matches the legend in your existing figures)
    def _identity_colour(pct: float) -> str:
        if pct == 0.0:
            return "#8B0000"   # dark red — no hit
        elif pct < 90:
            return "#FF0000"   # red
        elif pct < 91:
            return "#FF2200"
        elif pct < 92:
            return "#FF5500"
        elif pct < 93:
            return "#FF7700"
        elif pct < 94:
            return "#FF9900"
        elif pct < 95:
            return "#FFBB00"
        elif pct < 96:
            return "#DDDD00"
        elif pct < 97:
            return "#AADD00"
        elif pct < 98:
            return "#77DD00"
        elif pct < 99:
            return "#44CC00"
        elif pct < 100:
            return "#22BB00"
        else:
            return "#00AA00"   # dark green — 100% match

    legend_bins = [
        ("<90%",  "#FF0000"),
        ("90%",   "#FF2200"),
        ("91%",   "#FF5500"),
        ("92%",   "#FF7700"),
        ("93%",   "#FF9900"),
        ("94%",   "#FFBB00"),
        ("95%",   "#DDDD00"),
        ("96%",   "#AADD00"),
        ("97%",   "#77DD00"),
        ("98%",   "#44CC00"),
        ("99%",   "#22BB00"),
        ("100%",  "#00AA00"),
    ]

    id_colour_pairs = []
    try:
        with open(novelty_tsv) as fh:
            next(fh)  # skip header
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) < 2:
                    continue
                qid = parts[0]
                try:
                    pct = float(parts[1])
                except ValueError:
                    pct = 0.0
                id_colour_pairs.append((qid, _identity_colour(pct)))
    except FileNotFoundError:
        logger.warning("[NOVELTY] novelty.tsv not found: %s", novelty_tsv)
        return output

    with open(output, "w") as fh:
        fh.write("DATASET_COLORSTRIP\n")
        fh.write("SEPARATOR COMMA\n")
        fh.write("DATASET_LABEL,Novelty (identity)\n")
        fh.write("COLOR,#FF0000\n")
        fh.write("MARGIN,5\n")
        fh.write("SHOW_INTERNAL,0\n")
        fh.write("LEGEND_TITLE,Novelty (nearest identity)\n")
        fh.write("LEGEND_SHAPES," + ",".join(["1"] * len(legend_bins)) + "\n")
        fh.write("LEGEND_COLORS," + ",".join(c for _, c in legend_bins) + "\n")
        fh.write("LEGEND_LABELS," + ",".join(label for label, _ in legend_bins) + "\n")
        fh.write("DATA\n")
        for qid, colour in id_colour_pairs:
            fh.write(f"{qid},{colour}\n")

    logger.info("[NOVELTY] iTOL novelty strip → %s", output)
    return output


def compute_db_nearest_identities(mapped_derep: str, outdir: str, db, run_dataset: str, threads: Optional[int] = None):
    """Compute nearest-neighbour identities of run sequences against all non-run DB sequences.

    The comparison pool includes ALL datasets stored in the DB except *run_dataset*,
    which means every registered baseline dataset (no matter how many) is automatically included.
    """
    db_preset_fasta = Path(outdir) / 'project_collection_reference.fasta'
    try:
        records, pool_ds = _query_db_pool(
            db, run_dataset=run_dataset, evidence_pool='primary',
        )
        if not records:
            return {}
        write_fasta(records, str(db_preset_fasta))
        logger.info(
            "[NOVELTY] Nearest-identity primary marker pool: %d sequences from datasets: %s",
            len(records), pool_ds,
        )
    except Exception:
        return {}

    novel_results = {}
    if db_preset_fasta.exists():
        try:
            from branchmanager.utils.subprocess import run_cmd
            novelty_matches = Path(outdir) / 'novelty_matches.tsv'
            thread_flag = f" --threads {int(threads)}" if threads and int(threads) > 0 else ""
            cmd = (
                f"vsearch --usearch_global {shlex.quote(str(mapped_derep))} "
                f"--db {shlex.quote(str(db_preset_fasta))} --id 0.5 --strand both "
                f"--userout {shlex.quote(str(novelty_matches))} "
                f"--userfields {VSEARCH_USERFIELDS} "
                f"--maxaccepts 32 --maxrejects 256 --maxhits 32 "
                f"--query_cov {MIN_NEAREST_QUERY_COVERAGE / 100.0:.2f}{thread_flag}"
            )
            logger.info("[NOVELTY] Running vsearch against registered project collection to compute novelty")
            run_vsearch_audited(
                cmd, audit_path=Path(outdir) / 'vsearch_search_audit.tsv',
                stage='B', purpose='legacy_project_novelty', query_fasta=mapped_derep,
                reference_fasta=db_preset_fasta, runner=run_cmd,
            )
            if novelty_matches.exists():
                with open(novelty_matches) as nm:
                    for line in nm:
                        parts = line.strip().split('\t')
                        if len(parts) < 3:
                            continue
                        try:
                            identity = float(parts[2])
                        except Exception:
                            continue
                        if parts[0] not in novel_results or identity > novel_results[parts[0]]:
                            novel_results[parts[0]] = identity
        except Exception as e:
            logger.warning("[NOVELTY] vsearch novelty run failed: %s", e)
    return novel_results


def build_run_novelty_itol(outdir: str, run_ids, mapped_derep: str, db, run_dataset: str, orig_to_short: dict, threads: Optional[int] = None):
    """Build the run novelty iTOL strip from nearest-neighbour identities."""
    if not run_ids:
        return None

    from branchmanager.pipeline import itol

    nov_lines = [
        'DATASET_COLORSTRIP',
        'SEPARATOR COMMA',
        'DATASET_LABEL,Novelty (nearest identity)',
        'COLOR,#AAAAAA',
        'MARGIN,5',
        'SHOW_INTERNAL,0'
    ]
    try:
        percents = list(range(90, 101))
        colours = [itol._novelty_colour_for_pct('<90')]
        colours += [itol._identity_to_colour(p / 100.0, vmin=0.90, vmax=1.0) for p in percents]
        labels = ['<90'] + [f"{p}%" for p in percents]
        nov_lines.append('LEGEND_TITLE,Novelty (identity)')
        nov_lines.append('LEGEND_SHAPES,' + ','.join(['1'] * len(labels)))
        nov_lines.append('LEGEND_COLORS,' + ','.join(colours))
        nov_lines.append('LEGEND_LABELS,' + ','.join(labels))
    except Exception:
        logger.debug('Could not generate iTOL novelty legend; skipping')

    novel_results = compute_db_nearest_identities(mapped_derep, outdir, db, run_dataset, threads=threads)

    def _candidate_forms_for_short(x):
        forms = []
        try:
            s = orig_to_short.get(x)
            if s:
                forms.append(s)
        except Exception:
            pass  # orig_to_short lookup failed; continuing with original ID
        forms.append(x)
        if '|' in x:
            forms.append(x.split('|')[-1])
        forms.append(x.split()[-1] if x.split() else x)
        try:
            forms.append(forms[0].upper())
            forms.append(forms[0].lower())
        except Exception:
            pass  # case conversion failed; continuing with original forms
        import re as _re
        forms.append(_re.sub(r'_n?\d+$', '', forms[0] if forms else x))
        seen = set()
        out = []
        for f in forms:
            if not f:
                continue
            if f not in seen:
                seen.add(f)
                out.append(f)
        return out

    nov_lines.append('DATA')
    for iid in run_ids:
        short_id = None
        if novel_results:
            for cand in _candidate_forms_for_short(iid):
                if cand in novel_results:
                    short_id = cand
                    break
        if short_id is None:
            short_id = iid
        colour = '#cccccc'
        if novel_results and short_id in novel_results:
            v = novel_results.get(short_id)
            if v is not None:
                try:
                    colour = itol._identity_to_colour(float(v), vmin=0.90, vmax=1.0)
                except Exception:
                    pass  # could not parse identity value; using default colour
                    colour = '#cccccc'
            else:
                try:
                    with db.connect() as conn:
                        cur = conn.cursor()
                        cur.execute('SELECT identity FROM distances WHERE id = ? LIMIT 1', (iid,))
                        r = cur.fetchone()
                        if r:
                            colour = itol._identity_to_colour(r[0], vmin=0.90, vmax=1.0)
                except Exception:
                    pass  # could not query database for identity; using default colour
                    colour = '#cccccc'
        nov_lines.append(f'{iid},{colour}')

    nov_path = Path(outdir) / 'itol_novelty.itol'
    nov_path.write_text('\n'.join(nov_lines) + '\n')
    logger.info('[ITOL] Wrote novelty iTOL colour strip to %s', nov_path)
    return str(nov_path)


def _novelty_priority_from_metrics(nearest_identity: float, matches_ge_99: int, matches_ge_97: int):
    if nearest_identity < 97.0 and matches_ge_99 <= 1 and matches_ge_97 <= 3:
        return 'HIGH'
    if nearest_identity < 98.5 and matches_ge_99 <= 2 and matches_ge_97 <= 10:
        return 'MEDIUM'
    return 'LOW'


def _novelty_score_from_metrics(nearest_identity: float, matches_ge_99: int, matches_ge_97: int, matches_ge_95: int) -> float:
    distance_component = max(0.0, min(60.0, (100.0 - float(nearest_identity)) * 3.0))
    density_component = 0.0
    if matches_ge_99 <= 1:
        density_component += 20.0
    elif matches_ge_99 <= 5:
        density_component += 10.0
    if matches_ge_97 <= 3:
        density_component += 15.0
    elif matches_ge_97 <= 10:
        density_component += 7.5
    if matches_ge_95 <= 10:
        density_component += 5.0
    return round(min(100.0, distance_component + density_component), 2)


def _normalise_dataset_names(dataset_names) -> list[str]:
    if not dataset_names:
        return []
    if isinstance(dataset_names, str):
        raw = [dataset_names]
    else:
        raw = list(dataset_names)
    names: list[str] = []
    for item in raw:
        for part in str(item).split(','):
            name = part.strip()
            if name and name not in names:
                names.append(name)
    return names


def _infer_baseline_tier_name(dataset: str, explicit: str | None = None) -> str:
    tier = str(explicit or '').strip().lower()
    if tier in {'primary', 'hungate'}:
        tier = 'priority'
    if tier in {'priority', 'secondary'}:
        return tier
    return 'priority' if 'hungate' in str(dataset or '').strip().lower() else 'secondary'


def _discover_tiered_baseline_datasets(
    db,
    *,
    baseline_datasets=None,
    hungate_datasets=None,
    secondary_baseline_datasets=None,
) -> tuple[list[str], list[str], list[str]]:
    priority: list[str] = []
    secondary: list[str] = []

    def add(names, tier: str):
        target = priority if tier == 'priority' else secondary
        for name in _normalise_dataset_names(names):
            if name not in priority and name not in secondary:
                target.append(name)

    role_tiers: dict[str, str] = {}
    if db is not None:
        try:
            grouped = db.get_baseline_datasets_by_tier()
            add(grouped.get('priority', []), 'priority')
            add(grouped.get('secondary', []), 'secondary')
        except Exception:
            try:
                for dataset, metadata in db.get_dataset_roles().items():
                    if metadata.get('role') != 'baseline':
                        continue
                    tier = _infer_baseline_tier_name(dataset, metadata.get('baseline_tier'))
                    role_tiers[dataset] = tier
                    add([dataset], tier)
            except Exception:
                logger.debug('Could not query dataset roles for baseline tiers')
        else:
            try:
                for dataset, metadata in db.get_dataset_roles().items():
                    if metadata.get('role') == 'baseline':
                        role_tiers[dataset] = _infer_baseline_tier_name(dataset, metadata.get('baseline_tier'))
            except Exception:
                logger.debug('Could not query dataset roles for baseline tier names')

    for name in _normalise_dataset_names(baseline_datasets):
        add([name], role_tiers.get(name) or _infer_baseline_tier_name(name))
    add(hungate_datasets, 'priority')
    add(secondary_baseline_datasets, 'secondary')

    combined = list(dict.fromkeys(priority + secondary))
    return list(dict.fromkeys(priority)), list(dict.fromkeys(secondary)), combined


def _query_db_pool(
    db,
    *,
    run_dataset: Optional[str] = None,
    include_datasets=None,
    include_ids=None,
    exclude_datasets=None,
    evidence_pool: str = 'primary',
):
    # The evidence registry is the only source of pool eligibility. Direct writers
    # to ``sequences`` are reconciled conservatively first, so HIGH/MODERATE primary
    # markers are the only sequences that can enter primary pools; Barrnap-derived
    # observations are queried separately and never affect primary novelty scores.
    from branchmanager.marker_evidence import (
        backfill_legacy_marker_evidence, query_marker_evidence_pool, sync_registry_with_sequences,
    )
    backfill_legacy_marker_evidence(db)
    sync_registry_with_sequences(db)
    with db.connect() as conn:
        return query_marker_evidence_pool(
            conn, pool=evidence_pool, run_dataset=run_dataset,
            include_datasets=_normalise_dataset_names(include_datasets) or None,
            include_ids=include_ids,
            exclude_datasets=_normalise_dataset_names(exclude_datasets) or None,
        )


def _write_pool_fasta(records, path: Path) -> Optional[str]:
    if not records:
        return None
    path.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(records, str(path))
    return str(path)


def _build_db_pool_fasta(
    *,
    db,
    outdir: str,
    name: str,
    run_dataset: Optional[str] = None,
    include_datasets=None,
    include_ids=None,
    exclude_datasets=None,
    evidence_pool: str = 'primary',
) -> tuple[Optional[str], str, list[str], int]:
    if db is None:
        return None, 'none', [], 0
    try:
        records, datasets = _query_db_pool(
            db,
            run_dataset=run_dataset,
            include_datasets=include_datasets,
            include_ids=include_ids,
            exclude_datasets=exclude_datasets,
            evidence_pool=evidence_pool,
        )
    except Exception as e:
        raise RuntimeError(f'failed to query the {name} novelty pool from the project database: {e}') from e
    path = _write_pool_fasta(records, Path(outdir) / f'novelty_{name}_pool.fasta')
    if not path:
        return None, 'none', datasets, 0
    return path, name, datasets, len(records)


def _selected_partner_sequence_ids(db) -> list[str]:
    with db.connect() as conn:
        rows = conn.execute(
            "SELECT s.id FROM sequences s "
            "JOIN sequencing_metadata m ON m.id = s.id "
            "LEFT JOIN dataset_roles r ON r.dataset = s.dataset "
            "WHERE s.sequence IS NOT NULL AND s.sequence != '' "
            "AND COALESCE(r.role, '') != 'baseline' "
            "AND (COALESCE(m.selected_for_sequencing, 0) = 1 "
            "OR COALESCE(m.selected_for_wgs, 0) = 1)"
        ).fetchall()
    return sorted({str(row[0]) for row in rows if row[0]})


def _combine_project_target_with_selected_partners(
    target_fasta: str,
    db,
    outdir: str,
) -> tuple[str, str, list[str], int]:
    target_records = dict(read_fasta(target_fasta))
    if db is None:
        return target_fasta, 'target_fasta', [], len(target_records)

    selected_ids = _selected_partner_sequence_ids(db)
    if not selected_ids:
        return target_fasta, 'target_fasta', [], len(target_records)

    selected_records, selected_datasets = _query_db_pool(db, include_ids=selected_ids)
    for sequence_id, sequence in selected_records:
        prior = target_records.get(sequence_id)
        if prior is not None and prior.upper().replace('U', 'T') != sequence.upper().replace('U', 'T'):
            logger.warning(
                '[NOVELTY] Selected project sequence %s differs from the same ID in target FASTA; '
                'using the database sequence in the combined project pool',
                sequence_id,
            )
        target_records[sequence_id] = sequence

    merged_path = _write_pool_fasta(
        sorted(target_records.items()),
        Path(outdir) / 'novelty_project_collection_pool.fasta',
    )
    if not merged_path:
        return target_fasta, 'target_fasta', [], len(target_records)
    return (
        merged_path,
        'target_fasta+selected_project_isolates',
        selected_datasets,
        len(target_records),
    )


def _run_nearest_search(
    input_fasta: str,
    db_fasta: str,
    matches_path: Path,
    threads: Optional[int] = None,
    *,
    exclude_self: bool = False,
):
    """Find one best hit/query using 90→50% identity and 95→80% coverage tiers.

    Full-coverage (>=95%) matches are preferred. Only queries without a
    full-coverage hit enter the partial-coverage pass (>=80%). Partial hits
    are reported but do not satisfy the existing >=95% support criterion
    used by downstream novelty scoring.
    """
    import csv
    import tempfile

    identity_ladder = NEAREST_IDENTITY_LADDER
    coverage_tiers = (
        ("FULL", MIN_NEAREST_QUERY_COVERAGE),
        ("PARTIAL", PARTIAL_NEAREST_QUERY_COVERAGE),
    )
    matches_path = Path(matches_path)
    matches_path.parent.mkdir(parents=True, exist_ok=True)
    query_records = list(read_fasta(input_fasta))
    unresolved = dict(query_records)
    best = {}
    summary_rows = []
    total_raw_hits = 0
    thread_flag = f" --threads {int(threads)}" if threads and int(threads) > 0 else ""

    with tempfile.TemporaryDirectory(
        prefix=f".{matches_path.stem}_ladder_", dir=str(matches_path.parent)
    ) as temp_dir:
        root = Path(temp_dir)
        for coverage_label, coverage_floor in coverage_tiers:
            tier_unresolved = dict(unresolved)
            for round_number, identity_floor in enumerate(identity_ladder, 1):
                if not tier_unresolved:
                    break
                query_path = root / f"queries_{coverage_label}_{identity_floor}.fasta"
                raw_path = root / f"raw_{coverage_label}_{identity_floor}.tsv"
                write_fasta(list(tier_unresolved.items()), str(query_path))
                cmd = (
                    f"vsearch --usearch_global {shlex.quote(str(query_path))} "
                    f"--db {shlex.quote(str(db_fasta))} "
                    f"--id {identity_floor / 100.0:.2f} --strand both "
                    f"--userout {shlex.quote(str(raw_path))} "
                    f"--userfields {VSEARCH_USERFIELDS} "
                    f"--maxaccepts 32 --maxrejects 256 --maxhits 32 "
                    f"--query_cov {coverage_floor / 100.0:.2f}{thread_flag}"
                )
                logger.info(
                    "[NOVELTY][LADDER] coverage=%s(%.1f%%) identity_floor=%d%% "
                    "round=%d unresolved_queries=%d",
                    coverage_label, coverage_floor, identity_floor,
                    round_number, len(tier_unresolved),
                )
                run_vsearch_audited(
                    cmd, audit_path=matches_path.parent / 'vsearch_search_audit.tsv',
                    stage='B', purpose=f'nearest_hit_{coverage_label.lower()}',
                    query_fasta=query_path, reference_fasta=db_fasta, runner=run_cmd,
                )
                round_best = {}
                raw_hit_count = 0
                if raw_path.exists():
                    with open(raw_path) as handle:
                        for line in handle:
                            parts = line.rstrip("\n").split("\t")
                            if len(parts) < 3 or parts[0] not in tier_unresolved:
                                continue
                            if exclude_self and parts[0] == parts[1]:
                                continue
                            try:
                                pct = float(parts[2])
                            except (TypeError, ValueError):
                                continue
                            qcov, alnlen = _alignment_evidence(parts)
                            if qcov is None or qcov < coverage_floor:
                                continue
                            raw_hit_count += 1
                            candidate = (pct, parts[1], qcov, alnlen,
                                         identity_floor, round_number, coverage_label,
                                         parts[:14])
                            rank = (pct, qcov, alnlen if alnlen is not None else -1, parts[1])
                            current = round_best.get(parts[0])
                            current_rank = (
                                current[0], current[2],
                                current[3] if current[3] is not None else -1, current[1]
                            ) if current else None
                            if current_rank is None or rank > current_rank:
                                round_best[parts[0]] = candidate
                total_raw_hits += raw_hit_count
                for qid, hit in round_best.items():
                    best[qid] = hit
                tier_unresolved = {
                    qid: seq for qid, seq in tier_unresolved.items()
                    if qid not in round_best
                }
                summary_rows.append({
                    "coverage_tier": coverage_label,
                    "minimum_query_coverage": f"{coverage_floor:.1f}",
                    "identity_floor": identity_floor,
                    "round": round_number,
                    "queries_searched": len(tier_unresolved) + len(round_best),
                    "queries_resolved": len(round_best),
                    "queries_remaining_in_tier": len(tier_unresolved),
                    "raw_qualifying_hits": raw_hit_count,
                })
                logger.info(
                    "[NOVELTY][LADDER] coverage=%s identity_floor=%d%% "
                    "raw_hits=%d resolved=%d remaining=%d",
                    coverage_label, identity_floor, raw_hit_count,
                    len(round_best), len(tier_unresolved),
                )
            unresolved = tier_unresolved

    fields = [
        "query", "target", "id", "alnlen", "mism", "gaps", "qlo", "qhi",
        "tlo", "thi", "ql", "tl", "qcov", "tcov",
        "search_threshold", "search_round", "coverage_tier", "match_status",
    ]
    with open(matches_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        for qid, _seq in query_records:
            hit = best.get(qid)
            if hit is None:
                writer.writerow([qid, "NA", "NA", * (["NA"] * 11),
                                 "NA", "NA", "NONE", "NO_MATCH"])
            else:
                pct, target, qcov, alnlen, threshold, round_no, tier, raw_fields = hit
                writer.writerow(
                    raw_fields + [
                        threshold, round_no, tier,
                        "MATCH" if tier == "FULL" else "PARTIAL_COVERAGE",
                    ]
                )

    summary_name = (
        "novelty_reference_search_summary.tsv"
        if "reference" in matches_path.name
        else f"{matches_path.stem}_search_summary.tsv"
    )
    summary_path = matches_path.with_name(summary_name)
    with open(summary_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]) if summary_rows else [
            "coverage_tier", "minimum_query_coverage", "identity_floor", "round",
            "queries_searched", "queries_resolved", "queries_remaining_in_tier",
            "raw_qualifying_hits",
        ], delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    full = sum(hit[6] == "FULL" for hit in best.values())
    partial = sum(hit[6] == "PARTIAL" for hit in best.values())
    logger.info(
        "[NOVELTY][LADDER] Complete: queries=%d full=%d partial=%d unmatched=%d "
        "raw_hits=%d; matches=%s; summary=%s",
        len(query_records), full, partial, len(query_records) - len(best),
        total_raw_hits, matches_path, summary_path,
    )
    return best


def _run_density_search(
    input_fasta: str,
    db_fasta: str,
    matches_path: Path,
    threads: Optional[int] = None,
    *,
    exclude_self: bool = False,
    coverage_floor: float = MIN_NEAREST_QUERY_COVERAGE,
    retain_matches: bool = True,
):
    thread_flag = f" --threads {int(threads)}" if threads and int(threads) > 0 else ""
    # Density is a bounded prioritisation signal. Counts are explicitly
    # reported as lower bounds; nearest-hit evidence remains uncapped separately.
    cmd = (
        f"vsearch --usearch_global {shlex.quote(str(input_fasta))} "
        f"--db {shlex.quote(str(db_fasta))} "
        f"--id 0.95 "
        f"--strand both "
        f"--userout {shlex.quote(str(matches_path))} "
        f"--userfields {VSEARCH_USERFIELDS} "
        f"--maxaccepts {DENSITY_HIT_CAP} --maxrejects 256 --maxhits {DENSITY_HIT_CAP} "
        f"--query_cov {coverage_floor / 100.0:.2f}{thread_flag}"
    )
    run_vsearch_audited(
        cmd, audit_path=Path(matches_path).parent / 'vsearch_search_audit.tsv',
        stage='B', purpose=f'density_qcov{int(coverage_floor)}', query_fasta=input_fasta,
        reference_fasta=db_fasta, runner=run_cmd,
    )

    density_counts = {}
    if matches_path.exists():
        with open(matches_path) as fh:
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 3:
                    continue
                qid = parts[0]
                if exclude_self and parts[1] == qid:
                    continue
                try:
                    pct = float(parts[2])
                except Exception:
                    continue
                bucket = density_counts.setdefault(qid, {'ge_99': 0, 'ge_97': 0, 'ge_95': 0})
                if pct >= 95.0:
                    bucket['ge_95'] += 1
                if pct >= 97.0:
                    bucket['ge_97'] += 1
                if pct >= 99.0:
                    bucket['ge_99'] += 1
    if not retain_matches:
        Path(matches_path).unlink(missing_ok=True)
    return density_counts


def _run_density_search_coverage_pair(
    input_fasta: str,
    db_fasta: str,
    matches_path: Path,
    threads: Optional[int] = None,
    *,
    exclude_self: bool = False,
    retain_95_matches: bool = True,
    batch_size: int = 500,
    hit_cap: int = DENSITY_HIT_CAP,
) -> tuple[dict, dict]:
    """Calculate bounded qcov95/qcov80 density in restartable query batches.

    Each batch is searched against the *entire* reference pool, but VSEARCH
    stops after ``hit_cap`` accepted hits per query. Counts are therefore
    lower bounds, not exhaustive neighbour counts.
    A batch is only marked complete after VSEARCH exits successfully and its
    output has been audited. On restart,
    completed batches are checksum-verified and reused; an interrupted batch
    is rerun. The final qcov95 compatibility file is rebuilt atomically from
    the verified batch outputs.

    The identity/coverage criteria remain 95% identity, both strands, and an
    80% query-coverage search floor. Per-hit query coverage splits qcov80 and
    qcov95 counts. The accepted-hit cap bounds runtime and output size.

    Batching trades a small amount of repeated reference-database setup for
    bounded intermediate files and recovery without repeating completed
    batches. The default (500 queries) gives two batches for a 692-query
    cohort; callers/tests may override it.
    """
    import hashlib
    import json

    started = time.monotonic()
    matches_path = Path(matches_path)
    matches_path.parent.mkdir(parents=True, exist_ok=True)
    query_records = list(read_fasta(input_fasta))
    if batch_size < 1:
        raise ValueError('density batch_size must be at least 1')
    if int(hit_cap) < 1:
        raise ValueError('density hit_cap must be at least 1')

    def sha256(path: str | Path) -> str:
        digest = hashlib.sha256()
        with open(path, 'rb') as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b''):
                digest.update(block)
        return digest.hexdigest()

    query_sha = sha256(input_fasta)
    reference_sha = sha256(db_fasta)
    settings = {
        'schema_version': 1,
        'query_sha256': query_sha,
        'reference_sha256': reference_sha,
        'identity': 0.95,
        'strand': 'both',
        'query_coverage_floor': 0.80,
        'maxaccepts': int(hit_cap),
        'maxrejects': 256,
        'maxhits': int(hit_cap),
        'batch_size': int(batch_size),
        'query_count': len(query_records),
        'density_count_semantics': 'LOWER_BOUND',
    }
    run_key = hashlib.sha256(json.dumps(settings, sort_keys=True).encode()).hexdigest()[:12]
    checkpoint_dir = matches_path.parent / f'.{matches_path.stem}_density_checkpoints_{run_key}'
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    settings_path = checkpoint_dir / 'run.json'
    if settings_path.exists():
        existing = json.loads(settings_path.read_text(encoding='utf-8'))
        if existing != settings:
            raise RuntimeError(
                f'Density checkpoint settings do not match current inputs: {checkpoint_dir}. '
                'Do not reuse these checkpoints for a different query/reference set.'
            )
    else:
        temporary_settings = settings_path.with_suffix('.json.tmp')
        temporary_settings.write_text(json.dumps(settings, indent=2, sort_keys=True) + '\n', encoding='utf-8')
        os.replace(temporary_settings, settings_path)

    thread_flag = f" --threads {int(threads)}" if threads and int(threads) > 0 else ''
    batch_outputs: list[Path] = []
    completed_batches = 0
    total_batches = (len(query_records) + batch_size - 1) // batch_size
    logger.info(
        '[NOVELTY][DENSITY] resumable exhaustive search START qcov=80/95 queries=%d batches=%d '
        'batch_size=%d db=%s',
        len(query_records), total_batches, batch_size, db_fasta,
    )

    for batch_index, offset in enumerate(range(0, len(query_records), batch_size), start=1):
        batch_records = query_records[offset:offset + batch_size]
        batch_query = checkpoint_dir / f'queries_{batch_index:04d}.fasta'
        batch_output = checkpoint_dir / f'hits_{batch_index:04d}.tsv'
        batch_state = checkpoint_dir / f'batch_{batch_index:04d}.json'
        # Re-write the small deterministic query chunk so its content is
        # explicit and can be independently audited.
        write_fasta(batch_records, str(batch_query))
        batch_query_sha = sha256(batch_query)
        valid_checkpoint = False
        if batch_state.is_file() and batch_output.is_file():
            try:
                state = json.loads(batch_state.read_text(encoding='utf-8'))
                valid_checkpoint = (
                    state.get('query_sha256') == batch_query_sha
                    and state.get('reference_sha256') == reference_sha
                    and state.get('settings') == settings
                    and state.get('output_sha256') == sha256(batch_output)
                )
            except (OSError, ValueError, TypeError):
                valid_checkpoint = False
        if valid_checkpoint:
            completed_batches += 1
            batch_outputs.append(batch_output)
            logger.info(
                '[NOVELTY][DENSITY] batch %d/%d reused from verified checkpoint (%d queries)',
                batch_index, total_batches, len(batch_records),
            )
            continue

        # Never allow an old/incomplete output to be mistaken for a successful
        # run. The completion record is written only after VSEARCH and the audit
        # wrapper both return successfully; without it, any output is discarded
        # and the batch is rerun on the next invocation.
        batch_state.unlink(missing_ok=True)
        batch_output.unlink(missing_ok=True)
        cmd = (
            f"vsearch --usearch_global {shlex.quote(str(batch_query))} "
            f"--db {shlex.quote(str(db_fasta))} --id 0.95 --strand both "
            f"--userout {shlex.quote(str(batch_output))} "
            f"--userfields {VSEARCH_USERFIELDS} "
            f"--maxaccepts {int(hit_cap)} --maxrejects 256 --maxhits {int(hit_cap)} "
            f"--query_cov 0.80{thread_flag}"
        )
        logger.info(
            '[NOVELTY][DENSITY] batch %d/%d START queries=%d',
            batch_index, total_batches, len(batch_records),
        )
        run_vsearch_audited(
            cmd, audit_path=matches_path.parent / 'vsearch_search_audit.tsv',
            stage='B', purpose='density_qcov80_and_95_resumable_batch',
            query_fasta=batch_query, reference_fasta=db_fasta, runner=run_cmd,
        )
        if not batch_output.exists():
            # A successful zero-hit VSEARCH invocation may leave no userout
            # file on some builds; represent that valid result as an empty TSV.
            batch_output.touch()
        with batch_output.open(encoding='utf-8', errors='replace') as completed_output:
            raw_hit_count = sum(1 for line in completed_output if line.strip())
        state = {
            'batch': batch_index,
            'query_count': len(batch_records),
            'query_sha256': batch_query_sha,
            'reference_sha256': reference_sha,
            'settings': settings,
            'output_sha256': sha256(batch_output),
            'raw_hit_count': raw_hit_count,
        }
        temporary_state = batch_state.with_suffix('.json.tmp')
        temporary_state.write_text(json.dumps(state, indent=2, sort_keys=True) + '\n', encoding='utf-8')
        os.replace(temporary_state, batch_state)
        completed_batches += 1
        batch_outputs.append(batch_output)
        logger.info(
            '[NOVELTY][DENSITY] batch %d/%d COMPLETE queries=%d raw_hits=%d',
            batch_index, total_batches, len(batch_records), state['raw_hit_count'],
        )

    counts_80, counts_95 = {}, {}
    aggregate_tmp = matches_path.with_name(matches_path.name + '.partial')
    aggregate_tmp.unlink(missing_ok=True)
    output_handle = aggregate_tmp.open('w', encoding='utf-8') if retain_95_matches else None
    try:
        for batch_output in batch_outputs:
            with batch_output.open(encoding='utf-8', errors='replace') as source:
                for line in source:
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) < 14:
                        continue
                    qid, target = parts[0], parts[1]
                    if exclude_self and target == qid:
                        continue
                    try:
                        identity, qcov = float(parts[2]), float(parts[12])
                    except (TypeError, ValueError):
                        continue
                    if identity < 95.0 or qcov < 80.0:
                        continue
                    bucket80 = counts_80.setdefault(qid, {'ge_99': 0, 'ge_97': 0, 'ge_95': 0})
                    bucket80['ge_95'] += 1
                    if identity >= 97.0:
                        bucket80['ge_97'] += 1
                    if identity >= 99.0:
                        bucket80['ge_99'] += 1
                    if qcov >= 95.0:
                        bucket95 = counts_95.setdefault(qid, {'ge_99': 0, 'ge_97': 0, 'ge_95': 0})
                        bucket95['ge_95'] += 1
                        if identity >= 97.0:
                            bucket95['ge_97'] += 1
                        if identity >= 99.0:
                            bucket95['ge_99'] += 1
                        if output_handle is not None:
                            output_handle.write(line)
        if output_handle is not None:
            output_handle.close()
            output_handle = None
            os.replace(aggregate_tmp, matches_path)
    except BaseException:
        if output_handle is not None:
            output_handle.close()
        aggregate_tmp.unlink(missing_ok=True)
        raise

    logger.info(
        '[NOVELTY][DENSITY] resumable exhaustive search COMPLETE batches=%d/%d '
        'qcov80_queries=%d qcov95_queries=%d elapsed=%.2fs checkpoints=%s',
        completed_batches, total_batches, len(counts_80), len(counts_95),
        time.monotonic() - started, checkpoint_dir,
    )
    return counts_95, counts_80


def _crowding_from_counts(counts: dict) -> str:
    if counts['ge_99'] <= 1 and counts['ge_97'] <= 1:
        return 'isolated'
    if counts['ge_97'] <= 3:
        return 'sparse'
    if counts['ge_97'] <= 10:
        return 'moderate'
    return 'crowded'


def _metrics_from_nearest_and_density(qids, nearest_hits, density_counts, *, source: str, id_threshold: float):
    novel_threshold_pct = id_threshold * 100.0
    rows = {}
    for qid in qids:
        hit = nearest_hits.get(qid, (0.0, 'None', None, None))
        pct_id, hit_id = hit[0], hit[1]
        query_coverage = hit[2] if len(hit) > 2 else None
        has_supported_coverage = (
            query_coverage is not None
            and query_coverage >= MIN_NEAREST_QUERY_COVERAGE
        )
        decision_identity = pct_id if has_supported_coverage else 0.0
        counts = density_counts.get(qid, {'ge_99': 0, 'ge_97': 0, 'ge_95': 0})
        crowding = _crowding_from_counts(counts)
        score = _novelty_score_from_metrics(
            decision_identity, counts['ge_99'], counts['ge_97'], counts['ge_95'],
        )
        priority = _novelty_priority_from_metrics(decision_identity, counts['ge_99'], counts['ge_97'])
        rows[qid] = {
            'nearest_identity': pct_id,
            'nearest_hit': hit_id,
            'nearest_query_coverage': query_coverage,
            'nearest_alignment_length': hit[3] if len(hit) > 3 else None,
            'novel': str(not has_supported_coverage or decision_identity < novel_threshold_pct),
            'matches_ge_99': counts['ge_99'],
            'matches_ge_97': counts['ge_97'],
            'matches_ge_95': counts['ge_95'],
            'novelty_score': score,
            'crowding': crowding,
            'sequencing_priority': priority,
            'density_source': source,
            'density_count_method': f'LOWER_BOUND_MAX_{DENSITY_HIT_CAP}_HITS',
        }
    return rows


def _empty_pool_metrics(qids, *, source: str = 'none'):
    return {
        qid: {
            'nearest_identity': 0.0,
            'nearest_hit': 'None',
            'nearest_query_coverage': None,
            'nearest_alignment_length': None,
            'novel': 'NA',
            'matches_ge_99': 0,
            'matches_ge_97': 0,
            'matches_ge_95': 0,
            'novelty_score': 0.0,
            'crowding': 'unknown',
            'sequencing_priority': 'NA',
            'density_source': source,
            'density_count_method': 'NOT_AVAILABLE',
        }
        for qid in qids
    }


def _empty_sequencing_context(qids, pangenome_target: int = DEFAULT_PANGENOME_TARGET):
    target_count = max(1, int(pangenome_target or DEFAULT_PANGENOME_TARGET))
    return {
        qid: {
            'partner_id': 'NA',
            'selected_for_sequencing': 'NA',
            'already_sequenced': 'NA',
            'selected_for_wgs': 'NA',
            'nearest_genome_hit': 'None',
            'nearest_genome_identity': 0.0,
            'closest_project_isolate': 'None',
            'closest_project_isolate_identity': 0.0,
            'project_isolate_only_match': 'False',
            'genome_ge_99': 0,
            'genome_ge_97': 0,
            'genome_ge_95': 0,
            'related_genome_clade': 'False',
            'adjusted_novelty_score': 'NA',
            'adjusted_priority': 'NA',
            'source': 'none',
            'genome_same_species_committed': 0,
            'genome_same_species_hungate': 0,
            'genome_same_species_secondary_baseline': 0,
            'genome_same_species_available': 0,
            'genome_same_species_selected': 0,
            'genome_same_species_pending': 0,
            'pangenome_target': 'NA',
            'pangenome_gap': 'NA',
            'pangenome_eligibility_status': 'DEFERRED_TO_TREE_VALIDATED_SELECTION_SETS',
        }
        for qid in qids
    }



def _build_combined_baseline_pool_with_provenance(
    *, db, outdir: str, hungate_datasets, secondary_datasets,
) -> tuple[Optional[str], dict[str, str], list[str], int]:
    """Build one cultured-baseline FASTA and retain each target's tier provenance."""
    datasets_by_tier = {
        'HUNGATE': _normalise_dataset_names(hungate_datasets),
        'SECONDARY': _normalise_dataset_names(secondary_datasets),
    }
    all_datasets = list(dict.fromkeys(
        datasets_by_tier['HUNGATE'] + datasets_by_tier['SECONDARY']
    ))
    if db is None or not all_datasets:
        return None, {}, [], 0

    from branchmanager.marker_evidence import (
        backfill_legacy_marker_evidence, sync_registry_with_sequences,
    )
    backfill_legacy_marker_evidence(db)
    sync_registry_with_sequences(db)
    placeholders = ','.join('?' for _ in all_datasets)
    with db.connect() as conn:
        rows = conn.execute(
            "SELECT sequence_id, sequence, dataset FROM marker_evidence "
            "WHERE eligible_primary=1 AND active_for_selection=1 "
            "AND evidence_tier IN ('HIGH','MODERATE') "
            f"AND dataset IN ({placeholders}) "
            "ORDER BY sequence_id, created_at DESC, evidence_id DESC",
            tuple(all_datasets),
        ).fetchall()

    tier_by_dataset = {name: tier for tier, names in datasets_by_tier.items() for name in names}
    records_by_id = {}
    provenance = {}
    for sequence_id, sequence, dataset in rows:
        if not sequence_id or not sequence or not dataset:
            continue
        sid = str(sequence_id)
        # The evidence registry can contain historical versions. Keep the first
        # (latest) record per sequence ID, matching query_marker_evidence_pool.
        if sid in records_by_id:
            continue
        records_by_id[sid] = str(sequence)
        provenance[sid] = tier_by_dataset.get(str(dataset), 'SECONDARY')

    pool_path = _write_pool_fasta(
        sorted(records_by_id.items()),
        Path(outdir) / 'novelty_combined_baseline_pool.fasta',
    )
    return pool_path, provenance, all_datasets, len(records_by_id)


def _run_nearest_search_by_provenance(
    input_fasta: str, db_fasta: str, matches_prefix: Path, target_provenance: dict[str, str],
    threads: Optional[int] = None,
) -> dict[str, dict]:
    """Run one ladder search against a combined pool, retaining best hit per tier/query."""
    import csv
    import tempfile

    query_records = list(read_fasta(input_fasta))
    query_by_id = dict(query_records)
    groups = ('HUNGATE', 'SECONDARY')
    active_groups = tuple(group for group in groups if group in set(target_provenance.values()))
    best = {group: {} for group in groups}
    unresolved = {(qid, group) for qid in query_by_id for group in active_groups}
    summary_rows = []
    thread_flag = f" --threads {int(threads)}" if threads and int(threads) > 0 else ''

    with tempfile.TemporaryDirectory(prefix='.baseline_ladder_', dir=str(matches_prefix.parent)) as temp_dir:
        root = Path(temp_dir)
        for coverage_label, coverage_floor in (
            ('FULL', MIN_NEAREST_QUERY_COVERAGE),
            ('PARTIAL', PARTIAL_NEAREST_QUERY_COVERAGE),
        ):
            for round_number, identity_floor in enumerate(NEAREST_IDENTITY_LADDER, 1):
                active_qids = sorted({qid for qid, group in unresolved})
                if not active_qids:
                    break
                query_path = root / f'queries_{coverage_label}_{identity_floor}.fasta'
                raw_path = root / f'raw_{coverage_label}_{identity_floor}.tsv'
                write_fasta([(qid, query_by_id[qid]) for qid in active_qids], str(query_path))
                cmd = (
                    f"vsearch --usearch_global {shlex.quote(str(query_path))} "
                    f"--db {shlex.quote(str(db_fasta))} --id {identity_floor / 100.0:.2f} "
                    f"--strand both --userout {shlex.quote(str(raw_path))} "
                    f"--userfields {VSEARCH_USERFIELDS} "
                    f"--maxaccepts 32 --maxrejects 256 --maxhits 32 "
                    f"--query_cov {coverage_floor / 100.0:.2f}{thread_flag}"
                )
                run_vsearch_audited(
                    cmd, audit_path=matches_prefix.parent / 'vsearch_search_audit.tsv',
                    stage='B', purpose=f'baseline_nearest_{coverage_label.lower()}',
                    query_fasta=query_path, reference_fasta=db_fasta, runner=run_cmd,
                )
                round_best = {}
                raw_hits = 0
                if raw_path.exists():
                    with open(raw_path) as handle:
                        for line in handle:
                            parts = line.rstrip('\n').split('\t')
                            if len(parts) < 14 or parts[0] not in query_by_id:
                                continue
                            group = target_provenance.get(parts[1])
                            if group not in groups or (parts[0], group) not in unresolved:
                                continue
                            try:
                                identity = float(parts[2])
                            except (TypeError, ValueError):
                                continue
                            qcov, alnlen = _alignment_evidence(parts)
                            if qcov is None or qcov < coverage_floor:
                                continue
                            raw_hits += 1
                            rank = (identity, qcov, alnlen or -1, parts[1])
                            key = (parts[0], group)
                            previous = round_best.get(key)
                            if previous is None or rank > previous[0]:
                                round_best[key] = (rank, (identity, parts[1], qcov, alnlen,
                                    identity_floor, round_number, coverage_label, parts[:14]))
                for (qid, group), (_rank, hit) in round_best.items():
                    best[group][qid] = hit
                    unresolved.discard((qid, group))
                summary_rows.append({
                    'coverage_tier': coverage_label, 'minimum_query_coverage': coverage_floor,
                    'identity_floor': identity_floor, 'round': round_number,
                    'queries_searched': len(active_qids), 'queries_resolved': len(round_best),
                    'query_group_pairs_remaining': len(unresolved), 'raw_qualifying_hits': raw_hits,
                })

    fields = ['query', 'target', 'id', 'alnlen', 'mism', 'gaps', 'qlo', 'qhi', 'tlo', 'thi',
              'ql', 'tl', 'qcov', 'tcov', 'search_threshold', 'search_round', 'coverage_tier', 'match_status']
    output_by_group = {}
    for group in groups:
        path = matches_prefix.with_name(
            'novelty_hungate_baseline_nearest_matches.tsv' if group == 'HUNGATE'
            else 'novelty_secondary_baseline_nearest_matches.tsv'
        )
        output_by_group[group] = path
        with open(path, 'w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(fields)
            for qid, _sequence in query_records:
                hit = best[group].get(qid)
                if hit is None:
                    writer.writerow([qid, 'NA', 'NA'] + ['NA'] * 11 + ['NA', 'NA', 'NONE', 'NO_MATCH'])
                else:
                    identity, target, qcov, alnlen, threshold, round_no, tier, raw_fields = hit
                    writer.writerow(raw_fields + [threshold, round_no, tier,
                        'MATCH' if tier == 'FULL' else 'PARTIAL_COVERAGE'])

    summary_fields = list(summary_rows[0]) if summary_rows else [
        'coverage_tier','minimum_query_coverage','identity_floor','round','queries_searched',
        'queries_resolved','query_group_pairs_remaining','raw_qualifying_hits',
    ]
    summary_path = matches_prefix.with_name('novelty_baseline_nearest_search_summary.tsv')
    with open(summary_path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter='\t')
        writer.writeheader(); writer.writerows(summary_rows)
    # Preserve the historical per-pool summary filenames for downstream scripts.
    for group in groups:
        compatibility_path = output_by_group[group].with_name(output_by_group[group].stem + '_search_summary.tsv')
        with open(compatibility_path, 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter='\t')
            writer.writeheader(); writer.writerows(summary_rows)
    combined_summary = matches_prefix.with_name('novelty_baseline_nearest_matches_search_summary.tsv')
    with open(combined_summary, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter='\t')
        writer.writeheader(); writer.writerows(summary_rows)
    return best


def _run_combined_baseline_density(
    input_fasta: str, db_fasta: str, outdir: str, target_provenance: dict[str, str],
    threads: Optional[int] = None,
) -> dict[str, dict]:
    """Run bounded 95% and 80% density comparisons against the combined baseline pool.

    At most DENSITY_HIT_CAP accepted hits are counted per query; all resulting
    density values are lower bounds.
    """
    results = {}
    thread_flag = f" --threads {int(threads)}" if threads and int(threads) > 0 else ''
    # Bounded qcov95/qcov80 lower-bound counts are derived from one qcov80
    # superset scan. Provenance is retained while streaming; only historical
    # qcov95 match files are written to disk.
    started = time.monotonic()
    raw_path = Path(outdir) / 'novelty_baseline_density_all_qcov80_raw.tsv'
    cmd = (
        f"vsearch --usearch_global {shlex.quote(str(input_fasta))} "
        f"--db {shlex.quote(str(db_fasta))} --id 0.95 --strand both "
        f"--userout {shlex.quote(str(raw_path))} --userfields {VSEARCH_USERFIELDS} "
        f"--maxaccepts {DENSITY_HIT_CAP} --maxrejects 256 --maxhits {DENSITY_HIT_CAP} "
        f"--query_cov 0.80{thread_flag}"
    )
    raw_path.unlink(missing_ok=True)
    logger.info(
        '[NOVELTY][DENSITY][BASELINE] START bounded combined qcov=80/95 hit_cap=%d pool=%s',
        DENSITY_HIT_CAP, db_fasta,
    )
    run_vsearch_audited(
        cmd, audit_path=Path(outdir) / 'vsearch_search_audit.tsv',
        stage='B', purpose='baseline_density_qcov80_and_95_bounded_lower_bound',
        query_fasta=input_fasta, reference_fasta=db_fasta, runner=run_cmd,
    )
    counts = {coverage: {group: {} for group in ('HUNGATE', 'SECONDARY')} for coverage in (95, 80)}
    output_paths = {
        group: Path(outdir) / f'novelty_{label}_density_matches.tsv'
        for group, label in (('HUNGATE', 'hungate_baseline'), ('SECONDARY', 'secondary_baseline'), ('ALL', 'baseline'))
    }
    handles = {group: path.open('w', encoding='utf-8') for group, path in output_paths.items()}
    combined_qcov95_handle = (Path(outdir) / 'novelty_baseline_density_all_qcov95.tsv').open('w', encoding='utf-8')
    try:
        with raw_path.open(encoding='utf-8', errors='replace') as source:
            for line in source:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 14:
                    continue
                group = target_provenance.get(parts[1])
                if group not in {'HUNGATE', 'SECONDARY'}:
                    continue
                try:
                    identity, qcov = float(parts[2]), float(parts[12])
                except (TypeError, ValueError):
                    continue
                for coverage in (80, 95):
                    if qcov < coverage:
                        continue
                    bucket = counts[coverage][group].setdefault(parts[0], {'ge_99': 0, 'ge_97': 0, 'ge_95': 0})
                    if identity >= 95: bucket['ge_95'] += 1
                    if identity >= 97: bucket['ge_97'] += 1
                    if identity >= 99: bucket['ge_99'] += 1
                if qcov >= 95:
                    handles[group].write(line)
                    handles['ALL'].write(line)
                    combined_qcov95_handle.write(line)
    finally:
        for handle in handles.values():
            handle.close()
        combined_qcov95_handle.close()
        raw_path.unlink(missing_ok=True)
    results[95] = counts[95]
    results[80] = counts[80]
    logger.info('[NOVELTY][DENSITY][BASELINE] COMPLETE qcov80_queries=%d qcov95_queries=%d elapsed=%.2fs',
                sum(len(v) for v in counts[80].values()), sum(len(v) for v in counts[95].values()),
                time.monotonic() - started)
    return results

def _compute_db_pool_metrics(
    *,
    input_fasta: str,
    outdir: str,
    qids,
    db,
    run_dataset: Optional[str],
    datasets,
    pool_name: str,
    source_prefix: str,
    threads: Optional[int],
):
    dataset_names = _normalise_dataset_names(datasets)
    if not dataset_names or db is None:
        return None, _empty_pool_metrics(qids, source='none'), 'none', [], 0
    pool_path, _pool_source, datasets_found, count = _build_db_pool_fasta(
        db=db,
        outdir=outdir,
        name=pool_name,
        run_dataset=run_dataset,
        include_datasets=dataset_names,
    )
    if not pool_path:
        return None, _empty_pool_metrics(qids, source='none'), 'none', datasets_found, count
    source = source_prefix + ':' + ','.join(datasets_found or dataset_names)
    logger.info(
        '[NOVELTY] %s pool: %d sequences from dataset(s): %s',
        source_prefix.replace('_', ' '), count, datasets_found or dataset_names,
    )
    nearest = _run_nearest_search(
        input_fasta,
        pool_path,
        Path(outdir) / f'novelty_{pool_name}_nearest_matches.tsv',
        threads=threads,
    )
    density = _run_density_search(
        input_fasta,
        pool_path,
        Path(outdir) / f'novelty_{pool_name}_density_matches.tsv',
        threads=threads,
    )
    metrics = _metrics_from_nearest_and_density(
        qids,
        nearest,
        density,
        source=source,
        id_threshold=0.97,
    )
    return pool_path, metrics, source, datasets_found, count


def _compute_sequencing_context(
    input_fasta: str,
    outdir: str,
    db,
    qids,
    primary_metrics,
    threads: Optional[int] = None,
    pangenome_target: int = DEFAULT_PANGENOME_TARGET,
    project_metrics=None,
    reference_metrics=None,
    cultured_metrics=None,
):
    if db is None:
        return _empty_sequencing_context(qids, pangenome_target)

    from branchmanager.marker_evidence import sync_registry_with_sequences
    sync_registry_with_sequences(db)
    project_metrics = project_metrics or {}
    reference_metrics = reference_metrics or {}
    cultured_metrics = cultured_metrics or {}

    try:
        metadata = db.get_sequencing_metadata_for_ids(qids)
    except Exception:
        metadata = {}

    genome_records = []
    genome_metadata = {}
    collection_metadata = {}
    query_taxonomy = {}
    try:
        dataset_roles = db.get_dataset_roles()
        baseline_genome_datasets = {
            dataset for dataset, values in dataset_roles.items()
            if values.get('role') == 'baseline' and values.get('genomes_available')
        }
        with db.connect() as conn:
            cur = conn.cursor()
            cur.execute(
                "SELECT s.id, s.sequence, s.dataset, s.baseline_tier, t.taxonomy, "
                "COALESCE(m.selected_for_wgs, 0), COALESCE(m.selected_for_sequencing, 0) "
                "FROM sequences s "
                "LEFT JOIN taxonomy t ON s.id = t.id AND (t.dataset = s.dataset OR t.dataset IS NULL) "
                "LEFT JOIN sequencing_metadata m ON s.id = m.id "
                "WHERE s.sequence IS NOT NULL AND s.sequence != '' "
                "AND s.id IN (SELECT sequence_id FROM marker_evidence WHERE eligible_primary=1 "
                "AND active_for_selection=1 AND evidence_tier IN ('HIGH','MODERATE'))"
            )
            rows = cur.fetchall()
        by_id = {}
        for sid, sequence, dataset, baseline_tier, taxonomy, available, selected in rows:
            sid = str(sid)
            entry = by_id.setdefault(sid, {
                'sequence': str(sequence),
                'dataset': str(dataset or ''),
                'baseline_tier': str(baseline_tier or ''),
                'taxonomy': '',
                'available': bool(available),
                'selected': bool(selected),
            })
            if taxonomy and not entry['taxonomy']:
                entry['taxonomy'] = str(taxonomy)
        # Only taxonomically validated records may contribute to species-level
        # pangenome counts. Require confident evidence in every configured DB
        # and exact normalised species agreement across those assignments.
        confident_taxonomy_ids = set()
        try:
            classification_evidence = db.get_classification_evidence(set(by_id))
            evidence_by_id_db = {}
            for (sid, ref_db), evidence_row in classification_evidence.items():
                evidence_by_id_db.setdefault(str(sid), {})[str(ref_db)] = evidence_row
            alt_ref_dbs = list(db.get_alt_ref_dbs() or [])
            alt_taxonomies = db.get_taxonomy_alt_for_ids(set(by_id)) if alt_ref_dbs else {}
            for sid, evidence_rows in evidence_by_id_db.items():
                if not evidence_rows or any(
                    str(item.get('assignment_status') or '').upper() != 'CONFIDENT'
                    for item in evidence_rows.values()
                ):
                    continue
                primary_species = str(parse_taxon_string(by_id[sid].get('taxonomy', '')).get('s') or '').strip()
                if not primary_species:
                    continue
                concordant = True
                for ref_db in alt_ref_dbs:
                    alt_evidence = evidence_rows.get(str(ref_db))
                    alt_entry = (alt_taxonomies.get(sid) or {}).get(str(ref_db))
                    if not alt_evidence or str(alt_evidence.get('assignment_status') or '').upper() != 'CONFIDENT' or not alt_entry:
                        concordant = False
                        break
                    alt_species = str(parse_taxon_string(str(alt_entry[0] or '')).get('s') or '').strip()
                    if not alt_species or normalise_taxon_name(alt_species) != normalise_taxon_name(primary_species):
                        concordant = False
                        break
                if concordant:
                    confident_taxonomy_ids.add(sid)
        except Exception as exc:
            logger.warning('[NOVELTY] Could not validate cross-database taxonomy for pangenome context: %s', exc)
            confident_taxonomy_ids = set()

        for sid, entry in by_id.items():
            if sid in qids:
                query_taxonomy[sid] = entry.get('taxonomy', '')
            dataset = entry.get('dataset', '')
            role = dataset_roles.get(dataset, {})
            is_baseline_genome = dataset in baseline_genome_datasets
            tier = entry.get('baseline_tier') or role.get('baseline_tier') or ''
            if is_baseline_genome and tier not in {'priority', 'secondary'}:
                tier = _infer_baseline_tier_name(dataset)
            is_hungate_baseline = is_baseline_genome and tier == 'priority'
            is_secondary_baseline = is_baseline_genome and tier == 'secondary'
            # A candidate committed to genome sequencing counts as project-collection
            # evidence whether the genome is already in hand or only selected/pending —
            # the distinction does not matter for novelty/pangenome-gap purposes because
            # a selected isolate WILL eventually be sequenced.
            is_partner_genome = bool(entry.get('available') or entry.get('selected')) and not is_baseline_genome
            collection_metadata[sid] = {
                'taxonomy': entry.get('taxonomy', ''),
                'dataset': dataset,
                'baseline_tier': tier,
                'baseline_available': is_baseline_genome,
                'hungate_baseline_available': is_hungate_baseline,
                'secondary_baseline_available': is_secondary_baseline,
                'partner_available': is_partner_genome,
                # Retained for output-column compatibility; selected-but-not-yet-sequenced
                # isolates are now folded into 'partner_available' above rather than
                # tracked as a separate pending category.
                'selected_pending': False,
            }
            if (is_baseline_genome or is_partner_genome) and sid in confident_taxonomy_ids:
                genome_records.append((sid, entry['sequence']))
                genome_metadata[sid] = collection_metadata[sid]
    except Exception as e:
        logger.warning('[NOVELTY] Failed to load genome collection from DB: %s', e)
        genome_records = []
        genome_metadata = {}
        collection_metadata = {}

    selected_hits = {qid: [] for qid in qids}
    if genome_records:
        selected_pool = Path(outdir) / 'novelty_selected_for_wgs_pool.fasta'
        write_fasta(genome_records, str(selected_pool))
        matches_path = Path(outdir) / 'novelty_selected_for_wgs_matches.tsv'
        try:
            _run_density_search(input_fasta, str(selected_pool), matches_path, threads=threads)
            if matches_path.exists():
                with open(matches_path) as fh:
                    for line in fh:
                        parts = line.rstrip('\n').split('\t')
                        if len(parts) < 3:
                            continue
                        qid, hit_id = parts[0], parts[1]
                        if qid not in selected_hits or hit_id == qid:
                            continue
                        try:
                            pct = float(parts[2])
                        except Exception:
                            continue
                        selected_hits[qid].append((pct, hit_id))
        except Exception as e:
            logger.warning('[NOVELTY] Failed to compute selected-for-WGS neighbourhood: %s', e)

    def species_key(taxonomy):
        parsed = parse_taxon_string(str(taxonomy or ''))
        species = str(parsed.get('s') or '').strip()
        if not species or species.lower() in ('na', 'none', 'unclassified'):
            return None
        if species.startswith('s__') and not species[3:].strip():
            return None
        return species.lower()

    context = {}
    for qid in qids:
        hits = selected_hits.get(qid, [])
        ge_99 = sum(1 for pct, _ in hits if pct >= 99.0)
        ge_97 = sum(1 for pct, _ in hits if pct >= 97.0)
        ge_95 = sum(1 for pct, _ in hits if pct >= 95.0)
        nearest_pct, nearest_hit = (0.0, 'None')
        if hits:
            nearest_pct, nearest_hit = max(hits, key=lambda item: item[0])
        meta = metadata.get(qid, {})
        current_available = meta.get('selected_for_wgs')
        current_selected = meta.get('selected_for_sequencing')
        query_species = species_key(query_taxonomy.get(qid)) if qid in confident_taxonomy_ids else None
        same_species_ids = [
            sid for sid, values in collection_metadata.items()
            if sid in confident_taxonomy_ids
            and query_species and species_key(values.get('taxonomy')) == query_species
        ]
        same_species_hungate = sum(
            1 for sid in same_species_ids if collection_metadata[sid].get('hungate_baseline_available')
        )
        same_species_secondary_baseline = sum(
            1 for sid in same_species_ids if collection_metadata[sid].get('secondary_baseline_available')
        )
        same_species_baseline = same_species_hungate + same_species_secondary_baseline
        same_species_partner_available = sum(
            1 for sid in same_species_ids if collection_metadata[sid].get('partner_available')
        )
        same_species_pending = sum(
            1 for sid in same_species_ids if collection_metadata[sid].get('selected_pending')
        )
        same_species_committed = same_species_baseline + same_species_partner_available + same_species_pending
        # The authoritative pangenome gap is deferred until the tree and
        # local-clade context have been built and cross-checked in selection_sets.
        pangenome_gap = None

        # The name of the nearest OTHER submitted isolate, regardless of whether it has
        # been committed to genome sequencing. This is reported so a reviewer can always
        # see exactly which project isolate a candidate resembles, not just whether a
        # "genome collection" match exists.
        project_row = project_metrics.get(qid) or {}
        closest_project_isolate = project_row.get('nearest_hit') or 'None'
        try:
            closest_project_isolate_identity = float(project_row.get('nearest_identity'))
        except (TypeError, ValueError):
            closest_project_isolate_identity = 0.0

        # A candidate can look "not novel" purely because it resembles another
        # submitted isolate that has not itself been validated against any cultured
        # baseline or external reference. Flag that scenario explicitly so it is not
        # mistaken for genuine evidence that the organism is already known.
        reference_row = reference_metrics.get(qid) or {}
        cultured_row = cultured_metrics.get(qid) or {}
        baseline_close_hit = _pool_has_close_hit(cultured_row)
        reference_close_hit = _pool_has_close_hit(reference_row)
        project_close_hit = _pool_has_close_hit(project_row)
        project_isolate_only_match = bool(project_close_hit and not baseline_close_hit and not reference_close_hit)

        base = primary_metrics.get(qid, {})
        try:
            base_score = float(base.get('novelty_score', 0.0))
        except Exception:
            base_score = 0.0
        base_priority = str(base.get('sequencing_priority', 'NA'))
        clade_already_selected = ge_97 > 0
        if current_available:
            adjusted_priority = 'SEQUENCED'
            adjusted_score = base_score
        elif current_selected:
            adjusted_priority = 'SELECTED_PENDING'
            adjusted_score = base_score
        elif clade_already_selected:
            adjusted_priority = 'MEDIUM_RELATED_GENOME_CLADE'
            adjusted_score = max(0.0, base_score - 10.0)
        elif ge_95 > 0 and base_priority == 'HIGH':
            adjusted_priority = 'MEDIUM_SELECTED_NEARBY'
            adjusted_score = max(0.0, base_score - 10.0)
        else:
            adjusted_priority = base_priority
            adjusted_score = base_score
        context[qid] = {
            'partner_id': meta.get('partner_id', 'NA'),
            'selected_for_sequencing': 'True' if current_selected is True else ('False' if current_selected is False else 'NA'),
            'already_sequenced': 'True' if current_available is True else ('False' if current_available is False else 'NA'),
            'selected_for_wgs': 'True' if current_available is True else ('False' if current_available is False else 'NA'),
            'nearest_genome_hit': nearest_hit,
            'nearest_genome_identity': nearest_pct,
            'closest_project_isolate': closest_project_isolate,
            'closest_project_isolate_identity': closest_project_isolate_identity,
            'project_isolate_only_match': str(project_isolate_only_match),
            'genome_ge_99': ge_99,
            'genome_ge_97': ge_97,
            'genome_ge_95': ge_95,
            'related_genome_clade': str(clade_already_selected),
            'adjusted_novelty_score': adjusted_score,
            'adjusted_priority': adjusted_priority,
            'source': 'baseline_genomes_and_project_ledger' if metadata or genome_records else 'none',
            'genome_same_species_committed': same_species_committed,
            'genome_same_species_hungate': same_species_hungate,
            'genome_same_species_secondary_baseline': same_species_secondary_baseline,
            'genome_same_species_available': same_species_baseline,
            'genome_same_species_selected': same_species_partner_available,
            'genome_same_species_pending': same_species_pending,
            'pangenome_target': 'NA',
            'pangenome_gap': 'NA',
            'pangenome_eligibility_status': 'DEFERRED_TO_TREE_VALIDATED_SELECTION_SETS',
        }
    return context


def _metric_values(metrics, qid):
    row = metrics.get(qid) or _empty_pool_metrics([qid])[qid]
    return [
        f"{float(row['nearest_identity']):.2f}",
        row['nearest_hit'],
        row['novel'],
        str(row['matches_ge_99']),
        str(row['matches_ge_97']),
        str(row['matches_ge_95']),
        f"{float(row['novelty_score']):.2f}",
        row['crowding'],
        row['sequencing_priority'],
        'NA' if row.get('nearest_query_coverage') is None else f"{float(row['nearest_query_coverage']):.2f}",
        'NA' if row.get('nearest_alignment_length') is None else str(int(row['nearest_alignment_length'])),
        row['density_source'],
        row.get('density_count_method', 'NOT_AVAILABLE'),
    ]


def _pool_has_hit(row: dict) -> bool:
    if not row:
        return False
    source = str(row.get('density_source') or '')
    return source not in {'', 'none', 'NA'}


def _pool_has_close_hit(row: dict, identity_threshold: float = 97.0) -> bool:
    """Whether a pool comparison found a hit at/above the species-level identity boundary."""
    if not _pool_has_hit(row):
        return False
    try:
        identity = float(row.get('nearest_identity'))
        coverage = float(row.get('nearest_query_coverage'))
    except (TypeError, ValueError):
        return False
    return (
        identity >= float(identity_threshold)
        and coverage >= MIN_NEAREST_QUERY_COVERAGE
    )


def _passes_identity_coverage(row: dict, identity_threshold: float, coverage_threshold: float) -> bool:
    if not _pool_has_hit(row):
        return False
    try:
        identity = float(row.get('nearest_identity'))
        coverage = float(row.get('nearest_query_coverage'))
    except (TypeError, ValueError):
        return False
    return identity >= float(identity_threshold) and coverage >= float(coverage_threshold)


def _baseline_evidence_class(qid, hungate_metrics, secondary_metrics, cultured_metrics) -> str:
    hungate = hungate_metrics.get(qid, {})
    secondary = secondary_metrics.get(qid, {})
    cultured = cultured_metrics.get(qid, {})
    hungate_available = _pool_has_hit(hungate)
    secondary_available = _pool_has_hit(secondary)
    if not any(_pool_has_hit(row) for row in (hungate, secondary, cultured)):
        return 'NO CULTURED BASELINE AVAILABLE'
    if any(
        _passes_identity_coverage(row, DEFAULT_BASELINE_REDUNDANCY_IDENTITY, DEFAULT_BASELINE_REDUNDANCY_QUERY_COVERAGE)
        for row in (hungate, secondary)
    ):
        return 'CULTURED BASELINE REDUNDANT'
    hungate_covered = _passes_identity_coverage(
        hungate, DEFAULT_BASELINE_EXTENSION_MIN_IDENTITY, DEFAULT_BASELINE_EXTENSION_MIN_QUERY_COVERAGE,
    )
    secondary_covered = _passes_identity_coverage(
        secondary, DEFAULT_BASELINE_EXTENSION_MIN_IDENTITY, DEFAULT_BASELINE_EXTENSION_MIN_QUERY_COVERAGE,
    )
    if hungate_covered:
        return 'HUNGATE COVERED'
    if secondary_covered:
        return (
            'HUNGATE GAP / SECONDARY COVERED' if hungate_available
            else 'SECONDARY COVERED / HUNGATE NOT ASSESSED'
        )
    if hungate_available and secondary_available:
        return 'NOVEL TO BOTH BASELINES'
    if hungate_available:
        return 'NOVEL TO HUNGATE / SECONDARY NOT ASSESSED'
    if secondary_available:
        return 'NOVEL TO SECONDARY BASELINE / HUNGATE NOT ASSESSED'
    return 'NO CULTURED BASELINE AVAILABLE'


def _write_novelty_metrics_table(
    output: Path,
    qids,
    primary_metrics,
    hungate_metrics,
    secondary_metrics,
    cultured_metrics,
    project_metrics,
    reference_metrics,
    sequencing_context=None,
):
    metric_fields = [
        'NearestIdentity', 'NearestHit', 'Novel', 'MatchesGE99', 'MatchesGE97', 'MatchesGE95',
        'NoveltyScore', 'Crowding', 'SequencingPriority', 'NearestQueryCoverage',
        'NearestAlignmentLength', 'DensitySource', 'DensityCountMethod',
    ]
    headers = ['ID']
    headers.extend(metric_fields)
    for prefix in ('Hungate', 'SecondaryBaseline', 'CulturedProject', 'Baseline', 'Project', 'Reference'):
        headers.extend(prefix + field for field in metric_fields)
    headers.extend([
        'BaselineEvidenceClass',
        'PartnerID', 'AlreadySelectedForGenomeSequencing', 'GenomeAlreadySequenced',
        'NearestGenomeHit', 'NearestGenomeIdentity',
        'GenomeCollectionMatchesGE99', 'GenomeCollectionMatchesGE97', 'GenomeCollectionMatchesGE95',
        'RelatedGenomeCladeGE97',
        'CommittedGenomeCountSameAssessmentSpecies',
        'HungateGenomeCountSameAssessmentSpecies',
        'SecondaryBaselineGenomeCountSameAssessmentSpecies',
        'CulturedProjectGenomeCountSameAssessmentSpecies',
        'BaselineGenomeCountSameAssessmentSpecies',
        'SequencedPartnerGenomeCountSameAssessmentSpecies',
        'SelectedPendingGenomeCountSameAssessmentSpecies',
        'PangenomeTarget', 'PangenomeGap', 'PangenomeEligibilityStatus',
        'GenomeSequencingMetadataSource',
        'ClosestProjectIsolate', 'ClosestProjectIsolateIdentity', 'ProjectIsolateOnlyMatch',
        'ReferenceEvidenceInterpretation', 'ProjectEvidenceInterpretation',
        'GenomeRepresentationInterpretation',
    ])

    with open(output, 'w') as fh:
        fh.write('\t'.join(headers) + '\n')
        sequencing_context = sequencing_context or _empty_sequencing_context(qids)
        for qid in qids:
            primary = _metric_values(primary_metrics, qid)
            hungate = _metric_values(hungate_metrics, qid)
            secondary = _metric_values(secondary_metrics, qid)
            cultured = _metric_values(cultured_metrics, qid)
            baseline_alias = cultured
            project = _metric_values(project_metrics, qid)
            reference = _metric_values(reference_metrics, qid)
            evidence_class = _baseline_evidence_class(qid, hungate_metrics, secondary_metrics, cultured_metrics)
            seq_ctx = sequencing_context.get(qid) or _empty_sequencing_context([qid])[qid]
            cultured_count = int(seq_ctx.get('genome_same_species_available') or 0)
            sequencing = [
                evidence_class,
                seq_ctx['partner_id'],
                seq_ctx['selected_for_sequencing'],
                seq_ctx['already_sequenced'],
                seq_ctx['nearest_genome_hit'],
                f"{float(seq_ctx['nearest_genome_identity']):.2f}",
                str(seq_ctx['genome_ge_99']),
                str(seq_ctx['genome_ge_97']),
                str(seq_ctx['genome_ge_95']),
                seq_ctx['related_genome_clade'],
                str(seq_ctx['genome_same_species_committed']),
                str(seq_ctx.get('genome_same_species_hungate', 0)),
                str(seq_ctx.get('genome_same_species_secondary_baseline', 0)),
                str(cultured_count),
                str(cultured_count),
                str(seq_ctx['genome_same_species_selected']),
                str(seq_ctx['genome_same_species_pending']),
                str(seq_ctx['pangenome_target']),
                str(seq_ctx['pangenome_gap']),
                str(seq_ctx.get('pangenome_eligibility_status', 'DEFERRED_TO_TREE_VALIDATED_SELECTION_SETS')),
                seq_ctx['source'],
                seq_ctx.get('closest_project_isolate', 'None'),
                f"{float(seq_ctx.get('closest_project_isolate_identity') or 0.0):.2f}",
                seq_ctx.get('project_isolate_only_match', 'False'),
                'MARKER_SEQUENCE_RECOGNITION_ONLY',
                'PROJECT_MARKER_SEQUENCE_REPRESENTATION_NOT_PROOF_OF_GENOME_AVAILABILITY',
                'CULTURED_BASELINE_AND_COMPLETED_PARTNER_GENOME_EVIDENCE',
            ]
            fh.write('\t'.join([qid] + primary + hungate + secondary + cultured + baseline_alias + project + reference + sequencing) + '\n')
    return str(output)


def _run_supplementary_reference_checks(
    input_fasta: str, references, outdir: str, threads: Optional[int] = None,
) -> Optional[str]:
    """Assess marker recognition in supplementary 16S references without genome scoring."""
    import csv
    rows = []
    seen = set()
    for item in references or []:
        if not item:
            continue
        ref_path = item[0] if isinstance(item, (tuple, list)) else item
        name = item[2] if isinstance(item, (tuple, list)) and len(item) > 2 else Path(str(ref_path)).stem
        if not ref_path or not Path(ref_path).exists():
            logger.warning('[NOVELTY] Supplementary reference %s is unavailable; skipping', name)
            continue
        resolved = str(Path(ref_path).resolve())
        if resolved in seen:
            continue
        seen.add(resolved)
        slug = ''.join(ch if ch.isalnum() else '_' for ch in str(name)).strip('_')[:40] or 'reference'
        pool = _ensure_uncompressed(str(ref_path), outdir)
        matches_path = Path(outdir) / f'novelty_supplementary_{slug}_nearest_matches.tsv'
        nearest = _run_nearest_search(input_fasta, pool, matches_path, threads=threads)
        for qid, _sequence in read_fasta(input_fasta):
            hit = nearest.get(qid)
            rows.append({
                'sequence_id': qid, 'reference_database': str(name),
                'nearest_hit': hit[1] if hit else 'NA',
                'identity': f'{hit[0]:.3f}' if hit else '0.000',
                'query_coverage': f'{hit[2]:.2f}' if hit and hit[2] is not None else 'NA',
                'alignment_length': hit[3] if hit and hit[3] is not None else 'NA',
                'coverage_tier': hit[6] if hit else 'NONE',
                'interpretation': 'MARKER_SEQUENCE_RECOGNITION_ONLY',
            })
    if not rows:
        return None
    output = Path(outdir) / 'novelty_supplementary_reference_summary.tsv'
    fields = ['sequence_id','reference_database','nearest_hit','identity','query_coverage',
              'alignment_length','coverage_tier','interpretation']
    with output.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    return str(output)


def build_reference_novelty_metrics(
    input_fasta: str,
    ref_fasta: str,
    outdir: str,
    threads: Optional[int] = None,
    db=None,
    run_dataset: Optional[str] = None,
    target_fasta: Optional[str] = None,
    baseline_datasets=None,
    hungate_datasets=None,
    secondary_baseline_datasets=None,
    supplementary_references=None,
    pangenome_target: int = DEFAULT_PANGENOME_TARGET,
):
    """Write per-sequence novelty metrics comparing against project and baseline sequences.

        Project density is computed against the rolling partner-candidate collection,
        including selected partner isolates and the current batch but excluding each
        query's self-hit. Nearest-hit searches require at least 95% query coverage and
        rank eligible matches by percent identity. Cultured-project baselines and the
        external GTDB reference remain separate pools.

    If ``target_fasta`` is provided it replaces the DB candidate pool and the
    ``ref_fasta`` fallback. Selected partner isolates are added from the DB so
    committed project isolates remain visible in project-neighbour comparisons.

    Falls back to ``ref_fasta`` only when no DB sequences and no
    ``target_fasta`` are available.

    Parameters
    ----------
    input_fasta  : Query sequences (dereplicated, QC-passed).
    ref_fasta    : External reference FASTA — used as a last resort fallback only.
    outdir       : Output directory.
    threads      : VSEARCH thread count.
    db           : Database instance (used to pull previously submitted seqs).
    run_dataset  : Current run dataset name (excluded from density DB).
    target_fasta : Explicit project comparison FASTA; selected partner isolates are added when available.
    baseline_datasets : DB dataset label(s) to treat as the cultured/reference
                        baseline pool, e.g. Hungate. Can be a string,
                        comma-separated string, or iterable.

    If the DB contains ``sequencing_metadata`` rows, the output also reports
    whether each query or any close 16S neighbour (>=97% identity) has already
    been selected for WGS/full-genome sequencing.
    """
    output = Path(outdir) / 'novelty_metrics.tsv'
    # The metrics query FASTA may contain the entire rolling candidate
    # collection, not only the newest submission represented in novelty.tsv.
    qids = [
        h for h, _ in read_fasta(input_fasta)
        if not tree_pipeline.is_ref_anchor(h)
    ]
    if not qids:
        _write_novelty_metrics_table(
            output, [], {}, {}, {}, {}, {}, {}, sequencing_context={},
        )
        return str(output)

    hungate_names, secondary_names, baseline_names = _discover_tiered_baseline_datasets(
        db,
        baseline_datasets=baseline_datasets,
        hungate_datasets=hungate_datasets,
        secondary_baseline_datasets=secondary_baseline_datasets,
    )

    all_pool_path: Optional[str] = None
    all_source = 'none'
    all_datasets: list[str] = []
    all_count = 0
    if target_fasta and Path(target_fasta).exists():
        all_pool_path, all_source, all_datasets, all_count = _combine_project_target_with_selected_partners(
            target_fasta, db, outdir,
        )
        logger.info(
            '[NOVELTY] Using explicit target FASTA%s for project novelty metrics: %s',
            ' plus selected partner isolates' if all_source != 'target_fasta' else '',
            target_fasta,
        )
    elif db is not None:
        try:
            candidate_datasets = db.get_dataset_names_by_role('candidate')
        except Exception:
            candidate_datasets = []
        selected_partner_ids = _selected_partner_sequence_ids(db) if candidate_datasets else []
        all_pool_path, all_source, all_datasets, all_count = _build_db_pool_fasta(
            db=db,
            outdir=outdir,
            name='project_collection',
            include_datasets=candidate_datasets or None,
            include_ids=selected_partner_ids,
            exclude_datasets=None if candidate_datasets else baseline_names,
        )
        all_source = 'project_collection' if all_pool_path else all_source
        if all_pool_path:
            logger.info(
                '[NOVELTY] Project pool: %d sequences from %d dataset(s): %s',
                all_count, len(all_datasets), all_datasets,
            )

    if all_pool_path is None and ref_fasta:
        all_pool_path = _ensure_uncompressed(ref_fasta, outdir)
        all_source = 'reference_fasta'
        logger.info(
            '[NOVELTY] WARNING: No partner-candidate project pool available. '
            'Falling back to external reference (%s).',
            ref_fasta,
        )

    if all_pool_path:
        all_nearest = _run_nearest_search(
            input_fasta,
            all_pool_path,
            Path(outdir) / 'novelty_project_nearest_matches.tsv',
            threads=threads,
            exclude_self=True,
        )
        all_density, all_density_80 = _run_density_search_coverage_pair(
            input_fasta, all_pool_path,
            Path(outdir) / 'novelty_project_density_matches.tsv',
            threads=threads, exclude_self=True,
        )
        all_metrics = _metrics_from_nearest_and_density(
            qids,
            all_nearest,
            all_density,
            source=all_source,
            id_threshold=0.97,
        )
    else:
        logger.warning('[NOVELTY] No project comparison pool available; project metrics will be empty')
        all_metrics = _empty_pool_metrics(qids, source='none')

    reference_pool_path: Optional[str] = None
    reference_source = 'none'
    if ref_fasta:
        try:
            reference_pool_path = _ensure_uncompressed(ref_fasta, outdir)
            reference_source = 'reference_fasta'
        except Exception as e:
            raise RuntimeError(f'could not prepare external-reference novelty pool {ref_fasta}: {e}') from e

    if reference_pool_path:
        reference_nearest = _run_nearest_search(
            input_fasta,
            reference_pool_path,
            Path(outdir) / 'novelty_reference_nearest_matches.tsv',
            threads=threads,
        )
        reference_density, reference_density_80 = _run_density_search_coverage_pair(
            input_fasta, reference_pool_path,
            Path(outdir) / 'novelty_reference_density_matches.tsv',
            threads=threads,
        )
        reference_metrics = _metrics_from_nearest_and_density(
            qids,
            reference_nearest,
            reference_density,
            source=reference_source,
            id_threshold=0.97,
        )
    else:
        reference_metrics = _empty_pool_metrics(qids, source='none')
        reference_density = {}
        reference_density_80 = {}

    if all_pool_path or reference_pool_path:
        import csv
        marker_comparison = Path(outdir) / 'novelty_marker_density_coverage_comparison.tsv'
        with open(marker_comparison, 'w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(['sequence_id','collection','query_coverage','matches_ge_95','matches_ge_97','matches_ge_99'])
            for collection, density_sets in (
                ('PROJECT', {95: all_density if all_pool_path else {}, 80: all_density_80 if all_pool_path else {}}),
                ('GTDB_REFERENCE', {95: reference_density, 80: reference_density_80}),
            ):
                for qid in qids:
                    for coverage, counts in density_sets.items():
                        values = counts.get(qid, {})
                        writer.writerow([qid, collection, coverage, values.get('ge_95',0), values.get('ge_97',0), values.get('ge_99',0)])

    # Hungate and secondary cultured-genome baselines are searched together once.
    # Target provenance is retained so group-specific nearest hits and density
    # counts can still be reported independently.
    combined_pool_path, baseline_provenance, baseline_found, combined_count = (
        _build_combined_baseline_pool_with_provenance(
            db=db, outdir=outdir, hungate_datasets=hungate_names,
            secondary_datasets=secondary_names,
        )
    )
    hungate_pool_path = combined_pool_path if any(v == 'HUNGATE' for v in baseline_provenance.values()) else None
    secondary_pool_path = combined_pool_path if any(v == 'SECONDARY' for v in baseline_provenance.values()) else None
    cultured_pool_path = combined_pool_path
    hungate_source = 'hungate_baseline:' + ','.join(hungate_names) if hungate_pool_path else 'none'
    secondary_source = 'secondary_baseline:' + ','.join(secondary_names) if secondary_pool_path else 'none'
    cultured_source = 'cultured_project:' + ','.join(baseline_found) if combined_pool_path else 'none'
    hungate_nearest = {}
    secondary_nearest = {}
    if combined_pool_path:
        nearest_by_group = _run_nearest_search_by_provenance(
            input_fasta, combined_pool_path,
            Path(outdir) / 'novelty_baseline_nearest_matches.tsv',
            baseline_provenance, threads=threads,
        )
        hungate_nearest = nearest_by_group['HUNGATE']
        secondary_nearest = nearest_by_group['SECONDARY']
        # The combined nearest-hit output is the best hit across either baseline.
        combined_nearest = {}
        for qid in qids:
            candidates = [h[qid] for h in (hungate_nearest, secondary_nearest) if qid in h]
            if candidates:
                combined_nearest[qid] = max(candidates, key=lambda x: (x[0], x[2] or 0, x[3] or 0, x[1]))
        # Write the combined result using the same stable schema as the group files.
        import csv
        combined_path = Path(outdir) / 'novelty_baseline_nearest_matches.tsv'
        match_fields = ['query','target','id','alnlen','mism','gaps','qlo','qhi','tlo','thi','ql','tl','qcov','tcov',
                        'search_threshold','search_round','coverage_tier','match_status']
        with open(combined_path, 'w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n'); writer.writerow(match_fields)
            for qid in qids:
                hit = combined_nearest.get(qid)
                if hit is None: writer.writerow([qid,'NA','NA'] + ['NA'] * 11 + ['NA','NA','NONE','NO_MATCH'])
                else: writer.writerow(hit[7] + [hit[4],hit[5],hit[6], 'MATCH' if hit[6]=='FULL' else 'PARTIAL_COVERAGE'])

        density_by_coverage = _run_combined_baseline_density(
            input_fasta, combined_pool_path, outdir, baseline_provenance, threads=threads,
        )
        # Keep the existing 95%-coverage density outputs as the primary metrics;
        # the 80%-coverage outputs are an explicit sensitivity comparison.
        density95 = density_by_coverage[95]
        def _merge_density(groups):
            merged = {}
            for group_counts in groups:
                for qid, values in group_counts.items():
                    bucket = merged.setdefault(qid, {'ge_99':0,'ge_97':0,'ge_95':0})
                    for key in bucket: bucket[key] += values[key]
            return merged
        hungate_density = density95['HUNGATE']
        secondary_density = density95['SECONDARY']
        combined_density = _merge_density([hungate_density, secondary_density])
        # Save a compact comparison table so the 80% and 95% coverage effects are visible.
        comparison_path = Path(outdir) / 'novelty_baseline_density_coverage_comparison.tsv'
        with open(comparison_path, 'w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(['sequence_id','collection','query_coverage','matches_ge_95','matches_ge_97','matches_ge_99'])
            for collection in ('HUNGATE','SECONDARY'):
                for qid in qids:
                    for coverage in (95,80):
                        values = density_by_coverage[coverage][collection].get(qid, {})
                        writer.writerow([qid, collection, coverage, values.get('ge_95',0), values.get('ge_97',0), values.get('ge_99',0)])
        hungate_metrics = _metrics_from_nearest_and_density(qids, hungate_nearest, hungate_density, source=hungate_source, id_threshold=0.97)
        secondary_metrics = _metrics_from_nearest_and_density(qids, secondary_nearest, secondary_density, source=secondary_source, id_threshold=0.97)
        cultured_metrics = _metrics_from_nearest_and_density(qids, combined_nearest, combined_density, source=cultured_source, id_threshold=0.97)
    else:
        hungate_metrics = _empty_pool_metrics(qids, source='none')
        secondary_metrics = _empty_pool_metrics(qids, source='none')
        cultured_metrics = _empty_pool_metrics(qids, source='none')
        if baseline_names and db is not None:
            raise RuntimeError(f'no cultured-baseline sequences were found for registered dataset(s): {baseline_names}')
    cultured_found = baseline_found
    cultured_count = combined_count
    primary_metrics = cultured_metrics if cultured_pool_path else all_metrics
    sequencing_context = _compute_sequencing_context(
        input_fasta,
        outdir,
        db,
        qids,
        primary_metrics,
        threads=threads,
        pangenome_target=pangenome_target,
        project_metrics=all_metrics,
        reference_metrics=reference_metrics,
        cultured_metrics=cultured_metrics,
    )
    _write_novelty_metrics_table(
        output,
        qids,
        primary_metrics,
        hungate_metrics,
        secondary_metrics,
        cultured_metrics,
        all_metrics,
        reference_metrics,
        sequencing_context,
    )
    genome_derived_context = _write_genome_derived_context(
        input_fasta, outdir, db, qids, threads=threads,
    )
    supplementary_summary = _run_supplementary_reference_checks(
        input_fasta, supplementary_references, outdir, threads=threads,
    )
    logger.info(
        '[NOVELTY] Wrote novelty metrics to %s (primary=%s; hungate=%s; secondary=%s; cultured=%s; project=%s; reference=%s; genome-derived-context=%s; supplementary=%s)',
        output,
        'cultured_project' if cultured_pool_path else all_source,
        hungate_source,
        secondary_source,
        cultured_source,
        all_source,
        reference_source,
        genome_derived_context or 'none',
        supplementary_summary or 'none',
    )
    return str(output)