"""Local browser UI for inspecting Sanger chromatograms and reviewing markers.

Studio is a thin interface over the existing Paper Trail/Merge Meeting engine:
it does not implement a second base-caller, trimmer, or assembler. It binds to
localhost, can load an Assistant cohort-preparation directory (pairing partner
*_onboarding maps with *_paper_trail taxonomy outputs) or standalone sample maps,
and writes inspections/review decisions to a user-selected directory.
"""

from __future__ import annotations

import csv
import copy
import json
import mimetypes
import os
import re
import threading
import urllib.parse
import uuid
import webbrowser
import zipfile
from datetime import datetime, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

from branchmanager.pipeline import paper_trail
from branchmanager.run_manifest import RunManifest


def _raw_ab1_index_for_trimmed_position(
    trimmed_position: int,
    direction: str,
    trim_start: int,
    trim_end: int,
    raw_length: int,
) -> int | None:
    """Translate a 1-based oriented trimmed-read position to a 0-based AB1 index.

    Paper Trail stores trim bounds in original AB1 coordinates. Forward reads
    advance from trim_start; reverse reads advance backwards from trim_end.
    """
    try:
        pos, start, end, length = map(int, (trimmed_position, trim_start, trim_end, raw_length))
    except (TypeError, ValueError):
        return None
    if pos < 1 or length < 1 or end <= start or pos > (end - start):
        return None
    raw_index = end - pos if str(direction or '').strip().lower().startswith('rev') else start + pos - 1
    return raw_index if 0 <= raw_index < length else None



_MAP_ID_KEYS = ('sequence_id', 'seq_id', 'id', 'isolate_id', 'sample_id', 'sample')
_MAP_FILE_KEYS = ('read_file', 'file', 'filename', 'path', 'file_path', 'ab1_file')
_MAP_PRIMER_KEYS = ('primer', 'primer_name')
_MAP_DIRECTION_KEYS = ('direction', 'orientation', 'read_direction')
_MAP_MODE_KEYS = ('processing_mode', 'mode', 'read_mode')
_READ_EXTENSIONS = {'.ab1', '.abi'}
_MAP_EXTENSIONS = {'.tsv', '.csv', '.txt'}
_REVIEW_DECISIONS = {'approve', 'approved', 'accept', 'accepted', 'pass', 'yes', 'reject', 'rejected', 'fail', 'failed', 'no'}

# Conservative starting points only: amplicon length varies substantially by assay.
# Users can edit both thresholds in the QC panel after choosing a profile.
_SEQUENCE_TYPE_PROFILES = {
    '16S': {'min_read_length': 300, 'min_final_length': 800, 'label': '16S rRNA'},
    'ITS': {'min_read_length': 100, 'min_final_length': 150, 'label': 'ITS'},
    '18S': {'min_read_length': 200, 'min_final_length': 300, 'label': '18S rRNA'},
    '28S_LSU': {'min_read_length': 200, 'min_final_length': 300, 'label': '28S / LSU rRNA'},
    'CUSTOM': {'min_read_length': None, 'min_final_length': None, 'label': 'Other / custom target'},
}

_DIALOG_LOCK = threading.Lock()

def _native_dialog(kind: str, *, initialdir: str | Path | None = None):
    """Open a native local file/folder dialog from the localhost UI.

    The dialog returns paths only; the HTTP API never exposes arbitrary directory
    browsing or file contents.
    """
    try:
        import tkinter as tk
        from tkinter import filedialog
    except Exception as exc:  # pragma: no cover - depends on desktop installation
        raise RuntimeError(
            'Native file dialogs are unavailable. Install/enable Python tkinter and run Sanger Studio in a desktop session.'
        ) from exc

    with _DIALOG_LOCK:
        window = None
        try:
            window = tk.Tk()
            window.withdraw()
            try:
                window.attributes('-topmost', True)
            except Exception:
                pass  # -topmost not supported on all platforms
            window.update_idletasks()
            start = str(Path(initialdir).expanduser().resolve()) if initialdir else str(Path.cwd())
            if kind == 'map':
                result = filedialog.askopenfilename(
                    parent=window, title='Select Sanger sample map', initialdir=start,
                    filetypes=[('Sample maps', '*.tsv *.csv *.txt'), ('All files', '*.*')],
                )
            elif kind == 'review':
                result = filedialog.askopenfilename(
                    parent=window, title='Attach cohort marker review template', initialdir=start,
                    filetypes=[('TSV / CSV review files', '*.tsv *.csv *.txt'), ('All files', '*.*')],
                )
            elif kind == 'reads':
                result = filedialog.askopenfilenames(
                    parent=window, title='Select AB1 / ABI chromatograms', initialdir=start,
                    filetypes=[('Sanger chromatograms', '*.ab1 *.abi'), ('All files', '*.*')],
                )
            elif kind == 'taxonomy_ref':
                result = filedialog.askopenfilename(
                    parent=window, title='Select taxonomy reference FASTA', initialdir=start,
                    filetypes=[
                        ('FASTA reference databases', '*.fasta *.fa *.fna *.fas *.fasta.gz *.fa.gz *.fna.gz *.fas.gz'),
                        ('All files', '*.*'),
                    ],
                )
            elif kind == 'taxonomy_taxa':
                result = filedialog.askopenfilename(
                    parent=window, title='Select reference taxonomy mapping (optional)', initialdir=start,
                    filetypes=[('Taxonomy tables', '*.tsv *.csv *.txt *.tsv.gz *.csv.gz'), ('All files', '*.*')],
                )
            elif kind in {'input_folder', 'onboarding'}:
                result = filedialog.askdirectory(
                    parent=window,
                    title='Select Sanger input folder (Assistant cohort, sample map, or raw AB1/ABI files)',
                    initialdir=start,
                    mustexist=True,
                )
            elif kind == 'output':
                result = filedialog.askdirectory(
                    parent=window, title='Choose Sanger Studio output directory', initialdir=start,
                    mustexist=True,
                )
            else:
                raise ValueError(f'Unknown native dialog type: {kind}')
            return result
        except tk.TclError as exc:
            raise RuntimeError(
                'Could not open the native file dialog. Check that Sanger Studio is running in a graphical desktop session.'
            ) from exc
        finally:
            if window is not None:
                try:
                    window.destroy()
                except Exception:
                    pass  # window may already be destroyed


def _table_delimiter(path: str | Path, sample: str = '') -> str:
    """Choose a delimiter without allowing cell contents to override a TSV/CSV suffix.

    Paper Trail tables can contain comma-rich JSON fields (for example
    FallbackReadOptions). csv.Sniffer can misclassify those genuine TSV files
    as comma-separated and collapse the entire assembly row into one column.
    File extensions are authoritative for .tsv/.csv; sniff only ambiguous files.
    """
    name = str(path).lower()
    if name.endswith(('.tsv', '.tsv.gz', '.txt', '.txt.gz')):
        return '\t'
    if name.endswith(('.csv', '.csv.gz')):
        return ','
    try:
        return csv.Sniffer().sniff(sample, delimiters=',\t').delimiter
    except csv.Error:
        return '\t'


def _read_table_with_fields(path: str | Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read a TSV/CSV while preserving the original column order."""
    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise ValueError(f'File does not exist: {source}')
    with open(source, newline='', encoding='utf-8-sig') as handle:
        sample = handle.read(8192)
        handle.seek(0)
        delimiter = _table_delimiter(source, sample)
        reader = csv.DictReader(handle, delimiter=delimiter)
        fields = list(reader.fieldnames or [])
        if not fields:
            raise ValueError(f'File has no header row: {source}')
        rows = [dict(row) for row in reader]
    return fields, rows


def _parse_custom_primer_sequences(value: object) -> dict[str, str]:
    """Parse optional two-column primer-name / IUPAC-sequence text."""
    text = str(value or '').strip()
    if not text:
        return {}
    allowed = set(paper_trail.IUPAC_BASES)
    parsed: dict[str, str] = {}
    for line_number, line in enumerate(text.splitlines(), start=1):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        fields = re.split(r'[\t, ]+', line, maxsplit=1)
        if len(fields) != 2:
            raise ValueError(
                f'Custom primer line {line_number} must contain a primer name and sequence separated by a tab, comma or space.'
            )
        name, sequence = fields[0].strip(), re.sub(r'\s+', '', fields[1]).upper()
        if not name or not sequence or set(sequence) - allowed:
            raise ValueError(
                f'Custom primer line {line_number} has an empty name/sequence or non-IUPAC nucleotide characters.'
            )
        key = name.upper()
        if key in parsed:
            raise ValueError(f'Custom primer {name!r} is defined more than once.')
        parsed[key] = sequence
    return parsed


def _normalise_key(value: object) -> str:
    return re.sub(r'[^a-z0-9]+', '_', str(value or '').strip().lstrip('\ufeff').lower()).strip('_')


def _is_within(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
        return True
    except ValueError:
        return False


def _read_tsv(path: str | Path) -> list[dict[str, str]]:
    with open(path, newline='', encoding='utf-8-sig') as handle:
        sample = handle.read(4096)
        handle.seek(0)
        delimiter = _table_delimiter(path, sample)
        return [dict(row) for row in csv.DictReader(handle, delimiter=delimiter)]


def _write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    from branchmanager.utils.io import write_tsv
    write_tsv(path, fieldnames, rows)


def _read_id_from_path(value: object) -> str:
    """Convert a normalised-read-map file path to the Paper Trail read ID.

    Paper Trail's read ID is the filename stem (for AB1/ABI inputs), which is
    also the ID written in read_taxonomy_screen/taxonomy.tsv.
    """
    raw = str(value or '').strip()
    name = Path(raw).name
    lowered = name.lower()
    for suffix in ('.ab1.gz', '.abi.gz', '.ab1', '.abi'):
        if lowered.endswith(suffix):
            return name[:-len(suffix)]
    return Path(name).stem


def _read_normalised_read_map(path: str | Path) -> tuple[dict[str, dict[str, str]], dict[str, list[str]]]:
    """Load the onboarding read map without requiring the raw AB1 files to exist.

    The map is used only to associate isolate IDs with read IDs; raw reads are
    independently loaded through the active Studio sample map.
    """
    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise ValueError(f'Normalised read map does not exist: {source}')
    by_read: dict[str, dict[str, str]] = {}
    by_sequence: dict[str, list[str]] = {}
    for line, row in enumerate(_read_tsv(source), start=2):
        sequence_id = _row_get(row, ('sequence_id', 'seq_id', 'isolate_id', 'sample_id'))
        read_file = _row_get(row, ('file', 'read_file', 'filename', 'path', 'file_path'))
        if not sequence_id or not read_file:
            if any(str(value or '').strip() for value in row.values()):
                raise ValueError(
                    f'Normalised read-map row {line} requires sequence_id and file/read_file columns'
                )
            continue
        read_id = _read_id_from_path(read_file)
        if not read_id:
            continue
        read_path = Path(read_file).expanduser()
        if not read_path.is_absolute():
            read_path = source.parent / read_path
        mode = _row_get(row, _MAP_MODE_KEYS).lower()
        if mode not in {'assemble', 'best_read'}:
            mode = ''
        entry = {
            'sequence_id': sequence_id,
            'read_id': read_id,
            'file': str(read_path.resolve()),
            'primer': _row_get(row, ('primer', 'primer_name')),
            'direction': _row_get(row, ('direction', 'orientation', 'read_direction')),
            'processing_mode': mode,
            'partner_id': _row_get(row, ('partner_id', 'partner')),
            'dataset': _row_get(row, ('dataset', 'cohort_dataset')),
        }
        previous = by_read.get(read_id)
        if previous and previous['sequence_id'] != sequence_id:
            raise ValueError(
                f'Read ID {read_id!r} maps to both {previous["sequence_id"]!r} and '
                f'{sequence_id!r} in {source}. Select a single onboarding cohort/partner folder.'
            )
        by_read[read_id] = entry
        by_sequence.setdefault(sequence_id, [])
        if read_id not in by_sequence[sequence_id]:
            by_sequence[sequence_id].append(read_id)

    # processing_mode is optional in Assistant's normalised_read_map.tsv.
    # Still populate it on every row because the Studio inspection handler
    # persists primer/direction/mode corrections by indexing all three keys.
    # When the source map has no mode, retain the same default shown in the UI:
    # assemble multiple reads, otherwise select the best single read.
    for sequence_id, read_ids in by_sequence.items():
        explicit_modes = {
            by_read[read_id]['processing_mode']
            for read_id in read_ids
            if by_read[read_id].get('processing_mode')
        }
        default_mode = (
            next(iter(explicit_modes)) if len(explicit_modes) == 1
            else ('assemble' if len(read_ids) > 1 else 'best_read')
        )
        for read_id in read_ids:
            if not by_read[read_id].get('processing_mode'):
                by_read[read_id]['processing_mode'] = default_mode

    if not by_read:
        raise ValueError(f'Normalised read map contains no usable isolate/read records: {source}')
    return by_read, by_sequence


def _find_onboarding_taxonomy_files(selected: Path) -> tuple[list[Path], list[Path]]:
    """Find per-read taxonomy summaries and their normalised read maps.

    No read_taxonomy_concordance.tsv is required. The isolate/read relationship
    comes from normalised_read_map.tsv; taxon assignments and hit evidence come
    from taxonomy.tsv and read_taxonomy_hits.tsv.
    """
    root = selected.expanduser().resolve()
    if not root.is_dir():
        raise ValueError(f'Cohort folder selection is not a directory: {root}')

    # The preferred selection is the single 00_cohort_preparation directory
    # containing matching *_onboarding and *_paper_trail siblings. Be forgiving
    # when a user selects one partner's child folder: resolve to its parent so
    # both sides of the pair are discovered automatically.
    if root.name.endswith(('_onboarding', '_paper_trail')):
        root = root.parent
    elif root.name == 'read_taxonomy_screen':
        root = root.parent.parent

    # Use the bounded directory walker so selecting a project-level folder
    # does not recursively inspect unrelated hidden/output directories.
    taxonomy_paths = _discover_taxonomy_tables(root)
    # Also support a Paper Trail screen directory whose taxonomy.tsv is
    # directly beneath a selected directory with a different name.
    if not taxonomy_paths and (root / 'taxonomy.tsv').is_file():
        taxonomy_paths = [root / 'taxonomy.tsv']
    if not taxonomy_paths:
        raise ValueError(
            f'No read_taxonomy_screen/taxonomy.tsv found beneath {root}. '
            'Select the relevant Assistant onboarding/Paper Trail output folder.'
        )

    map_paths = _discover_normalised_maps(root)
    # In some onboarding layouts the normalised map is one directory above the
    # selected Paper Trail folder. Search a short, deterministic ancestor chain.
    if not map_paths:
        for ancestor in list(root.parents)[:6]:
            candidate = ancestor / 'normalised_read_map.tsv'
            if candidate.is_file():
                map_paths = [candidate]
                break
    if not map_paths:
        raise ValueError(
            f'No normalised_read_map.tsv found beneath {root} or in its parent folders. '
            'The read map is needed to link each taxonomic read assignment to its isolate.'
        )

    return taxonomy_paths, map_paths


def _normalise_input_directory(selected: str | Path) -> Path:
    """Resolve common selected child folders to their containing input directory."""
    root = Path(selected).expanduser().resolve()
    if not root.is_dir():
        raise ValueError(f'Input selection is not a directory: {root}')
    if root.name.endswith(('_onboarding', '_paper_trail')):
        return root.parent
    if root.name == 'read_taxonomy_screen':
        return root.parent.parent
    return root


def _iter_input_files(root: Path, suffixes: set[str], *, limit: int = 50000):
    """Walk an input tree while avoiding generated output and development folders."""
    ignored = {
        '.git', '.svn', '.hg', '__pycache__', '.pytest_cache', '.mypy_cache',
        '.venv', 'venv', 'node_modules', 'sanger_studio_runs', 'inspections',
    }
    found = 0
    visited_directories = 0
    for current, directories, filenames in os.walk(root):
        visited_directories += 1
        if visited_directories > 10000:
            return
        relative_depth = len(Path(current).relative_to(root).parts)
        directories[:] = sorted(
            name for name in directories
            if name not in ignored and not name.startswith('.')
        )
        if relative_depth >= 10:
            directories[:] = []
        for filename in sorted(filenames):
            path = Path(current) / filename
            lower = filename.lower()
            if any(lower.endswith(suffix) for suffix in suffixes):
                yield path.resolve()
                found += 1
                if found >= limit:
                    return


def _discover_normalised_maps(root: Path) -> list[Path]:
    return sorted({
        path for path in _iter_input_files(root, {'.tsv'})
        if path.name.lower() == 'normalised_read_map.tsv'
    })


def _discover_taxonomy_tables(root: Path) -> list[Path]:
    return sorted({
        path for path in _iter_input_files(root, {'.tsv'})
        if path.name.lower() == 'taxonomy.tsv' and path.parent.name == 'read_taxonomy_screen'
    })


def _discover_standalone_maps(root: Path) -> list[Path]:
    """Find TSV/CSV/TXT files that actually have sample-map columns.

    Merely having a .tsv extension is not enough: report tables, taxonomy files,
    and QC manifests are ignored unless they expose both an isolate ID and a
    physical read-file column.
    """
    candidates = []
    for path in _iter_input_files(root, {'.tsv', '.csv', '.txt'}):
        if path.name.lower() == 'normalised_read_map.tsv':
            continue
        try:
            with open(path, newline='', encoding='utf-8-sig') as handle:
                sample = handle.read(8192)
                handle.seek(0)
                delimiter = _table_delimiter(path, sample)
                reader = csv.DictReader(handle, delimiter=delimiter)
                keys = {_normalise_key(field) for field in (reader.fieldnames or [])}
                if not (keys & set(_MAP_ID_KEYS)) or not (keys & set(_MAP_FILE_KEYS)):
                    continue
                # Only inspect a few initial records while discovering files;
                # the complete map is parsed once a candidate is selected.
                looks_like_map = False
                for index, row in enumerate(reader):
                    if index >= 20:
                        break
                    if (
                        _row_get(row, _MAP_ID_KEYS)
                        and Path(_row_get(row, _MAP_FILE_KEYS)).name.lower().endswith(('.ab1', '.abi'))
                    ):
                        looks_like_map = True
                        break
        except (OSError, ValueError, csv.Error, UnicodeError):
            continue
        if looks_like_map:
            candidates.append(path)
    return sorted(set(candidates))


def _combine_normalised_read_maps(paths: list[Path]) -> tuple[list[dict[str, str]], list[str]]:
    """Combine partner normalised maps without silently accepting ID conflicts."""
    by_read: dict[str, dict[str, str]] = {}
    source_paths = []
    for path in sorted(paths):
        current, _by_sequence = _read_normalised_read_map(path)
        source_paths.append(str(path))
        for read_id, row in current.items():
            previous = by_read.get(read_id)
            if previous and previous != row:
                raise ValueError(
                    f'Normalised read ID {read_id!r} occurs in multiple input maps with conflicting metadata. '
                    f'Check {path} and {previous.get("file", "the earlier map")}. '
                    'Select the specific partner folder or correct the duplicate read ID.'
                )
            by_read[read_id] = row
    if not by_read:
        raise ValueError('No usable records were found in the normalised read maps.')
    return list(by_read.values()), source_paths


def _combine_standalone_sample_maps(paths: list[Path]) -> tuple[list[dict[str, str]], list[str]]:
    """Combine recognised standalone maps, deduplicating identical physical reads."""
    by_file: dict[str, dict[str, str]] = {}
    source_paths = []
    for path in sorted(paths):
        rows = _read_map(path)
        source_paths.append(str(path))
        for row in rows:
            key = str(Path(row['file']).resolve())
            previous = by_file.get(key)
            if previous and previous != row:
                raise ValueError(
                    f'Read file {key} is described differently in multiple sample maps. '
                    'Select one map or reconcile the duplicate rows before loading the folder.'
                )
            by_file[key] = row
    combined = list(by_file.values())
    if not combined:
        raise ValueError('The discovered standalone sample maps contain no usable read records.')
    return combined, source_paths


def _load_onboarding_taxonomy(selected: str | Path) -> dict[str, object]:
    """Load onboarding taxonomy outputs once into read- and isolate-indexed maps."""
    root = Path(selected).expanduser().resolve()
    taxonomy_paths, map_paths = _find_onboarding_taxonomy_files(root)

    read_map_by_id: dict[str, dict[str, str]] = {}
    read_ids_by_sequence: dict[str, list[str]] = {}
    used_map_paths: set[Path] = set()
    taxonomy_by_read: dict[str, dict[str, str]] = {}
    hits_by_read: dict[str, list[dict[str, str]]] = {}
    taxonomy_fields: list[str] = []
    hit_fields: list[str] = []
    policy_rows: list[dict[str, str]] = []
    policy_fields: list[str] = []
    loaded_files: list[str] = []
    trimmed_fasta_paths: list[str] = []
    trimmed_sequences_by_read: dict[str, str] = {}

    def map_for_screen(taxonomy_path: Path) -> Path | None:
        screen_dir = taxonomy_path.parent
        paper_trail_dir = screen_dir.parent if screen_dir.name == 'read_taxonomy_screen' else screen_dir

        # Assistant cohort PREPARE writes partner-specific sibling folders:
        #   01_<dataset>_onboarding/normalised_read_map.tsv
        #   01_<dataset>_paper_trail/read_taxonomy_screen/taxonomy.tsv
        # Pair these by the shared dataset prefix before considering a global
        # fallback. This is essential when several partners are selected.
        paper_name = paper_trail_dir.name
        if paper_name.endswith('_paper_trail'):
            onboarding_name = paper_name[:-len('_paper_trail')] + '_onboarding'
            sibling_map = paper_trail_dir.parent / onboarding_name / 'normalised_read_map.tsv'
            if sibling_map.is_file():
                return sibling_map.resolve()
        elif paper_name == '02_paper_trail_merge_meeting':
            sibling_map = paper_trail_dir.parent / '01_onboarding' / 'normalised_read_map.tsv'
            if sibling_map.is_file():
                return sibling_map.resolve()

        # Also support other layouts with exactly one onboarding sibling.
        sibling_maps = sorted(
            candidate / 'normalised_read_map.tsv'
            for candidate in paper_trail_dir.parent.iterdir()
            if candidate.is_dir() and candidate.name.endswith('_onboarding')
            and (candidate / 'normalised_read_map.tsv').is_file()
        ) if paper_trail_dir.parent.is_dir() else []
        if len(sibling_maps) == 1:
            return sibling_maps[0].resolve()

        for ancestor in (paper_trail_dir, *paper_trail_dir.parents):
            candidate = ancestor / 'normalised_read_map.tsv'
            if candidate.is_file():
                return candidate.resolve()
        return map_paths[0].resolve() if len(map_paths) == 1 else None

    # Associate each screen result directory with its nearest normalised map.
    screen_map_pairs = []
    for taxonomy_path in taxonomy_paths:
        mapped_path = map_for_screen(taxonomy_path)
        if mapped_path is None:
            raise ValueError(
                f'Could not associate {taxonomy_path} with a normalised_read_map.tsv. '
                'Select one partner/Paper Trail output folder at a time.'
            )
        screen_map_pairs.append((taxonomy_path.resolve(), mapped_path))

    for _taxonomy_path, map_path in screen_map_pairs:
        if map_path in used_map_paths:
            continue
        by_read, by_sequence = _read_normalised_read_map(map_path)
        for read_id, entry in by_read.items():
            previous = read_map_by_id.get(read_id)
            if previous and previous['sequence_id'] != entry['sequence_id']:
                raise ValueError(
                    f'Read ID {read_id!r} maps to multiple isolates across the selected onboarding outputs.'
                )
            read_map_by_id[read_id] = entry
        for sequence_id, ids in by_sequence.items():
            read_ids_by_sequence.setdefault(sequence_id, [])
            read_ids_by_sequence[sequence_id].extend(
                read_id for read_id in ids if read_id not in read_ids_by_sequence[sequence_id]
            )
        used_map_paths.add(map_path)
        loaded_files.append(str(map_path))

    for taxonomy_path, _map_path in screen_map_pairs:
        screen_dir = taxonomy_path.parent
        paper_trail_dir = screen_dir.parent if screen_dir.name == 'read_taxonomy_screen' else screen_dir
        fields, rows = _read_table_with_fields(taxonomy_path)
        if not taxonomy_fields:
            taxonomy_fields = fields
        elif taxonomy_fields != fields:
            raise ValueError('The selected onboarding taxonomy.tsv files have different column layouts.')
        for row in rows:
            read_id = _row_get(row, ('id', 'readid', 'read_id'))
            if not read_id:
                continue
            previous = taxonomy_by_read.get(read_id)
            if previous and previous != row:
                raise ValueError(
                    f'Multiple onboarding taxonomy assignments disagree for read {read_id!r}. '
                    'Select the intended onboarding run rather than combining different runs.'
                )
            taxonomy_by_read[read_id] = row
        loaded_files.append(str(taxonomy_path))

        hits_path = screen_dir / 'read_taxonomy_hits.tsv'
        if hits_path.is_file():
            fields, rows = _read_table_with_fields(hits_path)
            if not hit_fields:
                hit_fields = fields
            elif hit_fields != fields:
                raise ValueError('The selected onboarding read_taxonomy_hits.tsv files have different column layouts.')
            for row in rows:
                read_id = _row_get(row, ('readid', 'read_id'))
                if read_id:
                    hits_by_read.setdefault(read_id, []).append(row)
            loaded_files.append(str(hits_path))

        policy_path = screen_dir / 'taxonomy_screen_policy.tsv'
        if policy_path.is_file():
            fields, rows = _read_table_with_fields(policy_path)
            if not policy_fields:
                policy_fields = fields
            policy_rows.extend(rows)
            loaded_files.append(str(policy_path))

        trimmed_path = paper_trail_dir / 'trimmed_oriented_reads.fasta'
        if trimmed_path.is_file():
            trimmed_fasta_paths.append(str(trimmed_path.resolve()))
            try:
                for header, sequence in paper_trail.read_fasta(str(trimmed_path)):
                    read_id = str(header).split()[0]
                    sequence = ''.join(str(sequence).split()).upper()
                    previous = trimmed_sequences_by_read.get(read_id)
                    if previous and previous != sequence:
                        raise ValueError(
                            f'Conflicting onboarding trimmed sequences were found for read {read_id!r}. '
                            'Select the intended onboarding run.'
                        )
                    trimmed_sequences_by_read[read_id] = sequence
            except Exception as exc:
                raise ValueError(f'Could not read onboarding trimmed-read FASTA {trimmed_path}: {exc}') from exc
            loaded_files.append(str(trimmed_path.resolve()))

    mapped_taxonomy_reads = set(read_map_by_id) & set(taxonomy_by_read)
    if not mapped_taxonomy_reads:
        raise ValueError(
            'The normalised read map and taxonomy.tsv did not share any read IDs. '
            'Expected normalised_read_map.tsv file stems to match taxonomy.tsv ID values.'
        )

    # De-duplicate candidate hits while preserving their original table schema.
    for read_id, rows in list(hits_by_read.items()):
        seen = set()
        unique = []
        for row in rows:
            key = tuple(row.get(field, '') for field in hit_fields)
            if key not in seen:
                seen.add(key)
                unique.append(row)
        hits_by_read[read_id] = unique

    return {
        'root': str(root),
        'taxonomy_paths': [str(path) for path, _ in screen_map_pairs],
        'map_paths': sorted(str(path) for path in used_map_paths),
        'partner_names': sorted(
            Path(path).parent.name[:-len('_onboarding')]
            if Path(path).parent.name.endswith('_onboarding') else Path(path).parent.name
            for path in used_map_paths
        ),
        'onboarding_directory_count': len(used_map_paths),
        'paper_trail_directory_count': len(screen_map_pairs),
        'trimmed_fasta_paths': trimmed_fasta_paths,
        'trimmed_sequences_by_read': trimmed_sequences_by_read,
        'map_rows': list(read_map_by_id.values()),
        'loaded_files': loaded_files,
        'read_map_by_id': read_map_by_id,
        'read_ids_by_sequence': read_ids_by_sequence,
        'taxonomy_by_read': taxonomy_by_read,
        'hits_by_read': hits_by_read,
        'taxonomy_fields': taxonomy_fields,
        'hit_fields': hit_fields,
        'policy_rows': policy_rows,
        'policy_fields': policy_fields,
        'read_count': len(read_map_by_id),
        'taxonomy_read_count': len(mapped_taxonomy_reads),
        'hit_count': sum(len(rows) for rows in hits_by_read.values()),
        'isolate_count': len(read_ids_by_sequence),
    }


def _row_get(row: dict[str, str], names: tuple[str, ...]) -> str:
    values = {_normalise_key(key): str(value or '').strip() for key, value in row.items()}
    return next((values[name] for name in names if values.get(name)), '')


def _read_map(path: str | Path) -> list[dict[str, str]]:
    """Read a user-selected sample map; mapped AB1 paths may live anywhere."""
    map_path = Path(path).expanduser().resolve()
    if not map_path.is_file():
        raise ValueError(f'Sample map does not exist: {map_path}')
    if map_path.suffix.lower() not in _MAP_EXTENSIONS:
        raise ValueError(f'Sample map must be TSV, CSV or TXT: {map_path}')
    rows = []
    for line, source in enumerate(_read_tsv(map_path), start=2):
        sequence_id = _row_get(source, _MAP_ID_KEYS)
        read_file = _row_get(source, _MAP_FILE_KEYS)
        if not sequence_id or not read_file:
            if any(str(value or '').strip() for value in source.values()):
                raise ValueError(f'Sample map row {line} requires sequence_id and read_file/file/path columns')
            continue
        if any(separator in read_file for separator in (';', '|')):
            raise ValueError(
                f'Sample map row {line} contains multiple reads in one field. Use one row per physical read '
                'so each primer and direction can be inspected separately.'
            )
        read_path = Path(read_file).expanduser()
        if not read_path.is_absolute():
            read_path = map_path.parent / read_path
        read_path = read_path.resolve()
        if not read_path.is_file():
            raise ValueError(f'Read file on row {line} does not exist: {read_path}')
        if read_path.suffix.lower() not in _READ_EXTENSIONS:
            raise ValueError(f'Row {line} does not refer to an AB1/ABI file: {read_path}')
        primer = _row_get(source, _MAP_PRIMER_KEYS)
        direction = _row_get(source, _MAP_DIRECTION_KEYS).lower()
        if direction not in {'forward', 'reverse'}:
            direction = ''
        mode = _row_get(source, _MAP_MODE_KEYS).lower()
        if mode not in {'assemble', 'best_read'}:
            mode = ''
        rows.append({
            'sequence_id': sequence_id,
            'file': str(read_path),
            'primer': primer,
            'direction': direction,
            'processing_mode': mode,
            'dataset': _row_get(source, ('dataset',)),
            'partner_id': _row_get(source, ('partner_id', 'partner')),
        })
    if not rows:
        raise ValueError(f'Sample map contains no usable read rows: {map_path}')
    return rows


def _samples_from_rows(rows: list[dict[str, str]]) -> list[dict]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(row['sequence_id'], []).append(row)
    samples = []
    for sequence_id, reads in sorted(grouped.items()):
        modes = {row['processing_mode'] for row in reads if row.get('processing_mode')}
        mode = next(iter(modes)) if len(modes) == 1 else ('assemble' if len(reads) > 1 else 'best_read')
        samples.append({
            'sequence_id': sequence_id,
            'read_count': len(reads),
            'processing_mode': mode,
            'missing_files': sum(not Path(row['file']).is_file() for row in reads),
            'reads': reads,
        })
    return samples


def _write_map(path: Path, rows: list[dict[str, str]]) -> None:
    fields = ['sequence_id', 'file', 'primer', 'direction', 'processing_mode', 'dataset', 'partner_id']
    with open(path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def _load_review_decisions(path: Path) -> dict[str, dict[str, str]]:
    if not path.is_file():
        return {}
    decisions = {}
    for row in _read_tsv(path):
        sequence_id = _row_get(row, ('sequence_id', 'isolate_id', 'id'))
        if sequence_id:
            decisions[sequence_id] = {
                'sequence_id': sequence_id,
                'decision': _row_get(row, ('decision', 'review_status')),
                'reviewer': _row_get(row, ('reviewer',)),
                'notes': _row_get(row, ('notes', 'detail')),
            }
    return decisions


def _read_report_row(path: str | Path, sequence_id: str) -> dict[str, str]:
    """Return one isolate's assembly report row or fail loudly on malformed input."""
    rows = _read_tsv(path)
    if not rows:
        raise ValueError(f'Paper Trail assembly report contains no isolate rows: {path}')
    first = rows[0]
    if not any(_normalise_key(key) in {'sequenceid', 'id'} for key in first):
        raise ValueError(
            f'Could not parse the Paper Trail assembly report header in {path}. '
            'The file may have been read with the wrong delimiter; TSV files must remain tab-delimited.'
        )
    match = next((row for row in rows if _row_get(row, ('sequenceid', 'id')) == sequence_id), None)
    if match is None:
        available = [_row_get(row, ('sequenceid', 'id')) for row in rows[:8]]
        raise ValueError(
            f'Paper Trail assembly report {path} has no row for isolate {sequence_id!r}. '
            f'First available isolate IDs: {", ".join(value for value in available if value) or "none"}.'
        )

    # Compatibility for inspections written by the immediately preceding
    # Paper Trail build: its assembly_report.tsv header included SelectionWarning
    # but the row writer omitted that value, shifting every subsequent field by
    # one column. Detect only that exact signature and repair the row in memory.
    keys = list(match)
    selection_warning_index = keys.index('SelectionWarning') if 'SelectionWarning' in keys else -1
    fallback_status = str(match.get('FallbackReadOptions') or '').strip()
    if (
        selection_warning_index >= 0
        and fallback_status in {
            'assembled', 'partial_no_overlap', 'single_read', 'best_read',
            'failed_no_reads', 'filtered_output_length_lt_800', 'filtered_output_qc',
        }
        and str(match.get('AssemblyStatus') or '').strip().lower() in {'yes', 'no'}
    ):
        values = [match.get(key) for key in keys]
        repaired = values[:selection_warning_index] + [''] + values[selection_warning_index:-1]
        match = dict(zip(keys, repaired))
    return match


def _resolve_config_path(value: object, *, config_dir: Path, project_dir: Path, fallback_dir: Path) -> Path | None:
    raw = str(value or '').strip()
    if not raw:
        return None
    candidate = Path(raw).expanduser()
    if candidate.is_absolute():
        return candidate.resolve() if candidate.is_file() else None
    for base in (project_dir, config_dir, fallback_dir):
        resolved = (base / candidate).resolve()
        if resolved.is_file():
            return resolved
    return None


def _discover_taxonomy_screen_settings(root: Path) -> dict[str, str]:
    """Use the project's Assistant config as the default Studio taxonomy screen.

    Sanger Studio remains usable without a config. When the standard Rumen
    Gateway Assistant config is present in the launch directory, its primary
    reference and optional taxonomy table are reused so single-isolate
    inspections apply the same taxonomy screen as Assistant Paper Trail.
    """
    candidates = [
        root / 'Rumen_Gateway_Assistant.json',
        root / 'BranchManager.json',
        root / 'branchmanager.json',
    ]
    try:
        candidates.extend(sorted(root.glob('*Assistant*.json')))
    except OSError:
        pass
    seen = set()
    for config_path in candidates:
        config_path = config_path.resolve()
        if config_path in seen or not config_path.is_file():
            continue
        seen.add(config_path)
        try:
            config = json.loads(config_path.read_text(encoding='utf-8'))
        except (OSError, json.JSONDecodeError):
            continue
        defaults = config.get('defaults') if isinstance(config, dict) else {}
        assistant = defaults.get('assistant', {}) if isinstance(defaults, dict) else {}
        data = config.get('data', {}) if isinstance(config, dict) else {}
        if not isinstance(assistant, dict):
            assistant = {}
        if not isinstance(data, dict):
            data = {}
        project_raw = str(config.get('project_directory') or '').strip() if isinstance(config, dict) else ''
        project_dir = Path(project_raw).expanduser() if project_raw else config_path.parent
        if not project_dir.is_absolute():
            project_dir = (config_path.parent / project_dir).resolve()
        else:
            project_dir = project_dir.resolve()
        reference = _resolve_config_path(
            assistant.get('ref') or data.get('primary_taxonomy_fasta'),
            config_dir=config_path.parent, project_dir=project_dir, fallback_dir=root,
        )
        if not reference:
            continue
        taxa = _resolve_config_path(
            assistant.get('taxa') or data.get('primary_taxonomy_table'),
            config_dir=config_path.parent, project_dir=project_dir, fallback_dir=root,
        )
        return {
            'screen_ref': str(reference),
            'screen_taxa': str(taxa) if taxa else '',
            'config_path': str(config_path),
        }
    return {'screen_ref': '', 'screen_taxa': '', 'config_path': ''}


def _taxonomy_rank_value(taxonomy: str, rank: str) -> str:
    """Extract a GTDB-style rank value (for example ``g__Bacillus``)."""
    prefix = f'{rank}__'
    for item in str(taxonomy or '').split(';'):
        item = item.strip()
        if item.startswith(prefix):
            value = item[len(prefix):].strip()
            if value and value.lower() not in {'na', 'none', 'unknown', 'unclassified'}:
                return value
    return ''


def _derive_read_taxonomy_comparison(read_rows: list[dict[str, str]]) -> dict[str, str]:
    """Derive isolate-level concordance from per-read taxonomy.tsv rows.

    Assistant onboarding provides the per-read taxonomy table and read map;
    this calculation intentionally does not require a separate concordance TSV.
    It follows the same genus-first, family-fallback interpretation used by
    Paper Trail. A discordance with a single near-best hit is labelled potential.
    """
    genus = []
    family = []
    assignments = []
    for row in read_rows:
        read_id = _row_get(row, ('id', 'readid', 'read_id'))
        taxon = _row_get(row, ('taxon', 'taxonomy'))
        if not taxon or taxon.upper() == 'NA':
            taxon = ''
        if read_id:
            assignments.append(f'{read_id}={taxon or "unassigned"}')
        g = _taxonomy_rank_value(taxon, 'g')
        f = _taxonomy_rank_value(taxon, 'f')
        if read_id and g:
            genus.append((read_id, g))
        if read_id and f:
            family.append((read_id, f))

    genus_values = {taxon for _read_id, taxon in genus}
    family_values = {taxon for _read_id, taxon in family}
    if len(genus_values) > 1:
        compared = genus
        rank_name = 'genus'
        conflict = True
        reason = 'different genus assignments among pre-screened primer reads'
    elif len(family_values) > 1:
        compared = family
        rank_name = 'family'
        conflict = True
        reason = 'different family assignments among pre-screened primer reads'
    elif len(genus) >= 2:
        compared = genus
        rank_name = 'genus'
        conflict = False
        reason = f'concordant genus assignment: {genus[0][1]}'
    elif len(family) >= 2:
        compared = family
        rank_name = 'family'
        conflict = False
        reason = f'concordant family assignment: {family[0][1]}'
    else:
        compared = []
        rank_name = 'NA'
        conflict = False
        reason = 'fewer than two primer reads have genus/family assignments'

    def near_best_count(read_id: str) -> int:
        matching = next(
            (row for row in read_rows if _row_get(row, ('id', 'readid', 'read_id')) == read_id),
            {},
        )
        raw = _row_get(matching, ('nearbesthitcount', 'near_best_hit_count'))
        try:
            return int(float(raw or 0))
        except (TypeError, ValueError):
            return 0

    weak_support = bool(conflict) and any(near_best_count(read_id) < 2 for read_id, _taxon in compared)
    status = (
        'POTENTIAL_CONFLICT' if conflict and weak_support
        else 'CONFLICT' if conflict
        else 'CONCORDANT' if rank_name != 'NA'
        else 'INSUFFICIENT_CLASSIFICATION'
    )
    if conflict and weak_support:
        reason += '; at least one read has fewer than two near-best reference hits'
    return {
        'status': status,
        'compared_rank': rank_name,
        'assignments': ' | '.join(assignments),
        'reason': reason,
    }


def _attach_pre_screened_taxonomy(
    result: dict,
    sequence_id: str,
    job_dir: Path | None,
    output_root: Path | None,
    onboarding: dict[str, object],
    *,
    current_trimmed_fasta: str | Path | None = None,
) -> dict:
    """Attach cached Assistant onboarding taxonomy results to a fresh QC inspection.

    The only joins are:
      normalised_read_map.tsv: sequence_id + file basename -> read ID
      taxonomy.tsv: ID -> per-read assignment
      read_taxonomy_hits.tsv: ReadID -> candidate reference hits
    No concordance output file is needed or assumed.
    """
    taxonomy_by_read = onboarding.get('taxonomy_by_read') or {}
    hits_by_read = onboarding.get('hits_by_read') or {}
    read_map_by_id = onboarding.get('read_map_by_id') or {}
    expected_read_ids = list((onboarding.get('read_ids_by_sequence') or {}).get(sequence_id, []))
    current_read_ids = [
        str(row.get('ReadID') or '').strip()
        for row in result.get('reads', []) if str(row.get('ReadID') or '').strip()
    ]
    # Prefer actual IDs produced by the current inspection, but use the
    # normalised onboarding map to recover the intended linkage where useful.
    read_ids = current_read_ids or expected_read_ids
    if expected_read_ids and current_read_ids:
        expected = set(expected_read_ids)
        read_ids = [read_id for read_id in current_read_ids if read_id in expected]
        if not read_ids:
            read_ids = current_read_ids

    taxonomy_fields = list(onboarding.get('taxonomy_fields') or [])
    hit_fields = list(onboarding.get('hit_fields') or [])
    taxonomy_rows = []
    for read_id in read_ids:
        row = taxonomy_by_read.get(read_id)
        if row:
            taxonomy_rows.append(dict(row))
            continue
        # A read can be present in the normalised map but absent from the
        # taxonomy table (for example if it was omitted from an old screen run).
        if read_map_by_id.get(read_id):
            fallback = {field: '' for field in taxonomy_fields}
            fallback.update({
                'ID': read_id,
                'AssignedRank': 'NA',
                'Taxon': 'NA',
                'ReadClassificationStatus': 'NOT_FOUND_IN_ONBOARDING_TAXONOMY',
                'TaxonomySupport': 'READ_NOT_PRESENT_IN_TAXONOMY_TSV',
                'QualifiedHitCount': '0',
                'NearBestHitCount': '0',
            })
            taxonomy_rows.append(fallback)
    hit_rows = [
        dict(hit)
        for read_id in read_ids
        for hit in hits_by_read.get(read_id, [])
    ]
    comparison = _derive_read_taxonomy_comparison(taxonomy_rows)

    # Check the exact trimmed sequences screened during onboarding whenever its
    # trimmed FASTA is available. QC parameter changes can make cached taxonomy
    # stale even though the raw read IDs are unchanged.
    original_trimmed = onboarding.get('trimmed_sequences_by_read') or {}
    freshness = 'UNVERIFIED_TRIMMED_SEQUENCE_PROVENANCE'
    if current_trimmed_fasta and Path(current_trimmed_fasta).is_file() and original_trimmed:
        current_sequences = {}
        try:
            for header, sequence in paper_trail.read_fasta(str(current_trimmed_fasta)):
                current_sequences[str(header).split()[0]] = ''.join(str(sequence).split()).upper()
        except Exception:
            current_sequences = {}
        comparisons = []
        for read_id in read_ids:
            if read_id not in original_trimmed or read_id not in current_sequences:
                comparisons.append(False)
            else:
                comparisons.append(original_trimmed[read_id] == current_sequences[read_id])
        freshness = (
            'MATCHES_CURRENT_TRIMMED_READS'
            if comparisons and all(comparisons)
            else 'STALE_TRIMMED_SEQUENCE_CHANGED'
        )
    elif current_trimmed_fasta and not original_trimmed:
        freshness = 'UNVERIFIED_ONBOARDING_TRIMMED_FASTA_NOT_AVAILABLE'
    elif not current_trimmed_fasta:
        freshness = 'UNVERIFIED_CURRENT_TRIMMED_FASTA_NOT_AVAILABLE'

    # Materialise only the isolate-specific slice beside this inspection. This
    # keeps the review bundle self-contained without reloading large source TSVs.
    if job_dir is not None:
        screen_dir = job_dir / 'read_taxonomy_screen'
        screen_dir.mkdir(parents=True, exist_ok=True)
        if taxonomy_fields and taxonomy_rows:
            _write_tsv(screen_dir / 'taxonomy.tsv', taxonomy_fields, taxonomy_rows)
        if hit_fields:
            _write_tsv(screen_dir / 'read_taxonomy_hits.tsv', hit_fields, hit_rows)
        policy_fields = list(onboarding.get('policy_fields') or [])
        policy_rows = list(onboarding.get('policy_rows') or [])
        if policy_fields and policy_rows:
            _write_tsv(screen_dir / 'taxonomy_screen_policy.tsv', policy_fields, policy_rows[:1])
        source_manifest = screen_dir / 'taxonomy_source_manifest.tsv'
        _write_tsv(source_manifest, ['Key', 'Value'], [
            {'Key': 'Source', 'Value': 'Assistant onboarding pre-screened results'},
            {'Key': 'OnboardingRoot', 'Value': str(onboarding.get('root') or '')},
            {'Key': 'NormalisedReadMap', 'Value': ';'.join(onboarding.get('map_paths') or [])},
            {'Key': 'TaxonomyFiles', 'Value': ';'.join(onboarding.get('taxonomy_paths') or [])},
            {'Key': 'ReadIDsUsed', 'Value': ';'.join(read_ids)},
            {'Key': 'CurrentTrimmedFASTA', 'Value': str(current_trimmed_fasta or '')},
            {'Key': 'TaxonomyFreshness', 'Value': freshness},
            {'Key': 'VSEARCHRerun', 'Value': 'no'},
        ])
        if output_root is not None:
            known = {str(item.get('path') or '') for item in result.get('files', [])}
            for path, label in (
                (screen_dir / 'taxonomy.tsv', 'Pre-screened per-read taxonomy assignments'),
                (screen_dir / 'read_taxonomy_hits.tsv', 'Pre-screened candidate reference hits'),
                (screen_dir / 'taxonomy_screen_policy.tsv', 'Original onboarding taxonomy policy'),
                (source_manifest, 'Taxonomy source manifest'),
            ):
                if not path.is_file() or not _is_within(path, output_root):
                    continue
                relative = path.resolve().relative_to(output_root.resolve()).as_posix()
                if relative not in known:
                    result.setdefault('files', []).append({
                        'label': label,
                        'path': relative,
                        'url': '/artifact?path=' + urllib.parse.quote(relative),
                    })
                    known.add(relative)

    result['taxonomy_screen_requested'] = bool(taxonomy_rows)
    result['taxonomy_data_available'] = bool(taxonomy_rows)
    result['taxonomy_source'] = 'Assistant onboarding (pre-screened; no VSEARCH rerun)'
    result['taxonomy_source_root'] = str(onboarding.get('root') or '')
    result['taxonomy_freshness'] = freshness
    if freshness == 'STALE_TRIMMED_SEQUENCE_CHANGED':
        result['taxonomy_staleness_note'] = (
            'The current trimmed sequence differs from the sequence screened during onboarding. '
            'These assignments are historical evidence, not a fresh classification of the current reads. '
            'Enable the optional VSEARCH re-run to refresh them.'
        )
    elif freshness == 'MATCHES_CURRENT_TRIMMED_READS':
        result['taxonomy_staleness_note'] = 'The current trimmed read sequences match the sequences screened during Assistant onboarding.'
    else:
        result['taxonomy_staleness_note'] = (
            'The onboarding output did not contain enough trimmed-sequence provenance to verify freshness. '
            'Enable the optional VSEARCH re-run if QC/trimming settings have changed.'
        )
    result['taxonomy_concordance'] = comparison['status']
    result['read_taxonomy_assignments'] = comparison['assignments']
    result['taxonomy_comparison_reason'] = comparison['reason']
    result['taxonomy_reads'] = taxonomy_rows
    result['taxonomy_hits'] = hit_rows
    result['taxonomy_hit_count'] = len(hit_rows)
    result['taxonomy_read_ids'] = read_ids
    return result


def _parse_fallback_read_options(value: object) -> list[dict[str, object]]:
    """Parse the JSON-encoded QC-passing read candidates written by Paper Trail."""
    if isinstance(value, list):
        parsed = value
    else:
        raw = str(value or '').strip()
        if not raw:
            return []
        try:
            parsed = json.loads(raw)
        except (TypeError, ValueError, json.JSONDecodeError):
            return []
    if not isinstance(parsed, list):
        return []
    return [dict(item) for item in parsed if isinstance(item, dict) and str(item.get('read_id') or '').strip()]


def _serialise_inspection(output: dict, sequence_id: str, job_dir: Path, output_root: Path) -> dict:
    assembly = _read_report_row(output['assembly_tsv'], sequence_id)
    read_rows = [row for row in _read_tsv(output['read_qc_tsv']) if row.get('SequenceID') == sequence_id]
    placement_rows = [row for row in _read_tsv(output['assembly_placements_tsv']) if row.get('SequenceID') == sequence_id]
    overlap_rows = [row for row in _read_tsv(output['overlap_candidates_tsv']) if row.get('SequenceID') == sequence_id] if output.get('overlap_candidates_tsv') else []
    candidate_rows = [row for row in _read_tsv(output['assembly_candidates_tsv']) if row.get('SequenceID') == sequence_id] if output.get('assembly_candidates_tsv') else []
    consensus_evidence_rows = [row for row in _read_tsv(output['consensus_evidence_tsv']) if row.get('SequenceID') == sequence_id] if output.get('consensus_evidence_tsv') else []
    # Attach the exact candidate sequence from the FASTA artifact. The browser
    # must not reconstruct candidate sequences from summary metrics or tables.
    candidate_sequences: dict[str, str] = {}
    candidate_fasta = Path(str(output.get('candidate_assemblies_fasta') or ''))
    if candidate_fasta.is_file():
        current_id = ''
        chunks: list[str] = []
        with candidate_fasta.open() as handle:
            for line in handle:
                line = line.strip()
                if line.startswith('>'):
                    if current_id:
                        candidate_sequences[current_id] = ''.join(chunks)
                    current_id = line[1:].split()[0]
                    chunks = []
                elif current_id:
                    chunks.append(line)
        if current_id:
            candidate_sequences[current_id] = ''.join(chunks)
    for candidate in candidate_rows:
        candidate['Sequence'] = candidate_sequences.get(str(candidate.get('CandidateID') or ''), '')

    # Build an explicit, coordinate-preserving alignment for the browser.
    # Reads without a stored placement are deliberately not projected onto the
    # consensus: the UI must not imply support where the assembler recorded none.
    trimmed_sequences: dict[str, str] = {}
    trimmed_path = Path(str(output.get('trimmed_fasta') or ''))
    if trimmed_path.is_file():
        try:
            from branchmanager.utils.fasta import read_fasta
            trimmed_sequences = {str(header).split()[0]: str(seq).upper()
                                 for header, seq in read_fasta(str(trimmed_path))}
        except Exception:
            trimmed_sequences = {}
    read_alignment = []
    for placement in placement_rows:
        read_id = str(placement.get('ReadID') or '')
        sequence = trimmed_sequences.get(read_id, '')
        if not sequence:
            continue
        aligned = ['-'] * max(0, int(assembly.get('OutputLength') or 0))
        segment_text = str(placement.get('AlignedSegments') or '')
        for segment in segment_text.split(';'):
            match = re.fullmatch(r'(\d+)-(\d+):(\d+)-(\d+)', segment.strip())
            if not match:
                continue
            c_start, c_end, r_start, r_end = map(int, match.groups())
            c_start -= 1
            r_start -= 1
            span = min(c_end - c_start, r_end - r_start)
            for offset in range(max(0, span)):
                cp, rp = c_start + offset, r_start + offset
                if 0 <= cp < len(aligned) and 0 <= rp < len(sequence):
                    aligned[cp] = sequence[rp]
        if any(base != '-' for base in aligned):
            read_alignment.append({
                'read_id': read_id,
                'primer': placement.get('Primer', ''),
                'direction': placement.get('Direction', ''),
                'status': placement.get('ReadStatus', ''),
                'contributes': placement.get('ContributesToConsensus', ''),
                'aligned_sequence': ''.join(aligned),
            })
    # Decode compact, review-focused trace data for the reads in this isolate.
    # Full AB1 arrays are deliberately not embedded in the page: each channel is
    # max-pooled to a bounded number of points while preserving peak amplitudes.
    trace_review: dict[str, dict[str, object]] = {}
    try:
        for row in read_rows:
            read_id = str(row.get('ReadID') or '')
            source_file = str(row.get('SourceFile') or '').strip()
            if not read_id or not source_file:
                continue
            trace_path = Path(source_file).expanduser()
            if not trace_path.is_file():
                trace_review[read_id] = {'available': False, 'reason': 'Original AB1/ABI file is not accessible at its recorded path.'}
                continue
            try:
                details = paper_trail._read_ab1_details(trace_path)
                if not details.get('trace_available'):
                    trace_review[read_id] = {
                        'available': False,
                        'reason': 'AB1 trace channels or peak locations could not be decoded; base calls/qualities may still be available.',
                    }
                    continue
                channels = details.get('trace_channels') or {}
                trace_length = max((len(v) for v in channels.values()), default=0)
                max_points = 1800
                bin_size = max(1, (trace_length + max_points - 1) // max_points)
                pooled = {}
                for base in 'ACGT':
                    values = channels.get(base) or []
                    pooled[base] = [
                        max(values[i:i + bin_size], default=0)
                        for i in range(0, len(values), bin_size)
                    ]
                trace_review[read_id] = {
                    'available': True,
                    'trace_length': trace_length,
                    'bin_size': bin_size,
                    'channels': pooled,
                    'peak_locations': details.get('peak_locations') or [],
                    'trim_start': int(row.get('TrimStart') or 0),
                    'trim_end': int(row.get('TrimEnd') or 0),
                    'direction': str(row.get('Direction') or ''),
                    'raw_length': len(details.get('sequence') or ''),
                    'qualities': details.get('qualities') or [],
                    'secondary_peak_ratios': details.get('secondary_peak_ratios') or [],
                    'secondary_peak_bases': details.get('secondary_peak_bases') or [],
                    'sequence': details.get('sequence') or '',
                }
            except Exception as exc:
                trace_review[read_id] = {'available': False, 'reason': f'Trace parsing failed: {type(exc).__name__}'}
    except Exception:
        # Trace review is an optional evidence layer; it must never prevent
        # sequence/QC inspection from loading.
        trace_review = {}

    # Candidate-specific read projections use the same global alignment
    # implementation and scoring policy as the core assembler. This keeps the
    # browser's candidate review coordinates consistent with Paper Trail rather
    # than introducing a separate local-alignment interpretation.
    candidate_read_alignments: dict[str, list[dict[str, str]]] = {}
    try:
        read_metadata = {str(row.get('ReadID') or ''): row for row in read_rows}
        for candidate in candidate_rows:
            candidate_id = str(candidate.get('CandidateID') or '')
            candidate_seq = str(candidate.get('Sequence') or '').upper()
            if not candidate_id or not candidate_seq:
                continue
            selected_ids = [x for x in str(candidate.get('ReadIDs') or '').split(';') if x]
            projected = []
            for read_id in selected_ids:
                read_seq = trimmed_sequences.get(read_id, '')
                if not read_seq:
                    continue
                aligned = ['-'] * len(candidate_seq)
                read_positions = [-1] * len(candidate_seq)
                try:
                    # _align_read_to_consensus delegates to Paper Trail's shared
                    # _pairwise_aligner(): global mode, identical scores, free
                    # terminal gaps and N wildcard handling.
                    blocks = paper_trail._align_read_to_consensus(candidate_seq, read_seq)
                    for t0, t1, q0, q1 in blocks:
                        span = min(int(t1 - t0), int(q1 - q0))
                        for offset in range(span):
                            tp, qp = int(t0) + offset, int(q0) + offset
                            if 0 <= tp < len(aligned) and 0 <= qp < len(read_seq):
                                aligned[tp] = read_seq[qp]
                                read_positions[tp] = qp + 1
                except Exception:
                    continue
                if any(base != '-' for base in aligned):
                    meta = read_metadata.get(read_id, {})
                    projected.append({
                        'read_id': read_id,
                        'primer': str(meta.get('Primer') or ''),
                        'direction': str(meta.get('Direction') or ''),
                        'status': str(meta.get('QCClass') or ''),
                        'contributes': 'yes',
                        'aligned_sequence': ''.join(aligned),
                        'read_positions': read_positions,
                    })
            candidate_read_alignments[candidate_id] = projected
    except Exception:
        candidate_read_alignments = {}

    read_ids = {str(row.get('ReadID') or '') for row in read_rows}
    taxonomy_dir = job_dir / 'read_taxonomy_screen'
    taxonomy_summary_path = taxonomy_dir / 'taxonomy.tsv'
    taxonomy_hits_path = taxonomy_dir / 'read_taxonomy_hits.tsv'
    taxonomy_policy_path = taxonomy_dir / 'taxonomy_screen_policy.tsv'
    taxonomy_reads = [
        row for row in _read_tsv(taxonomy_summary_path)
        if str(row.get('ID') or '') in read_ids
    ] if taxonomy_summary_path.is_file() else []
    taxonomy_hits = [
        row for row in _read_tsv(taxonomy_hits_path)
        if str(row.get('ReadID') or '') in read_ids
    ] if taxonomy_hits_path.is_file() else []
    taxonomy_hits.sort(key=lambda row: (
        str(row.get('ReadID') or ''),
        str(row.get('NearBest') or '').lower() != 'yes',
        str(row.get('Qualified') or '').lower() != 'yes',
        -float(row.get('Identity') or 0.0),
    ))
    artifacts = []
    for key, label in (
        ('assembly_pngs', 'Assembly overview'),
        ('chromatogram_pngs', 'Chromatogram'),
        ('read_error_pngs', 'Read quality profile'),
    ):
        for item in output.get(key, []):
            path = Path(item).resolve()
            if _is_within(path, output_root):
                artifacts.append({'label': label, 'path': path.relative_to(output_root).as_posix(), 'url': '/artifact?path=' + urllib.parse.quote(path.relative_to(output_root).as_posix())})
    files = []
    for key, label in (
        ('assembly_tsv', 'Assembly report'),
        ('read_qc_tsv', 'Read QC table'),
        ('per_base_error_tsv', 'Per-base error table'),
        ('assembly_placements_tsv', 'Assembly placements'),
        ('overlap_candidates_tsv', 'Candidate read overlaps'),
        ('recommendations_tsv', 'Resequencing recommendations'),
        ('marker_review_template_tsv', 'Review template'),
        ('qc_policy_tsv', 'QC policy'),
        ('taxonomy_screen_tsv', 'Primer-read taxonomy concordance'),
        ('summary', 'Run summary'),
    ):
        raw_path = str(output.get(key) or '').strip()
        if not raw_path:
            continue
        path = Path(raw_path).resolve()
        if path.is_file() and _is_within(path, output_root):
            relative = path.relative_to(output_root).as_posix()
            files.append({'label': label, 'path': relative, 'url': '/artifact?path=' + urllib.parse.quote(relative)})
    for path, label in (
        (taxonomy_summary_path, 'Per-read taxonomy assignments'),
        (taxonomy_hits_path, 'Per-read taxonomy reference hits'),
        (taxonomy_policy_path, 'Taxonomy screen policy'),
        (taxonomy_dir / 'taxonomy_source_manifest.tsv', 'Taxonomy source manifest'),
    ):
        if path.is_file() and _is_within(path, output_root):
            relative = path.relative_to(output_root).as_posix()
            files.append({'label': label, 'path': relative, 'url': '/artifact?path=' + urllib.parse.quote(relative)})
    return {
        'sequence_id': sequence_id,
        'status': assembly.get('Status', ''),
        'qc_class': assembly.get('QCClass', ''),
        'recommendation': assembly.get('Recommendation', ''),
        'reasons': assembly.get('Reasons', ''),
        'output_length': assembly.get('OutputLength', ''),
        'overlap_length': assembly.get('OverlapLength', ''),
        'overlap_identity': assembly.get('OverlapIdentity', ''),
        'used_reads': assembly.get('UsedReads', ''),
        'read_count': assembly.get('ReadCount', ''),
        'processing_mode': assembly.get('ProcessingMode', ''),
        'assembly_status': assembly.get('AssemblyStatus') or assembly.get('Status', ''),
        'selected_read_id': assembly.get('SelectedReadID', ''),
        'read_selection_source': assembly.get('SelectionSource', ''),
        'fallback_reason': assembly.get('FallbackReason', ''),
        'recommended_read_id': assembly.get('RecommendedReadID', ''),
        'fallback_read_options': _parse_fallback_read_options(assembly.get('FallbackReadOptions')),
        'taxonomy_concordance': assembly.get('TaxonomyConcordance', ''),
        'read_taxonomy_assignments': assembly.get('ReadTaxonomyAssignments', ''),
        'taxonomy_screen_requested': bool(
            output.get('taxonomy_screen_tsv') or output.get('taxonomy_preloaded')
        ),
        'taxonomy_data_available': bool(taxonomy_reads),
        'taxonomy_source': output.get('taxonomy_source') or (
            'Fresh VSEARCH re-screen' if output.get('taxonomy_screen_tsv') else ''
        ),
        'taxonomy_staleness_note': output.get('taxonomy_staleness_note', ''),
        'taxonomy_comparison_reason': output.get('taxonomy_comparison_reason', ''),
        'taxonomy_reads': taxonomy_reads,
        'taxonomy_hits': taxonomy_hits,
        'taxonomy_hit_count': len(taxonomy_hits),
        'suggested_action': assembly.get('SuggestedAction', ''),
        'reads': read_rows,
        'placements': placement_rows,
        'overlap_candidates': overlap_rows,
        'assembly_candidates': candidate_rows,
        'consensus_evidence': consensus_evidence_rows,
        'read_alignment': read_alignment,
        'candidate_read_alignments': candidate_read_alignments,
        'trace_review': trace_review,
        'artifacts': artifacts,
        'files': files,
        'job_dir': str(job_dir),
    }


def _has_nonempty_manual_candidate(row: dict[str, str]) -> bool:
    """Determine whether a failed-QC review row describes a candidate sequence.

    New templates explicitly set eligible_for_manual_acceptance. Older saved
    Studio TSVs can contain a stale 'no' from the former length/N gate, so a
    positive output_length is also accepted as evidence that a candidate was
    generated. The downstream Assistant independently validates that the
    candidate is actually present in its prepared FASTA before resuming.
    """
    eligible = str(_row_get(row, ('eligible_for_manual_acceptance',))).strip().lower()
    if eligible in {'yes', 'true', '1'}:
        return True
    for length_field in ('candidate_length', 'output_length', 'length'):
        raw_length = _row_get(row, (length_field,))
        try:
            if float(raw_length) > 0:
                return True
        except (TypeError, ValueError):
            continue
    candidate_path = _row_get(row, ('candidate_fasta', 'fasta', 'sequence_fasta'))
    sequence_id = _row_get(row, ('sequence_id', 'id', 'isolate_id'))
    candidate = Path(candidate_path).expanduser() if candidate_path else None
    if candidate and candidate.is_file() and sequence_id:
        try:
            for header, sequence in paper_trail.read_fasta(str(candidate)):
                if str(header).split('|', 1)[0].split()[0] == sequence_id and str(sequence).strip():
                    return True
        except Exception:
            return False
    return False


def _review_decision_options(row: dict[str, str]) -> tuple[list[str], bool]:
    """Rebuild permitted decisions from the current QC/manual-review policy.

    In particular, do not trust stale allowed_decisions values stored by older
    Studio versions when eligible_for_manual_acceptance or output_length shows
    that a failed-QC candidate exists.
    """
    qc_class = str(_row_get(row, ('qc_class',))).upper()
    review_type = str(_row_get(row, ('review_type',))).upper()
    if qc_class == 'FAIL_QC':
        eligible = _has_nonempty_manual_candidate(row)
        options = (['APPROVE_FAIL_QC'] if eligible else []) + ['REJECT', 'RESEQUENCE']
        return options, eligible
    if qc_class == 'QUALITY_UNVERIFIED':
        eligible = str(_row_get(row, ('eligible_for_manual_acceptance',))).strip().lower() in {'yes', 'true', '1'}
        options = (['APPROVE_UNVERIFIED'] if eligible else []) + ['REJECT', 'RESEQUENCE']
        return options, eligible
    if qc_class == 'PASS_WITH_WARNINGS' or review_type == 'WARNING_REVIEW':
        return ['APPROVE', 'REJECT', 'RESEQUENCE'], True
    allowed = _row_get(row, ('allowed_decisions',))
    if allowed:
        return [value.strip().upper() for value in allowed.split('|') if value.strip()], False
    return ['REJECT', 'RESEQUENCE'], False


class SangerStudio:
    def __init__(self, root: str | Path | None = None, out: str | Path | None = None):
        # The launch directory is only a starting location for native file dialogs;
        # it is not a read-path security boundary.
        self.root = Path(root or Path.cwd()).expanduser().resolve()
        if not self.root.is_dir():
            raise ValueError(f'Studio starting directory is not an existing directory: {self.root}')
        self.out = Path(out or (self.root / 'sanger_studio_runs')).expanduser().resolve()
        self.out.mkdir(parents=True, exist_ok=True)
        taxonomy_defaults = _discover_taxonomy_screen_settings(self.root)
        self.screen_ref = taxonomy_defaults.get('screen_ref', '')
        self.screen_taxa = taxonomy_defaults.get('screen_taxa', '')
        self.taxonomy_config_path = taxonomy_defaults.get('config_path', '')
        self.jobs = self.out / 'inspections'
        self.jobs.mkdir(parents=True, exist_ok=True)
        self.review_path = self.out / 'marker_review_decisions.tsv'
        self.lock = threading.RLock()
        self.dialog_start = self.root
        self.map_path: Path | None = None
        self.map_rows: list[dict[str, str]] = []
        self.last_inspection: dict[str, dict] = {}
        self.last_batch: dict[str, dict] = {}
        self.decisions = _load_review_decisions(self.review_path)
        self.marker_review_source: Path | None = None
        self.marker_review_fields: list[str] = []
        self.marker_review_rows: list[dict[str, str]] = []
        self.marker_review_output_name = 'cohort_marker_review_decisions.tsv'
        self.onboarding_taxonomy: dict[str, object] | None = None
        self.pre_screen_taxonomy_by_isolate: dict[str, dict[str, str]] = {}
        self.input_source = 'none'
        self.input_directory = ''
        self.input_files: list[str] = []

    def choose_and_load_map(self) -> dict:
        selected = _native_dialog('map', initialdir=self.map_path.parent if self.map_path else self.dialog_start)
        if not selected:
            return {'cancelled': True}
        return self.load_map(selected)

    def choose_and_quick_map(self) -> dict:
        selected = _native_dialog('reads', initialdir=self.dialog_start)
        if not selected:
            return {'cancelled': True}
        return self.quick_map(list(selected))

    def choose_output_directory(self) -> dict:
        selected = _native_dialog('output', initialdir=self.out)
        if not selected:
            return {'cancelled': True}
        return self.set_output_directory(selected)

    def set_taxonomy_screen(self, screen_ref: str | Path | None, screen_taxa: str | Path | None = None) -> dict:
        reference_text = str(screen_ref or '').strip()
        taxa_text = str(screen_taxa or '').strip()
        reference = Path(reference_text).expanduser().resolve() if reference_text else None
        taxa = Path(taxa_text).expanduser().resolve() if taxa_text else None
        if reference_text and (reference is None or not reference.is_file()):
            raise ValueError(f'Taxonomy reference FASTA does not exist: {reference_text}')
        if taxa_text and (taxa is None or not taxa.is_file()):
            raise ValueError(f'Taxonomy mapping file does not exist: {taxa_text}')
        if not reference and taxa:
            raise ValueError('Select a taxonomy reference FASTA before supplying a taxonomy mapping file.')
        self.screen_ref = str(reference) if reference else ''
        self.screen_taxa = str(taxa) if taxa else ''
        return {
            'screen_ref': self.screen_ref,
            'screen_taxa': self.screen_taxa,
            'enabled': bool(self.screen_ref),
        }

    def choose_taxonomy_reference(self) -> dict:
        start = Path(self.screen_ref).parent if self.screen_ref else self.root
        selected = _native_dialog('taxonomy_ref', initialdir=start)
        if not selected:
            return {'cancelled': True}
        return self.set_taxonomy_screen(selected, self.screen_taxa)

    def choose_taxonomy_table(self) -> dict:
        start = Path(self.screen_taxa).parent if self.screen_taxa else (
            Path(self.screen_ref).parent if self.screen_ref else self.root
        )
        selected = _native_dialog('taxonomy_taxa', initialdir=start)
        if not selected:
            return {'cancelled': True}
        return self.set_taxonomy_screen(self.screen_ref, selected)

    def choose_and_load_input_folder(self) -> dict:
        """Choose a directory and automatically detect cohort, sample-map, or AB1/ABI input."""
        selected = _native_dialog('input_folder', initialdir=self.dialog_start)
        if not selected:
            return {'cancelled': True}
        return self.load_input_directory(selected)

    def choose_and_load_cohort_folder(self) -> dict:
        """Backward-compatible alias for the generic input-folder picker."""
        return self.choose_and_load_input_folder()

    def choose_and_load_onboarding_output(self) -> dict:
        """Backward-compatible alias for the generic input-folder picker."""
        return self.choose_and_load_input_folder()

    def load_uploaded_ab1_zip(self, archive_bytes: bytes, filename: str = 'sanger_reads.zip') -> dict:
        """Safely extract AB1/ABI files from a browser-uploaded ZIP and load them."""
        import io
        import shutil

        if not archive_bytes:
            raise ValueError('The uploaded ZIP file is empty.')
        if len(archive_bytes) > 500 * 1024 * 1024:
            raise ValueError('ZIP upload exceeds the 500 MB limit.')
        if not str(filename or '').lower().endswith('.zip'):
            raise ValueError('Please upload a .zip archive.')

        upload_root = self.out / 'input_uploads'
        upload_root.mkdir(parents=True, exist_ok=True)
        destination = upload_root / uuid.uuid4().hex
        destination.mkdir(parents=True, exist_ok=False)
        extracted_bytes = 0
        extracted_count = 0
        try:
            with zipfile.ZipFile(io.BytesIO(archive_bytes), 'r') as archive:
                members = archive.infolist()
                if len(members) > 20000:
                    raise ValueError('ZIP contains too many entries (maximum 20,000).')
                for member in members:
                    if member.is_dir():
                        continue
                    name = member.filename.replace('\\\\', '/')
                    parts = [part for part in name.split('/') if part not in {'', '.'}]
                    if not parts or any(part == '..' for part in parts):
                        raise ValueError(f'Unsafe path in ZIP archive: {member.filename!r}')
                    if Path(parts[-1]).suffix.lower() not in _READ_EXTENSIONS:
                        continue
                    extracted_bytes += int(member.file_size)
                    if extracted_bytes > 2 * 1024 * 1024 * 1024:
                        raise ValueError('Uncompressed chromatogram data exceeds the 2 GB limit.')
                    target = destination.joinpath(*parts)
                    target.parent.mkdir(parents=True, exist_ok=True)
                    with archive.open(member, 'r') as source, open(target, 'wb') as output:
                        while True:
                            chunk = source.read(1024 * 1024)
                            if not chunk:
                                break
                            output.write(chunk)
                    extracted_count += 1
                    if extracted_count > 20000:
                        raise ValueError('ZIP contains too many AB1/ABI files (maximum 20,000).')
            if not extracted_count:
                raise ValueError('No .ab1 or .abi chromatogram files were found in the ZIP archive.')
            result = self.load_input_directory(destination)
            result['input_kind'] = 'raw_chromatograms'
            result['uploaded_archive'] = Path(filename).name
            result['status_message'] = (
                f'Loaded {result.get("read_count", extracted_count)} AB1/ABI read(s) from '
                f'{Path(filename).name}. Check inferred isolate IDs, primers and directions before processing.'
            )
            return result
        except Exception:
            shutil.rmtree(destination, ignore_errors=True)
            raise

    def load_input_directory(self, path: str | Path) -> dict:
        """Detect and load the best Sanger input route in a selected directory.

        Priority:
          1. Assistant cohort/Paper Trail pair (read maps plus cached taxonomy).
          2. One or more normalised read maps without taxonomy.
          3. One or more standalone sample maps.
          4. Raw AB1/ABI chromatograms, converted to a quick map.

        No Assistant cohort mode or cached taxonomy is required.
        """
        selected = _normalise_input_directory(path)
        taxonomy_paths = _discover_taxonomy_tables(selected)
        normalised_maps = _discover_normalised_maps(selected)

        # A complete pre-screened cohort has priority, so map/taxonomy pairing
        # remains dataset-aware and is validated by the existing cohort loader.
        if taxonomy_paths:
            result = self.load_onboarding_output(selected)
            result.update({
                'input_kind': 'assistant_cohort',
                'input_directory': str(selected),
                'status_message': (
                    f'Assistant cohort detected: loaded {result.get("isolate_count", 0)} isolates, '
                    f'{result.get("mapped_read_count", 0)} reads and pre-screened taxonomy from '
                    f'{result.get("partner_count", 0)} paired partner dataset(s).'
                ),
            })
            return result

        # An Assistant-style read map can also be used without the cohort
        # taxonomy stage. This is a legitimate standalone route and remains
        # NOT_SCREENED until VSEARCH is deliberately requested.
        if normalised_maps:
            rows, source_paths = _combine_normalised_read_maps(normalised_maps)
            with self.lock:
                self.onboarding_taxonomy = None
                self.pre_screen_taxonomy_by_isolate = {}
                self.map_rows = [dict(row) for row in rows]
                self.map_path = Path(source_paths[0]).resolve() if len(source_paths) == 1 else None
                self.dialog_start = selected
                self.input_source = 'normalised_read_map_without_taxonomy'
                self.input_directory = str(selected)
                self.input_files = source_paths
            samples = _samples_from_rows(self.map_rows)
            return {
                'cancelled': False, 'input_kind': 'normalised_read_map',
                'input_directory': str(selected), 'path': str(self.map_path) if self.map_path else '',
                'map_paths': source_paths, 'samples': samples, 'read_count': len(rows),
                'taxonomy_loaded': False, 'taxonomy_status': 'NOT_SCREENED',
                'onboarding_taxonomy': self.onboarding_taxonomy_state(),
                'status_message': (
                    f'Loaded {len(samples)} isolates and {len(rows)} reads from '
                    f'{len(source_paths)} normalised read map(s). No pre-screened taxonomy was found '
                    '(NOT_SCREENED); optional VSEARCH can be configured below.'
                ),
            }

        standalone_maps = _discover_standalone_maps(selected)
        if standalone_maps:
            rows, source_paths = _combine_standalone_sample_maps(standalone_maps)
            with self.lock:
                self.onboarding_taxonomy = None
                self.pre_screen_taxonomy_by_isolate = {}
                self.map_rows = [dict(row) for row in rows]
                self.map_path = Path(source_paths[0]).resolve() if len(source_paths) == 1 else None
                self.dialog_start = selected
                self.input_source = 'standalone_sample_map'
                self.input_directory = str(selected)
                self.input_files = source_paths
            samples = _samples_from_rows(self.map_rows)
            return {
                'cancelled': False, 'input_kind': 'standalone_sample_map',
                'input_directory': str(selected), 'path': str(self.map_path) if self.map_path else '',
                'map_paths': source_paths, 'samples': samples, 'read_count': len(rows),
                'taxonomy_loaded': False, 'taxonomy_status': 'NOT_SCREENED',
                'onboarding_taxonomy': self.onboarding_taxonomy_state(),
                'status_message': (
                    f'Loaded {len(samples)} isolates and {len(rows)} reads from '
                    f'{len(source_paths)} standalone sample map(s). No Assistant onboarding taxonomy '
                    'was found (NOT_SCREENED); optional VSEARCH can be configured below.'
                ),
            }

        chromatograms = sorted(_iter_input_files(selected, {'.ab1', '.abi'}))
        if chromatograms:
            result = self.quick_map([str(item) for item in chromatograms])
            self.input_source = 'raw_chromatogram_directory'
            self.input_directory = str(selected)
            self.input_files = [str(item) for item in chromatograms]
            result.update({
                'cancelled': False, 'input_kind': 'raw_chromatograms',
                'input_directory': str(selected), 'taxonomy_loaded': False,
                'taxonomy_status': 'NOT_SCREENED',
                'status_message': (
                    f'No Assistant cohort or sample map was found. Created a quick map from '
                    f'{len(chromatograms)} AB1/ABI file(s) across {len(result.get("samples", []))} isolates. '
                    'Review inferred isolate IDs, primers and directions. Taxonomy is NOT_SCREENED; '
                    'optional VSEARCH can be configured below.'
                ),
            })
            return result

        raise ValueError(
            f'No recognised Sanger inputs were found beneath {selected}. Expected one of: '
            'an Assistant *_onboarding/*_paper_trail structure, normalised_read_map.tsv, '
            'a standalone TSV/CSV/TXT sample map with isolate and read-file columns, '
            'or AB1/ABI chromatograms.'
        )

    def load_onboarding_output(self, path: str | Path) -> dict:
        """Load pre-screened read taxonomy once from an Assistant output folder."""
        loaded = _load_onboarding_taxonomy(path)
        with self.lock:
            self.onboarding_taxonomy = loaded
            self.input_source = 'assistant_cohort'
            self.input_directory = str(Path(path).expanduser().resolve())
            self.input_files = list(loaded.get('loaded_files') or [])
            # The normalised onboarding map is also a valid Studio read map.
            # Populate it without rejecting records whose AB1 path is currently
            # unavailable; the sample list will mark missing files, while the
            # read-ID/taxonomy relationship remains available for review.
            self.map_rows = [dict(row) for row in (loaded.get('map_rows') or [])]
            map_paths = list(loaded.get('map_paths') or [])
            self.map_path = Path(map_paths[0]).resolve() if len(map_paths) == 1 else None
            self.dialog_start = Path(self.map_path).parent if self.map_path else Path(loaded['root'])
            self.pre_screen_taxonomy_by_isolate = {}
            for sequence_id, read_ids in (loaded.get('read_ids_by_sequence') or {}).items():
                rows = [
                    loaded.get('taxonomy_by_read', {}).get(read_id, {})
                    for read_id in read_ids
                    if loaded.get('taxonomy_by_read', {}).get(read_id)
                ]
                comparison = _derive_read_taxonomy_comparison(rows)
                self.pre_screen_taxonomy_by_isolate[sequence_id] = {
                    'TaxonomyConcordance': comparison['status'],
                    'TaxonomyComparedRank': comparison['compared_rank'],
                    'ReadTaxonomyAssignments': comparison['assignments'],
                    'Reason': comparison['reason'],
                    'ReadIDs': ';'.join(read_ids),
                }
            # Reuse the recorded reference path for optional re-screening only.
            # Loading/displaying cached assignments never opens the GTDB FASTA.
            policy_rows = loaded.get('policy_rows') or []
            if policy_rows:
                policy = policy_rows[0]
                reference_path = str(_row_get(policy, ('referencefasta', 'reference_fasta')) or '').strip()
                taxa_path = str(_row_get(policy, ('taxonomyfile', 'taxonomy_file')) or '').strip()
                if reference_path and Path(reference_path).is_file():
                    self.screen_ref = str(Path(reference_path).resolve())
                if taxa_path and taxa_path.lower() not in {'reference fasta headers', 'none', 'na'} and Path(taxa_path).is_file():
                    self.screen_taxa = str(Path(taxa_path).resolve())
                elif taxa_path.lower() == 'reference fasta headers':
                    self.screen_taxa = ''
        result = self.onboarding_taxonomy_state()
        result.update({
            'input_kind': 'assistant_cohort',
            'onboarding_taxonomy': self.onboarding_taxonomy_state(),
            'input_directory': str(Path(path).expanduser().resolve()),
            'path': str(self.map_path) if self.map_path else '',
            'samples': _samples_from_rows(self.map_rows),
            'read_count': len(self.map_rows),
            'taxonomy_loaded': True,
            'taxonomy_status': 'PRE_SCREENED_ONBOARDING',
        })
        return result

    # Backwards-compatible plural spelling used by earlier Studio development
    # builds and by saved local workflows.
    def load_onboarding_outputs(self, path: str | Path) -> dict:
        return self.load_onboarding_output(path)

    def _attach_prescreened_taxonomy(self, inspection: dict, sequence_id: str, output: dict) -> dict:
        """Attach cached onboarding assignments and verify trimmed-sequence freshness."""
        if not self.onboarding_taxonomy:
            return inspection
        job_value = output.get('job_dir') or output.get('_job_dir')
        job_dir = Path(job_value) if job_value else None
        return _attach_pre_screened_taxonomy(
            inspection,
            sequence_id,
            job_dir,
            self.out,
            self.onboarding_taxonomy,
            current_trimmed_fasta=output.get('trimmed_fasta'),
        )

    def onboarding_taxonomy_state(self) -> dict:
        loaded = self.onboarding_taxonomy
        if not loaded:
            return {
                'loaded': False, 'root': '', 'read_count': 0, 'taxonomy_read_count': 0,
                'hit_count': 0, 'isolate_count': 0, 'status': 'NOT_SCREENED',
                'message': (
                    'No pre-screened taxonomy is loaded. This is expected for standalone Sanger runs; '
                    'taxonomy can remain NOT_SCREENED or be generated with optional VSEARCH.'
                ),
            }
        return {
            'loaded': True,
            'root': str(loaded.get('root') or ''),
            'read_count': int(loaded.get('read_count') or 0),
            'taxonomy_read_count': int(loaded.get('taxonomy_read_count') or 0),
            'hit_count': int(loaded.get('hit_count') or 0),
            'taxonomy_hit_count': int(loaded.get('hit_count') or 0),
            'isolate_count': int(loaded.get('isolate_count') or 0),
            'mapped_read_count': int(loaded.get('read_count') or 0),
            'map_loaded': bool(self.map_rows),
            'map_path': str(self.map_path) if self.map_path else '',
            'input_source': self.input_source,
            'input_directory': self.input_directory,
            'input_files': list(self.input_files),
            'samples': _samples_from_rows(self.map_rows),
            'map_paths': list(loaded.get('map_paths') or []),
            'taxonomy_paths': list(loaded.get('taxonomy_paths') or []),
            'partner_count': len(loaded.get('map_paths') or []),
            'partner_names': list(loaded.get('partner_names') or []),
            'message': (
                f'Loaded {len(loaded.get("map_paths") or [])} paired onboarding/Paper Trail dataset(s). '
                'Read maps and pre-screened taxonomy are available together; routine inspections will not run VSEARCH.'
            ),
        }

    def choose_and_load_marker_review(self) -> dict:
        selected = _native_dialog('review', initialdir=self.marker_review_source.parent if self.marker_review_source else self.dialog_start)
        if not selected:
            return {'cancelled': True}
        return self.load_marker_review(selected)

    def set_output_directory(self, path: str | Path) -> dict:
        resolved = Path(path).expanduser().resolve()
        if not resolved.exists() or not resolved.is_dir():
            raise ValueError(f'Output directory does not exist or is not a directory: {resolved}')
        with self.lock:
            self.out = resolved
            self.jobs = self.out / 'inspections'
            self.jobs.mkdir(parents=True, exist_ok=True)
            self.review_path = self.out / 'marker_review_decisions.tsv'
            self.decisions = _load_review_decisions(self.review_path)
        return {'out': str(self.out), 'marker_review_output': str(self.out / self.marker_review_output_name)}

    def load_map(self, path: str | Path) -> dict:
        resolved = Path(path).expanduser().resolve()
        rows = _read_map(resolved)
        self.map_path = resolved
        self.map_rows = rows
        self.onboarding_taxonomy = None
        self.pre_screen_taxonomy_by_isolate = {}
        self.input_source = 'standalone_sample_map'
        self.input_directory = str(resolved.parent)
        self.input_files = [str(resolved)]
        self.dialog_start = resolved.parent
        return {
            'path': str(resolved), 'samples': _samples_from_rows(rows), 'read_count': len(rows),
            'input_kind': 'standalone_sample_map', 'taxonomy_loaded': False,
            'taxonomy_status': 'NOT_SCREENED', 'onboarding_taxonomy': self.onboarding_taxonomy_state(),
            'status_message': 'Standalone sample map loaded. No pre-screened taxonomy is available (NOT_SCREENED).',
        }

    def load_marker_review(self, path: str | Path) -> dict:
        source = Path(path).expanduser().resolve()
        fields, rows = _read_table_with_fields(source)
        required = {'sequence_id', 'decision', 'review_type', 'qc_class'}
        missing = sorted(required - {str(field).strip().lower() for field in fields})
        if missing:
            raise ValueError(
                'This is not a cohort marker review template. Missing columns: ' + ', '.join(missing)
            )
        seen = set()
        normalised = []
        for line, row in enumerate(rows, start=2):
            sid = _row_get(row, ('sequence_id',))
            if not sid:
                continue
            if sid in seen:
                raise ValueError(f'Duplicate sequence_id {sid!r} in marker review file (around line {line})')
            seen.add(sid)
            row['sequence_id'] = sid
            candidate_value = _row_get(row, ('candidate_fasta', 'fasta', 'sequence_fasta'))
            if candidate_value:
                candidate_path = Path(candidate_value).expanduser()
                if not candidate_path.is_absolute():
                    candidate_path = source.parent / candidate_path
                row['candidate_fasta'] = str(candidate_path.resolve())
            qc = str(_row_get(row, ('qc_class',))).upper()
            fallback_options = _parse_fallback_read_options(row.get('fallback_read_options'))
            row['fallback_read_options_list'] = fallback_options
            option_ids = {str(option.get('read_id') or '') for option in fallback_options}
            recommended_id = str(row.get('recommended_read_id') or (fallback_options[0].get('read_id') if fallback_options else '')).strip()
            selected_id = str(row.get('selected_read_id') or '').strip()
            selection_source = str(row.get('read_selection_source') or '').strip().upper()
            if selected_id and (selection_source != 'MANUAL' or selected_id not in option_ids):
                selected_id = ''
                selection_source = ''
            if not selected_id and recommended_id and option_ids:
                selected_id = recommended_id
                selection_source = 'AUTO_LONGEST_USABLE_READ'
            row['recommended_read_id'] = recommended_id
            row['selected_read_id'] = selected_id
            row['read_selection_source'] = selection_source
            decisions, eligible = _review_decision_options(row)
            if qc == 'FAIL_QC':
                # Upgrade stale saved review TSVs generated before the human
                # override policy. A positive candidate length is sufficient to
                # restore approval eligibility; Assistant resume still validates
                # the candidate FASTA independently.
                row['eligible_for_manual_acceptance'] = 'yes' if eligible else 'no'
                if eligible and 'manual_acceptance_basis' in row:
                    row['manual_acceptance_basis'] = 'HUMAN_REVIEW_OVERRIDE'
            row['allowed_decisions'] = '|'.join(decisions)
            row['_allowed_decisions'] = decisions
            row['_review_eligible'] = eligible
            row['_original_decision'] = str(row.get('decision') or '').strip().upper()
            normalised.append(row)
        if not normalised:
            raise ValueError(
                'Marker review file contains no sequence records. The file has a valid header but no review rows; '
                'check that Assistant PREPARE produced a populated cohort_marker_review_template.tsv, then rerun PREPARE.'
            )
        self.marker_review_source = source
        self.marker_review_fields = fields
        self.marker_review_rows = normalised
        self.dialog_start = source.parent
        return self.marker_review_state()

    def marker_review_state(self) -> dict:
        rows = []
        decided = 0
        for source_row in self.marker_review_rows:
            row = {key: value for key, value in source_row.items() if not key.startswith('_')}
            row['allowed_decisions_list'] = list(source_row.get('_allowed_decisions') or [])
            if str(row.get('decision') or '').strip():
                decided += 1
            rows.append(row)
        return {
            'source': str(self.marker_review_source) if self.marker_review_source else '',
            'output': str(self.out / self.marker_review_output_name),
            'total': len(rows),
            'decided': decided,
            'rows': rows,
        }

    def save_marker_review(self, data: dict) -> dict:
        """Validate and persist a complete review atomically.

        Work on a deep copy so an invalid row cannot partially mutate the
        in-memory review. The TSV is replaced only after every row validates
        and the temporary file has been written successfully.
        """
        if not isinstance(data, dict):
            raise ValueError('Invalid save request: expected a JSON object.')
        submitted = data.get('rows') or []
        if not isinstance(submitted, list):
            raise ValueError('Invalid save request: rows must be a list.')
        reviewer_default = str(data.get('reviewer') or '').strip()
        updates = {}
        for index, item in enumerate(submitted, start=1):
            if not isinstance(item, dict):
                raise ValueError(f'Invalid review row at submitted position {index}: expected an object.')
            sid = str(item.get('sequence_id') or '').strip()
            if not sid:
                raise ValueError(f'Invalid review row at submitted position {index}: sequence_id is missing.')
            if sid in updates:
                raise ValueError(f'Duplicate submitted review decision for {sid}')
            updates[sid] = item

        today = datetime.now().astimezone().date().isoformat()
        with self.lock:
            working_rows = copy.deepcopy(self.marker_review_rows)
            row_lookup = {str(row.get('sequence_id')): row for row in working_rows}
            unknown = sorted(set(updates) - set(row_lookup))
            if unknown:
                raise ValueError(
                    'Review contains sequence IDs not present in the attached template: '
                    + ', '.join(unknown[:10])
                )

            for sid, item in updates.items():
                row = row_lookup[sid]
                decision = str(item.get('decision') or '').strip().upper()
                reviewer = str(item.get('reviewer') or reviewer_default or row.get('reviewer') or '').strip()
                raw_notes = item.get('notes')
                notes = str(raw_notes if raw_notes is not None else (row.get('notes') or '')).strip()
                manual_selection_changed = False

                if 'selected_read_id' in item:
                    selected_read_id = str(item.get('selected_read_id') or '').strip()

                    # An inspection may have been run after the cohort review TSV
                    # was attached. In that workflow the browser receives the
                    # fallback options from the inspection response and adds them
                    # to its local review row, but the server's review row still
                    # reflects the original PREPARE template (which may not contain
                    # any fallback columns). Use the latest server-generated
                    # inspection as the authoritative source, and persist its
                    # evidence into the working review row before validating the
                    # submitted choice. Never trust candidate IDs supplied by the
                    # browser itself.
                    inspection = self.last_inspection.get(sid)
                    if inspection is not None:
                        fallback_options = _parse_fallback_read_options(
                            inspection.get('fallback_read_options')
                        )
                        recommended_read_id = str(
                            inspection.get('recommended_read_id') or ''
                        ).strip()
                        row['fallback_read_options_list'] = fallback_options
                        row['fallback_read_options'] = json.dumps(
                            fallback_options, ensure_ascii=False, separators=(',', ':')
                        ) if fallback_options else ''
                        row['recommended_read_id'] = recommended_read_id
                        row['fallback_reason'] = str(inspection.get('fallback_reason') or '')
                        row['processing_mode'] = str(inspection.get('processing_mode') or row.get('processing_mode') or '')
                        row['assembly_status'] = str(inspection.get('assembly_status') or row.get('assembly_status') or '')
                    else:
                        fallback_options = row.get('fallback_read_options_list') or _parse_fallback_read_options(
                            row.get('fallback_read_options')
                        )
                        recommended_read_id = str(row.get('recommended_read_id') or '').strip()

                    option_ids = {
                        str(option.get('read_id') or '').strip() for option in fallback_options
                    }
                    previous_selected = str(row.get('selected_read_id') or '').strip()
                    previous_source = str(row.get('read_selection_source') or '').strip().upper()
                    if selected_read_id and selected_read_id not in option_ids:
                        raise ValueError(
                            f'{sid}: selected fallback read {selected_read_id!r} is not among '
                            'the reads that passed individual QC.'
                        )
                    if selected_read_id and selected_read_id != recommended_read_id:
                        if not reviewer:
                            raise ValueError(
                                f'{sid}: reviewer identity is required when selecting a '
                                'non-default fallback read.'
                            )
                        row['selected_read_id'] = selected_read_id
                        row['read_selection_source'] = 'MANUAL'
                        row['reviewer'] = reviewer
                        manual_selection_changed = (
                            previous_selected != selected_read_id or previous_source != 'MANUAL'
                        )
                    elif selected_read_id:
                        row['selected_read_id'] = selected_read_id
                        row['read_selection_source'] = 'AUTO_LONGEST_USABLE_READ'
                        manual_selection_changed = previous_source == 'MANUAL'
                    else:
                        row['selected_read_id'] = (
                            recommended_read_id if recommended_read_id in option_ids else ''
                        )
                        row['read_selection_source'] = (
                            'AUTO_LONGEST_USABLE_READ' if row['selected_read_id'] else ''
                        )
                        manual_selection_changed = previous_source == 'MANUAL'

                if decision:
                    allowed = row.get('_allowed_decisions') or []
                    if decision not in allowed:
                        raise ValueError(
                            f'{sid}: {decision} is not permitted. Allowed: '
                            f'{", ".join(allowed) or "none"}'
                        )
                    qc_class = str(row.get('qc_class') or '').strip().upper()
                    expected_qc = {
                        'APPROVE': 'PASS_WITH_WARNINGS',
                        'APPROVE_FAIL_QC': 'FAIL_QC',
                        'APPROVE_UNVERIFIED': 'QUALITY_UNVERIFIED',
                    }.get(decision)
                    if expected_qc and qc_class != expected_qc:
                        raise ValueError(
                            f'{sid}: {decision} requires qc_class {expected_qc}, '
                            f'not {qc_class or "unknown"}.'
                        )
                    if not reviewer:
                        raise ValueError(
                            f'{sid}: reviewer is required when recording a decision. '
                            'Enter a reviewer name or initials above the review table.'
                        )
                    if decision.startswith('APPROVE') and not notes:
                        raise ValueError(
                            f'{sid}: explanatory notes are required for an approval. '
                            'Add a short rationale in the Reviewer notes column.'
                        )
                    if decision.startswith('APPROVE') and str(
                        row.get('eligible_for_manual_acceptance') or ''
                    ).strip().lower() not in {'yes', 'true', '1'}:
                        raise ValueError(f'{sid}: this record is not eligible for manual acceptance.')
                    row['decision'] = decision
                    row['reviewer'] = reviewer
                    row['notes'] = notes
                    if 'reviewed_at' in self.marker_review_fields and (
                        decision != str(row.get('_original_decision') or '').upper()
                        or manual_selection_changed
                    ):
                        row['reviewed_at'] = today
                    row['_original_decision'] = decision
                else:
                    # A blank decision may still carry reviewer notes as a saved draft.
                    row['decision'] = ''
                    row['reviewer'] = reviewer or str(row.get('reviewer') or '').strip()
                    row['notes'] = notes
                    if 'reviewed_at' in self.marker_review_fields:
                        if manual_selection_changed:
                            row['reviewed_at'] = (
                                today if row.get('read_selection_source') == 'MANUAL' else ''
                            )
                        elif row.get('read_selection_source') != 'MANUAL':
                            row['reviewed_at'] = ''
                    row['_original_decision'] = ''

            fields = list(self.marker_review_fields)
            for required_field in (
                'sequence_id', 'decision', 'reviewer', 'reviewed_at', 'notes',
                'eligible_for_manual_acceptance', 'allowed_decisions',
                'fallback_read_options', 'recommended_read_id', 'selected_read_id',
                'read_selection_source', 'fallback_reason', 'processing_mode',
                'assembly_status', 'output_type',
            ):
                if required_field not in fields:
                    fields.append(required_field)

            target = self.out / self.marker_review_output_name
            if self.marker_review_source and target.resolve() == self.marker_review_source.resolve():
                suffix = datetime.now().astimezone().strftime('%Y%m%dT%H%M%S')
                target = self.out / f'cohort_marker_review_decisions_{suffix}.tsv'
            target.parent.mkdir(parents=True, exist_ok=True)
            temp = target.with_name(target.name + f'.{threading.get_ident()}.tmp')
            try:
                with open(temp, 'w', newline='', encoding='utf-8') as handle:
                    writer = csv.DictWriter(
                        handle, fieldnames=fields, delimiter='\t',
                        extrasaction='ignore', lineterminator='\n',
                    )
                    writer.writeheader()
                    for row in working_rows:
                        row['allowed_decisions'] = '|'.join(row.get('_allowed_decisions') or [])
                        if str(row.get('qc_class') or '').upper() == 'FAIL_QC':
                            row['eligible_for_manual_acceptance'] = (
                                'yes' if row.get('_review_eligible') else 'no'
                            )
                            if row.get('_review_eligible') and 'manual_acceptance_basis' in row:
                                row['manual_acceptance_basis'] = 'HUMAN_REVIEW_OVERRIDE'
                        output_row = {field: row.get(field, '') for field in fields}
                        writer.writerow(output_row)
                temp.replace(target)
            except OSError as exc:
                try:
                    temp.unlink() if temp.exists() else None
                except OSError:
                    pass
                raise ValueError(
                    f'Could not write the review file to {target}: {exc}. '
                    'Check that the Sanger Studio output directory is writable and that the disk has space.'
                ) from exc
            except Exception:
                try:
                    temp.unlink() if temp.exists() else None
                except OSError:
                    pass  # cleanup of temp file failed
                raise

            # Commit the new in-memory state only after the TSV has been replaced.
            self.marker_review_rows = working_rows

        result = self.marker_review_state()
        result['saved'] = str(target)
        result['download_url'] = '/artifact?path=' + urllib.parse.quote(target.name)
        return result
    def state(self) -> dict:
        policy = paper_trail.DEFAULT_QC_POLICY
        return {
            'root': str(self.root),
            'out': str(self.out),
            'map_path': str(self.map_path) if self.map_path else '',
            'samples': _samples_from_rows(self.map_rows),
            'screen_ref': self.screen_ref,
            'screen_taxa': self.screen_taxa,
            'taxonomy_config_path': self.taxonomy_config_path,
            'taxonomy_screen_enabled': bool(self.screen_ref),
            'onboarding_taxonomy': self.onboarding_taxonomy_state(),
            'review_decisions': list(self.decisions.values()),
            'marker_review': self.marker_review_state() if self.marker_review_rows else None,
            'policy': {
                'min_quality': policy.min_quality,
                'min_read_length': policy.min_read_length,
                'min_final_length': policy.min_final_length,
                'min_mean_quality': policy.min_mean_quality,
                'max_output_expected_errors': policy.max_output_expected_errors,
                'warn_n_percent': policy.warn_n_percent,
                'max_n_percent': policy.max_n_percent,
                'min_overlap': policy.min_overlap,
                'min_overlap_identity': policy.min_overlap_identity,
            },
        }

    def browse(self, relative: str = '.') -> dict:
        path = (self.root / relative).resolve()
        if not _is_within(path, self.root) or not path.is_dir():
            raise ValueError('Browse path must be a directory inside the configured root')
        entries = []
        try:
            children = sorted(path.iterdir(), key=lambda item: (not item.is_dir(), item.name.lower()))
        except OSError as exc:
            raise ValueError(f'Could not read directory: {exc}') from exc
        for child in children[:5000]:
            if child.name.startswith('.'):
                continue
            if child.is_dir():
                entries.append({'name': child.name, 'path': child.relative_to(self.root).as_posix(), 'kind': 'directory'})
            elif child.suffix.lower() in _MAP_EXTENSIONS:
                entries.append({'name': child.name, 'path': child.relative_to(self.root).as_posix(), 'kind': 'map'})
            elif child.suffix.lower() in _READ_EXTENSIONS:
                entries.append({'name': child.name, 'path': child.relative_to(self.root).as_posix(), 'kind': 'read'})
        return {'current': path.relative_to(self.root).as_posix() or '.', 'entries': entries}

    def quick_map(self, paths: list[str]) -> dict:
        if not paths:
            raise ValueError('Select at least one AB1/ABI file')
        rows = []
        for raw_path in paths:
            path = Path(raw_path).expanduser().resolve()
            if not path.is_file() or path.suffix.lower() not in _READ_EXTENSIONS:
                raise ValueError(f'Not an existing AB1/ABI file: {path}')
            sequence_id, primer, direction = paper_trail.infer_read_info(path.name, paper_trail.DEFAULT_PRIMERS)
            rows.append({
                'sequence_id': sequence_id,
                'file': str(path),
                'primer': '' if primer == 'unknown' else primer,
                'direction': direction if primer != 'unknown' else '',
                'processing_mode': '',
                'dataset': 'SangerStudio',
                'partner_id': '',
            })
        grouped: dict[str, list[dict]] = {}
        for row in rows:
            grouped.setdefault(row['sequence_id'], []).append(row)
        for group in grouped.values():
            mode = 'assemble' if len(group) > 1 else 'best_read'
            for row in group:
                row['processing_mode'] = mode
        self.map_path = None
        self.map_rows = rows
        self.onboarding_taxonomy = None
        self.pre_screen_taxonomy_by_isolate = {}
        self.input_source = 'selected_raw_chromatograms'
        self.input_directory = str(Path(rows[0]['file']).parent) if rows else ''
        self.input_files = [str(row['file']) for row in rows]
        if rows:
            self.dialog_start = Path(rows[0]['file']).parent
        return {
            'path': '', 'samples': _samples_from_rows(rows), 'read_count': len(rows),
            'input_kind': 'raw_chromatograms', 'taxonomy_loaded': False,
            'taxonomy_status': 'NOT_SCREENED', 'onboarding_taxonomy': self.onboarding_taxonomy_state(),
            'status_message': 'Raw AB1/ABI files loaded. No pre-screened taxonomy is available (NOT_SCREENED).',
        }

    def _persist_review_decisions(self) -> None:
        fields = ['sequence_id', 'decision', 'reviewer', 'notes']
        with open(self.review_path, 'w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
            writer.writeheader()
            writer.writerows(self.decisions[key] for key in sorted(self.decisions))

    def save_review(self, data: dict) -> dict:
        sequence_id = str(data.get('sequence_id') or '').strip()
        decision = str(data.get('decision') or '').strip().lower()
        reviewer = str(data.get('reviewer') or '').strip()
        notes = str(data.get('notes') or '').strip()
        inspection = self.last_inspection.get(sequence_id)
        if not inspection:
            raise ValueError(f'Run QC for {sequence_id!r} before recording a review')
        if inspection.get('qc_class') != 'PASS_WITH_WARNINGS' or inspection.get('recommendation') != 'MANUAL_REVIEW':
            raise ValueError('Manual approval/rejection is only available for PASS_WITH_WARNINGS markers. FAIL_QC markers require resequencing.')
        if decision not in _REVIEW_DECISIONS:
            raise ValueError('Decision must be approve or reject')
        if decision in {'approve', 'accepted', 'accept', 'pass', 'yes'}:
            decision = 'approve'
        else:
            decision = 'reject'
        if not reviewer:
            raise ValueError('Enter the reviewer name or initials')
        self.decisions[sequence_id] = {
            'sequence_id': sequence_id,
            'decision': decision,
            'reviewer': reviewer,
            'notes': notes,
        }
        self._persist_review_decisions()
        return {'path': str(self.review_path), 'decisions': list(self.decisions.values())}

    def inspect(self, data: dict) -> dict:
        if not self.map_rows:
            raise ValueError('Load a Mailroom sample map or create a quick map from selected AB1 files first')
        sequence_id = str(data.get('sequence_id') or '').strip()
        sample_reads = [dict(row) for row in self.map_rows if row['sequence_id'] == sequence_id]
        if not sample_reads:
            raise ValueError(f'Isolate {sequence_id!r} is not in the loaded map')
        overrides = data.get('reads') or []
        override_by_file = {str(row.get('file') or ''): row for row in overrides}
        for row in sample_reads:
            update = override_by_file.get(row['file'], {})
            primer = str(update.get('primer') or '').strip()
            direction = str(update.get('direction') or '').strip().lower()
            if primer:
                row['primer'] = primer
            if direction:
                if direction not in {'forward', 'reverse'}:
                    raise ValueError(f'Invalid read direction for {row["file"]}: {direction}')
                row['direction'] = direction
        mode = str(data.get('processing_mode') or 'map').lower()
        if mode not in {'map', 'assemble', 'best_read'}:
            raise ValueError('processing_mode must be map, assemble, or best_read')
        if mode != 'map':
            for row in sample_reads:
                row['processing_mode'] = mode
        manual_read_selection = {}
        review_row = next(
            (row for row in self.marker_review_rows if str(row.get('sequence_id') or '') == sequence_id),
            None,
        ) if self.marker_review_rows else None
        if review_row and str(review_row.get('read_selection_source') or '').strip().upper() == 'MANUAL':
            selected_read_id = str(review_row.get('selected_read_id') or '').strip()
            if selected_read_id:
                manual_read_selection[sequence_id] = selected_read_id
        # Keep corrections in the active Studio session so batch processing
        # uses the same primer, direction, and mode choices as individual QC.
        for stored in self.map_rows:
            updated = next((row for row in sample_reads if row['file'] == stored['file']), None)
            if updated:
                stored.update({key: updated[key] for key in ('primer', 'direction', 'processing_mode')})
        if any(not row.get('primer') or row.get('direction') not in {'forward', 'reverse'} for row in sample_reads):
            raise ValueError('Every selected read needs a primer name and forward/reverse direction. Correct the mapping before QC.')

        allowed_settings = {
            'min_quality': int, 'min_read_length': int, 'min_final_length': int,
            'min_mean_quality': float, 'max_output_expected_errors': float,
            'warn_n_percent': float, 'max_n_percent': float,
            'min_overlap': int, 'min_overlap_identity': float,
        }
        defaults = self.state()['policy']
        submitted = data.get('settings') or {}
        settings = {}
        for name, kind in allowed_settings.items():
            value = submitted.get(name, defaults[name])
            try:
                settings[name] = kind(value)
            except (TypeError, ValueError) as exc:
                raise ValueError(f'Invalid QC setting {name}: {value!r}') from exc
        sequence_type = str(data.get('sequence_type') or '16S').strip().upper()
        if sequence_type not in _SEQUENCE_TYPE_PROFILES:
            raise ValueError(f'Unsupported sequence type profile: {sequence_type!r}')
        custom_primer_sequences = _parse_custom_primer_sequences(data.get('custom_primer_sequences'))
        if not (100 <= settings['min_read_length'] <= 5000):
            raise ValueError('Minimum usable read length must be between 100 and 5000 bp')
        if not (100 <= settings['min_final_length'] <= 10000):
            raise ValueError('Minimum final consensus length must be between 100 and 10000 bp')
        if not (0 <= settings['min_quality'] <= 60):
            raise ValueError('Phred trimming cutoff must be between 0 and 60')
        if not (0 <= settings['min_overlap_identity'] <= 1):
            raise ValueError('Minimum overlap identity must be between 0 and 1')
        if not (0 <= settings['warn_n_percent'] <= settings['max_n_percent'] <= 100):
            raise ValueError('N-percent warning/failure thresholds must satisfy 0 <= warning <= failure <= 100')
        if settings['min_overlap'] < 1:
            raise ValueError('Minimum overlap must be at least 1 bp')

        rerun_taxonomy = bool(data.get('rerun_taxonomy'))
        if rerun_taxonomy and not self.screen_ref:
            raise ValueError(
                'VSEARCH re-screening was requested, but no taxonomy reference FASTA is configured. '
                'Load onboarding results for the existing assignments, or configure a reference in the optional re-screening section.'
            )

        with self.lock:
            stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
            job_dir = self.jobs / f'{stamp}_{re.sub(r"[^A-Za-z0-9_.-]+", "_", sequence_id)}_{uuid.uuid4().hex[:6]}'
            job_dir.mkdir(parents=True, exist_ok=False)
            map_path = job_dir / 'selected_read_map.tsv'
            _write_map(map_path, sample_reads)
            manifest = RunManifest(job_dir, 'sanger_studio_inspection')
            manifest.data['settings'] = dict(
                settings, processing_mode=mode, sequence_id=sequence_id, sequence_type=sequence_type, custom_primer_sequences=custom_primer_sequences,
                taxonomy_screen_requested=rerun_taxonomy,
                manual_fallback_read_id=manual_read_selection.get(sequence_id, ''),
            )
            if self.map_path and self.map_path.is_file():
                manifest.add_input(self.map_path, role='source_sample_map')
            if self.review_path.is_file():
                manifest.add_input(self.review_path, role='manual_review_decisions')
            if rerun_taxonomy and self.screen_ref:
                manifest.add_input(self.screen_ref, role='taxonomy_screen_reference')
            if rerun_taxonomy and self.screen_taxa:
                manifest.add_input(self.screen_taxa, role='taxonomy_screen_mapping')
            manifest.data['taxonomy_screen'] = {
                'mode': 'fresh_vsearch' if rerun_taxonomy else (
                    'assistant_onboarding_cache' if self.onboarding_taxonomy else 'disabled'
                ),
                'enabled': rerun_taxonomy,
                'reference': self.screen_ref if rerun_taxonomy else None,
                'taxonomy_table': self.screen_taxa or None if rerun_taxonomy else None,
                'pre_screened_source': (self.onboarding_taxonomy or {}).get('root') if not rerun_taxonomy else None,
            }
            if not rerun_taxonomy and self.onboarding_taxonomy:
                for source_path in self.onboarding_taxonomy.get('loaded_files', []):
                    if Path(source_path).is_file():
                        manifest.add_input(source_path, role='pre_screened_onboarding_taxonomy')
            manifest.add_input(map_path, role='resolved_read_map')
            for row in sample_reads:
                manifest.add_input(row['file'], role=f'raw_trace:{row["sequence_id"]}')
            manifest.write()
            manifest.add_stage('paper_trail', 'RUNNING', detail='Per-read QC, trimming, primer orientation and read inspection')
            if rerun_taxonomy:
                manifest.add_stage(
                    'taxonomy_screen', 'RUNNING',
                    detail=f'Fresh VSEARCH screening of current trimmed reads against {Path(self.screen_ref).name}',
                )
            try:
                output = paper_trail.run_paper_trail(
                    [], job_dir, sample_map=map_path,
                    primers=paper_trail.DEFAULT_PRIMERS,
                    primer_sequences=custom_primer_sequences,
                    min_quality=settings['min_quality'],
                    min_read_length=settings['min_read_length'],
                    min_length=settings['min_final_length'],
                    min_mean_quality=settings['min_mean_quality'],
                    max_output_expected_errors=settings['max_output_expected_errors'],
                    warn_n_percent=settings['warn_n_percent'],
                    max_n_percent=settings['max_n_percent'],
                    min_overlap=settings['min_overlap'],
                    min_overlap_identity=settings['min_overlap_identity'],
                    screen_ref=self.screen_ref or None if rerun_taxonomy else None,
                    screen_taxa=self.screen_taxa or None if rerun_taxonomy else None,
                    threads=4,
                    selected_read_ids=manual_read_selection,
                    assemble=True,
                )
                taxonomy_summary_path = job_dir / 'read_taxonomy_screen' / 'taxonomy.tsv'
                taxonomy_hits_path = job_dir / 'read_taxonomy_screen' / 'read_taxonomy_hits.tsv'
                taxonomy_assignment_count = max(0, len(_read_tsv(taxonomy_summary_path))) if taxonomy_summary_path.is_file() else 0
                taxonomy_hit_count = max(0, len(_read_tsv(taxonomy_hits_path))) if taxonomy_hits_path.is_file() else 0
                if rerun_taxonomy:
                    manifest.add_stage(
                        'taxonomy_screen', 'COMPLETE',
                        detail=(f'Fresh VSEARCH completed against {Path(self.screen_ref).name}: '
                                f'{taxonomy_assignment_count} read assignment(s), {taxonomy_hit_count} candidate hit(s)'),
                    )
                for key in (
                    'assembled_fasta', 'read_qc_tsv', 'per_base_error_tsv', 'assembly_tsv',
                    'assembly_placements_tsv', 'recommendations_tsv', 'marker_review_template_tsv',
                    'failed_review_candidates_fasta', 'qc_policy_tsv', 'summary',
                ):
                    manifest.add_output(output[key], role=key)
                for key in ('read_error_pngs', 'chromatogram_pngs', 'assembly_pngs'):
                    for artifact in output.get(key, []):
                        manifest.add_output(artifact, role=key)
                manifest.add_stage('paper_trail', 'COMPLETE')
                manifest.add_stage('merge_meeting', 'COMPLETE', detail='Used existing quality-aware assembler and consensus-level QC')
                manifest.finish('COMPLETE')
            except Exception as exc:
                manifest.finish('FAILED', error=exc)
                raise

        result = _serialise_inspection(output, sequence_id, job_dir, self.out)
        if rerun_taxonomy:
            result['taxonomy_source'] = 'Fresh VSEARCH re-screen using the current trimmed reads'
            result['taxonomy_data_available'] = bool(result.get('taxonomy_reads'))
            result['taxonomy_freshness'] = 'FRESH_VSEARCH_RERUN'
            result['taxonomy_staleness_note'] = ''
            result['taxonomy_screen_message'] = (
                f'Fresh VSEARCH completed against {Path(self.screen_ref).name}: '
                f'{taxonomy_assignment_count} per-read assignment(s), {taxonomy_hit_count} candidate hit(s).'
            )
        elif self.onboarding_taxonomy:
            output['_job_dir'] = str(job_dir)
            result = self._attach_prescreened_taxonomy(result, sequence_id, output)
        else:
            result['taxonomy_source'] = ''
            result['taxonomy_data_available'] = False
        self.last_inspection[sequence_id] = result
        return result

    def process_all(self, data: dict) -> dict:
        """Run Paper Trail once across the loaded map and create a download bundle."""
        if not self.map_rows:
            raise ValueError('Load a Mailroom sample map or create a quick map from selected AB1 files first')
        missing = [row['file'] for row in self.map_rows if not Path(row['file']).is_file()]
        if missing:
            raise ValueError(f'{len(missing)} mapped read file(s) are missing; first missing file: {missing[0]}')
        if any(not row.get('primer') or row.get('direction') not in {'forward', 'reverse'} for row in self.map_rows):
            raise ValueError('Every read needs a primer name and forward/reverse direction before batch processing')

        active_id = str(data.get('sequence_id') or '').strip()
        if active_id and data.get('reads'):
            overrides = {str(row.get('file') or ''): row for row in data['reads']}
            for row in self.map_rows:
                if row['sequence_id'] != active_id:
                    continue
                update = overrides.get(row['file'], {})
                primer = str(update.get('primer') or '').strip()
                direction = str(update.get('direction') or '').strip().lower()
                if primer:
                    row['primer'] = primer
                if direction in {'forward', 'reverse'}:
                    row['direction'] = direction
            mode = str(data.get('processing_mode') or 'map').lower()
            if mode in {'assemble', 'best_read'}:
                for row in self.map_rows:
                    if row['sequence_id'] == active_id:
                        row['processing_mode'] = mode

        allowed_settings = {
            'min_quality': int, 'min_read_length': int, 'min_final_length': int,
            'min_mean_quality': float, 'max_output_expected_errors': float,
            'warn_n_percent': float, 'max_n_percent': float,
            'min_overlap': int, 'min_overlap_identity': float,
        }
        defaults = self.state()['policy']
        submitted = data.get('settings') or {}
        settings = {}
        for name, kind in allowed_settings.items():
            try:
                settings[name] = kind(submitted.get(name, defaults[name]))
            except (TypeError, ValueError) as exc:
                raise ValueError(f'Invalid QC setting {name}: {submitted.get(name)!r}') from exc
        sequence_type = str(data.get('sequence_type') or '16S').strip().upper()
        if sequence_type not in _SEQUENCE_TYPE_PROFILES:
            raise ValueError(f'Unsupported sequence type profile: {sequence_type!r}')
        custom_primer_sequences = _parse_custom_primer_sequences(data.get('custom_primer_sequences'))
        if not (100 <= settings['min_read_length'] <= 5000):
            raise ValueError('Minimum usable read length must be between 100 and 5000 bp')
        if not (100 <= settings['min_final_length'] <= 10000):
            raise ValueError('Minimum final sequence length must be between 100 and 10000 bp')
        if not (0 <= settings['min_quality'] <= 60):
            raise ValueError('Phred trimming cutoff must be between 0 and 60')
        if not (0 <= settings['min_overlap_identity'] <= 1):
            raise ValueError('Minimum overlap identity must be between 0 and 1')
        if settings['min_overlap'] < 1:
            raise ValueError('Minimum read overlap must be at least 1 bp')
        if not (0 <= settings['warn_n_percent'] <= settings['max_n_percent'] <= 100):
            raise ValueError('N-percent warning/failure thresholds must satisfy 0 <= warning <= failure <= 100')

        rerun_taxonomy = bool(data.get('rerun_taxonomy'))
        if rerun_taxonomy and not self.screen_ref:
            raise ValueError(
                'VSEARCH re-screening was requested, but no taxonomy reference FASTA is configured. '
                'Load onboarding results for the existing assignments, or configure a reference in the optional re-screening section.'
            )

        with self.lock:
            stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
            batch_name = f'{stamp}_batch_{uuid.uuid4().hex[:6]}'
            job_dir = self.jobs / batch_name
            job_dir.mkdir(parents=True, exist_ok=False)
            map_path = job_dir / 'resolved_read_map.tsv'
            _write_map(map_path, self.map_rows)
            sample_ids = {row['sequence_id'] for row in self.map_rows}
            batch_decisions = [self.decisions[key] for key in sorted(sample_ids) if key in self.decisions]
            batch_review_path = job_dir / 'marker_review_decisions.tsv'
            if batch_decisions:
                with open(batch_review_path, 'w', newline='', encoding='utf-8') as handle:
                    writer = csv.DictWriter(handle, fieldnames=['sequence_id', 'decision', 'reviewer', 'notes'], delimiter='\t')
                    writer.writeheader()
                    writer.writerows(batch_decisions)
            manual_read_selections = {
                str(row.get('sequence_id') or ''): str(row.get('selected_read_id') or '').strip()
                for row in self.marker_review_rows
                if str(row.get('read_selection_source') or '').strip().upper() == 'MANUAL'
                and str(row.get('selected_read_id') or '').strip()
            } if self.marker_review_rows else {}
            manifest = RunManifest(job_dir, 'sanger_studio_batch')
            manifest.data['settings'] = dict(
                settings, sequence_type=sequence_type, custom_primer_sequences=custom_primer_sequences,
                manual_fallback_read_selections=manual_read_selections,
            )
            manifest.data['sample_count'] = len({row['sequence_id'] for row in self.map_rows})
            if self.map_path and self.map_path.is_file():
                manifest.add_input(self.map_path, role='source_sample_map')
            if self.review_path.is_file():
                manifest.add_input(self.review_path, role='manual_review_decisions')
            if rerun_taxonomy and self.screen_ref:
                manifest.add_input(self.screen_ref, role='taxonomy_screen_reference')
            if rerun_taxonomy and self.screen_taxa:
                manifest.add_input(self.screen_taxa, role='taxonomy_screen_mapping')
            manifest.data['taxonomy_screen'] = {
                'mode': 'fresh_vsearch' if rerun_taxonomy else (
                    'assistant_onboarding_cache' if self.onboarding_taxonomy else 'disabled'
                ),
                'enabled': rerun_taxonomy,
                'reference': self.screen_ref if rerun_taxonomy else None,
                'taxonomy_table': (self.screen_taxa or None) if rerun_taxonomy else None,
                'pre_screened_source': (self.onboarding_taxonomy or {}).get('root') if not rerun_taxonomy else None,
            }
            if not rerun_taxonomy and self.onboarding_taxonomy:
                for source_path in self.onboarding_taxonomy.get('loaded_files', []):
                    if Path(source_path).is_file():
                        manifest.add_input(source_path, role='pre_screened_onboarding_taxonomy')
            manifest.add_input(map_path, role='resolved_read_map')
            for row in self.map_rows:
                manifest.add_input(row['file'], role=f'raw_trace:{row["sequence_id"]}')
            manifest.write()
            manifest.add_stage('paper_trail', 'RUNNING', detail='Batch base-calling, trimming, QC and assembly')
            if rerun_taxonomy:
                manifest.add_stage(
                    'taxonomy_screen', 'RUNNING',
                    detail=f'Fresh VSEARCH screening of current trimmed reads against {Path(self.screen_ref).name}',
                )
            try:
                output = paper_trail.run_paper_trail(
                    [], job_dir, sample_map=map_path,
                    primers=paper_trail.DEFAULT_PRIMERS,
                    primer_sequences=custom_primer_sequences,
                    min_quality=settings['min_quality'],
                    min_read_length=settings['min_read_length'],
                    min_length=settings['min_final_length'],
                    min_mean_quality=settings['min_mean_quality'],
                    max_output_expected_errors=settings['max_output_expected_errors'],
                    warn_n_percent=settings['warn_n_percent'],
                    max_n_percent=settings['max_n_percent'],
                    min_overlap=settings['min_overlap'],
                    min_overlap_identity=settings['min_overlap_identity'],
                    screen_ref=self.screen_ref or None if rerun_taxonomy else None,
                    screen_taxa=self.screen_taxa or None if rerun_taxonomy else None,
                    selected_read_ids=manual_read_selections,
                    threads=4,
                    assemble=True,
                )
                taxonomy_summary_path = job_dir / 'read_taxonomy_screen' / 'taxonomy.tsv'
                taxonomy_hits_path = job_dir / 'read_taxonomy_screen' / 'read_taxonomy_hits.tsv'
                taxonomy_assignment_count = max(0, len(_read_tsv(taxonomy_summary_path))) if taxonomy_summary_path.is_file() else 0
                taxonomy_hit_count = max(0, len(_read_tsv(taxonomy_hits_path))) if taxonomy_hits_path.is_file() else 0
                if rerun_taxonomy:
                    manifest.add_stage(
                        'taxonomy_screen', 'COMPLETE',
                        detail=(f'Fresh VSEARCH completed against {Path(self.screen_ref).name}: '
                                f'{taxonomy_assignment_count} read assignment(s), {taxonomy_hit_count} candidate hit(s)'),
                    )
                batch_taxonomy_freshness = 'NOT_USED'
                if not rerun_taxonomy and self.onboarding_taxonomy:
                    batch_taxonomy_freshness = self._write_pre_screened_batch_taxonomy(
                        job_dir, current_trimmed_fasta=output.get('trimmed_fasta'),
                    )
                if batch_decisions:
                    manifest.add_output(batch_review_path, role='marker_review_decisions')
                manifest.add_stage('paper_trail', 'COMPLETE')
                manifest.add_stage('merge_meeting', 'COMPLETE', detail='Batch assembly and consensus-level QC complete')
                manifest.finish('COMPLETE')
                bundle = self.out / f'{batch_name}.zip'
                manifest.write()
                with zipfile.ZipFile(bundle, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
                    for path in sorted(job_dir.rglob('*')):
                        if path.is_file():
                            archive.write(path, path.relative_to(job_dir).as_posix())
            except Exception as exc:
                manifest.finish('FAILED', error=exc)
                raise

        report_rows = _read_tsv(output['assembly_tsv'])
        qc_counts: dict[str, int] = {}
        for row in report_rows:
            qc_class = row.get('QCClass') or 'UNKNOWN'
            qc_counts[qc_class] = qc_counts.get(qc_class, 0) + 1
        result = {
            'job_dir': str(job_dir),
            'bundle': bundle.name,
            'bundle_url': '/artifact?path=' + urllib.parse.quote(bundle.name),
            'sample_count': len({row['sequence_id'] for row in self.map_rows}),
            'output': {key: output[key] for key in (
                'assembled_fasta', 'assembly_tsv', 'read_qc_tsv', 'per_base_error_tsv',
                'recommendations_tsv', 'marker_review_template_tsv', 'failed_review_candidates_fasta',
                'taxonomy_screen_tsv', 'summary',
            )},
            'qc_counts': qc_counts,
            'taxonomy_mode': 'fresh_vsearch' if rerun_taxonomy else (
                'assistant_onboarding_cache' if self.onboarding_taxonomy else 'disabled'
            ),
            'taxonomy_freshness': (
                'FRESH_VSEARCH_RERUN' if rerun_taxonomy else batch_taxonomy_freshness
            ),
            'taxonomy_screen_message': (
                f'Fresh VSEARCH completed against {Path(self.screen_ref).name}: '
                f'{taxonomy_assignment_count} read assignment(s), {taxonomy_hit_count} candidate hit(s).'
                if rerun_taxonomy else ''
            ),
        }
        self.last_batch[batch_name] = result
        return result

    def _write_pre_screened_batch_taxonomy(
        self, job_dir: Path, *, current_trimmed_fasta: str | Path | None = None,
    ) -> str:
        """Copy pre-screened taxonomy into the batch bundle and check sequence freshness."""
        if not self.onboarding_taxonomy:
            return 'NOT_LOADED'
        onboarding = self.onboarding_taxonomy
        read_ids = {
            _read_id_from_path(row.get('file'))
            for row in self.map_rows if row.get('file')
        }
        taxonomy_rows = [
            dict(onboarding['taxonomy_by_read'][read_id])
            for read_id in sorted(read_ids)
            if read_id in onboarding.get('taxonomy_by_read', {})
        ]
        hit_rows = [
            dict(hit)
            for read_id in sorted(read_ids)
            for hit in onboarding.get('hits_by_read', {}).get(read_id, [])
        ]
        original_trimmed = onboarding.get('trimmed_sequences_by_read') or {}
        current_trimmed = {}
        if current_trimmed_fasta and Path(current_trimmed_fasta).is_file():
            try:
                for header, sequence in paper_trail.read_fasta(str(current_trimmed_fasta)):
                    current_trimmed[str(header).split()[0]] = ''.join(str(sequence).split()).upper()
            except Exception:
                current_trimmed = {}
        if not original_trimmed or not current_trimmed:
            freshness = 'UNVERIFIED_TRIMMED_SEQUENCE_PROVENANCE'
        else:
            comparisons = [
                read_id in original_trimmed
                and read_id in current_trimmed
                and original_trimmed[read_id] == current_trimmed[read_id]
                for read_id in read_ids
            ]
            freshness = (
                'MATCHES_CURRENT_TRIMMED_READS'
                if comparisons and all(comparisons)
                else 'STALE_TRIMMED_SEQUENCE_CHANGED'
            )

        screen_dir = job_dir / 'read_taxonomy_screen'
        if taxonomy_rows and onboarding.get('taxonomy_fields'):
            _write_tsv(screen_dir / 'taxonomy.tsv', list(onboarding['taxonomy_fields']), taxonomy_rows)
        if onboarding.get('hit_fields'):
            _write_tsv(screen_dir / 'read_taxonomy_hits.tsv', list(onboarding['hit_fields']), hit_rows)
        if onboarding.get('policy_fields') and onboarding.get('policy_rows'):
            _write_tsv(screen_dir / 'taxonomy_screen_policy.tsv', list(onboarding['policy_fields']), list(onboarding['policy_rows'])[:1])
        _write_tsv(screen_dir / 'taxonomy_source_manifest.tsv', ['Key', 'Value'], [
            {'Key': 'Source', 'Value': 'Assistant onboarding pre-screened results'},
            {'Key': 'OnboardingRoot', 'Value': str(onboarding.get('root') or '')},
            {'Key': 'VSEARCHRerun', 'Value': 'no'},
            {'Key': 'ReadIDsIncluded', 'Value': str(len(taxonomy_rows))},
            {'Key': 'CurrentTrimmedFASTA', 'Value': str(current_trimmed_fasta or '')},
            {'Key': 'TaxonomyFreshness', 'Value': freshness},
        ])
        return freshness


def serve(*, port: int = 8765, open_browser: bool = True) -> None:
    """Start Sanger Studio; all input and output selection happens in the browser UI."""
    from importlib.resources import files
    start_dir = Path.cwd().resolve()
    studio = SangerStudio(start_dir, start_dir / 'sanger_studio_runs')
    html = files('branchmanager').joinpath('data/sanger_studio.html').read_text(encoding='utf-8')

    class Handler(BaseHTTPRequestHandler):
        server_version = 'BranchManagerSangerStudio/1'

        def _send(self, status: int, body: bytes, content_type: str) -> None:
            self.send_response(status)
            self.send_header('Content-Type', content_type)
            self.send_header('Content-Length', str(len(body)))
            self.send_header('Cache-Control', 'no-store')
            self.end_headers()
            self.wfile.write(body)

        def _json(self, status: int, value: dict) -> None:
            self._send(status, json.dumps(value, ensure_ascii=False).encode('utf-8'), 'application/json; charset=utf-8')

        def _body(self) -> dict:
            size = int(self.headers.get('Content-Length') or '0')
            if size > 2_000_000:
                raise ValueError('Request is too large')
            raw = self.rfile.read(size)
            return json.loads(raw.decode('utf-8') or '{}')

        def do_GET(self):  # noqa: N802 - stdlib handler name
            parsed = urllib.parse.urlparse(self.path)
            params = urllib.parse.parse_qs(parsed.query)
            if parsed.path == '/':
                self._send(200, html.encode('utf-8'), 'text/html; charset=utf-8')
                return
            try:
                if parsed.path == '/api/state':
                    self._json(200, studio.state())
                elif parsed.path == '/api/browse':
                    self._json(200, studio.browse((params.get('path') or ['.'])[0]))
                elif parsed.path == '/artifact':
                    relative = (params.get('path') or [''])[0]
                    path = (studio.out / relative).resolve()
                    if not _is_within(path, studio.out) or not path.is_file():
                        raise ValueError('Artifact not found or is outside the Studio output directory')
                    content_type = mimetypes.guess_type(path.name)[0] or 'application/octet-stream'
                    self._send(200, path.read_bytes(), content_type)
                else:
                    self._json(404, {'error': 'Not found'})
            except Exception as exc:
                self._json(400, {'error': str(exc)})

        def do_POST(self):  # noqa: N802 - stdlib handler name
            try:
                if self.path == '/api/upload-ab1-zip':
                    size = int(self.headers.get('Content-Length') or '0')
                    if size <= 0 or size > 500 * 1024 * 1024:
                        raise ValueError('ZIP upload must be between 1 byte and 500 MB.')
                    content_type = self.headers.get('Content-Type', '').split(';', 1)[0].strip().lower()
                    if content_type not in {'application/zip', 'application/x-zip-compressed', 'application/octet-stream'}:
                        raise ValueError('Upload must be a ZIP archive.')
                    archive_bytes = self.rfile.read(size)
                    if len(archive_bytes) != size:
                        raise ValueError('ZIP upload was incomplete.')
                    filename = urllib.parse.unquote(
                        self.headers.get('X-Filename', 'sanger_reads.zip')
                    )
                    self._json(200, studio.load_uploaded_ab1_zip(archive_bytes, filename))
                    return
                payload = self._body()
                if self.path == '/api/pick-map':
                    self._json(200, studio.choose_and_load_map())
                elif self.path == '/api/load-map':
                    self._json(200, studio.load_map(payload.get('path', '')))
                elif self.path == '/api/pick-reads':
                    self._json(200, studio.choose_and_quick_map())
                elif self.path == '/api/quick-map':
                    self._json(200, studio.quick_map(payload.get('paths') or []))
                elif self.path == '/api/pick-output':
                    self._json(200, studio.choose_output_directory())
                elif self.path == '/api/set-output':
                    self._json(200, studio.set_output_directory(payload.get('path', '')))
                elif self.path == '/api/pick-taxonomy-reference':
                    self._json(200, studio.choose_taxonomy_reference())
                elif self.path == '/api/pick-taxonomy-table':
                    self._json(200, studio.choose_taxonomy_table())
                elif self.path == '/api/set-taxonomy-screen':
                    self._json(200, studio.set_taxonomy_screen(
                        payload.get('screen_ref', ''), payload.get('screen_taxa', ''),
                    ))
                elif self.path == '/api/pick-input-folder':
                    self._json(200, studio.choose_and_load_input_folder())
                elif self.path == '/api/load-input-folder':
                    self._json(200, studio.load_input_directory(payload.get('path', '')))
                elif self.path == '/api/pick-cohort-folder':
                    self._json(200, studio.choose_and_load_input_folder())
                elif self.path == '/api/load-cohort-folder':
                    self._json(200, studio.load_onboarding_output(payload.get('path', '')))
                elif self.path == '/api/pick-onboarding-results':
                    self._json(200, studio.choose_and_load_cohort_folder())
                elif self.path == '/api/load-onboarding-results':
                    self._json(200, studio.load_onboarding_output(payload.get('path', '')))
                elif self.path == '/api/pick-marker-review':
                    self._json(200, studio.choose_and_load_marker_review())
                elif self.path == '/api/load-marker-review':
                    self._json(200, studio.load_marker_review(payload.get('path', '')))
                elif self.path == '/api/save-marker-review':
                    self._json(200, studio.save_marker_review(payload))
                elif self.path == '/api/inspect':
                    self._json(200, studio.inspect(payload))
                elif self.path == '/api/process-all':
                    self._json(200, studio.process_all(payload))
                elif self.path == '/api/review':
                    self._json(200, studio.save_review(payload))
                elif self.path == '/api/shutdown':
                    if self.client_address[0] not in {'127.0.0.1', '::1'}:
                        self._json(403, {'error': 'Shutdown is available only from the local machine.'})
                        return
                    self._json(200, {
                        'stopping': True,
                        'message': 'Sanger Studio is shutting down. Your output files and saved decisions remain on disk.',
                    })
                    # shutdown() must run outside the serve_forever thread.
                    threading.Thread(
                        target=self.server.shutdown,
                        name='sanger-studio-shutdown', daemon=True,
                    ).start()
                else:
                    self._json(404, {'error': 'Not found'})
            except Exception as exc:
                if self.path == '/api/save-marker-review':
                    print(
                        f'[sanger-studio] save-marker-review rejected: '
                        f'{type(exc).__name__}: {exc}',
                        flush=True,
                    )
                self._json(400, {'error': str(exc)})

        def log_message(self, format, *args):  # noqa: A003 - stdlib handler signature
            print(f'[sanger-studio] {self.address_string()} - {format % args}')

    try:
        server = ThreadingHTTPServer(('127.0.0.1', int(port)), Handler)
    except OSError as exc:
        raise SystemExit(f'[sanger-studio] Could not bind to 127.0.0.1:{port}: {exc}') from exc
    url = f'http://127.0.0.1:{server.server_port}/'
    print(f'[sanger-studio] Local GUI: {url}\n[sanger-studio] Starting directory: {studio.root}\n[sanger-studio] Choose an input folder, standalone map, or AB1/ABI files in the browser.\nUse Shut down Studio in the page or press Ctrl+C to stop.')
    if open_browser:
        try:
            webbrowser.open(url)
        except Exception:
            logging.getLogger(__name__).debug('Could not open browser at %s', url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print('\n[sanger-studio] Stopped from the browser or terminal.')
    finally:
        server.server_close()