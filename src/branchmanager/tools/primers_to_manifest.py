#!/usr/bin/env python3
"""Copy Sanger trace files from a primer-mapping Excel sheet and create a manifest.

Output trace names are standardised to:
    ISOLATE_ID_PRIMER.ab1
Repeated isolate/primer combinations receive unique _attemptNN suffixes.

The manifest retains source paths, submission-folder grouping, checksums and
the workbook cell/column provenance for every mapped chromatogram. Generic
"Other" columns use the primer designation inferred from the source filename
when possible; unresolved labels remain explicit for review.

The script accepts .ab1, .abi, .ab1.gz and .abi.gz inputs. Compressed files are
written as uncompressed .ab1 files in the output directory. Unsupported input
files are ignored. Unsafe isolate IDs are reported; use --sanitise-ids to apply
a consistent safe ID to both the manifest and the output filename.
"""

import argparse
import csv
import os
import re
import shutil
import sys
import unicodedata
import contextlib
from pathlib import Path

from branchmanager.utils.io import write_tsv, sha256sum, sha256sum_gzip


AB1_SUFFIXES = ('.ab1', '.abi', '.ab1.gz', '.abi.gz')
VALID_ID_RE = re.compile(r'^[A-Za-z0-9][A-Za-z0-9_-]*$')
VALID_COMPONENT_RE = re.compile(r'^[A-Za-z0-9_-]+$')


def _cell_text(value):
    """Convert an Excel cell value to clean text without turning 12.0 into '12.0'."""
    if value is None:
        return ''
    try:
        import pandas as pd
        if pd.isna(value):
            return ''
    except (ImportError, TypeError, ValueError):
        pass
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value).strip()


def _strip_ab1_suffix(value):
    name = Path(str(value)).name
    lower = name.lower()
    for suffix in sorted(AB1_SUFFIXES, key=len, reverse=True):
        if lower.endswith(suffix):
            return name[:-len(suffix)]
    return Path(name).stem


def _has_ab1_suffix(path):
    lower = str(path).lower()
    return any(lower.endswith(suffix) for suffix in AB1_SUFFIXES)


def _squash(key):
    """Make lookup tolerant of case and a stray space after a dot."""
    return re.sub(r'\.\s+', '.', str(key).strip()).casefold()


def sanitise_component(value):
    """Transliterate and replace disallowed characters with underscores.

    The partner guide allows only ASCII letters, digits, underscores and
    hyphens in the isolate ID and primer part of a filename. Periods are used
    only for the final .ab1 extension.
    """
    text = unicodedata.normalize('NFKD', str(value))
    text = text.encode('ascii', 'ignore').decode('ascii')
    text = re.sub(r'[^A-Za-z0-9_-]+', '_', text)
    text = re.sub(r'_+', '_', text)
    text = re.sub(r'-+', '-', text)
    text = text.strip('_-')
    return text


def _build_source_index(source_dir):
    """Index supported ABI trace files by their possible provided names."""
    index = {}
    for child in source_dir.rglob('*'):
        if not child.is_file() or not _has_ab1_suffix(child.name):
            continue
        stem = _strip_ab1_suffix(child)
        keys = {
            child.name.casefold(),
            child.stem.casefold(),
            stem.casefold(),
            Path(stem).name.casefold(),
        }
        # Also index every relative-path suffix. Partner exports sometimes add
        # a date/person prefix to the top-level submission folder while the
        # workbook retains only the stable submission identifier.
        try:
            relative_parts = child.resolve().relative_to(source_dir.resolve()).parts
            for offset in range(len(relative_parts)):
                suffix_path = Path(*relative_parts[offset:]).as_posix()
                keys.add(suffix_path.casefold())
                if offset == len(relative_parts) - 1:
                    continue
                suffix_stem = str(Path(suffix_path).with_suffix(''))
                keys.add(suffix_stem.casefold())
        except ValueError:
            pass
        for key in keys | {_squash(k) for k in keys}:
            if key:
                index.setdefault(key, []).append(child)
    return index


def _resolve_source_file(source_dir, value, source_index):
    """Resolve a workbook reference without crossing submission boundaries.

    A directory-qualified workbook reference is authoritative. If the export
    has prefixed the actual submission directory with date/partner text, allow
    that directory name to end with the supplied directory component(s), but
    never fall back to a matching basename in an unrelated directory.
    Basename-only lookup is reserved for references that contain no directory.
    """
    supplied_text = str(value).strip().replace('\\\\', '/')
    supplied = Path(supplied_text).expanduser()
    candidates = [supplied] if supplied.is_absolute() else [source_dir / supplied, supplied]

    # First honour the literal path, including AB1 extension variants.
    for candidate in candidates:
        if candidate.is_file() and _has_ab1_suffix(candidate.name):
            return candidate.resolve()
        for suffix in AB1_SUFFIXES:
            expanded = Path(str(candidate) + suffix)
            if expanded.is_file():
                return expanded.resolve()

    # A qualified relative path may correspond to an exported folder whose
    # name has a date/partner prefix, e.g. "2026-11-12 Partner EC00170083".
    has_directory = supplied.parent != Path('.')
    if has_directory and not supplied.is_absolute():
        requested_parts = [part.casefold() for part in supplied.parent.parts if part not in {'', '.'}]
        requested_stem = _strip_ab1_suffix(supplied.name).casefold()
        qualified_matches = []
        for path in source_dir.rglob('*'):
            if not path.is_file() or not _has_ab1_suffix(path.name):
                continue
            if _strip_ab1_suffix(path.name).casefold() != requested_stem:
                continue
            actual_parts = list(path.parent.relative_to(source_dir).parts)
            if len(actual_parts) < len(requested_parts):
                continue
            suffix_parts = actual_parts[-len(requested_parts):]
            # Exact component match, except the first component may have an
            # export-added prefix. This supports partner export folder naming.
            if all(
                actual.casefold() == requested
                or (i == 0 and actual.casefold().endswith(requested))
                for i, (actual, requested) in enumerate(zip(suffix_parts, requested_parts))
            ):
                qualified_matches.append(path.resolve())
        qualified_matches = sorted(set(qualified_matches))
        if len(qualified_matches) == 1:
            return qualified_matches[0]
        if len(qualified_matches) > 1:
            print(f"Warning: multiple submission-qualified files match '{value}':")
            for match in qualified_matches:
                print(f'  - {match}')
        return None

    # Only an unqualified basename may use global basename matching.
    stem = _strip_ab1_suffix(supplied_text)
    keys = {supplied_text.casefold(), supplied.name.casefold(),
            stem.casefold(), Path(stem).name.casefold()}
    keys.update(f'{stem}{suffix}'.casefold() for suffix in AB1_SUFFIXES)
    keys |= {_squash(key) for key in keys}
    matches, seen = [], set()
    for key in keys:
        for path in source_index.get(key, []):
            resolved = path.resolve()
            if resolved not in seen:
                matches.append(resolved)
                seen.add(resolved)
    if len(matches) > 1:
        print(f"Warning: multiple source files match unqualified reference '{value}':")
        for match in sorted(matches):
            print(f'  - {match}')
        return None
    return matches[0] if matches else None


def _unique_in_order(values):
    result = []
    seen = set()
    for value in values:
        if value not in seen:
            result.append(value)
            seen.add(value)
    return result


def _inspect_ids(ids):
    """Return ID issues using the naming guide's rules."""
    issues = []
    ids = _unique_in_order([_cell_text(value) for value in ids if _cell_text(value)])

    by_casefold = {}
    for raw_id in ids:
        if not VALID_ID_RE.fullmatch(raw_id):
            issues.append({
                'original_id': raw_id,
                'suggested_id': sanitise_component(raw_id),
                'issue': 'Contains disallowed characters or does not start with a letter/number.',
            })
        folded = raw_id.casefold()
        by_casefold.setdefault(folded, []).append(raw_id)

    for variants in by_casefold.values():
        if len(variants) > 1:
            for raw_id in variants:
                issues.append({
                    'original_id': raw_id,
                    'suggested_id': sanitise_component(raw_id),
                    'issue': 'ID differs from another ID only by capitalisation; IDs must be unique case-insensitively.',
                })

    # Also detect different originals that would become the same after cleaning.
    by_cleaned = {}
    for raw_id in ids:
        cleaned = sanitise_component(raw_id)
        if cleaned:
            by_cleaned.setdefault(cleaned.casefold(), []).append(raw_id)
    for variants in by_cleaned.values():
        distinct_casefolded = {v.casefold() for v in variants}
        if len(distinct_casefolded) > 1:
            for raw_id in variants:
                issues.append({
                    'original_id': raw_id,
                    'suggested_id': sanitise_component(raw_id),
                    'issue': 'Sanitising this ID would collide with another distinct ID; review the IDs carefully.',
                })

    # De-duplicate issue rows.
    deduped = []
    seen_issue_rows = set()
    for issue in issues:
        key = (issue['original_id'], issue['suggested_id'], issue['issue'])
        if key not in seen_issue_rows:
            deduped.append(issue)
            seen_issue_rows.add(key)
    return deduped


def _build_safe_id_map(ids):
    """Build a stable, safe ID map; preserve case-only variants as one ID."""
    id_map = {}
    casefold_map = {}
    used = set()

    for raw_id in _unique_in_order([_cell_text(value) for value in ids if _cell_text(value)]):
        folded_raw = raw_id.casefold()
        if folded_raw in casefold_map:
            id_map[raw_id] = casefold_map[folded_raw]
            continue

        candidate = sanitise_component(raw_id)
        if not candidate:
            candidate = 'Isolate'

        base = candidate
        counter = 2
        while candidate.casefold() in used:
            candidate = f'{base}_{counter}'
            counter += 1

        candidate = sanitise_component(candidate)
        casefold_map[folded_raw] = candidate
        id_map[raw_id] = candidate
        used.add(candidate.casefold())

    return id_map


def _write_tsv(path, rows, fieldnames):
    write_tsv(path, list(fieldnames), rows)


def _write_id_map(path, id_map):
    rows = [{'original_id': original, 'sanitised_id': sanitised}
            for original, sanitised in id_map.items() if original != sanitised]
    _write_tsv(path, rows, ['original_id', 'sanitised_id'])


def _sha256_file(path):
    return sha256sum(path)


def _sha256_trace_source(path):
    return sha256sum_gzip(path)


def _copy_as_ab1(source, destination):
    """Copy a trace to an uncompressed .ab1 destination without overwriting data silently."""
    source_hash = _sha256_trace_source(source)
    if destination.exists():
        if _sha256_file(destination) == source_hash:
            print(f'Already exists with identical contents; reusing: {destination.name}')
            return
        raise FileExistsError(
            f"Output file '{destination}' already exists with different contents. "
            'Move it away or choose a new output directory; it was not overwritten.'
        )

    temp_path = destination.with_name(destination.name + '.tmp')
    opener = gzip.open if str(source).lower().endswith('.gz') else open
    try:
        with opener(source, 'rb') as source_handle, open(temp_path, 'wb') as target_handle:
            shutil.copyfileobj(source_handle, target_handle)
        shutil.copystat(source, temp_path)
        # Avoid replacing a file that appeared while copying.
        if destination.exists():
            if _sha256_file(destination) == source_hash:
                temp_path.unlink(missing_ok=True)
                return
            raise FileExistsError(f"Output file appeared during copying: '{destination}'")
        os.replace(temp_path, destination)
    finally:
        if temp_path.exists():
            temp_path.unlink()


def _extract_primer_identifier(column_name):
    column_name = str(column_name).strip()
    if '(' in column_name and ')' in column_name and column_name.rfind(')') > column_name.rfind('('):
        return column_name[column_name.rfind('(') + 1:column_name.rfind(')')].strip()
    return column_name


def _primer_direction(primer_identifier, sanitised_primer):
    """Infer direction for common numeric Sanger primers; retain other primer labels."""
    match = re.fullmatch(r'\d+([FR])', str(primer_identifier).strip().upper())
    if match and match.group(1) == 'F':
        return 'forward'
    if match and match.group(1) == 'R':
        return 'reverse'
    return sanitised_primer


def _infer_primer_from_filename(path: Path) -> str:
    """Infer common bacterial 16S and other marker primers from an AB1 filename.

    Numeric primer labels accept optional separators (e.g. 907R, 907-R,
    1492_R); ITS/18S/28S/DS labels are also recognised. The returned label is
    canonicalised (907R, ITS1, etc.) for consistent grouping and filenames.
    """
    stem = _strip_ab1_suffix(path.name)
    patterns = (
        r'(?<![A-Za-z0-9])(\d{2,4})[_. -]*([FR])(?![A-Za-z0-9])',
        r'(?<![A-Za-z0-9])(ITS|DS)[_. -]*([12])(?![A-Za-z0-9])',
        r'(?<![A-Za-z0-9])(18S|28S)(?![A-Za-z0-9])',
    )
    matches = []
    for pattern in patterns:
        for match in re.finditer(pattern, stem, flags=re.IGNORECASE):
            if len(match.groups()) == 2:
                matches.append((match.start(), ''.join(match.groups()).upper()))
            else:
                matches.append((match.start(), match.group(1).upper()))
    return max(matches, default=(0, ''), key=lambda item: item[0])[1]


def _attempt_group(source_file: Path, source_dir: Path) -> str:
    """Use the source submission folder as provenance, not inferred chronology."""
    try:
        parent = source_file.resolve().relative_to(source_dir.resolve()).parent.as_posix()
    except ValueError:
        parent = source_file.parent.name
    return parent if parent not in {'', '.'} else 'source_root'


def _process_files(input_excel, input_dir, output_dir, manifest_name, sanitise_ids=False):
    try:
        import pandas as pd
    except ImportError:
        print('Error: this tool requires pandas and openpyxl (pip install pandas openpyxl).')
        return 1

    input_path = Path(input_excel)
    source_dir = Path(input_dir)
    dest_dir = Path(output_dir)
    manifest_path = dest_dir / manifest_name

    if not input_path.is_file():
        print(f"Error: Excel file '{input_path}' does not exist.")
        return 1
    if not source_dir.is_dir():
        print(f"Error: Source directory '{source_dir}' does not exist or is not a directory.")
        return 1

    dest_dir.mkdir(parents=True, exist_ok=True)
    print(f'Reading Excel: {input_path}')

    try:
        df = pd.read_excel(input_path, sheet_name='Primer file mapping')
        if len(df.columns) and 'Ensure that the ID' in str(df.columns[0]):
            df = pd.read_excel(input_path, sheet_name='Primer file mapping', skiprows=1)

        df.columns = [re.sub(r'\s+', ' ', str(column)).strip() for column in df.columns]
        if not df.empty and 'Ensure that the ID' in str(df.iloc[0, 0]):
            df = df.iloc[1:].reset_index(drop=True)

        id_col = 'Partner culture collection isolate ID'
        if id_col not in df.columns:
            print(f"Error: Could not find column '{id_col}' in the sheet.")
            print(f'Found columns: {df.columns.tolist()}')
            return 1

        source_index = _build_source_index(source_dir)
        records = []
        source_isolates = {}
        all_ids = []

        print('Reading mapping rows and resolving source files...')
        for _, row in df.iterrows():
            raw_sample_id = _cell_text(row[id_col])
            if not raw_sample_id or 'Ensure that the ID' in raw_sample_id:
                continue
            all_ids.append(raw_sample_id)

            for column in df.columns:
                if column == id_col:
                    continue
                filename = _cell_text(row[column])
                if not filename or 'Please, leave blank' in filename:
                    continue

                source_file = _resolve_source_file(source_dir, filename, source_index)
                if source_file is None:
                    print(f"Warning: File not found or ambiguous: {filename} (Sample: {raw_sample_id})")
                    continue
                if not _has_ab1_suffix(source_file.name):
                    print(f"Warning: Not a supported ABI trace file: {source_file}; row skipped.")
                    continue

                primer_identifier = _extract_primer_identifier(column)
                if primer_identifier.casefold().startswith('other'):
                    inferred = _infer_primer_from_filename(source_file)
                    if inferred:
                        primer_identifier = inferred
                    else:
                        print(f"Warning: Could not infer primer from '{source_file.name}' for {raw_sample_id} ({column}); retaining review label.")
                safe_primer = sanitise_component(primer_identifier)
                if not safe_primer:
                    print(f"Warning: Could not derive a valid primer name from column '{column}'; row skipped.")
                    continue

                record = {
                    'original_id': raw_sample_id,
                    'spreadsheet_filename': filename,
                    'mapping_column': column,
                    'mapping_row': str(int(_) + 2),
                    'source_file': source_file.resolve(),
                    'primer_original': primer_identifier,
                    'primer': safe_primer,
                    'direction': _primer_direction(primer_identifier, safe_primer),
                    'attempt_group': _attempt_group(source_file, source_dir),
                }
                records.append(record)
                source_isolates.setdefault(source_file.resolve(), set()).add(raw_sample_id)

        # A repeated workbook reference to the exact same physical file is a
        # duplicate mapping row, not a new sequencing attempt.
        deduplicated_records = []
        seen_mapping_records = set()
        for record in records:
            key = (
                record['original_id'].casefold(),
                record['primer'].casefold(),
                str(record['source_file'].resolve()).casefold(),
            )
            if key in seen_mapping_records:
                print(
                    f"Warning: duplicate mapping reference for {record['original_id']} "
                    f"/ {record['primer']} -> {record['source_file']}; skipped."
                )
                continue
            seen_mapping_records.add(key)
            deduplicated_records.append(record)
        records = deduplicated_records

        if not records:
            print('Error: No usable trace files were matched to the primer-mapping spreadsheet.')
            return 1

        raw_ids = _unique_in_order(all_ids)
        id_issues = _inspect_ids(raw_ids)
        report_path = dest_dir / 'id_sanitation_report.tsv'
        if id_issues:
            _write_tsv(report_path, id_issues, ['original_id', 'suggested_id', 'issue'])
            print(f'Warning: {len(id_issues)} isolate ID naming issue(s); see {report_path}')
            if not sanitise_ids:
                print('Error: output filenames must match the manifest IDs. Re-run with --sanitise-ids to sanitise IDs and filenames together.')
                return 1
            id_map = _build_safe_id_map(raw_ids)
            _write_id_map(dest_dir / 'id_map.tsv', id_map)
            print(f'Isolate IDs sanitised consistently; original-to-new mapping: {dest_dir / "id_map.tsv"}')
        else:
            id_map = {value: value for value in raw_ids}
            if report_path.exists():
                report_path.unlink()
            stale_map = dest_dir / 'id_map.tsv'
            if stale_map.exists():
                stale_map.unlink()
            if sanitise_ids:
                print('All spreadsheet isolate IDs already follow the naming rules; no id_map.tsv was needed.')
                print('Original-to-output AB1 filename changes will be recorded in ab1_filename_map.tsv.')

        # Repeated isolate/primer combinations receive unique attempt-aware names.
        combo_counts = {}
        for record in records:
            output_id = sanitise_component(id_map.get(record['original_id'], record['original_id']))
            key = (output_id.casefold(), record['primer'].casefold())
            combo_counts[key] = combo_counts.get(key, 0) + 1
        combo_seen, targets, planned_records = {}, {}, []
        for record in records:
            output_id = sanitise_component(id_map.get(record['original_id'], record['original_id']))
            if not output_id or not output_id[0].isalnum():
                print(f"Error: Could not create a valid output ID from '{record['original_id']}'.")
                return 1
            key = (output_id.casefold(), record['primer'].casefold())
            combo_seen[key] = combo_seen.get(key, 0) + 1
            if combo_counts[key] > 1:
                number = combo_seen[key]
                suffix, attempt_id = f'_attempt{number:02d}', f'attempt{number:02d}'
            else:
                suffix, attempt_id = '', 'attempt01'
            target_name = f"{output_id}_{record['primer']}{suffix}.ab1"
            if target_name.casefold() in targets:
                print(f'Error: output filename collision: {target_name}')
                return 1
            record.update({
                'output_id': output_id, 'attempt_id': attempt_id,
                'target_name': target_name,
                'source_sha256': _sha256_trace_source(record['source_file']),
            })
            targets[target_name.casefold()] = record
            planned_records.append(record)

        # Warn when the same physical source trace is assigned to multiple IDs.
        for source_file, isolate_ids in source_isolates.items():
            if len({value.casefold() for value in isolate_ids}) > 1:
                print(
                    f"Warning: AB1 file '{source_file}' is assigned to multiple isolates: "
                    f"{', '.join(sorted(isolate_ids))}"
                )

        # Preflight all existing destination names before copying anything. This prevents a
        # late collision from leaving a partly updated output batch.
        for record in planned_records:
            destination = dest_dir / record['target_name']
            if destination.exists():
                source_hash = _sha256_trace_source(record['source_file'])
                if _sha256_file(destination) != source_hash:
                    print(
                        f"Error: output file '{destination}' already exists with different contents. "
                        'No files have been copied. Move the existing file or choose a new output directory.'
                    )
                    return 1

        # Copy traces (decompressing .gz inputs) under their canonical output names.
        manifest_data = []
        try:
            for record in planned_records:
                destination = dest_dir / record['target_name']
                _copy_as_ab1(record['source_file'], destination)
                manifest_data.append([
                    record['output_id'], record['target_name'], record['direction'],
                    record['primer'], record['attempt_id'], record['attempt_group'],
                    str(record['source_file']), record['source_file'].name,
                    record['source_sha256'], record['spreadsheet_filename'],
                    record['mapping_column'], record['mapping_row'],
                ])
        except (OSError, FileExistsError) as error:
            print(f'Error while copying trace files: {error}')
            return 1

        with open(manifest_path, mode='w', newline='', encoding='utf-8') as manifest_handle:
            writer = csv.writer(manifest_handle)
            writer.writerow([
                'sequence_id', 'read_file', 'direction', 'primer', 'attempt_id',
                'attempt_group', 'source_path', 'source_file', 'source_sha256',
                'spreadsheet_filename', 'mapping_column', 'mapping_row',
            ])
            writer.writerows(manifest_data)

        # Record every original source filename -> canonical output filename mapping.
        # This is separate from id_map.tsv: source filenames may contain legacy names
        # (e.g. 'JJ A1.1492r.ab1') while the spreadsheet links that trace to a valid
        # isolate ID (e.g. 'BSI_0001').
        filename_map_path = dest_dir / 'ab1_filename_map.tsv'
        filename_map_rows = []
        for record in planned_records:
            filename_map_rows.append({
                'source_path': str(record['source_file']),
                'source_filename': record['source_file'].name,
                'spreadsheet_filename': record['spreadsheet_filename'],
                'original_isolate_id': record['original_id'],
                'output_isolate_id': record['output_id'],
                'original_primer': record['primer_original'],
                'output_primer': record['primer'],
                'direction': record['direction'],
                'output_filename': record['target_name'],
                'attempt_id': record['attempt_id'],
                'attempt_group': record['attempt_group'],
                'source_sha256': record['source_sha256'],
                'mapping_column': record['mapping_column'],
                'mapping_row': record['mapping_row'],
            })
        _write_tsv(
            filename_map_path,
            filename_map_rows,
            ['source_path', 'source_filename', 'spreadsheet_filename',
             'original_isolate_id', 'output_isolate_id', 'original_primer',
             'output_primer', 'direction', 'output_filename', 'attempt_id',
             'attempt_group', 'source_sha256', 'mapping_column', 'mapping_row'],
        )

        print(f'Success! Sanitised .ab1 files copied to: {dest_dir}')
        print(f'Manifest created at: {manifest_path}')
        print(f'Original-to-output filename mapping created at: {filename_map_path}')
        print('Output filenames use ISOLATE_ID_PRIMER.ab1 for unique reads and add _attemptNN for repeated isolate/primer reads.')
        return 0

    except Exception as error:
        print(f'An error occurred: {error}')
        return 1



class _Tee:
    """Write console output to both the original stream and a UTF-8 log."""
    def __init__(self, console_stream, log_handle):
        self.console_stream = console_stream
        self.log_handle = log_handle

    def write(self, text):
        self.console_stream.write(text)
        self.log_handle.write(text)
        self.log_handle.flush()
        return len(text)

    def flush(self):
        self.console_stream.flush()
        self.log_handle.flush()

    def isatty(self):
        return getattr(self.console_stream, 'isatty', lambda: False)()


def process_files(input_excel, input_dir, output_dir, manifest_name,
                  sanitise_ids=False, log_file=None):
    """Run the converter while recording all stdout/stderr in a log file."""
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    log_path = Path(log_file) if log_file else destination / 'primers_to_manifest.log'
    if not log_path.is_absolute():
        # Relative --log paths are interpreted relative to the current directory.
        log_path = log_path.resolve()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open('w', encoding='utf-8') as handle:
        stdout_tee = _Tee(sys.stdout, handle)
        stderr_tee = _Tee(sys.stderr, handle)
        with contextlib.redirect_stdout(stdout_tee), contextlib.redirect_stderr(stderr_tee):
            print(f'Run started: {__import__("datetime").datetime.now().astimezone().isoformat()}')
            print(f'Log file: {log_path}')
            result = _process_files(input_excel, input_dir, output_dir, manifest_name,
                                    sanitise_ids=sanitise_ids)
            print(f'Run finished with exit status: {result}')
    return result

def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument('-i', '--input', required=True, help='Path to the Excel (.xlsx) file')
    parser.add_argument('-s', '--source', required=True, help='Directory containing source .ab1/.abi files (may include .gz files)')
    parser.add_argument('-o', '--output', required=True, help='Directory to copy and standardise trace filenames into')
    parser.add_argument('-m', '--manifest', default='manifest.csv', help='Name of output manifest (default: manifest.csv)')
    parser.add_argument('--log', default=None, help='Log file path (default: OUTPUT/primers_to_manifest.log)')
    parser.add_argument(
        '--sanitise-ids', '--sanitize-ids', dest='sanitise_ids', action='store_true',
        help='Sanitise unsafe isolate IDs as well as output filenames; write id_map.tsv and id_sanitation_report.tsv. '
             'Required if the input spreadsheet contains unsafe or colliding IDs.'
    )


def run_cli(args: argparse.Namespace):
    return process_files(args.input, args.source, args.output, args.manifest, sanitise_ids=args.sanitise_ids, log_file=args.log)


def main():
    parser = argparse.ArgumentParser(
        description="Copy .ab1 files from a primer mapping workbook, preserving alternative/resequenced reads and source provenance."
    )
    add_arguments(parser)
    args = parser.parse_args()
    raise SystemExit(run_cli(args) or 0)


if __name__ == '__main__':
    main()