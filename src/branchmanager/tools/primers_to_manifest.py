#!/usr/bin/env python3
"""Copy Sanger trace files from a primer-mapping Excel sheet and create a manifest.

Output trace names are standardised to:
    ISOLATE_ID_PRIMER.ab1

The script accepts .ab1, .abi, .ab1.gz and .abi.gz inputs. Compressed files are
written as uncompressed .ab1 files in the output directory. Unsupported input
files are ignored. Unsafe isolate IDs are reported; use --sanitise-ids to apply
a consistent safe ID to both the manifest and the output filename.
"""

import argparse
import csv
import gzip
import hashlib
import os
import re
import shutil
import sys
import unicodedata
from pathlib import Path


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
        for key in keys | {_squash(k) for k in keys}:
            if key:
                index.setdefault(key, []).append(child)
    return index


def _resolve_source_file(source_dir, value, source_index):
    supplied_text = str(value).strip()
    supplied = Path(supplied_text).expanduser()
    if supplied.is_absolute() and supplied.is_file() and _has_ab1_suffix(supplied.name):
        return supplied

    direct = source_dir / supplied
    if direct.is_file() and _has_ab1_suffix(direct.name):
        return direct

    stem = _strip_ab1_suffix(supplied_text)
    candidates = {
        supplied_text.casefold(),
        Path(supplied_text).name.casefold(),
        supplied.name.casefold(),
        stem.casefold(),
        Path(stem).name.casefold(),
    }
    for suffix in AB1_SUFFIXES:
        candidates.add(f'{stem}{suffix}'.casefold())
    candidates |= {_squash(candidate) for candidate in candidates}

    matches = []
    seen = set()
    for key in candidates:
        for path in source_index.get(key, []):
            resolved = path.resolve()
            if resolved not in seen:
                matches.append(resolved)
                seen.add(resolved)
    # Do not silently choose between ambiguous matches with the same name.
    if len(matches) > 1:
        print(f"Warning: multiple source files match '{value}':")
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
    with open(path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def _write_id_map(path, id_map):
    rows = [{'original_id': original, 'sanitised_id': sanitised}
            for original, sanitised in id_map.items() if original != sanitised]
    _write_tsv(path, rows, ['original_id', 'sanitised_id'])


def _sha256_file(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _sha256_trace_source(path):
    digest = hashlib.sha256()
    opener = gzip.open if str(path).lower().endswith('.gz') else open
    with opener(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


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


def process_files(input_excel, input_dir, output_dir, manifest_name, sanitise_ids=False):
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

                primer_identifier = _extract_primer_identifier(column)
                safe_primer = sanitise_component(primer_identifier)
                if not safe_primer:
                    print(f"Warning: Could not derive a valid primer name from column '{column}'; row skipped.")
                    continue

                source_file = _resolve_source_file(source_dir, filename, source_index)
                if source_file is None:
                    print(f"Warning: File not found or ambiguous: {filename} (Sample: {raw_sample_id})")
                    continue
                if not _has_ab1_suffix(source_file.name):
                    print(f"Warning: Not a supported ABI trace file: {source_file}; row skipped.")
                    continue

                record = {
                    'original_id': raw_sample_id,
                    # Keep exactly what the mapping spreadsheet supplied, as well as
                    # the actual matched source file. This supports a separate audit
                    # trail for filename normalisation even when the isolate ID itself
                    # did not require sanitisation.
                    'spreadsheet_filename': filename,
                    'source_file': source_file.resolve(),
                    'primer_original': primer_identifier,
                    'primer': safe_primer,
                    'direction': _primer_direction(primer_identifier, safe_primer),
                }
                records.append(record)
                source_isolates.setdefault(source_file.resolve(), set()).add(raw_sample_id)

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

        # Assign canonical target names before copying anything, and detect collisions.
        targets = {}
        planned_records = []
        for record in records:
            output_id = id_map.get(record['original_id'], record['original_id'])
            output_id = sanitise_component(output_id)
            if not output_id or not output_id[0].isalnum():
                print(f"Error: Could not create a valid output ID from '{record['original_id']}'.")
                return 1

            # This follows the attached guide's explicit rules/examples: ID_PRIMER.ab1.
            target_name = f"{output_id}_{record['primer']}.ab1"
            if len(target_name) > 40:
                print(f"Warning: output filename is {len(target_name)} characters (>40 preferred): {target_name}")

            record_key = (target_name.casefold(), record['source_file'])
            prior = targets.get(target_name.casefold())
            if prior is not None:
                if prior['source_file'] == record['source_file'] and prior['output_id'].casefold() == output_id.casefold() and prior['primer'].casefold() == record['primer'].casefold():
                    print(f'Warning: duplicate mapping row for {target_name}; duplicate manifest entry skipped.')
                    continue
                print(f'Error: more than one source trace would be written as {target_name}:')
                print(f"  - {prior['source_file']}")
                print(f"  - {record['source_file']}")
                print('No files have been copied. Correct the primer map or resolve the duplicate naming first.')
                return 1

            record['output_id'] = output_id
            record['target_name'] = target_name
            targets[target_name.casefold()] = record
            planned_records.append(record)

        # Warn when the same physical source trace is assigned to multiple IDs.
        for source_file, isolate_ids in source_isolates.items():
            if len({value.casefold() for value in isolate_ids}) > 1:
                print(
                    f"Warning: source trace '{source_file}' is assigned to multiple isolate IDs: "
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
                    record['output_id'],
                    record['target_name'],
                    record['direction'],
                    record['primer'],
                ])
        except (OSError, FileExistsError) as error:
            print(f'Error while copying trace files: {error}')
            return 1

        with open(manifest_path, mode='w', newline='', encoding='utf-8') as manifest_handle:
            writer = csv.writer(manifest_handle)
            writer.writerow(['sequence_id', 'read_file', 'direction', 'primer'])
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
            })
        _write_tsv(
            filename_map_path,
            filename_map_rows,
            ['source_path', 'source_filename', 'spreadsheet_filename',
             'original_isolate_id', 'output_isolate_id', 'original_primer',
             'output_primer', 'direction', 'output_filename'],
        )

        print(f'Success! Sanitised .ab1 files copied to: {dest_dir}')
        print(f'Manifest created at: {manifest_path}')
        print(f'Original-to-output filename mapping created at: {filename_map_path}')
        print('Output filenames use the form ISOLATE_ID_PRIMER.ab1, e.g. SH_A1_27F.ab1.')
        return 0

    except Exception as error:
        print(f'An error occurred: {error}')
        return 1


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument('-i', '--input', required=True, help='Path to the Excel (.xlsx) file')
    parser.add_argument('-s', '--source', required=True, help='Directory containing source .ab1/.abi files (may include .gz files)')
    parser.add_argument('-o', '--output', required=True, help='Directory to copy and standardise trace filenames into')
    parser.add_argument('-m', '--manifest', default='manifest.csv', help='Name of output manifest (default: manifest.csv)')
    parser.add_argument(
        '--sanitise-ids', '--sanitize-ids', dest='sanitise_ids', action='store_true',
        help='Sanitise unsafe isolate IDs as well as output filenames; write id_map.tsv and id_sanitation_report.tsv. '
             'Required if the input spreadsheet contains unsafe or colliding IDs.'
    )


def run_cli(args: argparse.Namespace):
    return process_files(args.input, args.source, args.output, args.manifest, sanitise_ids=args.sanitise_ids)


def main():
    parser = argparse.ArgumentParser(
        description="Copy .ab1 files based on the 'Primer file mapping' Excel sheet, standardising names to ISOLATE_ID_PRIMER.ab1."
    )
    add_arguments(parser)
    args = parser.parse_args()
    raise SystemExit(run_cli(args) or 0)


if __name__ == '__main__':
    main()