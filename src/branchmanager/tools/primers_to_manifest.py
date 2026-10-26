import re
import shutil
import argparse
import csv
from pathlib import Path

from branchmanager.tools.id_sanitation import build_id_map, find_id_issues, write_id_map, write_issue_report


AB1_SUFFIXES = ('.ab1', '.abi', '.ab1.gz', '.abi.gz')


def _strip_ab1_suffix(value):
    name = Path(str(value)).name
    lower = name.lower()
    for suffix in sorted(AB1_SUFFIXES, key=len, reverse=True):
        if lower.endswith(suffix):
            return name[:-len(suffix)]
    return Path(name).stem


def _squash(key):
    """Casefold and drop spaces after dots ('JJ A1. 27f' == 'JJ A1.27f')."""
    return re.sub(r'\.\s+', '.', str(key).strip()).casefold()


def _build_source_index(source_dir):
    index = {}
    for child in source_dir.rglob('*'):
        if not child.is_file():
            continue
        stem = _strip_ab1_suffix(child)
        keys = {
            child.name.casefold(),
            child.stem.casefold(),
            stem.casefold(),
            Path(stem).name.casefold(),
        }
        for key in keys | {_squash(k) for k in keys}:
            if not key:
                continue
            index.setdefault(key, []).append(child)
    return index


def _resolve_source_file(source_dir, value, source_index):
    supplied = Path(str(value)).expanduser()
    if supplied.is_absolute() and supplied.is_file():
        return supplied

    direct = source_dir / supplied
    if direct.is_file():
        return direct

    stem = _strip_ab1_suffix(value)
    candidates = {
        str(value).strip().casefold(),
        Path(str(value)).name.casefold(),
        supplied.name.casefold(),
        stem.casefold(),
        Path(stem).name.casefold(),
    }

    for suffix in AB1_SUFFIXES:
        candidates.add(f"{stem}{suffix}".casefold())

    candidates |= {_squash(c) for c in candidates}
    matches = []
    seen = set()
    for key in candidates:
        for path in source_index.get(key, []):
            resolved = path.resolve()
            if resolved in seen:
                continue
            matches.append(resolved)
            seen.add(resolved)
    return sorted(matches)[0] if matches else None


def process_files(input_excel, input_dir, output_dir, manifest_name, sanitise_ids=False):
    try:
        import pandas as pd
    except ImportError:
        print("Error: this tool requires pandas and openpyxl (pip install pandas openpyxl).")
        return 1
    # Convert strings to Path objects
    input_path = Path(input_excel)
    source_dir = Path(input_dir)
    dest_dir = Path(output_dir)
    manifest_path = dest_dir / manifest_name

    # Create output directory if it doesn't exist
    dest_dir.mkdir(parents=True, exist_ok=True)

    if not source_dir.exists():
        print(f"Error: Source directory '{source_dir}' does not exist.")
        return

    print(f"Reading Excel: {input_path}")

    try:
        # Read the specific sheet. pandas treats Excel row 1 as the header by
        # default, so an instruction in A1 appears in df.columns, not df.iloc[0].
        df = pd.read_excel(input_path, sheet_name='Primer file mapping')

        if len(df.columns) and "Ensure that the ID" in str(df.columns[0]):
            df = pd.read_excel(input_path, sheet_name='Primer file mapping', skiprows=1)

        # Clean up column names (remove newlines and extra whitespace)
        df.columns = [str(c).replace('\n', ' ').replace('  ', ' ').strip() for c in df.columns]

        # Also handle instruction text in the first data row.
        if not df.empty and "Ensure that the ID" in str(df.iloc[0, 0]):
            df = df.iloc[1:].reset_index(drop=True)

        id_col = 'Partner culture collection isolate ID'
        if id_col not in df.columns:
            print(f"Error: Could not find column '{id_col}' in the sheet.")
            print(f"Found columns: {df.columns.tolist()}")
            return

        source_index = _build_source_index(source_dir)
        manifest_data = []
        source_isolates = {}

        print("Processing rows...")
        for _, row in df.iterrows():
            sample_id = row[id_col]

            # --- SKIP INSTRUCTION ROWS ---
            # Skip rows where the ID column contains the "Ensure that the ID..." instruction
            if pd.isna(sample_id) or "Ensure that the ID" in str(sample_id):
                continue

            # Iterate through all columns to find primer file entries
            for col in df.columns:
                if col == id_col:
                    continue

                filename = row[col]

                # --- SKIP EMPTY OR INSTRUCTION CELLS ---
                # Skip if cell is empty, NaN, or contains the "Please, leave blank" instruction
                if pd.isna(filename) or "Please, leave blank" in str(filename) or str(filename).strip() == "":
                    continue

                filename = str(filename).strip()

                # Extract the primer identifier (e.g., "27F") from the header.
                if '(' in col and ')' in col:
                    primer_identifier = col.split('(')[-1].replace(')', '').strip()
                else:
                    primer_identifier = col.strip()

                primer_name = primer_identifier

                # Determine direction based ONLY on the extracted identifier.
                if 'F' in primer_identifier.upper() and 'R' not in primer_identifier.upper():
                    direction = "forward"
                elif 'R' in primer_identifier.upper():
                    direction = "reverse"
                else:
                    # Fallback for non-standard primers like ITS1, DS1, 18S, etc.
                    direction = primer_identifier


                src_file = _resolve_source_file(source_dir, filename, source_index)

                if src_file:
                    dest_file = dest_dir / src_file.name
                    shutil.copy2(src_file, dest_file)

                    source_isolates.setdefault(src_file.resolve(), set()).add(str(sample_id).strip())

                    # Append to manifest list: [ID, filename, direction, clean_primer_name]
                    manifest_data.append([str(sample_id).strip(), src_file.name, direction, primer_name])
                else:
                    print(f"Warning: File not found: {filename} (Sample: {sample_id})")

        for src_file, isolates in source_isolates.items():
            if len(isolates) > 1:
                print(
                    f"Warning: AB1 file '{src_file}' is assigned to multiple isolates: "
                    f"{', '.join(sorted(isolates))}"
                )

        # ID sanitation: always report problems; rewrite IDs only with --sanitise-ids
        ids = [row[0] for row in manifest_data]
        issues = find_id_issues(ids)
        report_path = dest_dir / 'id_sanitation_report.tsv'
        if issues:
            write_issue_report(issues, report_path)
            print(f"Warning: {len({i['id'] for i in issues})} unsafe/colliding sequence ID(s); see {report_path}")
            if sanitise_ids:
                id_map = build_id_map(ids)
                write_id_map(id_map, dest_dir / 'id_map.tsv')
                for row in manifest_data:
                    row[0] = id_map[row[0]]
                print(f"Sequence IDs sanitised; original->new map: {dest_dir / 'id_map.tsv'}")
            else:
                print("Re-run with --sanitise-ids to rewrite them (IDs with spaces break vsearch/Newick).")
        elif report_path.exists():
            report_path.unlink()

        # Write the manifest file as a true CSV (comma delimited)
        with open(manifest_path, mode='w', newline='', encoding='utf-8') as m_file:
            writer = csv.writer(m_file)
            writer.writerow(['sequence_id', 'read_file', 'direction', 'primer'])
            writer.writerows(manifest_data)

        print(f"Success! Files copied to: {dest_dir}")
        print(f"Manifest created at: {manifest_path}")

    except Exception as e:
        print(f"An error occurred: {e}")


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("-i", "--input", required=True, help="Path to the Excel (.xlsx) file")
    parser.add_argument("-s", "--source", required=True, help="Directory containing the source .ab1 files")
    parser.add_argument("-o", "--output", required=True, help="Directory to copy processed files to")
    parser.add_argument("-m", "--manifest", default="manifest.csv",
                        help="Name of the output manifest (default: manifest.csv)")
    parser.add_argument("--sanitise-ids", "--sanitize-ids", dest="sanitise_ids", action="store_true",
                        help="Rewrite unsafe sequence IDs (spaces, punctuation, non-ASCII) and write id_map.tsv. "
                             "Without this flag unsafe IDs are only reported.")


def run_cli(args: argparse.Namespace):
    return process_files(args.input, args.source, args.output, args.manifest, sanitise_ids=args.sanitise_ids)


def main():
    parser = argparse.ArgumentParser(description="Process .ab1 files based on 'Primer file mapping' Excel sheet.")
    add_arguments(parser)
    run_cli(parser.parse_args())


if __name__ == "__main__":
    main()
