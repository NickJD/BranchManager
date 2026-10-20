"""Isolate/sequence ID sanitation checks.

Downstream tools (vsearch, FastTree, Newick/iTOL, SQLite keys, file names) treat
whitespace and punctuation in sequence IDs inconsistently, so IDs are checked
for safety up front and can optionally be rewritten with a recorded ID map.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import unicodedata
from collections import defaultdict
from pathlib import Path
from typing import Iterable

SAFE_ID = re.compile(r'^[A-Za-z0-9][A-Za-z0-9._-]*$')
ID_COLUMNS = ('sequence_id', 'seq_id', 'isolate_id', 'sample_id', 'id', 'unique_id')


def sanitise_id(value: object) -> str:
    """Return a filesystem/vsearch/Newick-safe version of *value*."""
    text = unicodedata.normalize('NFKD', str(value or '')).encode('ascii', 'ignore').decode('ascii')
    text = re.sub(r'[^A-Za-z0-9._-]+', '_', text.strip())
    return text.strip('_.-')


def find_id_issues(ids: Iterable[object]) -> list[dict[str, str]]:
    """Return one row per problem found: id, code, detail."""
    ids = [str(i) if i is not None else '' for i in ids]
    issues: list[dict[str, str]] = []

    def add(value, code, detail):
        issues.append({'id': value, 'code': code, 'detail': detail})

    for value in ids:
        if not value.strip():
            add(value, 'EMPTY_ID', 'ID is empty')
            continue
        if value != value.strip():
            add(value, 'LEADING_TRAILING_SPACE', 'ID has leading or trailing whitespace')
        if re.search(r'\s', value):
            add(value, 'WHITESPACE', 'ID contains whitespace; vsearch and Newick truncate or mangle it')
        if not value.isascii():
            add(value, 'NON_ASCII', 'ID contains non-ASCII characters')
        if not SAFE_ID.match(value) and not re.search(r'\s', value) and value.isascii():
            add(value, 'UNSAFE_CHARACTERS', 'ID contains characters outside A-Z a-z 0-9 . _ -')

    unique = sorted({i for i in ids if i.strip()})
    by_fold = defaultdict(list)
    by_clean = defaultdict(list)
    for value in unique:
        by_fold[value.casefold()].append(value)
        by_clean[sanitise_id(value).casefold()].append(value)
    for group in by_fold.values():
        if len(group) > 1:
            for value in group:
                add(value, 'CASE_COLLISION', f'differs only by case from: {"; ".join(g for g in group if g != value)}')
    for clean, group in by_clean.items():
        if len(group) > 1 and len({g.casefold() for g in group}) > 1:
            for value in group:
                add(value, 'SANITISED_COLLISION', f'sanitises to {clean!r}, same as: {"; ".join(g for g in group if g != value)}')
    return issues


def build_id_map(ids: Iterable[object]) -> dict[str, str]:
    """Map each distinct original ID to its sanitised form; raise on collisions or empties."""
    mapping: dict[str, str] = {}
    owner: dict[str, str] = {}
    for value in dict.fromkeys(str(i) for i in ids if i is not None):
        clean = sanitise_id(value)
        if not clean:
            raise ValueError(f'ID {value!r} is empty after sanitation')
        key = clean.casefold()
        if key in owner and owner[key] != value:
            raise ValueError(f'IDs {owner[key]!r} and {value!r} both sanitise to {clean!r}')
        owner[key] = value
        mapping[value] = clean
    return mapping


def write_issue_report(issues: list[dict[str, str]], path: str | Path) -> None:
    with open(path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=['id', 'code', 'detail'], delimiter='\t')
        writer.writeheader()
        writer.writerows(issues)


def write_id_map(mapping: dict[str, str], path: str | Path) -> None:
    with open(path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['original_id', 'sanitised_id'])
        for original, clean in mapping.items():
            writer.writerow([original, clean])


def _delimiter(path: Path) -> str:
    return '\t' if path.name.lower().endswith(('.tsv', '.txt', '.tsv.gz')) else ','


def _read_table_ids(path: Path, column: str | None):
    with open(path, newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle, delimiter=_delimiter(path))
        rows = list(reader)
        fields = reader.fieldnames or []
    if column:
        if column not in fields:
            raise SystemExit(f'Column {column!r} not found. Columns: {", ".join(fields)}')
        chosen = column
    else:
        lookup = {f.strip().lower().replace(' ', '_'): f for f in fields}
        chosen = next((lookup[c] for c in ID_COLUMNS if c in lookup), None)
        if chosen is None:
            raise SystemExit(f'No ID column found; pass --column. Columns: {", ".join(fields)}')
    return rows, fields, chosen


def check_file(path: str | Path, column: str | None = None, fix_out: str | Path | None = None,
               report: str | Path | None = None) -> list[dict[str, str]]:
    """Check IDs in a CSV/TSV table or FASTA; optionally write a sanitised copy and id_map."""
    path = Path(path)
    is_fasta = path.name.lower().endswith(('.fasta', '.fa', '.fna', '.fasta.gz', '.fa.gz', '.fna.gz'))
    if is_fasta:
        from branchmanager.utils.fasta import read_fasta, write_fasta
        records = list(read_fasta(path))
        # The whole header is the ID in BranchManager FASTAs.
        ids = [h for h, _ in records]
    else:
        rows, fields, chosen = _read_table_ids(path, column)
        ids = [row.get(chosen, '') for row in rows]

    issues = find_id_issues(ids)
    if report:
        write_issue_report(issues, report)
    if fix_out:
        try:
            mapping = build_id_map(ids)
        except ValueError as exc:
            raise SystemExit(f'[check-ids] cannot sanitise: {exc}')
        write_id_map(mapping, f'{fix_out}.id_map.tsv')
        if is_fasta:
            write_fasta([(mapping[h], s) for h, s in records], fix_out)
        else:
            with open(fix_out, 'w', newline='', encoding='utf-8') as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter=_delimiter(Path(fix_out)))
                writer.writeheader()
                for row in rows:
                    row[chosen] = mapping[row[chosen]]
                    writer.writerow(row)
    return issues


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument('-i', '--input', required=True, help='CSV/TSV table or FASTA whose IDs are checked.')
    parser.add_argument('--column', default=None, help='ID column for tables (default: auto-detect sequence_id/isolate_id/...).')
    parser.add_argument('--report', default=None, help='Write the issue list to this TSV.')
    parser.add_argument('--fix', default=None, metavar='OUT',
                        help='Write a sanitised copy to OUT (and OUT.id_map.tsv with original->sanitised IDs).')


def run(args: argparse.Namespace) -> int:
    issues = check_file(args.input, args.column, args.fix, args.report)
    for issue in issues[:50]:
        print(f"{issue['code']}\t{issue['id']}\t{issue['detail']}")
    if len(issues) > 50:
        print(f'... {len(issues) - 50} more (use --report for the full list)')
    print(f'[check-ids] {"FAIL" if issues else "PASS"}: {len(issues)} issue(s) in {args.input}')
    if args.fix:
        print(f'[check-ids] sanitised copy: {args.fix}')
    return 1 if issues and not args.fix else 0


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    sys.exit(run(parser.parse_args(argv)))


if __name__ == '__main__':
    main()
