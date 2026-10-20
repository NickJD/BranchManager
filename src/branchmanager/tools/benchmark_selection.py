#!/usr/bin/env python3
"""Evaluate isolate-ranking outputs against independently adjudicated labels.

Input columns: case_id, method, isolate_id, rank, relevant. `relevant` is 1
when the isolate met the predeclared external outcome (for example, selected by
an independent expert panel and successfully sequenced), otherwise 0.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
from collections import defaultdict
from pathlib import Path


def read_rows(path: str | Path) -> list[dict]:
    with open(path, newline='', encoding='utf-8-sig') as handle:
        sample = handle.read(4096)
        handle.seek(0)
        delimiter = '\t' if '\t' in sample and ',' not in sample else ','
        rows = list(csv.DictReader(handle, delimiter=delimiter))
    required = {'case_id', 'method', 'isolate_id', 'rank', 'relevant'}
    missing = required - set(rows[0]) if rows else required
    if missing:
        raise ValueError(f"Missing required column(s): {', '.join(sorted(missing))}")
    return rows


def _truth(value: str) -> int:
    normal = str(value).strip().lower()
    if normal in {'1', 'true', 'yes', 'y'}:
        return 1
    if normal in {'0', 'false', 'no', 'n'}:
        return 0
    raise ValueError(f"relevant must be 0/1 or yes/no, received {value!r}")


def score_case(rows: list[dict], k: int) -> dict:
    rows = sorted(rows, key=lambda row: (int(row['rank']), row['isolate_id']))
    labels = [_truth(row['relevant']) for row in rows]
    total_positive = sum(labels)
    top_labels = labels[:k]
    hits = sum(top_labels)
    precision = hits / len(top_labels) if top_labels else 0.0
    recall = hits / total_positive if total_positive else math.nan
    seen = 0
    precision_sum = 0.0
    for rank, label in enumerate(labels, start=1):
        if label:
            seen += 1
            precision_sum += seen / rank
    average_precision = precision_sum / total_positive if total_positive else math.nan
    return {'precision_at_k': precision, 'recall_at_k': recall, 'average_precision': average_precision,
            'positives': total_positive, 'candidates': len(rows), 'selected': len(top_labels)}


def run(input_path: str, output_path: str, k_values: list[int], bootstrap: int = 2000, seed: int = 1) -> None:
    grouped = defaultdict(list)
    seen = set()
    for row in read_rows(input_path):
        try:
            rank = int(row['rank'])
            if rank < 1:
                raise ValueError
            _truth(row['relevant'])
        except ValueError as exc:
            raise ValueError(f"Invalid rank or relevance value in row {row!r}") from exc
        key = (row['case_id'], row['method'], row['isolate_id'])
        if key in seen:
            raise ValueError(f"Duplicate case/method/isolate row: {key}")
        seen.add(key)
        grouped[(row['case_id'], row['method'])].append(row)
    summaries = defaultdict(lambda: defaultdict(list))
    for (case_id, method), rows in grouped.items():
        for k in k_values:
            values = score_case(rows, k)
            for metric in ('precision_at_k', 'recall_at_k', 'average_precision'):
                if not math.isnan(values[metric]):
                    summaries[(method, k)][metric].append(values[metric])
            summaries[(method, k)]['case_count'].append(case_id)
    fields = ['method', 'k', 'cases', 'mean_precision_at_k', 'precision_ci_low', 'precision_ci_high',
              'mean_recall_at_k', 'recall_ci_low', 'recall_ci_high',
              'mean_average_precision', 'map_ci_low', 'map_ci_high']
    with open(output_path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for method, k in sorted(summaries):
            metrics = summaries[(method, k)]
            mean = lambda name: sum(metrics[name]) / len(metrics[name]) if metrics[name] else ''
            output = {'method': method, 'k': k, 'cases': len(metrics['case_count']),
                      'mean_precision_at_k': mean('precision_at_k'),
                      'mean_recall_at_k': mean('recall_at_k'),
                      'mean_average_precision': mean('average_precision')}
            rng = random.Random(f'{seed}:{method}:{k}')
            for metric, low_key, high_key in (
                ('precision_at_k', 'precision_ci_low', 'precision_ci_high'),
                ('recall_at_k', 'recall_ci_low', 'recall_ci_high'),
                ('average_precision', 'map_ci_low', 'map_ci_high'),
            ):
                case_values = metrics[metric]
                if bootstrap > 0 and case_values:
                    estimates = sorted(
                        sum(rng.choice(case_values) for _ in case_values) / len(case_values)
                        for _ in range(bootstrap)
                    )
                    output[low_key] = estimates[int(0.025 * (bootstrap - 1))]
                    output[high_key] = estimates[int(0.975 * (bootstrap - 1))]
                else:
                    output[low_key] = output[high_key] = ''
            writer.writerow(output)


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument('--input', required=True, help='TSV/CSV ranking and adjudication table.')
    parser.add_argument('--out', required=True, help='Output summary TSV.')
    parser.add_argument('--k', type=int, nargs='+', default=[1, 5, 10], help='Selection cutoffs to evaluate.')
    parser.add_argument('--bootstrap', type=int, default=2000, help='Case-level bootstrap replicates for 95%% percentile intervals (0 disables).')
    parser.add_argument('--seed', type=int, default=1, help='Bootstrap random seed.')


def run_cli(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    if not args.k or min(args.k) < 1:
        parser.error('all --k values must be positive')
    if args.bootstrap < 0:
        parser.error('--bootstrap must be zero or greater')
    run(args.input, args.out, sorted(set(args.k)), bootstrap=args.bootstrap, seed=args.seed)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    run_cli(parser.parse_args(), parser)


if __name__ == '__main__':
    main()
