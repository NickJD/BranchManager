"""Prepare segmented PacBio HiFi 16S read sets for BranchManager.

SMRT Link is expected to have already segmented/deconcatenated and demultiplexed
the library molecules. Consequently each input FASTQ represents one isolate,
but contains many HiFi observations of its approximately 1.5 kb 16S amplicon.
This module quality filters, dereplicates, clusters those observations, and
selects an abundance-supported observed sequence as the marker representative.
"""
from __future__ import annotations

import csv
import gzip
import shutil
from collections import defaultdict
from pathlib import Path

from branchmanager.utils.fasta import write_fasta
from branchmanager.utils.subprocess import run_cmd


def _open_text(path: Path):
    return gzip.open(path, 'rt') if path.name.lower().endswith('.gz') else open(path)


def _fastq_records(path: Path):
    with _open_text(path) as handle:
        line_number = 0
        while True:
            header = handle.readline()
            if not header:
                return
            sequence, plus, quality = handle.readline(), handle.readline(), handle.readline()
            line_number += 4
            if not quality or not header.startswith('@') or not plus.startswith('+'):
                raise ValueError(f'malformed FASTQ near line {line_number - 3} in {path}')
            sequence = sequence.strip().upper()
            quality = quality.rstrip('\r\n')
            if len(sequence) != len(quality):
                raise ValueError(f'sequence/quality length mismatch near line {line_number - 3} in {path}')
            yield sequence, sum(max(0, ord(base) - 33) for base in quality) / max(1, len(quality))


def _read_map(path: str | Path) -> list[tuple[str, Path]]:
    path = Path(path).resolve()
    delimiter = '\t' if path.suffix.lower() in {'.tsv', '.tab'} else ','
    rows: list[tuple[str, Path]] = []
    seen: set[str] = set()
    with open(path, newline='') as handle:
        for row in csv.DictReader(handle, delimiter=delimiter):
            lowered = {str(key).strip().lower(): str(value or '').strip() for key, value in row.items() if key}
            sequence_id = (lowered.get('sequence_id') or lowered.get('sequenceid') or lowered.get('isolate_id')
                           or lowered.get('isolateid') or lowered.get('sample_id') or lowered.get('sampleid'))
            filename = lowered.get('file') or lowered.get('fastq') or lowered.get('fastq_file') or lowered.get('read_file')
            if not sequence_id or not filename:
                raise ValueError('PacBio map requires sequence_id (or isolate_id) and file (or fastq_file) columns')
            if sequence_id in seen:
                raise ValueError(f'duplicate PacBio sequence_id: {sequence_id}')
            file_path = Path(filename).expanduser()
            if not file_path.is_absolute():
                file_path = path.parent / file_path
            if not file_path.is_file():
                raise ValueError(f'PacBio FASTQ not found for {sequence_id}: {file_path}')
            seen.add(sequence_id)
            rows.append((sequence_id, file_path.resolve()))
    if not rows:
        raise ValueError('PacBio map contains no isolate FASTQ rows')
    return rows


def _cluster_representative(
    sequence_id: str,
    fastq: Path,
    workdir: Path,
    *,
    min_read_length: int,
    max_read_length: int,
    min_mean_quality: float,
    cluster_identity: float,
    library_prep: str,
) -> tuple[str | None, dict]:
    """Return the modal sequence from the largest vsearch cluster and its QC."""
    counts: dict[str, int] = defaultdict(int)
    quality_sums: dict[str, float] = defaultdict(float)
    total = passed = 0
    for sequence, mean_quality in _fastq_records(fastq):
        total += 1
        if min_read_length <= len(sequence) <= max_read_length and mean_quality >= min_mean_quality:
            counts[sequence] += 1
            quality_sums[sequence] += mean_quality
            passed += 1
    row = {
        'sequence_id': sequence_id, 'source_fastq': str(fastq), 'raw_read_count': total,
        'passed_read_count': passed, 'unique_sequence_count': len(counts),
        'cluster_identity': cluster_identity, 'dominant_cluster_read_count': 0,
        'dominant_cluster_fraction': 0.0, 'representative_length': '',
        'representative_mean_quality': '', 'qc_class': 'FAIL_QC',
        'recommendation': 'RESEQUENCE', 'reasons': '',
        'sequencing_platform': 'PacBio HiFi', 'library_prep': library_prep,
        'processing_workflow': 'SMRT Link Read Segmentation',
    }
    if not counts:
        row['reasons'] = 'no_reads_pass_length_and_quality_filter'
        return None, row

    derep = workdir / f'{sequence_id}.derep.fasta'
    labels: dict[str, str] = {}
    records = []
    for index, (sequence, count) in enumerate(counts.items(), start=1):
        label = f'read{index:08d};size={count};'
        labels[label] = sequence
        records.append((label, sequence))
    write_fasta(records, derep)
    uc = workdir / f'{sequence_id}.clusters.uc'
    if shutil.which('vsearch') is None:
        raise RuntimeError('vsearch is required for PacBio read clustering; install it before using --pacbio-map')
    run_cmd([
        'vsearch', '--cluster_fast', str(derep), '--id', str(cluster_identity),
        '--sizein', '--uc', str(uc), '--threads', '1', '--quiet',
    ])
    clusters: dict[str, list[str]] = defaultdict(list)
    with open(uc) as handle:
        for line in handle:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 9 and parts[0] in {'S', 'H'}:
                clusters[parts[1]].append(parts[8])
    if not clusters:
        raise RuntimeError(f'vsearch produced no clusters for {sequence_id}')
    winning = max(clusters.values(), key=lambda members: sum(counts[labels[name]] for name in members))
    winning_sequences = [labels[name] for name in winning]
    representative = max(
        winning_sequences,
        key=lambda sequence: (counts[sequence], quality_sums[sequence] / counts[sequence], sequence),
    )
    dominant_reads = sum(counts[sequence] for sequence in winning_sequences)
    row.update({
        'dominant_cluster_read_count': dominant_reads,
        'dominant_cluster_fraction': dominant_reads / passed,
        'representative_length': len(representative),
        'representative_mean_quality': round(quality_sums[representative] / counts[representative], 2),
    })
    return representative, row


def run_kinnex_import(
    map_path: str | Path,
    outdir: str | Path,
    *,
    min_read_length: int = 1000,
    max_read_length: int = 2000,
    min_mean_quality: float = 20.0,
    cluster_identity: float = 0.995,
    min_reads: int = 20,
    min_dominant_fraction: float = 0.80,
    library_prep: str = 'Kinnex 16S rRNA Kit',
) -> dict:
    """Build one reviewed PacBio 16S representative per isolate FASTQ."""
    if not 0 < cluster_identity <= 1 or not 0 < min_dominant_fraction <= 1:
        raise ValueError('cluster identity and dominant-cluster fraction must be between 0 and 1')
    output = Path(outdir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    workdir = output / 'clustering_work'
    workdir.mkdir(exist_ok=True)
    representatives, rows = [], []
    for sequence_id, fastq in _read_map(map_path):
        representative, row = _cluster_representative(
            sequence_id, fastq, workdir, min_read_length=min_read_length,
            max_read_length=max_read_length, min_mean_quality=min_mean_quality,
            cluster_identity=cluster_identity, library_prep=library_prep,
        )
        reasons = []
        if row['passed_read_count'] < min_reads:
            reasons.append('insufficient_passing_reads')
        if row['dominant_cluster_fraction'] < min_dominant_fraction:
            reasons.append('no_dominant_16s_cluster')
        if representative and not reasons:
            prep_reason = library_prep.replace(' ', '_')
            row.update(qc_class='PASS_HIGH_CONFIDENCE', recommendation='ACCEPT', reasons=f'abundance_supported_pacbio_representative;library_prep={prep_reason};processing=SMRT_Link_Read_Segmentation')
            representatives.append((sequence_id, representative))
        else:
            row['reasons'] = ';'.join(reasons) or row['reasons']
        rows.append(row)
    fasta = output / 'pacbio_16s_representatives.fasta'
    qc = output / 'marker_qc.tsv'
    report = output / 'pacbio_16s_qc.tsv'
    read_qc = output / 'read_qc.tsv'
    write_fasta(representatives, fasta)
    fields = list(rows[0]) if rows else []
    with open(report, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader(); writer.writerows(rows)
    with open(qc, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['sequence_id', 'qc_class', 'recommendation', 'reasons', 'read_ids'], delimiter='\t')
        writer.writeheader()
        for row in rows:
            writer.writerow({
                'sequence_id': row['sequence_id'], 'qc_class': row['qc_class'],
                'recommendation': row['recommendation'], 'reasons': row['reasons'],
                'read_ids': f"{row['passed_read_count']} HiFi reads; dominant cluster {row['dominant_cluster_read_count']}",
            })
    with open(read_qc, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['sequence_id', 'source_file'], delimiter='\t')
        writer.writeheader()
        writer.writerows({'sequence_id': row['sequence_id'], 'source_file': row['source_fastq']} for row in rows)
    return {
        'fasta': str(fasta), 'marker_qc': str(qc), 'report': str(report), 'read_qc': str(read_qc),
        'accepted': len(representatives), 'total': len(rows),
    }
