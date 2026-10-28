"""Local browser UI for inspecting and assembling Sanger chromatograms.

The Studio is deliberately a thin interface over Paper Trail/Merge Meeting:
it does not implement a second base-caller, trimmer, or assembler. It binds to
localhost, confines map/read browsing to the selected project root, and writes
every inspection to an auditable run directory.
"""

from __future__ import annotations

import csv
import json
import mimetypes
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


_MAP_ID_KEYS = ('sequence_id', 'seq_id', 'id', 'isolate_id', 'sample_id', 'sample')
_MAP_FILE_KEYS = ('read_file', 'file', 'filename', 'path', 'file_path', 'ab1_file')
_MAP_PRIMER_KEYS = ('primer', 'primer_name')
_MAP_DIRECTION_KEYS = ('direction', 'orientation', 'read_direction')
_MAP_MODE_KEYS = ('processing_mode', 'mode', 'read_mode')
_READ_EXTENSIONS = {'.ab1', '.abi'}
_MAP_EXTENSIONS = {'.tsv', '.csv', '.txt'}
_REVIEW_DECISIONS = {'approve', 'approved', 'accept', 'accepted', 'pass', 'yes', 'reject', 'rejected', 'fail', 'failed', 'no'}


def _normalise_key(value: object) -> str:
    return re.sub(r'[^a-z0-9]+', '_', str(value or '').strip().lstrip('\ufeff').lower()).strip('_')


def _is_within(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
        return True
    except ValueError:
        return False


def _require_within(path: str | Path, root: Path, label: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if not _is_within(resolved, root):
        raise ValueError(f'{label} must be inside the Studio project root ({root}): {resolved}')
    return resolved


def _read_tsv(path: str | Path) -> list[dict[str, str]]:
    with open(path, newline='', encoding='utf-8-sig') as handle:
        sample = handle.read(4096)
        handle.seek(0)
        delimiter = '\t' if Path(path).suffix.lower() in {'.tsv', '.txt'} else ','
        try:
            delimiter = csv.Sniffer().sniff(sample, delimiters=',\t').delimiter
        except csv.Error:
            pass
        return [dict(row) for row in csv.DictReader(handle, delimiter=delimiter)]


def _row_get(row: dict[str, str], names: tuple[str, ...]) -> str:
    values = {_normalise_key(key): str(value or '').strip() for key, value in row.items()}
    return next((values[name] for name in names if values.get(name)), '')


def _read_map(path: str | Path, root: Path) -> list[dict[str, str]]:
    map_path = _require_within(path, root, 'Sample map')
    if not map_path.is_file():
        raise ValueError(f'Sample map does not exist: {map_path}')
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
        read_path = _require_within(read_path, root, f'Read path on row {line}')
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
    return next((row for row in _read_tsv(path) if _row_get(row, ('sequenceid', 'id')) == sequence_id), {})


def _serialise_inspection(output: dict, sequence_id: str, job_dir: Path, output_root: Path) -> dict:
    assembly = _read_report_row(output['assembly_tsv'], sequence_id)
    read_rows = [row for row in _read_tsv(output['read_qc_tsv']) if row.get('SequenceID') == sequence_id]
    placement_rows = [row for row in _read_tsv(output['assembly_placements_tsv']) if row.get('SequenceID') == sequence_id]
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
        ('recommendations_tsv', 'Resequencing recommendations'),
        ('marker_review_template_tsv', 'Review template'),
        ('qc_policy_tsv', 'QC policy'),
        ('summary', 'Run summary'),
    ):
        path = Path(output[key]).resolve()
        if _is_within(path, output_root):
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
        'taxonomy_concordance': assembly.get('TaxonomyConcordance', ''),
        'read_taxonomy_assignments': assembly.get('ReadTaxonomyAssignments', ''),
        'suggested_action': assembly.get('SuggestedAction', ''),
        'reads': read_rows,
        'placements': placement_rows,
        'artifacts': artifacts,
        'files': files,
        'job_dir': str(job_dir),
    }


class SangerStudio:
    def __init__(self, root: str | Path, out: str | Path, sample_map: str | Path | None = None):
        self.root = Path(root).expanduser().resolve()
        if not self.root.is_dir():
            raise ValueError(f'Studio project root is not an existing directory: {self.root}')
        self.out = Path(out).expanduser().resolve()
        self.out.mkdir(parents=True, exist_ok=True)
        self.jobs = self.out / 'inspections'
        self.jobs.mkdir(parents=True, exist_ok=True)
        self.review_path = self.out / 'marker_review_decisions.tsv'
        self.lock = threading.Lock()
        self.map_path: Path | None = None
        self.map_rows: list[dict[str, str]] = []
        self.last_inspection: dict[str, dict] = {}
        self.last_batch: dict[str, dict] = {}
        self.decisions = _load_review_decisions(self.review_path)
        if sample_map:
            self.load_map(sample_map)

    def load_map(self, path: str | Path) -> dict:
        resolved = _require_within(path, self.root, 'Sample map')
        rows = _read_map(resolved, self.root)
        self.map_path = resolved
        self.map_rows = rows
        return {'path': str(resolved), 'samples': _samples_from_rows(rows), 'read_count': len(rows)}

    def state(self) -> dict:
        policy = paper_trail.DEFAULT_QC_POLICY
        return {
            'root': str(self.root),
            'out': str(self.out),
            'map_path': str(self.map_path) if self.map_path else '',
            'samples': _samples_from_rows(self.map_rows),
            'review_decisions': list(self.decisions.values()),
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
            path = _require_within(raw_path, self.root, 'Read path')
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
        return {'path': '', 'samples': _samples_from_rows(rows), 'read_count': len(rows)}

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

        with self.lock:
            stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
            job_dir = self.jobs / f'{stamp}_{re.sub(r"[^A-Za-z0-9_.-]+", "_", sequence_id)}_{uuid.uuid4().hex[:6]}'
            job_dir.mkdir(parents=True, exist_ok=False)
            map_path = job_dir / 'selected_read_map.tsv'
            _write_map(map_path, sample_reads)
            manifest = RunManifest(job_dir, 'sanger_studio_inspection')
            manifest.data['settings'] = dict(settings, processing_mode=mode, sequence_id=sequence_id)
            if self.map_path and self.map_path.is_file():
                manifest.add_input(self.map_path, role='source_sample_map')
            if self.review_path.is_file():
                manifest.add_input(self.review_path, role='manual_review_decisions')
            manifest.add_input(map_path, role='resolved_read_map')
            for row in sample_reads:
                manifest.add_input(row['file'], role=f'raw_trace:{row["sequence_id"]}')
            manifest.write()
            manifest.add_stage('paper_trail', 'RUNNING', detail='Per-read QC, trimming, primer orientation and read inspection')
            try:
                output = paper_trail.run_paper_trail(
                    [], job_dir, sample_map=map_path,
                    primers=paper_trail.DEFAULT_PRIMERS,
                    min_quality=settings['min_quality'],
                    min_read_length=settings['min_read_length'],
                    min_length=settings['min_final_length'],
                    min_mean_quality=settings['min_mean_quality'],
                    max_output_expected_errors=settings['max_output_expected_errors'],
                    warn_n_percent=settings['warn_n_percent'],
                    max_n_percent=settings['max_n_percent'],
                    min_overlap=settings['min_overlap'],
                    min_overlap_identity=settings['min_overlap_identity'],
                    assemble=True,
                )
                for key in (
                    'assembled_fasta', 'read_qc_tsv', 'per_base_error_tsv', 'assembly_tsv',
                    'assembly_placements_tsv', 'recommendations_tsv', 'marker_review_template_tsv',
                    'qc_policy_tsv', 'summary',
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
            manifest = RunManifest(job_dir, 'sanger_studio_batch')
            manifest.data['settings'] = settings
            manifest.data['sample_count'] = len({row['sequence_id'] for row in self.map_rows})
            if self.map_path and self.map_path.is_file():
                manifest.add_input(self.map_path, role='source_sample_map')
            if self.review_path.is_file():
                manifest.add_input(self.review_path, role='manual_review_decisions')
            manifest.add_input(map_path, role='resolved_read_map')
            for row in self.map_rows:
                manifest.add_input(row['file'], role=f'raw_trace:{row["sequence_id"]}')
            manifest.write()
            manifest.add_stage('paper_trail', 'RUNNING', detail='Batch base-calling, trimming, QC and assembly')
            try:
                output = paper_trail.run_paper_trail(
                    [], job_dir, sample_map=map_path,
                    primers=paper_trail.DEFAULT_PRIMERS,
                    min_quality=settings['min_quality'],
                    min_read_length=settings['min_read_length'],
                    min_length=settings['min_final_length'],
                    min_mean_quality=settings['min_mean_quality'],
                    max_output_expected_errors=settings['max_output_expected_errors'],
                    warn_n_percent=settings['warn_n_percent'],
                    max_n_percent=settings['max_n_percent'],
                    min_overlap=settings['min_overlap'],
                    min_overlap_identity=settings['min_overlap_identity'],
                    assemble=True,
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
                'recommendations_tsv', 'marker_review_template_tsv', 'summary',
            )},
            'qc_counts': qc_counts,
        }
        self.last_batch[batch_name] = result
        return result


def serve(root: str | Path, out: str | Path, sample_map: str | Path | None = None, *, port: int = 8765, open_browser: bool = True) -> None:
    """Start the local-only Sanger Studio server and optionally open its UI."""
    from importlib.resources import files

    studio = SangerStudio(root, out, sample_map)
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
                payload = self._body()
                if self.path == '/api/load-map':
                    self._json(200, studio.load_map(payload.get('path', '')))
                elif self.path == '/api/quick-map':
                    self._json(200, studio.quick_map(payload.get('paths') or []))
                elif self.path == '/api/inspect':
                    self._json(200, studio.inspect(payload))
                elif self.path == '/api/process-all':
                    self._json(200, studio.process_all(payload))
                elif self.path == '/api/review':
                    self._json(200, studio.save_review(payload))
                else:
                    self._json(404, {'error': 'Not found'})
            except Exception as exc:
                self._json(400, {'error': str(exc)})

        def log_message(self, format, *args):  # noqa: A003 - stdlib handler signature
            print(f'[sanger-studio] {self.address_string()} - {format % args}')

    try:
        server = ThreadingHTTPServer(('127.0.0.1', int(port)), Handler)
    except OSError as exc:
        raise SystemExit(f'[sanger-studio] Could not bind to 127.0.0.1:{port}: {exc}') from exc
    url = f'http://127.0.0.1:{server.server_port}/'
    print(f'[sanger-studio] Local GUI: {url}\n[sanger-studio] Project root: {studio.root}\n[sanger-studio] Outputs: {studio.out}\nPress Ctrl+C to stop.')
    if open_browser:
        try:
            webbrowser.open(url)
        except Exception:
            pass
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print('\n[sanger-studio] Stopped.')
    finally:
        server.server_close()
