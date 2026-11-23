"""Versioned, cumulative project continuation helpers."""
from __future__ import annotations

import json
import shutil
import sqlite3
from datetime import datetime, timezone
from pathlib import Path


def utc_stamp() -> str:
    return datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')


def sha256_file(path: str | Path) -> str:
    import hashlib
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def backup_sqlite(source: str | Path, destination: str | Path) -> Path:
    source = Path(source).expanduser().resolve()
    destination = Path(destination).expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.unlink(missing_ok=True)
    with sqlite3.connect(str(source)) as src, sqlite3.connect(str(destination)) as dst:
        src.backup(dst)
        check = dst.execute('PRAGMA integrity_check').fetchone()[0]
        if check != 'ok':
            raise RuntimeError(f'SQLite snapshot failed integrity_check: {check}')
    return destination


def create_version_dir(root: str | Path, stamp: str | None = None) -> Path:
    root = Path(root).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    stamp = stamp or utc_stamp()
    candidate = root / stamp
    suffix = 1
    while candidate.exists():
        candidate = root / f'{stamp}_{suffix:02d}'
        suffix += 1
    candidate.mkdir()
    return candidate


def copy_assessment_snapshot(manifest_path: str | Path | None, destination: str | Path) -> Path | None:
    """Copy the latest prior review's output tree when its manifest can be located."""
    if not manifest_path:
        return None
    manifest_path = Path(manifest_path).expanduser().resolve()
    if not manifest_path.is_file():
        return None
    try:
        data = json.loads(manifest_path.read_text(encoding='utf-8'))
    except (OSError, ValueError):
        return None
    assessment_file = next(
        (Path(item.get('path', '')) for item in data.get('outputs', [])
         if item.get('role') in {'full_assessment', 'board_summary'} and item.get('path')),
        None,
    )
    if assessment_file is None or not assessment_file.exists():
        return None
    assessment_root = assessment_file.parent if assessment_file.name == 'selection_summary.tsv' else assessment_file.parent
    # For a selection summary the directory is assessment/; copy its parent run folder
    # to retain the full tree, taxonomy, dashboard, and assessment evidence.
    if assessment_root.name == 'assessment':
        assessment_root = assessment_root.parent
    target = Path(destination) / 'previous_assessment'
    if assessment_root.is_dir():
        shutil.copytree(assessment_root, target, dirs_exist_ok=True)
        return target
    return None


def latest_performance_review_manifest(db_path: str | Path) -> Path | None:
    db_path = Path(db_path).expanduser().resolve()
    try:
        with sqlite3.connect(str(db_path)) as conn:
            row = conn.execute(
                "SELECT manifest_path FROM project_runs WHERE lower(workflow)='performance_review' "
                "AND lower(status)='complete' AND manifest_path IS NOT NULL AND manifest_path!='' "
                "ORDER BY completed_at DESC, started_at DESC LIMIT 1"
            ).fetchone()
    except sqlite3.Error:
        return None
    if not row:
        return None
    path = Path(row[0]).expanduser().resolve()
    return path if path.is_file() else None



def inherit_performance_review_args(manifest_path: str | Path):
    """Recover a prior Performance Review CLI namespace from its recorded command."""
    from branchmanager.cli.parser import build_parser
    data = json.loads(Path(manifest_path).read_text(encoding='utf-8'))
    command = data.get('command') or []
    if not isinstance(command, list) or not command:
        raise ValueError('Previous Performance Review manifest has no recorded command line.')
    # Recorded commands may include a Python/module wrapper; parse from the subcommand.
    try:
        idx = command.index('performance-review')
    except ValueError:
        raise ValueError('Previous Performance Review manifest command does not contain performance-review.')
    config_record = data.get('project_config') or {}
    config_path = str(config_record.get('path') or '') if isinstance(config_record, dict) else ''
    if config_path and not Path(config_path).is_file():
        raise ValueError(f'Previous Performance Review project profile is no longer available: {config_path}')
    prefix = ['--config', config_path] if config_path and Path(config_path).is_file() else []
    recovered = build_parser().parse_args(prefix + command[idx:])
    # Prefer absolute, manifest-recorded paths over historical relative command
    # arguments, which may have been run from a different working directory.
    inputs = data.get('inputs') or []
    by_role = {}
    for item in inputs:
        if item.get('path') and item.get('exists'):
            by_role.setdefault(item.get('role'), []).append(item['path'])
    for role, attribute in (
        ('primary_reference', 'ref'), ('reference_taxonomy', 'taxa'),
        ('baseline_fasta', 'baseline_fasta'), ('mwl', 'mwl'),
        ('baseline_taxonomy', 'baseline_taxa_assignments'),
        ('chimera_reference', 'chimera_ref'), ('target_collection', 'target'),
        ('anchors', 'anchors'),
    ):
        if by_role.get(role):
            setattr(recovered, attribute, by_role[role][-1])
    for role, attribute in (
        ('alternate_reference', 'alt_ref'),
        ('alternate_reference_taxonomy', 'alt_taxa'),
    ):
        if by_role.get(role):
            setattr(recovered, attribute, list(by_role[role]))
    return recovered, data


def export_candidate_marker_fasta(db, outpath: str | Path) -> Path:
    """Export the cumulative non-baseline marker collection from the staged project DB."""
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    roles = db.get_dataset_roles()
    candidate_datasets = [dataset for dataset, role in roles.items() if role.get('role') == 'candidate']
    with db.connect() as conn:
        if candidate_datasets:
            placeholders = ','.join('?' for _ in candidate_datasets)
            rows = conn.execute(
                f"SELECT id, sequence FROM sequences WHERE dataset IN ({placeholders}) "
                "AND sequence IS NOT NULL AND trim(sequence)!='' ORDER BY id",
                tuple(candidate_datasets),
            ).fetchall()
        else:
            rows = []
    # Fallback to non-baseline datasets for older DBs without dataset_roles metadata.
    if not rows:
        baseline_datasets = set(db.get_baseline_datasets())
        with db.connect() as conn:
            rows = conn.execute(
                "SELECT id, sequence, dataset FROM sequences WHERE sequence IS NOT NULL "
                "AND trim(sequence)!='' ORDER BY id"
            ).fetchall()
        rows = [(sid, seq) for sid, seq, dataset in rows if dataset not in baseline_datasets]
    from branchmanager.utils.fasta import write_fasta
    write_fasta([(str(sid), str(seq)) for sid, seq in rows], str(outpath))
    if not rows:
        raise RuntimeError('No cumulative candidate marker sequences were found for reassessment.')
    return outpath


def write_version_manifest(path: str | Path, payload: dict) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + '\n', encoding='utf-8')
    return path
