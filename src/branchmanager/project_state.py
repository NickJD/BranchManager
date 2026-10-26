"""Operational isolate and genome state for BranchManager's rolling project."""

from __future__ import annotations

import csv
import gzip
import hashlib
from pathlib import Path
from typing import Iterable, Optional


ISOLATE_STATES = (
    'RECEIVED',
    'TRACE_REVIEW',
    'MARKER_QC_FAILED',
    'MARKER_QC_PASSED',
    'PROPOSED',
    'SAB_APPROVED',
    'DNA_EXTRACTION_PENDING',
    'DNA_EXTRACTION_FAILED',
    'LIBRARY_PENDING',
    'SEQUENCED',
    'GENOME_QC_FAILED',
    'GENOME_QC_PASSED',
    'WITHDRAWN',
)

TRUE_VALUES = {'1', 'true', 't', 'yes', 'y', 'pass', 'passed'}
FALSE_VALUES = {'0', 'false', 'f', 'no', 'n', 'fail', 'failed'}


def _normalise(name: object) -> str:
    return str(name or '').strip().lstrip('\ufeff').lower().replace('-', '_').replace(' ', '_')


def _open(path: str | Path):
    return gzip.open(path, 'rt', newline='') if str(path).lower().endswith('.gz') else open(path, newline='')


def _read_rows(path: str | Path) -> list[dict[str, str]]:
    with _open(path) as handle:
        sample = handle.read(4096)
        handle.seek(0)
        delimiter = '\t' if str(path).lower().endswith(('.tsv', '.tsv.gz', '.txt', '.txt.gz')) else ','
        try:
            delimiter = csv.Sniffer().sniff(sample, delimiters=',\t').delimiter
        except Exception:
            pass
        return [
            {_normalise(key): str(value or '').strip() for key, value in row.items() if key is not None}
            for row in csv.DictReader(handle, delimiter=delimiter)
        ]


def _value(row: dict[str, str], *names: str) -> str:
    for name in names:
        value = row.get(_normalise(name), '')
        if value:
            return value
    return ''


def _optional_float(value: object) -> Optional[float]:
    if value in (None, '', 'NA', 'N/A', 'None'):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _optional_bool(value: object) -> Optional[bool]:
    text = str(value or '').strip().lower()
    if not text:
        return None
    if text in TRUE_VALUES:
        return True
    if text in FALSE_VALUES:
        return False
    return None


def _sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def normalise_state(value: object, *, default: str = 'RECEIVED') -> str:
    state = str(value or default).strip().upper().replace('-', '_').replace(' ', '_')
    aliases = {
        'APPROVED': 'SAB_APPROVED',
        'DNA_FAILED': 'DNA_EXTRACTION_FAILED',
        'DNA_PENDING': 'DNA_EXTRACTION_PENDING',
        'QC_PASSED': 'GENOME_QC_PASSED',
        'QC_FAILED': 'GENOME_QC_FAILED',
        'GENOME_AVAILABLE': 'GENOME_QC_PASSED',
    }
    state = aliases.get(state, state)
    if state not in ISOLATE_STATES:
        raise ValueError(f'Unknown isolate status {value!r}; expected one of {", ".join(ISOLATE_STATES)}')
    return state


def import_status_updates(db, path: str | Path) -> list[dict]:
    rows = _read_rows(path)
    known = set(db.get_all_ids())
    results = []
    for line_number, row in enumerate(rows, start=2):
        sequence_id = _value(row, 'sequence_id', 'isolate_id', 'sample_id', 'id')
        state_text = _value(row, 'status', 'operational_status', 'isolate_status')
        detail = _value(row, 'detail', 'status_detail', 'notes')
        if not sequence_id:
            results.append({'line': line_number, 'sequence_id': '', 'result': 'REJECTED', 'detail': 'missing sequence_id'})
            continue
        if sequence_id not in known:
            results.append({'line': line_number, 'sequence_id': sequence_id, 'result': 'REJECTED', 'detail': 'sequence ID is not in the project database'})
            continue
        try:
            state = normalise_state(state_text)
        except ValueError as exc:
            results.append({'line': line_number, 'sequence_id': sequence_id, 'result': 'REJECTED', 'detail': str(exc)})
            continue
        db.update_isolate_status(sequence_id, state, detail=detail, source_file=str(Path(path).resolve()))
        results.append({'line': line_number, 'sequence_id': sequence_id, 'result': 'IMPORTED', 'detail': state})
    return results


def import_genome_results(
    db, path: str | Path, *, min_completeness: float = 90.0,
    max_contamination: float = 5.0,
) -> list[dict]:
    """Import completed genome evidence without equating sequencing with QC success."""
    rows = _read_rows(path)
    known = set(db.get_all_ids())
    accepted = []
    results = []
    for line_number, row in enumerate(rows, start=2):
        sequence_id = _value(row, 'sequence_id', 'isolate_id', 'sample_id', 'id')
        genome_id = _value(row, 'genome_id', 'assembly_id', 'accession', 'genome_accession')
        if not sequence_id:
            results.append({'line': line_number, 'sequence_id': '', 'genome_id': genome_id, 'result': 'REJECTED', 'detail': 'missing sequence_id'})
            continue
        if sequence_id not in known:
            results.append({'line': line_number, 'sequence_id': sequence_id, 'genome_id': genome_id, 'result': 'REJECTED', 'detail': 'sequence ID is not in the project database'})
            continue
        if not genome_id:
            results.append({'line': line_number, 'sequence_id': sequence_id, 'genome_id': '', 'result': 'REJECTED', 'detail': 'missing genome_id/accession'})
            continue

        status_text = _value(row, 'genome_status', 'status') or 'SEQUENCED'
        try:
            status = normalise_state(status_text, default='SEQUENCED')
        except ValueError as exc:
            results.append({'line': line_number, 'sequence_id': sequence_id, 'genome_id': genome_id, 'result': 'REJECTED', 'detail': str(exc)})
            continue
        if status not in {'SEQUENCED', 'GENOME_QC_FAILED', 'GENOME_QC_PASSED'}:
            results.append({'line': line_number, 'sequence_id': sequence_id, 'genome_id': genome_id, 'result': 'REJECTED', 'detail': 'genome_status must be SEQUENCED, GENOME_QC_FAILED, or GENOME_QC_PASSED'})
            continue

        completeness = _optional_float(_value(row, 'completeness', 'genome_completeness'))
        contamination = _optional_float(_value(row, 'contamination', 'genome_contamination'))
        explicit_pass = _optional_bool(_value(row, 'genome_qc_pass', 'qc_pass', 'passed_qc'))
        metrics_available = completeness is not None and contamination is not None
        metrics_pass = bool(
            metrics_available
            and completeness >= float(min_completeness)
            and contamination <= float(max_contamination)
        )
        if explicit_pass is None and status == 'SEQUENCED' and metrics_available:
            qc_pass = metrics_pass
            status = 'GENOME_QC_PASSED' if metrics_pass else 'GENOME_QC_FAILED'
            qc_basis = f'automatic:{min_completeness:g}% completeness/{max_contamination:g}% contamination'
        else:
            qc_pass = status == 'GENOME_QC_PASSED' if explicit_pass is None else explicit_pass
            qc_basis = 'explicit_status_or_qc_flag'
        if status == 'GENOME_QC_FAILED' and qc_pass:
            results.append({'line': line_number, 'sequence_id': sequence_id, 'genome_id': genome_id, 'result': 'REJECTED', 'detail': 'GENOME_QC_FAILED conflicts with genome_qc_pass=yes'})
            continue
        if status == 'GENOME_QC_PASSED' and not qc_pass:
            results.append({'line': line_number, 'sequence_id': sequence_id, 'genome_id': genome_id, 'result': 'REJECTED', 'detail': 'GENOME_QC_PASSED conflicts with genome_qc_pass=no'})
            continue
        if qc_pass and metrics_available and not metrics_pass:
            results.append({
                'line': line_number, 'sequence_id': sequence_id, 'genome_id': genome_id,
                'result': 'REJECTED',
                'detail': (
                    f'QC pass conflicts with completeness={completeness}, contamination={contamination}; '
                    f'thresholds are >= {min_completeness} and <= {max_contamination}'
                ),
            })
            continue

        genome_path = _value(row, 'genome_path', 'assembly_path', 'fasta_path')
        genome_sha256 = ''
        if genome_path:
            candidate = Path(genome_path).expanduser()
            if not candidate.is_absolute():
                candidate = Path(path).resolve().parent / candidate
            if not candidate.is_file():
                results.append({'line': line_number, 'sequence_id': sequence_id, 'genome_id': genome_id, 'result': 'REJECTED', 'detail': f'genome file not found: {candidate}'})
                continue
            genome_path = str(candidate.resolve())
            genome_sha256 = _sha256(candidate)

        accepted.append({
            'genome_id': genome_id,
            'sequence_id': sequence_id,
            'accession': _value(row, 'accession', 'genome_accession'),
            'genome_status': status,
            'genome_qc_pass': bool(qc_pass),
            'completeness': completeness,
            'contamination': contamination,
            'gtdb_taxonomy': _value(row, 'gtdb_taxonomy', 'gtdbtk_taxonomy', 'taxonomy'),
            'ani_cluster': _value(row, 'ani_cluster', 'species_cluster', 'ani_species'),
            'genome_path': genome_path,
            'genome_sha256': genome_sha256,
            'notes': '; '.join(filter(None, (_value(row, 'notes', 'detail'), f'qc_basis={qc_basis}'))),
            'source_file': str(Path(path).resolve()),
        })
        results.append({'line': line_number, 'sequence_id': sequence_id, 'genome_id': genome_id, 'result': 'IMPORTED', 'detail': status})

    db.upsert_genome_records(accepted)
    for row in accepted:
        sequence_id = row['sequence_id']
        state = 'GENOME_QC_PASSED' if row['genome_qc_pass'] else row['genome_status']
        db.update_isolate_status(sequence_id, state, detail=row.get('notes', ''), source_file=str(Path(path).resolve()))
        if row['genome_qc_pass']:
            current = db.get_sequencing_metadata_for_ids([sequence_id]).get(sequence_id, {})
            db.upsert_sequencing_metadata([{
                'id': sequence_id,
                'partner_id': current.get('partner_id', ''),
                'dataset': current.get('dataset', ''),
                'selected_for_sequencing': current.get('selected_for_sequencing', False),
                'selected_for_wgs': True,
                'source_id': sequence_id,
                'source_file': str(Path(path).resolve()),
                'raw_selected_value': 'genome_qc_pass',
                'raw_commitment_value': current.get('raw_commitment_value', ''),
                'project_metadata': current.get('project_metadata', {}),
            }])
    return results


def import_external_genomes(
    db, path: str | Path, *, marker_fasta: str | Path | None = None,
    default_dataset: str = 'ExternalGenomeBaseline', default_tier: str = 'secondary',
) -> list[dict]:
    """Register genomes acquired outside the candidate-selection workflow.

    A marker sequence is optional. Genome-only records receive an empty marker
    slot, remain available for genome-coverage/taxonomy context, and are
    excluded from marker similarity searches and trees.
    """
    from branchmanager.taxonomy import canonicalise_sequence_id
    from branchmanager.utils.fasta import read_fasta

    rows = _read_rows(path)
    markers = {}
    if marker_fasta:
        for header, sequence in read_fasta(str(marker_fasta)):
            markers[canonicalise_sequence_id(header)] = sequence
    results = []
    sequences_by_dataset = {}
    tax_entries = []
    genome_entries = []
    metadata_entries = []
    role_entries = {}
    known = db.get_all_ids()
    existing_metadata = db.get_sequencing_metadata_for_ids()
    for line_number, row in enumerate(rows, start=2):
        sequence_id = _value(row, 'sequence_id', 'isolate_id', 'sample_id', 'id')
        genome_id = _value(row, 'genome_id', 'assembly_id', 'accession', 'genome_accession')
        if not sequence_id or not genome_id:
            results.append({'line': line_number, 'sequence_id': sequence_id, 'genome_id': genome_id,
                            'result': 'REJECTED', 'detail': 'sequence_id and genome_id/accession are required'})
            continue
        dataset = _value(row, 'dataset', 'collection', 'baseline_dataset') or default_dataset
        tier = (_value(row, 'baseline_tier', 'tier') or default_tier).strip().lower()
        if tier not in {'priority', 'secondary'}:
            results.append({'line': line_number, 'sequence_id': sequence_id, 'genome_id': genome_id,
                            'result': 'REJECTED', 'detail': 'baseline_tier must be priority or secondary'})
            continue
        marker = markers.get(canonicalise_sequence_id(sequence_id), '')
        if sequence_id not in known:
            sequences_by_dataset.setdefault(dataset, []).append((sequence_id, marker))
            role_entries[dataset] = tier
            known.add(sequence_id)
        elif marker:
            with db.connect() as conn:
                conn.execute(
                    "UPDATE sequences SET sequence=?, length=? WHERE id=? AND (sequence IS NULL OR sequence='')",
                    (marker, len(marker), sequence_id),
                )
        taxonomy = _value(row, 'gtdb_taxonomy', 'gtdbtk_taxonomy', 'taxonomy')
        old = existing_metadata.get(sequence_id, {})
        if taxonomy:
            confidence = _optional_float(_value(row, 'taxonomy_confidence', 'confidence'))
            tax_entries.append((sequence_id, taxonomy, confidence, old.get('dataset') or dataset))
        metadata_entries.append({
            'id': sequence_id, 'partner_id': _value(row, 'partner_id', 'partner') or old.get('partner_id', 'external'),
            'dataset': old.get('dataset') or dataset,
            'selected_for_sequencing': bool(old.get('selected_for_sequencing')),
            'selected_for_wgs': True, 'source_id': sequence_id, 'source_file': str(Path(path).resolve()),
            'raw_selected_value': 'external_genome_available',
            'raw_commitment_value': old.get('raw_commitment_value', ''),
            'project_metadata': {
                key: value for key, value in row.items()
                if key not in {
                    'sequence_id', 'isolate_id', 'sample_id', 'id', 'genome_id',
                    'assembly_id', 'accession', 'genome_accession', 'dataset',
                    'collection', 'baseline_dataset', 'baseline_tier', 'tier',
                    'gtdb_taxonomy', 'gtdbtk_taxonomy', 'taxonomy',
                    'taxonomy_confidence', 'confidence', 'genome_status',
                    'genome_qc_pass', 'qc_pass', 'completeness', 'genome_completeness',
                    'contamination', 'genome_contamination', 'ani_cluster',
                    'species_cluster', 'genome_path', 'assembly_path', 'notes', 'detail',
                }
            },
        })
        genome_status = _value(row, 'genome_status') or 'SEQUENCED'
        genome_qc = _optional_bool(_value(row, 'genome_qc_pass', 'qc_pass'))
        if genome_qc is None:
            genome_qc = genome_status.strip().upper() == 'GENOME_QC_PASSED'
        genome_entries.append({
            'genome_id': genome_id, 'sequence_id': sequence_id,
            'accession': _value(row, 'accession', 'genome_accession') or genome_id,
            'genome_status': genome_status,
            'genome_qc_pass': genome_qc,
            'completeness': _optional_float(_value(row, 'completeness', 'genome_completeness')),
            'contamination': _optional_float(_value(row, 'contamination', 'genome_contamination')),
            'gtdb_taxonomy': taxonomy,
            'ani_cluster': _value(row, 'ani_cluster', 'species_cluster'),
            'genome_path': _value(row, 'genome_path', 'assembly_path'),
            'notes': _value(row, 'notes', 'detail'), 'source_file': str(Path(path).resolve()),
        })
        results.append({'line': line_number, 'sequence_id': sequence_id, 'genome_id': genome_id,
                        'result': 'IMPORTED', 'detail': 'genome available; marker ' + ('attached' if marker else 'not supplied')})

    for dataset, tier in role_entries.items():
        db.upsert_dataset_role(dataset, 'baseline', genomes_available=True, baseline_tier=tier)
    if sequences_by_dataset:
        # Empty sequence values are intentional genome-only records.
        for dataset, subset in sequences_by_dataset.items():
            db.insert_sequences(subset, dataset=dataset, baseline_tier=role_entries[dataset])
    if tax_entries:
        db.insert_taxonomy(tax_entries)
    if metadata_entries:
        db.upsert_sequencing_metadata(metadata_entries)
    if genome_entries:
        db.upsert_genome_records(genome_entries)
    return results


def import_maldi_identifications(
    db, path: str | Path, *, default_dataset: str = 'MALDI_TOF_Partner_Isolates',
) -> list[dict]:
    """Register partner isolates with MALDI-TOF taxonomy and no marker sequence."""
    from branchmanager.taxonomy import parse_taxon_string

    rows = _read_rows(path)
    known = db.get_all_ids()
    roles = db.get_dataset_roles()
    existing_metadata = db.get_sequencing_metadata_for_ids()
    with db.connect() as conn:
        existing_datasets = dict(conn.execute('SELECT id, dataset FROM sequences').fetchall())
    accepted = []
    new_sequences: dict[str, list[tuple[str, str]]] = {}
    new_dataset_roles: set[str] = set()
    metadata_rows = []
    status_rows = []
    metadata_ids_added: set[str] = set()
    results = []
    seen_source_rows: set[tuple[str, str, str]] = set()
    rank_aliases = {
        'f': 'f', 'family': 'f', 'g': 'g', 'genus': 'g', 's': 's', 'species': 's',
    }

    for line_number, row in enumerate(rows, start=2):
        sequence_id = _value(row, 'sequence_id', 'isolate_id', 'sample_id', 'id')
        partner_id = _value(row, 'partner_id', 'partner', 'partner_acronym')
        taxon_id = _value(row, 'taxon_id', 'taxid', 'maldi_taxon_id')
        taxon_name = _value(row, 'taxon_name', 'identified_taxon', 'taxon', 'maldi_taxonomy')
        rank_text = _value(row, 'taxonomic_rank', 'taxon_rank', 'rank').strip().lower()
        rank = rank_aliases.get(rank_text, '') if rank_text else ''
        if rank_text and not rank:
            results.append({
                'line': line_number, 'sequence_id': sequence_id, 'result': 'REJECTED',
                'detail': 'taxonomic_rank must be species, genus, or family for soft same-taxon context',
            })
            continue
        if taxon_name:
            parsed = parse_taxon_string(taxon_name)
            if rank:
                taxon_name = parsed.get(rank, taxon_name)
            else:
                rank = next((candidate for candidate in ('s', 'g', 'f') if parsed.get(candidate)), '')
            if rank:
                taxon_name = parsed.get(rank, taxon_name)
        if not sequence_id or not partner_id or not (taxon_name or taxon_id):
            results.append({
                'line': line_number, 'sequence_id': sequence_id, 'result': 'REJECTED',
                'detail': 'sequence_id, partner_id, and taxon_name or taxon_id are required',
            })
            continue
        dataset = _value(row, 'dataset', 'collection') or default_dataset
        if roles.get(dataset, {}).get('role') == 'baseline':
            results.append({
                'line': line_number, 'sequence_id': sequence_id, 'result': 'REJECTED',
                'detail': f'dataset {dataset!r} is registered as a baseline; MALDI partner isolates cannot be added there',
            })
            continue
        signature = (sequence_id, taxon_id, taxon_name)
        if signature in seen_source_rows:
            results.append({
                'line': line_number, 'sequence_id': sequence_id, 'result': 'REJECTED',
                'detail': 'duplicate MALDI identification row for this isolate and taxon',
            })
            continue
        seen_source_rows.add(signature)

        if sequence_id not in known:
            new_sequences.setdefault(dataset, []).append((sequence_id, ''))
            new_dataset_roles.add(dataset)
            known.add(sequence_id)
            metadata = {
                'id': sequence_id, 'partner_id': partner_id, 'dataset': dataset,
                'selected_for_sequencing': False, 'selected_for_wgs': False,
                'source_id': sequence_id, 'source_file': str(Path(path).resolve()),
                'raw_selected_value': 'no', 'raw_commitment_value': 'no',
                'project_metadata': {
                    'marker_sequence_status': 'NOT_AVAILABLE',
                    'addition_method': 'maldi_tof',
                },
            }
            if sequence_id not in metadata_ids_added:
                metadata_rows.append(metadata)
                metadata_ids_added.add(sequence_id)
            status_rows.append(sequence_id)
        else:
            if roles.get(existing_datasets.get(sequence_id) or '', {}).get('role') == 'baseline':
                results.append({
                    'line': line_number, 'sequence_id': sequence_id, 'result': 'REJECTED',
                    'detail': 'existing baseline isolate cannot receive a partner MALDI identification',
                })
                continue
            if sequence_id not in existing_metadata and sequence_id not in metadata_ids_added:
                metadata_rows.append({
                    'id': sequence_id, 'partner_id': partner_id, 'dataset': dataset,
                    'selected_for_sequencing': False, 'selected_for_wgs': False,
                    'source_id': sequence_id, 'source_file': str(Path(path).resolve()),
                    'raw_selected_value': 'no', 'raw_commitment_value': 'no',
                    'project_metadata': {
                        'marker_sequence_status': 'NOT_AVAILABLE',
                        'addition_method': 'maldi_tof',
                    },
                })
                metadata_ids_added.add(sequence_id)
        accepted.append({
            'sequence_id': sequence_id, 'taxon_id': taxon_id or None,
            'taxon_name': taxon_name or None, 'taxonomic_rank': rank or None,
            'identification_score': _value(row, 'identification_score', 'score', 'maldi_score') or None,
            'instrument': _value(row, 'instrument', 'instrument_model') or None,
            'library_name': _value(row, 'library_name', 'reference_library') or None,
            'identification_date': _value(row, 'identification_date', 'date') or None,
            'notes': _value(row, 'notes', 'detail') or None,
            'source_file': str(Path(path).resolve()), 'source_row': line_number,
        })
        results.append({
            'line': line_number, 'sequence_id': sequence_id, 'result': 'IMPORTED',
            'taxon_id': taxon_id, 'taxon_name': taxon_name,
            'detail': 'MALDI-TOF identification stored as soft taxonomic context; not marker or genome evidence',
        })

    for dataset in new_dataset_roles:
        db.upsert_dataset_role(dataset, 'candidate', genomes_available=False)
    for dataset, records in new_sequences.items():
        db.insert_sequences(records, dataset=dataset)
    if metadata_rows:
        db.upsert_sequencing_metadata(metadata_rows)
    if status_rows:
        for sequence_id in status_rows:
            db.update_isolate_status(
                sequence_id, 'RECEIVED',
                detail='MALDI-TOF identification imported; marker sequence unavailable',
                source_file=str(Path(path).resolve()),
            )
    if accepted:
        db.insert_maldi_identifications(accepted)
    return results


def write_maldi_import_report(path: str | Path, rows: Iterable[dict]) -> str:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = ['line', 'sequence_id', 'result', 'taxon_id', 'taxon_name', 'detail']
    with open(output, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    return str(output)


def import_partner_workbook(db, path: str | Path, *, sheet_name='Bulk Upload Template', dataset='PartnerWorkbook') -> list[dict]:
    """Import source-labeled metadata and soft MALDI observations from a partner XLSM."""
    import json
    import re
    from branchmanager.partner_metadata import load_partner_workbook_rows

    source_file = str(Path(path).expanduser().resolve())
    source_rows = load_partner_workbook_rows(source_file, sheet_name=sheet_name)
    if not source_rows:
        return []

    dataset_roles = db.get_dataset_roles()
    if dataset_roles.get(dataset, {}).get('role') == 'baseline':
        raise ValueError(f'dataset {dataset!r} is registered as a baseline; partner workbook records cannot be added there')
    known_maldi = db.get_maldi_identifications()
    maldi_signatures = {
        (
            str(row.get('sequence_id') or ''), str(row.get('taxon_id') or ''),
            str(row.get('taxon_name') or ''), str(row.get('instrument') or ''),
            str(row.get('library_name') or ''), str(row.get('source_file') or ''),
            str(row.get('source_row') or ''),
        )
        for row in known_maldi
    }
    new_sequences = []
    new_metadata = []
    maldi_to_insert = []
    metadata_updates = []
    status_updates = []
    report = []

    for row in source_rows:
        source_id = str(row['source_id'])
        sequence_id = db.resolve_sequence_id(source_id) or source_id
        with db.connect() as conn:
            sequence_row = conn.execute(
                'SELECT sequence FROM sequences WHERE id=?', (sequence_id,),
            ).fetchone()
        is_new_sequence = sequence_row is None
        has_marker = bool(sequence_row and str(sequence_row[0] or '').strip())
        if is_new_sequence:
            new_sequences.append((sequence_id, ''))
            status_updates.append(sequence_id)

        existing_metadata = db.get_sequencing_metadata_for_ids([sequence_id]).get(sequence_id)
        project_metadata = dict((existing_metadata or {}).get('project_metadata') or {})
        project_metadata.update(row.get('project_metadata') or {})
        project_metadata.update({
            'partner_workbook_source': source_file,
            'partner_workbook_source_row': row['source_row'],
            'partner_workbook_dataset': dataset,
            'partner_workbook_partner_raw': row.get('partner_raw') or '',
            'amplicon_method_reported': bool(row.get('amplicon_method_reported')),
            'maldi_method_reported': bool(row.get('maldi_method_reported')),
        })
        has_maldi_record = bool(row.get('maldi_method_reported')) or any(
            str(item.get('sequence_id') or '') == sequence_id for item in known_maldi
        )
        if has_marker and has_maldi_record:
            project_metadata['addition_method'] = 'both'
        elif has_marker:
            project_metadata['addition_method'] = 'marker_sequence'
        elif has_maldi_record:
            project_metadata['addition_method'] = 'maldi_tof'

        if existing_metadata:
            metadata_updates.append((json.dumps(project_metadata, sort_keys=True, ensure_ascii=False), sequence_id))
        else:
            new_metadata.append({
                'id': sequence_id,
                'partner_id': row.get('partner_id') or source_id.split('_', 1)[0],
                'dataset': dataset,
                'selected_for_sequencing': False,
                'selected_for_wgs': False,
                'source_id': source_id,
                'source_file': source_file,
                'raw_selected_value': '',
                'raw_commitment_value': '',
                'project_metadata': project_metadata,
            })

        maldi_result = 'NOT_SELECTED'
        raw_taxon = str(row.get('maldi_result_raw') or '').strip()
        if row.get('maldi_method_reported'):
            if not raw_taxon:
                maldi_result = 'MALDI_SELECTED_NO_RESULT'
            else:
                taxon_id = raw_taxon if raw_taxon.isdigit() else None
                taxon_name = None if taxon_id else raw_taxon
                rank = None
                unresolved = False
                if taxon_name:
                    if re.fullmatch(r'[A-Za-z][A-Za-z0-9.-]*_\d+', taxon_name):
                        # These are partner culture IDs in this source field,
                        # not a taxonomic name or NCBI taxon ID.
                        taxon_name = None
                        unresolved = True
                    else:
                        tokens = taxon_name.split()
                        if len(tokens) == 1 and tokens[0][:1].isupper():
                            rank = 'g'
                        elif len(tokens) == 2 and tokens[0][:1].isupper() and tokens[1][:1].islower():
                            rank = 's'
                        else:
                            unresolved = True
                signature = (
                    sequence_id, str(taxon_id or ''), str(taxon_name or ''),
                    str(row.get('maldi_instrument') or ''), str(row.get('maldi_library') or ''),
                    source_file, str(row['source_row']),
                )
                if signature in maldi_signatures:
                    maldi_result = 'ALREADY_IMPORTED'
                else:
                    maldi_to_insert.append({
                        'sequence_id': sequence_id,
                        'taxon_id': taxon_id,
                        'taxon_name': taxon_name,
                        'taxonomic_rank': rank,
                        'instrument': row.get('maldi_instrument') or None,
                        'library_name': row.get('maldi_library') or None,
                        'notes': f"Unresolved source taxon value: {raw_taxon}" if unresolved else None,
                        'source_file': source_file,
                        'source_row': row['source_row'],
                    })
                    maldi_signatures.add(signature)
                    maldi_result = 'MALDI_RESULT_NEEDS_REVIEW' if unresolved else 'SOFT_TAXON_CONTEXT_IMPORTED'
        elif raw_taxon or row.get('maldi_instrument') or row.get('maldi_library'):
            maldi_result = 'MALDI_FIELDS_WITHOUT_SELECTION'

        report.append({
            'source_row': row['source_row'],
            'sequence_id': sequence_id,
            'partner_id': row.get('partner_id') or '',
            'metadata_result': 'INSERTED_ISOLATE' if is_new_sequence else 'METADATA_UPDATED',
            'maldi_result': maldi_result,
            'source_taxon_value': raw_taxon,
            'detail': (
                'MALDI data is stored separately from marker taxonomy; unresolved values do not contribute '
                'to same-taxon context or selection scores.'
            ),
        })

    if new_sequences:
        db.upsert_dataset_role(dataset, 'candidate', genomes_available=False)
        db.insert_sequences(new_sequences, dataset=dataset)
    if new_metadata:
        db.upsert_sequencing_metadata(new_metadata)
    if metadata_updates:
        with db.connect() as conn:
            conn.executemany(
                'UPDATE sequencing_metadata SET project_metadata_json=? WHERE id=?',
                metadata_updates,
            )
    if status_updates:
        for sequence_id in status_updates:
            db.update_isolate_status(
                sequence_id, 'RECEIVED',
                detail='Partner metadata imported; marker sequence or MALDI result may be pending',
                source_file=source_file,
            )
    if maldi_to_insert:
        db.insert_maldi_identifications(maldi_to_insert)
    return report


def write_partner_workbook_import_report(path: str | Path, rows: Iterable[dict]) -> str:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        'source_row', 'sequence_id', 'partner_id', 'metadata_result',
        'maldi_result', 'source_taxon_value', 'detail',
    ]
    with open(output, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    return str(output)


def write_import_report(path: str | Path, rows: Iterable[dict]) -> str:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = ['line', 'sequence_id', 'genome_id', 'result', 'detail']
    with open(output, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    return str(output)
