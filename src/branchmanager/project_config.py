"""Project profiles used to supply collection-specific BranchManager defaults.

Profiles are JSON so they can be loaded without an extra dependency and kept
under version control alongside the data they describe.
"""

from __future__ import annotations

import json
import os
import re
from pathlib import Path
from typing import Sequence


_ENV = re.compile(r"\$\{([A-Za-z_][A-Za-z0-9_]*)(?::-([^}]*))?\}")
_PATH_KEYS = {
    'ref', 'taxa', 'alt_ref', 'alt_taxa', 'mwl', 'partner_metadata', 'fasta', 'db',
    'baseline_fasta', 'baseline_taxa_assignments', 'baseline_colours',
    'taxa_assignments', 'previous_review', 'chimera_ref', 'anchors',
    'tree', 'alignment', 'screen_ref', 'read_dir', 'metadata', 'sample_map',
    'marker_fasta', 'input', 'out',
    'primary_taxonomy_fasta', 'secondary_taxonomy_fasta', 'secondary_taxonomy_table',
    'secondary_baseline_fasta', 'secondary_baseline_taxonomy',
    'baseline_fasta', 'most_wanted_list', 'anchor_fasta', 'anchor_metadata',
    'database', 'output_directory', 'pacbio_map',
}


def _expand(value, base: Path):
    if isinstance(value, list):
        return [_expand(item, base) for item in value]
    if isinstance(value, dict):
        return {key: _expand(item, base) for key, item in value.items()}
    if not isinstance(value, str):
        return value

    def replace(match):
        name, fallback = match.groups()
        result = os.environ.get(name, fallback)
        if result is None:
            raise ValueError(f"Project config requires environment variable {name}")
        return result

    expanded = _ENV.sub(replace, value)
    return expanded


def _normalise_assistant_databases(config: dict) -> None:
    """Convert named Assistant database entries to the CLI's legacy lists."""
    assistant = config.get('defaults', {}).get('assistant')
    if not isinstance(assistant, dict) or 'alt_databases' not in assistant:
        return
    if any(key in assistant for key in ('alt_ref', 'alt_taxa', 'alt_ref_name')):
        raise ValueError(
            "Assistant config cannot mix 'alt_databases' with 'alt_ref', "
            "'alt_taxa', or 'alt_ref_name'; use one format"
        )
    entries = assistant.pop('alt_databases')
    if not isinstance(entries, list):
        raise ValueError("Assistant 'alt_databases' must be a list of objects")
    names, refs, taxa = [], [], []
    for index, entry in enumerate(entries):
        if not isinstance(entry, dict):
            raise ValueError(f"Assistant alt_databases[{index}] must be an object")
        unknown = set(entry) - {'name', 'ref', 'taxa'}
        if unknown:
            raise ValueError(
                f"Unsupported key(s) in Assistant alt_databases[{index}]: "
                f"{', '.join(sorted(unknown))}"
            )
        name, ref = entry.get('name'), entry.get('ref')
        taxonomy = entry.get('taxa')
        if not isinstance(name, str) or not name.strip():
            raise ValueError(f"Assistant alt_databases[{index}] requires a non-empty 'name'")
        if not isinstance(ref, str) or not ref.strip():
            raise ValueError(f"Assistant alt_databases[{index}] requires a non-empty 'ref'")
        if taxonomy is not None and (not isinstance(taxonomy, str) or not taxonomy.strip()):
            raise ValueError(f"Assistant alt_databases[{index}].taxa must be a path or null")
        names.append(name)
        refs.append(ref)
        taxa.append(taxonomy)
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError(f"Assistant alt_databases names must be unique: {', '.join(duplicates)}")
    primary_name = assistant.get('ref_name', 'GTDB')
    if primary_name in names:
        raise ValueError(f"Assistant alternate database name duplicates primary ref_name {primary_name!r}")
    assistant['alt_ref_name'] = names
    assistant['alt_ref'] = refs
    assistant['alt_taxa'] = taxa


def _apply_assistant_run_data(config: dict) -> None:
    """Use Assistant run fields in data as defaults, without overriding CLI defaults."""
    data = config.get('data', {})
    if not isinstance(data, dict):
        return
    defaults = config.setdefault('defaults', {}).setdefault('assistant', {})
    if not isinstance(defaults, dict):
        return
    mappings = {
        'database': 'db',
        'dataset': 'dataset',
        'sample_map': 'sample_map',
        'output_directory': 'out',
    }
    for source, target in mappings.items():
        if source in data and data[source] is not None:
            defaults.setdefault(target, data[source])


def _resolve_config_path(value: str, base: Path) -> str:
    path = Path(value).expanduser()
    return str((path if path.is_absolute() else base / path).resolve())


def load_project_config(path: str | Path) -> dict:
    """Read and validate a profile, expanding environment variables and paths."""
    config_path = Path(path).expanduser().resolve()
    with config_path.open(encoding='utf-8') as handle:
        config = json.load(handle)
    if not isinstance(config, dict) or config.get('schema_version') != 1:
        raise ValueError('Project config must be a JSON object with schema_version: 1')
    defaults = config.get('defaults', {})
    if not isinstance(defaults, dict):
        raise ValueError("Project config 'defaults' must be an object keyed by CLI command")
    valid_commands = {'assistant', 'filing-cabinet', 'performance-review', 'quarterly-review', 'add-external-genomes'}
    unknown = set(defaults) - valid_commands
    if unknown:
        raise ValueError(f"Unsupported command(s) in project config: {', '.join(sorted(unknown))}")
    expanded = _expand(config, config_path.parent)
    configured_directory = expanded.get('project_directory')
    if configured_directory is None:
        project_directory = config_path.parent
    elif not isinstance(configured_directory, str) or not configured_directory.strip():
        raise ValueError("Project config 'project_directory' must be a non-empty path")
    else:
        project_directory = Path(configured_directory).expanduser()
        if not project_directory.is_absolute():
            project_directory = config_path.parent / project_directory
        project_directory = project_directory.resolve()
    expanded['project_directory'] = str(project_directory)
    _normalise_assistant_databases(expanded)
    _apply_assistant_run_data(expanded)
    baselines = expanded.get('baselines', [])
    if not isinstance(baselines, list) or any(not isinstance(item, dict) for item in baselines):
        raise ValueError("Project config 'baselines' must be a list of objects")
    for item in baselines:
        for key in _PATH_KEYS | {'taxonomy_assignments'}:
            raw = item.get(key)
            if isinstance(raw, str) and raw:
                item[key] = _resolve_config_path(raw, project_directory)
    for key in _PATH_KEYS & expanded.get('data', {}).keys():
        raw = expanded['data'][key]
        if isinstance(raw, str) and raw:
            expanded['data'][key] = _resolve_config_path(raw, project_directory)
    for command, values in expanded.get('defaults', {}).items():
        if not isinstance(values, dict):
            raise ValueError(f"Defaults for {command!r} must be an object")
        for key in _PATH_KEYS & values.keys():
            raw = values[key]
            if isinstance(raw, str) and raw:
                values[key] = _resolve_config_path(raw, project_directory)
            elif isinstance(raw, list):
                values[key] = [_resolve_config_path(item, project_directory) if isinstance(item, str) and item else item for item in raw]
    return expanded


def apply_project_config_defaults(parser, config: dict) -> None:
    """Apply profile values to command parsers; explicit CLI arguments win."""
    by_command = {action.dest: action for action in parser._actions if action.dest == 'command'}
    if not by_command:
        return
    subparsers = by_command['command'].choices
    configured_defaults = {command: dict(values) for command, values in config.get('defaults', {}).items()}
    budget = config.get('project', {}).get('selection_budget')
    if budget is not None:
        configured_defaults.setdefault('quarterly-review', {}).setdefault('genome_budget', int(budget))
    for command, defaults in configured_defaults.items():
        target = subparsers.get(command)
        if target is None:
            raise ValueError(f"Config provides defaults for unavailable command {command!r}")
        accepted = {action.dest for action in target._actions}
        unknown = set(defaults) - accepted
        if unknown:
            raise ValueError(f"Unsupported {command} config key(s): {', '.join(sorted(unknown))}")
        for action in target._actions:
            if action.dest in defaults and action.required:
                action.required = False
        target.set_defaults(**defaults)


def config_path_from_argv(argv: Sequence[str] | None) -> str | None:
    """Find --config whether it appears before or after the subcommand."""
    import sys
    words = list(sys.argv[1:] if argv is None else argv)
    for index, word in enumerate(words):
        if word == '--config' and index + 1 < len(words):
            return words[index + 1]
        if word.startswith('--config='):
            return word.split('=', 1)[1]
    return None
