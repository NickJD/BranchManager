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
    'baseline_fasta', 'most_wanted_list', 'anchor_fasta', 'anchor_metadata',
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
    baselines = expanded.get('baselines', [])
    if not isinstance(baselines, list) or any(not isinstance(item, dict) for item in baselines):
        raise ValueError("Project config 'baselines' must be a list of objects")
    for item in baselines:
        for key in _PATH_KEYS | {'taxonomy_assignments'}:
            raw = item.get(key)
            if isinstance(raw, str) and raw and not raw.startswith(('~', '/')):
                item[key] = str((config_path.parent / raw).resolve())
    for key in _PATH_KEYS & expanded.get('data', {}).keys():
        raw = expanded['data'][key]
        if isinstance(raw, str) and raw and not raw.startswith(('~', '/')):
            expanded['data'][key] = str((config_path.parent / raw).resolve())
    for command, values in expanded.get('defaults', {}).items():
        if not isinstance(values, dict):
            raise ValueError(f"Defaults for {command!r} must be an object")
        for key in _PATH_KEYS & values.keys():
            raw = values[key]
            if isinstance(raw, str) and raw and not raw.startswith(('~', '/')):
                values[key] = str((config_path.parent / raw).resolve())
            elif isinstance(raw, list):
                values[key] = [str((config_path.parent / item).resolve()) if isinstance(item, str) and item and not item.startswith(('~', '/')) else item for item in raw]
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
