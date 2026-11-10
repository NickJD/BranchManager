#!/usr/bin/env python3
"""Regenerate the public BranchManager dashboard demo from committed fixtures.

The TSV/JSON fixtures are deliberately synthetic and are committed under
docs/demo/. This script does not run the biological assessment pipeline; it
passes those fixtures to the same reporting function used by BranchManager.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parents[1]
DEMO_ROOT = REPO_ROOT / "docs" / "demo"

REQUIRED_INPUTS = (
    "assessment/selection_summary.tsv",
    "assessment/sequence_assessment.tsv",
    "assessment/pangenome_hierarchy.tsv",
    "assessment/novelty_metrics.tsv",
    "assessment/genome_collection_assessment.tsv",
    "run_manifest.json",
)
REQUIRED_TAXONOMY = (
    "taxonomy/taxonomy_GTDB_r232.tsv",
    "taxonomy/taxonomy_GG2_2024_09.tsv",
    "taxonomy/taxonomy_NCBI_16S_2026_08.tsv",
    "taxonomy/taxonomy_SILVA_144.tsv",
)

REQUIRED_HTML_IDS = (
    "overview",
    "quality",
    "taxonomy",
    "novelty",
    "genomes",
    "coverage",
    "phylogeny",
    "explorer",
    "isolateBody",
    "exportBtn",
    "toTop",
    "detail",
    "pangenomeDetail",
    "dashboard-data",
)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def _git_revision() -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=REPO_ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        return result.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def _check_inputs(demo_root: Path) -> None:
    missing = [
        str(Path(relative))
        for relative in (*REQUIRED_INPUTS, *REQUIRED_TAXONOMY)
        if not (demo_root / relative).is_file()
    ]
    if missing:
        raise FileNotFoundError(
            "Dashboard demo fixtures are incomplete. Missing: " + ", ".join(missing)
        )


def _prepare_output_assets(demo_root: Path, output_dir: Path) -> None:
    """Copy the static fixture assets when building into a temporary directory."""
    output_dir.mkdir(parents=True, exist_ok=True)
    for directory in ("assessment", "taxonomy"):
        source = demo_root / directory
        target = output_dir / directory
        if source.resolve() == target.resolve():
            continue
        if target.exists():
            shutil.rmtree(target)
        shutil.copytree(source, target)


def _validate_dashboard(html_path: Path) -> dict[str, int]:
    if not html_path.is_file() or html_path.stat().st_size < 20_000:
        raise RuntimeError(f"Dashboard HTML was not generated correctly: {html_path}")

    html = html_path.read_text(encoding="utf-8")
    missing_ids = [element_id for element_id in REQUIRED_HTML_IDS if f'id="{element_id}"' not in html]
    if missing_ids:
        raise RuntimeError(
            "Generated dashboard is missing expected sections or controls: "
            + ", ".join(missing_ids)
        )

    marker = '<script id="dashboard-data" type="application/json">'
    if marker not in html:
        raise RuntimeError("Dashboard does not contain its embedded data payload.")

    payload_start = html.index(marker) + len(marker)
    payload_end = html.find("</script>", payload_start)
    if payload_end < 0:
        raise RuntimeError("Dashboard data payload is not closed.")

    try:
        payload = json.loads(html[payload_start:payload_end])
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"Dashboard data payload is invalid JSON: {exc}") from exc

    records = payload.get("records", [])
    if not isinstance(records, list) or not records:
        raise RuntimeError("Dashboard data payload contains no isolate records.")

    return {"records": len(records), "html_bytes": html_path.stat().st_size}


def build_dashboard(
    demo_root: Path = DEMO_ROOT,
    output_dir: Path | None = None,
) -> Path:
    """Build and validate the demo dashboard; return the generated HTML path."""
    demo_root = demo_root.resolve()
    output_dir = (output_dir or demo_root).resolve()
    _check_inputs(demo_root)
    _prepare_output_assets(demo_root, output_dir)

    # Import only after the caller has configured the package environment.
    try:
        from branchmanager import __version__
        from branchmanager.reporting import write_performance_review_dashboard
    except ImportError as exc:
        raise RuntimeError(
            "BranchManager is not importable. Install it with `python -m pip install -e .` "
            "or set PYTHONPATH=src."
        ) from exc

    manifest_source = json.loads((demo_root / "run_manifest.json").read_text(encoding="utf-8"))
    timestamp = _utc_now()
    manifest = dict(manifest_source)
    manifest.update(
        {
            "workflow": "performance_review",
            "status": "DEMONSTRATION DATA",
            "started_at": timestamp,
            "completed_at": timestamp,
            "branchmanager_version": __version__,
            "git_revision": _git_revision(),
            "notes": (
                "All records and values are synthetic. The dashboard was regenerated "
                "from the committed demonstration fixtures using the current "
                "BranchManager reporting function."
            ),
        }
    )
    manifest_path = output_dir / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    assessment = demo_root / "assessment"
    html_path = Path(
        write_performance_review_dashboard(
            assessment / "selection_summary.tsv",
            output_dir,
            sequence_assessment=assessment / "sequence_assessment.tsv",
            pangenome_hierarchy=assessment / "pangenome_hierarchy.tsv",
            novelty_metrics=assessment / "novelty_metrics.tsv",
            run_manifest=manifest_path,
            genome_assessment=assessment / "genome_collection_assessment.tsv",
        )
    )
    stats = _validate_dashboard(html_path)
    print(
        f"Dashboard generated: {html_path}\n"
        f"Synthetic isolate records: {stats['records']}\n"
        f"HTML size: {stats['html_bytes']:,} bytes"
    )
    return html_path


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEMO_ROOT,
        help="Output directory (default: docs/demo).",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Build in a temporary directory and validate without changing docs/demo.",
    )
    args = parser.parse_args(argv)

    if args.check:
        with tempfile.TemporaryDirectory(prefix="branchmanager-demo-") as temporary:
            build_dashboard(DEMO_ROOT, Path(temporary) / "demo")
    else:
        build_dashboard(DEMO_ROOT, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
