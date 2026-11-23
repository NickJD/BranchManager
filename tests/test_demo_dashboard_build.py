from pathlib import Path
import json

from scripts.build_demo_dashboard import DEMO_ROOT, build_dashboard


def test_demo_dashboard_builds_from_committed_fixtures(tmp_path):
    output_dir = tmp_path / "demo"

    html_path = build_dashboard(DEMO_ROOT, output_dir)

    assert html_path.is_file()
    assert html_path.name == "performance_review_dashboard.html"
    html = html_path.read_text(encoding="utf-8")
    for element_id in (
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
    ):
        assert f'id="{element_id}"' in html

    marker = '<script id="dashboard-data" type="application/json">'
    payload_start = html.index(marker) + len(marker)
    payload_end = html.index("</script>", payload_start)
    payload = json.loads(html[payload_start:payload_end])
    assert payload["records"]
    assert len(payload["explorer_records"]) == 60
    assert payload["counts"]["registered_total"] == 60
    assert payload["counts"]["partners"] == {"North Farm": 15, "Rumen Research Unit": 15, "Field Station": 15, "Culture Collection": 15}
    catalogue_only = [row for row in payload["explorer_records"] if row["assessment_present"] is False]
    assert len(catalogue_only) == 36
    assert any(row["already_selected_for_genome_sequencing"] for row in catalogue_only)
    catalogue_members = [member for group in payload["pangenome"] for member in group["members"] if member["kind"] == "Registered catalogue isolate"]
    assert catalogue_members
    assert any(member["committed"] for member in catalogue_members)
    assert (output_dir / "taxonomy" / "taxonomy_GTDB_r232.tsv").is_file()
    assert (output_dir / "assessment" / "selection_summary.tsv").is_file()
    assert (output_dir / "assessment" / "registered_catalogue.tsv").is_file()
    assert (output_dir / "run_manifest.json").is_file()
    assert "Registered catalogue fixture TSV" in html


def test_demo_root_resolves_to_repository_docs_demo():
    assert DEMO_ROOT == Path(__file__).resolve().parents[1] / "docs" / "demo"
    assert (DEMO_ROOT / "assessment" / "selection_summary.tsv").is_file()
