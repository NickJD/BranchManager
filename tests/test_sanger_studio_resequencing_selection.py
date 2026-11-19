from pathlib import Path


def test_studio_exposes_per_primer_read_selection_before_assembly():
    html = (Path(__file__).parents[1] / "src/branchmanager/data/sanger_studio.html").read_text()
    assert "Select the read to retain for each primer" in html
    assert 'id="primerReadChoices"' in html
    assert "Only these selections are eligible for assembly" in html
    assert "Select the assembly or read to retain" in html


def test_paper_trail_excludes_same_primer_overlap_and_candidate_pairs():
    source = (Path(__file__).parents[1] / "src/branchmanager/pipeline/paper_trail.py").read_text()
    assert "Repeated reads of the same primer are resequencing" in source
    assert "itertools.product(*selected_groups)" in source
    assert "exactly one attempt from each included primer" in source


def test_read_qc_table_precedes_per_primer_selection_and_assembly():
    html = (Path(__file__).parents[1] / "src/branchmanager/data/sanger_studio.html").read_text()
    qc = html.index("Read quality and trimming evidence")
    selection = html.index("Select the read to retain for each primer")
    assembly = html.index("Select the assembly or read to retain")
    assert qc < selection < assembly



def test_trace_viewer_reads_trace_evidence_from_top_level_inspection_mapping():
    html = (Path(__file__).parents[1] / "src/branchmanager/data/sanger_studio.html").read_text()
    assert "const trace=(r.trace_review||{})[read.read_id];" in html
    assert "const trace=(read.trace_review||{})[read.read_id];" not in html


def test_primer_choice_updates_candidate_preview_not_only_visibility():
    html = (Path(__file__).parents[1] / "src/branchmanager/data/sanger_studio.html").read_text()
    apply_block = html[html.index("function applyPrimerReadChoices"):html.index("choiceBody?.addEventListener", html.index("function applyPrimerReadChoices"))]
    assert "activateCandidateAlignment(preview)" in apply_block
    assert "No assembly candidate is compatible with the current read selections" in apply_block


def test_repeated_input_basenames_are_given_unique_read_ids():
    source = (Path(__file__).parents[1] / "src/branchmanager/pipeline/paper_trail.py").read_text()
    assert "seen_read_ids: set[str]" in source
    assert "original_record_id in seen_read_ids" in source
    assert "candidate_id = f'{original_record_id}__{source_tag}'" in source



def test_active_candidate_metrics_replace_stale_default_assembly_metrics():
    html = (Path(__file__).parents[1] / "src/branchmanager/data/sanger_studio.html").read_text()
    assert "metricValues[0].textContent=activeConsensus.length.toLocaleString()" in html
    assert "metricValues[2].textContent=String(activeDisagreements.length)" in html
    assert "No assembly candidate is compatible with the current read selections." in html


def test_primer_selection_supports_explicit_deselection_and_single_active_read():
    html = (Path(__file__).parents[1] / "src/branchmanager/data/sanger_studio.html").read_text()
    assert 'type="checkbox" name="primer-choice-' in html
    assert 'selectedReadByPrimer.delete(primer)' in html
    assert 'only one read can be active per primer' in html.lower()
    assert 'selectedCandidate=null;' in html
    assert 'selectedReadByPrimer.size' in html


def test_candidate_taxonomy_is_explicitly_run_on_selected_candidate_sequence():
    source = (Path(__file__).parents[1] / "src/branchmanager/sanger_studio.py").read_text()
    html = (Path(__file__).parents[1] / "src/branchmanager/data/sanger_studio.html").read_text()
    assert "def assign_candidate_taxonomy(self, data: dict)" in source
    assert "paper_trail._screen_primer_read_taxonomy(" in source
    assert "'/api/assign-candidate-taxonomy'" in html
    assert "component-read assignments are not inherited" in html


def test_chromatogram_hover_reports_channel_support_and_peak_quality():
    html = (Path(__file__).parents[1] / "src/branchmanager/data/sanger_studio.html").read_text()
    assert "function attachTraceHover(svg, series, rangeForX)" in html
    assert "dominant channel" in html
    assert "secondary_peak_ratios" in html
    assert "trace-hover-readout" in html


def test_read_taxonomy_is_presented_as_per_read_expandable_evidence():
    html = (Path(__file__).parents[1] / "src/branchmanager/data/sanger_studio.html").read_text()
    assert 'class="read-taxonomy-card"' in html
    assert 'hits.filter(h=>h.ReadID===x.ID)' in html


def test_candidate_taxonomy_rejects_unknown_candidate_and_missing_reference(tmp_path, monkeypatch):
    from branchmanager.sanger_studio import SangerStudio

    studio = SangerStudio.__new__(SangerStudio)
    studio.jobs = tmp_path / "inspections"
    studio.jobs.mkdir()
    studio.screen_ref = ""
    studio.screen_taxa = ""
    studio.last_inspection = {"ss1": {"assembly_candidates": [
        {"CandidateID": "candidate_1", "ReadIDs": "read1", "Sequence": "ACGT"}
    ]}}
    import pytest
    with pytest.raises(ValueError, match="taxonomy reference FASTA"):
        studio.assign_candidate_taxonomy({"sequence_id": "ss1", "candidate_id": "candidate_1"})


def test_candidate_taxonomy_screens_exact_candidate_sequence(tmp_path, monkeypatch):
    import csv
    from branchmanager import sanger_studio
    from branchmanager.sanger_studio import SangerStudio

    reference = tmp_path / "reference.fasta"
    reference.write_text(">ref1 Bacteria;Firmicutes;Example\nACGT\n")
    studio = SangerStudio.__new__(SangerStudio)
    studio.jobs = tmp_path / "inspections"
    studio.jobs.mkdir()
    studio.screen_ref = str(reference)
    studio.screen_taxa = ""
    studio.last_inspection = {"ss1": {"assembly_candidates": [
        {"CandidateID": "candidate_1", "ReadIDs": "read1;read2", "Sequence": "ACGTACGT"}
    ]}}

    observed = {}
    def fake_screen(reads, trimmed_fasta, outdir, *, ref_fasta, taxa_tsv, threads):
        observed["sequence"] = Path(trimmed_fasta).read_text()
        observed["read_id"] = reads[0].read_id
        hit_dir = Path(outdir) / "read_taxonomy_screen"
        hit_dir.mkdir(parents=True, exist_ok=True)
        with (hit_dir / "read_taxonomy_hits.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["ReadID", "TargetID", "Taxonomy", "Identity", "QueryCoverage", "TargetCoverage", "Qualified", "NearBest"], delimiter="\t")
            writer.writeheader()
            writer.writerow({"ReadID": reads[0].read_id, "TargetID": "ref1", "Taxonomy": "Bacteria;Firmicutes;Example", "Identity": "99.0", "QueryCoverage": "100", "TargetCoverage": "100", "Qualified": "yes", "NearBest": "yes"})
        return ({reads[0].read_id: {"taxonomy": "Bacteria;Firmicutes;Example", "assigned_rank": "g", "status": "PROVISIONAL", "support": "SINGLE_NEAR_BEST_REFERENCE_HIT", "best_hit": "ref1", "best_identity": 99.0, "best_query_coverage": 100.0, "best_target_coverage": 100.0, "qualified_hit_count": 1, "near_best_hit_count": 1}}, str(Path(outdir) / "concordance.tsv"))

    monkeypatch.setattr(sanger_studio.paper_trail, "_screen_primer_read_taxonomy", fake_screen)
    result = studio.assign_candidate_taxonomy({"sequence_id": "ss1", "candidate_id": "candidate_1"})
    assert observed["sequence"].startswith(">candidate_")
    assert "ACGTACGT" in observed["sequence"]
    assert observed["read_id"].startswith("candidate_")
    assert result["candidate_type"] == "ASSEMBLY"
    assert result["taxonomy"] == "Bacteria;Firmicutes;Example"
    assert result["hits"][0]["TargetID"] == "ref1"
    assert (Path(result["run_directory"]) / "candidate_taxonomy_result.json").is_file()
    assert len(result["sequence_sha256"]) == 64


def test_no_candidate_disables_sequence_download_instead_of_falling_back_to_stale_consensus():
    html = (Path(__file__).parents[1] / "src/branchmanager/data/sanger_studio.html").read_text()
    assert 'id="downloadConsensus" disabled' in html
    assert 'const seq=activeConsensus;if(!seq)return;' in html
    assert "$('#downloadConsensus').disabled=true;" in html


def test_incompatible_candidate_cards_are_hidden_by_author_css():
    html = (Path(__file__).parents[1] / "src/branchmanager/data/sanger_studio.html").read_text()
    assert ".candidate-card[hidden]{display:none!important}" in html
    assert "ids.every((id,i)=>selectedSet.has(id)&&!!primers[i])" in html
    assert "new Set(primers).size===primers.length" in html


def test_selection_change_clears_or_rebuilds_active_candidate_from_compatible_set():
    html = (Path(__file__).parents[1] / "src/branchmanager/data/sanger_studio.html").read_text()
    block = html[html.index("function applyPrimerReadChoices"):html.index("choiceBody?.addEventListener", html.index("function applyPrimerReadChoices"))]
    assert "const compatibleCandidates=[]" in block
    assert "const preview=compatibleCandidates[0]||null" in block
    assert "No compatible candidate is available for chromatogram projection." in block
    assert "$('#downloadConsensus').disabled=true" in block


def test_chromatogram_hover_uses_pointer_events_and_reports_errors_in_readout():
    html = (Path(__file__).parents[1] / "src/branchmanager/data/sanger_studio.html").read_text()
    block = html[html.index("function attachTraceHover"):html.index("function renderTraceReview", html.index("function attachTraceHover"))]
    assert "addEventListener('pointermove',inspect)" in block
    assert "addEventListener('mousemove',inspect)" in block
    assert "getBoundingClientRect()" in block
    assert "Unable to inspect this trace position" in block
    assert "plotFraction" in block
