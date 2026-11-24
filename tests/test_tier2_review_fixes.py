"""Findings of the final review of the tier-2 sweep (2026-11-24)."""

import argparse
import json

import pytest


def test_the_interactive_heatmap_renders_an_unknown_gini():
    """I1: `render_primer_heatmap` computed `1 - p.gini`, which raised once a
    Gini could be None, so `report --level full --interactive` exited 1."""
    pytest.importorskip("plotly")
    from neoswga.core.report.metrics import PrimerMetrics
    from neoswga.core.report.visualizations import render_primer_heatmap

    primers = [
        PrimerMetrics.from_row({"primer": "ACGTACGTAC", "gini": "0.3"}),
        PrimerMetrics.from_row({"primer": "TTGACCATGA", "gini": "nan"}),
    ]
    assert render_primer_heatmap(primers)


@pytest.mark.usefixtures("primed_plasmid_example")
def test_analyze_coverage_without_primers_neither_scans_nor_refuses(tmp_path):
    """I2: with no primers (the BAM-only mode) `analyze-coverage` passed the
    placeholder "AAAAAAAA" into the scanning cache: every genome was scanned
    for it, and without an aligned genome list the command refused naming it
    and wrote design_failure.json, which then blocked export."""
    import shutil
    import subprocess
    import sys
    from pathlib import Path

    example = Path(__file__).resolve().parent.parent / "examples" / "plasmid_example"
    for entry in example.iterdir():
        if entry.is_file():
            shutil.copy2(entry, tmp_path)
    params = json.loads((tmp_path / "params.json").read_text())
    params["bg_genomes"] = []
    (tmp_path / "params.json").write_text(json.dumps(params))
    result = subprocess.run(
        [sys.executable, "-m", "neoswga.cli_unified", "analyze-coverage", "-j", "params.json",
         "-o", "out"],
        cwd=tmp_path, capture_output=True, text=True, timeout=300,
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert "AAAAAAAA" not in output and "Scanning" not in output
    assert not (tmp_path / "design_failure.json").exists()


def test_the_design_context_tm_window_follows_the_polymerase():
    """I3: the request defaulted the Tm window from the polymerase and
    `design_context_from_params` still used 20-50, so plan-pool optimized
    under a window its recorded request did not describe."""
    from neoswga.core.design_context import design_context_from_params
    from neoswga.core.parameter import default_tm_range

    params = {
        "fg_prefixes": ["t"],
        "fg_genomes": ["t.fna"],
        "fg_seq_lengths": [1000],
        "polymerase": "bst",
        "reaction_temp": 63.0,
        "min_k": 15,
        "max_k": 15,
        "num_primers": 6,
    }
    config = design_context_from_params(params).optimizer_config(verbose=False)
    assert (config.min_tm, config.max_tm) == default_tm_range("bst")


def test_the_manifest_records_a_seed_from_the_file(tmp_path):
    """Re-graded review Minor 2: the run was seeded from params.json and the
    manifest recorded seed None."""
    from neoswga.cli._common import manifest_seed

    path = tmp_path / "params.json"
    path.write_text(json.dumps({"seed": 7}))
    assert manifest_seed(argparse.Namespace(seed=None, json_file=str(path)), None) == 7


@pytest.mark.parametrize("value", [42.0, "42", True])
def test_a_seed_that_is_not_an_integer_is_refused(tmp_path, value):
    """Re-graded review Minor 2: a non-integer seed was hashed and silently
    seeded nothing, the defect Task 8 fixed for integers."""
    from neoswga.cli._common import resolve_seed

    path = tmp_path / "params.json"
    path.write_text(json.dumps({"seed": value}))
    with pytest.raises(ValueError, match="seed"):
        resolve_seed(argparse.Namespace(seed=None, json_file=str(path)))


def test_the_set_editing_loader_refuses_an_ambiguity_code():
    """Task 12 ruling, disputed by the review: a primer with N matches nothing
    in the index or an exact scan, so it scored as zero sites."""
    from neoswga.cli.iterate import _load_primer_list

    with pytest.raises(ValueError, match="ACGTNCGTAC"):
        _load_primer_list(["ACGTNCGTAC"], None)
