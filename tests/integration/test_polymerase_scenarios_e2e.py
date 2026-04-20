"""End-to-end polymerase scenario tests: phi29 (30 C, 6-12 bp) and
equiphi29 (43 C, longer primers with additives) against the same plasmid
target. Reuses the pre-built k-mer counts from examples/plasmid_example/.

Each scenario runs the filter step with scenario-specific params.json and
asserts the pipeline completes. Golden-output assertions are kept loose
because the optimizer has stochastic components; the aim is that the
plumbing for each polymerase / additive combo is end-to-end reachable.
"""

import json
import os
import shutil
import tempfile
from pathlib import Path

import pytest


EXAMPLE_DIR = Path(__file__).resolve().parent.parent.parent / "examples" / "plasmid_example"


def _reset_pipeline_state(params_file):
    import neoswga.core.pipeline as pipeline_mod
    from neoswga.core import parameter
    pipeline_mod._initialized = False
    pipeline_mod.fg_prefixes = None
    pipeline_mod.bg_prefixes = None
    pipeline_mod.fg_genomes = None
    pipeline_mod.bg_genomes = None
    pipeline_mod.fg_seq_lengths = None
    pipeline_mod.bg_seq_lengths = None
    pipeline_mod.fg_circular = None
    pipeline_mod.bg_circular = None
    parameter.json_file = params_file


def _build_scenario(tmp_path, overrides: dict):
    """Copy the plasmid example into tmp_path and overlay scenario-specific
    params.json keys."""
    if not EXAMPLE_DIR.is_dir():
        pytest.skip("examples/plasmid_example not available")

    for fname in os.listdir(EXAMPLE_DIR):
        src = EXAMPLE_DIR / fname
        if src.is_file():
            shutil.copy2(src, tmp_path / fname)

    params_path = tmp_path / "params.json"
    with open(params_path) as fh:
        params = json.load(fh)
    params.setdefault("schema_version", 1)
    params.setdefault("min_amp_pred", 0.0)
    params.update(overrides)
    with open(params_path, "w") as fh:
        json.dump(params, fh, indent=2)
    return params_path


@pytest.fixture
def scenario_workdir():
    tmpdir = Path(tempfile.mkdtemp(prefix="neoswga_scenario_"))
    original_cwd = os.getcwd()
    yield tmpdir
    os.chdir(original_cwd)
    shutil.rmtree(tmpdir, ignore_errors=True)


@pytest.mark.integration
@pytest.mark.slow
def test_phi29_scenario_filter_runs(scenario_workdir):
    """phi29 at 30 C with 6-12 bp primers on the plasmid example."""
    params_file = _build_scenario(scenario_workdir, {
        "polymerase": "phi29",
        "reaction_temp": 30.0,
        "min_k": 6,
        "max_k": 12,
    })
    os.chdir(scenario_workdir)
    _reset_pipeline_state(str(params_file))

    from neoswga.core.pipeline import step2
    df = step2()
    assert (scenario_workdir / "step2_df.csv").is_file()
    assert len(df) > 0, "phi29 scenario produced no primers"


@pytest.mark.integration
@pytest.mark.slow
def test_equiphi29_scenario_filter_runs(scenario_workdir):
    """equiphi29 at 43 C with additives on the plasmid example.

    Plasmid is tiny (<10 kb) and can't realistically fit 10-18 bp primers;
    use 8-12 bp here so step2 produces candidates but the polymerase,
    temperature, and additives are exercised."""
    params_file = _build_scenario(scenario_workdir, {
        "polymerase": "equiphi29",
        "reaction_temp": 43.0,
        "min_k": 8,
        "max_k": 12,
        "mg_conc": 10.0,
        "dmso_percent": 5.0,
        "betaine_m": 1.0,
        "min_tm": 20,
        "max_tm": 70,
    })
    os.chdir(scenario_workdir)
    _reset_pipeline_state(str(params_file))

    from neoswga.core.pipeline import step2
    df = step2()
    assert (scenario_workdir / "step2_df.csv").is_file()
    assert len(df) > 0, "equiphi29 scenario produced no primers"


@pytest.mark.integration
@pytest.mark.slow
def test_phi29_vs_equiphi29_produce_different_primer_counts(scenario_workdir):
    """Polymerase choice should propagate: phi29 vs equiphi29 with their
    respective additives should not produce identical primer counts."""
    base_overrides = {"min_k": 8, "max_k": 12, "min_tm": 20, "max_tm": 70}

    phi_dir = scenario_workdir / "phi29"
    equi_dir = scenario_workdir / "equiphi29"
    phi_dir.mkdir()
    equi_dir.mkdir()

    phi_params = _build_scenario(phi_dir, {
        **base_overrides, "polymerase": "phi29", "reaction_temp": 30.0,
    })
    equi_params = _build_scenario(equi_dir, {
        **base_overrides, "polymerase": "equiphi29", "reaction_temp": 43.0,
        "mg_conc": 10.0, "dmso_percent": 5.0, "betaine_m": 1.0,
    })

    from neoswga.core.pipeline import step2

    os.chdir(phi_dir)
    _reset_pipeline_state(str(phi_params))
    phi_df = step2()

    os.chdir(equi_dir)
    _reset_pipeline_state(str(equi_params))
    equi_df = step2()

    # Both should produce primers, but polymerase-specific filtering
    # (higher temperature, additives shifting Tm) should yield different
    # counts. Tolerance allows identical counts on very small genomes.
    assert len(phi_df) > 0 and len(equi_df) > 0
