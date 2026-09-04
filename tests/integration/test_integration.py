"""Scenario definitions for the polymerase integration cases.

WHAT THIS FILE ACTUALLY COVERS. Each scenario directory holds a `params.json`
and nothing else: the genomes its `fg_genomes` key names are not in the
repository. So these tests validate the scenario DEFINITIONS -- that the config
is well formed and internally consistent -- and cannot run the pipeline.

That distinction was previously lost. The module docstring said "integration
tests ... skips gracefully if test genomes are not available", CLAUDE.md
described the three directories as polymerase scenarios under Integration
tests, and `check_genome_available()` was defined and never called. Nothing
skipped, because nothing tried: the two tests were assertions about JSON keys,
and they passed for the same reason an empty test passes.

`test_the_pipeline_runs_when_the_genomes_are_present` closes that gap. It runs
the real pipeline when a scenario's genomes exist and skips with a stated
reason when they do not, so adding the genomes turns these into genuine
integration coverage without editing this file again.

One of the key assertions is also worth reading twice: `optimization_method` in
params.json is documented as inert (Known Issue 8), so asserting its presence
pins a key that changes nothing. It is kept because the schema declares it and
these are schema tests, but it is not evidence that the scenario runs the
method it names.
"""

import json
import pytest
from pathlib import Path

# Integration test scenarios
SCENARIOS = [
    pytest.param("phi29_baseline", id="phi29-standard"),
    pytest.param("equiphi29_baseline", id="equiphi29-standard"),
    pytest.param("phi29_with_bg", id="phi29-with-background"),
]

INTEGRATION_DIR = Path(__file__).parent


def get_scenario_path(scenario: str) -> Path:
    """Get the path to a scenario directory."""
    return INTEGRATION_DIR / scenario


def load_params(scenario: str) -> dict:
    """Load params.json for a scenario."""
    params_file = get_scenario_path(scenario) / "params.json"
    with open(params_file) as f:
        return json.load(f)


def check_genome_available(scenario: str) -> bool:
    """Check if the test genome is available for a scenario."""
    scenario_dir = get_scenario_path(scenario)
    params = load_params(scenario)

    for genome in params.get("fg_genomes", []):
        genome_path = scenario_dir / genome
        if not genome_path.exists():
            return False

    for genome in params.get("bg_genomes", []):
        if genome:  # Skip empty strings
            genome_path = scenario_dir / genome
            if not genome_path.exists():
                return False

    return True


@pytest.mark.integration
@pytest.mark.parametrize("scenario", SCENARIOS)
def test_params_valid(scenario: str):
    """Test that params.json exists and has required fields."""
    params = load_params(scenario)

    # Required fields
    assert "fg_genomes" in params
    assert "polymerase" in params
    assert "min_k" in params
    assert "max_k" in params
    assert "optimization_method" in params

    # Sanity checks
    assert params["min_k"] <= params["max_k"]
    assert params["polymerase"] in ["phi29", "equiphi29", "bst", "klenow"]


@pytest.mark.integration
@pytest.mark.slow
@pytest.mark.parametrize("scenario", SCENARIOS)
def test_scenario_description(scenario: str):
    """Verify each scenario has a description."""
    params = load_params(scenario)
    assert "_description" in params, f"Scenario {scenario} missing _description"


@pytest.mark.integration
@pytest.mark.slow
@pytest.mark.parametrize("scenario", SCENARIOS)
def test_the_pipeline_runs_when_the_genomes_are_present(scenario: str, tmp_path):
    """Run the scenario end to end, or say plainly why it was not run.

    `check_genome_available` existed for this and had no caller, so a scenario
    with no genomes was indistinguishable from one that passed. A skip with a
    reason is the honest outcome; silence is not.
    """
    if not check_genome_available(scenario):
        pytest.skip(
            f"{scenario}: genomes named in params.json are not in the repository, "
            "so this scenario validates its configuration only"
        )

    import os
    import shutil

    from neoswga.core.kmer_counter import check_jellyfish_available

    if not check_jellyfish_available():
        pytest.skip("jellyfish not on PATH")

    src = get_scenario_path(scenario)
    for item in src.iterdir():
        dest = tmp_path / item.name
        shutil.copytree(item, dest) if item.is_dir() else shutil.copy2(item, dest)

    cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        import neoswga.core.pipeline as pipeline_mod
        from neoswga.core import parameter
        from neoswga.core.unified_optimizer import optimize_step4

        pipeline_mod._initialized = False
        parameter.json_file = "params.json"
        pipeline_mod._initialize()
        pipeline_mod.step1()
        pipeline_mod._initialized = False
        assert len(pipeline_mod.step2()) > 0, f"{scenario}: filter produced no candidates"
        pipeline_mod._initialized = False
        assert len(pipeline_mod.step3()) > 0, f"{scenario}: score produced nothing"
        pipeline_mod._initialized = False
        primer_sets, _scores, _cache = optimize_step4(verbose=False, seed=42)
        assert primer_sets and primer_sets[0], f"{scenario}: optimize returned no primers"
    finally:
        os.chdir(cwd)
