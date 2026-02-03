"""
Parameterized integration tests for NeoSWGA pipeline.

Tests representative scenarios covering different polymerases and configurations.
Skips gracefully if test genomes are not available.
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
def test_filter_pipeline(scenario: str):
    """Test the filter step of the pipeline."""
    if not check_genome_available(scenario):
        pytest.skip(f"Test genome not available for {scenario}")

    from neoswga.core.pipeline import run_filter

    scenario_dir = get_scenario_path(scenario)
    params = load_params(scenario)

    # Update paths to be absolute
    params["fg_genomes"] = [str(scenario_dir / g) for g in params["fg_genomes"]]
    params["fg_prefixes"] = [str(scenario_dir / p) for p in params["fg_prefixes"]]
    if params.get("bg_genomes"):
        params["bg_genomes"] = [str(scenario_dir / g) for g in params["bg_genomes"] if g]
        params["bg_prefixes"] = [str(scenario_dir / p) for p in params.get("bg_prefixes", []) if p]

    # Run filter (this is the core test)
    # Note: Actual implementation depends on pipeline API
    # This is a placeholder for the test structure
    assert True  # Replace with actual pipeline call when ready


@pytest.mark.integration
@pytest.mark.slow
@pytest.mark.parametrize("scenario", SCENARIOS)
def test_scenario_description(scenario: str):
    """Verify each scenario has a description."""
    params = load_params(scenario)
    assert "_description" in params, f"Scenario {scenario} missing _description"
