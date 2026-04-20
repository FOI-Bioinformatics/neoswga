"""Same seed must yield same primers.

Covers the unified_optimizer path where --seed is forwarded through
`kwargs` into `OptimizerFactory.create()` and down to the underlying
optimizer. Before the Phase 13A fix, MOEA's inner MOEAConfig was built
without a seed, so two identical --seed runs could produce different
Pareto fronts.
"""

import json
import os
import shutil
import tempfile
from pathlib import Path

import pytest

EXAMPLE_DIR = Path(__file__).resolve().parent.parent / "examples" / "plasmid_example"


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


def _build_workdir(tmp_path):
    if not EXAMPLE_DIR.is_dir():
        pytest.skip("plasmid example not available")
    for fname in os.listdir(EXAMPLE_DIR):
        src = EXAMPLE_DIR / fname
        if src.is_file():
            shutil.copy2(src, tmp_path / fname)
    params_path = tmp_path / "params.json"
    with open(params_path) as fh:
        params = json.load(fh)
    params.setdefault("schema_version", 1)
    params.setdefault("min_amp_pred", 0.0)
    params["min_k"] = 8
    params["max_k"] = 10
    with open(params_path, "w") as fh:
        json.dump(params, fh, indent=2)
    return params_path


def _prepare(tmp_path):
    params_file = _build_workdir(tmp_path)
    os.chdir(tmp_path)
    _reset_pipeline_state(str(params_file))
    from neoswga.core.pipeline import step2, step3
    step2()
    step3()
    return params_file


@pytest.mark.integration
@pytest.mark.slow
def test_hybrid_same_seed_same_primers(tmp_path):
    """Two runs of hybrid with seed=42 should produce identical primers.

    Hybrid stage 2 is the historically non-deterministic one; seeding at
    the unified_optimizer entry point sets python random + numpy RNGs,
    which covers the bulk of randomness sources."""
    params_file = _prepare(tmp_path)

    from neoswga.core.unified_optimizer import optimize_step4

    _reset_pipeline_state(str(params_file))
    r1 = optimize_step4(optimization_method="hybrid", verbose=False, seed=42)

    _reset_pipeline_state(str(params_file))
    r2 = optimize_step4(optimization_method="hybrid", verbose=False, seed=42)

    primers1, _, _ = r1
    primers2, _, _ = r2
    # At least the top primer set must match. Accept list-of-sets comparison.
    if primers1 and primers2:
        assert set(primers1[0]) == set(primers2[0]), (
            f"Same seed should produce identical primer sets; "
            f"run1={primers1[0]}, run2={primers2[0]}"
        )


@pytest.mark.integration
@pytest.mark.slow
def test_moea_factory_config_receives_seed():
    """MOEABaseOptimizer wrapper must copy seed into MOEAConfig."""
    # Smoke test via introspection of the wrapper source to catch regression.
    import inspect
    try:
        from neoswga.core import moea_optimizer
    except ImportError:
        pytest.skip("moea_optimizer unavailable (pymoo not installed)")
    src = inspect.getsource(moea_optimizer)
    assert "seed=kwargs.get('seed')" in src, (
        "MOEABaseOptimizer must forward seed into MOEAConfig for reproducibility"
    )


@pytest.mark.integration
@pytest.mark.slow
def test_different_seeds_may_differ(tmp_path):
    """Different seeds may produce different sets — weaker but useful
    sanity. Tolerate identical output on the small plasmid example where
    the problem is nearly deterministic."""
    params_file = _prepare(tmp_path)
    from neoswga.core.unified_optimizer import optimize_step4

    _reset_pipeline_state(str(params_file))
    r1 = optimize_step4(optimization_method="hybrid", verbose=False, seed=1)
    _reset_pipeline_state(str(params_file))
    r7 = optimize_step4(optimization_method="hybrid", verbose=False, seed=7)

    p1, _, _ = r1
    p7, _, _ = r7
    # If deterministic on this example, both produce the same set; that's fine.
    # The assertion is only that neither crashed and both returned non-empty.
    assert p1 and p7
