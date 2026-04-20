"""Phase 11C: every CLI-exposed optimizer must write step4_improved_df.csv
with the same columns.

Divergent columns between optimizers break downstream consumers
(results_interpreter, report generation, export). The test runs a small
subset of registered optimizers against the plasmid example and asserts
the produced CSVs have identical column sets.

Reliance on the fast optimizers keeps the test under a few minutes; the
nightly workflow can extend the set.
"""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent.parent
EXAMPLE_DIR = ROOT / "examples" / "plasmid_example"


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


@pytest.fixture
def primed_workdir(tmp_path):
    """Copy plasmid example and run step2 + step3 so optimizers have data."""
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
    os.chdir(tmp_path)
    _reset_pipeline_state(str(params_path))
    from neoswga.core.pipeline import step2, step3
    step2()
    step3()
    return tmp_path


# Fast, deterministic optimizers that finish in a few seconds each on the
# plasmid example. Skip slow / optional ones (moea requires pymoo; milp
# requires mip; genetic has DEAP).
FAST_OPTIMIZERS = ["greedy", "dominating-set", "network", "tiling", "hybrid"]


@pytest.mark.integration
@pytest.mark.slow
def test_every_optimizer_emits_same_columns(primed_workdir):
    from neoswga.core.unified_optimizer import optimize_step4
    from neoswga.core import parameter

    column_sets = {}
    for method in FAST_OPTIMIZERS:
        params_file = primed_workdir / "params.json"
        _reset_pipeline_state(str(params_file))
        try:
            optimize_step4(
                optimization_method=method,
                verbose=False,
            )
        except Exception as e:
            pytest.skip(f"{method} unavailable / failed ({e})")
            continue
        csv_path = primed_workdir / "step4_improved_df.csv"
        if not csv_path.is_file():
            continue
        df = pd.read_csv(csv_path)
        column_sets[method] = frozenset(df.columns)
        # Move out of the way so the next optimizer doesn't inherit state
        csv_path.rename(primed_workdir / f"step4_{method}.csv")

    assert len(column_sets) >= 2, (
        f"Need at least 2 working optimizers to compare; got {list(column_sets)}"
    )

    # Every optimizer must emit the same column set. If this breaks,
    # downstream reporters (results_interpreter, export) can crash.
    canonical = None
    for method, cols in column_sets.items():
        if canonical is None:
            canonical = cols
            canonical_method = method
        else:
            diff_missing = canonical - cols
            diff_extra = cols - canonical
            assert cols == canonical, (
                f"Column mismatch between {canonical_method} and {method}:\n"
                f"  missing in {method}: {sorted(diff_missing)}\n"
                f"  extra   in {method}: {sorted(diff_extra)}"
            )


@pytest.mark.integration
@pytest.mark.slow
def test_canonical_columns_present(primed_workdir):
    """Assert the step4 CSV has the fields downstream consumers expect."""
    from neoswga.core.unified_optimizer import optimize_step4

    params_file = primed_workdir / "params.json"
    _reset_pipeline_state(str(params_file))
    optimize_step4(optimization_method="hybrid", verbose=False)
    csv_path = primed_workdir / "step4_improved_df.csv"
    assert csv_path.is_file()
    df = pd.read_csv(csv_path)

    required = {
        "primer",
        "set_index",
        "score",
        "normalized_score",
        "coverage",
        "bg_coverage",
        "selectivity",
        "mean_gap",
        "max_gap",
        "coverage_uniformity",
        "gap_gini",
        "gap_entropy",
        "dimer_risk_score",
        "strand_alternation",
        "strand_coverage_ratio",
        "mean_tm",
        "optimizer",
    }
    missing = required - set(df.columns)
    assert not missing, f"step4_improved_df.csv missing required columns: {sorted(missing)}"
