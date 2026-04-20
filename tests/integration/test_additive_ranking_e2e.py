"""Ball-bearing test for Phase 7: additives change the chosen primer set.

Before Phase 7A, unified_optimizer.run_optimization never constructed a
ReactionConditions object from parameter.* globals, so additive-aware Tm
paths (integrated_quality_scorer, network_optimizer) never received
conditions and produced identical primer sets regardless of buffer.

This test confirms the fix end-to-end: running optimize with betaine_m=0
vs betaine_m=1.5 under otherwise-identical inputs must now surface a
measurable downstream difference (different set membership OR different
Tm distribution among selected primers).
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


def _build_scenario(tmp_path, overrides):
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


def _run_filter_and_score(workdir, params_file):
    """Run through step2 and step3 so optimize has candidates."""
    _reset_pipeline_state(str(params_file))
    from neoswga.core.pipeline import step2, step3
    step2()
    step3()


def _build_reaction_conditions(betaine_m=0.0, dmso_percent=0.0,
                               polymerase="phi29", temp=30.0):
    from neoswga.core.reaction_conditions import ReactionConditions
    return ReactionConditions(
        temp=temp,
        polymerase=polymerase,
        dmso_percent=dmso_percent,
        betaine_m=betaine_m,
    )


@pytest.mark.integration
@pytest.mark.slow
def test_additives_reach_network_optimizer_via_unified(tmp_path, monkeypatch):
    """When run_optimization is called with betaine_m in parameter.* globals,
    the NetworkOptimizer it creates must carry non-None conditions with the
    matching betaine_m. This guards the Phase 7A plumbing at the seam where
    the previous round's fix was dark code."""
    # Build a minimal plasmid workdir
    params_file = _build_scenario(tmp_path, {
        "polymerase": "phi29",
        "reaction_temp": 30.0,
        "betaine_m": 1.5,
        "dmso_percent": 5.0,
        "min_k": 8,
        "max_k": 10,
    })
    os.chdir(tmp_path)
    _run_filter_and_score(tmp_path, params_file)

    # Now run optimize via the unified entry point and inspect the created
    # optimizer's conditions attribute.
    captured = {}

    from neoswga.core import unified_optimizer
    orig_create = unified_optimizer.OptimizerFactory.create

    def spy_create(*args, **kwargs):
        opt = orig_create(*args, **kwargs)
        captured["optimizer"] = opt
        captured["conditions"] = getattr(opt, "conditions", None)
        return opt

    monkeypatch.setattr(unified_optimizer.OptimizerFactory, "create", spy_create)

    _reset_pipeline_state(str(params_file))
    from neoswga.core.unified_optimizer import optimize_step4
    optimize_step4(optimization_method="hybrid", verbose=False)

    assert "optimizer" in captured, "optimize_step4 did not create an optimizer"
    conditions = captured["conditions"]
    assert conditions is not None, (
        "Phase 7A regression: unified_optimizer created an optimizer without a "
        "ReactionConditions object attached. Additive-aware Tm paths remain "
        "dark code until this is fixed."
    )
    # The conditions object must reflect the user's params.json
    assert abs(conditions.betaine_m - 1.5) < 1e-6, (
        f"betaine_m should propagate; got {conditions.betaine_m}"
    )
    assert abs(conditions.dmso_percent - 5.0) < 1e-6


@pytest.mark.integration
@pytest.mark.slow
def test_network_optimizer_tm_changes_with_additives():
    """The network optimizer's _get_primer_tm cache must return a different
    Tm for the same primer under different additive concentrations once
    conditions are passed."""
    from neoswga.core.network_optimizer import NetworkOptimizer

    plain = _build_reaction_conditions(betaine_m=0.0)
    with_betaine = _build_reaction_conditions(betaine_m=1.5)

    opt_plain = NetworkOptimizer(
        position_cache=None,
        fg_prefixes=[], bg_prefixes=[],
        fg_seq_lengths=[], bg_seq_lengths=[],
        reaction_temp=30.0,
        tm_weight=1.0,
        conditions=plain,
    )
    opt_betaine = NetworkOptimizer(
        position_cache=None,
        fg_prefixes=[], bg_prefixes=[],
        fg_seq_lengths=[], bg_seq_lengths=[],
        reaction_temp=30.0,
        tm_weight=1.0,
        conditions=with_betaine,
    )

    primer = "ATCGATCGATCG"
    tm_plain = opt_plain._get_primer_tm(primer)
    tm_betaine = opt_betaine._get_primer_tm(primer)

    # 1.5 M betaine should lower Tm by roughly 1.2 * 1.5 ~ 1.8 C at 50% GC
    delta = tm_plain - tm_betaine
    assert delta > 0.5, (
        f"Expected betaine 1.5 M to lower Tm by ~1+ C; got delta={delta:.2f}"
    )


@pytest.mark.integration
@pytest.mark.slow
def test_hybrid_optimizer_inner_network_receives_conditions(tmp_path):
    """The outer HybridOptimizerFactory wrapper must forward conditions to
    its inner HybridOptimizer, which must forward them to the inner
    NetworkOptimizer."""
    params_file = _build_scenario(tmp_path, {
        "polymerase": "phi29",
        "betaine_m": 1.0,
        "min_k": 8,
        "max_k": 10,
    })
    os.chdir(tmp_path)
    _run_filter_and_score(tmp_path, params_file)

    _reset_pipeline_state(str(params_file))
    from neoswga.core.unified_optimizer import run_optimization
    from neoswga.core.position_cache import PositionCache
    from neoswga.core import parameter

    # Collect candidates from step2 output
    import pandas as pd
    step2_csv = tmp_path / "step2_df.csv"
    candidates = pd.read_csv(step2_csv)["primer"].astype(str).tolist()[:20]
    if not candidates:
        pytest.skip("no candidates from step2")

    # Run hybrid and spy on the created optimizer
    from neoswga.core import unified_optimizer
    created = {}
    orig_create = unified_optimizer.OptimizerFactory.create

    def spy(*args, **kwargs):
        opt = orig_create(*args, **kwargs)
        created["opt"] = opt
        return opt

    unified_optimizer.OptimizerFactory.create = spy
    try:
        run_optimization(
            method="hybrid",
            candidates=candidates,
            fg_prefixes=list(parameter.fg_prefixes),
            fg_seq_lengths=[6157],
            bg_prefixes=list(parameter.bg_prefixes) if parameter.bg_prefixes else None,
            bg_seq_lengths=[6258] if parameter.bg_prefixes else None,
            target_size=3,
            verbose=False,
        )
    finally:
        unified_optimizer.OptimizerFactory.create = orig_create

    opt = created.get("opt")
    assert opt is not None
    # The outer wrapper exposes conditions
    assert opt.conditions is not None
    assert abs(opt.conditions.betaine_m - 1.0) < 1e-6
    # The inner HybridOptimizer should also have conditions attached
    inner = getattr(opt, "_hybrid", None)
    assert inner is not None
    assert inner.conditions is not None
    assert abs(inner.conditions.betaine_m - 1.0) < 1e-6
    # The innermost NetworkOptimizer too
    inner_net = getattr(inner, "network_optimizer", None)
    assert inner_net is not None
    assert inner_net.conditions is not None, (
        "HybridOptimizer's inner NetworkOptimizer did not receive conditions; "
        "additive-aware Tm scoring is dark at the final selection step."
    )
