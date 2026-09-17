"""Strict defaults and explicit relaxation reach the greedy selection stages."""

import pytest

from neoswga.core import parameter, unified_optimizer
from neoswga.core.base_optimizer import OptimizerConfig
from neoswga.core.optimizer_factory import OptimizerFactory


@pytest.mark.parametrize("method", ["dominating-set", "hybrid", "background-aware", "network"])
@pytest.mark.parametrize("allow", [False, True])
def test_policy_reaches_inner_optimizers(method, allow):
    unified_optimizer._ensure_optimizers_registered()
    optimizer = OptimizerFactory.create(
        method,
        position_cache=object(),
        fg_prefixes=["fg"],
        fg_seq_lengths=[10000],
        bg_prefixes=["bg"],
        bg_seq_lengths=[10000],
        config=OptimizerConfig(allow_dimer_relaxation=allow),
    )
    if method == "dominating-set":
        assert optimizer._optimizer.relax_dimer_constraint_when_stuck is allow
    elif method == "network":
        assert optimizer._network.allow_dimer_relaxation is allow
    else:
        assert optimizer._hybrid.dominating_optimizer.relax_dimer_constraint_when_stuck is allow
        assert optimizer._hybrid.network_optimizer.allow_dimer_relaxation is allow


def test_default_is_strict():
    assert OptimizerConfig().allow_dimer_relaxation is False


def test_kwarg_overrides_params_policy(monkeypatch):
    monkeypatch.setattr(parameter, "allow_dimer_relaxation", True, raising=False)
    build = unified_optimizer._build_optimizer_config
    assert build(6, False, 3000, False, {}).allow_dimer_relaxation is True
    assert (
        build(6, False, 3000, False, {"allow_dimer_relaxation": False}).allow_dimer_relaxation
        is False
    )


def test_invalid_policy_is_rejected():
    with pytest.raises(ValueError, match="must be a boolean"):
        OptimizerConfig(allow_dimer_relaxation="false").validate()


def test_cli_default_does_not_override_params():
    import argparse

    from neoswga.cli._optimize_parser import _add_optimize_option_groups

    parser = argparse.ArgumentParser()
    _add_optimize_option_groups(parser)
    assert parser.parse_args([]).allow_dimer_relaxation is None
    assert parser.parse_args(["--allow-dimer-relaxation"]).allow_dimer_relaxation is True


def test_strict_selection_rejects_incompatible_fixed_primers():
    from neoswga.core.dominating_set_optimizer import DominatingSetOptimizer

    optimizer = DominatingSetOptimizer(object(), ["fg"], [10000], max_dimer_bp=3)
    assert optimizer.relax_dimer_constraint_when_stuck is False
    with pytest.raises(ValueError, match="Fixed primers"):
        optimizer._build_dimer_matrix_for_greedy([], ["AAGGTGCGAATA", "TATTCGCACCTT"])


def test_unsupported_matrix_threshold_does_not_disable_screen():
    """A threshold the dense matrix cannot hold is enforced, not refused.

    It used to raise, which prevented silent disabling at the cost of making
    the threshold unusable. `lazy_dimer.dimer_screen` now routes it to the
    pairwise screen, which has no 4**8 code-space limit, so what the caller
    configured is what gets applied. The two sequences are exact reverse
    complements, a 12 bp run, so they exceed 8 and must still be flagged.
    """
    from neoswga.core.dominating_set_optimizer import DominatingSetOptimizer

    optimizer = DominatingSetOptimizer(object(), ["fg"], [10000], max_dimer_bp=8)
    screen = optimizer._build_dimer_matrix_for_greedy(["AAGGTGCGAATA", "TATTCGCACCTT"])

    assert screen.dimerises("AAGGTGCGAATA", ["TATTCGCACCTT"])


@pytest.mark.parametrize("method", ["hybrid", "background-aware"])
def test_swap_settings_reach_hybrid(method):
    unified_optimizer._ensure_optimizers_registered()
    optimizer = OptimizerFactory.create(
        method,
        position_cache=object(),
        fg_prefixes=["fg"],
        fg_seq_lengths=[10000],
        config=OptimizerConfig(
            refinement_method="swap", swap_max_evaluations=7, swap_max_seconds=0.5
        ),
    )
    assert optimizer._hybrid.refinement_method == "swap"
    assert optimizer._hybrid.swap_max_evaluations == 7
    assert optimizer._hybrid.swap_max_seconds == 0.5
