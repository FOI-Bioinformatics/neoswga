"""Chemistry and panel constraints survive every optimization stage."""

from dataclasses import replace
from types import SimpleNamespace

import pytest

from neoswga.core.base_optimizer import (
    OptimizationResult,
    OptimizationStatus,
    OptimizerConfig,
    PrimerSetMetrics,
)
from neoswga.core.design_context import design_context_from_params
from neoswga.core.network_optimizer import NetworkOptimizer
from neoswga.core.pool_objective import PoolConstraints, PoolObjective
from neoswga.core.unified_optimizer import _minimize_primer_count

A, B, C = "ACGGACGGACGG", "AGGAGGAGGAGG", "ACACACACACAC"


def test_planning_context_preserves_search_controls():
    context = design_context_from_params({"stage1_objective_width": 64, "max_frontier_refills": 7})
    config = context.optimizer_config(refinement_method="swap")
    assert config.stage1_objective_width == 64
    assert config.max_frontier_refills == 7
    assert context.optimizer_config(stage1_objective_width=12).stage1_objective_width == 12


def test_network_rejects_failed_configured_tm():
    def broken(sequence):
        raise ValueError("invalid chemistry")

    holder = SimpleNamespace(
        conditions=SimpleNamespace(fingerprint=lambda: "bad", calculate_effective_tm=broken),
        _tm_cache={},
    )
    with pytest.raises(ValueError, match="invalid chemistry"):
        NetworkOptimizer._get_primer_tm(holder, A)
    assert holder._tm_cache == {}


def test_invalid_conditions_abort_before_loading_candidates(monkeypatch):
    from neoswga.core import reaction_conditions, unified_optimizer

    def broken():
        raise ValueError("invalid chemistry")

    monkeypatch.setattr(reaction_conditions, "build_reaction_conditions", broken)
    monkeypatch.setattr(
        unified_optimizer,
        "_pool_for_this_run",
        lambda *a: pytest.fail("loaded candidates after chemistry failed"),
    )
    with pytest.raises(ValueError, match="invalid chemistry"):
        unified_optimizer.run_optimization(fg_prefixes=["fg"], fg_seq_lengths=[1000], verbose=False)


class Evaluator:
    conditions = object()
    bg_prefixes = ["bg"]
    bg_seq_lengths = [1000]
    config = OptimizerConfig(max_dimer_bp=3, max_self_dimer_bp=4, verbose=False)

    def __init__(self, singleton_effective=0.3, singleton_density=20):
        self.singleton_effective = singleton_effective
        self.singleton_density = singleton_density
        self.pool_objective = PoolObjective(
            self.compute_metrics, PoolConstraints(min_selectivity_density=10)
        )

    def compute_metrics(self, primers):
        return replace(
            PrimerSetMetrics.empty(),
            fg_coverage=0.99,
            effective_fg_coverage=self.singleton_effective if len(primers) == 1 else 0.95,
            selectivity_density=self.singleton_density if len(primers) == 1 else 20,
        )


def result(opt):
    return OptimizationResult(
        primers=(A, B),
        score=0.5,
        metrics=opt.compute_metrics([A, B]),
        status=OptimizationStatus.SUCCESS,
        optimizer_name="test",
        iterations=1,
    )


@pytest.mark.parametrize("effective,density", [(0.3, 20), (0.95, 2)])
def test_minimization_preserves_effective_coverage_and_limits(effective, density):
    opt = Evaluator(effective, density)
    original = result(opt)
    trimmed = _minimize_primer_count(original, opt, 0.9, False)
    assert trimmed.primers == original.primers


def test_minimization_still_removes_valid_redundancy():
    opt = Evaluator(0.95, 20)
    trimmed = _minimize_primer_count(result(opt), opt, 0.9, False)
    assert len(trimmed.primers) == 1
    assert not opt.pool_objective.violations(trimmed.primers)


def test_refinement_and_minimization_share_the_objective():
    from neoswga.core.panel_refinement import refine_result

    opt = Evaluator()

    def metrics(primers):
        return replace(
            PrimerSetMetrics.empty(),
            fg_coverage=0.99,
            effective_fg_coverage=0.97 if C in primers else 0.5,
            selectivity_density=20,
        )

    opt.compute_metrics = metrics
    opt.pool_objective = PoolObjective(metrics, PoolConstraints(min_selectivity_density=10))
    improved = refine_result(result(opt), opt, [A, B, C])
    assert C in improved.primers
    assert improved.metrics.effective_fg_coverage == 0.97


def test_dispatch_attaches_contract_before_selection_and_refines(monkeypatch):
    from neoswga.core import unified_optimizer as unified

    opt = Evaluator()
    del opt.pool_objective
    opt.name = "test"

    def metrics(primers):
        return replace(
            PrimerSetMetrics.empty(),
            fg_coverage=0.99,
            effective_fg_coverage=0.97 if C in primers else 0.5,
            selectivity_density=20,
        )

    opt.compute_metrics = metrics

    def select(candidates, target_size):
        assert opt.pool_objective.constraints.min_selectivity_density == 10
        return result(opt)

    opt.optimize = select
    monkeypatch.setattr(unified.OptimizerFactory, "create", lambda **kwargs: opt)
    monkeypatch.setattr(
        unified,
        "constraints_from_parameter",
        lambda source: PoolConstraints(min_selectivity_density=10),
    )
    delivered, evaluator = unified._dispatch_optimizer(
        "hybrid",
        None,
        [A, B, C],
        ["fg"],
        [1000],
        ["bg"],
        [1000],
        2,
        opt.config,
        opt.conditions,
        1,
        False,
        None,
        {},
    )
    assert evaluator is opt
    assert C in delivered.primers
    contract = opt.pool_objective
    _minimize_primer_count(delivered, opt, 0.9, False)
    assert opt.pool_objective is contract


def test_refinement_does_not_trade_a_feasible_panel_for_coverage():
    from neoswga.core.panel_refinement import refine_result

    opt = Evaluator()

    def metrics(primers):
        return replace(
            PrimerSetMetrics.empty(),
            fg_coverage=0.99,
            effective_fg_coverage=0.99 if C in primers else 0.9,
            selectivity_density=2 if C in primers else 20,
        )

    opt.compute_metrics = metrics
    opt.pool_objective = PoolObjective(metrics, PoolConstraints(min_selectivity_density=10))
    original = result(opt)
    assert refine_result(original, opt, [A, B, C]).primers == original.primers
