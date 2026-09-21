"""Behavioral checks for the shared panel stages."""

from dataclasses import replace

import pytest

from neoswga.core.base_optimizer import (
    OptimizationResult,
    OptimizationStatus,
    OptimizerConfig,
    PrimerSetMetrics,
)
from neoswga.core.optimization_service import OptimizationRequest, run_panel_search
from neoswga.core.pool_objective import PoolConstraints, PoolObjective
from neoswga.core.unified_optimizer import _select_ensemble_winner

A, B, C, D = "ACGGACGGACGG", "AGGAGGAGGAGG", "ACACACACACAC", "AAAACCCCAAAA"


class TableOptimizer:
    conditions = object()
    bg_prefixes = ["bg"]
    bg_seq_lengths = [1000]
    config = OptimizerConfig(swap_max_evaluations=100, swap_max_seconds=2)

    def __init__(self, table, initial):
        self.table = {frozenset(k): v for k, v in table.items()}
        self.initial = initial
        self.pool_objective = PoolObjective(
            self.compute_metrics, PoolConstraints(min_selectivity_density=10)
        )

    def compute_metrics(self, panel):
        coverage, density = self.table.get(frozenset(panel), (0.1, 20))
        return replace(
            PrimerSetMetrics.empty(),
            fg_coverage=0.99,
            effective_fg_coverage=coverage,
            selectivity_density=density,
            total_bg_sites=5,
        )

    def optimize(self, candidates, target_size, **kwargs):
        return OptimizationResult(
            tuple(self.initial),
            1,
            OptimizationStatus.SUCCESS,
            self.compute_metrics(self.initial),
            1,
            "table",
        )


def test_reduction_reopens_swaps_and_can_enable_another_deletion():
    # AB is the first smaller qualifying set; replacing B with C makes A redundant.
    opt = TableOptimizer(
        {(A, B, D): (0.98, 20), (A, B): (0.91, 20), (A, C): (0.97, 20), (C,): (0.95, 20)}, (A, B, D)
    )
    result = run_panel_search(
        OptimizationRequest(opt, (A, B, C, D), 3, minimize=True, target_coverage=0.9)
    )
    assert result.primers == (C,)
    assert result.metrics.effective_fg_coverage == 0.95
    changed = [s["stage"] for s in result.stage_history if s["changed"]]
    assert changed == ["selection", "reduction", "refinement", "reduction"]
    assert result.to_dict()["stage_history"]


def test_fixed_and_excluded_primers_survive_every_stage():
    opt = TableOptimizer({(A, B): (0.91, 20), (A, C): (0.97, 20), (A,): (0.95, 20)}, (A, B))
    result = run_panel_search(
        OptimizationRequest(
            opt,
            (A, B, C),
            2,
            fixed_primers=(B,),
            excluded_primers=(C,),
            minimize=True,
            target_coverage=0.9,
        )
    )
    assert set(result.primers) == {A, B}


def test_excluded_proposal_is_refused():
    opt = TableOptimizer({}, (A, B))
    with pytest.raises(ValueError, match="candidate/size bounds"):
        run_panel_search(OptimizationRequest(opt, (A, B), 2, excluded_primers=(B,)))


def test_ensemble_shared_objective_prefers_feasible_panel_over_composite_score():
    opt = TableOptimizer({(A,): (0.99, 2), (B,): (0.8, 20)}, (A,))
    high = opt.optimize([], 1)
    low = replace(high, primers=(B,), metrics=opt.compute_metrics([B]), score=-10)
    assert (
        _select_ensemble_winner({"high": high, "feasible": low}, "balanced", opt.pool_objective)
        == "feasible"
    )


def test_configured_expansion_passes_conditions_fixed_primers_and_limits(monkeypatch):
    from neoswga.core.design_context import design_context_from_params
    from neoswga.core.optimizer_factory import OptimizerFactory
    from neoswga.core.primer_expansion import PrimerExpander

    context = design_context_from_params(
        {"dmso_percent": 2, "min_selectivity_density": 10, "swap_max_evaluations": 27}
    )
    captured = {}
    opt = TableOptimizer({(A, B): (0.9, 20)}, (A, B))
    del opt.pool_objective

    def create(**kwargs):
        captured.update(kwargs)
        opt.config = kwargs["config"]
        opt.conditions = kwargs["conditions"]
        return opt

    monkeypatch.setattr(OptimizerFactory, "create", create)
    expander = PrimerExpander(None, ["fg"], [1000], ["bg"], [1000], context=context)
    result = expander._expand_hybrid([B], [A], 1, False)
    assert result["new_primers"] == [B]
    assert captured["conditions"] is context.conditions
    assert captured["config"].swap_max_evaluations == 27
    assert opt.pool_objective.constraints.min_selectivity_density == 10
    assert result["stage_history"]


def test_dominating_set_adapter_counts_new_slots_and_preserves_fixed(monkeypatch):
    import numpy as np

    from neoswga.core.dominating_set_adapter import DominatingSetAdapter

    class Cache:
        def get_positions(self, primer, prefix, strand):
            return np.array([100, 600], dtype=np.int64)

    opt = DominatingSetAdapter(Cache(), ["fg"], [1000], config=OptimizerConfig(verbose=False))
    captured = {}

    def greedy(**kwargs):
        captured.update(kwargs)
        return dict(primers=[A, B], coverage=0.9, covered_regions=2, total_regions=3)

    monkeypatch.setattr(opt._optimizer, "optimize_greedy", greedy)
    monkeypatch.setattr(opt, "compute_metrics", lambda panel: PrimerSetMetrics.empty())
    result = opt.optimize([A, B], target_size=2, fixed_primers=[A])
    assert captured["fixed_primers"] == [A]
    assert captured["max_primers"] == 1
    assert A in result.primers


def test_request_normalizes_candidate_and_exclusion_case():
    opt = TableOptimizer({(A,): (0.9, 20)}, (A,))
    request = OptimizationRequest(opt, (A.lower(), B.lower(), A), 1, excluded_primers=(B.lower(),))
    result = run_panel_search(request)
    assert result.primers == (A,)
    assert request.candidates == (A, B)
