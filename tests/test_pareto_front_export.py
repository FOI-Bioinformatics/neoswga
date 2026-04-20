"""Phase 14A/B: MOEA exports the full Pareto front; CLI can render it.

MOEA's NSGA-III computes a non-dominated set internally but the legacy
code collapsed it to a single best solution before returning. These
tests assert:

- OptimizationResult has a `pareto_front` field populated by MOEA (tuple
  of tuples of primer sequences).
- `pareto_metrics` carries objectives per solution.
- `to_dict` serialises pareto_front when present.
- Other optimizers leave pareto_front None.
"""

import pytest

from neoswga.core.base_optimizer import (
    OptimizationResult, OptimizationStatus, PrimerSetMetrics,
)


def test_optimization_result_accepts_pareto_fields():
    """OptimizationResult should take pareto_front and pareto_metrics."""
    r = OptimizationResult(
        primers=("ATCG", "GCTA"),
        score=0.8,
        status=OptimizationStatus.SUCCESS,
        metrics=PrimerSetMetrics.empty(),
        iterations=1,
        optimizer_name="test",
        pareto_front=(
            ("ATCG", "GCTA"),
            ("ATCG", "TACG"),
        ),
        pareto_metrics=(
            {"target_coverage": 0.9, "n_primers": 2},
            {"target_coverage": 0.85, "n_primers": 2},
        ),
    )
    assert r.pareto_front is not None
    assert len(r.pareto_front) == 2
    assert r.pareto_metrics[0]["target_coverage"] == 0.9


def test_to_dict_includes_pareto_front_when_set():
    r = OptimizationResult(
        primers=("ATCG",),
        score=0.9,
        status=OptimizationStatus.SUCCESS,
        metrics=PrimerSetMetrics.empty(),
        iterations=1,
        optimizer_name="test",
        pareto_front=(("ATCG",), ("GCTA",)),
        pareto_metrics=({"target_coverage": 0.9},) * 2,
    )
    d = r.to_dict()
    assert "pareto_front" in d
    assert d["pareto_front"] == [["ATCG"], ["GCTA"]]
    assert "pareto_metrics" in d


def test_to_dict_excludes_pareto_when_none():
    r = OptimizationResult(
        primers=("ATCG",),
        score=0.9,
        status=OptimizationStatus.SUCCESS,
        metrics=PrimerSetMetrics.empty(),
        iterations=1,
        optimizer_name="hybrid",
    )
    d = r.to_dict()
    assert "pareto_front" not in d


def test_non_moea_optimizers_leave_pareto_front_none():
    """Hybrid / greedy / etc. should not populate pareto_front."""
    r = OptimizationResult.failure("hybrid", "test")
    assert r.pareto_front is None


def test_moea_wrapper_source_populates_pareto_fields():
    """MOEABaseOptimizer.optimize must construct pareto_front and
    pareto_metrics when MOEA has a non-empty pareto_front result."""
    import inspect
    try:
        from neoswga.core import moea_optimizer
    except ImportError:
        pytest.skip("moea_optimizer requires pymoo")

    src = inspect.getsource(moea_optimizer)
    assert "pareto_front=pf_primers" in src, (
        "MOEABaseOptimizer.optimize must pass pareto_front into OptimizationResult"
    )
    assert "pareto_metrics=pf_metrics" in src
