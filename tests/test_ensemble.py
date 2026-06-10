"""Ensemble optimizer tests.

Key regression: the ensemble must pick the winner by `normalized_score`
(comparable across optimizers), NOT raw `score` (per-optimizer scale).
"""

import pytest

from neoswga.core.base_optimizer import (
    OptimizationResult,
    OptimizationStatus,
    PrimerSetMetrics,
)


def _metrics(fg_coverage, selectivity, dimer_risk, gap_gini, tm_range):
    return PrimerSetMetrics(
        fg_coverage=fg_coverage,
        bg_coverage=0.0,
        coverage_uniformity=gap_gini,
        total_fg_sites=100,
        total_bg_sites=1,
        selectivity_ratio=selectivity,
        mean_tm=35.0,
        tm_range=tm_range,
        dimer_risk_score=dimer_risk,
        mean_gap=1000.0,
        max_gap=2000.0,
        gap_gini=gap_gini,
        gap_entropy=1.0,
        strand_alternation_score=0.5,
        strand_coverage_ratio=0.5,
    )


def _result(name, raw_score, metrics):
    return OptimizationResult(
        primers=("ATCG", "GGGG"),
        score=raw_score,
        status=OptimizationStatus.SUCCESS,
        metrics=metrics,
        iterations=1,
        optimizer_name=name,
    )


class _StubOptimizer:
    def __init__(self, name, result):
        self.name = name
        self._result = result

    def optimize(self, candidates, target_size=None, **kwargs):
        return self._result


def test_ensemble_picks_by_normalized_score_not_raw(monkeypatch):
    """Method with high raw score but low normalized score must LOSE."""
    from neoswga.core import unified_optimizer as uo

    # 'a': huge raw score, poor real quality -> low normalized.
    a = _result("a", raw_score=9999.0,
                metrics=_metrics(0.10, 0.0, 1.0, 1.0, (10.0, 40.0)))
    # 'b': tiny raw score, excellent real quality -> high normalized.
    b = _result("b", raw_score=0.01,
                metrics=_metrics(0.97, 100.0, 0.0, 0.0, (34.0, 35.0)))
    stubs = {"a": _StubOptimizer("a", a), "b": _StubOptimizer("b", b)}

    def fake_create(name, **kwargs):
        return stubs[name]

    monkeypatch.setattr(uo.OptimizerFactory, "create", staticmethod(fake_create))

    winner = uo._run_ensemble(
        methods=["a", "b"], cache=object(), candidates=["ATCG", "GGGG"],
        fg_prefixes=["fg"], fg_seq_lengths=[1000],
        bg_prefixes=[], bg_seq_lengths=[],
        target_size=2, config=None, conditions=None,
        application="balanced", seed=1, verbose=False,
    )

    assert winner.optimizer_name == "b"          # normalized winner
    assert winner.score == pytest.approx(0.01)   # carried through
    # comparison table present with one row per method, exactly one selected.
    assert winner.ensemble_comparison is not None
    methods = {r["method"] for r in winner.ensemble_comparison}
    assert methods == {"a", "b"}
    selected = [r for r in winner.ensemble_comparison if r["selected"]]
    assert len(selected) == 1 and selected[0]["method"] == "b"


def test_ensemble_degrades_when_some_methods_fail(monkeypatch):
    from neoswga.core import unified_optimizer as uo

    good = _result("good", 1.0, _metrics(0.9, 50.0, 0.0, 0.1, (34.0, 36.0)))

    def fake_create(name, **kwargs):
        if name == "broken":
            raise RuntimeError("boom")
        return _StubOptimizer(name, good)

    monkeypatch.setattr(uo.OptimizerFactory, "create", staticmethod(fake_create))

    winner = uo._run_ensemble(
        methods=["broken", "good"], cache=object(), candidates=["ATCG"],
        fg_prefixes=["fg"], fg_seq_lengths=[1000],
        bg_prefixes=[], bg_seq_lengths=[],
        target_size=1, config=None, conditions=None,
        application="balanced", seed=None, verbose=False,
    )

    assert winner.optimizer_name == "good"
    rows = {r["method"]: r for r in winner.ensemble_comparison}
    assert rows["broken"]["status"] == "error"
    assert rows["broken"]["selected"] is False
    assert rows["good"]["selected"] is True


def test_ensemble_all_methods_fail_returns_error(monkeypatch):
    from neoswga.core import unified_optimizer as uo

    def fake_create(name, **kwargs):
        raise RuntimeError("nope")

    monkeypatch.setattr(uo.OptimizerFactory, "create", staticmethod(fake_create))

    res = uo._run_ensemble(
        methods=["a", "b"], cache=object(), candidates=["ATCG"],
        fg_prefixes=["fg"], fg_seq_lengths=[1000],
        bg_prefixes=[], bg_seq_lengths=[],
        target_size=1, config=None, conditions=None,
        application="balanced", seed=None, verbose=False,
    )
    assert res.status == OptimizationStatus.ERROR


def test_ensemble_is_in_cli_optimize_choices():
    """The optimize command must accept --optimization-method ensemble."""
    from neoswga.cli_unified import create_parser

    parser = create_parser()
    args = parser.parse_args(
        ["optimize", "-j", "params.json", "--optimization-method", "ensemble"]
    )
    assert args.optimization_method == "ensemble"
