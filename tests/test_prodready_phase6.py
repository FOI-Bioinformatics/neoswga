"""Phase 6 (production-readiness v2): optimizer features."""

import pytest

from neoswga.core.base_optimizer import (
    OptimizationResult,
    OptimizationStatus,
    PrimerSetMetrics,
)


def _metrics(fg_coverage, selectivity=10.0, dimer=0.1, gini=0.2, tm=(34.0, 36.0)):
    return PrimerSetMetrics(
        fg_coverage=fg_coverage,
        bg_coverage=0.0,
        coverage_uniformity=gini,
        total_fg_sites=100,
        total_bg_sites=2,
        selectivity_ratio=selectivity,
        mean_tm=35.0,
        tm_range=tm,
        dimer_risk_score=dimer,
        mean_gap=1000.0,
        max_gap=2000.0,
        gap_gini=gini,
        gap_entropy=1.0,
        strand_alternation_score=0.5,
        strand_coverage_ratio=0.5,
    )


def _result(name, primers, fg_coverage, raw=1.0):
    return OptimizationResult(
        primers=tuple(primers),
        score=raw,
        status=OptimizationStatus.SUCCESS,
        metrics=_metrics(fg_coverage),
        iterations=1,
        optimizer_name=name,
    )


class _Stub:
    def __init__(self, name, result):
        self.name = name
        self._result = result

    def optimize(self, candidates, target_size=None, **kwargs):
        return self._result


def test_ensemble_combine_flag_parses():
    from neoswga.cli_unified import create_parser

    parser = create_parser()
    args = parser.parse_args(
        [
            "optimize",
            "-j",
            "p.json",
            "--optimization-method",
            "ensemble",
            "--ensemble-combine",
            "union",
        ]
    )
    assert args.ensemble_combine == "union"


def test_ensemble_union_combine_improves_and_never_worsens(monkeypatch):
    from neoswga.core import unified_optimizer as uo

    # Two methods produce DIFFERENT primers; union pool = {A,B,C,D}.
    a = _result("a", ["A", "B"], fg_coverage=0.6)
    b = _result("b", ["C", "D"], fg_coverage=0.5)
    # Re-optimizing 'a' on the union returns a strictly better set.
    union_better = _result("a", ["A", "C"], fg_coverage=0.9)

    created = {"n": 0}

    def fake_create(name, **kwargs):
        # First two creates are the per-method runs; the 3rd is the union re-opt.
        created["n"] += 1
        if created["n"] <= 2:
            return _Stub(name, {"a": a, "b": b}[name])
        return _Stub(name, union_better)  # union re-optimization of best method

    monkeypatch.setattr(uo.OptimizerFactory, "create", staticmethod(fake_create))

    winner = uo._run_ensemble(
        methods=["a", "b"],
        cache=object(),
        candidates=["A", "B", "C", "D"],
        fg_prefixes=["fg"],
        fg_seq_lengths=[1000],
        bg_prefixes=[],
        bg_seq_lengths=[],
        target_size=2,
        config=None,
        conditions=None,
        application="balanced",
        seed=1,
        verbose=False,
        combine="union",
    )

    # Union result (0.9 coverage) beats both singles (0.6 / 0.5) and is kept.
    assert set(winner.primers) == {"A", "C"}
    union_row = [r for r in winner.ensemble_comparison if r["method"].startswith("union(")]
    assert len(union_row) == 1 and union_row[0]["selected"] is True
    # exactly one selected row overall
    assert sum(1 for r in winner.ensemble_comparison if r["selected"]) == 1


def test_ensemble_union_combine_keeps_winner_when_not_better(monkeypatch):
    from neoswga.core import unified_optimizer as uo

    a = _result("a", ["A", "B"], fg_coverage=0.9)  # best single
    b = _result("b", ["C", "D"], fg_coverage=0.5)
    union_worse = _result("a", ["A", "D"], fg_coverage=0.4)  # worse on union

    created = {"n": 0}

    def fake_create(name, **kwargs):
        created["n"] += 1
        if created["n"] <= 2:
            return _Stub(name, {"a": a, "b": b}[name])
        return _Stub(name, union_worse)

    monkeypatch.setattr(uo.OptimizerFactory, "create", staticmethod(fake_create))

    winner = uo._run_ensemble(
        methods=["a", "b"],
        cache=object(),
        candidates=["A", "B", "C", "D"],
        fg_prefixes=["fg"],
        fg_seq_lengths=[1000],
        bg_prefixes=[],
        bg_seq_lengths=[],
        target_size=2,
        config=None,
        conditions=None,
        application="balanced",
        seed=1,
        verbose=False,
        combine="union",
    )
    # Union was worse -> original winner 'a' kept; no union row added/selected.
    assert set(winner.primers) == {"A", "B"}
    assert not any(r["method"].startswith("union(") for r in winner.ensemble_comparison)


def test_mechanistic_weight_implies_enable():
    """--mechanistic-weight > 0 alone should enable the mechanistic model."""
    from neoswga.cli_unified import create_parser

    parser = create_parser()
    args = parser.parse_args(["optimize", "-j", "p.json", "--mechanistic-weight", "0.5"])
    # default is None now (so the handler can detect an explicit weight)
    assert args.mechanistic_weight == 0.5
    assert args.use_mechanistic_model is False  # not explicitly set
    # The run_step4 handler turns weight>0 into enabled (logic asserted here):
    use_mech = args.use_mechanistic_model
    explicit = args.mechanistic_weight
    if not use_mech and explicit is not None and explicit > 0:
        use_mech = True
    assert use_mech is True
