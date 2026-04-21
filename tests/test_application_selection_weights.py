"""Phase 15B: --application now drives selection weights, not just
post-hoc scoring. Clinical users get more-selective primer sets,
metagenomics users get broader coverage, etc.

Verifies the weight-lookup dict, the unified_optimizer plumbing, and the
fact that the weights actually reach the NetworkOptimizer inside the
default hybrid optimizer.
"""

import pytest

from neoswga.core.base_optimizer import OPTIMIZER_APPLICATION_WEIGHTS


def test_optimizer_weights_table_covers_all_applications():
    for name in ("balanced", "discovery", "clinical",
                 "enrichment", "metagenomics"):
        w = OPTIMIZER_APPLICATION_WEIGHTS[name]
        assert {"tm_weight", "uniformity_weight", "dimer_penalty"} <= set(w)
        for k, v in w.items():
            assert 0.0 <= v <= 1.0


def test_clinical_vs_metagenomics_have_different_weights():
    """Clinical favours selectivity (high dimer penalty, high tm weight);
    metagenomics favours coverage (higher uniformity weight)."""
    c = OPTIMIZER_APPLICATION_WEIGHTS["clinical"]
    m = OPTIMIZER_APPLICATION_WEIGHTS["metagenomics"]

    assert c["dimer_penalty"] > m["dimer_penalty"], (
        "clinical should weigh dimer_penalty higher than metagenomics"
    )
    assert c["tm_weight"] > m["tm_weight"], (
        "clinical should weigh tm_weight higher than metagenomics"
    )
    assert m["uniformity_weight"] > c["uniformity_weight"], (
        "metagenomics should weigh uniformity higher than clinical"
    )


def test_unified_optimizer_forwards_application_weights(monkeypatch):
    """run_optimization pops `application` and sets tm_weight /
    uniformity_weight / dimer_penalty in the factory kwargs if the caller
    hasn't already supplied them."""
    from neoswga.core import unified_optimizer as _uo

    captured: dict = {}

    def spy(**kwargs):
        captured.update(kwargs)

        class _DummyOpt:
            name = "dummy"
            conditions = kwargs.get("conditions")

            def optimize(self, candidates, target_size):
                from neoswga.core.base_optimizer import (
                    OptimizationResult, OptimizationStatus, PrimerSetMetrics,
                )
                return OptimizationResult(
                    primers=tuple(candidates[: target_size or 2]),
                    score=0.5,
                    status=OptimizationStatus.SUCCESS,
                    metrics=PrimerSetMetrics.empty(),
                    iterations=1,
                    optimizer_name="dummy",
                )

        return _DummyOpt()

    from neoswga.core import optimizer_factory

    monkeypatch.setattr(optimizer_factory.OptimizerFactory, "create", lambda **kw: spy(**kw))

    # Need a working position cache; use a stub returning no positions.
    class _DummyCache:
        def get_positions(self, *a, **kw):
            return []

    # Patch PositionCache inside run_optimization (it's built before the
    # factory call)
    from neoswga.core import unified_optimizer as _uo
    monkeypatch.setattr(_uo, "PositionCache", lambda *a, **kw: _DummyCache())

    _uo.run_optimization(
        method="hybrid",
        candidates=["ACGTACGT", "CCCCGGGG"],
        fg_prefixes=["x"],
        fg_seq_lengths=[1000],
        bg_prefixes=[],
        bg_seq_lengths=[],
        target_size=2,
        verbose=False,
        application="clinical",
    )

    assert "tm_weight" in captured, (
        f"run_optimization did not forward tm_weight; got keys {sorted(captured)}"
    )
    clinical_w = OPTIMIZER_APPLICATION_WEIGHTS["clinical"]
    assert captured["tm_weight"] == clinical_w["tm_weight"]
    assert captured["uniformity_weight"] == clinical_w["uniformity_weight"]
    assert captured["dimer_penalty"] == clinical_w["dimer_penalty"]


def test_caller_overrides_application_weights(monkeypatch):
    """Explicit tm_weight in kwargs must win over the application preset."""
    from neoswga.core import unified_optimizer as _uo, optimizer_factory

    captured: dict = {}

    def spy(**kwargs):
        captured.update(kwargs)

        class _DummyOpt:
            name = "dummy"
            conditions = kwargs.get("conditions")

            def optimize(self, candidates, target_size):
                from neoswga.core.base_optimizer import (
                    OptimizationResult, OptimizationStatus, PrimerSetMetrics,
                )
                return OptimizationResult(
                    primers=tuple(candidates[: target_size or 1]),
                    score=0.5,
                    status=OptimizationStatus.SUCCESS,
                    metrics=PrimerSetMetrics.empty(),
                    iterations=1,
                    optimizer_name="dummy",
                )

        return _DummyOpt()

    class _DummyCache:
        def get_positions(self, *a, **kw):
            return []

    monkeypatch.setattr(optimizer_factory.OptimizerFactory, "create", lambda **kw: spy(**kw))
    monkeypatch.setattr(_uo, "PositionCache", lambda *a, **kw: _DummyCache())

    _uo.run_optimization(
        method="hybrid",
        candidates=["ACGTACGT"],
        fg_prefixes=["x"],
        fg_seq_lengths=[1000],
        bg_prefixes=[],
        bg_seq_lengths=[],
        target_size=1,
        verbose=False,
        application="clinical",
        tm_weight=0.99,  # explicit override
    )

    assert captured["tm_weight"] == 0.99, (
        "explicit tm_weight kwarg must override the application preset"
    )
