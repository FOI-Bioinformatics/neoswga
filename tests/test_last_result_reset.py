"""Phase 16 critical gap #1: _LAST_RESULT must not leak stale data.

Before the fix, `unified_optimizer._LAST_RESULT` was set only on the happy
path. A second `run_optimization` that crashed after the first succeeded
left `_LAST_RESULT` pointing at the first run. CLI consumers
(`--show-frontier`, external audit tools) would silently pick up the
wrong primer set.
"""

import pytest


def test_last_result_starts_none_after_successful_run(monkeypatch):
    """Back-to-back runs must not leak state even when the second raises."""
    from neoswga.core import unified_optimizer as _uo
    from neoswga.core.base_optimizer import (
        OptimizationResult, OptimizationStatus, PrimerSetMetrics,
    )

    class DummyCache:
        def get_positions(self, *a, **kw):
            return []

    # Manually pre-set _LAST_RESULT to simulate a prior successful run
    fake_prior = OptimizationResult(
        primers=("OLD_PRIMER_1", "OLD_PRIMER_2"),
        score=0.99,
        status=OptimizationStatus.SUCCESS,
        metrics=PrimerSetMetrics.empty(),
        iterations=1,
        optimizer_name="prior_run",
    )
    _uo._LAST_RESULT = fake_prior

    # Make the factory raise so run_optimization fails before reaching its
    # post-optimization success path.
    def exploding_create(**kwargs):
        raise RuntimeError("simulated optimizer-construction failure")

    from neoswga.core import optimizer_factory
    monkeypatch.setattr(
        optimizer_factory.OptimizerFactory, "create", exploding_create
    )
    monkeypatch.setattr(_uo, "PositionCache", lambda *a, **kw: DummyCache())

    result = _uo.run_optimization(
        method="hybrid",
        candidates=["ACGTACGT"],
        fg_prefixes=["x"],
        fg_seq_lengths=[1000],
        bg_prefixes=[],
        bg_seq_lengths=[],
        target_size=1,
        verbose=False,
    )

    # The second call returned a failure. _LAST_RESULT must have been
    # reset at the top of run_optimization, not left pointing at the
    # prior success.
    assert result.status == OptimizationStatus.ERROR, (
        f"Expected ERROR from simulated failure, got {result.status}"
    )
    assert _uo._LAST_RESULT is None, (
        f"_LAST_RESULT leaked a stale prior-run result: "
        f"{_uo._LAST_RESULT!r} (was set to fake_prior before the failing call)"
    )


def test_last_result_clears_before_each_run(monkeypatch):
    """Even a clean successful run must clear _LAST_RESULT at entry,
    so external observers cannot see the OLD result while a new run
    is in progress."""
    from neoswga.core import unified_optimizer as _uo
    from neoswga.core.base_optimizer import (
        OptimizationResult, OptimizationStatus, PrimerSetMetrics,
    )

    class DummyCache:
        def get_positions(self, *a, **kw):
            return []

    # Pre-set a prior result
    _uo._LAST_RESULT = OptimizationResult(
        primers=("OLD",),
        score=1.0,
        status=OptimizationStatus.SUCCESS,
        metrics=PrimerSetMetrics.empty(),
        iterations=1,
        optimizer_name="prior",
    )

    observed_during_run = {}

    class _Dummy:
        name = "dummy"
        conditions = None

        def optimize(self, candidates, target_size):
            # Called from inside run_optimization — inspect _LAST_RESULT
            observed_during_run["value"] = _uo._LAST_RESULT
            return OptimizationResult(
                primers=tuple(candidates[: target_size or 1]),
                score=0.5,
                status=OptimizationStatus.SUCCESS,
                metrics=PrimerSetMetrics.empty(),
                iterations=1,
                optimizer_name="dummy",
            )

    from neoswga.core import optimizer_factory
    monkeypatch.setattr(
        optimizer_factory.OptimizerFactory,
        "create",
        lambda **kw: _Dummy(),
    )
    monkeypatch.setattr(_uo, "PositionCache", lambda *a, **kw: DummyCache())

    _uo.run_optimization(
        method="hybrid",
        candidates=["NEW_PRIMER_A", "NEW_PRIMER_B"],
        fg_prefixes=["x"],
        fg_seq_lengths=[1000],
        bg_prefixes=[],
        bg_seq_lengths=[],
        target_size=1,
        verbose=False,
    )

    # Inside optimize(), _LAST_RESULT should have been None — not the
    # prior "OLD" result we seeded.
    assert observed_during_run["value"] is None, (
        f"_LAST_RESULT was {observed_during_run.get('value')!r} inside the "
        f"optimize call; expected None (cleared at run_optimization entry)"
    )
