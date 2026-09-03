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
    for name in ("balanced", "discovery", "clinical", "enrichment", "metagenomics"):
        w = OPTIMIZER_APPLICATION_WEIGHTS[name]
        assert {"tm_weight", "uniformity_weight", "dimer_penalty"} <= set(w)
        for k, v in w.items():
            assert 0.0 <= v <= 1.0


def test_clinical_vs_metagenomics_have_different_weights():
    """Clinical favours selectivity (high dimer penalty, high tm weight);
    metagenomics favours coverage (higher uniformity weight)."""
    c = OPTIMIZER_APPLICATION_WEIGHTS["clinical"]
    m = OPTIMIZER_APPLICATION_WEIGHTS["metagenomics"]

    assert (
        c["dimer_penalty"] > m["dimer_penalty"]
    ), "clinical should weigh dimer_penalty higher than metagenomics"
    assert (
        c["tm_weight"] > m["tm_weight"]
    ), "clinical should weigh tm_weight higher than metagenomics"
    assert (
        m["uniformity_weight"] > c["uniformity_weight"]
    ), "metagenomics should weigh uniformity higher than clinical"


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
                    OptimizationResult,
                    OptimizationStatus,
                    PrimerSetMetrics,
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

    assert (
        "tm_weight" in captured
    ), f"run_optimization did not forward tm_weight; got keys {sorted(captured)}"
    clinical_w = OPTIMIZER_APPLICATION_WEIGHTS["clinical"]
    assert captured["tm_weight"] == clinical_w["tm_weight"]
    assert captured["uniformity_weight"] == clinical_w["uniformity_weight"]
    assert captured["dimer_penalty"] == clinical_w["dimer_penalty"]


def test_caller_overrides_application_weights(monkeypatch):
    """Explicit tm_weight in kwargs must win over the application preset."""
    from neoswga.core import optimizer_factory
    from neoswga.core import unified_optimizer as _uo

    captured: dict = {}

    def spy(**kwargs):
        captured.update(kwargs)

        class _DummyOpt:
            name = "dummy"
            conditions = kwargs.get("conditions")

            def optimize(self, candidates, target_size):
                from neoswga.core.base_optimizer import (
                    OptimizationResult,
                    OptimizationStatus,
                    PrimerSetMetrics,
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

    assert (
        captured["tm_weight"] == 0.99
    ), "explicit tm_weight kwarg must override the application preset"


# ----------------------------------------------------------------------
# The profile has to survive a caller that always forwards the weight.
#
# `optimize_step4` and the CLI forward `uniformity_weight` on every call, so
# the key was always present in kwargs and `setdefault` could never install
# the profile value. From the CLI, metagenomics (profile 0.25) and discovery
# (0.20) both ran at 0.0 while the log claimed the profile was applied.
# ----------------------------------------------------------------------


def _capture_factory_kwargs(monkeypatch) -> dict:
    """Stand in for the optimizer factory and record the kwargs it receives."""
    from neoswga.core import optimizer_factory
    from neoswga.core import unified_optimizer as _uo

    captured: dict = {}

    def spy(**kwargs):
        captured.update(kwargs)

        class _DummyOpt:
            name = "dummy"
            conditions = kwargs.get("conditions")

            def optimize(self, candidates, target_size):
                from neoswga.core.base_optimizer import (
                    OptimizationResult,
                    OptimizationStatus,
                    PrimerSetMetrics,
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

    return captured


def _run(monkeypatch, **kwargs) -> dict:
    from neoswga.core import unified_optimizer as _uo

    captured = _capture_factory_kwargs(monkeypatch)
    _uo.run_optimization(
        method="hybrid",
        candidates=["ACGTACGT", "CCCCGGGG"],
        fg_prefixes=["x"],
        fg_seq_lengths=[1000],
        bg_prefixes=[],
        bg_seq_lengths=[],
        target_size=2,
        verbose=kwargs.pop("verbose", False),
        **kwargs,
    )
    return captured


def test_metagenomics_uniformity_weight_reaches_the_optimizer(monkeypatch):
    """`--application metagenomics` with no weight flag: the CLI forwards the
    unset sentinel, and the profile's 0.25 must reach the optimizer."""
    captured = _run(monkeypatch, application="metagenomics", uniformity_weight=None)

    assert captured["uniformity_weight"] == pytest.approx(
        OPTIMIZER_APPLICATION_WEIGHTS["metagenomics"]["uniformity_weight"]
    )


def test_explicit_uniformity_weight_beats_the_profile(monkeypatch):
    captured = _run(monkeypatch, application="metagenomics", uniformity_weight=0.42)

    assert captured["uniformity_weight"] == pytest.approx(0.42)


def test_explicit_zero_uniformity_weight_is_a_choice(monkeypatch):
    """0.0 is a value a user can ask for. Only the sentinel means unset, so a
    fix that special-cases 0.0 would silently reinstate the profile here."""
    captured = _run(monkeypatch, application="metagenomics", uniformity_weight=0.0)

    assert captured["uniformity_weight"] == pytest.approx(0.0)


def test_sentinel_without_an_application_leaves_the_weight_off(monkeypatch):
    """No profile to fall back on: the sentinel must not reach the optimizer,
    where `uniformity_weight > 0` would raise on None."""
    captured = _run(monkeypatch, uniformity_weight=None)

    assert captured.get("uniformity_weight", 0.0) == pytest.approx(0.0)


def test_optimize_step4_does_not_overwrite_the_sentinel(monkeypatch):
    """The layer that caused the defect: optimize_step4 forwarded its own
    signature default on every call."""
    from neoswga.core import unified_optimizer as _uo

    captured: dict = {}

    def spy(**kwargs):
        captured.update(kwargs)
        from neoswga.core.base_optimizer import OptimizationResult

        # A failure result short-circuits the CSV writing below, which needs a
        # populated `parameter` module this test has no use for.
        return OptimizationResult.failure("hybrid", "stopped after capturing kwargs")

    monkeypatch.setattr(_uo, "run_optimization", spy)

    _uo.optimize_step4(optimization_method="hybrid", verbose=False, application="metagenomics")

    assert (
        captured.get("uniformity_weight") is None
    ), "optimize_step4 must forward 'unset' rather than a concrete 0.0"


def test_cli_uniformity_weight_defaults_to_the_unset_sentinel():
    """argparse is where the distinction between "not passed" and "passed 0.0"
    is made or lost."""
    from neoswga.cli_unified import create_parser

    args = create_parser().parse_args(["optimize", "-j", "test.json"])

    assert args.uniformity_weight is None


def test_application_log_reports_the_applied_weights(monkeypatch, caplog):
    """The log line named the profile, not what was applied, so an override
    was reported as the profile value."""
    import logging

    with caplog.at_level(logging.INFO, logger="neoswga.core.unified_optimizer"):
        _run(monkeypatch, application="metagenomics", uniformity_weight=0.42, verbose=True)

    assert "uniformity=0.42" in caplog.text
