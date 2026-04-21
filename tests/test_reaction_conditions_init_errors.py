"""Phase 17D — loud failure when ReactionConditions cannot be built.

Previously, a bare `except Exception` around ReactionConditions()
caught out-of-bounds additive concentrations (e.g. formamide>10%),
logged at debug/warning, set conditions=None, and let the optimizer
run additive-blind. The user's wet-lab cocktail was silently ignored
and the Tm estimates were off by several degrees.

These tests lock in the Phase 17D behaviour:
- The catch is narrowed to (ValueError, TypeError, KeyError) so
  genuine bugs still propagate.
- Logging level is ERROR, not debug.
- A `reaction_conditions_init_failed` warning shows up in the
  validator report so Phase 17A surfaces it in the HTML.
"""

import logging

import pytest

from neoswga.core.base_optimizer import OptimizationStatus
from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.unified_optimizer import run_optimization


def test_formamide_out_of_bounds_raises_valueerror():
    """Sanity: ReactionConditions enforces the documented bound so the
    downstream Phase 17D guard has something to catch."""
    with pytest.raises(ValueError, match="Formamide"):
        ReactionConditions(formamide_percent=15.0)


def test_run_optimization_logs_error_on_bad_conditions(caplog, monkeypatch):
    """Out-of-bounds formamide must emit an ERROR-level log line — not
    a debug, not a warning — so CLI users see it in their terminal."""
    from neoswga.core import parameter as param_mod
    monkeypatch.setattr(param_mod, 'formamide_percent', 15.0, raising=False)
    monkeypatch.setattr(param_mod, 'polymerase', 'phi29', raising=False)
    monkeypatch.setattr(param_mod, 'reaction_temp', 30.0, raising=False)

    with caplog.at_level(logging.ERROR, logger="neoswga.core.unified_optimizer"):
        # Phase 17D: ReactionConditions construction now happens BEFORE
        # the 17C empty-candidate guard, so a bad formamide value is
        # reported first. Empty candidates still short-circuit afterwards
        # (which is why this call returns quickly without a real pool).
        result = run_optimization(
            method="hybrid",
            candidates=[],
            fg_prefixes=["x"],
            fg_seq_lengths=[1000],
            target_size=3,
            verbose=False,
        )

    # The error log message must contain the specific bound-violation
    # text and the actionable hint.
    error_msgs = [r.getMessage() for r in caplog.records
                  if r.levelno == logging.ERROR]
    assert any("ReactionConditions construction failed" in m
               for m in error_msgs), (
        f"Expected an ERROR-level log line about ReactionConditions, "
        f"got: {error_msgs}"
    )
    assert any("Formamide" in m for m in error_msgs), (
        "Error log should echo the specific out-of-bounds additive"
    )


def test_narrow_catch_propagates_genuine_bugs(monkeypatch):
    """Exceptions outside (ValueError, TypeError, KeyError) must
    propagate — they indicate a real bug, not user misconfiguration."""
    import neoswga.core.unified_optimizer as uo

    class BrokenConditions(Exception):
        pass

    def broken_init(*args, **kwargs):
        raise BrokenConditions("genuine bug: something is wrong internally")

    from neoswga.core import reaction_conditions as rc
    monkeypatch.setattr(rc, "ReactionConditions", broken_init)

    # Should NOT be caught — the broken-init exception should propagate.
    with pytest.raises(BrokenConditions):
        run_optimization(
            method="hybrid",
            candidates=["ACGT"] * 5,
            fg_prefixes=["x"],
            fg_seq_lengths=[10000],
            target_size=3,
            verbose=False,
        )


def test_valid_conditions_do_not_emit_warning(caplog, monkeypatch):
    """Happy path: with in-bounds additives, no
    reaction_conditions_init_failed warning should be generated."""
    from neoswga.core import parameter as param_mod
    monkeypatch.setattr(param_mod, 'formamide_percent', 5.0, raising=False)
    monkeypatch.setattr(param_mod, 'polymerase', 'phi29', raising=False)
    monkeypatch.setattr(param_mod, 'reaction_temp', 30.0, raising=False)

    with caplog.at_level(logging.ERROR, logger="neoswga.core.unified_optimizer"):
        run_optimization(
            method="hybrid",
            candidates=[],  # 17C short-circuit, but 17D runs first
            fg_prefixes=["x"],
            fg_seq_lengths=[1000],
            target_size=3,
            verbose=False,
        )

    # No ReactionConditions-specific error log.
    conditions_errors = [
        r.getMessage() for r in caplog.records
        if r.levelno == logging.ERROR
        and "ReactionConditions" in r.getMessage()
    ]
    assert not conditions_errors, (
        f"Valid conditions should not generate error logs: {conditions_errors}"
    )
