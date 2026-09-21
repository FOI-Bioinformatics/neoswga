"""Invalid requested chemistry aborts before any design is attempted."""

import logging

import pytest

from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.unified_optimizer import run_optimization


def test_formamide_out_of_bounds_raises_valueerror():
    """Sanity: ReactionConditions enforces the documented bound so the
    downstream Phase 17D guard has something to catch."""
    with pytest.raises(ValueError, match="Formamide"):
        ReactionConditions(formamide_percent=15.0)


def test_run_optimization_rejects_bad_conditions(monkeypatch):
    from neoswga.core import parameter as param_mod

    monkeypatch.setattr(param_mod, "formamide_percent", 15.0, raising=False)
    monkeypatch.setattr(param_mod, "polymerase", "phi29", raising=False)
    monkeypatch.setattr(param_mod, "reaction_temp", 30.0, raising=False)
    with pytest.raises(ValueError, match="Formamide"):
        run_optimization(
            method="hybrid",
            candidates=[],
            fg_prefixes=["x"],
            fg_seq_lengths=[1000],
            target_size=3,
            verbose=False,
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

    monkeypatch.setattr(param_mod, "formamide_percent", 5.0, raising=False)
    monkeypatch.setattr(param_mod, "polymerase", "phi29", raising=False)
    monkeypatch.setattr(param_mod, "reaction_temp", 30.0, raising=False)

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
        r.getMessage()
        for r in caplog.records
        if r.levelno == logging.ERROR and "ReactionConditions" in r.getMessage()
    ]
    assert (
        not conditions_errors
    ), f"Valid conditions should not generate error logs: {conditions_errors}"
