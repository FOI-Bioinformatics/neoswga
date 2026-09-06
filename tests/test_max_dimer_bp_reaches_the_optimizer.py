"""The dimer limit the user configures is the one the optimizer screens on.

Finding A2. max_dimer_bp reached dimer.is_dimer and nothing else. The screen ran
on a hardcoded free-energy threshold, so the pipeline computed one criterion for
tens of minutes and reported the delivered pool against another. Measured
consequence: the shipped E. coli panel carries an 11 bp heterodimer against a
configured 3, and the fresh Wolbachia panel carries 10 bp.
"""

import pytest

from neoswga.core import hybrid_optimizer, parameter


def _optimizer(**kwargs):
    # position_cache is positional and first; None is enough for these tests.
    return hybrid_optimizer.HybridOptimizer(
        None, fg_prefixes=["x"], fg_seq_lengths=[10000], polymerase="equiphi29", **kwargs
    )


def test_constructor_argument_wins(monkeypatch):
    monkeypatch.setattr(parameter, "max_dimer_bp", 3, raising=False)
    assert _optimizer(max_dimer_bp=6).max_dimer_bp == 6


def test_params_value_is_used_when_no_argument(monkeypatch):
    monkeypatch.setattr(parameter, "max_dimer_bp", 5, raising=False)
    assert _optimizer().max_dimer_bp == 5


def test_default_is_three_when_nothing_is_configured(monkeypatch):
    monkeypatch.delattr(parameter, "max_dimer_bp", raising=False)
    assert _optimizer().max_dimer_bp == 3


def test_screen_receives_the_configured_threshold(monkeypatch):
    seen = {}

    class _Filter:
        def __init__(self, criteria):
            self.criteria = criteria

        def filter_candidates(self, candidates, **kwargs):
            seen.update(kwargs)
            return list(candidates), {}

    monkeypatch.setattr(
        "neoswga.core.thermodynamic_filter.ThermodynamicFilter", _Filter, raising=True
    )
    _optimizer(max_dimer_bp=4)._thermo_filter_candidates(["AAACCCGGGTTT", "ACCCGGGTTTAA"])
    assert seen.get("max_dimer_bp") == 4
