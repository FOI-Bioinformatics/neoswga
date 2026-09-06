"""The dimer limit the user configures is the one the optimizer screens on.

Finding A2. max_dimer_bp reached dimer.is_dimer and nothing else. The screen ran
on a hardcoded free-energy threshold, so the pipeline computed one criterion for
tens of minutes and reported the delivered pool against another. Measured
consequence: the shipped E. coli panel carries an 11 bp heterodimer against a
configured 3, and the fresh Wolbachia panel carries 10 bp.
"""

import logging

import numpy as np
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


def test_greedy_does_not_select_a_primer_that_dimerises_with_the_set():
    """Finding B6: no default path screened a primer against another primer."""
    from neoswga.core.dimer_matrix import build

    # AAAATTTTAAAA and TTTTAAAATTTT share a long complementary run.
    primers = ["AAAATTTTAAAA", "TTTTAAAATTTT", "CGCGATCGCGAT"]
    matrix = build(primers, 3)
    assert bool(matrix.pairs[0, 1]) is True
    assert bool(matrix.pairs[0, 2]) is False

    from neoswga.core.dominating_set_optimizer import DominatingSetOptimizer

    # The real signature is DominatingSetOptimizer(cache, fg_prefixes,
    # fg_seq_lengths, bin_size=10000, extension_reach=0). There is no
    # bg_prefixes parameter; the cache is positional and first.
    optimizer = DominatingSetOptimizer(
        None, fg_prefixes=["x"], fg_seq_lengths=[10000], max_dimer_bp=3
    )
    assert optimizer._would_dimerise("TTTTAAAATTTT", ["AAAATTTTAAAA"], matrix) is True
    assert optimizer._would_dimerise("CGCGATCGCGAT", ["AAAATTTTAAAA"], matrix) is False


class _StaticCache:
    """Minimal position-cache stand-in: one forward binding site per primer."""

    def __init__(self, positions_by_primer):
        self._by_primer = positions_by_primer

    def get_positions(self, prefix, primer, strand="both"):
        pos = self._by_primer.get(primer, [])
        if strand == "reverse":
            return np.array([], dtype=np.int64)
        return np.asarray(pos, dtype=np.int64)


# Three primers, non-overlapping single-site coverage windows, one dimerising
# pair (AAAACCCCGGGG x TTTTGGGGCCCC, same pair as the unit test above).
# CGCGATCGCGAT is dimer-free against both. A working greedy must eventually
# select all three: TTTTGGGGCCCC is not redundant, it is simply blocked until
# the relaxation admits it.
_GENOME = 1_000_000
_REACH = 3000
_DIMER_STALL_POSITIONS = {
    "AAAACCCCGGGG": [10_000],
    "TTTTGGGGCCCC": [50_000],
    "CGCGATCGCGAT": [90_000],
}


def _dimer_stall_optimizer():
    from neoswga.core.dominating_set_optimizer import DominatingSetOptimizer

    return DominatingSetOptimizer(
        cache=_StaticCache(_DIMER_STALL_POSITIONS),
        fg_prefixes=["fg"],
        fg_seq_lengths=[_GENOME],
        bin_size=_REACH // 4,
        extension_reach=_REACH,
        max_dimer_bp=3,
    )


def test_greedy_falls_back_rather_than_returning_an_undersized_set():
    """If every remaining candidate dimerises, the run must say so and continue
    rather than silently deliver fewer primers than asked for.

    Replaces a vacuous version of this test that asserted only
    `relax_dimer_constraint_when_stuck is True` without ever calling
    `optimize_greedy`, so it could not catch the bug below: the relaxation
    used to consume a loop iteration via `continue` inside
    `for iteration in range(max_primers)`, so a stall landing on the last
    iteration ended the run one primer short -- exactly the failure the
    relaxation exists to prevent. Reproduced before the fix: `max_primers=3`
    on this three-primer pool returned only 2 primers, while `max_primers=4`
    on the same pool correctly returned all 3 -- the difference was purely one
    spare iteration of budget.
    """
    optimizer = _dimer_stall_optimizer()

    result = optimizer.optimize_greedy(
        candidates=list(_DIMER_STALL_POSITIONS), max_primers=3, verbose=False
    )

    assert len(result["primers"]) == 3
    assert set(result["primers"]) == set(_DIMER_STALL_POSITIONS)


def test_relaxation_warning_names_the_primer_it_actually_admitted(caplog):
    """The warning used to fire before the retry resolved, so it promised a
    dimerising pair would be delivered even on a run where nothing ended up
    admitted. It must describe what happened, not what might: here a primer
    genuinely is admitted, so it must be named."""
    optimizer = _dimer_stall_optimizer()

    with caplog.at_level(logging.WARNING, logger="neoswga.core.dominating_set_optimizer"):
        result = optimizer.optimize_greedy(
            candidates=list(_DIMER_STALL_POSITIONS), max_primers=3, verbose=False
        )

    assert len(result["primers"]) == 3
    assert "TTTTGGGGCCCC" in caplog.text
    assert "at least one pair above max_dimer_bp=3" in caplog.text


def test_relaxation_warning_wording_when_nothing_was_admitted(caplog):
    """The other branch of the same wording fix: lifting the constraint does
    not by itself admit anything, because a candidate must still add
    coverage, and the warning must say so rather than claim a dimerising pair
    was delivered.

    This calls `_log_dimer_relaxation_outcome` directly rather than driving it
    through a live `optimize_greedy` run. Attempting the live version first
    (a redundant, dimerising second primer at the same binding position as
    the first) never logged anything: `graph.regions` is defined purely by
    what the candidates contribute, so if a redundant candidate is the only
    thing left, `covered_regions == graph.regions` already holds and the
    "Full coverage achieved" branch fires before selection is even attempted
    -- correctly, since nothing was actually left to gain. A 3,000-trial
    randomised search over synthetic pools (varying primer count, dimer
    pairs, and binding positions) found no live case where the stall fires,
    the graph is not yet fully covered, and the retry still admits nothing:
    whenever `covered_regions != graph.regions`, some not-yet-selected
    candidate accounts for the gap, and relaxation lifts the dimer check for
    every remaining candidate at once, so that candidate always scores
    positively on retry. The report that this branch fired live was against
    the round-1 code, where the *other* bug (the relaxation consuming a loop
    iteration) could end the run before the retry ran at all -- which looked
    like "nothing admitted" but was actually that bug, now fixed above. The
    branch is still real, reachable code (a future caller could set
    `relax_dimer_constraint_when_stuck` False mid-run, or a subclass could
    change what "stuck" means), so its wording is pinned directly here.
    """
    optimizer = _dimer_stall_optimizer()

    with caplog.at_level(logging.WARNING, logger="neoswga.core.dominating_set_optimizer"):
        optimizer._log_dimer_relaxation_outcome(None, n_selected=2, requested=4)

    assert "no further candidate could add coverage" in caplog.text.lower()
    assert "2 primers selected" in caplog.text
    assert "4 primers requested" in caplog.text
    assert "TTTTGGGGCCCC" not in caplog.text
