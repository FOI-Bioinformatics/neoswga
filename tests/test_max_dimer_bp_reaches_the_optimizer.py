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
    assert "unscreened against the already-selected set" in caplog.text
    assert "above max_dimer_bp=3" in caplog.text


def test_relaxation_warning_wording_when_nothing_was_admitted(caplog):
    """The other branch of the same wording fix: lifting the constraint does
    not by itself admit anything, because a candidate must still add
    coverage, and the warning must say so rather than claim a dimerising pair
    was delivered.

    This calls `_log_dimer_relaxation_outcome` directly rather than driving it
    through a live `optimize_greedy` run, because the branch is unreachable
    there by construction, not merely unobserved: `graph.regions` is built
    entirely from the candidates' own binding positions, so no region in it is
    orphaned -- every region exists because some candidate covers it. Combined
    with the full-coverage break at the top of the selection loop, the case is
    closed. If coverage is complete, that break fires before any scan, so
    there is no stall and nothing to relax. If coverage is incomplete, some
    region is uncovered, and by construction some not-yet-selected candidate
    covers it: either it was skipped for dimer reasons, in which case relaxing
    admits it and it scores on retry, or it was not skipped, in which case
    there was no dimer stall to relax in the first place. Either way, whenever
    the stall condition holds, the retry has something to admit. (Attempting a
    live repro confirms this: a redundant, dimerising second primer at the
    same binding position as the first never logs anything, because the
    "Full coverage achieved" branch fires first, before selection is even
    attempted.) The round-1 report of this branch firing live was almost
    certainly finding 1 itself: under the old code, the relaxation's
    `continue` could exhaust the iteration budget before the retry ever ran,
    which presents identically to "nothing admitted" without being that
    branch. The branch stays as real, reachable code regardless
    (`relax_dimer_constraint_when_stuck` could be set False mid-run by a
    subclass or a future caller, or `graph.regions` could stop being built
    this way), so its wording is pinned directly here.
    """
    optimizer = _dimer_stall_optimizer()

    with caplog.at_level(logging.WARNING, logger="neoswga.core.dominating_set_optimizer"):
        optimizer._log_dimer_relaxation_outcome(None, n_selected=2, requested=4)

    assert "no further candidate could add coverage" in caplog.text.lower()
    assert "2 primers selected" in caplog.text
    assert "4 primers requested" in caplog.text
    assert "TTTTGGGGCCCC" not in caplog.text


# Four primers, all four mutually dimerising at max_dimer_bp=3, each with one
# binding site in its own coverage window. Greedy therefore stalls on every
# pick after the first: whichever primer it takes, the whole remainder of the
# pool dimerises with the set. Three unscreened admissions, three warnings.
_MULTI_STALL_POSITIONS = {
    "AAAACCCCGGGG": [10_000],
    "TTTTGGGGCCCC": [50_000],
    "GGGGCCCCTTTT": [90_000],
    "CCCCGGGGAAAA": [130_000],
}


def _multi_stall_optimizer():
    from neoswga.core.dominating_set_optimizer import DominatingSetOptimizer

    return DominatingSetOptimizer(
        cache=_StaticCache(_MULTI_STALL_POSITIONS),
        fg_prefixes=["fg"],
        fg_seq_lengths=[_GENOME],
        bin_size=_REACH // 4,
        extension_reach=_REACH,
        max_dimer_bp=3,
    )


def _violating_pairs(primers, max_dimer_bp):
    """Delivered pairs above the threshold, by the exact pairwise relation."""
    from neoswga.core.dimer import is_dimer_fast

    return [
        (a, b)
        for i, a in enumerate(primers)
        for b in primers[i + 1 :]
        if is_dimer_fast(a, b, max_dimer_bp)
    ]


def _relaxation_warnings(caplog):
    return [r for r in caplog.records if "unscreened against the already-selected set" in r.message]


def test_every_unscreened_admission_is_logged_not_just_the_first(caplog):
    """The relaxation must be re-armed after each admission it makes.

    Setting the local `dimers` to None and clearing only the `relaxed` flag
    turned the dimer screen off for the whole remainder of the run while the
    log still described a single admission. On the 200-primer S. aureus panel
    that delivered 9172 of 19900 pairs above a configured max_dimer_bp of 3
    with five primers named. The three older relaxation tests above cannot see
    it: their pool admits exactly one relaxed pick, so one warning is both the
    correct and the incorrect answer.

    Here the pool stalls on every pick after the first, so the number of
    warnings is the number of primers admitted without a screen -- and the
    property an operator actually needs holds: every delivered pair above the
    threshold involves a primer the log named, so the log is a complete list
    of what to swap by hand.
    """
    optimizer = _multi_stall_optimizer()

    with caplog.at_level(logging.WARNING, logger="neoswga.core.dominating_set_optimizer"):
        result = optimizer.optimize_greedy(
            candidates=list(_MULTI_STALL_POSITIONS), max_primers=4, verbose=False
        )

    delivered = result["primers"]
    assert set(delivered) == set(_MULTI_STALL_POSITIONS)

    warnings = _relaxation_warnings(caplog)
    named = {p for p in delivered if any(p in r.getMessage() for r in warnings)}

    # One warning per primer admitted unscreened: the first pick is screened
    # (nothing is selected yet), the other three are not.
    assert len(warnings) == len(delivered) - 1 == 3
    assert len(named) == len(warnings)

    # And no delivered pair above the threshold is between two primers the log
    # left unnamed. Under the unrestored constraint this fails: only the first
    # relaxed admission was named, so the pairs among the later ones were
    # silent.
    violations = _violating_pairs(delivered, 3)
    assert violations, "the pool must actually deliver violating pairs for this to test anything"
    unnamed_pairs = [(a, b) for a, b in violations if a not in named and b not in named]
    assert unnamed_pairs == []
