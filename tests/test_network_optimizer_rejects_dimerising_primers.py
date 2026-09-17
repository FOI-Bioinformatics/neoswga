"""The `network` method had no dimer rejection at all.

`DominatingSetOptimizer.optimize_greedy` screens candidates against the
selected set with a hard guard. `NetworkOptimizer.optimize_greedy` has its own
selection loop and never reached it. Its only dimer-aware term is a soft
multiplier gated on `dimer_penalty`, which defaults to 0.0, so by default there
was no dimer term at all -- and even at 1.0 a dimerising primer is downweighted
rather than excluded.

This matters because the shipped M. tuberculosis design comes from `network`,
and `ensemble` picks its winner on normalized score alone, so an unguarded
`network` set can win the ensemble and its union re-optimization is unguarded
too.

Deliberately NOT ported from the dominating-set guard: the relaxation that
admits an unscreened primer when the pool is exhausted. Plan 2's sweep measured
that relaxation admitting 171, 129 and 10 unscreened primers at set 0 on the
three shipped designs, which is what produces their 11 bp worst heterodimers.
Here the loop stops instead, and says why.
"""

import numpy as np
import pytest

from neoswga.core import dimer_matrix
from neoswga.core.network_optimizer import NetworkOptimizer

# A dimerising pair: each is the reverse complement of the other, so they carry
# a 12 bp complementary run, far above any usable max_dimer_bp.
PAIR_A = "AAGGTGCGAATA"
PAIR_B = "TATTCGCACCTT"
# Three primers with no long complementary run against the pair or each other.
SAFE = ("ACACACACACAA", "AGAGAGAGAGAA", "ATATATATATAA")


class _Cache:
    """Positions for every primer on one foreground prefix."""

    def __init__(self, mapping):
        self._mapping = mapping

    def get_positions(self, prefix, primer, strand="both"):
        if prefix != "fg":
            return np.array([], dtype=np.int64)
        return np.asarray(self._mapping.get(primer, []), dtype=np.int64)


def _optimizer(primers, max_dimer_bp=3):
    # Spread the sites so the objective can tell candidates apart.
    mapping = {p: [1000 + 4000 * i, 30000 + 4000 * i] for i, p in enumerate(primers)}
    return NetworkOptimizer(
        position_cache=_Cache(mapping),
        fg_prefixes=["fg"],
        bg_prefixes=[],
        fg_seq_lengths=[100_000],
        bg_seq_lengths=[],
        max_dimer_bp=max_dimer_bp,
    )


def test_the_pair_really_does_dimerise_at_the_threshold():
    """Guard the guard: if the fixture pair stopped dimerising, every
    assertion below would pass for the wrong reason."""
    matrix = dimer_matrix.build([PAIR_A, PAIR_B], 3)
    assert matrix.dimerises(PAIR_B, [PAIR_A])


def test_the_greedy_does_not_select_both_members_of_a_dimerising_pair():
    """The regression. Both used to be selectable."""
    primers = [PAIR_A, PAIR_B, *SAFE]
    selected = _optimizer(primers).optimize_greedy(primers, num_primers=5)

    assert not (PAIR_A in selected and PAIR_B in selected), selected


def test_the_safe_primers_are_still_selected():
    """A guard that rejects everything would pass the test above."""
    primers = [PAIR_A, PAIR_B, *SAFE]
    selected = _optimizer(primers).optimize_greedy(primers, num_primers=5)

    assert set(SAFE) <= set(selected), selected


def test_a_loose_threshold_admits_the_pair():
    """The guard reads max_dimer_bp rather than rejecting on some fixed rule."""
    primers = [PAIR_A, PAIR_B]
    # 7 is the loosest value the matrix representation supports.
    selected = _optimizer(primers, max_dimer_bp=7).optimize_greedy(primers, num_primers=2)
    assert len(selected) == 1, "a 12 bp complementary run exceeds 7 as well"


def test_the_loop_stops_rather_than_admitting_an_unscreened_primer(caplog):
    """The dominating-set guard relaxes when the pool is exhausted, which is
    what produces the 11 bp heterodimers in the shipped panels. This one stops.
    """
    import logging

    primers = [PAIR_A, PAIR_B]
    with caplog.at_level(logging.WARNING):
        selected = _optimizer(primers).optimize_greedy(primers, num_primers=2)

    assert len(selected) == 1
    assert "max_dimer_bp" in caplog.text, caplog.text


def test_an_empty_candidate_list_is_not_an_error():
    assert _optimizer([]).optimize_greedy([], num_primers=3) == []


def test_a_threshold_the_matrix_cannot_represent_is_still_enforced():
    """The property this has always protected, asserted directly.

    `dimer_matrix.build` cannot represent a threshold above 7: it codes t-mers
    in a 4**8 space. This used to assert that the greedy raised, which kept it
    from quietly continuing with no screening but also meant a configured
    threshold of 9 could not be used at all.

    Since 2026-09-17 `lazy_dimer.dimer_screen` sends such a threshold to the
    pairwise screen, which has no code-space limit, so it is enforced rather
    than refused. PAIR_A and PAIR_B are exact reverse complements and carry a
    12 bp complementary run, so they exceed 9 and must still be rejected.
    """
    primers = [PAIR_A, PAIR_B, *SAFE]

    selected = _optimizer(primers, max_dimer_bp=9).optimize_greedy(primers, num_primers=3)

    assert not (PAIR_A in selected and PAIR_B in selected), (
        "both halves of a 12 bp complementary pair were selected at max_dimer_bp=9, "
        "so the screen was disabled rather than enforced"
    )


def test_explicit_relaxation_can_admit_the_pair(caplog):
    optimizer = _optimizer([PAIR_A, PAIR_B])
    optimizer.allow_dimer_relaxation = True
    selected = optimizer.optimize_greedy([PAIR_A, PAIR_B], num_primers=2)
    assert set(selected) == {PAIR_A, PAIR_B}
    assert "unscreened against the already-selected set" in caplog.text
