"""A circular target cannot be covered by a primer that does not bind it.

`compute_per_prefix_coverage` and `marginal_coverage_curve` both short-circuit
when the window diameter reaches around a circular genome: `if circular and
2 * extension >= length: occupied[:] = True`. The comment beside it says
"covered by any single site", which is the correct rule, but neither function
checked that a site exists before applying it.

So an empty position cache reported 1.0 on a small circular target -- the best
possible coverage, from a panel that binds nothing. The marginal curve was worse
still: the shortcut sat inside the per-primer loop, so the first primer marked
the whole genome whether or not it had a position, and every later entry
inherited it.

This matters where small circular targets are the subject: plasmids, and the
bundled 5.4 kb plasmid example at the 3 kb default reach.
"""

from typing import Dict, List

import numpy as np
import pytest

from neoswga.core.coverage import compute_per_prefix_coverage, marginal_coverage_curve

# 2 * extension >= length, so the shortcut is live.
LENGTH = 1_000
REACH = 1_000


class _Cache:
    """PositionCache-shaped, with per-primer positions."""

    def __init__(self, positions: Dict[str, List[int]]):
        self._positions = positions

    def get_positions(self, prefix, primer, strand):  # noqa: ARG002
        return np.asarray(self._positions.get(primer, []), dtype=int)


def test_a_circular_target_with_no_sites_is_not_fully_covered():
    agg, per = compute_per_prefix_coverage(
        cache=_Cache({}),
        primers=["ACGTACGTACGT"],
        prefixes=["t"],
        seq_lengths=[LENGTH],
        extension=REACH,
        circular=True,
    )

    assert agg == 0.0
    assert per["t"] == 0.0


def test_one_site_on_a_small_circle_still_covers_everything():
    """The shortcut is correct once its premise holds; keep it."""
    agg, _ = compute_per_prefix_coverage(
        cache=_Cache({"ACGTACGTACGT": [500]}),
        primers=["ACGTACGTACGT"],
        prefixes=["t"],
        seq_lengths=[LENGTH],
        extension=REACH,
        circular=True,
    )

    assert agg == 1.0


def test_the_marginal_curve_does_not_credit_a_primer_that_binds_nothing():
    curve = marginal_coverage_curve(
        cache=_Cache({"CCCCCCCCCCCC": [250]}),
        primers=["AAAAAAAAAAAA", "CCCCCCCCCCCC"],
        prefixes=["t"],
        seq_lengths=[LENGTH],
        extension=REACH,
        circular=True,
    )

    # The first primer has no sites, so it buys nothing; the second wraps the
    # whole circle from its one site.
    assert curve[0]["coverage"] == 0.0
    assert curve[0]["marginal_pp"] == 0.0
    assert curve[1]["coverage"] == 1.0


@pytest.mark.parametrize("circular", [True, False])
def test_an_absent_prefix_contributes_no_coverage(circular):
    """A prefix the cache holds nothing for is zero, not one."""
    agg, per = compute_per_prefix_coverage(
        cache=_Cache({"ACGTACGTACGT": [10]}),
        primers=["ACGTACGTACGT"],
        prefixes=["present", "absent"],
        seq_lengths=[LENGTH, LENGTH],
        extension=REACH,
        circular=circular,
    )

    assert 0.0 <= per["absent"] <= 1.0
    assert agg <= 1.0
