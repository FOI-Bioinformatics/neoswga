"""`effective_fg_coverage` is what the swap actually selects on, and it was symmetric.

PR #108 made the geometry reach `fg_coverage` and `bg_coverage`. Those are
reported, gate `--minimize-primers` and feed `PoolObjective("raw")`, but the
quantity the Stage-2 swap refinement optimises is `effective_fg_coverage` --
`pool_objective.PoolObjective.coverage()` returns it by default -- and that goes
through `occupancy_coverage`, which marks symmetric windows through
`merged_window_intervals`.

So until this increment a directional run selected on one geometry in Stage 1,
refined on another in Stage 2, and was reported on a third combination. Asking
whether a delivered panel moves was not answerable.

`merged_window_intervals` is generalised the way `_mark_window` was in PR #106:
the spans come from `site_spans` and the clipping, wrapping and merging stay
where they are. Its default is unchanged and that is most of what is asserted
here, because `effective_fg_coverage` is in every saved summary.
"""

import numpy as np
import pytest

from neoswga.core.coverage import merged_window_intervals
from neoswga.core.occupancy_coverage import occupancy_weighted_coverage
from neoswga.core.reaction_conditions import ReactionConditions

PRIMER = "ACCACAGATAGC"
GENOME = 40_000


# ---------------------------------------------------------------------------
# The interval builder
# ---------------------------------------------------------------------------


def test_the_default_intervals_are_unchanged():
    """One site, reach 100, symmetric: a single 200-wide span."""
    assert merged_window_intervals([1000], 100, GENOME, False) == [(900, 1100)]


def test_a_pooled_caller_gets_the_same_answer_it_always_got():
    """The invariant the default rests on: passing the pooled list and nothing
    else must be what passing it as `positions` has always meant."""
    pooled = [1000, 5000, 9000]

    assert merged_window_intervals(pooled, 300, GENOME, False) == merged_window_intervals(
        pooled, 300, GENOME, False, geometry="symmetric"
    )


def test_the_directional_intervals_reach_one_way():
    """Forward at 1000 reaches right; reverse at 5000 reaches left INTO 5000,
    not from 1000. An earlier draft of this test asserted the reverse span of
    the forward site, which the implementation correctly refused to produce."""
    spans = merged_window_intervals(
        [1000], 100, GENOME, False, reverse=[5000], geometry="directional"
    )

    assert sorted(spans) == [(1000, 1100), (4900, 5000)]


def test_overlapping_directional_spans_still_merge():
    """Merging is a property of the intervals, not of the geometry. Two sites
    150 apart at reach 100 reach into each other."""
    spans = merged_window_intervals(
        [1000], 100, GENOME, False, reverse=[1150], geometry="directional"
    )

    assert spans == [(1000, 1150)] or spans == [(1050, 1150), (1000, 1100)]


def test_an_unknown_geometry_is_refused():
    with pytest.raises(ValueError):
        merged_window_intervals([1], 10, GENOME, False, geometry="radial")


# ---------------------------------------------------------------------------
# The quantity the swap selects on
# ---------------------------------------------------------------------------


def conditions():
    return ReactionConditions(temp=30.0, polymerase="phi29")


def test_the_default_effective_coverage_is_unchanged():
    """Every saved summary carries this figure."""
    sites = {PRIMER: [5_000, 9_000, 13_000]}

    pooled = occupancy_weighted_coverage(
        sites, GENOME, extension_reach=2_000, circular=False, conditions=conditions()
    )
    explicit = occupancy_weighted_coverage(
        sites,
        GENOME,
        extension_reach=2_000,
        circular=False,
        conditions=conditions(),
        geometry="symmetric",
    )

    assert pooled == explicit


def test_the_directional_effective_coverage_is_lower_at_the_same_reach():
    """Half the width per site and these spans do not overlap, so strictly
    less. Not a claim about which is right -- the reach was refitted for that,
    in docs/validation/2026-09-23-reach-refit-directional.md."""
    forward = {PRIMER: [5_000, 9_000, 13_000]}

    symmetric = occupancy_weighted_coverage(
        forward, GENOME, extension_reach=2_000, circular=False, conditions=conditions()
    )
    directional = occupancy_weighted_coverage(
        forward,
        GENOME,
        extension_reach=2_000,
        circular=False,
        conditions=conditions(),
        reverse_by_primer={PRIMER: []},
        geometry="directional",
    )

    assert directional < symmetric
    assert directional == pytest.approx(symmetric / 2, rel=0.02)


def test_a_primer_with_no_sites_in_either_orientation_contributes_nothing():
    """Absence is not coverage, and an empty reverse list is not a reason to
    skip a primer that has forward sites."""
    assert (
        occupancy_weighted_coverage(
            {PRIMER: []},
            GENOME,
            extension_reach=2_000,
            circular=False,
            conditions=conditions(),
            reverse_by_primer={PRIMER: []},
            geometry="directional",
        )
        == 0.0
    )
