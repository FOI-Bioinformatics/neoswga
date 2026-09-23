"""Where a site reaches, derived from what the polymerase does.

A primer P occurring literally at `i` anneals to the MINUS strand there. The
nascent strand is plus-sense and extension runs toward increasing coordinates.
An occurrence of `rc(P)` at `j` is the mirror: P anneals to the PLUS strand and
extension runs toward decreasing coordinates. One direction per occurrence.

`PositionCache` already separates the two -- `"forward"` is `db[primer]` and
`"reverse"` is `db[rc(primer)]`, both plus-strand offsets
(`position_cache.py:374-386`) -- so nothing needs re-indexing.

This holds both models in one function so a caller chooses rather than inherits.
Nothing selects `directional` yet: the three stages of a run currently disagree
about both direction and width, and
`docs/validation/2026-09-23-directional-coverage.md` records which is which.
Picking one is a later increment and needs the reach refitted under it.

**The oracle below is written from the two sentences at the top and calls
nothing in `coverage.py`.** `tests/test_coverage_independent_oracle.py` cannot
serve here: its own oracle hard-codes `low, high = pos - reach, pos + reach` at
`:43`, so it is independent of the IMPLEMENTATION and not of the CONVENTION,
and adapting it would make the new model agree with itself.
"""

import pytest

from neoswga.core.coverage import COVERAGE_GEOMETRIES, site_spans


def oracle(forward, reverse, reach, geometry):
    """The mechanism, written out. Deliberately slow and deliberately naive."""
    if geometry == "symmetric":
        return sorted((pos - reach, pos + reach) for pos in set(list(forward) + list(reverse)))
    spans = [(i, i + reach) for i in forward]
    spans += [(j - reach, j) for j in reverse]
    return sorted(spans)


# ---------------------------------------------------------------------------
# The mechanism
# ---------------------------------------------------------------------------


def test_a_forward_site_reaches_only_forward():
    """P occurs literally, so it primes the minus strand and copies toward
    increasing coordinates. The oligo's own footprint is neglected; the
    docstring says why."""
    assert site_spans([1000], [], 100, "directional") == [(1000, 1100)]


def test_a_reverse_site_reaches_only_backward():
    """`rc(P)` occurs, so P primes the plus strand and copies toward decreasing
    coordinates, ending at the site."""
    assert site_spans([], [1000], 100, "directional") == [(900, 1000)]


def test_the_two_directions_do_not_overlap_at_the_same_coordinate():
    """The headline. A site at one coordinate reaches one way, not both, and a
    forward and a reverse site at the SAME coordinate reach opposite ways."""
    spans = site_spans([1000], [1000], 100, "directional")

    assert sorted(spans) == [(900, 1000), (1000, 1100)]


# ---------------------------------------------------------------------------
# The symmetric mode is what every recorded figure was produced under
# ---------------------------------------------------------------------------


def test_the_symmetric_mode_is_unchanged():
    assert site_spans([1000], [], 100, "symmetric") == [(900, 1100)]


def test_a_palindrome_gains_no_span_in_symmetric_mode():
    """`get_positions(..., "both")` returns `np.unique(...)`, so a
    self-reverse-complementary oligo has one site there. The symmetric mode has
    to dedup or it credits that site twice and inflates coverage for exactly
    the oligos a dimer screen is already suspicious of."""
    assert site_spans([500], [500], 100, "symmetric") == [(400, 600)]


def test_the_directional_mode_does_not_dedup():
    """The mirror of the test above, and not an oversight. A forward and a
    reverse occurrence at one coordinate are two different priming events
    reaching opposite ways, so collapsing them would delete a real span."""
    assert len(site_spans([500], [500], 100, "directional")) == 2


# ---------------------------------------------------------------------------
# Against the oracle
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("geometry", COVERAGE_GEOMETRIES)
@pytest.mark.parametrize(
    "forward,reverse",
    [
        ([], []),
        ([0], []),
        ([], [0]),
        ([5, 5000, 120000], [17, 4999]),
        ([100, 100, 100], [100]),
        ([2**31 + 5], [2**31 + 9]),
    ],
    ids=["empty", "one-forward", "one-reverse", "scattered", "repeated", "past-int32"],
)
def test_the_spans_match_an_oracle_written_from_the_mechanism(forward, reverse, geometry):
    assert sorted(site_spans(forward, reverse, 3000, geometry)) == oracle(
        forward, reverse, 3000, geometry
    )


def test_spans_may_fall_outside_the_genome_and_the_caller_clips():
    """`_mark_window` and `merged_window_intervals` already clip, and doing it
    twice would mean two places deciding what a genome edge is. A circular
    reference wraps rather than clipping, which only the caller knows."""
    assert site_spans([10], [], 100, "directional") == [(10, 110)]
    assert site_spans([], [10], 100, "directional") == [(-90, 10)]


# ---------------------------------------------------------------------------
# Refusals
# ---------------------------------------------------------------------------


def test_an_unknown_geometry_is_refused():
    """A typo must not silently pick a model. `fg_coverage` is the
    authoritative figure in every saved summary, and a default here would be a
    silent answer to the question this function exists to make explicit."""
    with pytest.raises(ValueError) as caught:
        site_spans([1], [], 100, "radial")

    assert "radial" in str(caught.value)
    assert "symmetric" in str(caught.value)


def test_both_geometries_are_named_in_the_public_tuple():
    assert COVERAGE_GEOMETRIES == ("symmetric", "directional")
