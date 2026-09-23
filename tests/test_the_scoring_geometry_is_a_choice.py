"""`fg_coverage` can be asked for either geometry, and asks for the old one.

`_union_coverage` took one pooled position list and marked a symmetric window
around each site. The list is pooled twice before it arrives -- once by
`get_positions(..., "both")`, once by `compute_metrics` unioning across primers
into a set -- so the orientation a site had is gone by then, and a directional
span cannot be formed.

It now takes the two orientations and a geometry. `symmetric` is the default
and must answer exactly what it answered before, which is what most of this
file checks: the recorded figures were all produced under it, and
`docs/validation/2026-09-23-reach-refit-directional.md` shows the reach itself
was fitted alongside it.

The one behaviour that is NOT a simple pass-through is the circular
short-circuit. A symmetric window is `2 * reach` wide and a directional span is
`reach`, so the test for "one site covers the whole circle" has to move with
the geometry. Left alone it would report 1.0 for a directional model that
covers half the sequence.
"""

import numpy as np
import pytest

from neoswga.core.base_optimizer import _union_coverage


def test_splitting_the_input_does_not_change_the_symmetric_answer():
    """The invariant the default rests on. `compute_metrics` will pass the two
    orientations where it used to pass their union, and the symmetric mode has
    to be blind to the difference."""
    forward, reverse = [100, 5000], [2500, 5000]
    pooled = sorted(set(forward) | set(reverse))

    assert _union_coverage(pooled, 20000, 300, False) == _union_coverage(
        forward, 20000, 300, False, reverse=reverse, geometry="symmetric"
    )


def test_a_site_in_both_orientations_is_counted_once_when_symmetric():
    """`get_positions(..., "both")` dedups with `np.unique`, so the symmetric
    mode must too, or a palindromic oligo gains coverage it does not have."""
    assert _union_coverage([500], 10000, 100, False) == _union_coverage(
        [500], 10000, 100, False, reverse=[500], geometry="symmetric"
    )


def test_the_directional_mode_covers_less_at_the_same_reach():
    """Half the width per site, so strictly less unless the spans happen to
    tile. Not a claim about which is right -- the reach was refitted for that."""
    forward, reverse = [2000], [8000]

    symmetric = _union_coverage(forward, 20000, 300, False, reverse=reverse)
    directional = _union_coverage(
        forward, 20000, 300, False, reverse=reverse, geometry="directional"
    )

    assert directional < symmetric


def test_the_circular_short_circuit_moves_with_the_geometry():
    """The one place the change is not a pass-through.

    A symmetric window of `2 * reach` covers a 500 bp circle from any single
    site when reach is 300. A directional span of 300 does not. Leaving the
    test at `2 * reach >= length` would report 1.0 for a model covering 60%.
    """
    assert _union_coverage([100], 500, 300, True) == 1.0

    directional = _union_coverage([100], 500, 300, True, geometry="directional")
    assert directional < 1.0
    assert directional == pytest.approx(300 / 500)


def test_a_directional_span_still_wraps_on_a_circular_reference():
    """Wrapping is a property of the reference, not of the geometry."""
    covered = _union_coverage([50], 1000, 300, True, reverse=[], geometry="directional")

    assert covered == pytest.approx(300 / 1000)


def test_no_sites_is_zero_in_both_geometries():
    """Absence is not a full circle. An empty cache reporting 1.0 on a small
    circular target is a defect `compute_per_prefix_coverage` already carries a
    guard against."""
    for geometry in ("symmetric", "directional"):
        assert _union_coverage([], 500, 300, True, reverse=[], geometry=geometry) == 0.0


def test_an_unknown_geometry_is_refused():
    with pytest.raises(ValueError):
        _union_coverage([1], 500, 10, False, geometry="radial")
