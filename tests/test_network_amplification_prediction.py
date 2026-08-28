"""Tests for the amplification prediction that Stage-2 optimization scores on.

`AmplificationNetwork.predict_amplification_fold` turns network structure into a
predicted fold-amplification. It is reported to users as
`stage2_predicted_amplification`, and `background_aware_optimizer` ranks
candidate sets on it directly -- so its shape decides which primer sets ship,
not just what gets printed.

The property that matters is monotonicity. A better-connected network cannot
predict less amplification than a worse-connected one; if it does, the score
actively prefers the wrong set. That is what these tests pin, alongside the
network construction the prediction reads from.
"""

import math

import numpy as np
import pytest

nx = pytest.importorskip("networkx")

from neoswga.core.network_optimizer import (
    _LINEAR_COMPONENT_LIMIT,
    _LINEAR_FOLD_PER_SITE,
    AmplificationNetwork,
)


def make_network(n_pairs, spacing=1_000, reach=70_000):
    """A network with `n_pairs` forward and `n_pairs` reverse sites.

    Opposite strands within reach of each other, so every pair can form an
    amplification edge.
    """
    net = AmplificationNetwork(max_extension=reach)
    net.add_primer_sites("P", np.array([i * spacing for i in range(n_pairs)]), "+")
    net.add_primer_sites("P", np.array([i * spacing + spacing // 2 for i in range(n_pairs)]), "-")
    net.build_edges()
    return net


# ----------------------------------------------------------------------
# Monotonicity -- the property the score depends on
# ----------------------------------------------------------------------


def test_prediction_is_monotone_in_component_size():
    """A bigger connected component must never predict less amplification.

    The two branches used to be independent: `largest * 5` below a component of
    10, `2 ** (largest / 10)` at or above it. They did not meet. Nine sites
    predicted 45-fold and ten predicted 2-fold, a 22x drop, and the exponential
    did not regain 45 until a component of 55 -- so every network with a largest
    component between 10 and 54 was predicted to amplify LESS than one with 9.

    `background_aware_optimizer` ranks sets on this value, so a better-connected
    network scored worse over that whole range.
    """
    predictions = [make_network(k).predict_amplification_fold() for k in range(1, 40)]

    for smaller, larger in zip(predictions, predictions[1:]):
        assert larger >= smaller, "a larger network predicted less amplification"


def test_the_two_branches_meet_at_the_crossover():
    """Continuity is what makes the function usable as a score."""
    below = make_network((_LINEAR_COMPONENT_LIMIT - 1) // 2).predict_amplification_fold()
    at = make_network(_LINEAR_COMPONENT_LIMIT // 2 + 1).predict_amplification_fold()

    assert at >= below
    # The jump across the crossover should be small, not an order of magnitude.
    assert at < below * 2


def test_linear_branch_below_the_crossover():
    """Small components cannot branch, so growth is linear in their size."""
    net = make_network(2)  # 4 sites, one component of 4
    assert net.largest_component_size() < _LINEAR_COMPONENT_LIMIT
    assert net.predict_amplification_fold() == pytest.approx(
        net.largest_component_size() * _LINEAR_FOLD_PER_SITE
    )


def test_growth_is_faster_than_linear_above_the_crossover():
    """The exponential branch has to actually be exponential, or the
    distinction between connected and isolated primers -- the premise of this
    whole module -- carries no weight."""
    small = make_network(10)
    large = make_network(30)

    ratio_sites = large.largest_component_size() / small.largest_component_size()
    ratio_pred = large.predict_amplification_fold() / small.predict_amplification_fold()

    assert ratio_pred > ratio_sites


def test_prediction_saturates():
    """A cap keeps an enormous network from reporting an absurd fold."""
    huge = make_network(400)
    assert math.isfinite(huge.predict_amplification_fold())


def test_empty_network_predicts_nothing():
    net = AmplificationNetwork(max_extension=70_000)
    assert net.predict_amplification_fold() == 0


# ----------------------------------------------------------------------
# Uncertainty bounds
# ----------------------------------------------------------------------


def test_confidence_interval_brackets_the_estimate():
    net = make_network(12)
    point, lower, upper = net.predict_amplification_with_uncertainty()

    assert lower <= point <= upper


def test_confidence_bounds_are_finite():
    """These are reported to users; inf or nan would render as a measurement."""
    for k in (1, 5, 12, 40):
        for value in make_network(k).predict_amplification_with_uncertainty():
            assert math.isfinite(value)


def test_more_sites_narrows_the_relative_interval():
    """Uncertainty is documented to fall with more binding sites (sqrt(N))."""

    def relative_width(k):
        point, lower, upper = make_network(k).predict_amplification_with_uncertainty()
        return (upper - lower) / point if point else 0.0

    assert relative_width(40) < relative_width(3)


def test_a_network_with_no_amplification_reports_unit_bounds():
    net = AmplificationNetwork(max_extension=70_000)
    assert net.predict_amplification_with_uncertainty() == (1.0, 1.0, 1.0)


# ----------------------------------------------------------------------
# The network the prediction reads from
# ----------------------------------------------------------------------


def test_edges_only_join_opposite_strands():
    """An amplicon needs a forward and a reverse primer.

    Two same-strand sites cannot prime toward each other, so an edge between
    them would inflate connectivity with pairs that cannot amplify.
    """
    net = make_network(6)
    for a, b in net.graph.edges():
        assert a.strand != b.strand


def test_sites_beyond_the_extension_reach_are_not_joined():
    """Reach is what makes this a model of phi29 rather than of adjacency."""
    net = AmplificationNetwork(max_extension=1_000)
    net.add_primer_sites("P", np.array([0]), "+")
    net.add_primer_sites("P", np.array([50_000]), "-")
    net.build_edges()

    assert net.graph.number_of_edges() == 0


def test_a_single_strand_forms_no_edges():
    """The failure mode found in AmpliconNetwork: only forward sites loaded.

    It produces nodes and no edges, which reads as a real result rather than as
    missing data, so it is worth having a test that names the shape.
    """
    net = AmplificationNetwork(max_extension=70_000)
    net.add_primer_sites("P", np.array([0, 1_000, 2_000]), "+")
    net.build_edges()

    assert net.graph.number_of_nodes() == 3
    assert net.graph.number_of_edges() == 0
    assert net.predict_amplification_fold() < _LINEAR_FOLD_PER_SITE * 2


def test_incremental_edge_building_matches_a_full_rebuild():
    """`add_edges_for_sites` is an optimisation of `build_edges`.

    If the two disagree, an optimizer that adds primers one at a time gets a
    different network from one that builds in a single pass.
    """
    full = make_network(8)

    incremental = AmplificationNetwork(max_extension=70_000)
    fwd = incremental.add_primer_sites("P", np.array([i * 1_000 for i in range(8)]), "+")
    incremental.add_edges_for_sites(fwd)
    rev = incremental.add_primer_sites("P", np.array([i * 1_000 + 500 for i in range(8)]), "-")
    incremental.add_edges_for_sites(rev)

    assert incremental.graph.number_of_edges() == full.graph.number_of_edges()
    assert incremental.largest_component_size() == full.largest_component_size()


def test_statistics_are_all_finite():
    """`get_statistics` feeds reports and the optimizer summary."""
    stats = make_network(12).get_statistics()

    for name, value in stats.items():
        if isinstance(value, float):
            assert math.isfinite(value), f"{name} is not finite"


def test_statistics_on_an_empty_network():
    stats = AmplificationNetwork(max_extension=70_000).get_statistics()
    assert stats["num_sites"] == 0
    assert stats["predicted_amplification"] == 0
