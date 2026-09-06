"""`coverage_uniformity` must reward even spacing, not punish it.

It returned the coefficient of variation of the component CENTROID POSITIONS,
with lower treated as better. That is a measure of how tightly the components
sit together, not of how evenly they are spread:

    eight binding islands spread across a 4 Mb genome   0.6546
    the same eight packed inside one 5 kb island        0.0000

so the term scored a set that amplifies one small region far better than a set
that covers the whole target, and `uniformity_weight` -- which `--application`
sets, highest for `metagenomics` and `discovery` -- pushed selection toward
concentrating every binding site in one place.

Until the selection weights were forwarded into the inner optimizer, this term
only reached the `network` method. It now reaches `hybrid`, the default, so the
sign matters much more than it did.

The quantity wanted is the evenness of the GAPS between components, including
the run-up from the start of the sequence and the run-out to its end. Evenly
spaced components give equal gaps and a CV of zero; everything in one cluster
gives one enormous gap and a CV near the maximum. Lower stays better, so no
call site changes.
"""

import numpy as np
import pytest

pytest.importorskip("networkx")

from neoswga.core.network_optimizer import AmplificationNetwork

GENOME = 4_000_000
MAX_EXTENSION = 10_000


def _network(island_starts):
    """One two-site component per island, each island its own component."""
    network = AmplificationNetwork(max_extension=MAX_EXTENSION)
    for start in island_starts:
        network.add_primer_sites("P", np.array([start], dtype=np.int64), "+")
        network.add_primer_sites("P", np.array([start + 400], dtype=np.int64), "-")
    network.build_edges()
    return network


#: Eight components either way, so the comparison is about WHERE they sit
#: and not how many there are. The clustered islands are spaced 25 kb apart
#: -- above `max_extension`, so they stay separate components -- but all of
#: them sit inside a 200 kb stretch of a 4 Mb genome.
SPREAD = [i * (GENOME // 8) + 50_000 for i in range(8)]
CLUSTERED = [2_000_000 + i * 25_000 for i in range(8)]


def test_the_two_layouts_really_are_different():
    """Guard the guard: both must produce eight separate components, or the
    comparison is measuring component counts rather than their spacing."""
    assert len(list(__import__("networkx").connected_components(_network(SPREAD).graph))) == 8
    assert len(list(__import__("networkx").connected_components(_network(CLUSTERED).graph))) == 8


def test_evenly_spread_components_score_better_than_a_single_cluster():
    spread = _network(SPREAD).coverage_uniformity(GENOME)
    clustered = _network(CLUSTERED).coverage_uniformity(GENOME)

    assert spread < clustered, (
        f"components spread across the genome scored {spread:.4f} and the same "
        f"number packed into one 5 kb island scored {clustered:.4f}. Lower is "
        "treated as better, so this rewards clustering."
    )


def test_perfectly_even_spacing_scores_zero():
    """The optimum of the measure, placed exactly.

    N components make N+1 gaps once the run-up from the start and the run-out
    to the end are counted, so equal gaps means placing them at i/(N+1) of the
    sequence -- not at (i+0.5)/N, which leaves the two end gaps half the size of
    the interior ones and scores a little above zero.
    """
    n = 8
    even = [int((i + 1) * GENOME / (n + 1)) for i in range(n)]
    assert _network(even).coverage_uniformity(GENOME) == pytest.approx(0.0, abs=1e-3)


def test_an_empty_network_is_the_worst_score():
    empty = AmplificationNetwork(max_extension=MAX_EXTENSION)
    empty.build_edges()
    assert empty.coverage_uniformity(GENOME) == float("inf")
