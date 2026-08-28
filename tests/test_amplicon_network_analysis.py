"""Tests for the amplicon-network analysis methods.

`AmpliconNetwork` models which pairs of binding sites could produce an amplicon:
a forward primer and a downstream reverse primer within the amplifiable length
range. Everything it reports -- coverage, hubs, dead ends, amplicon statistics,
critical primers -- is read off that graph.

Until the loader was fixed it read no reverse-strand sites, so the graph had
nodes and zero edges for every input, and every one of these methods was being
exercised against an empty graph. Fixing that exposed a second bug they had been
hiding: the component analysis used STRONGLY connected components, and amplicon
edges only run forward -> downstream reverse, so the graph is acyclic and every
such component is a single node.

The fixture below builds a network with known, hand-checkable geometry so the
numbers can be verified rather than merely observed.
"""

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")
nx = pytest.importorskip("networkx")

from neoswga.core.amplicon_network import AmpliconNetwork
from neoswga.core.thermodynamics import reverse_complement

PRIMER = "TTGACCATGA"
GENOME_LENGTH = 20_000

# Two clusters, deliberately far enough apart that no amplicon spans them.
FORWARD_SITES = [1_000, 4_000, 8_000]
REVERSE_SITES = [2_500, 6_000, 9_500]


@pytest.fixture
def network(tmp_path):
    prefix = str(tmp_path / "g")
    with h5py.File(f"{prefix}_{len(PRIMER)}mer_positions.h5", "w") as f:
        f.create_dataset(PRIMER, data=np.array(FORWARD_SITES, dtype=np.int32))
        f.create_dataset(reverse_complement(PRIMER), data=np.array(REVERSE_SITES, dtype=np.int32))

    net = AmpliconNetwork(primers=[PRIMER], genome_length=GENOME_LENGTH)
    net.load_positions_from_hdf5(prefix)
    net.build_network()
    return net


@pytest.fixture
def empty_network():
    return AmpliconNetwork(primers=[PRIMER], genome_length=GENOME_LENGTH)


# ----------------------------------------------------------------------
# The graph itself
# ----------------------------------------------------------------------


def test_the_graph_is_acyclic(network):
    """This is why strongly connected components were the wrong tool.

    An amplicon edge runs forward -> downstream reverse, so it always advances
    along the genome and can never close a cycle.
    """
    assert nx.is_directed_acyclic_graph(network.G)


def test_every_edge_points_downstream(network):
    """A reverse primer upstream of the forward one cannot form an amplicon."""
    for fwd, rev in network.G.edges():
        assert network.G.nodes[rev]["position"] > network.G.nodes[fwd]["position"]


def test_edges_respect_the_amplicon_length_limits(network):
    for _f, _r, data in network.G.edges(data=True):
        assert network.min_amplicon_length <= data["length"] <= network.max_amplicon_length


def test_edges_run_from_forward_to_reverse(network):
    for fwd, rev in network.G.edges():
        assert network.G.nodes[fwd]["strand"] == "forward"
        assert network.G.nodes[rev]["strand"] == "reverse"


# ----------------------------------------------------------------------
# Coverage
# ----------------------------------------------------------------------


def test_coverage_is_non_zero_for_a_network_with_amplicons(network):
    """It returned 0.0 for every network.

    Coverage summed over components of size > 1, but every strongly connected
    component of an acyclic graph has size 1, so the branch never ran.
    """
    assert network.calculate_coverage() > 0


def test_coverage_is_a_fraction(network):
    assert 0.0 < network.calculate_coverage() <= 1.0


def test_coverage_spans_the_amplicon_clusters(network):
    """Covered bases are the span of each connected cluster, end inclusive.

    With the default 200-5000 bp amplicon window the four edges are
    1000->2500, 1000->6000, 4000->6000 and 8000->9500, which gives two
    clusters: 1000..6000 (5001 bases) and 8000..9500 (1501). The 4000->9500
    pair is 5500 bp and exceeds the window, which is what keeps them separate.
    """
    expected = ((6_000 - 1_000 + 1) + (9_500 - 8_000 + 1)) / GENOME_LENGTH
    assert network.calculate_coverage() == pytest.approx(expected)


def test_empty_network_has_zero_coverage(empty_network):
    """Zero here is a real answer, not the bug -- there are no nodes at all."""
    assert empty_network.calculate_coverage() == 0.0


# ----------------------------------------------------------------------
# Components
# ----------------------------------------------------------------------


def test_components_are_amplicon_clusters_not_single_nodes(network):
    """`largest_component_size` was always 1, which told the reader nothing."""
    metrics = network.analyze()

    assert metrics.largest_component_size > 1
    assert metrics.connected_components < network.G.number_of_nodes()


def test_component_count_matches_networkx(network):
    metrics = network.analyze()
    expected = len(list(nx.weakly_connected_components(network.G)))
    assert metrics.connected_components == expected


def test_deprecated_alias_still_reads(network):
    """Kept so existing readers do not break on the rename."""
    metrics = network.analyze()
    assert metrics.strongly_connected_components == metrics.connected_components


# ----------------------------------------------------------------------
# Hubs and dead ends
# ----------------------------------------------------------------------


def test_hubs_respect_the_degree_threshold(network):
    for node in network.find_hubs(min_degree=2):
        degree = network.G.in_degree(node) + network.G.out_degree(node)
        assert degree >= 2


def test_a_higher_threshold_finds_no_more_hubs(network):
    assert len(network.find_hubs(min_degree=4)) <= len(network.find_hubs(min_degree=1))


def test_dead_ends_have_no_outgoing_edges(network):
    for node in network.find_dead_ends():
        assert network.G.out_degree(node) == 0


def test_every_reverse_node_is_a_dead_end(network):
    """Edges only leave forward nodes, so reverse nodes terminate every chain.

    That is what makes them dead ends, and it is a structural property worth
    pinning: if a reverse node ever gained an outgoing edge, the graph would no
    longer be acyclic and the component analysis would change meaning.
    """
    dead_ends = set(network.find_dead_ends())
    reverse_nodes = {n for n, d in network.G.nodes(data=True) if d["strand"] == "reverse"}
    assert reverse_nodes <= dead_ends


# ----------------------------------------------------------------------
# Amplicon statistics
# ----------------------------------------------------------------------


def test_amplicon_statistics_describe_the_edges(network):
    stats = network.calculate_amplicon_statistics()

    assert stats["count"] == network.G.number_of_edges()
    assert stats["min"] <= stats["mean"] <= stats["max"]
    assert stats["min"] >= network.min_amplicon_length


def test_amplicon_statistics_on_an_empty_graph(empty_network):
    stats = empty_network.calculate_amplicon_statistics()
    assert stats["count"] == 0
    assert stats["mean"] == 0


def test_amplicon_lengths_match_the_planted_geometry(network):
    """Hand-checkable: the shortest amplicon here is 1000 -> 2500."""
    stats = network.calculate_amplicon_statistics()
    assert stats["min"] == 1_500


# ----------------------------------------------------------------------
# Critical primers
# ----------------------------------------------------------------------


def test_critical_primers_are_scored_and_sorted(network):
    ranked = network.find_critical_primers()

    assert ranked
    scores = [score for _primer, score in ranked]
    assert scores == sorted(scores, reverse=True)


def test_critical_primers_only_names_primers_in_the_network(network):
    for primer, _score in network.find_critical_primers():
        assert primer == PRIMER


# ----------------------------------------------------------------------
# The full analysis
# ----------------------------------------------------------------------


def test_analyze_reports_finite_metrics(network):
    import math

    metrics = network.analyze()
    for name, value in vars(metrics).items():
        if isinstance(value, float):
            assert math.isfinite(value), f"{name} is not finite"


def test_analyze_on_an_empty_network_does_not_crash(empty_network):
    """A design that yields no binding sites still has to produce a report."""
    metrics = empty_network.analyze()

    assert metrics.coverage_fraction == 0.0
    assert metrics.num_amplicons == 0
    assert metrics.largest_component_size == 0


def test_subgraph_extraction_is_bounded(network):
    """Visualisation caps the node count; an unbounded graph would hang a plot."""
    sub = network.visualize_subgraph(max_nodes=3)
    assert sub.number_of_nodes() <= 3
