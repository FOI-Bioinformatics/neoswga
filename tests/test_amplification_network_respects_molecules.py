"""An amplification edge must not join sites on different molecules.

`_build_network` added every prefix's binding positions into ONE
`AmplificationNetwork`, and `PositionCache` returns genome-local coordinates
with nothing offsetting them. A forward site at 500 on contig C1 and a reverse
site at 900 on contig C2 therefore looked 400 bp apart and were joined by an
edge, describing an amplicon that spans two separate pieces of DNA.

This changes selection rather than only reporting. `largest_component_size` is
the entire foreground term of the greedy objective
(`_evaluate_primer_addition`), so phantom cross-contig connectivity moves which
primers get picked on any draft assembly, any plasmid-plus-chromosome design,
and any multi-target run.
"""

import numpy as np
import pytest

pytest.importorskip("networkx")

from neoswga.core.network_optimizer import NetworkOptimizer

CONTIG_LENGTH = 50_000
MAX_EXTENSION = 10_000


class _Cache:
    """Positions keyed by (prefix, primer, strand)."""

    def __init__(self, positions):
        self._positions = positions

    def get_positions(self, prefix, primer, strand="both"):
        return np.array(self._positions.get((prefix, primer, strand), []), dtype=np.int64)


def _optimizer(positions, prefixes):
    return NetworkOptimizer(
        position_cache=_Cache(positions),
        fg_prefixes=prefixes,
        bg_prefixes=[],
        fg_seq_lengths=[CONTIG_LENGTH] * len(prefixes),
        bg_seq_lengths=[],
        max_extension=MAX_EXTENSION,
    )


def test_two_sites_on_the_same_contig_are_joined():
    """Guard the guard: the pair below must be joinable in principle, or the
    cross-contig test proves only that nothing connects at all."""
    positions = {
        ("C1", "P", "forward"): [500],
        ("C1", "P", "reverse"): [900],
    }
    network = _optimizer(positions, ["C1"])._build_network(["P"], ["C1"])

    assert network.largest_component_size() == 2, (
        "a forward site at 500 and a reverse site at 900 on ONE contig are "
        f"400 bp apart and must form one amplicon; got "
        f"{network.largest_component_size()}"
    )


def test_sites_on_different_contigs_are_not_joined():
    """The same two coordinates, on two different molecules."""
    positions = {
        ("C1", "P", "forward"): [500],
        ("C2", "P", "reverse"): [900],
    }
    network = _optimizer(positions, ["C1", "C2"])._build_network(["P"], ["C1", "C2"])

    assert network.largest_component_size() == 1, (
        "a forward site on C1 and a reverse site on C2 were joined into one "
        "amplification component; they are on separate pieces of DNA and no "
        "extension can bridge them"
    )


def test_each_contig_keeps_its_own_connectivity():
    """Real structure within a contig survives the separation."""
    positions = {
        ("C1", "P", "forward"): [500, 5_000],
        ("C1", "P", "reverse"): [900, 5_400],
        ("C2", "P", "forward"): [500],
        ("C2", "P", "reverse"): [900],
    }
    network = _optimizer(positions, ["C1", "C2"])._build_network(["P"], ["C1", "C2"])

    # C1 holds four mutually reachable sites; C2 holds two. Neither may merge.
    assert network.largest_component_size() == 4, (
        f"expected C1's four sites to form the largest component, got "
        f"{network.largest_component_size()}"
    )
    assert (
        network.num_components() == 2
    ), f"expected one component per contig, got {network.num_components()}"


def test_the_simulated_addition_uses_the_same_separation():
    """`_simulate_add_primer` is the scoring path, and it queries the cache
    prefix by prefix exactly as `_build_network` does. If only one of them
    separates the molecules, the score disagrees with the network it claims to
    describe."""
    positions = {
        ("C1", "SEED", "forward"): [500],
        ("C1", "SEED", "reverse"): [900],
        ("C1", "CAND", "forward"): [40_000],
        ("C2", "CAND", "reverse"): [40_400],
    }
    optimizer = _optimizer(positions, ["C1", "C2"])
    network = optimizer._build_network(["SEED"], ["C1", "C2"])

    simulated = optimizer._simulate_add_primer(network, "CAND", ["C1", "C2"])

    assert simulated.largest_component_size() == 2, (
        "the candidate's two sites are on different contigs, so adding it must "
        "not create a component larger than the seed's; got "
        f"{simulated.largest_component_size()}"
    )
