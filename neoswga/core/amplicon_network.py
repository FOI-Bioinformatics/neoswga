"""
Amplicon network analysis for SWGA using graph theory.

Models primer binding sites and potential amplicons as a directed graph:
- Nodes: Primer binding positions
- Edges: Potential amplicons (forward -> reverse primer within distance)

Enables:
- Coverage prediction via graph connectivity
- Amplification hotspot identification
- Dead-end detection
- Optimal primer set selection via network properties
"""

import numpy as np
import networkx as nx
from typing import List, Dict, Tuple, Set, Optional
from dataclasses import dataclass
from collections import defaultdict
import h5py

import neoswga.core.thermodynamics as thermo
import neoswga.core.reaction_conditions as rc


@dataclass
class NetworkMetrics:
    """Metrics from amplicon network analysis."""
    coverage_fraction: float
    num_amplicons: int
    num_hubs: int
    mean_amplicon_length: float
    median_amplicon_length: float
    max_amplicon_length: float
    strongly_connected_components: int
    largest_component_size: int
    average_clustering: float
    network_density: float
    mean_degree: float


class AmpliconNetwork:
    """
    Graph-based model of SWGA amplification.

    Constructs and analyzes network of potential amplicons
    to predict coverage and identify optimal primer sets.
    """

    def __init__(self,
                 primers: List[str],
                 genome_length: int,
                 max_amplicon_length: int = 5000,
                 min_amplicon_length: int = 200):
        """
        Initialize amplicon network.

        Args:
            primers: Primer sequences
            genome_length: Length of target genome
            max_amplicon_length: Maximum amplicon length to consider
            min_amplicon_length: Minimum amplicon length
        """
        self.primers = primers
        self.genome_length = genome_length
        self.max_amplicon_length = max_amplicon_length
        self.min_amplicon_length = min_amplicon_length

        self.G = nx.DiGraph()
        self.positions_forward = {}  # primer -> positions
        self.positions_reverse = {}  # primer -> positions

    def load_positions_from_hdf5(self, hdf5_prefix: str):
        """
        Load primer positions from HDF5 files.

        Args:
            hdf5_prefix: Prefix for HDF5 position files
        """
        for primer in self.primers:
            k = len(primer)
            hdf5_file = f"{hdf5_prefix}_{k}mer_positions.h5"

            try:
                with h5py.File(hdf5_file, 'r') as f:
                    if primer in f:
                        positions = f[primer][:]

                        # Separate forward and reverse
                        # Assuming positions stored as: [pos1_fwd, pos2_fwd, ..., -pos1_rev, -pos2_rev, ...]
                        # Negative positions indicate reverse strand
                        forward = positions[positions >= 0]
                        reverse = -positions[positions < 0]

                        self.positions_forward[primer] = sorted(forward)
                        self.positions_reverse[primer] = sorted(reverse)
            except (FileNotFoundError, KeyError):
                warnings.warn(f"Could not load positions for {primer}")
                self.positions_forward[primer] = []
                self.positions_reverse[primer] = []

    def build_network(self):
        """
        Construct amplicon network.

        Creates nodes for binding sites and edges for valid amplicons.
        """
        print("Building amplicon network...")

        node_id = 0
        primer_nodes = defaultdict(list)  # primer -> node IDs

        # Add nodes for all binding sites
        for primer in self.primers:
            # Forward strand sites
            for pos in self.positions_forward[primer]:
                self.G.add_node(
                    node_id,
                    primer=primer,
                    position=pos,
                    strand='forward'
                )
                primer_nodes[(primer, 'forward')].append(node_id)
                node_id += 1

            # Reverse strand sites
            for pos in self.positions_reverse[primer]:
                self.G.add_node(
                    node_id,
                    primer=primer,
                    position=pos,
                    strand='reverse'
                )
                primer_nodes[(primer, 'reverse')].append(node_id)
                node_id += 1

        print(f"  Nodes: {self.G.number_of_nodes():,}")

        # Add edges for valid amplicons
        # Amplicon: forward primer -> downstream reverse primer
        edge_count = 0

        for primer_fwd in self.primers:
            fwd_nodes = primer_nodes[(primer_fwd, 'forward')]

            for primer_rev in self.primers:
                rev_nodes = primer_nodes[(primer_rev, 'reverse')]

                # Connect forward -> reverse if within distance limits
                for fwd_node in fwd_nodes:
                    fwd_pos = self.G.nodes[fwd_node]['position']

                    for rev_node in rev_nodes:
                        rev_pos = self.G.nodes[rev_node]['position']

                        amplicon_length = rev_pos - fwd_pos

                        if self.min_amplicon_length <= amplicon_length <= self.max_amplicon_length:
                            self.G.add_edge(
                                fwd_node, rev_node,
                                length=amplicon_length,
                                primers=(primer_fwd, primer_rev)
                            )
                            edge_count += 1

        print(f"  Edges: {edge_count:,}")
        print(f"  Potential amplicons: {edge_count:,}")

    def calculate_coverage(self) -> float:
        """
        Calculate genome coverage from network.

        Uses connected components to estimate covered regions.

        Returns:
            Coverage fraction (0-1)
        """
        if self.G.number_of_nodes() == 0:
            return 0.0

        # Find all strongly connected components
        components = list(nx.strongly_connected_components(self.G))

        # Calculate covered bases
        covered = set()

        for component in components:
            if len(component) > 1:
                positions = [self.G.nodes[n]['position'] for n in component]
                start = min(positions)
                end = max(positions)

                # Add all bases in this region
                covered.update(range(start, end + 1))

        coverage = len(covered) / self.genome_length

        return coverage

    def find_hubs(self, min_degree: int = 5) -> List[int]:
        """
        Find hub nodes (high-degree vertices).

        Hubs drive amplification due to multiple incoming/outgoing edges.

        Args:
            min_degree: Minimum degree to be considered hub

        Returns:
            List of hub node IDs
        """
        hubs = []

        for node in self.G.nodes():
            degree = self.G.in_degree(node) + self.G.out_degree(node)
            if degree >= min_degree:
                hubs.append(node)

        return hubs

    def find_dead_ends(self) -> List[int]:
        """
        Find dead-end nodes (out-degree 0).

        Dead ends terminate amplification chains.

        Returns:
            List of dead-end node IDs
        """
        dead_ends = []

        for node in self.G.nodes():
            if self.G.out_degree(node) == 0:
                dead_ends.append(node)

        return dead_ends

    def calculate_amplicon_statistics(self) -> Dict:
        """
        Calculate statistics on amplicon lengths.

        Returns:
            Dictionary with length statistics
        """
        if self.G.number_of_edges() == 0:
            return {
                'count': 0,
                'mean': 0,
                'median': 0,
                'std': 0,
                'min': 0,
                'max': 0
            }

        lengths = [data['length'] for _, _, data in self.G.edges(data=True)]

        return {
            'count': len(lengths),
            'mean': np.mean(lengths),
            'median': np.median(lengths),
            'std': np.std(lengths),
            'min': np.min(lengths),
            'max': np.max(lengths)
        }

    def analyze(self) -> NetworkMetrics:
        """
        Comprehensive network analysis.

        Returns:
            NetworkMetrics object
        """
        print("Analyzing amplicon network...")

        coverage = self.calculate_coverage()
        print(f"  Coverage: {coverage:.1%}")

        hubs = self.find_hubs()
        print(f"  Hubs: {len(hubs)}")

        amplicon_stats = self.calculate_amplicon_statistics()
        print(f"  Mean amplicon length: {amplicon_stats['mean']:.0f} bp")

        # Connected components
        components = list(nx.strongly_connected_components(self.G))
        largest_component = max(components, key=len) if components else set()
        print(f"  Connected components: {len(components)}")
        print(f"  Largest component: {len(largest_component)} nodes")

        # Graph metrics
        if self.G.number_of_nodes() > 0:
            density = nx.density(self.G)
            mean_degree = sum(dict(self.G.degree()).values()) / self.G.number_of_nodes()

            # Clustering (convert to undirected for calculation)
            G_undirected = self.G.to_undirected()
            clustering = nx.average_clustering(G_undirected)
        else:
            density = 0.0
            mean_degree = 0.0
            clustering = 0.0

        return NetworkMetrics(
            coverage_fraction=coverage,
            num_amplicons=amplicon_stats['count'],
            num_hubs=len(hubs),
            mean_amplicon_length=amplicon_stats['mean'],
            median_amplicon_length=amplicon_stats['median'],
            max_amplicon_length=amplicon_stats['max'],
            strongly_connected_components=len(components),
            largest_component_size=len(largest_component),
            average_clustering=clustering,
            network_density=density,
            mean_degree=mean_degree
        )

    def visualize_subgraph(self, max_nodes: int = 100) -> nx.DiGraph:
        """
        Extract subgraph for visualization.

        Args:
            max_nodes: Maximum nodes to include

        Returns:
            Subgraph for plotting
        """
        if self.G.number_of_nodes() <= max_nodes:
            return self.G

        # Sample largest connected component
        components = list(nx.strongly_connected_components(self.G))
        largest = max(components, key=len)

        if len(largest) > max_nodes:
            sampled = list(largest)[:max_nodes]
        else:
            sampled = list(largest)

        return self.G.subgraph(sampled)

    def find_critical_primers(self) -> List[Tuple[str, float]]:
        """
        Identify primers critical for network connectivity.

        Uses node betweenness centrality to find primers
        that bridge different regions.

        Returns:
            List of (primer, centrality) sorted by importance
        """
        print("Finding critical primers...")

        # Calculate betweenness centrality
        betweenness = nx.betweenness_centrality(self.G)

        # Aggregate by primer
        primer_centrality = defaultdict(list)

        for node, centrality in betweenness.items():
            primer = self.G.nodes[node]['primer']
            primer_centrality[primer].append(centrality)

        # Average centrality per primer
        primer_importance = []
        for primer, centralities in primer_centrality.items():
            avg_centrality = np.mean(centralities)
            primer_importance.append((primer, avg_centrality))

        # Sort by importance
        primer_importance.sort(key=lambda x: x[1], reverse=True)

        return primer_importance


def network_based_primer_selection(
    primer_candidates: List[str],
    genome_length: int,
    hdf5_prefix: str,
    target_set_size: int = 6,
    max_amplicon_length: int = 5000
) -> List[str]:
    """
    Select primer set that maximizes network connectivity.

    Greedy algorithm:
    1. Start with empty set
    2. Add primer that maximally improves coverage
    3. Repeat until target size

    Args:
        primer_candidates: Available primers
        genome_length: Target genome length
        hdf5_prefix: HDF5 file prefix for positions
        target_set_size: Desired number of primers
        max_amplicon_length: Max amplicon length

    Returns:
        Selected primer set
    """
    selected = []
    remaining = set(primer_candidates)

    print(f"Network-based primer selection (target: {target_set_size} primers)")

    for iteration in range(target_set_size):
        best_primer = None
        best_coverage = 0.0

        print(f"\nIteration {iteration + 1}/{target_set_size}")

        for candidate in remaining:
            # Test adding this primer
            test_set = selected + [candidate]

            # Build network
            network = AmpliconNetwork(
                test_set, genome_length, max_amplicon_length
            )
            network.load_positions_from_hdf5(hdf5_prefix)
            network.build_network()

            # Calculate coverage
            coverage = network.calculate_coverage()

            if coverage > best_coverage:
                best_coverage = coverage
                best_primer = candidate

        if best_primer:
            selected.append(best_primer)
            remaining.remove(best_primer)
            print(f"  Selected: {best_primer}")
            print(f"  Coverage: {best_coverage:.1%}")
        else:
            print("  No improvement possible")
            break

    return selected


if __name__ == "__main__":
    print("Amplicon Network Analysis")
    print("=" * 60)
    print("\nGraph-based modeling of SWGA amplification:")
    print("  - Nodes: Primer binding sites")
    print("  - Edges: Potential amplicons")
    print("  - Analysis: Coverage, hubs, connectivity")
    print("\nKey insights:")
    print("  - High-degree nodes drive exponential amplification")
    print("  - Connected components predict covered regions")
    print("  - Betweenness centrality identifies critical primers")
    print("\nExample usage:")
    print("""
    from neoswga.core import amplicon_network

    # Build network
    network = amplicon_network.AmpliconNetwork(
        primers=['ATCGAT', 'GCTAGC', 'TTAACC'],
        genome_length=4400000,
        max_amplicon_length=5000
    )

    # Load positions from HDF5
    network.load_positions_from_hdf5('target_kmers')

    # Construct network
    network.build_network()

    # Analyze
    metrics = network.analyze()
    print(f"Coverage: {metrics.coverage_fraction:.1%}")
    print(f"Amplicons: {metrics.num_amplicons}")
    print(f"Hubs: {metrics.num_hubs}")

    # Find critical primers
    critical = network.find_critical_primers()
    print("Most critical primers:")
    for primer, centrality in critical[:5]:
        print(f"  {primer}: {centrality:.4f}")
    """)
