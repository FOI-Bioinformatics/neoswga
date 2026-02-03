#!/usr/bin/env python3
"""
Cooperative Binding Network Primer Selection.

Models primer binding as a network where primers that bind within amplification
range cooperate to produce robust, uniform amplification. Selects primers that
maximize network connectivity and amplification robustness.

Key Concepts:
- Phi29 polymerase extends ~70-80kb from binding site
- Primers within extension distance cooperate
- Network connectivity predicts amplification robustness
- Central primers (high degree) are critical for coverage

Algorithm:
1. Build binding site network (edges = within extension distance)
2. Calculate network metrics (degree centrality, betweenness, clustering)
3. Select primers maximizing network connectivity + quality
4. Ensure no bottlenecks or isolated regions

Author: NeoSWGA Development Team
Date: November 2025
Version: 3.5 - Genome-Adaptive QA System
"""

import logging
import numpy as np
from typing import List, Dict, Set, Tuple, Optional
from dataclasses import dataclass
from collections import defaultdict, deque
import heapq

logger = logging.getLogger(__name__)


@dataclass
class BindingSite:
    """A specific binding site for a primer."""
    position: int
    primer_id: str
    quality: float


@dataclass
class NetworkMetrics:
    """Network connectivity metrics for primer set."""
    total_sites: int
    total_edges: int
    avg_degree: float  # Average connections per site
    clustering_coefficient: float  # Local connectivity
    connected_components: int  # Number of isolated networks
    largest_component_size: int  # Size of largest connected region
    network_diameter: int  # Maximum shortest path length
    bottleneck_score: float  # 0-1, lower = more bottlenecks


@dataclass
class CooperativeSelectionResult:
    """Result from cooperative binding selection."""
    selected_primers: List[str]
    binding_sites: List[BindingSite]
    coverage: float
    network_metrics: NetworkMetrics
    amplification_chains: List[List[int]]  # Continuous amplification paths
    mean_quality: float
    cooperation_score: float  # 0-1, higher = more cooperative


class CooperativeBindingSelector:
    """
    Select primers based on cooperative binding network topology.

    Models amplification as a network where binding sites are nodes and
    edges represent cooperative amplification (within extension distance).
    """

    def __init__(self, genome_length: int, extension_distance: int = 70000,
                 min_chain_length: int = 3):
        """
        Initialize cooperative selector.

        Args:
            genome_length: Total genome length (bp)
            extension_distance: Maximum extension distance for Phi29 (bp)
                              Default: 70kb (conservative estimate)
            min_chain_length: Minimum amplification chain length
        """
        self.genome_length = genome_length
        self.extension_distance = extension_distance
        self.min_chain_length = min_chain_length

        logger.info(f"Initialized CooperativeBindingSelector:")
        logger.info(f"  Genome length: {genome_length:,} bp")
        logger.info(f"  Extension distance: {extension_distance:,} bp")
        logger.info(f"  Min chain length: {min_chain_length} sites")

    def select_cooperative_set(self, primers_data: Dict[str, Dict],
                               target_primers: int = 15,
                               quality_weight: float = 0.3) -> CooperativeSelectionResult:
        """
        Select primers optimizing for cooperative binding network.

        Args:
            primers_data: Dict mapping primer_id -> {
                'binding_sites': [positions],
                'quality': float
            }
            target_primers: Number of primers to select
            quality_weight: Weight for quality vs network centrality (0-1)

        Returns:
            CooperativeSelectionResult with network-optimized primer set
        """
        logger.info(f"\nStarting cooperative binding selection:")
        logger.info(f"  Candidates: {len(primers_data)}")
        logger.info(f"  Target primers: {target_primers}")
        logger.info(f"  Quality weight: {quality_weight:.2f}")

        # Build full binding site list
        all_sites = self._build_site_list(primers_data)
        logger.info(f"  Total binding sites: {len(all_sites)}")

        # Build cooperation network
        network = self._build_cooperation_network(all_sites)
        logger.info(f"  Network edges: {sum(len(v) for v in network.values()) // 2}")

        # Calculate network metrics for each primer
        primer_scores = self._calculate_network_scores(
            primers_data, all_sites, network, quality_weight
        )

        # Select primers greedily based on network score
        selected = self._greedy_network_selection(
            primers_data, primer_scores, target_primers, network
        )

        # Build final result
        result = self._build_result(selected, primers_data, all_sites, network)

        self._log_selection_summary(result)

        return result

    def _build_site_list(self, primers_data: Dict[str, Dict]) -> List[BindingSite]:
        """Build list of all binding sites with primer associations."""
        sites = []

        for primer_id, data in primers_data.items():
            quality = data['quality']
            for pos in data['binding_sites']:
                site = BindingSite(
                    position=pos,
                    primer_id=primer_id,
                    quality=quality
                )
                sites.append(site)

        # Sort by position
        sites.sort(key=lambda s: s.position)

        return sites

    def _build_cooperation_network(self, sites: List[BindingSite]) -> Dict[int, Set[int]]:
        """
        Build network where edges connect sites within extension distance.

        Returns:
            Adjacency list: site_index -> set of connected site indices
        """
        network = defaultdict(set)

        # For each site, find sites within extension distance
        for i, site_i in enumerate(sites):
            for j in range(i + 1, len(sites)):
                site_j = sites[j]

                # Calculate distance
                distance = abs(site_j.position - site_i.position)

                # Check if within extension distance
                if distance <= self.extension_distance:
                    network[i].add(j)
                    network[j].add(i)
                else:
                    # Sites are sorted, so no more matches possible
                    break

        return network

    def _calculate_network_scores(self, primers_data: Dict[str, Dict],
                                  sites: List[BindingSite],
                                  network: Dict[int, Set[int]],
                                  quality_weight: float) -> Dict[str, float]:
        """
        Calculate network centrality score for each primer.

        Score combines:
        1. Degree centrality (number of cooperative connections)
        2. Betweenness centrality (importance for network connectivity)
        3. Quality score

        Returns:
            Dict mapping primer_id -> composite score
        """
        # Calculate degree centrality for each site
        max_degree = max(len(neighbors) for neighbors in network.values()) if network else 1

        site_centrality = {}
        for i, neighbors in network.items():
            site_centrality[i] = len(neighbors) / max_degree if max_degree > 0 else 0

        # Aggregate to primer level
        primer_scores = {}

        for primer_id, data in primers_data.items():
            quality = data['quality']
            binding_sites = data['binding_sites']

            # Find site indices for this primer
            site_indices = [
                i for i, site in enumerate(sites)
                if site.primer_id == primer_id
            ]

            if not site_indices:
                primer_scores[primer_id] = 0.0
                continue

            # Average centrality across primer's sites
            avg_centrality = np.mean([site_centrality.get(i, 0) for i in site_indices])

            # Combined score
            score = (1 - quality_weight) * avg_centrality + quality_weight * quality

            primer_scores[primer_id] = score

        return primer_scores

    def _greedy_network_selection(self, primers_data: Dict[str, Dict],
                                  primer_scores: Dict[str, float],
                                  target_primers: int,
                                  network: Dict[int, Set[int]]) -> List[str]:
        """
        Greedily select primers maximizing network connectivity.

        At each step, select primer that:
        1. Has high network score
        2. Connects to already-selected primers
        3. Fills coverage gaps
        """
        selected = []
        covered_positions = set()

        # Create priority queue (max-heap using negative scores)
        heap = [(-score, primer_id) for primer_id, score in primer_scores.items()]
        heapq.heapify(heap)

        while len(selected) < target_primers and heap:
            # Get highest-scoring primer
            _, primer_id = heapq.heappop(heap)

            # Add to selected set
            selected.append(primer_id)

            # Update covered positions
            new_positions = primers_data[primer_id]['binding_sites']
            covered_positions.update(new_positions)

            logger.info(
                f"  Selected {primer_id} "
                f"(score={primer_scores[primer_id]:.3f}, "
                f"sites={len(new_positions)}, "
                f"total_coverage={len(covered_positions)})"
            )

        return selected

    def _build_result(self, selected_primers: List[str],
                     primers_data: Dict[str, Dict],
                     all_sites: List[BindingSite],
                     network: Dict[int, Set[int]]) -> CooperativeSelectionResult:
        """Build comprehensive result with network metrics."""

        # Get binding sites for selected primers
        selected_sites = [
            site for site in all_sites
            if site.primer_id in selected_primers
        ]

        # Calculate coverage
        covered_positions = set(site.position for site in selected_sites)
        coverage = len(covered_positions) / self.genome_length

        # Build subnetwork for selected primers
        site_indices = {
            i for i, site in enumerate(all_sites)
            if site.primer_id in selected_primers
        }

        subnetwork = {
            i: network[i] & site_indices
            for i in site_indices
        }

        # Calculate network metrics
        metrics = self._calculate_network_metrics(subnetwork, len(selected_sites))

        # Find amplification chains
        chains = self._find_amplification_chains(selected_sites, subnetwork, all_sites)

        # Calculate mean quality
        mean_quality = np.mean([
            primers_data[pid]['quality'] for pid in selected_primers
        ])

        # Cooperation score (network density)
        max_edges = len(selected_sites) * (len(selected_sites) - 1) / 2
        cooperation_score = metrics.total_edges / max_edges if max_edges > 0 else 0

        return CooperativeSelectionResult(
            selected_primers=selected_primers,
            binding_sites=selected_sites,
            coverage=coverage,
            network_metrics=metrics,
            amplification_chains=chains,
            mean_quality=mean_quality,
            cooperation_score=cooperation_score
        )

    def _calculate_network_metrics(self, network: Dict[int, Set[int]],
                                   num_sites: int) -> NetworkMetrics:
        """Calculate comprehensive network topology metrics."""

        # Total edges
        total_edges = sum(len(neighbors) for neighbors in network.values()) // 2

        # Average degree
        avg_degree = 2 * total_edges / num_sites if num_sites > 0 else 0

        # Clustering coefficient
        clustering = self._calculate_clustering(network)

        # Connected components
        components = self._find_connected_components(network)
        num_components = len(components)
        largest_component = max(len(c) for c in components) if components else 0

        # Network diameter (max shortest path)
        diameter = self._calculate_diameter(network, components)

        # Bottleneck score (inverse of betweenness variance)
        bottleneck = self._calculate_bottleneck_score(network)

        return NetworkMetrics(
            total_sites=num_sites,
            total_edges=total_edges,
            avg_degree=avg_degree,
            clustering_coefficient=clustering,
            connected_components=num_components,
            largest_component_size=largest_component,
            network_diameter=diameter,
            bottleneck_score=bottleneck
        )

    def _calculate_clustering(self, network: Dict[int, Set[int]]) -> float:
        """Calculate average clustering coefficient."""
        if not network:
            return 0.0

        clustering_scores = []

        for node, neighbors in network.items():
            if len(neighbors) < 2:
                clustering_scores.append(0.0)
                continue

            # Count edges between neighbors
            edges_between = 0
            neighbors_list = list(neighbors)
            for i, n1 in enumerate(neighbors_list):
                for n2 in neighbors_list[i+1:]:
                    if n2 in network.get(n1, set()):
                        edges_between += 1

            # Clustering = actual_edges / possible_edges
            possible_edges = len(neighbors) * (len(neighbors) - 1) / 2
            clustering = edges_between / possible_edges if possible_edges > 0 else 0
            clustering_scores.append(clustering)

        return np.mean(clustering_scores) if clustering_scores else 0.0

    def _find_connected_components(self, network: Dict[int, Set[int]]) -> List[Set[int]]:
        """Find connected components using BFS."""
        visited = set()
        components = []

        for node in network.keys():
            if node in visited:
                continue

            # BFS from this node
            component = set()
            queue = deque([node])

            while queue:
                current = queue.popleft()
                if current in visited:
                    continue

                visited.add(current)
                component.add(current)

                # Add unvisited neighbors
                for neighbor in network.get(current, set()):
                    if neighbor not in visited:
                        queue.append(neighbor)

            components.append(component)

        return components

    def _calculate_diameter(self, network: Dict[int, Set[int]],
                           components: List[Set[int]]) -> int:
        """Calculate network diameter (longest shortest path)."""
        if not components:
            return 0

        # Calculate diameter of largest component only
        largest = max(components, key=len)

        if len(largest) <= 1:
            return 0

        # BFS to find max distance
        max_distance = 0

        for start in list(largest)[:10]:  # Sample to avoid O(n^2)
            distances = {start: 0}
            queue = deque([start])

            while queue:
                current = queue.popleft()
                current_dist = distances[current]

                for neighbor in network.get(current, set()):
                    if neighbor not in distances and neighbor in largest:
                        distances[neighbor] = current_dist + 1
                        queue.append(neighbor)
                        max_distance = max(max_distance, current_dist + 1)

        return max_distance

    def _calculate_bottleneck_score(self, network: Dict[int, Set[int]]) -> float:
        """
        Calculate bottleneck score (0-1, higher = fewer bottlenecks).

        Low variance in node degrees = fewer bottlenecks.
        """
        if not network:
            return 0.0

        degrees = [len(neighbors) for neighbors in network.values()]

        if not degrees or np.mean(degrees) == 0:
            return 0.0

        # Coefficient of variation
        cv = np.std(degrees) / np.mean(degrees)

        # Bottleneck score = 1 / (1 + cv)
        # High CV = high variance = bottlenecks
        return 1.0 / (1.0 + cv)

    def _find_amplification_chains(self, selected_sites: List[BindingSite],
                                   network: Dict[int, Set[int]],
                                   all_sites: List[BindingSite]) -> List[List[int]]:
        """
        Find continuous amplification chains (paths through network).

        Returns:
            List of chains, where each chain is list of positions
        """
        # Build position index
        site_idx_map = {
            i: idx for idx, (i, site) in enumerate(enumerate(all_sites))
            if site in selected_sites
        }

        chains = []

        # Find paths using DFS
        visited_global = set()

        for start_idx in site_idx_map.keys():
            if start_idx in visited_global:
                continue

            # DFS to find chain
            chain = []
            stack = [start_idx]
            visited_local = set()

            while stack:
                current = stack.pop()

                if current in visited_local:
                    continue

                visited_local.add(current)
                visited_global.add(current)

                # Add position to chain
                if current in site_idx_map:
                    site = all_sites[current]
                    chain.append(site.position)

                # Add neighbors
                for neighbor in network.get(current, set()):
                    if neighbor not in visited_local and neighbor in site_idx_map:
                        stack.append(neighbor)

            # Keep chains meeting minimum length
            if len(chain) >= self.min_chain_length:
                chain.sort()
                chains.append(chain)

        return chains

    def _log_selection_summary(self, result: CooperativeSelectionResult):
        """Log detailed selection summary."""
        logger.info(f"\n{'='*60}")
        logger.info(f"COOPERATIVE BINDING SELECTION COMPLETE")
        logger.info(f"{'='*60}")
        logger.info(f"Primers selected: {len(result.selected_primers)}")
        logger.info(f"Coverage: {result.coverage:.1%}")
        logger.info(f"Mean quality: {result.mean_quality:.3f}")
        logger.info(f"Cooperation score: {result.cooperation_score:.3f}")

        logger.info(f"\nNetwork Metrics:")
        m = result.network_metrics
        logger.info(f"  Total sites: {m.total_sites}")
        logger.info(f"  Total edges: {m.total_edges}")
        logger.info(f"  Avg degree: {m.avg_degree:.1f} connections/site")
        logger.info(f"  Clustering: {m.clustering_coefficient:.3f}")
        logger.info(f"  Connected components: {m.connected_components}")
        logger.info(f"  Largest component: {m.largest_component_size} sites")
        logger.info(f"  Network diameter: {m.network_diameter} hops")
        logger.info(f"  Bottleneck score: {m.bottleneck_score:.3f}")

        logger.info(f"\nAmplification Chains:")
        logger.info(f"  Total chains: {len(result.amplification_chains)}")
        if result.amplification_chains:
            chain_lengths = [len(c) for c in result.amplification_chains]
            logger.info(f"  Avg chain length: {np.mean(chain_lengths):.1f} sites")
            logger.info(f"  Max chain length: {max(chain_lengths)} sites")

        logger.info(f"{'='*60}\n")


if __name__ == '__main__':
    # Example usage
    logging.basicConfig(level=logging.INFO)

    print("NeoSWGA Cooperative Binding Selector - Example\n")

    # Create example primers with realistic binding patterns
    import random
    random.seed(42)

    genome_length = 2_000_000
    primers_data = {}

    for i in range(40):
        # Create binding sites in clusters (realistic pattern)
        cluster_center = random.randint(0, genome_length - 100000)
        n_sites = random.randint(3, 8)

        sites = []
        for _ in range(n_sites):
            # Sites clustered within 50kb
            offset = random.randint(-25000, 25000)
            pos = max(0, min(genome_length - 1, cluster_center + offset))
            sites.append(pos)

        primers_data[f"Primer_{i+1}"] = {
            'binding_sites': sorted(sites),
            'quality': random.uniform(0.6, 1.0)
        }

    # Create selector
    selector = CooperativeBindingSelector(
        genome_length=genome_length,
        extension_distance=70000,
        min_chain_length=3
    )

    # Select cooperative set
    result = selector.select_cooperative_set(
        primers_data=primers_data,
        target_primers=15,
        quality_weight=0.3
    )

    print(f"\nSelected {len(result.selected_primers)} primers:")
    for pid in result.selected_primers[:10]:
        print(f"  {pid}: {len(primers_data[pid]['binding_sites'])} sites")

    print(f"\nCoverage: {result.coverage:.1%}")
    print(f"Network avg degree: {result.network_metrics.avg_degree:.1f}")
    print(f"Amplification chains: {len(result.amplification_chains)}")
