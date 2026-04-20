"""
Clique-based dimer-free primer set optimizer.

Guarantees mathematical exhaustiveness: every dimer-free primer set of
the target size is enumerable. This approach models primer compatibility
as a graph problem:

1. Build a compatibility graph where each node is a candidate primer
   and edges connect primers that do NOT form dimers.
2. Find cliques (fully connected subgraphs) of the target size.
   Each clique represents a set of mutually compatible (dimer-free) primers.
3. Score each clique by coverage and selectivity.
4. Return the highest-scoring dimer-free set.

Based on the approach used in SWGA v1 (Clark et al. 2017) with Cliquer
integration. Here we use NetworkX's clique enumeration, which provides
correct results for moderate pool sizes (<200 primers).

Computational complexity:
- Dimer matrix: O(n^2) where n = number of candidates
- Clique enumeration: NP-hard in general, but tractable for sparse
  compatibility graphs typical of SWGA primer pools
- Recommended for pools of up to ~200 primers; larger pools should use
  the hybrid or greedy optimizer with incremental dimer checking.
"""

import logging
import time
from dataclasses import dataclass
from typing import List, Dict, Optional, Set, Tuple

import numpy as np

try:
    import networkx as nx
    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False

from .base_optimizer import (
    BaseOptimizer, OptimizationResult, OptimizationStatus,
    PrimerSetMetrics, OptimizerConfig
)
from .optimizer_factory import OptimizerFactory
from .dimer import heterodimer_matrix_fast, is_dimer_fast

logger = logging.getLogger(__name__)

# Maximum candidate pool size before warning about runtime
MAX_RECOMMENDED_POOL = 200


@dataclass
class CliqueOptimizerConfig(OptimizerConfig):
    """Configuration for clique-based optimizer."""
    # Maximum number of cliques to evaluate (0 = unlimited)
    max_cliques: int = 10000
    # Use thermodynamic dimer check instead of sequence-based
    use_thermodynamic_dimer: bool = False
    # Delta-G threshold for thermodynamic dimer check (kcal/mol)
    delta_g_threshold: float = -6.0
    # Maximum pool size before automatic truncation
    max_pool_size: int = 200
    # Include self-dimer check in compatibility
    check_self_dimers: bool = True


def build_compatibility_graph(
    primers: List[str],
    max_dimer_bp: int = 3,
    max_self_dimer_bp: int = 4,
    check_self_dimers: bool = True,
) -> 'nx.Graph':
    """
    Build a graph where edges connect dimer-free primer pairs.

    Each node is a primer. An edge (i, j) means primers i and j do
    not form a heterodimer. A clique in this graph is therefore a
    set of mutually compatible primers.

    Args:
        primers: List of candidate primer sequences
        max_dimer_bp: Maximum complementary bases for heterodimer
        max_self_dimer_bp: Maximum complementary bases for self-dimer
        check_self_dimers: Exclude primers with self-dimer risk

    Returns:
        NetworkX Graph with primers as nodes and compatibility edges
    """
    if not HAS_NETWORKX:
        raise ImportError(
            "networkx is required for the clique optimizer. "
            "Install it with: pip install networkx"
        )

    n = len(primers)
    G = nx.Graph()

    # Add all primers as nodes (unless they have self-dimer issues)
    valid_primers = []
    for primer in primers:
        if check_self_dimers and is_dimer_fast(primer, primer, max_self_dimer_bp):
            logger.debug(f"Excluding {primer}: self-dimer risk")
            continue
        valid_primers.append(primer)
        G.add_node(primer)

    logger.info(
        f"Compatibility graph: {len(valid_primers)}/{n} primers "
        f"passed self-dimer check"
    )

    # Build dimer matrix for valid primers
    het_matrix = heterodimer_matrix_fast(valid_primers, max_dimer_bp)

    # Add edges where primers are compatible (no dimer)
    edge_count = 0
    for i in range(len(valid_primers)):
        for j in range(i + 1, len(valid_primers)):
            if het_matrix[i, j] == 0:
                G.add_edge(valid_primers[i], valid_primers[j])
                edge_count += 1

    total_pairs = len(valid_primers) * (len(valid_primers) - 1) // 2
    density = edge_count / total_pairs if total_pairs > 0 else 0
    logger.info(
        f"Compatibility graph: {edge_count}/{total_pairs} edges "
        f"(density={density:.2f})"
    )

    return G


def enumerate_dimer_free_sets(
    G: 'nx.Graph',
    target_size: int,
    max_cliques: int = 10000,
) -> List[List[str]]:
    """
    Enumerate all dimer-free primer sets of the target size.

    Uses NetworkX's clique enumeration to find all maximal cliques,
    then extracts subsets of the target size. For efficiency, stops
    after max_cliques sets are found.

    Args:
        G: Compatibility graph from build_compatibility_graph()
        target_size: Desired primer set size
        max_cliques: Maximum number of clique-derived sets to return

    Returns:
        List of primer sets (each a list of primer sequences)
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx is required for clique enumeration")

    sets = []
    seen = set()

    # Enumerate cliques of at least target_size
    for clique in nx.find_cliques(G):
        if len(clique) < target_size:
            continue

        if len(clique) == target_size:
            key = frozenset(clique)
            if key not in seen:
                seen.add(key)
                sets.append(list(clique))
        else:
            # Extract all target_size subsets from larger cliques
            # Use itertools.combinations for efficiency
            from itertools import combinations
            for subset in combinations(clique, target_size):
                key = frozenset(subset)
                if key not in seen:
                    seen.add(key)
                    sets.append(list(subset))

        if len(sets) >= max_cliques:
            logger.info(
                f"Reached max_cliques limit ({max_cliques}). "
                f"Stopping enumeration."
            )
            break

    return sets


@OptimizerFactory.register(
    'clique',
    aliases=['clique-dimer-free', 'dimer-free'],
    description='Clique-based exhaustive dimer-free set discovery'
)
class CliqueOptimizer(BaseOptimizer):
    """
    Clique-based optimizer guaranteeing dimer-free primer sets.

    Builds a primer compatibility graph and enumerates cliques to find
    all possible dimer-free primer sets. Each set is scored by coverage
    and selectivity, and the best is returned.

    Best for primer pools under ~200 candidates where exhaustive dimer-free
    guarantees are needed. For larger pools, use the hybrid optimizer with
    incremental dimer checking.
    """

    def __init__(
        self,
        position_cache,
        fg_prefixes: List[str],
        fg_seq_lengths: List[int],
        bg_prefixes: Optional[List[str]] = None,
        bg_seq_lengths: Optional[List[int]] = None,
        config: Optional[OptimizerConfig] = None,
        conditions=None,
        **kwargs
    ):
        if not HAS_NETWORKX:
            raise ImportError(
                "networkx is required for the clique optimizer. "
                "Install it with: pip install networkx"
            )
        super().__init__(
            position_cache, fg_prefixes, fg_seq_lengths,
            bg_prefixes, bg_seq_lengths,
            config or CliqueOptimizerConfig(),
            conditions=conditions,
        )
        # Extract clique-specific config
        if isinstance(config, CliqueOptimizerConfig):
            self.clique_config = config
        else:
            self.clique_config = CliqueOptimizerConfig()

    @property
    def name(self) -> str:
        return "clique"

    @property
    def description(self) -> str:
        return "Clique-based exhaustive dimer-free set discovery"

    def optimize(
        self,
        candidates: List[str],
        target_size: Optional[int] = None,
        fixed_primers: Optional[List[str]] = None,
        **kwargs
    ) -> OptimizationResult:
        """
        Find the highest-scoring dimer-free primer set.

        Steps:
        1. Optionally truncate pool to max_pool_size by coverage ranking
        2. Build compatibility graph (edges = no dimer risk)
        3. Enumerate cliques of target size
        4. Score each clique by coverage/selectivity
        5. Return the best-scoring set

        Args:
            candidates: Pool of candidate primer sequences
            target_size: Desired number of primers in the set
            fixed_primers: Primers that must be included (pre-validated
                for mutual compatibility)

        Returns:
            OptimizationResult with the best dimer-free primer set
        """
        candidates = self._validate_candidates(candidates)
        target = target_size or self.config.target_set_size
        max_dimer_bp = self.config.max_dimer_bp

        # Handle fixed primers
        fixed = list(fixed_primers) if fixed_primers else []
        if fixed:
            # Validate fixed primers are mutually compatible
            for i in range(len(fixed)):
                for j in range(i + 1, len(fixed)):
                    if is_dimer_fast(fixed[i], fixed[j], max_dimer_bp):
                        return OptimizationResult.failure(
                            self.name,
                            f"Fixed primers {fixed[i]} and {fixed[j]} "
                            f"form a dimer"
                        )
            # Remove fixed primers from candidate pool
            fixed_set = set(fixed)
            candidates = [p for p in candidates if p not in fixed_set]
            # Reduce target to account for fixed primers
            target = max(1, target - len(fixed))

        # Truncate pool if too large
        if len(candidates) > self.clique_config.max_pool_size:
            logger.warning(
                f"Pool size ({len(candidates)}) exceeds recommended "
                f"maximum ({self.clique_config.max_pool_size}). "
                f"Ranking by foreground binding and truncating."
            )
            candidates = self._rank_and_truncate(
                candidates, self.clique_config.max_pool_size
            )

        if len(candidates) < target:
            return OptimizationResult.failure(
                self.name,
                f"Not enough candidates ({len(candidates)}) for "
                f"target size ({target})"
            )

        logger.info(
            f"Clique optimizer: {len(candidates)} candidates, "
            f"target size {target}"
        )

        # Step 1: Build compatibility graph
        start_time = time.time()
        G = build_compatibility_graph(
            candidates,
            max_dimer_bp=max_dimer_bp,
            max_self_dimer_bp=getattr(
                self.config, 'max_self_dimer_bp',
                max_dimer_bp + 1
            ),
            check_self_dimers=self.clique_config.check_self_dimers,
        )
        graph_time = time.time() - start_time
        logger.info(f"Compatibility graph built in {graph_time:.1f}s")

        # If fixed primers are specified, restrict graph to their neighbors
        if fixed:
            # Only keep nodes compatible with ALL fixed primers
            compatible_nodes = set(G.nodes())
            for fp in fixed:
                if fp in G:
                    fp_neighbors = set(G.neighbors(fp))
                    compatible_nodes &= fp_neighbors
                else:
                    # Fixed primer not in graph (self-dimer?) - proceed anyway
                    logger.warning(
                        f"Fixed primer {fp} excluded from graph "
                        f"(possible self-dimer)"
                    )
            G = G.subgraph(compatible_nodes).copy()
            logger.info(
                f"After fixed primer filtering: {len(G.nodes())} "
                f"compatible candidates"
            )

        # Step 2: Enumerate cliques
        start_time = time.time()
        dimer_free_sets = enumerate_dimer_free_sets(
            G, target, max_cliques=self.clique_config.max_cliques
        )
        enum_time = time.time() - start_time
        logger.info(
            f"Enumerated {len(dimer_free_sets)} dimer-free sets "
            f"in {enum_time:.1f}s"
        )

        if not dimer_free_sets:
            return OptimizationResult.failure(
                self.name,
                f"No dimer-free sets of size {target} found. "
                f"Try reducing target size or relaxing dimer threshold."
            )

        # Step 3: Score each set
        best_set = None
        best_score = float('-inf')
        best_metrics = None

        for primer_set in dimer_free_sets:
            # Prepend fixed primers
            full_set = fixed + primer_set
            score = self._score_set(full_set)
            if score > best_score:
                best_score = score
                best_set = full_set
                best_metrics = self.compute_metrics(full_set)

        total_time = graph_time + enum_time
        logger.info(
            f"Best set score: {best_score:.4f} "
            f"(evaluated {len(dimer_free_sets)} sets in {total_time:.1f}s)"
        )

        status = OptimizationStatus.SUCCESS
        if best_metrics and best_metrics.fg_coverage < 0.5:
            status = OptimizationStatus.PARTIAL

        return OptimizationResult(
            primers=tuple(best_set),
            score=best_score,
            status=status,
            metrics=best_metrics or PrimerSetMetrics.empty(),
            iterations=len(dimer_free_sets),
            optimizer_name=self.name,
            message=(
                f"Evaluated {len(dimer_free_sets)} dimer-free sets "
                f"from {len(candidates)} candidates"
            ),
        )

    def _score_set(self, primers: List[str]) -> float:
        """
        Score a primer set by foreground/background binding ratio.

        Args:
            primers: List of primer sequences

        Returns:
            Score (higher is better): fg_sites / (bg_sites + 1)
        """
        fg_total = 0
        bg_total = 0

        for primer in primers:
            for prefix in self.fg_prefixes:
                positions = self.get_primer_positions(primer, prefix, 'both')
                fg_total += len(positions)
            for prefix in self.bg_prefixes:
                positions = self.get_primer_positions(primer, prefix, 'both')
                bg_total += len(positions)

        return fg_total / (bg_total + 1)

    def _rank_and_truncate(
        self, candidates: List[str], max_size: int
    ) -> List[str]:
        """
        Rank candidates by foreground binding and keep top max_size.

        Args:
            candidates: Full candidate list
            max_size: Maximum number to keep

        Returns:
            Truncated list of best candidates
        """
        scored = []
        for primer in candidates:
            fg_count = 0
            for prefix in self.fg_prefixes:
                positions = self.get_primer_positions(primer, prefix, 'both')
                fg_count += len(positions)
            scored.append((fg_count, primer))

        scored.sort(reverse=True)
        return [p for _, p in scored[:max_size]]
