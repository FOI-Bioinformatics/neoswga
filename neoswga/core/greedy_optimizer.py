"""
Greedy breadth-first search optimizer for primer set selection.

Refactored implementation of the original optimize.py algorithm using:
- BaseOptimizer interface for consistent API
- Typed dataclasses instead of **kwargs
- Clear separation of concerns
- Proper error handling and logging

The algorithm builds primer sets incrementally:
1. Start with empty sets (or initial sets if provided)
2. For each iteration:
   a. Try adding each candidate primer to each current set
   b. Score all resulting sets
   c. Keep top-k sets for next iteration
3. Return best set found

This greedy approach provides good results quickly but may miss
globally optimal solutions due to local optima.
"""

import numpy as np
import logging
from collections import OrderedDict
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional, Set, FrozenSet
from functools import partial
from concurrent.futures import ProcessPoolExecutor

from .base_optimizer import (
    BaseOptimizer, OptimizationResult, OptimizationStatus,
    PrimerSetMetrics, OptimizerConfig
)
from .optimizer_factory import OptimizerFactory
from .search_context import (
    BFSSearchContext, BFSConfig, GenomeInfo, DimerConstraints,
    SearchState, PositionData, SearchResult
)
from .exceptions import NoCandidatesError, OptimizerConvergenceError

logger = logging.getLogger(__name__)


class BoundedScoreCache(OrderedDict):
    """
    LRU cache with maximum size for primer set scores.

    Prevents unbounded memory growth during optimization.
    Uses OrderedDict for O(1) insertion and eviction.
    """

    def __init__(self, max_size: int = 100000):
        super().__init__()
        self.max_size = max_size
        self.hits = 0
        self.misses = 0

    def get(self, key: FrozenSet[str], default=None):
        """Get value, moving key to end (most recently used)."""
        if key in self:
            self.hits += 1
            self.move_to_end(key)
            return self[key]
        self.misses += 1
        return default

    def set(self, key: FrozenSet[str], value: float) -> None:
        """Set value, evicting oldest if at capacity."""
        if key in self:
            self.move_to_end(key)
        else:
            if len(self) >= self.max_size:
                self.popitem(last=False)  # Remove oldest (LRU)
        self[key] = value

    def hit_rate(self) -> float:
        """Return cache hit rate."""
        total = self.hits + self.misses
        return self.hits / total if total > 0 else 0.0

    def clear_stats(self) -> None:
        """Reset hit/miss counters."""
        self.hits = 0
        self.misses = 0


@dataclass
class GreedyConfig(OptimizerConfig):
    """
    Configuration for greedy optimizer.

    Extends base config with greedy-specific parameters.
    """
    max_sets: int = 10  # Number of parallel sets to maintain
    selection_method: str = 'deterministic'  # 'deterministic', 'softmax', 'normalized'
    drop_out_iteration: int = 4  # Iteration for dropout layer
    enable_dropout: bool = True  # Whether to use dropout for escaping local optima
    score_cache_size: int = 100000  # Maximum entries in score cache


@OptimizerFactory.register('greedy', aliases=['bfs', 'greedy-bfs', 'greedy-dropout', 'dropout'])
class GreedyOptimizer(BaseOptimizer):
    """
    Greedy breadth-first search optimizer.

    Maintains multiple candidate sets in parallel and greedily adds
    the best primer to each set at every iteration. Uses score caching
    to avoid redundant evaluations.

    Advantages:
    - Fast convergence for well-separated primer pools
    - Deterministic results (with same inputs)
    - Low memory usage

    Disadvantages:
    - Can get stuck in local optima
    - May miss synergistic primer combinations
    - No guarantee of global optimum

    Use when:
    - You need fast results
    - Primer pool is well-curated
    - Approximate solutions are acceptable
    """

    def __init__(
        self,
        position_cache,
        fg_prefixes: List[str],
        fg_seq_lengths: List[int],
        bg_prefixes: Optional[List[str]] = None,
        bg_seq_lengths: Optional[List[int]] = None,
        config: Optional[GreedyConfig] = None,
        dimer_matrix: Optional[np.ndarray] = None,
        **kwargs
    ):
        """
        Initialize greedy optimizer.

        Args:
            position_cache: PositionCache for primer position lookups
            fg_prefixes: Foreground genome HDF5 prefixes
            fg_seq_lengths: Foreground genome lengths
            bg_prefixes: Background genome HDF5 prefixes
            bg_seq_lengths: Background genome lengths
            config: Greedy-specific configuration
            dimer_matrix: Pre-computed dimer compatibility matrix
        """
        super().__init__(
            position_cache, fg_prefixes, fg_seq_lengths,
            bg_prefixes, bg_seq_lengths, config or GreedyConfig()
        )
        self.greedy_config = config if isinstance(config, GreedyConfig) else GreedyConfig()
        self.dimer_matrix = dimer_matrix
        # Use bounded LRU cache to prevent memory leaks
        self._score_cache = BoundedScoreCache(
            max_size=self.greedy_config.score_cache_size
        )

    @property
    def name(self) -> str:
        """Optimizer identifier for logging and factory registration."""
        return "greedy"

    @property
    def description(self) -> str:
        """One-line summary of the optimization strategy."""
        return "Greedy breadth-first search optimizer with parallel set tracking"

    def optimize(
        self,
        candidates: List[str],
        target_size: Optional[int] = None,
        initial_sets: Optional[List[List[str]]] = None,
        banned_primers: Optional[Set[str]] = None,
        **kwargs
    ) -> OptimizationResult:
        """
        Find optimal primer set using greedy BFS.

        Args:
            candidates: Pool of candidate primers
            target_size: Desired primer set size (default from config)
            initial_sets: Optional starting sets to build from
            banned_primers: Primers to exclude from search

        Returns:
            OptimizationResult with best primer set found
        """
        # Validate and prepare
        candidates = self._validate_candidates(candidates)
        target = target_size or self.config.target_set_size
        banned = banned_primers or set()

        if len(candidates) == 0:
            return OptimizationResult.failure(self.name, "No valid candidates")

        if self.greedy_config.verbose:
            logger.info(f"Starting greedy optimization: {len(candidates)} candidates, target={target}")

        # Build dimer constraints if not provided
        dimer_constraints = self._build_dimer_constraints(candidates)

        # Create search context
        fg_genome = GenomeInfo.from_lists(
            self.fg_prefixes, self.fg_seq_lengths, circular=True
        )
        bg_genome = GenomeInfo.from_lists(
            self.bg_prefixes, self.bg_seq_lengths, circular=False
        ) if self.bg_prefixes else None

        bfs_config = BFSConfig(
            max_sets=self.greedy_config.max_sets,
            iterations=self.config.max_iterations,
            selection_method=self.greedy_config.selection_method,
            verbose=self.config.verbose,
        )

        ctx = BFSSearchContext(
            primer_pool=candidates,
            fg_genome=fg_genome,
            bg_genome=bg_genome,
            dimer_constraints=dimer_constraints,
            config=bfs_config,
            banned_primers=banned,
        )

        # Initialize states
        ctx.initialize_states(initial_sets)

        # Run BFS iterations
        try:
            result = self._run_bfs(ctx, target)
        except Exception as e:
            logger.error(f"BFS optimization failed: {e}")
            return OptimizationResult.failure(self.name, str(e))

        # Convert to OptimizationResult
        best_primers = result.best_primers
        metrics = self.compute_metrics(best_primers)

        status = OptimizationStatus.SUCCESS
        if len(best_primers) < target:
            status = OptimizationStatus.PARTIAL

        return OptimizationResult(
            primers=tuple(best_primers),
            score=result.best_score,
            status=status,
            metrics=metrics,
            iterations=result.iterations_completed,
            optimizer_name=self.name,
            message=result.message,
        )

    def _build_dimer_constraints(self, candidates: List[str]) -> DimerConstraints:
        """Build dimer constraints from candidates."""
        if self.dimer_matrix is not None:
            primer_to_idx = {p: i for i, p in enumerate(candidates)}
            return DimerConstraints(
                matrix=self.dimer_matrix,
                primer_to_index=primer_to_idx,
                max_dimer_bp=self.config.max_dimer_bp,
            )

        # Build dimer matrix using optimized function
        from . import dimer

        n = len(candidates)
        if self.config.verbose:
            logger.info(f"Computing dimer matrix for {n} candidates...")

        # Use fast sequential for small sets, parallel for large sets
        if n > 200:
            # Parallel computation for large primer pools
            matrix = dimer.heterodimer_matrix_parallel(
                candidates,
                max_dimer_bp=self.config.max_dimer_bp
            )
        else:
            # Fast sequential with early termination
            matrix = dimer.heterodimer_matrix_fast(
                candidates,
                max_dimer_bp=self.config.max_dimer_bp
            )

        # Convert to bool
        matrix = matrix.astype(bool)

        primer_to_idx = {p: i for i, p in enumerate(candidates)}
        return DimerConstraints(
            matrix=matrix,
            primer_to_index=primer_to_idx,
            max_dimer_bp=self.config.max_dimer_bp,
        )

    def _run_bfs(self, ctx: BFSSearchContext, target_size: int) -> SearchResult:
        """
        Run BFS optimization iterations.

        Args:
            ctx: Search context with all state
            target_size: Target primer set size

        Returns:
            SearchResult with best states found
        """
        best_score = float('-inf')
        finished_states = []
        iteration = 0  # Initialize to handle case when iterations is 0

        for iteration in range(1, ctx.config.iterations + 1):
            if ctx.config.verbose:
                logger.info(f"Iteration {iteration}/{ctx.config.iterations}")

            # Expand each current state
            all_new_states = []
            for state in ctx.top_states:
                if len(state.primers) >= target_size:
                    # Already at target size, keep as-is
                    finished_states.append(state)
                    continue

                # Get available primers for this state
                available = ctx.get_available_primers(state.primers)
                if not available:
                    finished_states.append(state)
                    continue

                # Try adding each available primer
                new_states = self._expand_state(ctx, state, available)
                all_new_states.extend(new_states)

            if not all_new_states:
                # No more expansion possible
                break

            # Select top states for next iteration
            all_new_states.sort(key=lambda s: s.score, reverse=True)
            ctx.top_states = all_new_states[:ctx.config.max_sets]

            # Only check convergence after sets have reached target size
            # (don't stop building sets just because score plateaued temporarily)
            best_set_size = max(len(s.primers) for s in ctx.top_states)
            if best_set_size >= target_size and self._check_convergence(ctx.top_states, best_score):
                if ctx.config.verbose:
                    logger.info("Converged - no improvement at target size")
                break

            # Track best (after convergence check)
            if ctx.top_states[0].score > best_score:
                best_score = ctx.top_states[0].score

        # Combine finished and current states
        all_states = finished_states + ctx.top_states
        all_states.sort(key=lambda s: s.score, reverse=True)

        base_result = SearchResult(
            states=all_states,
            best_score=best_score,
            iterations_completed=iteration,
            converged=True,
            message=f"Completed {iteration} iterations",
        )

        # Apply dropout if enabled and we have enough iterations
        if (self.greedy_config.enable_dropout and
                base_result.iterations_completed >= self.greedy_config.drop_out_iteration):

            if ctx.config.verbose:
                logger.info("Applying dropout layer")

            dropped_states = self._apply_dropout(ctx, base_result.states)
            combined = base_result.states + dropped_states
            combined.sort(key=lambda s: s.score, reverse=True)

            return SearchResult(
                states=combined[:ctx.config.max_sets],
                best_score=max(s.score for s in combined),
                iterations_completed=base_result.iterations_completed,
                converged=base_result.converged,
                message=base_result.message + " (with dropout)",
            )

        return base_result

    def _expand_state(
        self,
        ctx: BFSSearchContext,
        state: SearchState,
        available: List[str]
    ) -> List[SearchState]:
        """
        Expand a state by trying to add each available primer.

        Args:
            ctx: Search context
            state: Current state to expand
            available: Primers that can be added

        Returns:
            List of new states, sorted by score
        """
        new_states = []

        for primer in available:
            # Use frozenset for O(1) hashing instead of O(n log n) sort+join
            new_primers_list = state.primers + [primer]
            cache_key = frozenset(new_primers_list)

            cached_score = self._score_cache.get(cache_key)
            if cached_score is not None:
                score = cached_score
            else:
                # Compute new positions and score
                fg_positions = self._get_combined_positions(
                    state.fg_positions, primer, ctx.fg_genome
                )
                bg_positions = self._get_combined_positions(
                    state.bg_positions, primer, ctx.bg_genome
                ) if ctx.has_background else PositionData.empty(0)

                score = self._evaluate(fg_positions, bg_positions, ctx)
                self._score_cache.set(cache_key, score)

            # Create new state (keep sorted for consistent output)
            new_state = SearchState(
                primers=sorted(new_primers_list),
                score=score,
                fg_positions=state.fg_positions.copy(),
                bg_positions=state.bg_positions.copy(),
            )
            new_states.append(new_state)

        # Sort by score descending
        new_states.sort(key=lambda s: s.score, reverse=True)

        # Return top 2 for this state
        return new_states[:2]

    def _get_combined_positions(
        self,
        current: PositionData,
        primer: str,
        genome: Optional[GenomeInfo]
    ) -> PositionData:
        """Get positions combining current with new primer."""
        if genome is None:
            return PositionData.empty(0)

        result = current.copy()

        for i, prefix in enumerate(genome.prefixes):
            fwd = self.get_primer_positions(primer, prefix, 'forward')
            rev = self.get_primer_positions(primer, prefix, 'reverse')
            result.add_positions(i, fwd, rev)

        return result

    def _evaluate(
        self,
        fg_positions: PositionData,
        bg_positions: PositionData,
        ctx: BFSSearchContext
    ) -> float:
        """
        Evaluate a primer set based on positions.

        Scoring balances foreground binding density with background avoidance.
        When background binding is zero or negligible, scoring focuses on
        foreground density to ensure differentiation between candidates
        (important for small genomes where selectivity is saturated).
        """
        # Calculate foreground coverage
        fg_sites = fg_positions.total_sites()
        fg_density = fg_sites / max(ctx.fg_genome.total_length, 1)

        if fg_sites == 0:
            return float('-inf')

        if not ctx.has_background:
            return np.log10(fg_density * 1e6)

        # With background: compute selectivity
        bg_sites = bg_positions.total_sites()
        bg_density = bg_sites / max(ctx.bg_genome.total_length, 1) if ctx.bg_genome else 0

        # When background binding is zero or negligible, focus on foreground
        # density to maintain differentiation between candidates
        if bg_density < 1e-8:
            # No background binding — score purely on foreground density
            # Add small bonus for having zero background
            return np.log10(fg_density * 1e6) + 5.0

        # Selectivity ratio
        selectivity = fg_density / bg_density

        # Combined score: selectivity + density (both in log space)
        score = np.log10(selectivity) + np.log10(fg_density * 1e6)

        return score

    def _check_convergence(
        self,
        states: List[SearchState],
        best_score: float
    ) -> bool:
        """Check if optimization has converged."""
        if not states:
            return True

        current_best = max(s.score for s in states)
        improvement = current_best - best_score

        return improvement < self.config.convergence_threshold

    def _apply_dropout(
        self,
        ctx: BFSSearchContext,
        states: List[SearchState]
    ) -> List[SearchState]:
        """
        Apply dropout by removing each primer from top states.

        This allows the algorithm to escape local optima by exploring
        subsets of existing solutions.
        """
        dropped_states = []

        for state in states[:ctx.config.max_sets]:
            if len(state.primers) < 2:
                continue

            for i, _ in enumerate(state.primers):
                # Create subset without primer i
                new_primers = state.primers[:i] + state.primers[i+1:]

                # Score the subset using frozenset key
                cache_key = frozenset(new_primers)
                cached_score = self._score_cache.get(cache_key)
                if cached_score is None:
                    # Heuristic penalty for dropped primer
                    new_score = state.score * 0.95
                    self._score_cache.set(cache_key, new_score)
                else:
                    new_score = cached_score

                dropped_states.append(SearchState(
                    primers=sorted(new_primers),
                    score=new_score,
                    fg_positions=PositionData.empty(ctx.fg_genome.num_sequences),
                    bg_positions=PositionData.empty(
                        ctx.bg_genome.num_sequences if ctx.bg_genome else 0
                    ),
                ))

        return dropped_states
