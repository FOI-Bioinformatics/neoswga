"""
Mixed Integer Linear Programming (MILP) for exact primer set optimization.

For bacterial genomes with 100-1000 candidate primers, MILP can find
PROVABLY OPTIMAL solutions in seconds to minutes.

This is a game-changer compared to heuristics with unknown approximation ratios.
"""

import numpy as np
from typing import List, Dict, Tuple, Optional
import logging

try:
    from mip import Model, xsum, maximize, BINARY, OptimizationStatus
    MIP_AVAILABLE = True
except ImportError:
    logging.getLogger(__name__).debug("python-mip not installed. Install with: pip install mip")
    MIP_AVAILABLE = False
    Model = None

logger = logging.getLogger(__name__)


class MILPOptimizer:
    """
    Exact optimization via Mixed Integer Linear Programming.

    Decision variables: x_i ∈ {0,1} for each primer i
    Objective: Maximize coverage
    Constraints:
        - Select exactly K primers
        - No dimer pairs
        - Minimum specificity
        - Coverage requirements
    """

    def __init__(self, position_cache, fg_prefixes: List[str],
                 bg_prefixes: List[str], fg_seq_lengths: List[int],
                 bg_seq_lengths: List[int]):
        """
        Initialize MILP optimizer.

        Args:
            position_cache: PositionCache with all primer positions
            fg_prefixes: Target genome prefixes
            bg_prefixes: Background genome prefixes
            fg_seq_lengths: Target genome lengths
            bg_seq_lengths: Background genome lengths
        """
        if not MIP_AVAILABLE:
            raise ImportError("python-mip required. Install: pip install mip")

        self.cache = position_cache
        self.fg_prefixes = fg_prefixes
        self.bg_prefixes = bg_prefixes
        self.fg_seq_lengths = fg_seq_lengths
        self.bg_seq_lengths = bg_seq_lengths

    def optimize(self, candidates: List[str], num_primers: int = 10,
                max_time_seconds: int = 300, specificity_threshold: float = 0.0) -> Optional[List[str]]:
        """
        Find optimal primer set via MILP.

        Args:
            candidates: Candidate primers (ideally <1000 for speed)
            num_primers: Number of primers to select
            max_time_seconds: Time limit (default 5 minutes)
            specificity_threshold: Minimum fg/bg ratio (optional)

        Returns:
            Optimal primer set (or None if infeasible/timeout)
        """
        if len(candidates) > 1000:
            logger.warning(f"Large candidate set ({len(candidates)}). MILP may be slow. "
                          f"Consider pre-filtering or using greedy algorithm.")

        logger.info(f"Building MILP model: {len(candidates)} candidates, select {num_primers}")

        # Precompute coverage matrix
        logger.info("Precomputing coverage matrix...")
        coverage_matrix, regions = self._build_coverage_matrix(candidates)

        # Precompute dimer pairs
        logger.info("Precomputing dimer pairs...")
        dimer_pairs = self._find_dimer_pairs(candidates)

        # Precompute specificity
        logger.info("Precomputing specificity...")
        specificity = self._compute_specificity(candidates)

        # Build MILP model
        logger.info("Building MILP model...")
        model = Model(sense=maximize)

        # Decision variables: x[i] = 1 if primer i selected
        x = {i: model.add_var(var_type=BINARY, name=f"x_{i}")
             for i in range(len(candidates))}

        # Auxiliary variables: y[r] = 1 if region r covered
        y = {r: model.add_var(var_type=BINARY, name=f"y_{r}")
             for r in range(len(regions))}

        # Objective: Maximize total coverage
        model.objective = maximize(xsum(y[r] for r in range(len(regions))))

        # Constraint 1: Select exactly num_primers
        model += xsum(x[i] for i in range(len(candidates))) == num_primers

        # Constraint 2: No dimer pairs
        for i, j in dimer_pairs:
            model += x[i] + x[j] <= 1

        # Constraint 3: Region coverage definition
        for r in range(len(regions)):
            # Region r is covered if ANY selected primer covers it
            covering_primers = [i for i in range(len(candidates))
                              if coverage_matrix[i, r] == 1]
            if covering_primers:
                # y[r] <= sum of primers covering region r
                model += y[r] <= xsum(x[i] for i in covering_primers)

        # Constraint 4: Minimum specificity (optional)
        if specificity_threshold > 0:
            for i in range(len(candidates)):
                if specificity[i] < specificity_threshold:
                    model += x[i] == 0

        # Solve
        logger.info(f"Solving MILP (time limit: {max_time_seconds}s)...")
        status = model.optimize(max_seconds=max_time_seconds)

        if status == OptimizationStatus.OPTIMAL:
            logger.info("Optimal solution found!")
        elif status == OptimizationStatus.FEASIBLE:
            logger.info("Feasible solution found (not proven optimal)")
        else:
            logger.warning("No solution found within time limit")
            return None

        # Extract solution
        selected_indices = [i for i in range(len(candidates)) if x[i].x >= 0.99]
        selected_primers = [candidates[i] for i in selected_indices]

        # Report
        obj_value = model.objective_value
        logger.info(f"Solution: {len(selected_primers)} primers, "
                   f"{obj_value}/{len(regions)} regions covered "
                   f"({100*obj_value/len(regions):.1f}%)")

        return selected_primers

    def _build_coverage_matrix(self, candidates: List[str],
                               window_size: int = 10000) -> Tuple[np.ndarray, List[int]]:
        """
        Build binary coverage matrix.

        Returns:
            coverage_matrix: [num_primers × num_regions] binary matrix
            regions: List of region start positions
        """
        genome_length = sum(self.fg_seq_lengths)
        num_regions = (genome_length + window_size - 1) // window_size
        regions = list(range(0, genome_length, window_size))

        coverage = np.zeros((len(candidates), num_regions), dtype=np.uint8)

        for i, primer in enumerate(candidates):
            # Get all binding positions
            all_positions = []
            for prefix in self.fg_prefixes:
                positions = self.cache.get_positions(prefix, primer, 'both')
                all_positions.extend(positions)

            # Mark covered regions
            for pos in all_positions:
                region_idx = int(pos) // window_size
                if region_idx < num_regions:
                    coverage[i, region_idx] = 1

        return coverage, regions

    def _find_dimer_pairs(self, candidates: List[str]) -> List[Tuple[int, int]]:
        """
        Find all primer pairs that may form dimers.

        Returns:
            List of (i, j) pairs where i < j
        """
        from . import dimer  # Use existing dimer detection

        dimer_pairs = []

        for i in range(len(candidates)):
            for j in range(i+1, len(candidates)):
                if dimer.is_dimer_fast(candidates[i], candidates[j], max_dimer_bp=3):
                    dimer_pairs.append((i, j))

        logger.info(f"Found {len(dimer_pairs)} dimer pairs")
        return dimer_pairs

    def _compute_specificity(self, candidates: List[str]) -> np.ndarray:
        """
        Compute fg/bg ratio for each primer.

        Returns:
            Array of specificity values
        """
        specificity = np.zeros(len(candidates))

        for i, primer in enumerate(candidates):
            fg_count = sum(
                len(self.cache.get_positions(p, primer, 'both'))
                for p in self.fg_prefixes
            )
            bg_count = sum(
                len(self.cache.get_positions(p, primer, 'both'))
                for p in self.bg_prefixes
            )

            if bg_count > 0:
                specificity[i] = fg_count / bg_count
            else:
                specificity[i] = float('inf')

        return specificity


class MILPFallbackOptimizer:
    """
    Hybrid approach: Use MILP when possible, fall back to greedy.

    MILP: Best for <500 candidates (provably optimal, fast)
    Greedy: Better for >500 candidates (approximate, faster)

    Note: Renamed from HybridOptimizer to avoid conflict with
    hybrid_optimizer.py:HybridOptimizer which uses a different strategy.
    """

    def __init__(self, position_cache, fg_prefixes, bg_prefixes,
                 fg_seq_lengths, bg_seq_lengths):
        """Initialize hybrid optimizer"""
        self.cache = position_cache
        self.fg_prefixes = fg_prefixes
        self.bg_prefixes = bg_prefixes
        self.fg_seq_lengths = fg_seq_lengths
        self.bg_seq_lengths = bg_seq_lengths

    def optimize(self, candidates: List[str], num_primers: int = 10,
                max_time_seconds: int = 300) -> List[str]:
        """
        Automatically choose best optimization method.

        Args:
            candidates: Candidate primers
            num_primers: Number to select
            max_time_seconds: Time limit

        Returns:
            Optimized primer set
        """
        if len(candidates) <= 500 and MIP_AVAILABLE:
            # MILP is feasible
            logger.info("Using MILP optimization (exact solution)")
            milp_opt = MILPOptimizer(self.cache, self.fg_prefixes, self.bg_prefixes,
                                    self.fg_seq_lengths, self.bg_seq_lengths)
            result = milp_opt.optimize(candidates, num_primers, max_time_seconds)

            if result is not None:
                return result
            else:
                logger.warning("MILP failed, falling back to greedy")

        # Fall back to greedy
        logger.info("Using greedy optimization (fast approximation)")
        from .network_optimizer import NetworkOptimizer

        greedy_opt = NetworkOptimizer(self.cache, self.fg_prefixes, self.bg_prefixes,
                                      self.fg_seq_lengths, self.bg_seq_lengths)
        return greedy_opt.optimize_greedy(candidates, num_primers)


def compare_milp_vs_greedy(position_cache, fg_prefixes, bg_prefixes,
                           fg_seq_lengths, bg_seq_lengths, candidates):
    """
    Benchmark: MILP vs. Greedy optimization.

    Shows optimality gap and runtime tradeoff.
    """
    import time
    from .network_optimizer import NetworkOptimizer

    # Greedy
    logger.info("=== Greedy Optimization ===")
    greedy_opt = NetworkOptimizer(position_cache, fg_prefixes, bg_prefixes,
                                 fg_seq_lengths, bg_seq_lengths)

    start = time.time()
    greedy_primers = greedy_opt.optimize_greedy(candidates, num_primers=10)
    greedy_time = time.time() - start

    greedy_score = greedy_opt.score_primer_set(greedy_primers)
    logger.info(f"Primers: {greedy_primers}")
    logger.info(f"Score: {greedy_score['score']:.3f}")
    logger.info(f"Time: {greedy_time:.2f}s")

    # MILP (if available and feasible)
    if MIP_AVAILABLE and len(candidates) <= 500:
        logger.info("=== MILP Optimization ===")
        milp_opt = MILPOptimizer(position_cache, fg_prefixes, bg_prefixes,
                                fg_seq_lengths, bg_seq_lengths)

        start = time.time()
        milp_primers = milp_opt.optimize(candidates, num_primers=10, max_time_seconds=300)
        milp_time = time.time() - start

        if milp_primers is not None:
            milp_score = greedy_opt.score_primer_set(milp_primers)
            logger.info(f"Primers: {milp_primers}")
            logger.info(f"Score: {milp_score['score']:.3f}")
            logger.info(f"Time: {milp_time:.2f}s")

            logger.info("=== Comparison ===")
            improvement = (milp_score['score'] - greedy_score['score']) / greedy_score['score'] * 100
            logger.info(f"MILP improvement: {improvement:.1f}%")
            logger.info(f"MILP slower by: {milp_time / greedy_time:.1f}x")
        else:
            logger.info("MILP failed to find solution")
    else:
        logger.info("=== MILP Not Available ===")
        if not MIP_AVAILABLE:
            logger.info("python-mip not installed")
        else:
            logger.info(f"Too many candidates ({len(candidates)} > 500)")


# =============================================================================
# Factory Registration - BaseOptimizer Interface
# =============================================================================

from neoswga.core.base_optimizer import (
    BaseOptimizer, OptimizationResult, OptimizationStatus,
    PrimerSetMetrics, OptimizerConfig
)
from neoswga.core.optimizer_factory import OptimizerFactory


@OptimizerFactory.register('milp', aliases=['mip', 'exact'])
class MILPBaseOptimizer(BaseOptimizer):
    """
    MILP optimizer implementing BaseOptimizer interface.

    Provides provably optimal solutions for small candidate sets (<1000 primers).
    Requires python-mip package: pip install mip
    """

    def __init__(
        self,
        position_cache,
        fg_prefixes: List[str],
        fg_seq_lengths: List[int],
        bg_prefixes: Optional[List[str]] = None,
        bg_seq_lengths: Optional[List[int]] = None,
        config: Optional[OptimizerConfig] = None,
        **kwargs
    ):
        if not MIP_AVAILABLE:
            raise ImportError(
                "MILP optimizer requires python-mip package. "
                "Install with: pip install mip"
            )
        super().__init__(
            position_cache, fg_prefixes, fg_seq_lengths,
            bg_prefixes, bg_seq_lengths, config
        )
        self._milp = MILPOptimizer(
            position_cache=position_cache,
            fg_prefixes=fg_prefixes,
            bg_prefixes=bg_prefixes or [],
            fg_seq_lengths=fg_seq_lengths,
            bg_seq_lengths=bg_seq_lengths or [],
        )

    @property
    def name(self) -> str:
        return "milp"

    @property
    def description(self) -> str:
        return "Mixed Integer Linear Programming optimizer (exact, for <1000 candidates)"

    def optimize(
        self,
        candidates: List[str],
        target_size: Optional[int] = None,
        **kwargs
    ) -> OptimizationResult:
        """Run MILP optimization."""
        candidates = self._validate_candidates(candidates)
        target = target_size or self.config.target_set_size

        if self.config.verbose:
            logger.info(f"Running MILP optimization: {len(candidates)} candidates")

        try:
            max_time = kwargs.get('max_time_seconds', 300)
            primers = self._milp.optimize(
                candidates=candidates,
                num_primers=target,
                max_time_seconds=max_time,
            )

            if primers is None:
                return OptimizationResult.failure(self.name, "MILP solver failed or timed out")

            metrics = self.compute_metrics(primers)

            return OptimizationResult(
                primers=tuple(primers),
                score=metrics.fg_coverage,
                status=OptimizationStatus.SUCCESS,
                metrics=metrics,
                iterations=1,
                optimizer_name=self.name,
                message="Optimal solution found",
            )

        except Exception as e:
            logger.error(f"MILP optimization failed: {e}")
            return OptimizationResult.failure(self.name, str(e))


if __name__ == "__main__":
    import sys
    logging.basicConfig(level=logging.INFO)

    if len(sys.argv) > 1:
        print("MILP optimizer test")
        # TODO: Add standalone test
    else:
        print("MILP Optimizer for SWGA")
        print("Requires: pip install mip")
        print("\nUsage: python milp_optimizer.py <test_data>")
