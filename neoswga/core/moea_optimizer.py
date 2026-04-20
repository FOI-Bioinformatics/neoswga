"""
Multi-Objective Evolutionary Algorithm (MOEA) for primer set optimization.

Uses NSGA-III to find Pareto-optimal solutions across multiple objectives:
1. Maximize target genome coverage
2. Minimize number of primers (cost)
3. Minimize off-target binding
4. Maximize binding uniformity

Returns Pareto front of solutions instead of single optimum.
This allows users to choose based on their priorities.
"""

import numpy as np
import logging
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
import time

logger = logging.getLogger(__name__)

try:
    from pymoo.algorithms.moo.nsga3 import NSGA3
    from pymoo.core.problem import Problem
    from pymoo.optimize import minimize
    from pymoo.util.ref_dirs import get_reference_directions
    HAS_PYMOO = True
except ImportError:
    HAS_PYMOO = False
    logger.warning("pymoo not installed. MOEA optimizer unavailable. Install: pip install pymoo")


@dataclass
class MOEAConfig:
    """Configuration for multi-objective optimization"""
    pop_size: int = 100  # Population size
    n_generations: int = 100  # Number of generations
    n_objectives: int = 4  # Number of objectives
    crossover_prob: float = 0.9
    mutation_prob: float = 0.1
    seed: Optional[int] = None


# Define stub base class if pymoo is not available
if not HAS_PYMOO:
    class Problem:
        """Stub Problem class when pymoo is not installed."""
        def __init__(self, *args, **kwargs):
            pass


class PrimerSetProblem(Problem):
    """
    Multi-objective optimization problem for primer selection.

    Decision variables: Binary vector (0/1 for each candidate primer)
    Objectives:
        1. Maximize target coverage (minimize -coverage)
        2. Minimize number of primers
        3. Minimize off-target binding
        4. Maximize uniformity (minimize -uniformity)
    """

    def __init__(self,
                 candidates: List[str],
                 cache,
                 fg_prefixes: List[str],
                 bg_prefixes: List[str],
                 fg_seq_lengths: List[int],
                 bg_seq_lengths: List[int],
                 max_primers: int = 20):
        """
        Initialize problem.

        Args:
            candidates: List of candidate primers
            cache: PositionCache with primer positions
            fg_prefixes: Target genome prefixes
            bg_prefixes: Background genome prefixes
            fg_seq_lengths: Target genome lengths
            bg_seq_lengths: Background genome lengths
            max_primers: Maximum number of primers to select
        """
        self.candidates = candidates
        self.cache = cache
        self.fg_prefixes = fg_prefixes
        self.bg_prefixes = bg_prefixes
        self.fg_seq_lengths = fg_seq_lengths
        self.bg_seq_lengths = bg_seq_lengths
        self.max_primers = max_primers

        # Problem definition
        n_var = len(candidates)  # Binary decision for each candidate
        n_obj = 4  # Four objectives
        n_constr = 1  # Constraint: select between 5 and max_primers

        super().__init__(
            n_var=n_var,
            n_obj=n_obj,
            n_constr=n_constr,
            xl=0,  # Lower bound (binary)
            xu=1,  # Upper bound (binary)
            type_var=int,  # Integer (binary) variables
        )

    def _evaluate(self, X, out, *args, **kwargs):
        """
        Evaluate objectives for population X.

        Args:
            X: Population matrix (n_individuals × n_variables)
            out: Output dictionary
        """
        n_individuals = X.shape[0]

        # Objectives matrix
        F = np.zeros((n_individuals, 4))

        # Constraint violations
        G = np.zeros((n_individuals, 1))

        for i in range(n_individuals):
            # Individual is binary vector
            individual = X[i]

            # Selected primers
            selected_primers = [self.candidates[j] for j in range(len(self.candidates))
                              if individual[j] == 1]

            n_selected = len(selected_primers)

            # Constraint: select between 5 and max_primers
            if n_selected < 5:
                G[i, 0] = 5 - n_selected
            elif n_selected > self.max_primers:
                G[i, 0] = n_selected - self.max_primers
            else:
                G[i, 0] = 0

            if n_selected == 0:
                # Invalid solution
                F[i, :] = [1e6, 1e6, 1e6, 1e6]
                continue

            # Objective 1: Maximize target coverage (minimize -coverage)
            target_coverage = self._calculate_coverage(selected_primers,
                                                      self.fg_prefixes,
                                                      self.fg_seq_lengths)
            F[i, 0] = -target_coverage  # Minimize negative (maximize positive)

            # Objective 2: Minimize number of primers
            F[i, 1] = n_selected

            # Objective 3: Minimize off-target binding
            bg_binding = self._calculate_binding(selected_primers,
                                                self.bg_prefixes)
            F[i, 2] = bg_binding

            # Objective 4: Maximize uniformity (minimize -uniformity)
            uniformity = self._calculate_uniformity(selected_primers,
                                                   self.fg_prefixes,
                                                   self.fg_seq_lengths)
            F[i, 3] = -uniformity  # Minimize negative (maximize positive)

        out["F"] = F
        out["G"] = G

    def _calculate_coverage(self, primers: List[str], prefixes: List[str],
                           seq_lengths: List[int]) -> float:
        """
        Calculate genome coverage.

        Returns fraction of genome covered by primers.
        """
        if not primers:
            return 0.0

        total_length = sum(seq_lengths)
        covered_bases = 0

        for prefix, length in zip(prefixes, seq_lengths):
            # Get all binding positions
            positions = []
            for primer in primers:
                pos_fwd = self.cache.get_positions(prefix, primer, 'forward')
                pos_rev = self.cache.get_positions(prefix, primer, 'reverse')
                positions.extend(pos_fwd)
                positions.extend(pos_rev)

            if positions:
                # Count unique positions (simplified coverage)
                covered_bases += len(set(positions))

        return covered_bases / total_length if total_length > 0 else 0.0

    def _calculate_binding(self, primers: List[str], prefixes: List[str]) -> float:
        """
        Calculate total background binding.

        Returns total number of binding sites in background.
        """
        if not primers:
            return 0.0

        total_binding = 0

        for prefix in prefixes:
            for primer in primers:
                pos_fwd = self.cache.get_positions(prefix, primer, 'forward')
                pos_rev = self.cache.get_positions(prefix, primer, 'reverse')
                total_binding += len(pos_fwd) + len(pos_rev)

        return total_binding

    def _calculate_uniformity(self, primers: List[str], prefixes: List[str],
                             seq_lengths: List[int]) -> float:
        """
        Calculate binding uniformity across genome.

        Returns 1 - CV (coefficient of variation) of binding density.
        Higher is better (more uniform).
        """
        if not primers:
            return 0.0

        # Divide genome into bins
        n_bins = 100

        densities = []
        for prefix, length in zip(prefixes, seq_lengths):
            bin_size = length // n_bins
            bin_counts = np.zeros(n_bins)

            # Count bindings per bin
            for primer in primers:
                pos_fwd = self.cache.get_positions(prefix, primer, 'forward')
                pos_rev = self.cache.get_positions(prefix, primer, 'reverse')
                all_pos = np.concatenate([pos_fwd, pos_rev])

                for pos in all_pos:
                    bin_idx = min(int(pos / bin_size), n_bins - 1)
                    bin_counts[bin_idx] += 1

            densities.extend(bin_counts)

        densities = np.array(densities)

        if densities.sum() == 0:
            return 0.0

        # Coefficient of variation with epsilon guard for numerical stability
        mean = densities.mean()
        std = densities.std()

        # Use epsilon to prevent division by very small numbers
        EPSILON = 1e-10
        if mean < EPSILON:
            return 0.0

        cv = std / max(mean, EPSILON)

        # Return 1 - CV (higher is more uniform)
        return 1.0 / (1.0 + cv)


class MOEAOptimizer:
    """
    Multi-objective evolutionary algorithm optimizer.

    Uses NSGA-III to find Pareto-optimal primer sets.
    Returns multiple solutions representing different tradeoffs.
    """

    def __init__(self, cache, fg_prefixes: List[str], bg_prefixes: List[str],
                 fg_seq_lengths: List[int], bg_seq_lengths: List[int],
                 config: Optional[MOEAConfig] = None):
        """
        Initialize MOEA optimizer.

        Args:
            cache: PositionCache with primer positions
            fg_prefixes: Target genome prefixes
            bg_prefixes: Background genome prefixes
            fg_seq_lengths: Target genome lengths
            bg_seq_lengths: Background genome lengths
            config: MOEA configuration
        """
        if not HAS_PYMOO:
            raise ImportError(
                "pymoo required for MOEA optimizer. "
                "Install: pip install pymoo"
            )

        self.cache = cache
        self.fg_prefixes = fg_prefixes
        self.bg_prefixes = bg_prefixes
        self.fg_seq_lengths = fg_seq_lengths
        self.bg_seq_lengths = bg_seq_lengths
        self.config = config or MOEAConfig()

    def optimize(self, candidates: List[str], max_primers: int = 20,
                verbose: bool = True) -> Dict:
        """
        Run MOEA optimization.

        Args:
            candidates: List of candidate primers
            max_primers: Maximum primers to select
            verbose: Print progress

        Returns:
            Dictionary with Pareto front solutions
        """
        if verbose:
            logger.info(f"Running MOEA optimization: {len(candidates)} candidates")
            logger.info(f"Population: {self.config.pop_size}, "
                       f"Generations: {self.config.n_generations}")

        # Create problem
        problem = PrimerSetProblem(
            candidates=candidates,
            cache=self.cache,
            fg_prefixes=self.fg_prefixes,
            bg_prefixes=self.bg_prefixes,
            fg_seq_lengths=self.fg_seq_lengths,
            bg_seq_lengths=self.bg_seq_lengths,
            max_primers=max_primers
        )

        # Reference directions for NSGA-III
        ref_dirs = get_reference_directions(
            "energy",
            self.config.n_objectives,
            self.config.pop_size
        )

        # Create algorithm
        algorithm = NSGA3(
            ref_dirs=ref_dirs,
            pop_size=self.config.pop_size,
            seed=self.config.seed
        )

        # Run optimization
        start_time = time.time()

        res = minimize(
            problem,
            algorithm,
            ('n_gen', self.config.n_generations),
            seed=self.config.seed,
            verbose=verbose
        )

        runtime = time.time() - start_time

        if verbose:
            logger.info(f"MOEA complete: {runtime:.1f}s")

        # Handle case where no feasible solutions were found
        if res.X is None:
            logger.warning("MOEA found no feasible solutions")
            return {
                'pareto_front': [],
                'best_solution': {'primers': [], 'score': 0.0},
                'runtime': runtime,
                'n_solutions': 0
            }

        if verbose:
            logger.info(f"Found {len(res.X)} Pareto-optimal solutions")

        # Extract solutions
        pareto_solutions = []

        for i in range(len(res.X)):
            individual = res.X[i]
            objectives = res.F[i]

            # Selected primers
            selected = [candidates[j] for j in range(len(candidates))
                       if individual[j] == 1]

            solution = {
                'primers': selected,
                'objectives': {
                    'target_coverage': -objectives[0],  # Negated back
                    'n_primers': int(objectives[1]),
                    'background_binding': objectives[2],
                    'uniformity': -objectives[3]  # Negated back
                },
                'score': self._compute_score(objectives)
            }

            pareto_solutions.append(solution)

        # Sort by score (best first)
        pareto_solutions.sort(key=lambda x: x['score'], reverse=True)

        return {
            'pareto_front': pareto_solutions,
            'best_solution': pareto_solutions[0],
            'runtime': runtime,
            'n_solutions': len(pareto_solutions)
        }

    def _compute_score(self, objectives: np.ndarray) -> float:
        """
        Compute scalar score from multi-objective values.

        Weighted combination (can be customized).
        """
        # Weights (can be tuned)
        w_coverage = 2.0
        w_primers = -0.5  # Fewer is better
        w_background = -1.0  # Less is better
        w_uniformity = 1.0

        score = (w_coverage * (-objectives[0]) +  # Coverage (negated)
                w_primers * objectives[1] +
                w_background * objectives[2] / 1000.0 +  # Scale down
                w_uniformity * (-objectives[3]))  # Uniformity (negated)

        return score


# =============================================================================
# Factory Registration - BaseOptimizer Interface
# =============================================================================

if HAS_PYMOO:
    from neoswga.core.base_optimizer import (
        BaseOptimizer, OptimizationResult, OptimizationStatus,
        PrimerSetMetrics, OptimizerConfig
    )
    from neoswga.core.optimizer_factory import OptimizerFactory

    @OptimizerFactory.register('moea', aliases=['nsga3', 'multi-objective', 'pareto'])
    class MOEABaseOptimizer(BaseOptimizer):
        """
        MOEA optimizer implementing BaseOptimizer interface.

        Multi-objective optimization using NSGA-III algorithm.
        Returns the best solution from the Pareto front.
        Requires pymoo package.
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
            super().__init__(
                position_cache, fg_prefixes, fg_seq_lengths,
                bg_prefixes, bg_seq_lengths, config,
            conditions=conditions,
        )
            moea_config = MOEAConfig(
                pop_size=kwargs.get('pop_size', 100),
                n_generations=kwargs.get('n_generations', 100),
                # Forward --seed so same-seed runs produce identical Pareto
                # fronts. Without this MOEA was non-deterministic even when
                # the user passed --seed; see tests/test_reproducibility_moea.
                seed=kwargs.get('seed'),
            )
            self._moea = MOEAOptimizer(
                cache=position_cache,
                fg_prefixes=fg_prefixes,
                bg_prefixes=bg_prefixes or [],
                fg_seq_lengths=fg_seq_lengths,
                bg_seq_lengths=bg_seq_lengths or [],
                config=moea_config,
            )

        @property
        def name(self) -> str:
            return "moea"

        @property
        def description(self) -> str:
            return "Multi-Objective Evolutionary Algorithm (NSGA-III, Pareto optimization)"

        def optimize(
            self,
            candidates: List[str],
            target_size: Optional[int] = None,
            **kwargs
        ) -> OptimizationResult:
            """Run MOEA optimization."""
            candidates = self._validate_candidates(candidates)
            target = target_size or self.config.target_set_size

            if self.config.verbose:
                logger.info(f"Running MOEA optimization: {len(candidates)} candidates")

            try:
                result = self._moea.optimize(
                    candidates=candidates,
                    max_primers=target,
                    verbose=self.config.verbose,
                )

                if not result['pareto_front']:
                    return OptimizationResult.failure(self.name, "No Pareto-optimal solutions found")

                # Get best solution
                best = result['best_solution']
                primers = best['primers']
                metrics = self.compute_metrics(primers)

                return OptimizationResult(
                    primers=tuple(primers),
                    score=best['score'],
                    status=OptimizationStatus.SUCCESS,
                    metrics=metrics,
                    iterations=result.get('runtime', 1),
                    optimizer_name=self.name,
                    message=f"Found {result['n_solutions']} Pareto-optimal solutions",
                )

            except Exception as e:
                logger.error(f"MOEA optimization failed: {e}")
                return OptimizationResult.failure(self.name, str(e))


if __name__ == "__main__":
    print("MOEA Optimizer for Primer Set Selection")
    print("\nRequires:")
    print("  - pymoo library: pip install pymoo")
    print("  - PositionCache with primer data")
    print("\nFeatures:")
    print("  - Multi-objective optimization (NSGA-III)")
    print("  - Returns Pareto front of solutions")
    print("  - Objectives: coverage, primers, off-target, uniformity")
    print("\nUsage:")
    print("  from neoswga.core.moea_optimizer import MOEAOptimizer")
    print("  optimizer = MOEAOptimizer(cache, fg_prefixes, bg_prefixes,")
    print("                           fg_seq_lengths, bg_seq_lengths)")
    print("  result = optimizer.optimize(candidates, max_primers=15)")
    print("  print(f'Found {result[\"n_solutions\"]} Pareto-optimal solutions')")
