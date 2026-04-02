"""
Genetic Algorithm optimization for SWGA primer set selection.

Implements evolutionary optimization to find optimal primer sets that:
- Maximize target genome coverage
- Minimize off-target amplification
- Avoid primer-dimer formation
- Maintain thermodynamic favorability

Superior to greedy search for escaping local optima and exploring
diverse primer combinations.
"""

import numpy as np
import random
import time
import warnings
from typing import List, Dict, Tuple, Optional, Callable
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

from neoswga.core import thermodynamics as thermo
from neoswga.core import reaction_conditions as rc
from neoswga.core import secondary_structure as ss


@dataclass
class GAConfig:
    """Configuration for genetic algorithm."""
    population_size: int = 200
    generations: int = 100
    mutation_rate: float = 0.15
    crossover_rate: float = 0.8
    elitism_fraction: float = 0.10
    tournament_size: int = 5
    min_set_size: int = 4
    max_set_size: int = 8
    max_dimer_severity: float = 0.5
    n_processes: int = None  # None = use all CPUs
    seed: Optional[int] = None  # Random seed for reproducibility


@dataclass
class Individual:
    """Represents a primer set (individual in population)."""
    primers: List[str]
    fitness: float = None
    metrics: Dict = None


class PrimerSetGA:
    """
    Genetic Algorithm for primer set optimization.

    Evolutionary operators:
    - Selection: Tournament selection
    - Crossover: Uniform crossover with dimer checking
    - Mutation: Add/remove/replace primers
    - Elitism: Preserve top performers
    """

    def __init__(self,
                 primer_pool: List[str],
                 fg_prefixes: List[str],
                 bg_prefixes: List[str],
                 fg_lengths: List[int],
                 bg_lengths: List[int],
                 conditions: rc.ReactionConditions,
                 config: Optional[GAConfig] = None,
                 position_cache=None):
        """
        Initialize genetic algorithm.

        Args:
            primer_pool: Available primers to choose from
            fg_prefixes: Foreground HDF5 file prefixes
            bg_prefixes: Background HDF5 file prefixes
            fg_lengths: Foreground genome lengths
            bg_lengths: Background genome lengths
            conditions: Reaction conditions
            config: GA configuration
            position_cache: Optional PositionCache for position lookups
        """
        self.fg_prefixes = fg_prefixes
        self.bg_prefixes = bg_prefixes if bg_prefixes else []
        self.fg_lengths = fg_lengths
        self.bg_lengths = bg_lengths if bg_lengths else []
        self.fg_total_length = sum(fg_lengths) if fg_lengths else 0
        self.bg_total_length = sum(self.bg_lengths) if self.bg_lengths else 0
        self.conditions = conditions
        self.config = config if config else GAConfig()
        self.position_cache = position_cache

        if self.config.n_processes is None:
            self.config.n_processes = multiprocessing.cpu_count()

        # Pre-filter large candidate pools to limit O(n^2) dimer matrix cost
        _MAX_DIMER_POOL = 200
        if len(primer_pool) > _MAX_DIMER_POOL and position_cache is not None:
            primer_pool = self._prefilter_by_selectivity(
                primer_pool, _MAX_DIMER_POOL
            )
        elif len(primer_pool) > _MAX_DIMER_POOL:
            primer_pool = primer_pool[:_MAX_DIMER_POOL]

        self.primer_pool = primer_pool

        # Precompute dimer matrix for primer pool
        print(f"Precomputing dimer matrix for {len(primer_pool)} primers...")
        self.dimer_matrix = ss.calculate_dimer_matrix(
            primer_pool, conditions
        )
        print(f"Dimer matrix: {self.dimer_matrix.shape}")

        # Map primers to indices
        self.primer_to_idx = {p: i for i, p in enumerate(primer_pool)}

        # Statistics tracking
        self.generation_stats = []
        self.best_individual_history = []

    def _prefilter_by_selectivity(self, primers: List[str], max_count: int) -> List[str]:
        """Pre-filter primers by selectivity ratio (fg/bg binding) to limit pool size."""
        scores = []
        for p in primers:
            fg_count = 0
            bg_count = 0
            for prefix in self.fg_prefixes:
                try:
                    fw = self.position_cache.get_positions(prefix, p, 'forward')
                    rv = self.position_cache.get_positions(prefix, p, 'reverse')
                    fg_count += len(fw) + len(rv)
                except (KeyError, Exception) as e:
                    logger.debug(f"Ignored error getting fg positions for {p}: {e}")
            for prefix in self.bg_prefixes:
                try:
                    fw = self.position_cache.get_positions(prefix, p, 'forward')
                    rv = self.position_cache.get_positions(prefix, p, 'reverse')
                    bg_count += len(fw) + len(rv)
                except (KeyError, Exception) as e:
                    logger.debug(f"Ignored error getting bg positions for {p}: {e}")
            selectivity = fg_count / (bg_count + 1) if fg_count > 0 else 0
            scores.append((selectivity, p))
        scores.sort(reverse=True)
        filtered = [p for _, p in scores[:max_count]]
        print(f"Pre-filtered {len(primers)} -> {len(filtered)} primers by selectivity")
        return filtered

    def evolve(self, verbose: bool = True, timeout: float = 0) -> Individual:
        """
        Run genetic algorithm evolution.

        Args:
            verbose: Print progress
            timeout: Maximum wall-clock seconds (0 = unlimited)

        Returns:
            Best individual found
        """
        start_time = time.time()

        # Set random seeds for reproducibility
        if self.config.seed is not None:
            random.seed(self.config.seed)
            np.random.seed(self.config.seed)
            if verbose:
                print(f"Random seed set to {self.config.seed} for reproducibility")

        # Initialize population
        population = self._initialize_population()

        if verbose:
            print(f"GA Configuration:")
            print(f"  Population: {self.config.population_size}")
            print(f"  Generations: {self.config.generations}")
            print(f"  Mutation rate: {self.config.mutation_rate:.2%}")
            print(f"  Crossover rate: {self.config.crossover_rate:.2%}")
            print(f"  Elitism: {self.config.elitism_fraction:.1%}")
            print(f"  Processes: {self.config.n_processes}")
            if timeout > 0:
                print(f"  Timeout: {timeout:.0f}s")
            print()

        # Evaluate initial population
        population = self._evaluate_population(population)

        best_overall = max(population, key=lambda x: x.fitness)
        timed_out = False

        # Evolution loop
        for generation in range(self.config.generations):
            # Check timeout
            if timeout > 0 and (time.time() - start_time) > timeout:
                if verbose:
                    print(f"Timeout ({timeout:.0f}s) reached at generation {generation}")
                timed_out = True
                break

            # Selection
            parents = self._select_parents(population)

            # Generate offspring
            offspring = []

            for i in range(0, len(parents) - 1, 2):
                parent1, parent2 = parents[i], parents[i + 1]

                # Crossover
                if random.random() < self.config.crossover_rate:
                    child1, child2 = self._crossover(parent1, parent2)
                else:
                    child1, child2 = parent1, parent2

                # Mutation
                if random.random() < self.config.mutation_rate:
                    child1 = self._mutate(child1)
                if random.random() < self.config.mutation_rate:
                    child2 = self._mutate(child2)

                offspring.extend([child1, child2])

            # Evaluate offspring
            offspring = self._evaluate_population(offspring)

            # Elitism: preserve best individuals
            n_elite = int(self.config.population_size * self.config.elitism_fraction)
            population.sort(key=lambda x: x.fitness, reverse=True)
            elite = population[:n_elite]

            # Combine elite + offspring
            population = elite + offspring[:(self.config.population_size - n_elite)]

            # Track best
            generation_best = max(population, key=lambda x: x.fitness)
            if generation_best.fitness > best_overall.fitness:
                best_overall = generation_best

            # Statistics
            fitnesses = [ind.fitness for ind in population]
            stats = {
                'generation': generation,
                'best_fitness': max(fitnesses),
                'mean_fitness': np.mean(fitnesses),
                'std_fitness': np.std(fitnesses),
                'median_fitness': np.median(fitnesses)
            }
            self.generation_stats.append(stats)
            self.best_individual_history.append(best_overall)

            if verbose and generation % 10 == 0:
                print(f"Generation {generation:3d}: "
                      f"Best={stats['best_fitness']:.4f}, "
                      f"Mean={stats['mean_fitness']:.4f}±{stats['std_fitness']:.4f}, "
                      f"Primers={len(generation_best.primers)}")

        elapsed = time.time() - start_time
        if verbose:
            status = "timed out" if timed_out else "complete"
            print(f"\nEvolution {status} in {elapsed:.1f}s")
            print(f"Best fitness: {best_overall.fitness:.4f}")
            print(f"Best primer set: {best_overall.primers}")
            print(f"Metrics: {best_overall.metrics}")

        best_overall.timed_out = timed_out
        return best_overall

    def _initialize_population(self) -> List[Individual]:
        """Create initial random population."""
        population = []

        for _ in range(self.config.population_size):
            # Random set size
            set_size = random.randint(
                self.config.min_set_size,
                self.config.max_set_size
            )

            # Random primers (ensuring no severe dimers)
            primer_set = self._random_compatible_set(set_size)

            population.append(Individual(primers=primer_set))

        return population

    def _random_compatible_set(self, size: int, max_attempts: int = 100) -> List[str]:
        """
        Generate random primer set with dimer constraints.

        Args:
            size: Desired set size
            max_attempts: Max attempts to find compatible set

        Returns:
            List of primer sequences
        """
        for attempt in range(max_attempts):
            primer_set = random.sample(self.primer_pool, size)

            # Check dimer compatibility
            if self._check_dimer_compatibility(primer_set):
                return primer_set

        # If failed, return best-effort set
        warnings.warn(f"Could not find fully compatible set of size {size}")
        return random.sample(self.primer_pool, min(size, self.config.min_set_size))

    def _check_dimer_compatibility(self, primer_set: List[str]) -> bool:
        """
        Check if primer set has acceptable dimer interactions.

        Args:
            primer_set: List of primers

        Returns:
            True if no severe dimers
        """
        indices = [self.primer_to_idx[p] for p in primer_set]

        for i in range(len(indices)):
            for j in range(i, len(indices)):
                severity = self.dimer_matrix[indices[i], indices[j]]
                if severity > self.config.max_dimer_severity:
                    return False

        return True

    def _evaluate_population(self, population: List[Individual]) -> List[Individual]:
        """
        Evaluate fitness for all individuals in population.

        Uses multiprocessing for parallel evaluation.

        Args:
            population: List of individuals

        Returns:
            Population with fitness scores
        """
        # Parallel evaluation
        with ProcessPoolExecutor(max_workers=self.config.n_processes) as executor:
            futures = {
                executor.submit(self._evaluate_individual, ind): ind
                for ind in population if ind.fitness is None
            }

            for future in as_completed(futures):
                ind = futures[future]
                fitness, metrics = future.result()
                ind.fitness = fitness
                ind.metrics = metrics

        return population

    def _evaluate_individual(self, individual: Individual) -> Tuple[float, Dict]:
        """
        Evaluate fitness of single individual.

        Fitness components:
        1. Coverage (binding site distribution)
        2. Specificity (fg/bg ratio)
        3. Evenness (Gini index)
        4. Thermodynamic favorability
        5. Lack of dimers

        Args:
            individual: Individual to evaluate

        Returns:
            (fitness, metrics)
        """
        if not individual.primers:
            return 0.0, {}

        # Calculate position-based metrics using position cache
        fg_total = 0
        bg_total = 0
        all_fg_positions = []

        if self.position_cache is not None:
            # Use actual position data from cache
            for primer in individual.primers:
                for prefix in self.fg_prefixes:
                    try:
                        fw = self.position_cache.get_positions(prefix, primer, 'forward')
                        rv = self.position_cache.get_positions(prefix, primer, 'reverse')
                        fg_count = len(fw) + len(rv)
                        fg_total += fg_count
                        all_fg_positions.extend(fw.tolist() if hasattr(fw, 'tolist') else list(fw))
                        all_fg_positions.extend(rv.tolist() if hasattr(rv, 'tolist') else list(rv))
                    except KeyError:
                        # Primer not found in cache - expected for some primers
                        pass
                    except Exception as e:
                        logger.debug(f"Error getting fg positions for {primer}: {e}")

                for prefix in self.bg_prefixes:
                    try:
                        fw = self.position_cache.get_positions(prefix, primer, 'forward')
                        rv = self.position_cache.get_positions(prefix, primer, 'reverse')
                        bg_total += len(fw) + len(rv)
                    except KeyError:
                        # Primer not found in cache - expected for some primers
                        pass
                    except Exception as e:
                        logger.debug(f"Error getting bg positions for {primer}: {e}")

            # Component 1: Coverage score (0-1)
            # Normalized binding frequency in target genome
            coverage_score = min(1.0, fg_total / (self.fg_total_length * 1e-5))

            # Component 2: Specificity score (0-1)
            # Ratio-based: higher fg/bg is better
            fg_freq = fg_total / self.fg_total_length if self.fg_total_length > 0 else 0
            bg_freq = bg_total / self.bg_total_length if self.bg_total_length > 0 else 1e-10
            ratio = fg_freq / (bg_freq + 1e-10)
            specificity_score = min(1.0, ratio / 100.0)  # Normalize: 100× ratio = 1.0

            # Component 3: Evenness score (0-1)
            # Lower Gini = more even distribution = higher score
            if len(all_fg_positions) > 1:
                all_fg_positions = sorted(all_fg_positions)
                gaps = np.diff(all_fg_positions)
                if len(gaps) > 0:
                    gini = self._calculate_gini(gaps)
                    evenness_score = 1.0 - gini  # Invert: low gini = high evenness
                else:
                    evenness_score = 0.5
            else:
                evenness_score = 0.0

        else:
            # Fallback to random if no cache (backward compatibility)
            coverage_score = random.random() * 0.5 + 0.25  # 0.25-0.75 range
            specificity_score = random.random() * 0.5 + 0.25
            evenness_score = random.random() * 0.5 + 0.25

        # Component 4: Thermodynamic score (0-1)
        # Higher if primers have favorable Tm and delta-G
        tm_scores = []
        dg_scores = []
        for primer in individual.primers:
            tm = self.conditions.calculate_effective_tm(primer)
            tm_midpoint = np.mean(self.conditions.get_polymerase_range())
            tm_score = np.exp(-((tm - tm_midpoint) / 5)**2)
            tm_scores.append(tm_score)

            dg = thermo.calculate_free_energy(primer, self.conditions.temp)
            dg_score = 1 / (1 + np.exp((dg + 10) / 3))
            dg_scores.append(dg_score)

        thermo_score = (np.mean(tm_scores) + np.mean(dg_scores)) / 2

        # Component 5: Dimer avoidance score (0-1)
        # Higher if no severe dimers
        dimer_score = 1 - self._calculate_max_dimer_severity(individual.primers)

        # Weighted combination
        fitness = (
            0.35 * coverage_score +
            0.30 * specificity_score +
            0.15 * evenness_score +
            0.10 * thermo_score +
            0.10 * dimer_score
        )

        metrics = {
            'coverage': coverage_score,
            'specificity': specificity_score,
            'evenness': evenness_score,
            'thermodynamic': thermo_score,
            'dimer': dimer_score,
            'fg_sites': fg_total,
            'bg_sites': bg_total
        }

        return fitness, metrics

    def _calculate_gini(self, values: np.ndarray) -> float:
        """Calculate Gini coefficient for evenness measurement."""
        if len(values) == 0:
            return 0.0
        sorted_vals = np.sort(values)
        n = len(sorted_vals)
        total_sum = np.sum(sorted_vals)
        # Guard against division by zero when all values are 0
        if total_sum == 0:
            return 0.0
        return (2 * np.sum((np.arange(1, n + 1) * sorted_vals)) / (n * total_sum)) - (n + 1) / n

    def _calculate_max_dimer_severity(self, primer_set: List[str]) -> float:
        """Calculate maximum dimer severity in set."""
        if len(primer_set) <= 1:
            return 0.0

        indices = [self.primer_to_idx[p] for p in primer_set]
        max_severity = 0.0

        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                severity = self.dimer_matrix[indices[i], indices[j]]
                max_severity = max(max_severity, severity)

        return max_severity

    def _select_parents(self, population: List[Individual]) -> List[Individual]:
        """
        Select parents using tournament selection.

        Args:
            population: Current population

        Returns:
            Selected parents
        """
        parents = []

        for _ in range(self.config.population_size):
            # Tournament
            tournament = random.sample(population, self.config.tournament_size)
            winner = max(tournament, key=lambda x: x.fitness)
            parents.append(winner)

        return parents

    def _crossover(self, parent1: Individual, parent2: Individual) -> Tuple[Individual, Individual]:
        """
        Uniform crossover: randomly mix primers from both parents.

        Args:
            parent1: First parent
            parent2: Second parent

        Returns:
            Two offspring
        """
        # Combine primer pools from both parents
        all_primers = list(set(parent1.primers + parent2.primers))

        # Random sizes for offspring
        size1 = random.randint(self.config.min_set_size, self.config.max_set_size)
        size2 = random.randint(self.config.min_set_size, self.config.max_set_size)

        # Select from combined parent primers, with dimer compatibility check
        # Falls back to full pool only if parent primers are insufficient
        child1_primers = self._select_from_pool(all_primers, min(size1, len(all_primers)))
        child2_primers = self._select_from_pool(all_primers, min(size2, len(all_primers)))

        return Individual(primers=child1_primers), Individual(primers=child2_primers)

    def _select_from_pool(self, pool: List[str], size: int, max_attempts: int = 50) -> List[str]:
        """Select a dimer-compatible primer set from a given pool.

        Tries to build a set from the provided pool first. Only falls back
        to the full primer_pool if the given pool is too small.
        """
        source = pool if len(pool) >= size else self.primer_pool
        for _ in range(max_attempts):
            primer_set = random.sample(source, min(size, len(source)))
            if self._check_dimer_compatibility(primer_set):
                return primer_set
        # Fallback: return best-effort sample from the source pool
        return random.sample(source, min(size, len(source)))

    def _mutate(self, individual: Individual) -> Individual:
        """
        Mutate individual by adding/removing/replacing primers.

        Mutation types (equal probability):
        1. Add random primer
        2. Remove random primer
        3. Replace random primer

        Args:
            individual: Individual to mutate

        Returns:
            Mutated individual
        """
        mutation_type = random.choice(['add', 'remove', 'replace'])

        primers = individual.primers.copy()

        if mutation_type == 'add' and len(primers) < self.config.max_set_size:
            # Add random primer from pool (avoiding dimers)
            available = [p for p in self.primer_pool if p not in primers]
            if available:
                for _ in range(10):  # Try 10 times to find compatible primer
                    candidate = random.choice(available)
                    test_set = primers + [candidate]
                    if self._check_dimer_compatibility(test_set):
                        primers.append(candidate)
                        break

        elif mutation_type == 'remove' and len(primers) > self.config.min_set_size:
            # Remove random primer
            primers.remove(random.choice(primers))

        elif mutation_type == 'replace' and primers:
            # Replace random primer
            old_primer = random.choice(primers)
            primers.remove(old_primer)

            available = [p for p in self.primer_pool if p not in primers]
            if available:
                for _ in range(10):
                    candidate = random.choice(available)
                    test_set = primers + [candidate]
                    if self._check_dimer_compatibility(test_set):
                        primers.append(candidate)
                        break

        return Individual(primers=primers)


def optimize_primer_set_ga(primer_pool: List[str],
                           fg_prefixes: List[str],
                           bg_prefixes: List[str],
                           fg_lengths: List[int],
                           bg_lengths: List[int],
                           conditions: Optional[rc.ReactionConditions] = None,
                           config: Optional[GAConfig] = None) -> Individual:
    """
    High-level function for GA optimization.

    Args:
        primer_pool: Available primers
        fg_prefixes: Foreground file prefixes
        bg_prefixes: Background file prefixes
        fg_lengths: Foreground genome lengths
        bg_lengths: Background genome lengths
        conditions: Reaction conditions
        config: GA configuration

    Returns:
        Best primer set found
    """
    if conditions is None:
        conditions = rc.get_enhanced_conditions()

    ga = PrimerSetGA(
        primer_pool, fg_prefixes, bg_prefixes,
        fg_lengths, bg_lengths,
        conditions, config
    )

    return ga.evolve(verbose=True)


# =============================================================================
# Factory Registration - BaseOptimizer Interface
# =============================================================================

import logging
from neoswga.core.base_optimizer import (
    BaseOptimizer, OptimizationResult, OptimizationStatus,
    PrimerSetMetrics, OptimizerConfig
)
from neoswga.core.optimizer_factory import OptimizerFactory

logger = logging.getLogger(__name__)


@OptimizerFactory.register('genetic', aliases=['ga', 'genetic-algorithm'])
class GeneticBaseOptimizer(BaseOptimizer):
    """
    Genetic algorithm optimizer implementing BaseOptimizer interface.

    Uses evolutionary optimization for multi-objective primer selection.
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
            bg_prefixes, bg_seq_lengths, config
        )
        self.conditions = conditions
        max_gen = min(50, config.max_iterations) if config else 50
        self.ga_config = GAConfig(
            population_size=kwargs.get('population_size', 50),
            generations=max_gen,
            min_set_size=max(4, (config.target_set_size if config else 6) - 2),
            max_set_size=(config.target_set_size if config else 6) + 2,
            seed=kwargs.get('seed', 42),
        )

    @property
    def name(self) -> str:
        """Optimizer identifier for logging and factory registration."""
        return "genetic"

    @property
    def description(self) -> str:
        """One-line summary of the optimization strategy."""
        return "Genetic algorithm for multi-objective primer optimization"

    def optimize(
        self,
        candidates: List[str],
        target_size: Optional[int] = None,
        **kwargs
    ) -> OptimizationResult:
        """Run genetic algorithm optimization."""
        candidates = self._validate_candidates(candidates)
        target = target_size or self.config.target_set_size

        if self.config.verbose:
            logger.info(f"Running genetic algorithm: {len(candidates)} candidates")

        try:
            # Get reaction conditions
            conditions = self.conditions
            if conditions is None:
                conditions = rc.get_standard_conditions()

            # Create GA optimizer
            ga = PrimerSetGA(
                primer_pool=candidates,
                fg_prefixes=self.fg_prefixes,
                bg_prefixes=self.bg_prefixes,
                fg_lengths=self.fg_seq_lengths,
                bg_lengths=self.bg_seq_lengths,
                conditions=conditions,
                config=self.ga_config,
                position_cache=self.cache,
            )

            # Run evolution with 60s timeout
            _GA_TIMEOUT = 60
            best = ga.evolve(verbose=self.config.verbose, timeout=_GA_TIMEOUT)

            primers = best.primers
            metrics = self.compute_metrics(primers)

            timed_out = getattr(best, 'timed_out', False)
            if not primers:
                status = OptimizationStatus.NO_CONVERGENCE
            elif timed_out:
                status = OptimizationStatus.PARTIAL
            else:
                status = OptimizationStatus.SUCCESS

            return OptimizationResult(
                primers=tuple(primers),
                score=best.fitness,
                status=status,
                metrics=metrics,
                iterations=len(ga.generation_stats),
                optimizer_name=self.name,
                message="Timeout: returning best result found" if timed_out else "",
            )

        except Exception as e:
            logger.error(f"Genetic algorithm failed: {e}")
            return OptimizationResult.failure(self.name, str(e))
