#!/usr/bin/env python3
"""
Simulation-based fitness evaluation for SWGA primer sets.

Uses Phi29Simulator to predict actual amplification performance rather than
relying solely on network-based heuristics. Provides:
1. Actual coverage predictions
2. Amplification uniformity
3. Fork dynamics validation
4. Experimental time estimates

This is computationally expensive but provides ground truth for validation.

Author: NeoSWGA Development Team
Date: November 2025
Version: 3.0 - Phase 2.2
"""

import logging
import time
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np

from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.replication_simulator import Phi29Simulator, SimulationConfig, simulate_primer_set

logger = logging.getLogger(__name__)


@dataclass
class SimulationFitness:
    """Fitness metrics from simulation-based evaluation"""

    # Coverage metrics
    mean_coverage: float  # 0-1, fraction of genome covered
    std_coverage: float  # Variation across replicates

    # Amplification metrics
    mean_forks_created: float
    mean_fork_travel: float  # Average distance traveled (bp)

    # Uniformity metrics
    coverage_uniformity: float  # 0-1, how even is coverage

    # Composite score
    fitness_score: float  # Overall fitness (0-1)

    # Metadata
    simulation_time: float  # Runtime (seconds)
    n_replicates: int


class SimulationBasedEvaluator:
    """
    Evaluate primer sets using simulation.

    Provides ground-truth predictions of amplification performance.
    Much slower than network-based heuristics but more accurate.

    Use cases:
    - Final validation of selected primer sets
    - Comparing optimization methods
    - Parameter tuning and benchmarking
    """

    def __init__(
        self,
        genome_sequence: str,
        genome_length: int,
        position_cache,
        conditions: Optional[ReactionConditions] = None,
        n_replicates: int = 3,
        simulation_duration: float = 3600.0,
    ):
        """
        Initialize simulation-based evaluator.

        Args:
            genome_sequence: Full genome sequence
            genome_length: Genome length (bp)
            position_cache: PositionCache with primer positions
            conditions: Reaction conditions (uses defaults if None)
            n_replicates: Number of simulation replicates
            simulation_duration: Simulation time (seconds)
        """
        self.genome_sequence = genome_sequence
        self.genome_length = genome_length
        self.position_cache = position_cache
        self.conditions = conditions
        self.n_replicates = n_replicates
        self.simulation_duration = simulation_duration

        logger.info("Simulation-based evaluator initialized")
        logger.info(f"  Genome length: {genome_length:,} bp")
        logger.info(f"  Replicates: {n_replicates}")
        logger.info(f"  Duration: {simulation_duration:.0f}s")

    def evaluate(self, primers: List[str], verbose: bool = True) -> SimulationFitness:
        """
        Evaluate primer set using simulation.

        Args:
            primers: Primer sequences to evaluate
            verbose: Print progress

        Returns:
            SimulationFitness with comprehensive metrics
        """
        if verbose:
            logger.info(f"\nSimulating {len(primers)} primers...")

        start_time = time.time()

        # Build primer positions dict for simulator
        primer_positions = self._build_position_dict(primers)

        # Check if we have any positions
        total_positions = sum(
            len(positions.get("forward", [])) + len(positions.get("reverse", []))
            for positions in primer_positions.values()
        )

        if total_positions == 0:
            logger.warning("No primer positions found - returning zero fitness")
            return SimulationFitness(
                mean_coverage=0.0,
                std_coverage=0.0,
                mean_forks_created=0.0,
                mean_fork_travel=0.0,
                coverage_uniformity=0.0,
                fitness_score=0.0,
                simulation_time=time.time() - start_time,
                n_replicates=0,
            )

        # Run simulation with multiple replicates
        if verbose:
            logger.info(f"Running {self.n_replicates} replicates...")

        results = []

        # Resolve conditions ONCE and derive the polymerase from them, so the
        # simulated enzyme and the simulated buffer cannot disagree. They did:
        # with conditions=None the config named phi29 while the default
        # conditions were equiphi29 at 42 C, so phi29-length primers (Tm around
        # 26-31 C) could never prime and every set came back with zero forks and
        # zero coverage -- for a reason that had nothing to do with the primers.
        conditions = self.conditions or self._get_default_conditions()
        polymerase_type = self._get_polymerase_type(conditions)

        for i in range(self.n_replicates):
            if verbose:
                logger.info(f"  Replicate {i+1}/{self.n_replicates}")

            config = SimulationConfig(
                duration=self.simulation_duration,
                time_step=1.0,
                polymerase_type=polymerase_type,
            )

            simulator = Phi29Simulator(
                primers=primers,
                primer_positions=primer_positions,
                genome_length=self.genome_length,
                genome_sequence=self.genome_sequence,
                conditions=conditions,
                config=config,
            )

            # Run simulation
            result = simulator.run(verbose=False)
            results.append(result)

        # Aggregate results
        coverages = [r.final_coverage_fraction for r in results]
        mean_coverage = np.mean(coverages)
        std_coverage = np.std(coverages)

        mean_forks = np.mean([r.num_forks_created for r in results])
        mean_travel = np.mean([r.mean_fork_travel for r in results])

        # Calculate coverage uniformity from first replicate
        # High uniformity = even coverage across genome
        coverage_array = results[0].coverage
        if len(coverage_array) > 0 and np.max(coverage_array) > 0:
            # Coefficient of variation (lower = more uniform)
            cv = np.std(coverage_array) / (np.mean(coverage_array) + 1e-10)
            # Convert to 0-1 score (higher = more uniform)
            coverage_uniformity = 1.0 / (1.0 + cv)
        else:
            coverage_uniformity = 0.0

        # Calculate composite fitness score
        # Weights: coverage (60%), uniformity (30%), efficiency (10%)
        efficiency = min(mean_forks / (len(primers) * 100), 1.0)  # Normalized fork creation
        fitness_score = 0.6 * mean_coverage + 0.3 * coverage_uniformity + 0.1 * efficiency

        simulation_time = time.time() - start_time

        if verbose:
            logger.info(f"\nSimulation results:")
            logger.info(f"  Coverage: {mean_coverage:.1%} ± {std_coverage:.1%}")
            logger.info(f"  Uniformity: {coverage_uniformity:.2f}")
            logger.info(f"  Mean forks: {mean_forks:.0f}")
            logger.info(f"  Mean travel: {mean_travel:.0f} bp")
            logger.info(f"  Fitness score: {fitness_score:.3f}")
            logger.info(f"  Runtime: {simulation_time:.1f}s")

        return SimulationFitness(
            mean_coverage=mean_coverage,
            std_coverage=std_coverage,
            mean_forks_created=mean_forks,
            mean_fork_travel=mean_travel,
            coverage_uniformity=coverage_uniformity,
            fitness_score=fitness_score,
            simulation_time=simulation_time,
            n_replicates=self.n_replicates,
        )

    def _build_position_dict(self, primers: List[str]) -> Dict[str, Dict[str, List[int]]]:
        """Build position dictionary for simulator"""
        # PositionCache is keyed by (fname_prefix, primer, strand), where the
        # prefix is the real path stem the HDF5 was written under. This asked
        # for the empty string, which matches no key, so every primer came back
        # with no binding sites -- and `get_positions` returns an empty array
        # rather than raising, so nothing surfaced. Every fitness score this
        # evaluator produced was therefore computed from an empty binding map,
        # which for its stated purpose (comparing optimization methods) means
        # every set scored identically.
        prefixes = list(getattr(self.position_cache, "fname_prefixes", []) or [])
        if not prefixes:
            logger.warning(
                "Position cache exposes no genome prefixes; simulation fitness "
                "will see no binding sites and every set will score the same."
            )

        positions = {}
        total_sites = 0

        for primer in primers:
            forward_positions: List[int] = []
            reverse_positions: List[int] = []

            for prefix in prefixes:
                try:
                    fwd = self.position_cache.get_positions(prefix, primer, "forward")
                    rev = self.position_cache.get_positions(prefix, primer, "reverse")
                except (KeyError, AttributeError, TypeError) as e:
                    logger.debug(f"Cache access failed for primer {primer}: {e}")
                    continue

                forward_positions.extend(int(p) for p in fwd)
                reverse_positions.extend(int(p) for p in rev)

            total_sites += len(forward_positions) + len(reverse_positions)
            positions[primer] = {"forward": forward_positions, "reverse": reverse_positions}

        if primers and total_sites == 0:
            logger.warning(
                "No binding sites found for any of the %d primers. Simulated "
                "fitness will be zero for reasons unrelated to the primers -- "
                "check that the position cache was built for these prefixes.",
                len(primers),
            )

        return positions

    def _get_polymerase_type(self, conditions=None) -> str:
        """Determine the polymerase to simulate, from the conditions in use.

        Takes the conditions explicitly so the caller can pass the RESOLVED
        ones. Reading self.conditions here meant that when it was None the
        answer was "phi29" while the defaults actually applied were equiphi29 at
        42 C -- an enzyme and a buffer that do not go together.
        """
        if conditions is None:
            conditions = self.conditions
        if conditions is None:
            return "phi29"

        polymerase = getattr(conditions, "polymerase", None)
        if polymerase:
            return polymerase

        # Fall back to inferring from temperature.
        temp = getattr(conditions, "temp", None)
        if temp is not None and temp >= 42:
            return "equiphi29"
        return "phi29"

    def _get_default_conditions(self) -> ReactionConditions:
        """Default reaction conditions: standard phi29 at 30 C.

        This returned the ENHANCED preset (equiphi29, 42 C). That is not the
        project default -- `parameter.polymerase` is phi29 -- so a primer set
        designed for phi29 was being evaluated in a buffer 12 C hotter than the
        one it was designed for, where its primers do not bind.
        """
        from neoswga.core.reaction_conditions import get_standard_conditions

        return get_standard_conditions()

    def compare_sets(
        self, primer_sets: List[Tuple[str, List[str]]], verbose: bool = True
    ) -> List[Tuple[str, SimulationFitness]]:
        """
        Compare multiple primer sets using simulation.

        Args:
            primer_sets: List of (name, primers) tuples
            verbose: Print progress

        Returns:
            List of (name, fitness) sorted by fitness score (descending)
        """
        logger.info(f"\nComparing {len(primer_sets)} primer sets using simulation...")

        results = []

        for name, primers in primer_sets:
            logger.info(f"\n{'='*60}")
            logger.info(f"Evaluating: {name}")
            logger.info(f"{'='*60}")

            fitness = self.evaluate(primers, verbose=verbose)
            results.append((name, fitness))

        # Sort by fitness score (descending)
        results.sort(key=lambda x: x[1].fitness_score, reverse=True)

        # Print comparison
        logger.info(f"\n{'='*60}")
        logger.info("COMPARISON RESULTS")
        logger.info(f"{'='*60}")

        for i, (name, fitness) in enumerate(results, 1):
            logger.info(f"{i}. {name}")
            logger.info(f"   Score: {fitness.fitness_score:.3f}")
            logger.info(f"   Coverage: {fitness.mean_coverage:.1%} ± {fitness.std_coverage:.1%}")
            logger.info(f"   Uniformity: {fitness.coverage_uniformity:.2f}")

        return results


if __name__ == "__main__":
    # Example usage
    logging.basicConfig(level=logging.INFO)

    print("Simulation-Based Fitness Evaluator - Example Usage\n")
    print("This module integrates SWGA simulation for ground-truth fitness evaluation.")
    print("\nExample:")
    print("""
    from neoswga.core.simulation_fitness import SimulationBasedEvaluator

    evaluator = SimulationBasedEvaluator(
        genome_sequence=genome_seq,
        genome_length=4658411,
        position_cache=cache,
        n_replicates=5
    )

    # Evaluate single primer set
    fitness = evaluator.evaluate(primers)
    print(f"Fitness: {fitness.fitness_score:.3f}")
    print(f"Coverage: {fitness.mean_coverage:.1%}")

    # Compare multiple sets
    results = evaluator.compare_sets([
        ('hybrid', hybrid_primers),
        ('network', network_primers),
        ('dominating_set', domset_primers)
    ])
    """)
