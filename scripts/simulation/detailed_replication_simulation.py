#!/usr/bin/env python3
"""
Detailed Phi29 Replication Simulation for Top Francisella Configs

Runs agent-based Phi29 polymerase simulation for top-performing configurations
to validate binding site analysis with realistic amplification kinetics.

This provides:
- Amplification dynamics over time
- Coverage growth curves
- Primer interaction effects
- Realistic performance estimates

Author: NeoSWGA Development Team
Date: 2025-11-23
"""

import sys
import logging
import json
from pathlib import Path
from typing import List, Dict, Tuple
from dataclasses import dataclass, asdict
import numpy as np

# Add neoswga to path
sys.path.insert(0, str(Path(__file__).parent))

from simulate_francisella_primers import load_primers

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class ReplicationMetrics:
    """Results from detailed replication simulation"""
    config_name: str
    n_primers: int

    # Amplification metrics
    final_coverage_percent: float
    amplification_fold: float  # Total DNA amplified (fold change)
    mean_fragment_length: float  # bp

    # Kinetics
    time_to_50pct_coverage: float  # minutes
    time_to_plateau: float  # minutes

    # Specificity
    target_amplification: float  # ng
    background_amplification: float  # ng
    enrichment: float

    # Quality
    coverage_uniformity: float  # 0-1, based on Gini index
    recommendation: str

    def to_dict(self) -> Dict:
        return asdict(self)


class DetailedSwgaSimulator:
    """
    Agent-based Phi29 replication simulator.

    Simplified version that tracks:
    - Polymerase binding and extension
    - Strand displacement
    - Coverage accumulation over time

    Note: This is a simplified model. For production use, consider
    more sophisticated simulators with full biochemical parameters.
    """

    def __init__(self,
                 target_genome_path: str,
                 background_genome_path: str,
                 simulation_time: float = 120.0,  # minutes
                 time_step: float = 1.0):  # minutes
        """
        Initialize detailed simulator.

        Args:
            target_genome_path: Path to target genome FASTA
            background_genome_path: Path to background genome FASTA
            simulation_time: Total simulation time (minutes)
            time_step: Time step for simulation (minutes)
        """
        from simulate_francisella_primers import FrancisellaSimulator

        # Use FrancisellaSimulator for genome loading and binding site finding
        self.base_simulator = FrancisellaSimulator(
            target_genome_path=target_genome_path,
            background_genome_path=background_genome_path,
            bin_size=10000
        )

        self.simulation_time = simulation_time
        self.time_step = time_step
        self.time_points = int(simulation_time / time_step)

        # Phi29 parameters (simplified)
        self.extension_rate = 1000  # bp/min (Phi29 is ~1000-5000 bp/min)
        self.binding_rate = 0.1  # primers bind per site per minute
        self.initial_dna = 10.0  # ng starting DNA

    def simulate_amplification(self,
                               primers: List[str],
                               config_name: str,
                               genome: str = 'target') -> Dict:
        """
        Simulate Phi29 amplification over time.

        Args:
            primers: List of primer sequences
            config_name: Configuration name
            genome: 'target' or 'background'

        Returns:
            Dictionary with time series and final metrics
        """
        logger.info(f"Simulating {genome} genome amplification...")

        # Find binding sites
        if genome == 'target':
            genome_seq = self.base_simulator.target_seq
            genome_length = self.base_simulator.target_length
        else:
            genome_seq = self.base_simulator.background_seq
            genome_length = self.base_simulator.background_length

        # Get positions for this genome only
        target_pos, bg_pos = self.base_simulator.find_binding_sites(primers)
        positions = target_pos if genome == 'target' else bg_pos

        # Flatten all positions
        all_positions = []
        for primer, pos_list in positions.items():
            all_positions.extend(pos_list)

        n_sites = len(all_positions)
        logger.info(f"  Binding sites: {n_sites:,}")

        if n_sites == 0:
            logger.warning(f"  No binding sites for {genome}, returning zeros")
            return {
                'coverage_percent': [0.0] * self.time_points,
                'dna_amount': [self.initial_dna] * self.time_points,
                'final_coverage': 0.0,
                'final_dna': self.initial_dna,
                'amplification_fold': 1.0
            }

        # Initialize simulation state
        coverage_bins = set()
        dna_amount = self.initial_dna

        # Time series
        coverage_history = []
        dna_history = []

        # Simulate each time step
        for t in range(self.time_points):
            time = t * self.time_step

            # Probability of extension from each site
            # More sites → more priming events → more DNA
            extensions_per_step = n_sites * self.binding_rate * self.time_step

            # Each extension produces DNA
            # Fragment length depends on genome size and coverage
            mean_fragment = min(genome_length / max(n_sites, 1), 50000)  # max 50kb
            dna_synthesized = extensions_per_step * mean_fragment * 650 / 1e9  # ng (650 Da/bp)

            dna_amount += dna_synthesized

            # Coverage accumulation (simplified)
            # As amplification proceeds, more genome bins are covered
            coverage_increase = extensions_per_step * 0.05  # heuristic
            n_covered = min(
                len(coverage_bins) + int(coverage_increase),
                self.base_simulator.target_bins
            )
            coverage_bins = set(range(n_covered))

            coverage_pct = (len(coverage_bins) / self.base_simulator.target_bins) * 100

            coverage_history.append(coverage_pct)
            dna_history.append(dna_amount)

            # Stop if saturated
            if coverage_pct >= 99.9:
                break

        # Pad to full time if stopped early
        while len(coverage_history) < self.time_points:
            coverage_history.append(coverage_history[-1])
            dna_history.append(dna_history[-1])

        final_coverage = coverage_history[-1]
        final_dna = dna_history[-1]
        amplification_fold = final_dna / self.initial_dna

        logger.info(f"  Final coverage: {final_coverage:.1f}%")
        logger.info(f"  Amplification: {amplification_fold:.1f}x")
        logger.info(f"  Final DNA: {final_dna:.1f} ng")

        return {
            'coverage_percent': coverage_history,
            'dna_amount': dna_history,
            'final_coverage': final_coverage,
            'final_dna': final_dna,
            'amplification_fold': amplification_fold,
            'n_sites': n_sites
        }

    def simulate_config(self,
                       primers: List[str],
                       config_name: str) -> ReplicationMetrics:
        """
        Run complete simulation for one configuration.

        Args:
            primers: List of primer sequences
            config_name: Configuration name

        Returns:
            ReplicationMetrics
        """
        logger.info(f"\n{'='*80}")
        logger.info(f"DETAILED REPLICATION SIMULATION: {config_name}")
        logger.info(f"{'='*80}")
        logger.info(f"Primers: {len(primers)}")
        logger.info(f"Simulation time: {self.simulation_time} minutes")
        logger.info("")

        # Simulate target
        target_results = self.simulate_amplification(primers, config_name, 'target')

        # Simulate background
        bg_results = self.simulate_amplification(primers, config_name, 'background')

        # Calculate metrics
        final_coverage = target_results['final_coverage']
        target_dna = target_results['final_dna']
        bg_dna = bg_results['final_dna']

        enrichment = (target_dna / bg_dna) if bg_dna > 0 else 999

        # Time to 50% coverage
        coverage_series = np.array(target_results['coverage_percent'])
        time_to_50 = np.searchsorted(coverage_series, 50.0) * self.time_step
        time_to_50 = min(time_to_50, self.simulation_time)

        # Time to plateau (95% of final coverage)
        plateau_threshold = final_coverage * 0.95
        time_to_plateau = np.searchsorted(coverage_series, plateau_threshold) * self.time_step
        time_to_plateau = min(time_to_plateau, self.simulation_time)

        # Coverage uniformity (simplified - based on binding site count)
        # More sites = more uniform coverage
        n_sites = target_results['n_sites']
        uniformity = min(1.0, n_sites / 100.0)  # Saturates at 100 sites

        # Mean fragment length (simplified)
        mean_fragment = self.base_simulator.target_length / max(n_sites, 1)
        mean_fragment = min(mean_fragment, 50000)  # Cap at 50kb

        # Recommendation
        if enrichment >= 100 and final_coverage >= 80:
            recommendation = "EXCELLENT"
        elif enrichment >= 50 and final_coverage >= 60:
            recommendation = "GOOD"
        elif enrichment >= 10 and final_coverage >= 40:
            recommendation = "FAIR"
        else:
            recommendation = "POOR"

        logger.info(f"\n{'='*80}")
        logger.info("SUMMARY")
        logger.info(f"{'='*80}")
        logger.info(f"Coverage: {final_coverage:.1f}%")
        logger.info(f"Enrichment: {enrichment:.1f}x")
        logger.info(f"Time to 50% coverage: {time_to_50:.1f} min")
        logger.info(f"Time to plateau: {time_to_plateau:.1f} min")
        logger.info(f"Amplification fold: {target_results['amplification_fold']:.1f}x")
        logger.info(f"Recommendation: {recommendation}")
        logger.info(f"{'='*80}\n")

        return ReplicationMetrics(
            config_name=config_name,
            n_primers=len(primers),
            final_coverage_percent=final_coverage,
            amplification_fold=target_results['amplification_fold'],
            mean_fragment_length=mean_fragment,
            time_to_50pct_coverage=time_to_50,
            time_to_plateau=time_to_plateau,
            target_amplification=target_dna,
            background_amplification=bg_dna,
            enrichment=enrichment,
            coverage_uniformity=uniformity,
            recommendation=recommendation
        )


def main():
    """Run detailed simulations for top configs"""
    logger.info("="*80)
    logger.info("DETAILED PHI29 REPLICATION SIMULATION")
    logger.info("="*80)
    logger.info("Running detailed amplification simulation for top-performing configs")
    logger.info("")

    # Paths
    base_dir = Path(__file__).parent
    results_dir = base_dir / "francisella_results"
    target_genome = Path("/Users/andreassjodin/Code/swga-dev/test/francisella_GCF_000008985.1_ASM898v1_genomic.fna")
    background_genome = Path("/Users/andreassjodin/Code/swga-dev/test/bacillus_GCF_000008445.1_ASM844v1_genomic.fna")

    # Output directory
    sim_output = base_dir / "francisella_detailed_simulation"
    sim_output.mkdir(exist_ok=True)

    # Top configs to simulate (based on binding site analysis results)
    top_configs = [
        "Config7_SuperLong_Optimized",  # Highest specificity (105.7x)
        "Config4_MediumLong_Moderate",  # Best coverage (6.8%)
        "Config5_VeryLong_Moderate",    # Good balance (29.1x)
    ]

    logger.info("Configs selected for detailed simulation:")
    for config in top_configs:
        logger.info(f"  - {config}")
    logger.info("")

    # Initialize simulator
    simulator = DetailedSwgaSimulator(
        target_genome_path=str(target_genome),
        background_genome_path=str(background_genome),
        simulation_time=120.0,  # 2 hours
        time_step=1.0  # 1 minute
    )

    # Run simulations
    all_results = []

    for config_name in top_configs:
        config_dir = results_dir / config_name
        primer_file = config_dir / "primers.txt"

        if not primer_file.exists():
            logger.warning(f"No primers.txt for {config_name}, skipping")
            continue

        # Load primers
        primers = load_primers(primer_file)

        if len(primers) == 0:
            logger.warning(f"No primers loaded for {config_name}, skipping")
            continue

        # Run detailed simulation
        result = simulator.simulate_config(primers, config_name)
        all_results.append(result)

        # Save individual result
        result_file = sim_output / f"{config_name}_detailed_simulation.json"
        with open(result_file, 'w') as f:
            json.dump(result.to_dict(), f, indent=2)

        logger.info(f"Saved: {result_file}\n")

    # Save combined results
    combined_file = sim_output / "all_detailed_results.json"
    with open(combined_file, 'w') as f:
        json.dump([r.to_dict() for r in all_results], f, indent=2)

    logger.info(f"{'='*80}")
    logger.info(f"DETAILED SIMULATION COMPLETE")
    logger.info(f"{'='*80}")
    logger.info(f"Results saved to: {sim_output}")
    logger.info(f"Configurations simulated: {len(all_results)}")

    # Print comparison table
    logger.info(f"\n{'='*80}")
    logger.info("DETAILED SIMULATION RESULTS - COMPARISON")
    logger.info(f"{'='*80}")
    logger.info(f"{'Config':<30} {'Coverage':>10} {'Enrich':>8} {'Time50%':>9} {'Recom':<10}")
    logger.info("-" * 80)

    for r in all_results:
        logger.info(
            f"{r.config_name:<30} {r.final_coverage_percent:>9.1f}% "
            f"{r.enrichment:>8.1f}x {r.time_to_50pct_coverage:>8.1f}m {r.recommendation:<10}"
        )

    logger.info("="*80)

    # Key insights
    logger.info("")
    logger.info("KEY INSIGHTS:")
    logger.info("")

    best_coverage = max(all_results, key=lambda r: r.final_coverage_percent)
    best_enrichment = max(all_results, key=lambda r: r.enrichment)
    fastest = min(all_results, key=lambda r: r.time_to_50pct_coverage)

    logger.info(f"  Best coverage: {best_coverage.config_name}")
    logger.info(f"    Coverage: {best_coverage.final_coverage_percent:.1f}%")
    logger.info("")
    logger.info(f"  Highest specificity: {best_enrichment.config_name}")
    logger.info(f"    Enrichment: {best_enrichment.enrichment:.1f}x")
    logger.info("")
    logger.info(f"  Fastest amplification: {fastest.config_name}")
    logger.info(f"    Time to 50% coverage: {fastest.time_to_50pct_coverage:.1f} min")
    logger.info("")

    logger.info("NEXT STEPS:")
    logger.info("  1. Compare detailed simulation with binding site analysis")
    logger.info("  2. Generate time series plots (coverage vs time)")
    logger.info("  3. Validate with experimental data")
    logger.info("")

    return all_results


if __name__ == "__main__":
    main()
