#!/usr/bin/env python3
"""
Francisella vs Human Background Simulation

Re-runs the Francisella primer simulation using human genome as background
instead of Bacillus. This provides more realistic specificity estimates for
clinical applications where Francisella would be extracted from human samples.

Author: NeoSWGA Development Team
Date: 2025-11-23
"""

import sys
import logging
import json
from pathlib import Path
from typing import List, Dict
from simulate_francisella_primers import (
    FrancisellaSimulator,
    SimulationMetrics,
    load_primers
)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def main():
    """Run Francisella simulation against human genome background"""
    logger.info("="*80)
    logger.info("FRANCISELLA vs HUMAN GENOME SIMULATION")
    logger.info("="*80)
    logger.info("Comparing primer specificity against realistic clinical background")
    logger.info("")

    # Paths
    base_dir = Path(__file__).parent
    results_dir = base_dir / "francisella_results"
    target_genome = Path("/Users/andreassjodin/Code/swga-dev/test/francisella_GCF_000008985.1_ASM898v1_genomic.fna")

    # Check for human genome
    human_genome_paths = [
        Path("/Users/andreassjodin/Code/swga-dev/test/human_bg/hg38_chr1.fna"),
        Path("/Users/andreassjodin/Code/swga-dev/human_bg/hg38_chr1.fna"),
        Path("~/Code/swga-dev/test/human_bg/hg38_chr1.fna"),
        Path("../test/human_bg/hg38_chr1.fna")
    ]

    human_genome = None
    for path in human_genome_paths:
        if path.expanduser().exists():
            human_genome = path.expanduser()
            break

    if human_genome is None:
        logger.warning("="*80)
        logger.warning("HUMAN GENOME NOT FOUND")
        logger.warning("="*80)
        logger.warning("Searched paths:")
        for path in human_genome_paths:
            logger.warning(f"  - {path}")
        logger.warning("")
        logger.warning("Human genome is required for realistic clinical validation.")
        logger.warning("Please download hg38 chromosome 1 (or full genome) and place at:")
        logger.warning("  /Users/andreassjodin/Code/swga-dev/test/human_bg/hg38_chr1.fna")
        logger.warning("")
        logger.warning("Download from:")
        logger.warning("  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/")
        logger.warning("")
        logger.warning("Falling back to Bacillus background...")
        human_genome = Path("/Users/andreassjodin/Code/swga-dev/test/bacillus_GCF_000008445.1_ASM844v1_genomic.fna")

    # Output directory
    sim_output = base_dir / "francisella_simulation_human"
    sim_output.mkdir(exist_ok=True)

    # Initialize simulator with HUMAN background
    logger.info(f"Target genome: {target_genome.name}")
    logger.info(f"Background genome: {human_genome.name}")
    logger.info("")

    simulator = FrancisellaSimulator(
        target_genome_path=str(target_genome),
        background_genome_path=str(human_genome),
        bin_size=10000
    )

    # Find all config directories
    config_dirs = sorted([d for d in results_dir.iterdir() if d.is_dir()])
    logger.info(f"Found {len(config_dirs)} configurations\n")

    # Run simulations
    all_results = []

    for config_dir in config_dirs:
        config_name = config_dir.name
        primer_file = config_dir / "primers.txt"

        if not primer_file.exists():
            logger.warning(f"No primers.txt for {config_name}, skipping")
            continue

        # Load primers
        primers = load_primers(primer_file)

        if len(primers) == 0:
            logger.warning(f"No primers loaded for {config_name}, skipping")
            continue

        # Run simulation
        result = simulator.simulate_config(primers, config_name)
        all_results.append(result)

        # Save individual result
        result_file = sim_output / f"{config_name}_simulation_human.json"
        with open(result_file, 'w') as f:
            json.dump(result.to_dict(), f, indent=2)

        logger.info(f"Saved: {result_file}")

    # Save combined results
    combined_file = sim_output / "all_results_human.json"
    with open(combined_file, 'w') as f:
        json.dump([r.to_dict() for r in all_results], f, indent=2)

    logger.info(f"\n{'='*80}")
    logger.info(f"SIMULATION COMPLETE")
    logger.info(f"{'='*80}")
    logger.info(f"Results saved to: {sim_output}")
    logger.info(f"Configurations simulated: {len(all_results)}")

    # Print comparison table
    logger.info(f"\n{'='*80}")
    logger.info("FRANCISELLA vs HUMAN BACKGROUND - SUMMARY")
    logger.info(f"{'='*80}")
    logger.info(f"{'Config':<25} {'Primers':>8} {'Coverage':>10} {'Enrich':>8} {'Recom':<10}")
    logger.info("-" * 80)

    for r in all_results:
        logger.info(
            f"{r.config_name:<25} {r.n_primers:>8} {r.target_coverage_percent:>9.1f}% "
            f"{r.enrichment:>8.1f}x {r.recommendation:<10}"
        )

    logger.info("="*80)

    # Highlight key findings
    logger.info("")
    logger.info("KEY FINDINGS:")
    logger.info("")

    best_enrichment = max(all_results, key=lambda r: r.enrichment)
    best_coverage = max(all_results, key=lambda r: r.target_coverage_percent)

    logger.info(f"  Highest specificity: {best_enrichment.config_name}")
    logger.info(f"    Enrichment: {best_enrichment.enrichment:.1f}x")
    logger.info(f"    Mean primer length: {best_enrichment.mean_primer_length:.1f}bp")
    logger.info("")
    logger.info(f"  Best coverage: {best_coverage.config_name}")
    logger.info(f"    Coverage: {best_coverage.target_coverage_percent:.1f}%")
    logger.info(f"    Enrichment: {best_coverage.enrichment:.1f}x")
    logger.info("")

    logger.info("NEXT STEPS:")
    logger.info("  1. Generate comparison report (Bacillus vs Human background)")
    logger.info("  2. Run detailed Phi29 replication simulation for top configs")
    logger.info("  3. Experimental wet-lab validation")
    logger.info("")

    return all_results


if __name__ == "__main__":
    main()
