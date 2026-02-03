#!/usr/bin/env python3
"""
Francisella Primer Performance Simulation

Simulates theoretical performance of all 8 Francisella primer configurations
to validate the terminal Tm bug fix and assess additive chemistry benefits.

Author: NeoSWGA Development Team
Date: 2025-11-23
"""

import sys
import logging
import json
import time
from pathlib import Path
from typing import List, Dict, Tuple
from dataclasses import dataclass, asdict
from Bio import SeqIO
import numpy as np
import h5py

# Add neoswga to path
sys.path.insert(0, str(Path(__file__).parent))

from src import string_search, utility

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class ConfigInfo:
    """Information about a Francisella test configuration"""
    name: str
    description: str
    primer_length_range: Tuple[int, int]
    qa_stringency: str
    swga_mode: bool
    additive_chemistry: Dict
    expected_primers: int


@dataclass
class SimulationMetrics:
    """Simulation results for one configuration"""
    config_name: str
    n_primers: int

    # Binding metrics
    target_binding_sites: int
    background_binding_sites: int
    target_frequency: float  # Sites per kb
    background_frequency: float

    # Coverage metrics (simplified - bins covered)
    target_coverage_bins: int  # Number of 10kb bins with ≥1 binding site
    target_coverage_percent: float

    # Specificity metrics
    enrichment: float  # target_freq / bg_freq
    specificity_score: float  # 0-1

    # Primer characteristics
    mean_primer_length: float
    mean_sites_per_primer_target: float
    mean_sites_per_primer_background: float

    # Quality summary
    recommendation: str  # EXCELLENT, GOOD, FAIR, POOR

    def to_dict(self) -> Dict:
        return asdict(self)


class FrancisellaSimulator:
    """
    Simplified simulator for Francisella primers.

    Uses binding site analysis to estimate coverage and specificity
    without full replication simulation (much faster).
    """

    def __init__(self,
                 target_genome_path: str,
                 background_genome_path: str,
                 bin_size: int = 10000):
        """
        Initialize simulator.

        Args:
            target_genome_path: Path to Francisella genome
            background_genome_path: Path to background genome
            bin_size: Bin size for coverage analysis (bp)
        """
        self.target_genome_path = target_genome_path
        self.background_genome_path = background_genome_path
        self.bin_size = bin_size

        # Load genome sequences
        logger.info("Loading genomes...")
        self.target_seq = self._load_genome(target_genome_path)
        self.background_seq = self._load_genome(background_genome_path)

        self.target_length = len(self.target_seq)
        self.background_length = len(self.background_seq)

        logger.info(f"Target genome: {self.target_length:,} bp")
        logger.info(f"Background genome: {self.background_length:,} bp")

        self.target_bins = (self.target_length // bin_size) + 1

    def _load_genome(self, genome_path: str) -> str:
        """Load genome sequence from FASTA"""
        seq = ""
        with open(genome_path) as f:
            for record in SeqIO.parse(f, "fasta"):
                seq += str(record.seq).upper()
        return seq

    def find_binding_sites(self, primers: List[str]) -> Tuple[Dict, Dict]:
        """
        Find all binding sites for primers in target and background.

        Returns:
            (target_positions, background_positions)
            where each is {primer: [pos1, pos2, ...]}
        """
        logger.info(f"Finding binding sites for {len(primers)} primers...")

        # Use existing string_search functionality
        target_positions = {}
        background_positions = {}

        # Group primers by length for efficiency
        primer_by_length = {}
        for p in primers:
            length = len(p)
            if length not in primer_by_length:
                primer_by_length[length] = []
            primer_by_length[length].append(p)

        # Find positions for each length group
        for length, primer_list in primer_by_length.items():
            logger.info(f"  Searching {len(primer_list)} primers of length {length}...")

            # Search target
            target_dict = self._find_positions_in_sequence(
                primer_list, self.target_seq
            )
            target_positions.update(target_dict)

            # Search background
            bg_dict = self._find_positions_in_sequence(
                primer_list, self.background_seq
            )
            background_positions.update(bg_dict)

        return target_positions, background_positions

    def _find_positions_in_sequence(self, primers: List[str], sequence: str) -> Dict:
        """Find all positions of primers in a sequence"""
        positions = {}
        for primer in primers:
            positions[primer] = []
            # Search forward strand
            pos = 0
            while True:
                pos = sequence.find(primer, pos)
                if pos == -1:
                    break
                positions[primer].append(pos)
                pos += 1

            # Search reverse complement
            rev_comp = str(utility.reverse_complement(primer))
            pos = 0
            while True:
                pos = sequence.find(rev_comp, pos)
                if pos == -1:
                    break
                positions[primer].append(pos)
                pos += 1

        return positions

    def calculate_coverage(self, positions: Dict[str, List[int]]) -> Tuple[int, float]:
        """
        Calculate coverage in terms of bins covered.

        Args:
            positions: {primer: [pos1, pos2, ...]}

        Returns:
            (bins_covered, coverage_percent)
        """
        covered_bins = set()

        for primer, pos_list in positions.items():
            for pos in pos_list:
                bin_idx = pos // self.bin_size
                covered_bins.add(bin_idx)

        coverage_percent = (len(covered_bins) / self.target_bins) * 100

        return len(covered_bins), coverage_percent

    def simulate_config(self,
                       primers: List[str],
                       config_name: str) -> SimulationMetrics:
        """
        Simulate one configuration.

        Args:
            primers: List of primer sequences
            config_name: Configuration name

        Returns:
            SimulationMetrics
        """
        logger.info(f"\n{'='*80}")
        logger.info(f"Simulating: {config_name}")
        logger.info(f"{'='*80}")
        logger.info(f"Primers: {len(primers)}")

        # Find binding sites
        target_pos, bg_pos = self.find_binding_sites(primers)

        # Count total binding sites
        target_sites = sum(len(pos) for pos in target_pos.values())
        bg_sites = sum(len(pos) for pos in bg_pos.values())

        logger.info(f"Target binding sites: {target_sites:,}")
        logger.info(f"Background binding sites: {bg_sites:,}")

        # Calculate frequencies (per kb)
        target_freq = (target_sites / self.target_length) * 1000
        bg_freq = (bg_sites / self.background_length) * 1000 if bg_sites > 0 else 0.0001

        # Calculate enrichment
        enrichment = target_freq / bg_freq if bg_freq > 0 else 999

        logger.info(f"Target frequency: {target_freq:.2f} sites/kb")
        logger.info(f"Background frequency: {bg_freq:.2f} sites/kb")
        logger.info(f"Enrichment: {enrichment:.1f}x")

        # Calculate coverage
        bins_covered, coverage_pct = self.calculate_coverage(target_pos)
        logger.info(f"Coverage: {coverage_pct:.1f}% ({bins_covered}/{self.target_bins} bins)")

        # Calculate mean sites per primer
        mean_sites_target = target_sites / len(primers) if len(primers) > 0 else 0
        mean_sites_bg = bg_sites / len(primers) if len(primers) > 0 else 0

        # Calculate mean primer length
        mean_length = np.mean([len(p) for p in primers])

        # Calculate specificity score (0-1)
        # Higher enrichment = higher score, scaled logarithmically
        specificity_score = min(1.0, np.log10(enrichment + 1) / 3.0)  # log10(1000) = 3

        # Determine recommendation
        if enrichment >= 100 and coverage_pct >= 60:
            recommendation = "EXCELLENT"
        elif enrichment >= 50 and coverage_pct >= 40:
            recommendation = "GOOD"
        elif enrichment >= 10 and coverage_pct >= 20:
            recommendation = "FAIR"
        else:
            recommendation = "POOR"

        logger.info(f"Recommendation: {recommendation}")

        return SimulationMetrics(
            config_name=config_name,
            n_primers=len(primers),
            target_binding_sites=target_sites,
            background_binding_sites=bg_sites,
            target_frequency=target_freq,
            background_frequency=bg_freq,
            target_coverage_bins=bins_covered,
            target_coverage_percent=coverage_pct,
            enrichment=enrichment,
            specificity_score=specificity_score,
            mean_primer_length=mean_length,
            mean_sites_per_primer_target=mean_sites_target,
            mean_sites_per_primer_background=mean_sites_bg,
            recommendation=recommendation
        )


def load_primers(primer_file: Path) -> List[str]:
    """Load primers from config output file"""
    primers = []
    with open(primer_file) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                # Parse: primer_sequence<tab>qa_score
                parts = line.split('\t')
                if parts:
                    primers.append(parts[0])
    return primers


def main():
    """Main simulation orchestration"""
    logger.info("="*80)
    logger.info("FRANCISELLA PRIMER PERFORMANCE SIMULATION")
    logger.info("="*80)

    # Paths
    base_dir = Path(__file__).parent
    results_dir = base_dir / "francisella_results"
    target_genome = Path("/Users/andreassjodin/Code/swga-dev/test/francisella_GCF_000008985.1_ASM898v1_genomic.fna")
    background_genome = Path("/Users/andreassjodin/Code/swga-dev/test/bacillus_GCF_000008445.1_ASM844v1_genomic.fna")

    # Output directory
    sim_output = base_dir / "francisella_simulation"
    sim_output.mkdir(exist_ok=True)

    # Initialize simulator
    simulator = FrancisellaSimulator(
        target_genome_path=str(target_genome),
        background_genome_path=str(background_genome),
        bin_size=10000
    )

    # Find all config directories
    config_dirs = sorted([d for d in results_dir.iterdir() if d.is_dir()])
    logger.info(f"\nFound {len(config_dirs)} configurations")

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
        result_file = sim_output / f"{config_name}_simulation.json"
        with open(result_file, 'w') as f:
            json.dump(result.to_dict(), f, indent=2)

        logger.info(f"Saved: {result_file}")

    # Save combined results
    combined_file = sim_output / "all_results.json"
    with open(combined_file, 'w') as f:
        json.dump([r.to_dict() for r in all_results], f, indent=2)

    logger.info(f"\n{'='*80}")
    logger.info(f"SIMULATION COMPLETE")
    logger.info(f"{'='*80}")
    logger.info(f"Results saved to: {sim_output}")
    logger.info(f"Configurations simulated: {len(all_results)}")

    # Print summary table
    logger.info(f"\n{'='*80}")
    logger.info("SUMMARY TABLE")
    logger.info(f"{'='*80}")
    logger.info(f"{'Config':<25} {'Primers':>8} {'Coverage':>10} {'Enrich':>8} {'Recom':<10}")
    logger.info("-" * 80)

    for r in all_results:
        logger.info(
            f"{r.config_name:<25} {r.n_primers:>8} {r.target_coverage_percent:>9.1f}% "
            f"{r.enrichment:>8.1f}x {r.recommendation:<10}"
        )

    logger.info("="*80)

    return all_results


if __name__ == "__main__":
    main()
