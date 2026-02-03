#!/usr/bin/env python3
"""
Generate synthetic training data for enhanced primer scoring model.

Uses the stochastic simulator (Gillespie algorithm) to generate
labeled training data with actual amplification outcomes, rather
than relying solely on empirical data.

This enables:
1. Training on genomes without experimental data
2. Exploring edge cases (extreme GC, large genomes)
3. Generating enough data to train 120+ feature models

Usage:
    python scripts/generate_training_data.py \
        --fg-genome target.fasta \
        --bg-genome background.fasta \
        --output training_data.csv \
        --num-primers 500 \
        --simulation-time 3600
"""

import argparse
import json
import logging
import os
import sys
import time
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import random

import numpy as np
import pandas as pd
from Bio import SeqIO

# Add neoswga to path if running as script
script_dir = Path(__file__).parent.absolute()
sys.path.insert(0, str(script_dir.parent))

from neoswga.core.advanced_features import AdvancedFeatureEngineer
from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.stochastic_simulator import GillespieSimulator, ReactionParameters
from neoswga.core.amplicon_network import AmpliconNetwork
from neoswga.core.thermodynamics import gc_content

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def generate_random_primers(genome: str, k: int, num_primers: int,
                           min_gc: float = 0.3, max_gc: float = 0.7) -> List[str]:
    """
    Generate random primers from genome with GC filtering.

    Args:
        genome: Genome sequence
        k: Primer length
        num_primers: Number of primers to generate
        min_gc: Minimum GC content
        max_gc: Maximum GC content

    Returns:
        List of primer sequences
    """
    primers = set()
    attempts = 0
    max_attempts = num_primers * 100

    genome_len = len(genome)

    while len(primers) < num_primers and attempts < max_attempts:
        pos = random.randint(0, genome_len - k)
        seq = genome[pos:pos + k].upper()

        if 'N' in seq:
            attempts += 1
            continue

        gc = gc_content(seq)
        if min_gc <= gc <= max_gc:
            primers.add(seq)

        attempts += 1

    return list(primers)


def find_primer_positions(genome: str, primer: str) -> Dict[str, List[int]]:
    """
    Find all occurrences of primer in genome (both strands).

    Args:
        genome: Genome sequence
        primer: Primer sequence

    Returns:
        Dict with 'forward' and 'reverse' position lists
    """
    from neoswga.core.thermodynamics import reverse_complement

    forward_positions = []
    reverse_positions = []

    primer_rc = reverse_complement(primer)

    # Forward strand
    start = 0
    while True:
        pos = genome.find(primer, start)
        if pos == -1:
            break
        forward_positions.append(pos)
        start = pos + 1

    # Reverse strand (search for reverse complement)
    start = 0
    while True:
        pos = genome.find(primer_rc, start)
        if pos == -1:
            break
        reverse_positions.append(pos)
        start = pos + 1

    return {
        'forward': forward_positions,
        'reverse': reverse_positions
    }


class SimplifiedNetwork:
    """Simplified amplification network for training data generation."""

    def __init__(self, primers: List[str], primer_positions: Dict[str, Dict],
                 genome_length: int, max_extension: int = 70000):
        self.primers = primers
        self.primer_positions = primer_positions
        self.genome_length = genome_length
        self.max_extension = max_extension

        # Collect all binding sites
        self.binding_sites = set()
        for primer in primers:
            positions = primer_positions.get(primer, {})
            for pos in positions.get('forward', []):
                self.binding_sites.add(pos)
            for pos in positions.get('reverse', []):
                self.binding_sites.add(pos)

        self.graph = self._build_graph()

    def _build_graph(self):
        """Build simple adjacency representation."""
        import networkx as nx
        G = nx.DiGraph()

        sites = sorted(self.binding_sites)
        for site in sites:
            G.add_node(site)

        # Add edges for sites within extension distance
        for i, site1 in enumerate(sites):
            for site2 in sites[i+1:]:
                distance = site2 - site1
                if distance <= self.max_extension:
                    G.add_edge(site1, site2, distance=distance)
                else:
                    break

        return G

    def largest_component_size(self) -> int:
        """Get size of largest connected component."""
        import networkx as nx
        if len(self.graph.nodes()) == 0:
            return 0
        components = list(nx.weakly_connected_components(self.graph))
        if not components:
            return 0
        return max(len(c) for c in components)

    def average_component_size(self) -> float:
        """Get average component size."""
        import networkx as nx
        if len(self.graph.nodes()) == 0:
            return 0.0
        components = list(nx.weakly_connected_components(self.graph))
        if not components:
            return 0.0
        return np.mean([len(c) for c in components])

    def predict_amplification_fold(self) -> float:
        """Predict amplification fold based on network structure."""
        largest = self.largest_component_size()
        if largest == 0:
            return 1.0
        return 2 ** (np.log2(largest + 1))


def simulate_primer_set(primers: List[str],
                       fg_genome: str, bg_genome: str,
                       fg_positions: Dict, bg_positions: Dict,
                       simulation_time: float = 1800.0) -> Dict:
    """
    Simulate amplification of a primer set and return metrics.

    Args:
        primers: List of primers
        fg_genome: Target genome sequence
        bg_genome: Background genome sequence
        fg_positions: Position dict for target
        bg_positions: Position dict for background
        simulation_time: Simulation time in seconds

    Returns:
        Dict with amplification metrics
    """
    # Build networks
    fg_network = SimplifiedNetwork(primers, fg_positions, len(fg_genome))
    bg_network = SimplifiedNetwork(primers, bg_positions, len(bg_genome))

    # Run simulation
    try:
        simulator = GillespieSimulator(
            primers, fg_network, bg_network,
            ReactionParameters(temperature=30.0)
        )
        history = simulator.simulate(
            max_time=simulation_time,
            sample_interval=simulation_time / 10
        )

        if history:
            final = history[-1]
            return {
                'fg_amplification': final.get('fg_amplification', 1.0),
                'bg_amplification': final.get('bg_amplification', 1.0),
                'enrichment': final.get('enrichment', 1.0),
                'simulation_success': True
            }
    except Exception as e:
        logger.warning(f"Simulation failed: {e}")

    # Fallback to network-based prediction
    fg_amp = fg_network.predict_amplification_fold()
    bg_amp = bg_network.predict_amplification_fold()

    return {
        'fg_amplification': fg_amp,
        'bg_amplification': bg_amp,
        'enrichment': fg_amp / (bg_amp + 0.001),
        'simulation_success': False
    }


def generate_training_data(fg_genome_path: str,
                          bg_genome_path: str,
                          output_path: str,
                          num_primers: int = 500,
                          primer_length: int = 10,
                          simulation_time: float = 1800.0,
                          num_sets: int = 100,
                          set_size: int = 6,
                          conditions: Optional[ReactionConditions] = None) -> pd.DataFrame:
    """
    Generate complete training dataset with features and labels.

    Args:
        fg_genome_path: Path to target genome FASTA
        bg_genome_path: Path to background genome FASTA
        output_path: Output CSV path
        num_primers: Number of candidate primers to generate
        primer_length: Length of primers
        simulation_time: Simulation time per set (seconds)
        num_sets: Number of primer sets to simulate
        set_size: Size of each primer set
        conditions: Reaction conditions

    Returns:
        DataFrame with features and labels
    """
    logger.info("Loading genomes...")

    # Load genomes
    fg_record = next(SeqIO.parse(fg_genome_path, "fasta"))
    fg_genome = str(fg_record.seq).upper()
    logger.info(f"Target genome: {len(fg_genome):,} bp")

    bg_record = next(SeqIO.parse(bg_genome_path, "fasta"))
    bg_genome = str(bg_record.seq).upper()
    logger.info(f"Background genome: {len(bg_genome):,} bp")

    if conditions is None:
        conditions = ReactionConditions(temp=30.0, na_conc=50.0)

    # Generate candidate primers
    logger.info(f"Generating {num_primers} candidate primers...")
    primers = generate_random_primers(fg_genome, primer_length, num_primers)
    logger.info(f"Generated {len(primers)} unique primers")

    # Find positions for all primers
    logger.info("Finding primer positions...")
    fg_positions = {}
    bg_positions = {}

    for primer in primers:
        fg_positions[primer] = find_primer_positions(fg_genome, primer)
        bg_positions[primer] = find_primer_positions(bg_genome, primer)

    # Calculate features for all primers
    logger.info("Engineering features...")
    feature_engineer = AdvancedFeatureEngineer(
        fg_genome, conditions, fg_positions
    )
    features_df = feature_engineer.engineer_features(primers)
    logger.info(f"Features shape: {features_df.shape}")

    # Generate training samples by simulating primer sets
    logger.info(f"Simulating {num_sets} primer sets...")
    all_results = []

    for i in range(num_sets):
        if (i + 1) % 10 == 0:
            logger.info(f"  Progress: {i + 1}/{num_sets}")

        # Sample a random set of primers
        selected_primers = random.sample(primers, min(set_size, len(primers)))

        # Simulate this set
        result = simulate_primer_set(
            selected_primers, fg_genome, bg_genome,
            fg_positions, bg_positions, simulation_time
        )

        # Store results for each primer in the set
        for primer in selected_primers:
            all_results.append({
                'primer': primer,
                'set_id': i,
                'set_size': len(selected_primers),
                'fg_amplification': result['fg_amplification'],
                'bg_amplification': result['bg_amplification'],
                'enrichment': result['enrichment'],
                'simulation_success': result['simulation_success']
            })

    # Create results DataFrame
    results_df = pd.DataFrame(all_results)

    # Aggregate results per primer (average across sets)
    primer_labels = results_df.groupby('primer').agg({
        'fg_amplification': 'mean',
        'bg_amplification': 'mean',
        'enrichment': 'mean',
        'set_id': 'count'
    }).rename(columns={
        'set_id': 'num_sets',
        'fg_amplification': 'mean_fg_amplification',
        'bg_amplification': 'mean_bg_amplification',
        'enrichment': 'mean_enrichment'
    }).reset_index()

    # Merge features with labels
    training_df = features_df.merge(primer_labels, on='primer', how='left')

    # Add derived labels
    training_df['log_enrichment'] = np.log10(training_df['mean_enrichment'].clip(lower=0.001))
    training_df['quality_score'] = (
        training_df['mean_enrichment'].rank(pct=True) * 0.5 +
        (1 / training_df['mean_bg_amplification'].clip(lower=0.001)).rank(pct=True) * 0.5
    )

    # Save
    training_df.to_csv(output_path, index=False)
    logger.info(f"Saved training data to {output_path}")
    logger.info(f"Final dataset: {len(training_df)} primers x {len(training_df.columns)} features/labels")

    return training_df


def main():
    parser = argparse.ArgumentParser(
        description='Generate synthetic training data for enhanced primer scoring'
    )

    parser.add_argument('--fg-genome', required=True,
                       help='Path to target genome FASTA')
    parser.add_argument('--bg-genome', required=True,
                       help='Path to background genome FASTA')
    parser.add_argument('--output', default='training_data.csv',
                       help='Output CSV path')
    parser.add_argument('--num-primers', type=int, default=500,
                       help='Number of candidate primers')
    parser.add_argument('--primer-length', type=int, default=10,
                       help='Primer length')
    parser.add_argument('--simulation-time', type=float, default=1800.0,
                       help='Simulation time per set (seconds)')
    parser.add_argument('--num-sets', type=int, default=100,
                       help='Number of primer sets to simulate')
    parser.add_argument('--set-size', type=int, default=6,
                       help='Primers per set')
    parser.add_argument('--temperature', type=float, default=30.0,
                       help='Reaction temperature (C)')
    parser.add_argument('--na-conc', type=float, default=50.0,
                       help='Sodium concentration (mM)')

    args = parser.parse_args()

    conditions = ReactionConditions(
        temp=args.temperature,
        na_conc=args.na_conc
    )

    start_time = time.time()

    generate_training_data(
        fg_genome_path=args.fg_genome,
        bg_genome_path=args.bg_genome,
        output_path=args.output,
        num_primers=args.num_primers,
        primer_length=args.primer_length,
        simulation_time=args.simulation_time,
        num_sets=args.num_sets,
        set_size=args.set_size,
        conditions=conditions
    )

    elapsed = time.time() - start_time
    logger.info(f"Total time: {elapsed/60:.1f} minutes")


if __name__ == '__main__':
    main()
