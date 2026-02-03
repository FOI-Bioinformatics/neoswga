#!/usr/bin/env python3
"""
SWGA Comprehensive Simulator

Validates primer sets by simulating SWGA amplification on target and background genomes.

Three simulation modes:
1. Fast (~1 min): Network-based prediction
2. Detailed (~10 min): Agent-based replication simulation
3. Validation (~30 min): Stochastic kinetics validation

Usage:
    from neoswga.core.swga_simulator import SwgaSimulator

    sim = SwgaSimulator(
        primers=['ATCGATCG', 'GCTAGCTA', ...],
        fg_genome='target.fasta',
        bg_genome='background.fasta',
        fg_positions_h5='target_pos.h5',
        bg_positions_h5='background_pos.h5'
    )

    result = sim.simulate_detailed()
    print(f"Coverage: {result['target']['coverage']:.1%}")
    print(f"Specificity: {result['specificity']['enrichment']:.0f}×")
    print(f"Recommendation: {result['recommendation']}")
"""

import numpy as np
import h5py
import time
import logging
from typing import List, Dict, Optional, Tuple
from pathlib import Path
from Bio import SeqIO
from dataclasses import dataclass, asdict

logger = logging.getLogger(__name__)


@dataclass
class SimulationResult:
    """Results from SWGA simulation"""
    mode: str  # 'fast', 'detailed', or 'validation'
    runtime: float  # seconds

    # Target genome results
    target_coverage: float  # Fraction of genome covered
    target_uniformity: float  # 1 - Gini coefficient
    target_amplification: float  # Predicted fold amplification
    target_gaps: List[Dict]  # Under-covered regions

    # Background genome results
    background_coverage: float
    background_amplification: float

    # Specificity metrics
    enrichment: float  # Target/background ratio
    specificity_score: float  # 0-1 normalized score

    # Quality metrics
    composite_score: float  # Overall quality (0-1)
    recommendation: str  # 'EXCELLENT', 'GOOD', 'FAIR', 'POOR'
    confidence: float  # Confidence in prediction (0-1)

    # Primer analysis
    primer_contributions: Dict  # Contribution of each primer

    # Detailed results (mode-dependent)
    details: Dict  # Additional mode-specific data

    def to_dict(self) -> Dict:
        """Convert to dictionary"""
        return asdict(self)


class SwgaSimulator:
    """
    Comprehensive SWGA simulation tool.

    Simulates primer set performance on target and background genomes using
    network analysis, agent-based replication, or stochastic kinetics.
    """

    def __init__(self,
                 primers: List[str],
                 fg_genome: str,
                 bg_genome: str,
                 fg_positions_h5: str,
                 bg_positions_h5: str,
                 reaction_conditions: Optional[Dict] = None,
                 bin_size: int = 10000,
                 max_extension: int = 70000):
        """
        Initialize simulator.

        Args:
            primers: List of primer sequences
            fg_genome: Path to target genome FASTA
            bg_genome: Path to background genome FASTA
            fg_positions_h5: Path to target positions HDF5
            bg_positions_h5: Path to background positions HDF5
            reaction_conditions: Optional dict with temp, polymerase, etc.
            bin_size: Size of bins for coverage analysis (default 10kb)
            max_extension: Max Phi29 extension distance (default 70kb)
        """
        self.primers = primers
        self.fg_genome_path = Path(fg_genome)
        self.bg_genome_path = Path(bg_genome)
        self.fg_positions_h5 = Path(fg_positions_h5)
        self.bg_positions_h5 = Path(bg_positions_h5)
        self.bin_size = bin_size
        self.max_extension = max_extension

        # Reaction conditions
        self.conditions = reaction_conditions or {
            'temperature': 30,
            'polymerase': 'phi29',
            'duration': 3600
        }

        # Load genomes
        logger.info("Loading genomes...")
        self.fg_genome, self.fg_length, self.fg_gc = self._load_genome(fg_genome)
        self.bg_genome, self.bg_length, self.bg_gc = self._load_genome(bg_genome)

        # Extract genome names from file paths for HDF5 lookup
        # Assumes genome files have format: genusspecies_*.fna or genus_*.fna
        fg_genome_name = self._extract_genome_name(fg_genome)
        bg_genome_name = self._extract_genome_name(bg_genome)

        # Load positions
        logger.info("Loading primer positions...")
        logger.info(f"  Target genome: {fg_genome_name}")
        logger.info(f"  Background genome: {bg_genome_name}")
        self.fg_positions = self._load_positions(fg_positions_h5, primers, fg_genome_name)
        self.bg_positions = self._load_positions(bg_positions_h5, primers, bg_genome_name)

        logger.info(f"Initialized simulator:")
        logger.info(f"  Primers: {len(primers)}")
        logger.info(f"  Target: {self.fg_length:,} bp ({self.fg_gc:.1%} GC)")
        logger.info(f"  Background: {self.bg_length:,} bp ({self.bg_gc:.1%} GC)")

    def _extract_genome_name(self, genome_path: str) -> str:
        """
        Extract genome name from file path.

        Examples:
            yersinia_GCF_000222975.1_ASM22297v1_genomic.fna → yersinia
            GCF_000001405.40_GRCh38.p14_genomic.fna → human (special case)
        """
        filename = Path(genome_path).name
        # Remove extension
        name = filename.split('.')[0].split('_')[0].lower()

        # Special cases
        if 'grch' in filename.lower() or 'human' in filename.lower():
            return 'human'

        return name

    def _load_genome(self, fasta_path: str) -> Tuple[str, int, float]:
        """Load genome sequence and calculate stats"""
        records = list(SeqIO.parse(fasta_path, 'fasta'))
        sequence = ''.join(str(rec.seq).upper() for rec in records)
        length = len(sequence)

        gc_count = sequence.count('G') + sequence.count('C')
        gc_content = gc_count / length if length > 0 else 0

        return sequence, length, gc_content

    def _load_positions(self, h5_path: Path, primers: List[str], genome_name: str = None) -> Dict:
        """
        Load primer positions from HDF5.

        Args:
            h5_path: Path to HDF5 file
            primers: List of primer sequences
            genome_name: Specific genome group to load (if None, uses first group)
        """
        positions = {}

        with h5py.File(h5_path, 'r') as f:
            # Get genome group
            if genome_name is None:
                genome_name = list(f.keys())[0]

            if genome_name not in f:
                logger.warning(f"Genome '{genome_name}' not found in HDF5, using first genome")
                genome_name = list(f.keys())[0]

            genome_group = f[genome_name]

            for primer in primers:
                if primer in genome_group:
                    # Load forward strand positions
                    pos_fwd = np.array(genome_group[primer]['+'])
                    pos_rev = np.array(genome_group[primer]['-'])

                    positions[primer] = {
                        '+': pos_fwd,
                        '-': pos_rev,
                        'forward': pos_fwd,
                        'reverse': pos_rev
                    }
                else:
                    # Primer not found
                    positions[primer] = {
                        '+': np.array([]),
                        '-': np.array([]),
                        'forward': np.array([]),
                        'reverse': np.array([])
                    }

        return positions

    def simulate_fast(self) -> SimulationResult:
        """
        Fast coverage-based simulation (~1 second).

        Uses genome binning to estimate coverage, similar to dominating-set method.
        Good for quick validation.

        Returns:
            SimulationResult with coverage-based predictions
        """
        logger.info("Running fast coverage-based simulation...")
        start_time = time.time()

        # Calculate target genome coverage
        fg_coverage, fg_covered_bins = self._calculate_bin_coverage(
            self.fg_positions, self.fg_length
        )
        fg_uniformity = self._estimate_uniformity(fg_covered_bins, self.fg_length)

        # Calculate background genome coverage
        bg_coverage, bg_covered_bins = self._calculate_bin_coverage(
            self.bg_positions, self.bg_length
        )

        # Count total binding sites for amplification estimate
        fg_total_sites = sum(len(pos['+']) + len(pos['-']) for pos in self.fg_positions.values())
        bg_total_sites = sum(len(pos['+']) + len(pos['-']) for pos in self.bg_positions.values())

        # Estimate amplification based on sites and coverage
        # Higher coverage + more sites = better amplification
        target_amplification = fg_total_sites * fg_coverage * 10  # Scaling factor
        background_amplification = bg_total_sites * bg_coverage * 10

        # Calculate enrichment (accounting for genome size difference)
        genome_size_ratio = self.bg_length / self.fg_length
        enrichment = (target_amplification / (background_amplification + 1)) * genome_size_ratio

        # Identify gaps
        gaps = self._identify_coverage_gaps(fg_covered_bins, self.fg_length)

        # Calculate composite score
        result = self._calculate_composite_score(
            target_coverage=fg_coverage,
            target_uniformity=fg_uniformity,
            enrichment=enrichment,
            mode='fast'
        )

        runtime = time.time() - start_time

        return SimulationResult(
            mode='fast',
            runtime=runtime,
            target_coverage=fg_coverage,
            target_uniformity=fg_uniformity,
            target_amplification=target_amplification,
            target_gaps=gaps,
            background_coverage=bg_coverage,
            background_amplification=background_amplification,
            enrichment=enrichment,
            specificity_score=result['specificity_score'],
            composite_score=result['composite_score'],
            recommendation=result['recommendation'],
            confidence=0.8,  # Good confidence for coverage-based approach
            primer_contributions={},
            details={
                'target_sites': fg_total_sites,
                'background_sites': bg_total_sites,
                'target_bins_covered': len(fg_covered_bins),
                'background_bins_covered': len(bg_covered_bins)
            }
        )

    def simulate_detailed(self, num_replicates: int = 5) -> SimulationResult:
        """
        Detailed agent-based simulation (~10 minutes).

        Simulates replication fork dynamics with GC-dependent speed.
        Best balance of accuracy and speed.

        Args:
            num_replicates: Number of simulation replicates for statistics

        Returns:
            SimulationResult with detailed replication dynamics
        """
        logger.info(f"Running detailed agent-based simulation ({num_replicates} replicates)...")
        start_time = time.time()

        from neoswga.core.replication_simulator import Phi29Simulator, ReactionConditions

        # Set up reaction conditions
        conditions = ReactionConditions(
            temperature=self.conditions.get('temperature', 30),
            polymerase=self.conditions.get('polymerase', 'phi29')
        )

        # Run target genome replicates
        fg_results = []
        for i in range(num_replicates):
            logger.info(f"  Target replicate {i+1}/{num_replicates}...")
            sim = Phi29Simulator(
                primers=self.primers,
                primer_positions=self.fg_positions,
                genome_sequence=self.fg_genome,
                genome_name='target',
                conditions=conditions
            )
            result = sim.run(duration=self.conditions.get('duration', 3600))
            fg_results.append(result)

        # Run background genome replicates
        bg_results = []
        for i in range(num_replicates):
            logger.info(f"  Background replicate {i+1}/{num_replicates}...")
            sim = Phi29Simulator(
                primers=self.primers,
                primer_positions=self.bg_positions,
                genome_sequence=self.bg_genome,
                genome_name='background',
                conditions=conditions
            )
            result = sim.run(duration=self.conditions.get('duration', 3600))
            bg_results.append(result)

        # Aggregate statistics
        target_coverage = np.mean([r.final_coverage_fraction for r in fg_results])
        target_coverage_std = np.std([r.final_coverage_fraction for r in fg_results])

        background_coverage = np.mean([r.final_coverage_fraction for r in bg_results])

        # Estimate amplification from coverage
        target_amplification = target_coverage * self.fg_length
        background_amplification = background_coverage * self.bg_length

        enrichment = target_amplification / (background_amplification + 1e-10)

        # Calculate uniformity (from coverage array)
        target_uniformity = 1.0 - self._calculate_gini(fg_results[0].coverage_array)

        # Identify gaps
        gaps = self._identify_gaps(fg_results[0].coverage_array, self.fg_length)

        # Calculate composite score
        result = self._calculate_composite_score(
            target_coverage=target_coverage,
            target_uniformity=target_uniformity,
            enrichment=enrichment,
            mode='detailed'
        )

        runtime = time.time() - start_time

        return SimulationResult(
            mode='detailed',
            runtime=runtime,
            target_coverage=target_coverage,
            target_uniformity=target_uniformity,
            target_amplification=target_amplification,
            target_gaps=gaps,
            background_coverage=background_coverage,
            background_amplification=background_amplification,
            enrichment=enrichment,
            specificity_score=result['specificity_score'],
            composite_score=result['composite_score'],
            recommendation=result['recommendation'],
            confidence=0.9,  # High confidence for detailed mode
            primer_contributions={},
            details={
                'fg_results': fg_results,
                'bg_results': bg_results,
                'target_coverage_std': target_coverage_std
            }
        )

    def simulate_validation(self) -> SimulationResult:
        """
        Validation with stochastic kinetics (~30 minutes).

        Uses Gillespie algorithm for exact stochastic simulation.
        Highest accuracy but computationally expensive.

        Returns:
            SimulationResult with stochastic kinetics validation
        """
        logger.info("Running validation with stochastic kinetics...")
        logger.warning("Stochastic simulation not yet implemented - using detailed mode")
        return self.simulate_detailed(num_replicates=10)

    def _calculate_bin_coverage(self, positions_dict: Dict, genome_length: int) -> Tuple[float, set]:
        """
        Calculate bin coverage for a genome.

        Returns:
            (coverage_fraction, set_of_covered_bins)
        """
        n_bins = (genome_length + self.bin_size - 1) // self.bin_size
        covered_bins = set()

        # Mark bins as covered based on primer binding positions
        for primer, pos_data in positions_dict.items():
            for strand in ['+', '-']:
                positions = pos_data.get(strand, np.array([]))
                for pos in positions:
                    bin_idx = int(pos) // self.bin_size
                    if 0 <= bin_idx < n_bins:
                        covered_bins.add(bin_idx)

        coverage = len(covered_bins) / n_bins if n_bins > 0 else 0.0
        return coverage, covered_bins

    def _estimate_uniformity(self, covered_bins: set, genome_length: int) -> float:
        """
        Estimate coverage uniformity.

        Uniformity = 1 - Gini coefficient of bin spacing
        Higher values = more uniform distribution
        """
        if len(covered_bins) < 2:
            return 0.0

        # Calculate gaps between covered bins
        sorted_bins = sorted(covered_bins)
        gaps = [sorted_bins[i+1] - sorted_bins[i] for i in range(len(sorted_bins)-1)]

        if len(gaps) == 0:
            return 1.0

        # Calculate Gini coefficient of gaps
        gaps_array = np.array(sorted(gaps))
        n = len(gaps_array)
        index = np.arange(1, n + 1)
        gini = (2 * np.sum(index * gaps_array)) / (n * np.sum(gaps_array)) - (n + 1) / n

        # Uniformity = 1 - Gini (higher is better)
        uniformity = 1.0 - gini
        return np.clip(uniformity, 0.0, 1.0)

    def _identify_coverage_gaps(self, covered_bins: set, genome_length: int) -> List[Dict]:
        """Identify under-covered regions as gaps between covered bins"""
        if len(covered_bins) == 0:
            return []

        sorted_bins = sorted(covered_bins)
        gaps = []

        # Find gaps larger than 5 bins (50kb default)
        min_gap_size = 5

        for i in range(len(sorted_bins) - 1):
            gap_size = sorted_bins[i+1] - sorted_bins[i] - 1
            if gap_size >= min_gap_size:
                start_pos = (sorted_bins[i] + 1) * self.bin_size
                end_pos = sorted_bins[i+1] * self.bin_size
                gaps.append({
                    'start': start_pos,
                    'end': end_pos,
                    'length': end_pos - start_pos,
                    'mean_coverage': 0.0  # No coverage in gaps
                })

        return gaps

    def _calculate_gini(self, coverage_array: np.ndarray) -> float:
        """Calculate Gini coefficient for coverage uniformity"""
        if len(coverage_array) == 0:
            return 0.0

        # Bin coverage
        n_bins = (len(coverage_array) + self.bin_size - 1) // self.bin_size
        bin_coverage = []

        for i in range(n_bins):
            start = i * self.bin_size
            end = min((i+1) * self.bin_size, len(coverage_array))
            bin_coverage.append(coverage_array[start:end].mean())

        # Calculate Gini
        bin_coverage = np.array(sorted(bin_coverage))
        n = len(bin_coverage)
        index = np.arange(1, n + 1)
        gini = (2 * np.sum(index * bin_coverage)) / (n * np.sum(bin_coverage)) - (n + 1) / n

        return gini

    def _identify_gaps(self, coverage_array: np.ndarray, genome_length: int,
                      threshold_percentile: float = 25) -> List[Dict]:
        """Identify under-covered regions"""
        if len(coverage_array) == 0:
            return []

        # Calculate threshold
        threshold = np.percentile(coverage_array[coverage_array > 0], threshold_percentile)

        # Find gaps
        gaps = []
        in_gap = False
        gap_start = 0

        for i, cov in enumerate(coverage_array):
            if cov < threshold and not in_gap:
                gap_start = i
                in_gap = True
            elif cov >= threshold and in_gap:
                gaps.append({
                    'start': gap_start,
                    'end': i,
                    'length': i - gap_start,
                    'mean_coverage': coverage_array[gap_start:i].mean()
                })
                in_gap = False

        # Close final gap if needed
        if in_gap:
            gaps.append({
                'start': gap_start,
                'end': len(coverage_array),
                'length': len(coverage_array) - gap_start,
                'mean_coverage': coverage_array[gap_start:].mean()
            })

        return gaps

    def _calculate_composite_score(self, target_coverage: float,
                                   target_uniformity: float,
                                   enrichment: float,
                                   mode: str) -> Dict:
        """
        Calculate composite quality score.

        Args:
            target_coverage: Fraction of target genome covered (0-1)
            target_uniformity: Coverage uniformity (0-1, higher is better)
            enrichment: Target/background amplification ratio
            mode: Simulation mode

        Returns:
            Dict with scores and recommendation
        """
        # Normalize coverage (0-1)
        coverage_score = np.clip(target_coverage, 0, 1)

        # Uniformity score (already 0-1)
        uniformity_score = np.clip(target_uniformity, 0, 1)

        # Specificity score (log-scale, normalized to 0-1)
        # 10× enrichment → 0.33, 100× → 0.67, 1000× → 1.0
        specificity_score = np.log10(max(enrichment, 1)) / 3
        specificity_score = np.clip(specificity_score, 0, 1)

        # Composite score (weighted average)
        composite = (
            0.40 * coverage_score +      # Coverage most important
            0.25 * uniformity_score +    # Uniformity important for sequencing
            0.35 * specificity_score     # Specificity critical for mixed samples
        )

        # Recommendation thresholds
        if composite >= 0.80:
            recommendation = "EXCELLENT"
        elif composite >= 0.65:
            recommendation = "GOOD"
        elif composite >= 0.45:
            recommendation = "FAIR"
        else:
            recommendation = "POOR"

        return {
            'coverage_score': coverage_score,
            'uniformity_score': uniformity_score,
            'specificity_score': specificity_score,
            'composite_score': composite,
            'recommendation': recommendation
        }

    def generate_report(self, result: SimulationResult, output_path: Optional[str] = None):
        """
        Generate human-readable report.

        Args:
            result: SimulationResult to report
            output_path: Optional path to save report (default: stdout)
        """
        report = []
        report.append("=" * 80)
        report.append("SWGA SIMULATION REPORT")
        report.append("=" * 80)
        report.append("")
        report.append(f"Simulation mode: {result.mode.upper()}")
        report.append(f"Runtime: {result.runtime:.1f} seconds")
        report.append(f"Confidence: {result.confidence:.0%}")
        report.append("")

        report.append("-" * 80)
        report.append("RECOMMENDATION: " + result.recommendation)
        report.append(f"Composite Score: {result.composite_score:.2f} / 1.00")
        report.append("-" * 80)
        report.append("")

        report.append("TARGET GENOME PERFORMANCE")
        report.append(f"  Coverage: {result.target_coverage:.1%}")
        report.append(f"  Uniformity: {result.target_uniformity:.1%}")
        report.append(f"  Amplification: {result.target_amplification:.2e}×")
        report.append(f"  Number of gaps: {len(result.target_gaps)}")
        report.append("")

        report.append("BACKGROUND GENOME")
        report.append(f"  Coverage: {result.background_coverage:.1%}")
        report.append(f"  Amplification: {result.background_amplification:.2e}×")
        report.append("")

        report.append("SPECIFICITY")
        report.append(f"  Enrichment: {result.enrichment:.1f}×")
        report.append(f"  Specificity score: {result.specificity_score:.2f}")
        report.append("")

        if result.target_gaps:
            report.append("LARGEST GAPS (Under-covered regions):")
            sorted_gaps = sorted(result.target_gaps, key=lambda g: g['length'], reverse=True)
            for i, gap in enumerate(sorted_gaps[:5], 1):
                report.append(f"  {i}. Position {gap['start']:,}-{gap['end']:,} "
                            f"({gap['length']:,} bp, coverage={gap['mean_coverage']:.2f})")
            report.append("")

        report.append("=" * 80)

        report_text = "\n".join(report)

        if output_path:
            with open(output_path, 'w') as f:
                f.write(report_text)
            logger.info(f"Report saved to: {output_path}")
        else:
            print(report_text)
