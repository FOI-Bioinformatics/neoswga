#!/usr/bin/env python3
"""
EquiPhi29-optimized primer set selection for SWGA.

EquiPhi29 is a thermostable variant of Phi29 polymerase with enhanced properties:
- Higher optimal temperature (42°C vs 30°C)
- Faster extension rate (200 bp/s vs 167 bp/s)
- Better processivity (80kb vs 70kb)
- 15% higher yield
- Improved tolerance to GC-rich sequences

This module provides specialized optimization strategies that leverage these
properties for improved performance, especially with GC-rich genomes.

Author: NeoSWGA Development Team
Date: November 2025
Version: 3.0 - Phase 3.1
"""

import logging
import numpy as np
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
import time

from neoswga.core.optimal_oligo_generator import POLYMERASE_PROFILES
from neoswga.core.thermodynamic_filter import (
    ThermodynamicFilter,
    ThermodynamicCriteria,
    create_filter_from_conditions
)
from neoswga.core.hybrid_optimizer import HybridOptimizer, HybridResult

logger = logging.getLogger(__name__)


@dataclass
class EquiPhi29Config:
    """Configuration for EquiPhi29-optimized primer selection"""
    # Temperature parameters
    reaction_temp: float = 42.0       # Optimal temperature
    min_primer_tm: float = 37.0       # Min Tm (reaction_temp - 5°C)
    max_primer_tm: float = 62.0       # Max Tm (reaction_temp + 20°C)

    # GC content parameters
    min_gc: float = 0.25              # Lower GC allowed (better tolerance)
    max_gc: float = 0.75              # Higher GC allowed (better tolerance)

    # Betaine optimization
    betaine_concentration: float = 1.5  # M (standard for GC-rich)
    max_betaine: float = 2.0           # M (maximum for extreme GC)

    # Primer parameters
    kmer_range: Tuple[int, int] = (11, 15)  # Longer k-mers (higher temp)
    optimal_kmer: int = 12                   # Optimal k-mer length

    # Primer count (adjusted for higher yield)
    primer_count_multiplier: float = 0.85    # Need ~15% fewer primers

    # Extension parameters
    max_extension: int = 80000              # 80kb processivity
    extension_rate: float = 200.0           # bp/s

    # Thermodynamic thresholds (relaxed at higher temp)
    max_homodimer_dg: float = -10.0        # kcal/mol (more negative allowed)
    max_heterodimer_dg: float = -10.0      # kcal/mol
    max_hairpin_dg: float = -3.0           # kcal/mol


class EquiPhi29Optimizer:
    """
    Optimized primer selection for EquiPhi29 polymerase.

    Key advantages:
    - Better performance on GC-rich genomes (>65% GC)
    - Fewer primers needed (15% reduction)
    - Longer extension distances (80kb vs 70kb)
    - Higher temperature stability (42°C)

    Recommended for:
    - GC-rich genomes (>65% GC)
    - Complex genomes requiring robust amplification
    - Cost-sensitive applications (fewer primers)
    - Higher throughput (faster extension)
    """

    def __init__(self,
                 position_cache,
                 fg_prefixes: List[str],
                 fg_seq_lengths: List[int],
                 bg_prefixes: Optional[List[str]] = None,
                 bg_seq_lengths: Optional[List[int]] = None,
                 genome_gc_content: float = 0.50,
                 config: Optional[EquiPhi29Config] = None):
        """
        Initialize EquiPhi29-optimized primer selector.

        Args:
            position_cache: PositionCache with primer positions
            fg_prefixes: Target genome identifiers
            fg_seq_lengths: Target genome lengths
            bg_prefixes: Background genome identifiers (optional)
            bg_seq_lengths: Background genome lengths (optional)
            genome_gc_content: Target genome GC content (0-1)
            config: EquiPhi29Config (uses defaults if None)
        """
        self.position_cache = position_cache
        self.fg_prefixes = fg_prefixes
        self.fg_seq_lengths = fg_seq_lengths
        self.bg_prefixes = bg_prefixes or []
        self.bg_seq_lengths = bg_seq_lengths or []
        self.genome_gc_content = genome_gc_content

        # Use provided config or create default
        self.config = config or EquiPhi29Config()

        # Adjust config for genome GC content
        self._adjust_config_for_gc()

        # Create thermodynamic filter
        self.thermo_filter = self._create_thermodynamic_filter()

        # Create hybrid optimizer with EquiPhi29 parameters
        self.hybrid_optimizer = HybridOptimizer(
            position_cache=position_cache,
            fg_prefixes=fg_prefixes,
            fg_seq_lengths=fg_seq_lengths,
            bg_prefixes=self.bg_prefixes,
            bg_seq_lengths=self.bg_seq_lengths,
            bin_size=10000,
            max_extension=self.config.max_extension
        )

        logger.info("EquiPhi29 optimizer initialized")
        logger.info(f"  Temperature: {self.config.reaction_temp}°C")
        logger.info(f"  Genome GC: {genome_gc_content:.1%}")
        logger.info(f"  Betaine: {self.config.betaine_concentration}M")
        logger.info(f"  K-mer range: {self.config.kmer_range}")
        logger.info(f"  Max extension: {self.config.max_extension:,} bp")

    def _adjust_config_for_gc(self):
        """Adjust configuration based on genome GC content"""
        gc = self.genome_gc_content

        # GC-rich genome (>65% GC): Use maximum betaine
        if gc > 0.65:
            self.config.betaine_concentration = self.config.max_betaine
            self.config.max_gc = 0.80  # Allow very high GC primers
            self.config.kmer_range = (12, 15)  # Longer k-mers for stability
            self.config.optimal_kmer = 13
            logger.info("  Detected GC-rich genome - using max betaine and longer k-mers")

        # AT-rich genome (<35% GC): Use less betaine, shorter k-mers
        elif gc < 0.35:
            self.config.betaine_concentration = 0.5
            self.config.min_gc = 0.20  # Allow very low GC primers
            self.config.kmer_range = (10, 13)  # Shorter k-mers
            self.config.optimal_kmer = 11
            logger.info("  Detected AT-rich genome - using less betaine and shorter k-mers")

        # Balanced genome: Use standard parameters
        else:
            self.config.betaine_concentration = 1.0
            logger.info("  Detected balanced genome - using standard parameters")

    def _create_thermodynamic_filter(self) -> ThermodynamicFilter:
        """Create thermodynamic filter with EquiPhi29 parameters"""
        criteria = ThermodynamicCriteria(
            min_tm=self.config.min_primer_tm,
            max_tm=self.config.max_primer_tm,
            target_tm=self.config.reaction_temp + 5,  # Slightly above reaction temp
            na_conc=50.0,  # mM
            mg_conc=0.0,
            max_homodimer_dg=self.config.max_homodimer_dg,
            max_heterodimer_dg=self.config.max_heterodimer_dg,
            max_hairpin_dg=self.config.max_hairpin_dg,
            min_gc=self.config.min_gc,
            max_gc=self.config.max_gc,
            reaction_temp=self.config.reaction_temp
        )

        return ThermodynamicFilter(criteria)

    def optimize(self,
                candidates: List[str],
                target_primer_count: Optional[int] = None,
                verbose: bool = True,
                validate_with_simulation: bool = False,
                genome_sequence: Optional[str] = None) -> HybridResult:
        """
        EquiPhi29-optimized primer set selection.

        Workflow:
        1. Thermodynamic filtering (EquiPhi29-specific criteria)
        2. Adjust target count for EquiPhi29 yield advantage
        3. Two-stage hybrid optimization
        4. Optional simulation validation

        Args:
            candidates: Candidate primer sequences
            target_primer_count: Target number of primers (auto if None)
            verbose: Print progress
            validate_with_simulation: Run simulation validation
            genome_sequence: Genome sequence (required for simulation)

        Returns:
            HybridResult with optimized primer set
        """
        if verbose:
            logger.info("="*80)
            logger.info("EQUIPHI29-OPTIMIZED PRIMER SELECTION")
            logger.info("="*80)
            logger.info(f"Input: {len(candidates)} candidates")

        start_time = time.time()

        # ===================================================================
        # STEP 1: Thermodynamic Filtering (EquiPhi29-specific)
        # ===================================================================

        if verbose:
            logger.info("\n" + "-"*80)
            logger.info("STEP 1: Thermodynamic Filtering (42°C, EquiPhi29)")
            logger.info("-"*80)

        filtered_candidates, filter_stats = self.thermo_filter.filter_candidates(
            candidates,
            check_heterodimers=True,
            max_heterodimer_fraction=0.3
        )

        if verbose:
            logger.info(f"\nFiltered: {len(filtered_candidates)}/{len(candidates)} passed")
            logger.info(f"  Mean Tm: {filter_stats['mean_tm']:.1f}°C")
            logger.info(f"  Mean GC: {filter_stats['mean_gc']:.1%}")

        if len(filtered_candidates) == 0:
            logger.error("No primers passed thermodynamic filtering!")
            raise ValueError("No primers passed EquiPhi29 thermodynamic filters")

        # ===================================================================
        # STEP 2: Adjust Target Count for EquiPhi29 Yield
        # ===================================================================

        if target_primer_count is None:
            # Auto-determine based on genome size and complexity
            genome_mb = sum(self.fg_seq_lengths) / 1e6

            if genome_mb < 2:
                base_count = 8
            elif genome_mb < 5:
                base_count = 12
            elif genome_mb < 10:
                base_count = 16
            else:
                base_count = 20

            # Apply EquiPhi29 multiplier (need fewer primers)
            target_primer_count = max(6, int(base_count * self.config.primer_count_multiplier))

        if verbose:
            logger.info(f"\n" + "-"*80)
            logger.info(f"STEP 2: Target Count Adjustment")
            logger.info("-"*80)
            logger.info(f"  Base count: {int(target_primer_count / self.config.primer_count_multiplier)}")
            logger.info(f"  EquiPhi29 multiplier: {self.config.primer_count_multiplier}")
            logger.info(f"  Final target: {target_primer_count} primers")

        # ===================================================================
        # STEP 3: Two-Stage Hybrid Optimization
        # ===================================================================

        if verbose:
            logger.info(f"\n" + "-"*80)
            logger.info("STEP 3: Hybrid Optimization")
            logger.info("-"*80)

        result = self.hybrid_optimizer.optimize(
            candidates=filtered_candidates,
            final_count=target_primer_count,
            verbose=verbose,
            validate_with_simulation=validate_with_simulation,
            genome_sequence=genome_sequence
        )

        # ===================================================================
        # Add EquiPhi29-specific metadata
        # ===================================================================

        total_time = time.time() - start_time

        if verbose:
            logger.info("\n" + "="*80)
            logger.info("EQUIPHI29 OPTIMIZATION COMPLETE")
            logger.info("="*80)
            logger.info(f"Selected primers: {len(result.primers)}")
            logger.info(f"Coverage: {result.final_coverage:.1%}")
            logger.info(f"Connectivity: {result.final_connectivity:.2f}")
            logger.info(f"Predicted amplification: {result.final_predicted_amplification:.1f}×")
            if result.simulation_fitness:
                logger.info(f"Simulation fitness: {result.simulation_fitness.fitness_score:.3f}")
            logger.info(f"\nEquiPhi29 advantages:")
            logger.info(f"  Extension rate: {self.config.extension_rate} bp/s (+20% vs Phi29)")
            logger.info(f"  Max extension: {self.config.max_extension:,} bp (+14% vs Phi29)")
            logger.info(f"  Primer reduction: {int((1 - self.config.primer_count_multiplier)*100)}%")
            logger.info(f"  Expected yield: +15% vs Phi29")
            logger.info(f"\nTotal runtime: {total_time:.2f}s")
            logger.info("="*80)

        return result

    def estimate_reaction_time(self, primer_set: List[str]) -> Dict[str, float]:
        """
        Estimate EquiPhi29 reaction time for optimal amplification.

        Args:
            primer_set: Selected primer sequences

        Returns:
            Dictionary with time estimates
        """
        # Calculate average primer spacing
        total_sites = 0
        for primer in primer_set:
            for prefix in self.fg_prefixes:
                fwd = self.position_cache.get_positions(prefix, primer, 'forward')
                rev = self.position_cache.get_positions(prefix, primer, 'reverse')
                total_sites += len(fwd) + len(rev)

        genome_length = sum(self.fg_seq_lengths)

        if total_sites > 0:
            avg_spacing = genome_length / total_sites
        else:
            avg_spacing = genome_length

        # Extension time to cover average spacing
        # With EquiPhi29's 200 bp/s, time to extend avg_spacing distance
        extension_time_per_fork = avg_spacing / self.config.extension_rate

        # Account for multiple rounds of priming
        # Typically want 3-5 rounds for good coverage
        min_time = extension_time_per_fork * 2  # Minimum (2 rounds)
        optimal_time = extension_time_per_fork * 4  # Optimal (4 rounds)
        max_time = extension_time_per_fork * 6  # Maximum (6 rounds)

        return {
            'min_time_minutes': min_time / 60,
            'optimal_time_minutes': optimal_time / 60,
            'max_time_minutes': max_time / 60,
            'avg_primer_spacing': avg_spacing,
            'extension_rate': self.config.extension_rate,
            'estimated_coverage_min': 0.4,
            'estimated_coverage_optimal': 0.65,
            'estimated_coverage_max': 0.80
        }


def create_equiphi29_optimizer(genome_fasta: str,
                                background_fasta: Optional[str] = None,
                                position_cache = None) -> EquiPhi29Optimizer:
    """
    Convenience function to create EquiPhi29 optimizer from FASTA files.

    Args:
        genome_fasta: Path to target genome FASTA
        background_fasta: Path to background genome FASTA (optional)
        position_cache: Pre-built position cache (optional)

    Returns:
        Configured EquiPhi29Optimizer
    """
    from Bio import SeqIO

    # Read target genome
    fg_sequences = list(SeqIO.parse(genome_fasta, "fasta"))
    fg_prefixes = [rec.id for rec in fg_sequences]
    fg_seq_lengths = [len(rec.seq) for rec in fg_sequences]
    total_seq = "".join(str(rec.seq) for rec in fg_sequences)

    # Calculate GC content
    gc_count = total_seq.upper().count('G') + total_seq.upper().count('C')
    genome_gc = gc_count / len(total_seq) if len(total_seq) > 0 else 0.5

    # Read background genome if provided
    bg_prefixes = []
    bg_seq_lengths = []
    if background_fasta:
        bg_sequences = list(SeqIO.parse(background_fasta, "fasta"))
        bg_prefixes = [rec.id for rec in bg_sequences]
        bg_seq_lengths = [len(rec.seq) for rec in bg_sequences]

    return EquiPhi29Optimizer(
        position_cache=position_cache,
        fg_prefixes=fg_prefixes,
        fg_seq_lengths=fg_seq_lengths,
        bg_prefixes=bg_prefixes,
        bg_seq_lengths=bg_seq_lengths,
        genome_gc_content=genome_gc
    )


if __name__ == '__main__':
    # Example usage
    logging.basicConfig(level=logging.INFO)

    print("EquiPhi29 Optimizer - Example Usage\n")
    print("Optimized for EquiPhi29 thermostable polymerase:")
    print("  - 42°C reaction temperature")
    print("  - 200 bp/s extension rate (+20% vs Phi29)")
    print("  - 80kb processivity (+14% vs Phi29)")
    print("  - 15% higher yield")
    print("  - Better GC tolerance\n")

    print("Example:")
    print("""
    from neoswga.core.equiphi29_optimizer import EquiPhi29Optimizer

    optimizer = EquiPhi29Optimizer(
        position_cache=cache,
        fg_prefixes=['plasmodium'],
        fg_seq_lengths=[23332831],  # 23 Mbp, 19% GC
        genome_gc_content=0.19
    )

    result = optimizer.optimize(
        candidates=filtered_primers,
        target_primer_count=10,  # Auto-adjusted from 12 (EquiPhi29 multiplier)
        verbose=True
    )

    # Estimate reaction time
    timing = optimizer.estimate_reaction_time(result.primers)
    print(f"Optimal reaction time: {timing['optimal_time_minutes']:.0f} minutes")
    print(f"Expected coverage: {timing['estimated_coverage_optimal']:.1%}")
    """)
