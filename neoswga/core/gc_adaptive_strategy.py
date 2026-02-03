#!/usr/bin/env python3
"""
GC-adaptive primer selection strategy for SWGA.

Automatically adjusts optimization parameters based on genome GC content:
- AT-rich genomes (<35% GC): Shorter primers, lower Tm, less betaine
- Balanced genomes (35-65% GC): Standard parameters
- GC-rich genomes (>65% GC): Longer primers, higher Tm, more betaine

Also recommends optimal polymerase (Phi29 vs EquiPhi29) based on genome
characteristics.

Author: NeoSWGA Development Team
Date: November 2025
Version: 3.0 - Phase 3.2
"""

import logging
import numpy as np
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
from enum import Enum

logger = logging.getLogger(__name__)


class GenomeClass(Enum):
    """Genome classification by GC content"""
    AT_RICH = "at_rich"      # <35% GC
    BALANCED = "balanced"     # 35-65% GC
    GC_RICH = "gc_rich"       # >65% GC


@dataclass
class GCAdaptiveParameters:
    """
    GC-content-adaptive parameters for primer selection.

    All parameters automatically adjusted based on genome GC content.
    """
    # Genome characteristics
    genome_class: GenomeClass
    gc_content: float

    # Recommended polymerase
    recommended_polymerase: str  # 'phi29' or 'equiphi29'

    # K-mer parameters
    kmer_range: Tuple[int, int]
    optimal_kmer: int

    # Temperature parameters
    reaction_temp: float
    min_primer_tm: float
    max_primer_tm: float

    # GC content constraints
    min_gc: float
    max_gc: float

    # Chemical additives
    betaine_concentration: float  # M
    dmso_concentration: float     # % v/v

    # Thermodynamic thresholds
    max_homodimer_dg: float    # kcal/mol
    max_heterodimer_dg: float  # kcal/mol
    max_hairpin_dg: float      # kcal/mol

    # Primer count adjustment
    primer_count_multiplier: float

    # Extension parameters
    max_extension: int  # bp
    extension_rate: float  # bp/s

    # Confidence score (0-1)
    confidence: float

    def __str__(self):
        return f"""GC-Adaptive Parameters:
  Genome: {self.gc_content:.1%} GC ({self.genome_class.value})
  Polymerase: {self.recommended_polymerase.upper()}
  Temperature: {self.reaction_temp}°C
  K-mer range: {self.kmer_range[0]}-{self.kmer_range[1]} bp
  Betaine: {self.betaine_concentration}M
  Confidence: {self.confidence:.2f}"""


class GCAdaptiveStrategy:
    """
    Automatically selects and configures optimization strategy based on
    genome GC content.

    Key decisions:
    1. Genome classification (AT-rich/Balanced/GC-rich)
    2. Polymerase selection (Phi29 vs EquiPhi29)
    3. Parameter optimization (k-mer, Tm, betaine, etc.)
    4. Confidence scoring

    Usage:
        strategy = GCAdaptiveStrategy(genome_sequence)
        params = strategy.get_parameters()
        print(f"Use {params.recommended_polymerase} at {params.reaction_temp}°C")
    """

    def __init__(self, genome_sequence: Optional[str] = None,
                 genome_gc_content: Optional[float] = None):
        """
        Initialize GC-adaptive strategy.

        Args:
            genome_sequence: Full genome sequence (optional)
            genome_gc_content: Pre-calculated GC content 0-1 (optional)

        Note: Provide either genome_sequence OR genome_gc_content
        """
        if genome_sequence is not None:
            self.gc_content = self._calculate_gc(genome_sequence)
            self.genome_length = len(genome_sequence)
        elif genome_gc_content is not None:
            self.gc_content = genome_gc_content
            self.genome_length = None
        else:
            raise ValueError("Must provide either genome_sequence or genome_gc_content")

        self.genome_class = self._classify_genome()

        logger.info(f"GC-adaptive strategy initialized")
        logger.info(f"  GC content: {self.gc_content:.1%}")
        logger.info(f"  Classification: {self.genome_class.value}")

    def _calculate_gc(self, sequence: str) -> float:
        """Calculate GC content from sequence"""
        upper_seq = sequence.upper()
        gc_count = upper_seq.count('G') + upper_seq.count('C')
        total = len(upper_seq)
        return gc_count / total if total > 0 else 0.5

    def _classify_genome(self) -> GenomeClass:
        """Classify genome by GC content"""
        if self.gc_content < 0.35:
            return GenomeClass.AT_RICH
        elif self.gc_content > 0.65:
            return GenomeClass.GC_RICH
        else:
            return GenomeClass.BALANCED

    def get_parameters(self, preferred_polymerase: Optional[str] = None) -> GCAdaptiveParameters:
        """
        Get GC-adaptive parameters for primer selection.

        Args:
            preferred_polymerase: Override polymerase selection ('phi29' or 'equiphi29')

        Returns:
            GCAdaptiveParameters configured for this genome
        """
        # Select polymerase
        if preferred_polymerase:
            polymerase = preferred_polymerase.lower()
        else:
            polymerase = self._select_polymerase()

        # Get base parameters for polymerase
        if polymerase == 'equiphi29':
            params = self._get_equiphi29_parameters()
        else:
            params = self._get_phi29_parameters()

        # Adjust for GC content
        params = self._adjust_for_gc(params, polymerase)

        return params

    def _select_polymerase(self) -> str:
        """
        Select optimal polymerase based on genome characteristics.

        Decision tree:
        - GC-rich (>65%): EquiPhi29 (better GC tolerance)
        - AT-rich (<35%): Phi29 (lower temp, simpler)
        - Balanced (35-65%): EquiPhi29 (slightly better overall)
        """
        if self.genome_class == GenomeClass.GC_RICH:
            logger.info("  Selected EquiPhi29 (GC-rich genome)")
            return 'equiphi29'
        elif self.genome_class == GenomeClass.AT_RICH:
            logger.info("  Selected Phi29 (AT-rich genome)")
            return 'phi29'
        else:
            # Balanced: EquiPhi29 slightly preferred for overall performance
            logger.info("  Selected EquiPhi29 (balanced genome)")
            return 'equiphi29'

    def _get_phi29_parameters(self) -> GCAdaptiveParameters:
        """Get base Phi29 parameters"""
        return GCAdaptiveParameters(
            genome_class=self.genome_class,
            gc_content=self.gc_content,
            recommended_polymerase='phi29',
            kmer_range=(8, 12),
            optimal_kmer=10,
            reaction_temp=30.0,
            min_primer_tm=25.0,
            max_primer_tm=50.0,
            min_gc=0.30,
            max_gc=0.70,
            betaine_concentration=0.5,
            dmso_concentration=0.0,
            max_homodimer_dg=-9.0,
            max_heterodimer_dg=-9.0,
            max_hairpin_dg=-2.0,
            primer_count_multiplier=1.0,
            max_extension=70000,
            extension_rate=167.0,
            confidence=0.85
        )

    def _get_equiphi29_parameters(self) -> GCAdaptiveParameters:
        """Get base EquiPhi29 parameters"""
        return GCAdaptiveParameters(
            genome_class=self.genome_class,
            gc_content=self.gc_content,
            recommended_polymerase='equiphi29',
            kmer_range=(11, 15),
            optimal_kmer=12,
            reaction_temp=42.0,
            min_primer_tm=37.0,
            max_primer_tm=62.0,
            min_gc=0.25,
            max_gc=0.75,
            betaine_concentration=1.5,
            dmso_concentration=0.0,
            max_homodimer_dg=-10.0,
            max_heterodimer_dg=-10.0,
            max_hairpin_dg=-3.0,
            primer_count_multiplier=0.85,
            max_extension=80000,
            extension_rate=200.0,
            confidence=0.90
        )

    def _adjust_for_gc(self, params: GCAdaptiveParameters, polymerase: str) -> GCAdaptiveParameters:
        """Adjust parameters based on GC content"""

        if self.genome_class == GenomeClass.GC_RICH:
            # GC-rich genome adjustments
            logger.info("  Applying GC-rich adjustments:")

            # Use maximum betaine
            if polymerase == 'equiphi29':
                params.betaine_concentration = 2.0
                params.max_gc = 0.80
                params.kmer_range = (12, 15)
                params.optimal_kmer = 13
            else:
                params.betaine_concentration = 1.5
                params.max_gc = 0.75
                params.kmer_range = (10, 13)
                params.optimal_kmer = 11

            # Add DMSO for extreme GC
            if self.gc_content > 0.70:
                params.dmso_concentration = 5.0  # 5% v/v
                logger.info(f"    Added DMSO: {params.dmso_concentration}% v/v")

            # Higher confidence for EquiPhi29 on GC-rich
            if polymerase == 'equiphi29':
                params.confidence = 0.95

            logger.info(f"    Betaine: {params.betaine_concentration}M")
            logger.info(f"    Max GC: {params.max_gc:.0%}")
            logger.info(f"    K-mer: {params.kmer_range}")

        elif self.genome_class == GenomeClass.AT_RICH:
            # AT-rich genome adjustments
            logger.info("  Applying AT-rich adjustments:")

            # Use less betaine
            if polymerase == 'equiphi29':
                params.betaine_concentration = 0.5
                params.min_gc = 0.20
                params.kmer_range = (10, 13)
                params.optimal_kmer = 11
            else:
                params.betaine_concentration = 0.0
                params.min_gc = 0.25
                params.kmer_range = (8, 11)
                params.optimal_kmer = 9

            # Lower Tm for AT-rich
            if self.gc_content < 0.25:
                params.min_primer_tm = params.reaction_temp - 8  # Wider margin
                logger.info(f"    Lowered min Tm: {params.min_primer_tm}°C")

            # Higher confidence for Phi29 on AT-rich
            if polymerase == 'phi29':
                params.confidence = 0.92

            logger.info(f"    Betaine: {params.betaine_concentration}M")
            logger.info(f"    Min GC: {params.min_gc:.0%}")
            logger.info(f"    K-mer: {params.kmer_range}")

        else:
            # Balanced genome - use standard parameters
            logger.info("  Using standard parameters (balanced genome)")
            if polymerase == 'equiphi29':
                params.betaine_concentration = 1.0

        return params

    def get_strategy_report(self) -> str:
        """Generate detailed strategy report"""
        params = self.get_parameters()

        report = f"""
{'='*80}
GC-ADAPTIVE STRATEGY REPORT
{'='*80}

Genome Characteristics:
  GC content: {self.gc_content:.1%}
  Classification: {self.genome_class.value.upper()}
  {"(AT-rich: challenging for amplification)" if self.genome_class == GenomeClass.AT_RICH else ""}
  {"(GC-rich: requires robust polymerase)" if self.genome_class == GenomeClass.GC_RICH else ""}
  {"(Balanced: standard SWGA conditions)" if self.genome_class == GenomeClass.BALANCED else ""}

Recommended Configuration:
  Polymerase: {params.recommended_polymerase.upper()}
  Temperature: {params.reaction_temp}°C
  K-mer range: {params.kmer_range[0]}-{params.kmer_range[1]} bp (optimal: {params.optimal_kmer})
  Tm range: {params.min_primer_tm:.0f}-{params.max_primer_tm:.0f}°C
  GC range: {params.min_gc:.0%}-{params.max_gc:.0%}

Chemical Additives:
  Betaine: {params.betaine_concentration}M {"(max)" if params.betaine_concentration >= 2.0 else ""}
  DMSO: {params.dmso_concentration}% v/v

Performance Expectations:
  Extension rate: {params.extension_rate} bp/s
  Max processivity: {params.max_extension:,} bp
  Primer count multiplier: {params.primer_count_multiplier}
  Confidence: {params.confidence:.0%}

Rationale:
"""
        if self.genome_class == GenomeClass.GC_RICH:
            report += """  GC-rich genomes benefit from:
  - Higher temperature (denatures GC-rich regions)
  - Betaine (equalizes AT/GC melting)
  - EquiPhi29 (better GC tolerance)
  - Longer k-mers (more stable binding)
"""
        elif self.genome_class == GenomeClass.AT_RICH:
            report += """  AT-rich genomes benefit from:
  - Lower temperature (prevents primer dimer)
  - Minimal betaine (AT binding already favorable)
  - Shorter k-mers (adequate stability)
  - Phi29 (simpler, cost-effective)
"""
        else:
            report += """  Balanced genomes:
  - Standard SWGA conditions work well
  - EquiPhi29 provides slight performance edge
  - Wide parameter tolerance
"""

        report += f"\n{'='*80}\n"
        return report


if __name__ == '__main__':
    # Example usage
    logging.basicConfig(level=logging.INFO)

    print("GC-Adaptive Strategy - Example Usage\n")
    print("Automatically optimizes parameters for any genome GC content\n")

    # Example 1: GC-rich genome (Plasmodium falciparum, 19% GC)
    print("\n" + "="*80)
    print("Example 1: AT-rich genome (Plasmodium, 19% GC)")
    print("="*80)
    strategy1 = GCAdaptiveStrategy(genome_gc_content=0.19)
    params1 = strategy1.get_parameters()
    print(params1)

    # Example 2: Balanced genome (E. coli, 50% GC)
    print("\n" + "="*80)
    print("Example 2: Balanced genome (E. coli, 50% GC)")
    print("="*80)
    strategy2 = GCAdaptiveStrategy(genome_gc_content=0.50)
    params2 = strategy2.get_parameters()
    print(params2)

    # Example 3: GC-rich genome (Streptomyces, 72% GC)
    print("\n" + "="*80)
    print("Example 3: GC-rich genome (Streptomyces, 72% GC)")
    print("="*80)
    strategy3 = GCAdaptiveStrategy(genome_gc_content=0.72)
    params3 = strategy3.get_parameters()
    print(params3)

    # Generate full report
    print("\n" + strategy3.get_strategy_report())
