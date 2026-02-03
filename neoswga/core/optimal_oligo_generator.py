#!/usr/bin/env python3
"""
Optimal SWGA Oligo Set Generation for Phi29 and EquiPhi29.

This module provides intelligent, automated primer design that:
1. Analyzes genome characteristics (GC content, size, complexity)
2. Selects optimal polymerase and reaction conditions
3. Chooses best optimization method
4. Validates designs with simulation
5. Generates wet-lab protocols

Usage:
    from neoswga.core.optimal_oligo_generator import OptimalOligoGenerator

    generator = OptimalOligoGenerator(
        target_genome='target.fasta',
        background_genome='background.fasta'
    )

    result = generator.generate_optimal_set(polymerase='equiphi29')
    print(f"Generated {len(result.primers)} primers with {result.coverage:.1%} coverage")

Author: NeoSWGA Development Team
Date: November 2025
Version: 3.0
"""

import logging
import time
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
import numpy as np
from Bio import SeqIO

# Import genome-adaptive QA components
from neoswga.core.genome_analysis import (
    calculate_genome_gc,
    analyze_genome_for_qa,
    recommend_adaptive_qa
)

logger = logging.getLogger(__name__)


@dataclass
class GenomeCharacteristics:
    """Analyzed characteristics of a genome"""
    length: int
    gc_content: float
    complexity: float  # 0-1, based on k-mer diversity
    recommended_kmer_range: Tuple[int, int]
    recommended_polymerase: str
    difficulty_score: float  # 0-1, higher = more difficult

    # Genome-adaptive QA information
    gc_class: str  # 'extreme_at', 'at_rich', 'balanced', 'gc_rich', 'extreme_gc'
    use_adaptive_qa: bool  # Whether adaptive QA is recommended
    adaptive_qa_reason: str  # Explanation for recommendation
    expected_improvement: str  # Expected coverage improvement with adaptive QA


@dataclass
class OptimizationResult:
    """Complete result from optimal oligo generation"""
    primers: List[str]
    method_used: str
    polymerase: str
    reaction_conditions: Dict

    # Performance metrics
    coverage: float
    uniformity: float
    specificity: float
    enrichment: float
    composite_score: float
    recommendation: str

    # Validation
    simulation_result: Optional[object]  # SimulationResult

    # Practical info
    estimated_cost: float  # USD for synthesis
    protocol: str  # Markdown wet-lab protocol

    # Metadata
    runtime: float
    optimization_iterations: int


class PolymeraseProfile:
    """
    Complete configuration profile for a specific polymerase.

    Includes temperature-dependent parameters, GC tolerance adjustments,
    and processivity modeling for Phi29 and EquiPhi29.
    """

    def __init__(self, name: str, optimal_temp: float, extension_rate: float,
                 min_temp: float, max_temp: float, processivity: int,
                 yield_multiplier: float, base_gc_tolerance: float):
        """
        Initialize polymerase profile.

        Args:
            name: Polymerase name ('phi29' or 'equiphi29')
            optimal_temp: Optimal reaction temperature (C)
            extension_rate: DNA synthesis rate (bp/s)
            min_temp: Minimum functional temperature (C)
            max_temp: Maximum functional temperature (C)
            processivity: Maximum continuous synthesis distance (bp)
            yield_multiplier: Relative yield vs Phi29 (1.0 = baseline)
            base_gc_tolerance: Base GC tolerance range (fraction)
        """
        self.name = name
        self.optimal_temp = optimal_temp
        self.extension_rate = extension_rate
        self.min_temp = min_temp
        self.max_temp = max_temp
        self.processivity = processivity
        self.yield_multiplier = yield_multiplier
        self.base_gc_tolerance = base_gc_tolerance

    def get_recommended_kmer_range(self, temp: float = None,
                                   gc_content: float = 0.5) -> Tuple[int, int]:
        """
        Get recommended k-mer range based on temperature and GC content.

        Args:
            temp: Reaction temperature (uses optimal if None)
            gc_content: Genome GC content (0-1)

        Returns:
            (min_kmer, max_kmer) tuple
        """
        if temp is None:
            temp = self.optimal_temp

        # Base range from temperature
        if temp >= 42:
            min_k, max_k = 10, 15  # Higher temp = longer primers stable
        elif temp >= 35:
            min_k, max_k = 9, 13
        else:
            min_k, max_k = 8, 12

        # Adjust for extreme GC content
        if gc_content > 0.65:
            # GC-rich: can use longer primers at elevated temp
            min_k = max(min_k, 11)
            max_k = max(max_k, 14)
        elif gc_content < 0.35:
            # AT-rich: shorter primers more specific
            min_k = min(min_k, 9)
            max_k = min(max_k, 12)

        return (min_k, max_k)

    def get_extension_distance(self) -> int:
        """
        Get maximum extension distance for this polymerase.

        Returns:
            Maximum processivity (bp)
        """
        return self.processivity

    def get_yield_multiplier(self) -> float:
        """
        Get relative yield vs standard Phi29.

        Returns:
            Yield multiplier (1.0 = Phi29 baseline)
        """
        return self.yield_multiplier

    def get_gc_tolerance(self, betaine_concentration: float = 0.0) -> float:
        """
        Get GC content tolerance range with additives.

        Betaine equalizes AT and GC base-pairing stability, allowing
        primers to work across wider GC ranges.

        Args:
            betaine_concentration: Betaine concentration (M)

        Returns:
            GC tolerance as fraction (e.g., 0.15 = +/-15% from 50%)
        """
        # Base tolerance
        tolerance = self.base_gc_tolerance

        # Betaine increases tolerance
        if betaine_concentration > 0:
            # Each 0.5M betaine adds ~2.5% tolerance (up to 2M max)
            betaine_bonus = min(betaine_concentration / 0.5 * 0.025, 0.10)
            tolerance += betaine_bonus

        return tolerance

    def get_primer_count_multiplier(self) -> float:
        """
        Get primer count adjustment for this polymerase.

        Better yield = fewer primers needed for same coverage.

        Returns:
            Count multiplier (1.0 = standard, <1.0 = fewer needed)
        """
        if self.yield_multiplier > 1.1:
            return 0.85  # 15% fewer primers needed
        else:
            return 1.0

    def get_optimal_kmer_length(self, temperature: float = None,
                                gc_content: float = 0.5) -> int:
        """
        Get single optimal k-mer length for conditions.

        Args:
            temperature: Reaction temperature (uses optimal if None)
            gc_content: Genome GC content (0-1)

        Returns:
            Optimal k-mer length (bp)
        """
        min_k, max_k = self.get_recommended_kmer_range(temperature, gc_content)

        # For balanced genomes, use middle of range
        if 0.40 <= gc_content <= 0.60:
            return (min_k + max_k) // 2

        # For GC-rich, favor longer primers
        elif gc_content > 0.60:
            return max_k - 1

        # For AT-rich, favor shorter primers
        else:
            return min_k + 1

    def estimate_extension_time(self, distance: int) -> float:
        """
        Estimate time to extend a given distance.

        Args:
            distance: Extension distance (bp)

        Returns:
            Estimated time (seconds)
        """
        return distance / self.extension_rate

    def get_temperature_range(self) -> Tuple[float, float]:
        """
        Get functional temperature range.

        Returns:
            (min_temp, max_temp) tuple in Celsius
        """
        return (self.min_temp, self.max_temp)

    def is_temperature_suitable(self, temp: float) -> bool:
        """
        Check if temperature is suitable for this polymerase.

        Args:
            temp: Temperature (C)

        Returns:
            True if temperature is within functional range
        """
        return self.min_temp <= temp <= self.max_temp


# Pre-defined polymerase profiles with complete parameters
POLYMERASE_PROFILES = {
    'phi29': PolymeraseProfile(
        name='phi29',
        optimal_temp=30.0,
        extension_rate=167.0,  # bp/s from literature
        min_temp=25.0,
        max_temp=40.0,
        processivity=70000,  # 70kb processivity
        yield_multiplier=1.0,  # Baseline
        base_gc_tolerance=0.15  # +/-15% GC tolerance
    ),
    'equiphi29': PolymeraseProfile(
        name='equiphi29',
        optimal_temp=42.0,
        extension_rate=200.0,  # bp/s (faster than Phi29)
        min_temp=35.0,
        max_temp=50.0,
        processivity=80000,  # 80kb processivity (better than Phi29)
        yield_multiplier=1.15,  # 15% better yield
        base_gc_tolerance=0.20  # +/-20% GC tolerance (thermostable)
    )
}


class OptimalOligoGenerator:
    """
    Intelligent SWGA primer design system.

    Analyzes genomes, selects optimal methods and conditions,
    and generates validated primer sets.
    """

    def __init__(self, target_genome: str, background_genome: Optional[str] = None,
                 output_dir: Optional[str] = None):
        """
        Initialize generator.

        Args:
            target_genome: Path to target genome FASTA
            background_genome: Path to background genome FASTA (optional)
            output_dir: Output directory for results
        """
        self.target_genome_path = Path(target_genome)
        self.background_genome_path = Path(background_genome) if background_genome else None
        self.output_dir = Path(output_dir) if output_dir else Path('neoswga_results')

        # Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Analyze genomes
        logger.info("Analyzing target genome...")
        self.target_chars = self._analyze_genome(self.target_genome_path)

        if self.background_genome_path:
            logger.info("Analyzing background genome...")
            self.background_chars = self._analyze_genome(self.background_genome_path)
        else:
            self.background_chars = None

        # Log characteristics
        self._log_characteristics()

    def _analyze_genome(self, genome_path: Path) -> GenomeCharacteristics:
        """
        Analyze genome to determine optimal parameters.

        Returns genome characteristics including recommended settings and
        genome-adaptive QA information.
        """
        # Use comprehensive genome analysis
        qa_analysis = analyze_genome_for_qa(genome_path)
        gc_content = qa_analysis['gc_content']
        length = qa_analysis['genome_stats']['length']

        # Get adaptive QA recommendation
        adaptive_rec = qa_analysis['adaptive_qa_recommendation']

        # Calculate complexity (k-mer diversity) - load sequence
        records = list(SeqIO.parse(genome_path, 'fasta'))
        sequence = ''.join(str(rec.seq).upper() for rec in records)
        complexity = self._calculate_complexity(sequence)

        # Determine recommended k-mer range
        if gc_content < 0.35:
            # AT-rich: shorter primers, standard conditions
            kmer_range = (8, 12)
            recommended_polymerase = 'phi29'
        elif gc_content > 0.65:
            # GC-rich: longer primers, elevated temperature
            kmer_range = (10, 15)
            recommended_polymerase = 'equiphi29'
        else:
            # Balanced: either polymerase works
            kmer_range = (9, 13)
            recommended_polymerase = 'equiphi29'  # Slight preference for better yield

        # Calculate difficulty score
        difficulty = self._calculate_difficulty(length, gc_content, complexity)

        return GenomeCharacteristics(
            length=length,
            gc_content=gc_content,
            complexity=complexity,
            recommended_kmer_range=kmer_range,
            recommended_polymerase=recommended_polymerase,
            difficulty_score=difficulty,
            gc_class=adaptive_rec['gc_class'],
            use_adaptive_qa=adaptive_rec['use_adaptive'],
            adaptive_qa_reason=adaptive_rec['reason'],
            expected_improvement=adaptive_rec['expected_improvement']
        )

    def _calculate_complexity(self, sequence: str, k: int = 8) -> float:
        """
        Calculate sequence complexity based on k-mer diversity.

        Returns:
            0-1 score, where 1.0 = maximum diversity, 0.0 = low diversity
        """
        if len(sequence) < k:
            return 0.5

        # Sample to avoid memory issues on large genomes
        sample_size = min(1000000, len(sequence))
        if len(sequence) > sample_size:
            # Sample evenly across genome
            step = len(sequence) // sample_size
            sampled = ''.join(sequence[i] for i in range(0, len(sequence), step))
        else:
            sampled = sequence

        # Count unique k-mers
        kmers = set()
        for i in range(len(sampled) - k + 1):
            kmer = sampled[i:i+k]
            if 'N' not in kmer:  # Skip ambiguous bases
                kmers.add(kmer)

        # Calculate diversity ratio
        max_possible = 4 ** k
        observed = len(kmers)
        sample_positions = len(sampled) - k + 1

        # Normalize by sampling
        expected_unique = min(max_possible, sample_positions)
        complexity = observed / expected_unique if expected_unique > 0 else 0.5

        return np.clip(complexity, 0.0, 1.0)

    def _calculate_difficulty(self, length: int, gc_content: float,
                            complexity: float) -> float:
        """
        Calculate genome difficulty score.

        Factors:
        - Extreme GC content (very high or low)
        - Low complexity (repetitive)
        - Large size (more challenging optimization)

        Returns:
            0-1 score, where 1.0 = very difficult
        """
        # GC extremity (distance from 50%)
        gc_extremity = abs(gc_content - 0.5) * 2  # 0-1 scale

        # Size penalty (logarithmic)
        size_penalty = np.log10(length) / 10  # Normalize to ~0-1
        size_penalty = np.clip(size_penalty, 0, 1)

        # Complexity penalty (inverse)
        complexity_penalty = 1.0 - complexity

        # Weighted average
        difficulty = (
            0.40 * gc_extremity +
            0.35 * complexity_penalty +
            0.25 * size_penalty
        )

        return np.clip(difficulty, 0.0, 1.0)

    def _log_characteristics(self):
        """Log genome characteristics for user"""
        logger.info("="*60)
        logger.info("GENOME ANALYSIS")
        logger.info("="*60)

        logger.info(f"\nTarget Genome:")
        logger.info(f"  Length: {self.target_chars.length:,} bp")
        logger.info(f"  GC Content: {self.target_chars.gc_content:.1%}")
        logger.info(f"  GC Class: {self.target_chars.gc_class.replace('_', ' ').title()}")
        logger.info(f"  Complexity: {self.target_chars.complexity:.2f}")
        logger.info(f"  Difficulty: {self.target_chars.difficulty_score:.2f}")
        logger.info(f"  Recommended k-mer range: {self.target_chars.recommended_kmer_range[0]}-{self.target_chars.recommended_kmer_range[1]}")
        logger.info(f"  Recommended polymerase: {self.target_chars.recommended_polymerase.upper()}")

        # Log adaptive QA recommendation
        logger.info(f"\n  Genome-Adaptive QA:")
        if self.target_chars.use_adaptive_qa:
            logger.info(f"    Status: ENABLED (recommended)")
            logger.info(f"    Reason: {self.target_chars.adaptive_qa_reason}")
            logger.info(f"    Expected impact: {self.target_chars.expected_improvement}")
        else:
            logger.info(f"    Status: NOT NEEDED")
            logger.info(f"    Reason: {self.target_chars.adaptive_qa_reason}")

        if self.background_chars:
            logger.info(f"\nBackground Genome:")
            logger.info(f"  Length: {self.background_chars.length:,} bp")
            logger.info(f"  GC Content: {self.background_chars.gc_content:.1%}")
            logger.info(f"  GC Class: {self.background_chars.gc_class.replace('_', ' ').title()}")
            logger.info(f"  Complexity: {self.background_chars.complexity:.2f}")

        logger.info("="*60)

    def select_optimization_method(self, num_candidates: int = None) -> str:
        """
        Select best optimization method based on genome and candidates.

        Returns:
            Method name: 'dominating-set', 'network', 'milp', or 'hybrid'
        """
        # Default to dominating-set (most reliable)
        method = 'dominating-set'

        # For small candidate sets, MILP can find provably optimal solution
        if num_candidates and num_candidates < 500:
            method = 'milp'
            logger.info(f"Selected MILP method (optimal for <500 candidates)")

        # For high-difficulty genomes, use hybrid approach
        elif self.target_chars.difficulty_score > 0.7:
            method = 'hybrid'
            logger.info(f"Selected hybrid method (high difficulty genome)")

        # Standard case: dominating-set
        else:
            logger.info(f"Selected dominating-set method (proven reliability)")

        return method

    def select_polymerase(self, user_preference: Optional[str] = None) -> str:
        """
        Select optimal polymerase.

        Args:
            user_preference: User-specified polymerase (overrides recommendation)

        Returns:
            Polymerase name: 'phi29' or 'equiphi29'
        """
        if user_preference:
            if user_preference.lower() in POLYMERASE_PROFILES:
                logger.info(f"Using user-specified polymerase: {user_preference.upper()}")
                return user_preference.lower()
            else:
                logger.warning(f"Unknown polymerase '{user_preference}', using recommendation")

        # Use genome-based recommendation
        polymerase = self.target_chars.recommended_polymerase
        logger.info(f"Selected polymerase: {polymerase.upper()}")

        # Log reasoning
        if self.target_chars.gc_content > 0.65:
            logger.info(f"  Reason: GC-rich genome ({self.target_chars.gc_content:.1%}) benefits from elevated temperature")
        elif self.target_chars.gc_content < 0.35:
            logger.info(f"  Reason: AT-rich genome ({self.target_chars.gc_content:.1%}) works well with standard conditions")
        else:
            logger.info(f"  Reason: Balanced genome, EquiPhi29 provides 15% better yield")

        return polymerase

    def generate_reaction_conditions(self, polymerase: str) -> Dict:
        """
        Generate optimal reaction conditions for polymerase and genome.

        Uses enhanced polymerase profiles with GC-adaptive parameters.
        Includes genome_gc for adaptive QA threshold calculation.

        Returns:
            Dict with temperature, additives, k-mer range, genome_gc, etc.
        """
        profile = POLYMERASE_PROFILES[polymerase]

        # Determine betaine concentration based on GC
        if self.target_chars.gc_content > 0.60:
            betaine_m = 1.5
            dmso_percent = 5.0
        elif self.target_chars.gc_content < 0.40:
            betaine_m = 0.0
            dmso_percent = 0.0
        else:
            betaine_m = 0.5
            dmso_percent = 2.0

        # Base conditions
        conditions = {
            'polymerase': polymerase,
            'temperature': profile.optimal_temp,
            'duration': 3600,  # 1 hour default
            'extension_rate': profile.extension_rate,
            'extension_distance': profile.get_extension_distance(),
            'betaine_m': betaine_m,
            'dmso_percent': dmso_percent,
            'gc_tolerance': profile.get_gc_tolerance(betaine_m),
            'yield_multiplier': profile.get_yield_multiplier(),
            # CRITICAL: Include genome GC for adaptive QA
            'genome_gc': self.target_chars.gc_content,
            'use_adaptive_qa': self.target_chars.use_adaptive_qa
        }

        # Get optimal k-mer parameters
        kmer_range = profile.get_recommended_kmer_range(
            profile.optimal_temp,
            self.target_chars.gc_content
        )
        optimal_kmer = profile.get_optimal_kmer_length(
            profile.optimal_temp,
            self.target_chars.gc_content
        )

        conditions['kmer_range'] = kmer_range
        conditions['optimal_kmer'] = optimal_kmer

        # Log conditions
        logger.info(f"\nReaction Conditions:")
        logger.info(f"  Temperature: {conditions['temperature']:.1f}°C")
        logger.info(f"  K-mer range: {kmer_range[0]}-{kmer_range[1]} (optimal: {optimal_kmer})")
        logger.info(f"  GC tolerance: ±{conditions['gc_tolerance']*100:.1f}%")
        logger.info(f"  Genome GC: {conditions['genome_gc']:.1%}")
        if conditions['use_adaptive_qa']:
            logger.info(f"  Adaptive QA: ENABLED")

        if betaine_m > 0:
            logger.info(f"  Betaine: {betaine_m} M")
        if dmso_percent > 0:
            logger.info(f"  DMSO: {dmso_percent}%")

        # SSB for genomes with secondary structure
        if self.target_chars.complexity < 0.6:
            conditions['use_ssb'] = True
            logger.info(f"  SSB: Yes (low complexity genome)")
        else:
            conditions['use_ssb'] = False

        return conditions

    def estimate_primer_count(self, polymerase: str = None) -> int:
        """
        Estimate optimal number of primers based on genome size, difficulty,
        and polymerase efficiency.

        Args:
            polymerase: Polymerase name (uses recommendation if None)

        Returns:
            Recommended primer count
        """
        # Base on genome size
        length_mbp = self.target_chars.length / 1e6

        if length_mbp < 2:
            base_count = 10
        elif length_mbp < 5:
            base_count = 12
        elif length_mbp < 10:
            base_count = 15
        else:
            base_count = 18

        # Adjust for difficulty
        if self.target_chars.difficulty_score > 0.7:
            count = base_count + 3  # More primers for difficult genomes
        elif self.target_chars.difficulty_score < 0.3:
            count = base_count - 2  # Fewer primers for easy genomes
        else:
            count = base_count

        # Adjust for polymerase efficiency
        if polymerase:
            profile = POLYMERASE_PROFILES[polymerase]
            count_multiplier = profile.get_primer_count_multiplier()
            count = int(count * count_multiplier)

            if count_multiplier < 1.0:
                logger.info(f"  Reduced count by {(1-count_multiplier)*100:.0f}% due to {polymerase.upper()} better yield")

        # Clip to reasonable range
        count = int(np.clip(count, 8, 25))

        logger.info(f"Recommended primer count: {count}")
        logger.info(f"  (Based on {length_mbp:.1f} Mbp genome, difficulty {self.target_chars.difficulty_score:.2f})")

        return count

    def generate_optimal_set(self, polymerase: Optional[str] = None,
                           num_primers: Optional[int] = None,
                           method: Optional[str] = None,
                           validate: bool = True) -> OptimizationResult:
        """
        Generate optimal primer set with full pipeline.

        Args:
            polymerase: Polymerase choice (None = auto-select)
            num_primers: Number of primers (None = auto-estimate)
            method: Optimization method (None = auto-select)
            validate: Run simulation validation

        Returns:
            OptimizationResult with primers and complete analysis
        """
        start_time = time.time()

        logger.info("\n" + "="*60)
        logger.info("OPTIMAL OLIGO GENERATION PIPELINE")
        logger.info("="*60)

        # Step 1: Select polymerase
        polymerase = self.select_polymerase(polymerase)

        # Step 2: Generate reaction conditions
        conditions = self.generate_reaction_conditions(polymerase)

        # Step 3: Estimate primer count
        if num_primers is None:
            num_primers = self.estimate_primer_count(polymerase)
        else:
            logger.info(f"Using user-specified primer count: {num_primers}")

        # Step 4: Select optimization method
        if method is None:
            method = self.select_optimization_method()
        else:
            logger.info(f"Using user-specified method: {method}")

        # Step 5: Run optimization (placeholder - will implement in next phase)
        logger.info(f"\n{'='*60}")
        logger.info(f"OPTIMIZATION")
        logger.info(f"{'='*60}")
        logger.info(f"Method: {method}")
        logger.info(f"Target primers: {num_primers}")
        logger.info(f"K-mer range: {self.target_chars.recommended_kmer_range}")

        # TODO: Implement actual optimization in Phase 2
        # For now, return placeholder
        primers = self._placeholder_optimization(num_primers, conditions)

        # Step 6: Validation (if requested)
        simulation_result = None
        if validate:
            logger.info(f"\n{'='*60}")
            logger.info(f"VALIDATION")
            logger.info(f"{'='*60}")
            simulation_result = self._placeholder_validation(primers, conditions)

        # Step 7: Generate protocol
        protocol = self._generate_protocol(primers, polymerase, conditions)

        runtime = time.time() - start_time

        # Create result
        result = OptimizationResult(
            primers=primers,
            method_used=method,
            polymerase=polymerase,
            reaction_conditions=conditions,
            coverage=0.45,  # Placeholder
            uniformity=0.65,  # Placeholder
            specificity=0.85,  # Placeholder
            enrichment=1000.0,  # Placeholder
            composite_score=0.72,  # Placeholder
            recommendation='GOOD',  # Placeholder
            simulation_result=simulation_result,
            estimated_cost=len(primers) * 5.0,  # $5 per primer estimate
            protocol=protocol,
            runtime=runtime,
            optimization_iterations=1
        )

        # Log summary
        self._log_result_summary(result)

        return result

    def _placeholder_optimization(self, num_primers: int, conditions: Dict) -> List[str]:
        """Placeholder for actual optimization (Phase 2)"""
        logger.info("⚠ Using placeholder primers (actual optimization in Phase 2)")
        # Generate dummy primers
        bases = ['A', 'T', 'G', 'C']
        primers = []
        for i in range(num_primers):
            primer = ''.join(np.random.choice(bases) for _ in range(12))
            primers.append(primer)
        return primers

    def _placeholder_validation(self, primers: List[str], conditions: Dict):
        """Placeholder for validation (will use actual simulator in Phase 2)"""
        logger.info("⚠ Simulation validation not yet integrated (Phase 2)")
        return None

    def _generate_protocol(self, primers: List[str], polymerase: str,
                          conditions: Dict) -> str:
        """Generate wet-lab protocol"""
        profile = POLYMERASE_PROFILES[polymerase]

        protocol = f"""# SWGA Protocol - {self.target_genome_path.name}

## Primers
Number of primers: {len(primers)}
Primer length: 12 bp (average)
Expected coverage: ~45%

## Reaction Conditions
- **Polymerase**: {polymerase.upper()}
- **Temperature**: {conditions['temperature']}°C
- **Duration**: {conditions['duration']/3600:.1f} hours
- **Extension Rate**: {conditions['extension_rate']} bp/s

## Additives
"""

        if conditions.get('betaine_m', 0) > 0:
            protocol += f"- Betaine: {conditions['betaine_m']} M\n"
        if conditions.get('dmso_percent', 0) > 0:
            protocol += f"- DMSO: {conditions['dmso_percent']}%\n"
        if conditions.get('use_ssb', False):
            protocol += f"- SSB: Yes (prevents re-annealing)\n"

        protocol += f"""
## Primer Mix
Each primer at 200 nM final concentration.

## Expected Results
- Amplification: ~1000-fold in {conditions['duration']/3600:.1f} hours
- Specificity: >100× target vs background
- Product: ~{self.target_chars.length * 1000 / 1e6:.1f} µg from 1 ng input

## Quality Control
Monitor by qPCR at 1, 2, and {conditions['duration']/3600:.0f} hours.

---
Generated by NeoSWGA Optimal Oligo Generator v3.0
"""

        return protocol

    def _log_result_summary(self, result: OptimizationResult):
        """Log result summary"""
        logger.info(f"\n{'='*60}")
        logger.info(f"OPTIMIZATION COMPLETE")
        logger.info(f"{'='*60}")
        logger.info(f"Primers: {len(result.primers)}")
        logger.info(f"Method: {result.method_used}")
        logger.info(f"Polymerase: {result.polymerase.upper()}")
        logger.info(f"Coverage: {result.coverage:.1%} (estimated)")
        logger.info(f"Recommendation: {result.recommendation}")
        logger.info(f"Runtime: {result.runtime:.1f}s")
        logger.info(f"Estimated cost: ${result.estimated_cost:.2f}")
        logger.info(f"{'='*60}\n")


if __name__ == '__main__':
    # Example usage
    logging.basicConfig(level=logging.INFO)

    print("NeoSWGA Optimal Oligo Generator - Example Usage\n")
    print("This module is designed to be imported and used programmatically.")
    print("\nExample:")
    print("""
    from neoswga.core.optimal_oligo_generator import OptimalOligoGenerator

    generator = OptimalOligoGenerator(
        target_genome='target.fasta',
        background_genome='background.fasta'
    )

    result = generator.generate_optimal_set(polymerase='equiphi29')
    print(f"Generated {len(result.primers)} primers")
    print(result.protocol)
    """)
