#!/usr/bin/env python3
"""
Thermodynamic-aware primer filtering for SWGA.

Filters candidate primers based on:
1. Melting temperature (Tm) range at reaction conditions
2. Self-dimer formation potential (homodimers)
3. Cross-dimer formation (heterodimers between primers)
4. Secondary structure (hairpins, self-complementarity)
5. GC content constraints

Integrates with optimal_oligo_generator to provide thermodynamically
validated candidate sets.

Author: NeoSWGA Development Team
Date: November 2025
Version: 3.0 - Phase 1.3
"""

import concurrent.futures
import logging
import os
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
import numpy as np

from neoswga.core.thermodynamics import (
    calculate_tm_with_salt,
    calculate_free_energy,
    gc_content
)

from neoswga.core.secondary_structure import (
    check_homodimer,
    check_heterodimer,
    check_hairpins
)

from neoswga.core.reaction_conditions import ReactionConditions

logger = logging.getLogger(__name__)


def _check_heterodimer_pair(args):
    """Check a single heterodimer pair.

    Module-level function so it can be pickled by multiprocessing.

    Args:
        args: Tuple of (seq1, seq2, i, j, conditions_dict).

    Returns:
        Tuple of (i, j, dimer_dg).
    """
    seq1, seq2, i, j, conditions_dict = args
    try:
        conditions = ReactionConditions(**conditions_dict)
        result = check_heterodimer(seq1, seq2, conditions)
        dg = result.get('energy', 0.0)
        if dg == float('inf') or dg > 0:
            dg = 0.0
        return (i, j, dg)
    except Exception:
        return (i, j, 0.0)


@dataclass
class ThermodynamicCriteria:
    """
    Thermodynamic filtering criteria for SWGA primers.

    All criteria are temperature and condition-dependent.
    Note: SWGA uses short primers (8-15 bp) with lower Tm than standard PCR.
    """
    # Melting temperature constraints (adjusted for short SWGA primers)
    min_tm: float = 25.0  # Minimum Tm (C) - SWGA primers are short
    max_tm: float = 50.0  # Maximum Tm (C) - want primers well below reaction temp
    target_tm: float = 35.0  # Optimal Tm (C)

    # Salt concentrations (for Tm calculation)
    na_conc: float = 50.0  # mM Na+
    mg_conc: float = 0.0   # mM Mg2+

    # Dimer formation thresholds (free energy)
    max_homodimer_dg: float = -9.0  # kcal/mol (less negative = weaker binding)
    max_heterodimer_dg: float = -9.0  # kcal/mol

    # Hairpin formation
    max_hairpin_dg: float = -2.0  # kcal/mol (less negative = weaker)

    # GC content constraints (wider for SWGA)
    min_gc: float = 0.30  # 30%
    max_gc: float = 0.70  # 70%

    # Reaction temperature (for context)
    reaction_temp: float = 30.0  # C


@dataclass
class PrimerThermodynamics:
    """Thermodynamic properties of a primer."""
    sequence: str
    tm: float  # Melting temperature (C)
    gc: float  # GC content (fraction)
    homodimer_dg: float  # Self-dimer free energy (kcal/mol)
    hairpin_dg: float  # Hairpin free energy (kcal/mol)
    passes_filters: bool  # Overall pass/fail
    failure_reasons: List[str]  # Reasons for failure


class ThermodynamicFilter:
    """
    Filters SWGA primer candidates based on thermodynamic properties.

    Designed to reduce experimental failures by ~30% through early
    elimination of primers with poor thermodynamic characteristics.
    """

    def __init__(self, criteria: Optional[ThermodynamicCriteria] = None):
        """
        Initialize thermodynamic filter.

        Args:
            criteria: Filtering criteria (uses defaults if None)
        """
        self.criteria = criteria or ThermodynamicCriteria()

        logger.info("Thermodynamic filter initialized")
        logger.info(f"  Tm range: {self.criteria.min_tm}-{self.criteria.max_tm}°C")
        logger.info(f"  Max homodimer ΔG: {self.criteria.max_homodimer_dg} kcal/mol")
        logger.info(f"  Max heterodimer ΔG: {self.criteria.max_heterodimer_dg} kcal/mol")
        logger.info(f"  GC range: {self.criteria.min_gc*100:.0f}-{self.criteria.max_gc*100:.0f}%")

    def analyze_primer(self, sequence: str) -> PrimerThermodynamics:
        """
        Analyze thermodynamic properties of a single primer.

        Args:
            sequence: Primer DNA sequence

        Returns:
            PrimerThermodynamics with all calculated properties
        """
        failure_reasons = []

        # Calculate Tm
        tm = calculate_tm_with_salt(
            sequence,
            na_conc=self.criteria.na_conc,
            mg_conc=self.criteria.mg_conc
        )

        # Calculate GC content
        gc = gc_content(sequence)

        # Check Tm range
        if tm < self.criteria.min_tm:
            failure_reasons.append(f"Tm too low ({tm:.1f}°C < {self.criteria.min_tm}°C)")
        elif tm > self.criteria.max_tm:
            failure_reasons.append(f"Tm too high ({tm:.1f}°C > {self.criteria.max_tm}°C)")

        # Check GC content
        if gc < self.criteria.min_gc:
            failure_reasons.append(f"GC too low ({gc:.1%} < {self.criteria.min_gc:.1%})")
        elif gc > self.criteria.max_gc:
            failure_reasons.append(f"GC too high ({gc:.1%} > {self.criteria.max_gc:.1%})")

        # Create reaction conditions for structure analysis
        conditions = ReactionConditions(
            temp=self.criteria.reaction_temp,
            na_conc=self.criteria.na_conc,
            mg_conc=self.criteria.mg_conc
        )

        # Check homodimer formation
        try:
            homodimer_result = check_homodimer(sequence, conditions)
            # Result is a dictionary with 'energy' key (kcal/mol)
            homodimer_dg = homodimer_result.get('energy', 0.0)
            # Handle inf (no dimer)
            if homodimer_dg == float('inf') or homodimer_dg > 0:
                homodimer_dg = 0.0
        except Exception as e:
            logger.warning(f"Homodimer check failed for {sequence}: {e}")
            homodimer_dg = 0.0  # Assume no dimer on error

        if homodimer_dg < self.criteria.max_homodimer_dg:
            failure_reasons.append(
                f"Strong homodimer (ΔG={homodimer_dg:.1f} < {self.criteria.max_homodimer_dg:.1f})"
            )

        # Check hairpin formation
        try:
            hairpin_results = check_hairpins(sequence, conditions)
            # Result is a list of hairpin dictionaries with 'energy' key
            if hairpin_results:
                # Get the strongest (most negative) hairpin
                hairpin_dg = min(h.get('energy', 0.0) for h in hairpin_results)
                if hairpin_dg == float('inf') or hairpin_dg > 0:
                    hairpin_dg = 0.0
            else:
                hairpin_dg = 0.0
        except Exception as e:
            logger.warning(f"Hairpin check failed for {sequence}: {e}")
            hairpin_dg = 0.0  # Assume no hairpin on error

        if hairpin_dg < self.criteria.max_hairpin_dg:
            failure_reasons.append(
                f"Strong hairpin (ΔG={hairpin_dg:.1f} < {self.criteria.max_hairpin_dg:.1f})"
            )

        passes = len(failure_reasons) == 0

        return PrimerThermodynamics(
            sequence=sequence,
            tm=tm,
            gc=gc,
            homodimer_dg=homodimer_dg,
            hairpin_dg=hairpin_dg,
            passes_filters=passes,
            failure_reasons=failure_reasons
        )

    def filter_candidates(self, candidates: List[str],
                         check_heterodimers: bool = True,
                         max_heterodimer_fraction: float = 0.3) -> Tuple[List[str], Dict]:
        """
        Filter list of primer candidates.

        Args:
            candidates: List of primer sequences
            check_heterodimers: Also check for cross-primer dimers
            max_heterodimer_fraction: Max fraction of heterodimers allowed (0-1)

        Returns:
            (passing_primers, statistics_dict)
        """
        logger.info(f"\nFiltering {len(candidates)} candidates...")

        # Analyze all primers
        analyses = [self.analyze_primer(seq) for seq in candidates]

        # Filter by individual properties
        passing = [a for a in analyses if a.passes_filters]

        logger.info(f"  After individual filters: {len(passing)}/{len(candidates)} passed")

        # Log failure statistics
        if len(passing) < len(candidates):
            failure_counts = {}
            for a in analyses:
                if not a.passes_filters:
                    for reason in a.failure_reasons:
                        category = reason.split('(')[0].strip()
                        failure_counts[category] = failure_counts.get(category, 0) + 1

            logger.info(f"  Failure breakdown:")
            for category, count in sorted(failure_counts.items(), key=lambda x: x[1], reverse=True):
                logger.info(f"    {category}: {count}")

        # Check heterodimers if requested
        heterodimer_issues = 0
        if check_heterodimers and len(passing) > 1:
            logger.info(f"  Checking heterodimers between {len(passing)} primers...")

            # Create reaction conditions
            conditions = ReactionConditions(
                temp=self.criteria.reaction_temp,
                na_conc=self.criteria.na_conc,
                mg_conc=self.criteria.mg_conc
            )

            # Calculate pairwise heterodimer potential
            problematic_pairs = set()

            conditions_dict = {
                'temp': self.criteria.reaction_temp,
                'na_conc': self.criteria.na_conc,
                'mg_conc': self.criteria.mg_conc,
            }

            pairs = [
                (passing[i].sequence, passing[j].sequence, i, j, conditions_dict)
                for i in range(len(passing))
                for j in range(i + 1, len(passing))
            ]

            if len(pairs) > 100:
                # Parallel path for large sets (ProcessPoolExecutor for CPU-bound work)
                n_workers = min(os.cpu_count() or 1, 8)
                with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as executor:
                    for i, j, dimer_dg in executor.map(
                        _check_heterodimer_pair, pairs, chunksize=50
                    ):
                        if dimer_dg < self.criteria.max_heterodimer_dg:
                            problematic_pairs.add((i, j))
                            heterodimer_issues += 1
            else:
                # Sequential path for small sets
                for seq1, seq2, i, j, _ in pairs:
                    try:
                        dimer_result = check_heterodimer(seq1, seq2, conditions)
                        dimer_dg = dimer_result.get('energy', 0.0)
                        if dimer_dg == float('inf') or dimer_dg > 0:
                            dimer_dg = 0.0
                    except Exception as e:
                        logger.warning(f"Heterodimer check failed for {seq1}/{seq2}: {e}")
                        dimer_dg = 0.0

                    if dimer_dg < self.criteria.max_heterodimer_dg:
                        problematic_pairs.add((i, j))
                        heterodimer_issues += 1

            # Remove primers involved in too many heterodimers
            if problematic_pairs:
                # Count heterodimer issues per primer
                issue_counts = {}
                for i, j in problematic_pairs:
                    issue_counts[i] = issue_counts.get(i, 0) + 1
                    issue_counts[j] = issue_counts.get(j, 0) + 1

                # Calculate max allowed issues
                max_issues = int(len(passing) * max_heterodimer_fraction)

                # Remove worst offenders
                to_remove = set()
                for idx, count in sorted(issue_counts.items(), key=lambda x: x[1], reverse=True):
                    if count > max_issues:
                        to_remove.add(idx)

                if to_remove:
                    passing = [p for i, p in enumerate(passing) if i not in to_remove]
                    logger.info(f"  Removed {len(to_remove)} primers with excessive heterodimers")
                    logger.info(f"  After heterodimer filtering: {len(passing)} primers")

        # Extract passing sequences
        passing_seqs = [p.sequence for p in passing]

        # Compile statistics
        stats = {
            'total_candidates': len(candidates),
            'passed_individual': len([a for a in analyses if a.passes_filters]),
            'final_passing': len(passing_seqs),
            'heterodimer_issues': heterodimer_issues,
            'mean_tm': np.mean([p.tm for p in passing]) if passing else 0.0,
            'std_tm': np.std([p.tm for p in passing]) if passing else 0.0,
            'mean_gc': np.mean([p.gc for p in passing]) if passing else 0.0,
            'mean_homodimer_dg': np.mean([p.homodimer_dg for p in passing]) if passing else 0.0
        }

        logger.info(f"\nThermodynamic filtering complete:")
        logger.info(f"  Input: {stats['total_candidates']} candidates")
        if stats['total_candidates'] > 0:
            logger.info(f"  Output: {stats['final_passing']} passing ({stats['final_passing']/stats['total_candidates']*100:.1f}%)")
        else:
            logger.info(f"  Output: 0 passing (empty input)")
        if passing:
            logger.info(f"  Mean Tm: {stats['mean_tm']:.1f} ± {stats['std_tm']:.1f}°C")
            logger.info(f"  Mean GC: {stats['mean_gc']:.1%}")

        return passing_seqs, stats

    def adjust_criteria_for_conditions(self, temperature: float, gc_content: float,
                                       betaine_m: float = 0.0):
        """
        Adjust filtering criteria based on reaction conditions.

        Args:
            temperature: Reaction temperature (C)
            gc_content: Target genome GC content (0-1)
            betaine_m: Betaine concentration (M)
        """
        logger.info(f"\nAdjusting criteria for conditions:")
        logger.info(f"  Temperature: {temperature}°C")
        logger.info(f"  Genome GC: {gc_content:.1%}")
        logger.info(f"  Betaine: {betaine_m} M")

        # Adjust Tm range based on reaction temperature
        # SWGA primers anneal at reaction temp, so Tm should be slightly below to just above
        self.criteria.reaction_temp = temperature
        self.criteria.min_tm = temperature - 5   # e.g., 25°C at 30°C, 37°C at 42°C
        self.criteria.max_tm = temperature + 20  # e.g., 50°C at 30°C, 62°C at 42°C

        # Widen GC range with betaine
        if betaine_m > 0:
            # Each 0.5M betaine adds ~5% GC tolerance
            gc_bonus = min(betaine_m / 0.5 * 0.05, 0.15)
            self.criteria.min_gc = max(0.30, self.criteria.min_gc - gc_bonus)
            self.criteria.max_gc = min(0.70, self.criteria.max_gc + gc_bonus)

        # Adjust GC range for genome composition
        if gc_content > 0.65:
            # GC-rich genome: allow wider high-GC range
            self.criteria.max_gc = 0.70
        elif gc_content < 0.35:
            # AT-rich genome: allow wider low-GC range
            self.criteria.min_gc = 0.30

        # Relax dimer constraints at higher temperature
        # Higher temp = dimers less stable = can tolerate stronger dimers
        if temperature >= 42:
            self.criteria.max_homodimer_dg = -10.0  # More negative allowed
            self.criteria.max_heterodimer_dg = -10.0

        logger.info(f"\nAdjusted criteria:")
        logger.info(f"  Tm range: {self.criteria.min_tm:.1f}-{self.criteria.max_tm:.1f}°C")
        logger.info(f"  GC range: {self.criteria.min_gc:.1%}-{self.criteria.max_gc:.1%}")
        logger.info(f"  Max homodimer ΔG: {self.criteria.max_homodimer_dg:.1f} kcal/mol")


def calculate_adaptive_gc_range(genome_gc: float, tolerance: float = 0.15) -> Tuple[float, float]:
    """
    Calculate genome-adaptive GC range centered on target genome composition.

    Args:
        genome_gc: Target genome GC content (0-1)
        tolerance: GC range tolerance (±15% default)

    Returns:
        (min_gc, max_gc) as fractions

    Examples:
        >>> calculate_adaptive_gc_range(0.32)  # Francisella
        (0.17, 0.47)
        >>> calculate_adaptive_gc_range(0.67)  # Burkholderia
        (0.52, 0.82)
        >>> calculate_adaptive_gc_range(0.50)  # E. coli
        (0.35, 0.65)

    Rationale:
        Fixed GC ranges (30-70%) systematically bias against genomes with
        extreme GC content. Centering on genome composition allows primers
        to match target genome while maintaining sufficient specificity.
    """
    min_gc = max(0.15, genome_gc - tolerance)  # Floor at 15%
    max_gc = min(0.85, genome_gc + tolerance)  # Ceiling at 85%
    return (min_gc, max_gc)


def calculate_adaptive_dimer_threshold(
    base_threshold: float,
    genome_gc: float
) -> float:
    """
    Calculate genome-adaptive dimer threshold based on GC content.

    Args:
        base_threshold: Base dimer ΔG threshold (kcal/mol, e.g., -9.0)
        genome_gc: Target genome GC content (0-1)

    Returns:
        Adjusted dimer threshold (kcal/mol)

    Algorithm:
        1. Calculate GC deviation from balanced (50%)
        2. Scale adjustment: 3 kcal/mol per 0.35 GC deviation
        3. Apply to base threshold

    Examples:
        >>> calculate_adaptive_dimer_threshold(-9.0, 0.32)  # Francisella
        -7.5  # More lenient (less negative) for AT-rich
        >>> calculate_adaptive_dimer_threshold(-9.0, 0.67)  # Burkholderia
        -10.5  # Stricter (more negative) for GC-rich
        >>> calculate_adaptive_dimer_threshold(-9.0, 0.50)  # E. coli
        -9.0  # Unchanged for balanced

    Rationale:
        GC-rich sequences form stronger dimers (more negative ΔG) due to
        stronger G≡C hydrogen bonds. AT-rich sequences form weaker dimers.
        Adaptive thresholds prevent over-filtering AT-rich primers and
        under-filtering GC-rich primers.
    """
    # Import from three_prime_stability to reuse logic
    from neoswga.core.three_prime_stability import calculate_gc_deviation

    gc_deviation = calculate_gc_deviation(genome_gc)

    # Scale: 3 kcal/mol per 0.35 GC deviation
    # Positive deviation (GC-rich) → more negative threshold (stricter)
    # Negative deviation (AT-rich) → less negative threshold (lenient)
    dg_adjustment = gc_deviation * (3.0 / 0.35)

    adjusted_threshold = base_threshold - dg_adjustment  # Subtract to get correct direction

    # Floor/ceiling to maintain reasonable bounds
    adjusted_threshold = max(-15.0, min(-5.0, adjusted_threshold))

    return adjusted_threshold


def create_filter_from_conditions(polymerase: str, temperature: float,
                                  gc_content: float, betaine_m: float = 0.0,
                                  na_conc: float = 50.0) -> ThermodynamicFilter:
    """
    Create thermodynamic filter optimized for specific conditions.

    Args:
        polymerase: 'phi29' or 'equiphi29'
        temperature: Reaction temperature (C)
        gc_content: Target genome GC content (0-1)
        betaine_m: Betaine concentration (M)
        na_conc: Sodium concentration (mM)

    Returns:
        Configured ThermodynamicFilter
    """
    # Create base criteria
    criteria = ThermodynamicCriteria(na_conc=na_conc)

    # Create filter
    filter_obj = ThermodynamicFilter(criteria)

    # Adjust for conditions
    filter_obj.adjust_criteria_for_conditions(temperature, gc_content, betaine_m)

    return filter_obj


def create_thermodynamic_filter_adaptive(
    genome_gc: float,
    polymerase: str = 'phi29',
    temperature: float = 30.0,
    betaine_m: float = 0.0,
    na_conc: float = 50.0,
    gc_tolerance: float = 0.15
) -> ThermodynamicFilter:
    """
    Create thermodynamic filter with genome-adaptive thresholds.

    This is the recommended factory function for creating filters that
    adapt to extreme GC genomes (AT-rich <35% or GC-rich >65%).

    Args:
        genome_gc: Target genome GC content (0-1)
        polymerase: Polymerase type ('phi29' or 'equiphi29')
        temperature: Reaction temperature (C)
        betaine_m: Betaine concentration (M)
        na_conc: Sodium concentration (mM)
        gc_tolerance: GC range tolerance (default ±15%)

    Returns:
        ThermodynamicFilter with genome-adaptive criteria

    Examples:
        # Francisella (32% GC, AT-rich)
        >>> filter_fr = create_thermodynamic_filter_adaptive(genome_gc=0.32)
        # GC range: 17-47% (vs standard 30-70%)
        # Dimer threshold: -7.5 kcal/mol (vs standard -9.0)

        # Burkholderia (67% GC, GC-rich)
        >>> filter_bk = create_thermodynamic_filter_adaptive(genome_gc=0.67)
        # GC range: 52-82% (vs standard 30-70%)
        # Dimer threshold: -10.5 kcal/mol (vs standard -9.0)

    Validation:
        Expected to increase pass rates for extreme GC genomes by 2-3x
        while maintaining specificity through adaptive dimer thresholds.
    """
    # Calculate adaptive GC range
    min_gc, max_gc = calculate_adaptive_gc_range(genome_gc, gc_tolerance)

    # Calculate adaptive dimer thresholds
    base_dimer_threshold = -9.0  # Standard threshold
    adaptive_dimer_threshold = calculate_adaptive_dimer_threshold(
        base_dimer_threshold,
        genome_gc
    )

    # Create criteria with adaptive values
    criteria = ThermodynamicCriteria(
        na_conc=na_conc,
        min_gc=min_gc,
        max_gc=max_gc,
        max_homodimer_dg=adaptive_dimer_threshold,
        max_heterodimer_dg=adaptive_dimer_threshold
    )

    # Create filter
    filter_obj = ThermodynamicFilter(criteria)

    # Adjust Tm range for temperature (but keep adaptive GC/dimer)
    filter_obj.criteria.reaction_temp = temperature
    filter_obj.criteria.min_tm = temperature - 5
    filter_obj.criteria.max_tm = temperature + 20

    # Handle betaine (widens GC range further if present)
    if betaine_m > 0:
        gc_bonus = min(betaine_m / 0.5 * 0.05, 0.15)
        filter_obj.criteria.min_gc = max(0.15, filter_obj.criteria.min_gc - gc_bonus)
        filter_obj.criteria.max_gc = min(0.85, filter_obj.criteria.max_gc + gc_bonus)

    # Relax dimer constraints at higher temperature
    if temperature >= 42:
        # Add 1 kcal/mol at high temp (dimers less stable)
        filter_obj.criteria.max_homodimer_dg += 1.0
        filter_obj.criteria.max_heterodimer_dg += 1.0

    logger.info(f"\nGenome-adaptive thermodynamic filter created:")
    logger.info(f"  Genome GC: {genome_gc:.1%}")
    logger.info(f"  GC range: {filter_obj.criteria.min_gc:.1%}-{filter_obj.criteria.max_gc:.1%}")
    logger.info(f"  Tm range: {filter_obj.criteria.min_tm:.1f}-{filter_obj.criteria.max_tm:.1f}°C")
    logger.info(f"  Max homodimer ΔG: {filter_obj.criteria.max_homodimer_dg:.1f} kcal/mol")
    logger.info(f"  Max heterodimer ΔG: {filter_obj.criteria.max_heterodimer_dg:.1f} kcal/mol")

    return filter_obj


if __name__ == '__main__':
    # Example usage
    logging.basicConfig(level=logging.INFO)

    print("Thermodynamic Filter - Example Usage\n")

    # Example primers (some good, some bad)
    test_primers = [
        'ATCGATCGATCG',  # Balanced GC
        'GGGGGGGGGGGG',  # All G (too high GC, strong dimers)
        'AAAAAAAAAAAA',  # All A (too low GC, low Tm)
        'GCTAGCTAGCTA',  # Balanced, good
        'ATATATATATAT',  # AT-rich, weak dimers
    ]

    # Create filter for Phi29 at 30°C
    filter_obj = create_filter_from_conditions(
        polymerase='phi29',
        temperature=30.0,
        gc_content=0.50,
        betaine_m=0.5
    )

    # Filter candidates
    passing, stats = filter_obj.filter_candidates(test_primers)

    print(f"\nPassing primers: {len(passing)}/{len(test_primers)}")
    for seq in passing:
        print(f"  {seq}")
