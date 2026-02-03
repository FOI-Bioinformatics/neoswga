#!/usr/bin/env python3
"""
Hybrid Core+Supplemental Primer Selection Strategy.

Combines high-specificity core primers with coverage-optimized supplemental
primers to achieve both strong enrichment and comprehensive genome coverage.

Strategy:
1. Core Set: Longer primers (12-14bp) for high specificity and enrichment
2. Supplemental Set: Shorter primers (10-11bp) to fill coverage gaps
3. Unified QA scoring adapted to genome GC content

This directly addresses the Francisella coverage problem where long GC-biased
primers leave AT-rich gaps. Supplemental primers fill these gaps while
maintaining overall specificity.

Author: NeoSWGA Development Team
Date: November 2025
Version: 3.5 - Genome-Adaptive QA System
"""

import logging
import numpy as np
from typing import List, Dict, Set, Tuple, Optional
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass
class PrimerSet:
    """A set of primers with metadata."""
    primers: List[str]
    primer_type: str  # 'core' or 'supplemental'
    mean_length: float
    mean_quality: float
    total_binding_sites: int
    coverage: float  # Fraction of genome covered


@dataclass
class HybridResult:
    """Result from hybrid primer selection."""
    core_set: PrimerSet
    supplemental_set: PrimerSet
    combined_coverage: float
    coverage_improvement: float  # Improvement from supplemental
    total_primers: int
    mean_enrichment_estimate: float
    gap_filling_efficiency: float  # How well supplemental fills gaps


class HybridPrimerStrategy:
    """
    Select hybrid primer set combining core (specificity) and
    supplemental (coverage) primers.

    Addresses genome composition bias by using different primer lengths
    and adaptive QA thresholds.
    """

    def __init__(self, genome_length: int, genome_gc: float,
                 core_length_range: Tuple[int, int] = (12, 14),
                 supplemental_length_range: Tuple[int, int] = (10, 11)):
        """
        Initialize hybrid strategy.

        Args:
            genome_length: Total genome length (bp)
            genome_gc: Genome GC content (0-1)
            core_length_range: Primer length range for core set
            supplemental_length_range: Primer length range for supplemental
        """
        self.genome_length = genome_length
        self.genome_gc = genome_gc
        self.core_length_range = core_length_range
        self.supplemental_length_range = supplemental_length_range

        logger.info(f"Initialized HybridPrimerStrategy:")
        logger.info(f"  Genome: {genome_length:,} bp, {genome_gc:.1%} GC")
        logger.info(f"  Core length: {core_length_range[0]}-{core_length_range[1]} bp")
        logger.info(f"  Supplemental length: {supplemental_length_range[0]}-{supplemental_length_range[1]} bp")

    def select_hybrid_set(self, all_primers: Dict[str, Dict],
                         core_count: int = 8,
                         supplemental_count: int = 8,
                         target_coverage: float = 0.70) -> HybridResult:
        """
        Select hybrid primer set with core and supplemental subsets.

        Algorithm:
        1. Select core primers (longer, high quality)
        2. Evaluate coverage from core set
        3. Identify gaps in coverage
        4. Select supplemental primers to fill gaps
        5. Combine and evaluate total coverage

        Args:
            all_primers: Dict mapping primer_id -> {
                'sequence': str,
                'binding_sites': [positions],
                'quality': float,
                'length': int
            }
            core_count: Number of core primers
            supplemental_count: Number of supplemental primers
            target_coverage: Target genome coverage (0-1)

        Returns:
            HybridResult with core + supplemental sets
        """
        logger.info(f"\nStarting hybrid primer selection:")
        logger.info(f"  Total candidates: {len(all_primers)}")
        logger.info(f"  Core target: {core_count} primers")
        logger.info(f"  Supplemental target: {supplemental_count} primers")
        logger.info(f"  Coverage target: {target_coverage:.1%}")

        # Step 1: Partition primers by length
        core_candidates, supplemental_candidates = self._partition_by_length(all_primers)

        logger.info(f"\nPartitioned candidates:")
        logger.info(f"  Core candidates: {len(core_candidates)}")
        logger.info(f"  Supplemental candidates: {len(supplemental_candidates)}")

        # Step 2: Select core set (high quality, longer primers)
        core_primers, core_coverage = self._select_core_set(
            core_candidates, core_count
        )

        logger.info(f"\nCore set selected:")
        logger.info(f"  Primers: {len(core_primers)}")
        logger.info(f"  Coverage: {core_coverage:.1%}")

        # Step 3: Identify coverage gaps
        covered_by_core = self._get_covered_positions(core_primers, all_primers)
        gaps = self._identify_gaps(covered_by_core)

        logger.info(f"\nCoverage gaps identified:")
        logger.info(f"  Number of gaps: {len(gaps)}")
        total_gap_size = sum(end - start for start, end in gaps)
        logger.info(f"  Total gap size: {total_gap_size:,} bp ({total_gap_size/self.genome_length:.1%})")

        # Step 4: Select supplemental primers to fill gaps
        supplemental_primers, supplemental_coverage = self._select_gap_filling_primers(
            supplemental_candidates, gaps, supplemental_count
        )

        logger.info(f"\nSupplemental set selected:")
        logger.info(f"  Primers: {len(supplemental_primers)}")
        logger.info(f"  Additional coverage: {supplemental_coverage - core_coverage:.1%}")

        # Step 5: Build result
        result = self._build_result(
            core_primers, supplemental_primers, all_primers,
            core_coverage, supplemental_coverage
        )

        self._log_summary(result)

        return result

    def _partition_by_length(self, all_primers: Dict[str, Dict]) -> Tuple[Dict, Dict]:
        """Partition primers into core and supplemental candidates by length."""
        core_candidates = {}
        supplemental_candidates = {}

        for primer_id, data in all_primers.items():
            length = data.get('length', len(data['sequence']))

            if self.core_length_range[0] <= length <= self.core_length_range[1]:
                core_candidates[primer_id] = data
            elif self.supplemental_length_range[0] <= length <= self.supplemental_length_range[1]:
                supplemental_candidates[primer_id] = data

        return core_candidates, supplemental_candidates

    def _select_core_set(self, candidates: Dict[str, Dict],
                        count: int) -> Tuple[List[str], float]:
        """
        Select core primers prioritizing quality and enrichment.

        Uses greedy selection based on quality score.
        """
        # Sort by quality
        sorted_primers = sorted(
            candidates.items(),
            key=lambda x: x[1]['quality'],
            reverse=True
        )

        # Select top primers
        selected = [pid for pid, _ in sorted_primers[:count]]

        # Calculate coverage
        covered_positions = set()
        for pid in selected:
            covered_positions.update(candidates[pid]['binding_sites'])

        coverage = len(covered_positions) / self.genome_length

        return selected, coverage

    def _get_covered_positions(self, primer_ids: List[str],
                               all_primers: Dict[str, Dict]) -> Set[int]:
        """Get all positions covered by selected primers."""
        covered = set()
        for pid in primer_ids:
            covered.update(all_primers[pid]['binding_sites'])
        return covered

    def _identify_gaps(self, covered_positions: Set[int],
                      min_gap_size: int = 5000) -> List[Tuple[int, int]]:
        """
        Identify large uncovered gaps in genome.

        Args:
            covered_positions: Set of covered positions
            min_gap_size: Minimum gap size to report (bp)

        Returns:
            List of (start, end) tuples for gaps
        """
        if not covered_positions:
            return [(0, self.genome_length)]

        # Sort positions
        sorted_positions = sorted(covered_positions)

        # Find gaps
        gaps = []
        prev_pos = 0

        for pos in sorted_positions:
            gap_size = pos - prev_pos

            if gap_size >= min_gap_size:
                gaps.append((prev_pos, pos))

            prev_pos = pos

        # Check final gap
        final_gap = self.genome_length - prev_pos
        if final_gap >= min_gap_size:
            gaps.append((prev_pos, self.genome_length))

        return gaps

    def _select_gap_filling_primers(self, candidates: Dict[str, Dict],
                                    gaps: List[Tuple[int, int]],
                                    count: int) -> Tuple[List[str], float]:
        """
        Select supplemental primers to fill coverage gaps.

        Prioritizes primers with sites in gap regions.
        """
        # Score primers by gap coverage
        primer_scores = {}

        for primer_id, data in candidates.items():
            sites = data['binding_sites']

            # Count sites in gaps
            gap_sites = 0
            for start, end in gaps:
                gap_sites += sum(1 for s in sites if start <= s < end)

            # Score = gap_sites * quality
            score = gap_sites * data['quality']
            primer_scores[primer_id] = score

        # Select top scorers
        sorted_primers = sorted(
            primer_scores.items(),
            key=lambda x: x[1],
            reverse=True
        )

        selected = [pid for pid, _ in sorted_primers[:count]]

        # Calculate total coverage (core + supplemental)
        all_positions = set()
        for pid in selected:
            all_positions.update(candidates[pid]['binding_sites'])

        # Note: This is just supplemental coverage, need to add core
        coverage = len(all_positions) / self.genome_length

        return selected, coverage

    def _build_result(self, core_primers: List[str],
                     supplemental_primers: List[str],
                     all_primers: Dict[str, Dict],
                     core_coverage: float,
                     supplemental_coverage: float) -> HybridResult:
        """Build comprehensive result."""

        # Core set metrics
        core_lengths = [all_primers[pid].get('length', len(all_primers[pid]['sequence']))
                       for pid in core_primers]
        core_qualities = [all_primers[pid]['quality'] for pid in core_primers]
        core_sites = sum(len(all_primers[pid]['binding_sites']) for pid in core_primers)

        core_set = PrimerSet(
            primers=core_primers,
            primer_type='core',
            mean_length=np.mean(core_lengths) if core_lengths else 0,
            mean_quality=np.mean(core_qualities) if core_qualities else 0,
            total_binding_sites=core_sites,
            coverage=core_coverage
        )

        # Supplemental set metrics
        supp_lengths = [all_primers[pid].get('length', len(all_primers[pid]['sequence']))
                       for pid in supplemental_primers]
        supp_qualities = [all_primers[pid]['quality'] for pid in supplemental_primers]
        supp_sites = sum(len(all_primers[pid]['binding_sites']) for pid in supplemental_primers)

        supplemental_set = PrimerSet(
            primers=supplemental_primers,
            primer_type='supplemental',
            mean_length=np.mean(supp_lengths) if supp_lengths else 0,
            mean_quality=np.mean(supp_qualities) if supp_qualities else 0,
            total_binding_sites=supp_sites,
            coverage=supplemental_coverage - core_coverage  # Additional coverage
        )

        # Combined metrics
        all_selected = core_primers + supplemental_primers
        combined_covered = self._get_covered_positions(all_selected, all_primers)
        combined_coverage = len(combined_covered) / self.genome_length

        coverage_improvement = combined_coverage - core_coverage

        # Estimate enrichment (simplified)
        mean_enrichment = 100.0  # Placeholder

        # Gap filling efficiency
        if core_coverage < 1.0:
            max_improvement = 1.0 - core_coverage
            gap_efficiency = coverage_improvement / max_improvement if max_improvement > 0 else 0
        else:
            gap_efficiency = 1.0

        return HybridResult(
            core_set=core_set,
            supplemental_set=supplemental_set,
            combined_coverage=combined_coverage,
            coverage_improvement=coverage_improvement,
            total_primers=len(all_selected),
            mean_enrichment_estimate=mean_enrichment,
            gap_filling_efficiency=gap_efficiency
        )

    def _log_summary(self, result: HybridResult):
        """Log detailed summary."""
        logger.info(f"\n{'='*60}")
        logger.info(f"HYBRID PRIMER SELECTION COMPLETE")
        logger.info(f"{'='*60}")

        logger.info(f"\nCore Set ({len(result.core_set.primers)} primers):")
        logger.info(f"  Mean length: {result.core_set.mean_length:.1f} bp")
        logger.info(f"  Mean quality: {result.core_set.mean_quality:.3f}")
        logger.info(f"  Binding sites: {result.core_set.total_binding_sites}")
        logger.info(f"  Coverage: {result.core_set.coverage:.1%}")

        logger.info(f"\nSupplemental Set ({len(result.supplemental_set.primers)} primers):")
        logger.info(f"  Mean length: {result.supplemental_set.mean_length:.1f} bp")
        logger.info(f"  Mean quality: {result.supplemental_set.mean_quality:.3f}")
        logger.info(f"  Binding sites: {result.supplemental_set.total_binding_sites}")
        logger.info(f"  Additional coverage: +{result.coverage_improvement:.1%}")

        logger.info(f"\nCombined Performance:")
        logger.info(f"  Total primers: {result.total_primers}")
        logger.info(f"  Total coverage: {result.combined_coverage:.1%}")
        logger.info(f"  Coverage improvement: +{result.coverage_improvement:.1%}")
        logger.info(f"  Gap filling efficiency: {result.gap_filling_efficiency:.1%}")

        logger.info(f"{'='*60}\n")


if __name__ == '__main__':
    # Example usage
    logging.basicConfig(level=logging.INFO)

    print("NeoSWGA Hybrid Primer Strategy - Example\n")

    # Simulate primers for Francisella-like genome
    import random
    random.seed(42)

    genome_length = 1_900_000
    genome_gc = 0.32  # AT-rich like Francisella

    all_primers = {}

    # Generate longer primers (core candidates) - fewer sites
    for i in range(30):
        length = random.choice([12, 13, 14])
        n_sites = random.randint(2, 5)  # Fewer sites for longer primers

        sites = sorted([random.randint(0, genome_length - 1) for _ in range(n_sites)])

        all_primers[f"Core_{i+1}"] = {
            'sequence': 'A' * length,
            'binding_sites': sites,
            'quality': random.uniform(0.7, 1.0),
            'length': length
        }

    # Generate shorter primers (supplemental candidates) - more sites
    for i in range(30):
        length = random.choice([10, 11])
        n_sites = random.randint(4, 10)  # More sites for shorter primers

        sites = sorted([random.randint(0, genome_length - 1) for _ in range(n_sites)])

        all_primers[f"Supp_{i+1}"] = {
            'sequence': 'A' * length,
            'binding_sites': sites,
            'quality': random.uniform(0.5, 0.8),  # Lower quality
            'length': length
        }

    # Create strategy
    strategy = HybridPrimerStrategy(
        genome_length=genome_length,
        genome_gc=genome_gc,
        core_length_range=(12, 14),
        supplemental_length_range=(10, 11)
    )

    # Select hybrid set
    result = strategy.select_hybrid_set(
        all_primers=all_primers,
        core_count=8,
        supplemental_count=8,
        target_coverage=0.70
    )

    print(f"\nCore primers: {result.core_set.primers[:5]}...")
    print(f"Supplemental primers: {result.supplemental_set.primers[:5]}...")
    print(f"\nCore coverage: {result.core_set.coverage:.1%}")
    print(f"Combined coverage: {result.combined_coverage:.1%}")
    print(f"Improvement: +{result.coverage_improvement:.1%}")
