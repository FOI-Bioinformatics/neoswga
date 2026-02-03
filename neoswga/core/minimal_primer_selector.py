#!/usr/bin/env python3
"""
Minimal Primer Set Selection via Dominating Set Algorithm.

This module implements graph-based algorithms for selecting the minimum number
of primers needed to achieve target genome coverage. Uses weighted set cover
and dominating set approximations.

Key Features:
- Greedy weighted set cover (log(n) approximation)
- Quality-weighted primer selection
- Coverage uniformity optimization
- Cooperative binding consideration
- Gap analysis and filling

Author: NeoSWGA Development Team
Date: November 2025
Version: 3.5 - Genome-Adaptive QA System
"""

import logging
import numpy as np
from typing import List, Dict, Set, Tuple, Optional
from dataclasses import dataclass
from collections import defaultdict

logger = logging.getLogger(__name__)


@dataclass
class Primer:
    """Primer with binding sites and quality metrics."""
    sequence: str
    binding_sites: List[int]  # Genome positions
    quality_score: float  # 0-1, higher = better
    gc_content: float
    primer_id: Optional[str] = None

    def __post_init__(self):
        if self.primer_id is None:
            self.primer_id = self.sequence


@dataclass
class SelectionResult:
    """Result from minimal primer selection."""
    selected_primers: List[Primer]
    coverage: float  # Fraction of genome covered
    coverage_uniformity: float  # 0-1, higher = more uniform
    total_binding_sites: int
    covered_positions: Set[int]
    gaps: List[Tuple[int, int]]  # (start, end) of uncovered regions
    mean_quality: float
    selection_iterations: int


class MinimalPrimerSelector:
    """
    Select minimal primer set using dominating set approximation.

    Uses weighted greedy set cover algorithm with quality scores
    and coverage uniformity optimization.
    """

    def __init__(self, genome_length: int, min_spacing: int = 1000,
                 quality_weight: float = 0.3):
        """
        Initialize selector.

        Args:
            genome_length: Total genome length (bp)
            min_spacing: Minimum spacing between binding sites (bp)
            quality_weight: Weight for quality vs coverage (0-1)
                          0 = pure coverage maximization
                          1 = pure quality maximization
                          0.3 = balanced (default)
        """
        self.genome_length = genome_length
        self.min_spacing = min_spacing
        self.quality_weight = quality_weight

        logger.info(f"Initialized MinimalPrimerSelector:")
        logger.info(f"  Genome length: {genome_length:,} bp")
        logger.info(f"  Min spacing: {min_spacing:,} bp")
        logger.info(f"  Quality weight: {quality_weight:.2f}")

    def select_minimal_set(self, primers: List[Primer],
                          target_coverage: float = 0.70,
                          max_primers: int = 20) -> SelectionResult:
        """
        Select minimal primer set using weighted greedy algorithm.

        Algorithm:
        1. Start with empty set
        2. Repeatedly select primer with best coverage/quality ratio
        3. Continue until target coverage reached or max primers used

        Args:
            primers: List of candidate primers with binding sites
            target_coverage: Target genome coverage (0-1)
            max_primers: Maximum primers to select

        Returns:
            SelectionResult with selected primers and metrics
        """
        logger.info(f"\nStarting minimal primer selection:")
        logger.info(f"  Candidates: {len(primers)}")
        logger.info(f"  Target coverage: {target_coverage:.1%}")
        logger.info(f"  Max primers: {max_primers}")

        # Initialize
        selected = []
        covered_positions = set()
        remaining_primers = primers.copy()
        iteration = 0

        # Greedy selection loop
        while len(selected) < max_primers and len(remaining_primers) > 0:
            iteration += 1

            # Calculate current coverage
            current_coverage = len(covered_positions) / self.genome_length

            # Check if target reached
            if current_coverage >= target_coverage:
                logger.info(f"  Target coverage {target_coverage:.1%} reached at iteration {iteration}")
                break

            # Find best primer (coverage + quality weighted)
            best_primer = self._select_best_primer(
                remaining_primers,
                covered_positions
            )

            if best_primer is None:
                logger.warning("  No more primers provide additional coverage")
                break

            # Add primer
            selected.append(best_primer)
            remaining_primers.remove(best_primer)

            # Update covered positions
            new_positions = set(best_primer.binding_sites) - covered_positions
            covered_positions.update(new_positions)

            # Log progress
            coverage = len(covered_positions) / self.genome_length
            logger.info(
                f"  Iteration {iteration}: Selected primer {best_primer.primer_id} "
                f"(+{len(new_positions)} sites, coverage={coverage:.1%})"
            )

        # Calculate final metrics
        result = self._calculate_metrics(selected, covered_positions, iteration)

        # Log summary
        self._log_selection_summary(result)

        return result

    def _select_best_primer(self, primers: List[Primer],
                           covered_positions: Set[int]) -> Optional[Primer]:
        """
        Select primer with best coverage/quality score.

        Score = (1 - quality_weight) * coverage_gain + quality_weight * quality

        Args:
            primers: Available primers
            covered_positions: Already covered positions

        Returns:
            Best primer or None if no primers provide coverage
        """
        best_primer = None
        best_score = -1.0

        for primer in primers:
            # Calculate new positions covered
            new_positions = set(primer.binding_sites) - covered_positions
            coverage_gain = len(new_positions)

            if coverage_gain == 0:
                continue  # Skip primers with no new coverage

            # Normalize coverage gain by genome length
            normalized_coverage = coverage_gain / self.genome_length

            # Combined score (coverage + quality)
            score = (
                (1 - self.quality_weight) * normalized_coverage +
                self.quality_weight * primer.quality_score
            )

            if score > best_score:
                best_score = score
                best_primer = primer

        return best_primer

    def _calculate_metrics(self, selected: List[Primer],
                          covered_positions: Set[int],
                          iterations: int) -> SelectionResult:
        """Calculate comprehensive metrics for selected primer set."""

        # Coverage
        coverage = len(covered_positions) / self.genome_length

        # Total binding sites
        total_sites = sum(len(p.binding_sites) for p in selected)

        # Coverage uniformity (Gini coefficient)
        uniformity = self._calculate_coverage_uniformity(covered_positions)

        # Find gaps
        gaps = self._find_coverage_gaps(covered_positions)

        # Mean quality
        mean_quality = np.mean([p.quality_score for p in selected]) if selected else 0.0

        return SelectionResult(
            selected_primers=selected,
            coverage=coverage,
            coverage_uniformity=uniformity,
            total_binding_sites=total_sites,
            covered_positions=covered_positions,
            gaps=gaps,
            mean_quality=mean_quality,
            selection_iterations=iterations
        )

    def _calculate_coverage_uniformity(self, covered_positions: Set[int]) -> float:
        """
        Calculate coverage uniformity using Gini coefficient.

        Divides genome into windows and measures evenness of coverage.

        Returns:
            Uniformity score (0-1, higher = more uniform)
        """
        if not covered_positions:
            return 0.0

        # Divide genome into 100 windows
        n_windows = 100
        window_size = self.genome_length // n_windows
        window_counts = np.zeros(n_windows)

        # Count positions in each window
        for pos in covered_positions:
            window_idx = min(pos // window_size, n_windows - 1)
            window_counts[window_idx] += 1

        # Calculate Gini coefficient
        if window_counts.sum() == 0:
            return 0.0

        # Sort counts
        sorted_counts = np.sort(window_counts)

        # Gini = 1 - 2 * area under Lorenz curve
        n = len(sorted_counts)
        cumsum = np.cumsum(sorted_counts)
        total = cumsum[-1]

        if total == 0:
            return 0.0

        # Lorenz curve area
        lorenz_area = np.sum(cumsum) / (n * total)

        # Gini coefficient
        gini = 1 - 2 * lorenz_area

        # Uniformity = 1 - Gini (higher = more uniform)
        uniformity = 1 - gini

        return uniformity

    def _find_coverage_gaps(self, covered_positions: Set[int],
                           min_gap_size: int = 10000) -> List[Tuple[int, int]]:
        """
        Find large uncovered regions in genome.

        Args:
            covered_positions: Set of covered positions
            min_gap_size: Minimum gap size to report (bp)

        Returns:
            List of (start, end) tuples for large gaps
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

    def _log_selection_summary(self, result: SelectionResult):
        """Log detailed selection summary."""
        logger.info(f"\n{'='*60}")
        logger.info(f"MINIMAL PRIMER SELECTION COMPLETE")
        logger.info(f"{'='*60}")
        logger.info(f"Primers selected: {len(result.selected_primers)}")
        logger.info(f"Coverage: {result.coverage:.1%}")
        logger.info(f"Coverage uniformity: {result.coverage_uniformity:.2f}")
        logger.info(f"Total binding sites: {result.total_binding_sites}")
        logger.info(f"Mean primer quality: {result.mean_quality:.3f}")
        logger.info(f"Selection iterations: {result.selection_iterations}")

        if result.gaps:
            logger.info(f"\nCoverage gaps (>{10000} bp):")
            for i, (start, end) in enumerate(result.gaps[:5], 1):
                gap_size = end - start
                logger.info(f"  Gap {i}: {start:,} - {end:,} ({gap_size:,} bp)")
            if len(result.gaps) > 5:
                logger.info(f"  ... and {len(result.gaps) - 5} more gaps")
        else:
            logger.info(f"\nNo large coverage gaps (excellent coverage!)")

        logger.info(f"{'='*60}\n")

    def optimize_for_uniformity(self, primers: List[Primer],
                               target_primers: int = 15) -> SelectionResult:
        """
        Select primers optimized for coverage uniformity.

        Uses a different strategy that prioritizes filling gaps.

        Args:
            primers: Candidate primers
            target_primers: Number of primers to select

        Returns:
            SelectionResult with uniform coverage
        """
        logger.info(f"\nOptimizing primer selection for uniformity:")
        logger.info(f"  Target primers: {target_primers}")

        selected = []
        covered_positions = set()

        # Divide genome into regions
        n_regions = target_primers
        region_size = self.genome_length // n_regions

        # Select best primer for each region
        for region_idx in range(n_regions):
            region_start = region_idx * region_size
            region_end = region_start + region_size

            # Find primers with sites in this region
            region_primers = [
                p for p in primers
                if any(region_start <= site < region_end for site in p.binding_sites)
            ]

            if not region_primers:
                logger.warning(f"  No primers for region {region_idx} ({region_start:,}-{region_end:,})")
                continue

            # Select best quality primer in region
            best_primer = max(region_primers, key=lambda p: p.quality_score)
            selected.append(best_primer)
            covered_positions.update(best_primer.binding_sites)

            logger.info(
                f"  Region {region_idx}: Selected {best_primer.primer_id} "
                f"(quality={best_primer.quality_score:.3f})"
            )

        # Calculate metrics
        result = self._calculate_metrics(selected, covered_positions, n_regions)
        self._log_selection_summary(result)

        return result


def create_primers_from_sequences(sequences: List[str],
                                  binding_sites_map: Dict[str, List[int]],
                                  quality_scores: Dict[str, float]) -> List[Primer]:
    """
    Create Primer objects from sequences and binding data.

    Args:
        sequences: List of primer sequences
        binding_sites_map: Dict mapping sequence -> list of binding positions
        quality_scores: Dict mapping sequence -> quality score (0-1)

    Returns:
        List of Primer objects
    """
    primers = []

    for seq in sequences:
        sites = binding_sites_map.get(seq, [])
        quality = quality_scores.get(seq, 0.5)  # Default 0.5

        # Calculate GC content
        gc_count = seq.count('G') + seq.count('C')
        gc_content = gc_count / len(seq) if len(seq) > 0 else 0.5

        primer = Primer(
            sequence=seq,
            binding_sites=sites,
            quality_score=quality,
            gc_content=gc_content,
            primer_id=seq
        )

        primers.append(primer)

    return primers


if __name__ == '__main__':
    # Example usage
    logging.basicConfig(level=logging.INFO)

    print("NeoSWGA Minimal Primer Selector - Example\n")

    # Create example primers
    genome_length = 2_000_000  # 2 Mbp

    # Simulate 50 primers with random binding sites
    import random
    random.seed(42)

    primers = []
    for i in range(50):
        # Random binding sites
        n_sites = random.randint(2, 8)
        sites = sorted([random.randint(0, genome_length - 1) for _ in range(n_sites)])

        # Random quality
        quality = random.uniform(0.5, 1.0)

        # Random sequence
        bases = ['A', 'T', 'G', 'C']
        seq = ''.join(random.choice(bases) for _ in range(12))

        primer = Primer(
            sequence=seq,
            binding_sites=sites,
            quality_score=quality,
            gc_content=random.uniform(0.3, 0.7),
            primer_id=f"Primer_{i+1}"
        )
        primers.append(primer)

    # Create selector
    selector = MinimalPrimerSelector(
        genome_length=genome_length,
        min_spacing=1000,
        quality_weight=0.3
    )

    # Select minimal set
    result = selector.select_minimal_set(
        primers=primers,
        target_coverage=0.60,
        max_primers=15
    )

    print(f"\nSelected {len(result.selected_primers)} primers:")
    for primer in result.selected_primers[:10]:
        print(f"  {primer.primer_id}: {len(primer.binding_sites)} sites, quality={primer.quality_score:.3f}")

    print(f"\nCoverage: {result.coverage:.1%}")
    print(f"Uniformity: {result.coverage_uniformity:.2f}")
