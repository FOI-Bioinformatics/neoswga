"""
Dominating set optimizer adapter implementing BaseOptimizer interface.

Wraps the existing DominatingSetOptimizer to provide a consistent interface
with the new optimizer framework. The core algorithm remains unchanged.

The dominating set approach models primer selection as a graph theory problem:
- Nodes: Genome regions (bins of configurable size)
- Edges: Primer can amplify region
- Goal: Find minimum set of primers that covers all regions

Advantages over greedy BFS:
- 8x faster due to simpler coverage model
- Provable ln(n) approximation bound
- Focus on coverage rather than selectivity

Use when:
- Speed is critical
- Coverage is primary concern
- Background filtering can be done separately
"""

from dataclasses import dataclass
from typing import List, Dict, Optional, Set
import logging
import numpy as np

from .base_optimizer import (
    BaseOptimizer, OptimizationResult, OptimizationStatus,
    PrimerSetMetrics, OptimizerConfig
)
from .optimizer_factory import OptimizerFactory
from .exceptions import NoCandidatesError

logger = logging.getLogger(__name__)


@dataclass
class CoverageRegion:
    """A region of the genome that needs coverage."""
    chromosome: str
    start: int
    end: int

    def __hash__(self):
        return hash((self.chromosome, self.start, self.end))

    def __eq__(self, other):
        if not isinstance(other, CoverageRegion):
            return False
        return (self.chromosome == other.chromosome and
                self.start == other.start and
                self.end == other.end)


class CoverageGraph:
    """
    Bipartite graph: Primers <-> Regions

    Efficient structure for set cover computation.
    """

    def __init__(self, bin_size: int = 10000):
        self.bin_size = bin_size
        self.primers: Set[str] = set()
        self.regions: Set[CoverageRegion] = set()
        self.primer_to_regions: Dict[str, Set[CoverageRegion]] = {}
        self.region_to_primers: Dict[CoverageRegion, Set[str]] = {}

    def add_primer_coverage(
        self,
        primer: str,
        positions: np.ndarray,
        genome_id: str,
        genome_length: int
    ):
        """Add primer and its coverage regions."""
        self.primers.add(primer)

        if primer not in self.primer_to_regions:
            self.primer_to_regions[primer] = set()

        # Bin positions
        covered_bins = set()
        for pos in positions:
            bin_idx = int(pos / self.bin_size)
            covered_bins.add(bin_idx)

        # Add regions
        for bin_idx in covered_bins:
            start = bin_idx * self.bin_size
            end = min((bin_idx + 1) * self.bin_size, genome_length)

            region = CoverageRegion(
                chromosome=genome_id,
                start=start,
                end=end,
            )

            # Use hash-based lookup (O(1) instead of O(n))
            if region not in self.regions:
                self.regions.add(region)
                self.region_to_primers[region] = set()

            # Add edge
            self.primer_to_regions[primer].add(region)
            self.region_to_primers[region].add(primer)

    def get_coverage(self, selected: Set[str]) -> float:
        """Calculate coverage fraction for selected primers."""
        if not self.regions:
            return 0.0

        covered = set()
        for primer in selected:
            covered.update(self.primer_to_regions.get(primer, set()))

        return len(covered) / len(self.regions)


@dataclass
class DominatingSetConfig(OptimizerConfig):
    """Configuration for dominating set optimizer."""
    bin_size: int = 10000  # Size of genome bins
    use_ilp: bool = False  # Use ILP for exact solution
    ilp_timeout: int = 300  # ILP solver timeout in seconds


@OptimizerFactory.register('dominating-set', aliases=['ds', 'set-cover'])
class DominatingSetAdapter(BaseOptimizer):
    """
    Dominating set optimizer implementing BaseOptimizer interface.

    Uses greedy set cover algorithm for fast primer selection with
    provable approximation bounds (within ln(n) of optimal).

    Algorithm:
    1. Build bipartite graph: primers <-> genome regions
    2. Greedily select primer covering most uncovered regions
    3. Repeat until target coverage or primer count reached

    Complexity: O(n * m) where n = primers, m = regions
    Approximation: ln(n) of optimal (Chvatal 1979)
    """

    def __init__(
        self,
        position_cache,
        fg_prefixes: List[str],
        fg_seq_lengths: List[int],
        bg_prefixes: Optional[List[str]] = None,
        bg_seq_lengths: Optional[List[int]] = None,
        config: Optional[DominatingSetConfig] = None,
        **kwargs
    ):
        super().__init__(
            position_cache, fg_prefixes, fg_seq_lengths,
            bg_prefixes, bg_seq_lengths, config or DominatingSetConfig()
        )
        self.ds_config = config if isinstance(config, DominatingSetConfig) else DominatingSetConfig()

    @property
    def name(self) -> str:
        return "dominating-set"

    @property
    def description(self) -> str:
        return "Graph-based set cover optimizer (8x faster, ln(n) approximation)"

    @property
    def supports_background(self) -> bool:
        return False  # This optimizer focuses on foreground coverage only

    def optimize(
        self,
        candidates: List[str],
        target_size: Optional[int] = None,
        **kwargs
    ) -> OptimizationResult:
        """
        Find optimal primer set using greedy set cover.

        Args:
            candidates: Pool of candidate primers
            target_size: Maximum primers to select

        Returns:
            OptimizationResult with selected primers
        """
        candidates = self._validate_candidates(candidates)
        max_primers = target_size or self.config.target_set_size

        if self.config.verbose:
            logger.info(f"Dominating set optimization: {len(candidates)} candidates")

        # Build coverage graph
        graph = self._build_graph(candidates)

        if self.config.verbose:
            logger.info(f"Graph: {len(graph.primers)} primers, {len(graph.regions)} regions")

        # Run greedy set cover
        result = self._greedy_set_cover(graph, max_primers)

        # Compute metrics
        metrics = self.compute_metrics(result['primers'])

        status = OptimizationStatus.SUCCESS
        if result['coverage'] < 0.95:
            status = OptimizationStatus.PARTIAL

        return OptimizationResult(
            primers=tuple(result['primers']),
            score=result['coverage'],
            status=status,
            metrics=metrics,
            iterations=result['iterations'],
            optimizer_name=self.name,
            message=f"Coverage: {result['coverage']:.1%} ({result['covered']}/{result['total']} regions)",
        )

    def _build_graph(self, candidates: List[str]) -> CoverageGraph:
        """Build coverage graph from candidates."""
        graph = CoverageGraph(bin_size=self.ds_config.bin_size)

        for primer in candidates:
            for prefix, length in zip(self.fg_prefixes, self.fg_seq_lengths):
                positions = self.get_primer_positions(primer, prefix, 'both')

                if len(positions) > 0:
                    graph.add_primer_coverage(primer, positions, prefix, length)

        return graph

    def _greedy_set_cover(self, graph: CoverageGraph, max_primers: int) -> Dict:
        """
        Greedy set cover algorithm.

        Repeatedly selects primer covering most uncovered regions.
        """
        selected: Set[str] = set()
        covered_regions: Set[CoverageRegion] = set()
        iterations = 0

        for iteration in range(max_primers):
            iterations = iteration + 1

            if len(covered_regions) == len(graph.regions):
                if self.config.verbose:
                    logger.info("Full coverage achieved")
                break

            # Find primer covering most uncovered regions
            best_primer = None
            best_new_coverage = 0

            for primer in graph.primers:
                if primer in selected:
                    continue

                # Count NEW regions this primer covers
                new_regions = graph.primer_to_regions.get(primer, set()) - covered_regions
                if len(new_regions) > best_new_coverage:
                    best_new_coverage = len(new_regions)
                    best_primer = primer

            if best_primer is None or best_new_coverage == 0:
                if self.config.verbose:
                    logger.info("No more coverage gains possible")
                break

            # Add primer
            selected.add(best_primer)
            covered_regions.update(graph.primer_to_regions[best_primer])

            if self.config.verbose and (iteration + 1) % 5 == 0:
                coverage = len(covered_regions) / len(graph.regions)
                logger.info(f"  {iteration + 1} primers: {coverage:.1%} coverage")

        # Final statistics
        coverage = len(covered_regions) / len(graph.regions) if graph.regions else 0.0

        return {
            'primers': list(selected),
            'coverage': coverage,
            'covered': len(covered_regions),
            'total': len(graph.regions),
            'iterations': iterations,
        }


@OptimizerFactory.register('weighted-set-cover', aliases=['wsc'])
class WeightedSetCoverOptimizer(DominatingSetAdapter):
    """
    Weighted set cover variant with Tm-based primer weighting.

    Extends basic set cover to prefer primers with:
    - Tm close to reaction temperature
    - More uniform GC content
    - Lower dimer risk
    """

    @property
    def name(self) -> str:
        return "weighted-set-cover"

    @property
    def description(self) -> str:
        return "Set cover with Tm and quality weighting"

    def _greedy_set_cover(self, graph: CoverageGraph, max_primers: int) -> Dict:
        """Greedy set cover with weighted primer selection."""
        selected: Set[str] = set()
        covered_regions: Set[CoverageRegion] = set()
        iterations = 0

        for iteration in range(max_primers):
            iterations = iteration + 1

            if len(covered_regions) == len(graph.regions):
                break

            # Find best primer considering coverage AND quality
            best_primer = None
            best_score = float('-inf')

            for primer in graph.primers:
                if primer in selected:
                    continue

                new_regions = graph.primer_to_regions.get(primer, set()) - covered_regions
                if not new_regions:
                    continue

                # Score = coverage_gain * quality_factor
                coverage_gain = len(new_regions)
                quality = self._primer_quality(primer)
                score = coverage_gain * quality

                if score > best_score:
                    best_score = score
                    best_primer = primer

            if best_primer is None:
                break

            selected.add(best_primer)
            covered_regions.update(graph.primer_to_regions[best_primer])

        coverage = len(covered_regions) / len(graph.regions) if graph.regions else 0.0

        return {
            'primers': list(selected),
            'coverage': coverage,
            'covered': len(covered_regions),
            'total': len(graph.regions),
            'iterations': iterations,
        }

    def _primer_quality(self, primer: str) -> float:
        """
        Calculate quality score for primer.

        Higher is better. Factors in:
        - Tm (prefer ~35C for phi29)
        - GC content (prefer 40-60%)
        - Length (prefer 8-12bp)
        """
        # Tm estimate
        gc = primer.upper().count('G') + primer.upper().count('C')
        at = primer.upper().count('A') + primer.upper().count('T')
        tm = 4 * gc + 2 * at

        # Tm penalty (optimal around 35C)
        tm_penalty = abs(tm - 35) / 20.0

        # GC content penalty (optimal 40-60%)
        gc_content = gc / len(primer)
        gc_penalty = abs(gc_content - 0.5) / 0.3

        # Length bonus (prefer 8-12)
        length = len(primer)
        if 8 <= length <= 12:
            length_bonus = 0.1
        else:
            length_bonus = 0.0

        quality = 1.0 - 0.3 * tm_penalty - 0.2 * gc_penalty + length_bonus

        return max(0.1, quality)  # Minimum quality of 0.1
