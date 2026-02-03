"""
Minimum Dominating Set algorithms for primer coverage optimization.

Models primer selection as a graph theory problem:
- Nodes: Genome regions (bins)
- Edges: Primer can amplify region
- Goal: Find minimum set of primers that covers all regions

This provides:
1. Provable approximation bounds
2. Faster than exhaustive search
3. Optimal coverage guarantees
"""

import numpy as np
import networkx as nx
from typing import List, Dict, Set, Tuple, Optional
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)


@dataclass
class CoverageRegion:
    """A region of the genome that needs coverage"""
    chromosome: str
    start: int
    end: int
    covered_by: Set[str]  # Primers that cover this region

    def __hash__(self):
        return hash((self.chromosome, self.start, self.end))


class BipartiteGraph:
    """
    Bipartite graph: Primers ↔ Regions

    Left nodes: Primers
    Right nodes: Genome regions
    Edges: Primer can amplify region
    """

    def __init__(self, bin_size: int = 10000):
        """
        Initialize bipartite graph.

        Args:
            bin_size: Size of genome bins (default: 10kb)
        """
        self.graph = nx.Graph()
        self.bin_size = bin_size

        # Node sets
        self.primers: Set[str] = set()
        self.regions: Set[CoverageRegion] = set()

        # Coverage map
        self.primer_to_regions: Dict[str, Set[CoverageRegion]] = {}
        self.region_to_primers: Dict[CoverageRegion, Set[str]] = {}

        # Fast lookup for existing regions: (chromosome, start, end) -> CoverageRegion
        self._region_lookup: Dict[Tuple[str, int, int], CoverageRegion] = {}

    def add_primer_coverage(self, primer: str, positions: np.ndarray,
                           genome_id: str, genome_length: int):
        """
        Add primer and its coverage regions.

        Args:
            primer: Primer sequence
            positions: Array of binding positions
            genome_id: Genome identifier
            genome_length: Total genome length
        """
        self.primers.add(primer)

        if primer not in self.primer_to_regions:
            self.primer_to_regions[primer] = set()

        # Create regions based on binning
        n_bins = (genome_length + self.bin_size - 1) // self.bin_size

        covered_bins = set()
        for pos in positions:
            bin_idx = int(pos / self.bin_size)
            covered_bins.add(bin_idx)

        # Add regions
        for bin_idx in covered_bins:
            start = bin_idx * self.bin_size
            end = min((bin_idx + 1) * self.bin_size, genome_length)

            # O(1) lookup for existing region instead of O(n) search
            region_key = (genome_id, start, end)
            existing = self._region_lookup.get(region_key)

            if existing:
                region = existing
            else:
                region = CoverageRegion(
                    chromosome=genome_id,
                    start=start,
                    end=end,
                    covered_by=set()
                )
                self.regions.add(region)
                self.region_to_primers[region] = set()
                self._region_lookup[region_key] = region

            # Add edge
            region.covered_by.add(primer)
            self.primer_to_regions[primer].add(region)
            self.region_to_primers[region].add(primer)

            # Add to graph
            self.graph.add_edge(f"primer_{primer}", f"region_{region.chromosome}_{region.start}")

    def get_uncovered_regions(self, selected_primers: Set[str]) -> Set[CoverageRegion]:
        """Get regions not covered by selected primers"""
        uncovered = set()

        for region in self.regions:
            if not any(p in selected_primers for p in region.covered_by):
                uncovered.add(region)

        return uncovered

    def get_coverage_score(self, selected_primers: Set[str]) -> float:
        """Calculate fraction of regions covered"""
        if not self.regions:
            return 0.0

        covered = len(self.regions) - len(self.get_uncovered_regions(selected_primers))
        return covered / len(self.regions)


class DominatingSetOptimizer:
    """
    Minimum dominating set optimizer for primer selection.

    Uses graph algorithms to find near-optimal primer sets.
    """

    def __init__(self, cache, fg_prefixes: List[str], fg_seq_lengths: List[int],
                 bin_size: int = 10000):
        """
        Initialize optimizer.

        Args:
            cache: PositionCache with primer positions
            fg_prefixes: Target genome prefixes
            fg_seq_lengths: Target genome lengths
            bin_size: Size of genome bins for coverage
        """
        self.cache = cache
        self.fg_prefixes = fg_prefixes
        self.fg_seq_lengths = fg_seq_lengths
        self.bin_size = bin_size

    def optimize_greedy(self, candidates: List[str], max_primers: int = 20,
                       fixed_primers: Optional[List[str]] = None,
                       verbose: bool = True) -> Dict:
        """
        Greedy set cover algorithm.

        Repeatedly selects primer that covers most uncovered regions.
        Approximation ratio: ln(n) where n = number of regions.

        Args:
            candidates: List of candidate primers
            max_primers: Maximum primers to select (new primers only,
                excludes fixed_primers)
            fixed_primers: Optional list of primers that are already selected
                and must be included. Their coverage is pre-computed and
                the algorithm selects additional primers to complement them.
            verbose: Print progress

        Returns:
            Dictionary with selected primers and coverage stats
        """
        fixed_primers = fixed_primers or []
        n_fixed = len(fixed_primers)

        if verbose:
            logger.info(f"Dominating set optimization: {len(candidates)} candidates")
            if n_fixed > 0:
                logger.info(f"  Fixed primers: {n_fixed} (pre-selected)")

        # Build bipartite graph including both candidates and fixed primers
        graph = BipartiteGraph(bin_size=self.bin_size)

        # Add fixed primers first
        for primer in fixed_primers:
            for prefix, length in zip(self.fg_prefixes, self.fg_seq_lengths):
                positions = self.cache.get_positions(prefix, primer, 'both')
                if len(positions) > 0:
                    graph.add_primer_coverage(primer, positions, prefix, length)

        # Add candidate primers
        for primer in candidates:
            for prefix, length in zip(self.fg_prefixes, self.fg_seq_lengths):
                # PositionCache uses 'forward', 'reverse', 'both' for strand
                positions = self.cache.get_positions(prefix, primer, 'both')

                if len(positions) > 0:
                    graph.add_primer_coverage(primer, positions, prefix, length)

        if verbose:
            logger.info(f"Graph: {len(graph.primers)} primers, {len(graph.regions)} regions")

        # Greedy set cover - start with fixed primers' coverage
        selected = set(fixed_primers)
        covered_regions = set()

        # Pre-fill coverage from fixed primers
        for primer in fixed_primers:
            covered_regions.update(graph.primer_to_regions.get(primer, set()))

        if verbose and n_fixed > 0:
            coverage_so_far = len(covered_regions) / len(graph.regions) if graph.regions else 0.0
            logger.info(f"  Fixed primer coverage: {coverage_so_far:.1%}")

        for iteration in range(max_primers):
            if len(covered_regions) == len(graph.regions):
                if verbose:
                    logger.info("Full coverage achieved")
                break

            # Find primer that covers most uncovered regions
            best_primer = None
            best_new_coverage = 0

            for primer in graph.primers:
                if primer in selected:
                    continue

                # Count how many NEW regions this primer covers
                new_regions = graph.primer_to_regions.get(primer, set()) - covered_regions
                if len(new_regions) > best_new_coverage:
                    best_new_coverage = len(new_regions)
                    best_primer = primer

            if best_primer is None or best_new_coverage == 0:
                if verbose:
                    logger.info("No more coverage gains possible")
                break

            # Add primer
            selected.add(best_primer)
            covered_regions.update(graph.primer_to_regions[best_primer])

            if verbose and (iteration + 1) % 5 == 0:
                coverage = len(covered_regions) / len(graph.regions)
                logger.info(f"  {iteration + 1} primers: {coverage:.1%} coverage")

        # Final statistics
        coverage = len(covered_regions) / len(graph.regions) if graph.regions else 0.0
        uncovered = graph.regions - covered_regions

        # Separate fixed and new primers in result
        new_primers = [p for p in selected if p not in fixed_primers]

        result = {
            'primers': list(selected),
            'new_primers': new_primers,
            'fixed_primers': fixed_primers,
            'n_primers': len(selected),
            'n_new': len(new_primers),
            'n_fixed': n_fixed,
            'coverage': coverage,
            'covered_regions': len(covered_regions),
            'total_regions': len(graph.regions),
            'uncovered_regions': len(uncovered),
            'graph': graph
        }

        if verbose:
            if n_fixed > 0:
                logger.info(f"Selected {len(new_primers)} new primers + {n_fixed} fixed = {len(selected)} total")
            else:
                logger.info(f"Selected {len(selected)} primers")
            logger.info(f"Coverage: {coverage:.1%} ({len(covered_regions)}/{len(graph.regions)} regions)")

        return result

    def optimize_ilp(self, candidates: List[str], max_primers: int = 20,
                    verbose: bool = True) -> Optional[Dict]:
        """
        Integer Linear Programming for exact minimum dominating set.

        Finds optimal solution (if feasible within time limit).

        Args:
            candidates: List of candidate primers
            max_primers: Maximum primers to select
            verbose: Print progress

        Returns:
            Dictionary with optimal solution, or None if infeasible
        """
        try:
            from mip import Model, BINARY, minimize, xsum
        except ImportError:
            logger.error("python-mip required for ILP. Install: pip install mip")
            return None

        if verbose:
            logger.info("Building ILP model for minimum dominating set...")

        # Build bipartite graph
        graph = BipartiteGraph(bin_size=self.bin_size)

        for primer in candidates:
            for prefix, length in zip(self.fg_prefixes, self.fg_seq_lengths):
                # PositionCache uses 'forward', 'reverse', 'both' for strand
                positions = self.cache.get_positions(prefix, primer, 'both')

                if len(positions) > 0:
                    graph.add_primer_coverage(primer, positions, prefix, length)

        primers_list = list(graph.primers)
        regions_list = list(graph.regions)

        # Create model
        model = Model(sense=minimize)

        # Decision variables: x[i] = 1 if primer i selected
        x = {primer: model.add_var(var_type=BINARY, name=f"x_{primer}")
             for primer in primers_list}

        # Objective: Minimize number of primers
        model.objective = xsum(x[primer] for primer in primers_list)

        # Constraints: Each region must be covered by at least one primer
        for region in regions_list:
            covering_primers = graph.region_to_primers[region]
            if covering_primers:
                model += xsum(x[primer] for primer in covering_primers
                            if primer in x) >= 1

        # Constraint: Select at most max_primers
        model += xsum(x[primer] for primer in primers_list) <= max_primers

        if verbose:
            logger.info(f"ILP model: {len(primers_list)} primers, {len(regions_list)} regions")
            logger.info("Solving...")

        # Solve
        status = model.optimize(max_seconds=300)

        if status.value == 0:  # Optimal
            selected = [primer for primer in primers_list if x[primer].x >= 0.99]

            covered_regions = set()
            for primer in selected:
                covered_regions.update(graph.primer_to_regions.get(primer, set()))

            coverage = len(covered_regions) / len(graph.regions) if graph.regions else 0.0

            result = {
                'primers': selected,
                'n_primers': len(selected),
                'coverage': coverage,
                'covered_regions': len(covered_regions),
                'total_regions': len(graph.regions),
                'optimal': True,
                'graph': graph
            }

            if verbose:
                logger.info(f"Optimal solution: {len(selected)} primers, {coverage:.1%} coverage")

            return result
        else:
            if verbose:
                logger.warning("ILP did not find optimal solution")
            return None


def optimize(verbose: bool = True, max_time: int = 300) -> Tuple[List[List[str]], List[float]]:
    """
    Standalone optimize function for CLI integration.

    Uses greedy set cover algorithm for fast primer selection with
    provable approximation bounds (within ln(n) of optimal).

    Args:
        verbose: Print detailed progress information
        max_time: Maximum optimization time in seconds (not directly used, for API compatibility)

    Returns:
        Tuple of (primer_sets, scores) where:
        - primer_sets: List containing one list of selected primers
        - scores: List containing the coverage score
    """
    import os
    import pandas as pd
    from neoswga.core import parameter
    from neoswga.core import pipeline as core_pipeline
    from neoswga.core.position_cache import PositionCache

    # Initialize pipeline (loads parameters)
    core_pipeline._initialize()

    # Get parameters from pipeline
    fg_prefixes = core_pipeline.fg_prefixes
    fg_seq_lengths = core_pipeline.fg_seq_lengths

    # Load candidates from step3 output
    step3_path = os.path.join(parameter.data_dir, "step3_df.csv")
    if not os.path.exists(step3_path):
        raise FileNotFoundError(f"Step 3 output not found: {step3_path}. Run 'neoswga score' first.")

    step3_df = pd.read_csv(step3_path)
    candidates = step3_df['primer'].tolist()

    if verbose:
        logger.info(f"Loaded {len(candidates)} candidate primers from step3_df.csv")

    # Initialize position cache for foreground genomes
    # PositionCache(fname_prefixes, primers) loads positions automatically
    if verbose:
        logger.info("Loading position data...")

    cache = PositionCache(fg_prefixes, candidates)

    # Get number of primers to select
    num_primers = getattr(parameter, 'num_primers', 10)
    target_set_size = getattr(parameter, 'target_set_size', num_primers)

    # Create dominating set optimizer
    optimizer = DominatingSetOptimizer(
        cache=cache,
        fg_prefixes=fg_prefixes,
        fg_seq_lengths=fg_seq_lengths,
        bin_size=10000  # 10 kb bins
    )

    # Run greedy optimization (fast, within ln(n) of optimal)
    if verbose:
        logger.info(f"Running dominating set optimization for {target_set_size} primers...")
        logger.info("Algorithm: Greedy set cover (within ln(n) of optimal)")

    result = optimizer.optimize_greedy(
        candidates=candidates,
        max_primers=target_set_size,
        verbose=verbose
    )

    primers = result['primers']
    coverage = result['coverage']

    # Save results to CSV
    output_csv = os.path.join(parameter.data_dir, "step4_improved_df.csv")
    results_df = pd.DataFrame({
        'primer': primers,
        'score': [coverage] * len(primers),
        'set_index': [0] * len(primers),
        'coverage': [coverage] * len(primers),
        'regions_covered': [result['covered_regions']] * len(primers),
        'total_regions': [result['total_regions']] * len(primers)
    })
    results_df.to_csv(output_csv, index=False)

    if verbose:
        logger.info(f"Results saved to {output_csv}")
        logger.info(f"Coverage: {coverage:.1%} ({result['covered_regions']}/{result['total_regions']} regions)")

    # Return in standard format (list of primer sets, list of scores)
    return [primers], [coverage]


if __name__ == "__main__":
    print("Minimum Dominating Set Optimizer")
    print("\nModels primer selection as graph theory problem:")
    print("  - Nodes: Genome regions (bins)")
    print("  - Edges: Primer covers region")
    print("  - Goal: Minimum primers to cover all regions")
    print("\nFeatures:")
    print("  - Greedy set cover (ln(n) approximation)")
    print("  - ILP for exact solution (when feasible)")
    print("  - Provable bounds on solution quality")
    print("\nUsage:")
    print("  from neoswga.core.dominating_set_optimizer import DominatingSetOptimizer")
    print("  optimizer = DominatingSetOptimizer(cache, fg_prefixes, fg_seq_lengths)")
    print("  result = optimizer.optimize_greedy(candidates, max_primers=15)")
    print("  print(f'Coverage: {result[\"coverage\"]:.1%}')")
