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

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

import networkx as nx
import numpy as np

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


def _deterministic_scan_order(fixed_primers, candidates, graph_primers) -> List[str]:
    """The order the greedy loop considers primers in.

    `graph.primers` is a set, and greedy selection keeps the first primer seen
    at the best score -- so ties were broken by set iteration order, which
    Python varies per process for strings. Four runs of the bundled example at
    seed=42 gave four different sets: same first primer, different tie-break on
    the second. No RNG is involved, which is why seeding never fixed it.

    Input order rather than `sorted()`: the caller's list arrives ranked by the
    scoring step, so preferring the earlier candidate on a tie respects that
    ranking instead of substituting an alphabetical one.
    """
    seen = set()
    order = []
    for primer in list(fixed_primers) + list(candidates):
        if primer in graph_primers and primer not in seen:
            seen.add(primer)
            order.append(primer)
    return order


def coverage_bin_size(bin_size: int, extension_reach: int) -> int:
    """Bin size that does not overstate what the polymerase reaches.

    `BipartiteGraph.add_primer_coverage` marks a bin covered when any part of
    it is within reach of a binding site. That approximation is only sound
    while a bin is no larger than the reach: at the shipped defaults --
    10 kb bins against a realistic ~3 kb per-primer reach -- a primer covering
    30% of a bin claimed all of it.

    Two things followed, both silent. Greedy set cover saw the genome fully
    covered after two primers and stopped, so every request between 2 and 12
    primers returned the same two; and the coverage reported for that set was
    1.000 where the honest figure was 0.433. A coarse bin rounds its own error
    up to a perfect score, which is the one direction a reader will not
    question.

    A bin equal to the reach is not enough. A site at `pos` reaches
    `[pos, pos + reach)` but marks every bin it touches, so it claims up to
    `reach + bin_size` bases: the relative overstatement is `bin_size / reach`.
    Binning at a quarter of the reach holds that under 25%, which is the
    accuracy this approximation is used at -- progress reporting and the
    background-pruning coverage floor. The authoritative figure remains
    `BaseOptimizer._compute_coverage`, which counts bases rather than bins.

    Callers that model no extension (`extension_reach=0`) keep their bin size:
    with no reach there is no reach to overstate.
    """
    if extension_reach and extension_reach > 0:
        return max(1, min(bin_size, extension_reach // 4 or 1))
    return bin_size


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

    def add_primer_coverage(
        self,
        primer: str,
        positions: np.ndarray,
        genome_id: str,
        genome_length: int,
        extension_reach: int = 0,
        forward_positions: np.ndarray = None,
        reverse_positions: np.ndarray = None,
    ):
        """
        Add primer and its coverage regions.

        When forward_positions and reverse_positions are provided separately,
        uses strand-direction-aware extension: forward sites extend downstream
        only, reverse sites extend upstream only. This models the biological
        reality that phi29 polymerase extends from the 3' end of the primer
        in one direction along the template.

        Args:
            primer: Primer sequence
            positions: Array of binding positions (used when strand info
                is not available, extends bidirectionally)
            genome_id: Genome identifier
            genome_length: Total genome length
            extension_reach: Polymerase extension distance in bp.
            forward_positions: Positions on forward strand (extend downstream)
            reverse_positions: Positions on reverse strand (extend upstream)
        """
        self.primers.add(primer)

        if primer not in self.primer_to_regions:
            self.primer_to_regions[primer] = set()

        n_bins = (genome_length + self.bin_size - 1) // self.bin_size

        covered_bins = set()

        if extension_reach > 0 and (forward_positions is not None or reverse_positions is not None):
            # Strand-direction-aware extension (realistic model)
            if forward_positions is not None:
                for pos in forward_positions:
                    # Forward strand: extend DOWNSTREAM only
                    pos = int(min(pos, genome_length - 1))
                    end_pos = min(genome_length, pos + extension_reach)
                    start_bin = pos // self.bin_size
                    end_bin = min(end_pos // self.bin_size, n_bins - 1)
                    for b in range(start_bin, end_bin + 1):
                        covered_bins.add(b)

            if reverse_positions is not None:
                for pos in reverse_positions:
                    # Reverse strand: extend UPSTREAM only
                    pos = int(min(pos, genome_length))
                    start_pos = max(0, pos - extension_reach)
                    start_bin = start_pos // self.bin_size
                    end_bin = min(pos // self.bin_size, n_bins - 1)
                    for b in range(start_bin, end_bin + 1):
                        covered_bins.add(b)

        elif extension_reach > 0:
            # Bidirectional fallback (when strand info not available)
            for pos in positions:
                start_pos = max(0, int(pos) - extension_reach)
                end_pos = min(genome_length, int(pos) + extension_reach)
                start_bin = start_pos // self.bin_size
                end_bin = min(end_pos // self.bin_size, n_bins - 1)
                for b in range(start_bin, end_bin + 1):
                    covered_bins.add(b)
        else:
            # No extension: original bin-only coverage
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
                    chromosome=genome_id, start=start, end=end, covered_by=set()
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

    def __init__(
        self,
        cache,
        fg_prefixes: List[str],
        fg_seq_lengths: List[int],
        bin_size: int = 10000,
        extension_reach: int = 0,
    ):
        """
        Initialize optimizer.

        Args:
            cache: PositionCache with primer positions
            fg_prefixes: Target genome prefixes
            fg_seq_lengths: Target genome lengths
            bin_size: Size of genome bins for coverage
            extension_reach: Polymerase extension distance in bp. When > 0,
                each binding site covers all bins within this distance.
                Use 70000 for phi29, 80000 for equiphi29.
        """
        self.cache = cache
        self.fg_prefixes = fg_prefixes
        self.fg_seq_lengths = fg_seq_lengths
        self.extension_reach = extension_reach
        self.bin_size = coverage_bin_size(bin_size, extension_reach)
        if self.bin_size != bin_size:
            logger.info(
                f"Bin size reduced {bin_size:,} -> {self.bin_size:,} bp to match the "
                f"{extension_reach:,} bp extension reach"
            )

    def _total_bins(self) -> int:
        """Bins in the genome, not bins some candidate happens to reach.

        `BipartiteGraph.regions` only gains a bin when a candidate covers it, so
        `len(covered) / len(graph.regions)` is the fraction of the *reachable*
        genome -- a number that rises as the pool narrows. A pool touching a
        tenth of the target and covering all of it scored 1.000, which became
        `score` and passed a 0.95 SUCCESS gate. Dividing by this instead makes
        coverage a property of the genome, which is what every consumer of the
        number assumes it already was.
        """
        return sum((length + self.bin_size - 1) // self.bin_size for length in self.fg_seq_lengths)

    def _genome_fraction(self, covered_regions) -> float:
        """Covered bins as a fraction of the genome. See `_total_bins`."""
        total = self._total_bins()
        return len(covered_regions) / total if total else 0.0

    def optimize_greedy(
        self,
        candidates: List[str],
        max_primers: int = 20,
        fixed_primers: Optional[List[str]] = None,
        min_coverage: Optional[float] = None,
        verbose: bool = True,
    ) -> Dict:
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
            min_coverage: Optional minimum coverage fraction (0.0-1.0).
                When set, the algorithm stops adding primers once the
                coverage target is met, even if max_primers is not reached.
            verbose: Print progress

        Returns:
            Dictionary with selected primers and coverage stats.
            Includes 'coverage_target' and 'target_met' when min_coverage is set.
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
                fw = self.cache.get_positions(prefix, primer, "forward")
                rv = self.cache.get_positions(prefix, primer, "reverse")
                both = self.cache.get_positions(prefix, primer, "both")
                if len(both) > 0:
                    graph.add_primer_coverage(
                        primer,
                        both,
                        prefix,
                        length,
                        extension_reach=self.extension_reach,
                        forward_positions=fw if len(fw) > 0 else None,
                        reverse_positions=rv if len(rv) > 0 else None,
                    )

        # Add candidate primers
        for primer in candidates:
            for prefix, length in zip(self.fg_prefixes, self.fg_seq_lengths):
                fw = self.cache.get_positions(prefix, primer, "forward")
                rv = self.cache.get_positions(prefix, primer, "reverse")
                both = self.cache.get_positions(prefix, primer, "both")

                if len(both) > 0:
                    graph.add_primer_coverage(
                        primer,
                        both,
                        prefix,
                        length,
                        extension_reach=self.extension_reach,
                        forward_positions=fw if len(fw) > 0 else None,
                        reverse_positions=rv if len(rv) > 0 else None,
                    )

        total_bins = self._total_bins()
        if verbose:
            logger.info(
                f"Graph: {len(graph.primers)} primers, {len(graph.regions)} regions "
                f"of {total_bins} genome bins"
            )

        # Greedy set cover - start with fixed primers' coverage.
        #
        # `order` records the sequence selection happened in, which `selected`
        # (a set) destroys. It matters downstream: greedy picks the primer with
        # the largest marginal coverage at each step, so its first N picks are
        # the same N whatever `max_primers` was, and coverage over that prefix
        # is non-decreasing in N. That makes the prefix a monotone baseline the
        # hybrid optimizer's Stage-2 refinement is floored against.
        selected = set(fixed_primers)
        order = list(fixed_primers)
        covered_regions = set()

        # Pre-fill coverage from fixed primers
        for primer in fixed_primers:
            covered_regions.update(graph.primer_to_regions.get(primer, set()))

        scan_order = _deterministic_scan_order(fixed_primers, candidates, graph.primers)

        if verbose and n_fixed > 0:
            coverage_so_far = self._genome_fraction(covered_regions)
            logger.info(f"  Fixed primer coverage: {coverage_so_far:.1%}")

        for iteration in range(max_primers):
            if len(covered_regions) == len(graph.regions):
                if verbose:
                    logger.info("Full coverage achieved")
                break

            # Find primer that covers most uncovered regions
            best_primer = None
            best_new_coverage = 0

            for primer in scan_order:
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
            order.append(best_primer)
            covered_regions.update(graph.primer_to_regions[best_primer])

            if verbose and (iteration + 1) % 5 == 0:
                coverage = self._genome_fraction(covered_regions)
                logger.info(f"  {iteration + 1} primers: {coverage:.1%} coverage")

            # Check if coverage target is met
            if min_coverage is not None and graph.regions:
                current_coverage = self._genome_fraction(covered_regions)
                if current_coverage >= min_coverage:
                    if verbose:
                        logger.info(
                            f"Coverage target {min_coverage:.1%} met "
                            f"({current_coverage:.1%}) with "
                            f"{len(selected) - n_fixed} new primers"
                        )
                    break

        # Final statistics
        coverage = self._genome_fraction(covered_regions)
        reachable_coverage = len(covered_regions) / len(graph.regions) if graph.regions else 0.0
        uncovered = graph.regions - covered_regions

        # Both lists come from `order`, not from `selected`.
        #
        # Returning `list(selected)` handed callers a set-ordered list, and
        # Python randomises string hashing per process, so the order differed
        # between runs of the same command. Stage-2 refinement breaks ties by
        # position, so a different order gave a different primer set: the
        # bundled example produced five different answers in five runs at
        # seed=42. No RNG was involved, which is why seeding never helped --
        # greedy selection is deterministic, and only the container forgot the
        # order it happened in.
        fixed_set = set(fixed_primers)
        new_primers = [p for p in order if p not in fixed_set]

        result = {
            "primers": list(order),
            "ordered_primers": order,
            "new_primers": new_primers,
            "fixed_primers": fixed_primers,
            "n_primers": len(selected),
            "n_new": len(new_primers),
            "n_fixed": n_fixed,
            "coverage": coverage,
            # Well above `coverage`: the pool cannot reach the rest of the genome.
            "reachable_coverage": reachable_coverage,
            "covered_regions": len(covered_regions),
            "total_regions": len(graph.regions),
            "uncovered_regions": len(uncovered),
            "graph": graph,
            "coverage_target": min_coverage,
            "target_met": min_coverage is not None and coverage >= min_coverage,
        }

        if verbose:
            if n_fixed > 0:
                logger.info(
                    f"Selected {len(new_primers)} new primers + {n_fixed} fixed = {len(selected)} total"
                )
            else:
                logger.info(f"Selected {len(selected)} primers")
            logger.info(
                f"Coverage: {coverage:.1%} ({len(covered_regions)}/{len(graph.regions)} regions)"
            )

        return result

    def optimize_ilp(
        self, candidates: List[str], max_primers: int = 20, verbose: bool = True
    ) -> Optional[Dict]:
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
            from mip import BINARY, Model, minimize, xsum
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
                positions = self.cache.get_positions(prefix, primer, "both")

                if len(positions) > 0:
                    graph.add_primer_coverage(primer, positions, prefix, length)

        primers_list = list(graph.primers)
        regions_list = list(graph.regions)

        # Create model
        model = Model(sense=minimize)

        # Decision variables: x[i] = 1 if primer i selected
        x = {primer: model.add_var(var_type=BINARY, name=f"x_{primer}") for primer in primers_list}

        # Objective: Minimize number of primers
        model.objective = xsum(x[primer] for primer in primers_list)

        # Constraints: Each region must be covered by at least one primer
        for region in regions_list:
            covering_primers = graph.region_to_primers[region]
            if covering_primers:
                model += xsum(x[primer] for primer in covering_primers if primer in x) >= 1

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

            total_bins = self._total_bins()
            coverage = len(covered_regions) / total_bins if total_bins else 0.0
            reachable = len(covered_regions) / len(graph.regions) if graph.regions else 0.0

            result = {
                "primers": selected,
                "n_primers": len(selected),
                "coverage": coverage,
                "reachable_coverage": reachable,
                "covered_regions": len(covered_regions),
                "total_regions": len(graph.regions),
                "optimal": True,
                "graph": graph,
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
        raise FileNotFoundError(
            f"Step 3 output not found: {step3_path}. Run 'neoswga score' first."
        )

    step3_df = pd.read_csv(step3_path)
    candidates = step3_df["primer"].tolist()

    if verbose:
        logger.info(f"Loaded {len(candidates)} candidate primers from step3_df.csv")

    # Initialize position cache for foreground genomes
    # PositionCache(fname_prefixes, primers) loads positions automatically
    if verbose:
        logger.info("Loading position data...")

    cache = PositionCache(fg_prefixes, candidates)

    # Get number of primers to select
    num_primers = getattr(parameter, "num_primers", 10)
    target_set_size = getattr(parameter, "target_set_size", num_primers)

    # Per-primer reach for coverage. Use realistic phi29 reach (Phase 16,
    # default 3 kb) so this entry's coverage matches the dominating-set
    # adapter and base_optimizer rather than reporting bin-only coverage.
    polymerase = getattr(parameter, "polymerase", "phi29")
    try:
        from neoswga.core.coverage import polymerase_extension_reach

        extension_reach = polymerase_extension_reach(polymerase, coverage_metric="realistic")
    except (ImportError, ValueError):
        extension_reach = 3000

    # Create dominating set optimizer
    optimizer = DominatingSetOptimizer(
        cache=cache,
        fg_prefixes=fg_prefixes,
        fg_seq_lengths=fg_seq_lengths,
        bin_size=10000,  # 10 kb bins
        extension_reach=extension_reach,
    )

    # Run greedy optimization (fast, within ln(n) of optimal)
    if verbose:
        logger.info(f"Running dominating set optimization for {target_set_size} primers...")
        logger.info("Algorithm: Greedy set cover (within ln(n) of optimal)")

    result = optimizer.optimize_greedy(
        candidates=candidates, max_primers=target_set_size, verbose=verbose
    )

    primers = result["primers"]
    coverage = result["coverage"]

    # Save results to CSV
    output_csv = os.path.join(parameter.data_dir, "step4_improved_df.csv")
    results_df = pd.DataFrame(
        {
            "primer": primers,
            "score": [coverage] * len(primers),
            "set_index": [0] * len(primers),
            "coverage": [coverage] * len(primers),
            "regions_covered": [result["covered_regions"]] * len(primers),
            "total_regions": [result["total_regions"]] * len(primers),
        }
    )
    results_df.to_csv(output_csv, index=False)

    if verbose:
        logger.info(f"Results saved to {output_csv}")
        logger.info(
            f"Coverage: {coverage:.1%} ({result['covered_regions']}/{result['total_regions']} regions)"
        )

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
