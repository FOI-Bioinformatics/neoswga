"""
Primer set expansion for iterative SWGA design.

Allows users to provide existing primers that have been experimentally
validated and design additional primers to fill coverage gaps.

This supports iterative wet-lab workflows where:
1. Initial primer set is synthesized and tested
2. Some primers work well, others fail
3. User wants to design additional primers to improve coverage

Usage:
    from neoswga.core.primer_expansion import PrimerExpander

    expander = PrimerExpander(
        position_cache=cache,
        fg_prefixes=fg_prefixes,
        fg_seq_lengths=fg_seq_lengths,
    )

    # Identify coverage gaps in current primer set
    gaps = expander.identify_gaps(
        primers=['ATCGATCG', 'GCTAGCTA'],
        min_gap_size=10000
    )

    # Expand primer set with additional primers
    result = expander.expand(
        candidates=candidates_df['primer'].tolist(),
        fixed_primers=['ATCGATCG', 'GCTAGCTA'],
        failed_primers=['BADPRIMER'],  # Exclude these
        target_new=6,
        optimization_method='hybrid'
    )
"""

import logging
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple, Set
import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class CoverageGap:
    """A region of the genome not covered by current primers."""
    chromosome: str
    start: int
    end: int
    size: int

    def __post_init__(self):
        if self.size == 0:
            self.size = self.end - self.start


@dataclass
class ExpansionInput:
    """
    Input for primer set expansion.

    Attributes:
        fixed_primers: Primers that must be kept (already validated)
        failed_primers: Primers that failed and should be excluded
        coverage_gaps: Optional list of specific regions needing coverage
    """
    fixed_primers: List[str] = field(default_factory=list)
    failed_primers: List[str] = field(default_factory=list)
    coverage_gaps: List[CoverageGap] = field(default_factory=list)


@dataclass
class ExpansionResult:
    """
    Result of primer set expansion.

    Attributes:
        new_primers: Newly designed primers
        combined_set: Fixed + new primers
        fixed_primers: Original fixed primers
        coverage_before: Coverage with fixed primers only
        coverage_after: Coverage with combined set
        gap_coverage: Fraction of gap regions now covered
        predicted_improvement: Estimated improvement in amplification
        gaps_remaining: Number of significant gaps still present
        optimization_method: Method used for optimization
    """
    new_primers: List[str]
    combined_set: List[str]
    fixed_primers: List[str]
    coverage_before: float
    coverage_after: float
    gap_coverage: float
    predicted_improvement: float
    gaps_remaining: int
    optimization_method: str
    message: str = ""

    @property
    def n_new(self) -> int:
        """Number of new primers designed."""
        return len(self.new_primers)

    @property
    def n_fixed(self) -> int:
        """Number of fixed primers."""
        return len(self.fixed_primers)

    @property
    def n_total(self) -> int:
        """Total primers in combined set."""
        return len(self.combined_set)

    def to_dict(self) -> Dict:
        """Convert to dictionary for JSON serialization."""
        return {
            'new_primers': self.new_primers,
            'combined_set': self.combined_set,
            'fixed_primers': self.fixed_primers,
            'n_new': self.n_new,
            'n_fixed': self.n_fixed,
            'n_total': self.n_total,
            'coverage_before': self.coverage_before,
            'coverage_after': self.coverage_after,
            'gap_coverage': self.gap_coverage,
            'predicted_improvement': self.predicted_improvement,
            'gaps_remaining': self.gaps_remaining,
            'optimization_method': self.optimization_method,
            'message': self.message,
        }

    def __str__(self) -> str:
        lines = [
            "Primer Expansion Result:",
            f"  Fixed primers: {self.n_fixed}",
            f"  New primers: {self.n_new}",
            f"  Total set: {self.n_total}",
            f"  Coverage: {self.coverage_before:.1%} -> {self.coverage_after:.1%}",
            f"  Gap coverage: {self.gap_coverage:.1%}",
            f"  Predicted improvement: {self.predicted_improvement:.1f}x",
            f"  Gaps remaining: {self.gaps_remaining}",
        ]
        if self.message:
            lines.append(f"  Note: {self.message}")
        return "\n".join(lines)


class PrimerExpander:
    """
    Expands existing primer sets with additional primers.

    Supports iterative design workflows where users have validated
    some primers experimentally and want to add more.
    """

    def __init__(
        self,
        position_cache,
        fg_prefixes: List[str],
        fg_seq_lengths: List[int],
        bg_prefixes: Optional[List[str]] = None,
        bg_seq_lengths: Optional[List[int]] = None,
        bin_size: int = 10000,
        max_extension: int = 70000,
    ):
        """
        Initialize primer expander.

        Args:
            position_cache: PositionCache with primer binding positions
            fg_prefixes: Target genome identifiers
            fg_seq_lengths: Target genome lengths
            bg_prefixes: Background genome identifiers (optional)
            bg_seq_lengths: Background genome lengths (optional)
            bin_size: Bin size for coverage analysis (bp)
            max_extension: Maximum polymerase extension distance (bp)
        """
        self.cache = position_cache
        self.fg_prefixes = fg_prefixes
        self.fg_seq_lengths = fg_seq_lengths
        self.bg_prefixes = bg_prefixes or []
        self.bg_seq_lengths = bg_seq_lengths or []
        self.bin_size = bin_size
        self.max_extension = max_extension
        self.total_length = sum(fg_seq_lengths)

    def identify_gaps(
        self,
        primers: List[str],
        min_gap_size: int = 10000,
    ) -> List[CoverageGap]:
        """
        Identify coverage gaps in a primer set.

        Args:
            primers: Current primer sequences
            min_gap_size: Minimum gap size to report (bp)

        Returns:
            List of CoverageGap objects sorted by size (largest first)
        """
        if not primers:
            # Entire genome is a gap
            gaps = []
            for prefix, length in zip(self.fg_prefixes, self.fg_seq_lengths):
                gaps.append(CoverageGap(
                    chromosome=prefix,
                    start=0,
                    end=length,
                    size=length
                ))
            return gaps

        # Collect all binding positions
        all_positions = {}  # prefix -> sorted positions

        for prefix, length in zip(self.fg_prefixes, self.fg_seq_lengths):
            positions = []
            for primer in primers:
                pos = self.cache.get_positions(prefix, primer, 'both')
                positions.extend(pos.tolist())
            all_positions[prefix] = sorted(set(positions))

        # Find gaps
        gaps = []

        for prefix, length in zip(self.fg_prefixes, self.fg_seq_lengths):
            positions = all_positions.get(prefix, [])

            if not positions:
                # Entire chromosome is a gap
                gaps.append(CoverageGap(
                    chromosome=prefix,
                    start=0,
                    end=length,
                    size=length
                ))
                continue

            # Check gap at start
            if positions[0] > min_gap_size:
                gaps.append(CoverageGap(
                    chromosome=prefix,
                    start=0,
                    end=positions[0],
                    size=positions[0]
                ))

            # Check internal gaps
            for i in range(1, len(positions)):
                gap_size = positions[i] - positions[i-1]
                if gap_size > min_gap_size:
                    gaps.append(CoverageGap(
                        chromosome=prefix,
                        start=positions[i-1],
                        end=positions[i],
                        size=gap_size
                    ))

            # Check gap at end
            if length - positions[-1] > min_gap_size:
                gaps.append(CoverageGap(
                    chromosome=prefix,
                    start=positions[-1],
                    end=length,
                    size=length - positions[-1]
                ))

        # Sort by size (largest first)
        gaps.sort(key=lambda g: g.size, reverse=True)

        return gaps

    def _calculate_coverage(self, primers: List[str]) -> float:
        """Calculate genome coverage fraction for a primer set."""
        if not primers:
            return 0.0

        from neoswga.core.dominating_set_optimizer import BipartiteGraph

        graph = BipartiteGraph(bin_size=self.bin_size)

        for primer in primers:
            for prefix, length in zip(self.fg_prefixes, self.fg_seq_lengths):
                positions = self.cache.get_positions(prefix, primer, 'both')
                if len(positions) > 0:
                    graph.add_primer_coverage(primer, positions, prefix, length)

        if not graph.regions:
            return 0.0

        total_bins = sum(
            (length + self.bin_size - 1) // self.bin_size
            for length in self.fg_seq_lengths
        )

        return len(graph.regions) / total_bins if total_bins > 0 else 0.0

    def expand(
        self,
        candidates: List[str],
        fixed_primers: List[str],
        failed_primers: Optional[List[str]] = None,
        target_new: int = 6,
        optimization_method: str = 'hybrid',
        verbose: bool = True,
    ) -> ExpansionResult:
        """
        Expand primer set with additional primers.

        Args:
            candidates: Pool of candidate primers to choose from
            fixed_primers: Primers that must be included (validated)
            failed_primers: Primers to exclude (failed in wet lab)
            target_new: Number of new primers to select
            optimization_method: 'hybrid', 'dominating-set', or 'network'
            verbose: Print progress

        Returns:
            ExpansionResult with new and combined primer sets
        """
        failed_primers = failed_primers or []
        fixed_set = set(p.upper() for p in fixed_primers)
        failed_set = set(p.upper() for p in failed_primers)

        # Filter candidates
        candidates_filtered = [
            c for c in candidates
            if c.upper() not in fixed_set and c.upper() not in failed_set
        ]

        if verbose:
            logger.info("="*60)
            logger.info("PRIMER SET EXPANSION")
            logger.info("="*60)
            logger.info(f"Fixed primers: {len(fixed_primers)}")
            logger.info(f"Failed primers (excluded): {len(failed_primers)}")
            logger.info(f"Candidate pool: {len(candidates_filtered)}")
            logger.info(f"Target new primers: {target_new}")
            logger.info(f"Optimization method: {optimization_method}")

        # Calculate initial coverage
        coverage_before = self._calculate_coverage(list(fixed_set))

        if verbose:
            logger.info(f"\nCurrent coverage: {coverage_before:.1%}")

        # Identify gaps before expansion
        gaps_before = self.identify_gaps(list(fixed_set), min_gap_size=self.bin_size)

        if verbose:
            logger.info(f"Significant gaps: {len(gaps_before)}")
            if gaps_before and len(gaps_before) <= 5:
                for gap in gaps_before[:5]:
                    logger.info(f"  {gap.chromosome}: {gap.start:,}-{gap.end:,} ({gap.size:,} bp)")

        # Select optimization method
        if optimization_method in ('hybrid', 'two-stage'):
            result = self._expand_hybrid(
                candidates_filtered, fixed_primers, target_new, verbose
            )
        elif optimization_method in ('dominating-set', 'dominating_set', 'ds'):
            result = self._expand_dominating_set(
                candidates_filtered, fixed_primers, target_new, verbose
            )
        else:
            # Default to hybrid
            logger.warning(f"Unknown method '{optimization_method}', using hybrid")
            result = self._expand_hybrid(
                candidates_filtered, fixed_primers, target_new, verbose
            )

        # Calculate coverage improvement
        combined_set = list(fixed_set) + result['new_primers']
        coverage_after = self._calculate_coverage(combined_set)

        # Count remaining gaps
        gaps_after = self.identify_gaps(combined_set, min_gap_size=self.bin_size)

        # Calculate gap coverage (fraction of original gaps now covered)
        gap_coverage = 1.0 - (len(gaps_after) / max(len(gaps_before), 1))

        # Estimate improvement factor
        if coverage_before > 0:
            predicted_improvement = coverage_after / coverage_before
        else:
            predicted_improvement = float('inf') if coverage_after > 0 else 1.0

        if verbose:
            logger.info("\n" + "="*60)
            logger.info("EXPANSION COMPLETE")
            logger.info("="*60)
            logger.info(f"New primers selected: {len(result['new_primers'])}")
            logger.info(f"Coverage: {coverage_before:.1%} -> {coverage_after:.1%}")
            logger.info(f"Improvement: {predicted_improvement:.2f}x")
            logger.info(f"Gaps remaining: {len(gaps_after)}")

        return ExpansionResult(
            new_primers=result['new_primers'],
            combined_set=combined_set,
            fixed_primers=fixed_primers,
            coverage_before=coverage_before,
            coverage_after=coverage_after,
            gap_coverage=gap_coverage,
            predicted_improvement=predicted_improvement,
            gaps_remaining=len(gaps_after),
            optimization_method=optimization_method,
            message=result.get('message', ''),
        )

    def _expand_hybrid(
        self,
        candidates: List[str],
        fixed_primers: List[str],
        target_new: int,
        verbose: bool,
    ) -> Dict:
        """Use hybrid optimizer for expansion."""
        from neoswga.core.hybrid_optimizer import HybridOptimizer

        optimizer = HybridOptimizer(
            position_cache=self.cache,
            fg_prefixes=self.fg_prefixes,
            fg_seq_lengths=self.fg_seq_lengths,
            bg_prefixes=self.bg_prefixes,
            bg_seq_lengths=self.bg_seq_lengths,
            bin_size=self.bin_size,
            max_extension=self.max_extension,
        )

        # Target total = fixed + new
        total_target = len(fixed_primers) + target_new

        result = optimizer.optimize(
            candidates=candidates,
            final_count=total_target,
            fixed_primers=fixed_primers,
            verbose=verbose,
        )

        # Extract new primers (those not in fixed)
        fixed_set = set(fixed_primers)
        new_primers = [p for p in result.primers if p not in fixed_set]

        return {
            'new_primers': new_primers,
            'message': f"Connectivity: {result.final_connectivity:.2f}",
        }

    def _expand_dominating_set(
        self,
        candidates: List[str],
        fixed_primers: List[str],
        target_new: int,
        verbose: bool,
    ) -> Dict:
        """Use dominating set optimizer for expansion."""
        from neoswga.core.dominating_set_optimizer import DominatingSetOptimizer

        optimizer = DominatingSetOptimizer(
            cache=self.cache,
            fg_prefixes=self.fg_prefixes,
            fg_seq_lengths=self.fg_seq_lengths,
            bin_size=self.bin_size,
        )

        result = optimizer.optimize_greedy(
            candidates=candidates,
            max_primers=target_new,
            fixed_primers=fixed_primers,
            verbose=verbose,
        )

        return {
            'new_primers': result['new_primers'],
            'message': f"Coverage: {result['coverage']:.1%}",
        }


def expand_primers(
    params_path: str,
    fixed_primers: List[str],
    failed_primers: Optional[List[str]] = None,
    target_new: int = 6,
    optimization_method: str = 'hybrid',
    output_dir: Optional[str] = None,
    verbose: bool = True,
) -> ExpansionResult:
    """
    Convenience function to expand primer set using params.json.

    Args:
        params_path: Path to params.json
        fixed_primers: Primers to keep (validated)
        failed_primers: Primers to exclude (failed)
        target_new: Number of new primers to add
        optimization_method: Optimization method to use
        output_dir: Directory to save results (optional)
        verbose: Print progress

    Returns:
        ExpansionResult with new and combined primer sets
    """
    import os
    import json
    import pandas as pd
    from neoswga.core.position_cache import PositionCache
    from neoswga.core import parameter

    # Load parameters
    with open(params_path) as f:
        params = json.load(f)

    # Get prefixes and lengths
    fg_prefixes = params.get('fg_prefixes', [])
    bg_prefixes = params.get('bg_prefixes', [])
    data_dir = params.get('data_dir', './')

    # Load genome lengths (from step2 or calculate)
    # This is a simplified version - in practice, get from genome files
    fg_seq_lengths = params.get('fg_seq_lengths', [])
    bg_seq_lengths = params.get('bg_seq_lengths', [])

    # Load candidates from step3
    step3_path = os.path.join(data_dir, 'step3_df.csv')
    if not os.path.exists(step3_path):
        raise FileNotFoundError(
            f"Step 3 output not found: {step3_path}. Run 'neoswga score' first."
        )

    step3_df = pd.read_csv(step3_path)
    candidates = step3_df['primer'].tolist()

    # Initialize position cache
    cache = PositionCache(fg_prefixes, candidates + fixed_primers)

    # Create expander
    expander = PrimerExpander(
        position_cache=cache,
        fg_prefixes=fg_prefixes,
        fg_seq_lengths=fg_seq_lengths,
        bg_prefixes=bg_prefixes,
        bg_seq_lengths=bg_seq_lengths,
    )

    # Run expansion
    result = expander.expand(
        candidates=candidates,
        fixed_primers=fixed_primers,
        failed_primers=failed_primers,
        target_new=target_new,
        optimization_method=optimization_method,
        verbose=verbose,
    )

    # Save results if output_dir specified
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

        # Save result as JSON
        result_path = os.path.join(output_dir, 'expansion_result.json')
        with open(result_path, 'w') as f:
            json.dump(result.to_dict(), f, indent=2)

        # Save combined primers as CSV
        primers_path = os.path.join(output_dir, 'expanded_primers.csv')
        pd.DataFrame({
            'primer': result.combined_set,
            'type': ['fixed'] * result.n_fixed + ['new'] * result.n_new,
        }).to_csv(primers_path, index=False)

        if verbose:
            logger.info(f"\nResults saved to {output_dir}")

    return result
