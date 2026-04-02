"""
Standalone background pruning pre-filter.

Extracted from BackgroundAwareOptimizer Stage 2. Removes primers with
the worst background-to-coverage-loss ratio, keeping only those that
contribute foreground coverage without excessive background binding.

Can be used as a stage in serial cascades or as a standalone filter
before any downstream optimizer.

Usage:
    neoswga optimize -j params.json --optimization-method bg-prefilter-hybrid

Programmatic:
    prefilter = BackgroundPrefilter(cache, fg_prefixes, fg_seq_lengths,
                                    bg_prefixes, bg_seq_lengths)
    result = prefilter.optimize(candidates, target_size=20)
    # result.primers contains the filtered candidates
"""

import logging
import numpy as np
from typing import List, Optional, Tuple

from .base_optimizer import (
    BaseOptimizer,
    OptimizerConfig,
    OptimizationResult,
    OptimizationStatus,
)
from .optimizer_factory import OptimizerFactory

logger = logging.getLogger(__name__)


@OptimizerFactory.register(
    'bg-prefilter',
    aliases=['background-prefilter', 'bg-prune'],
    description='Background pruning pre-filter (removes high-background primers)',
)
class BackgroundPrefilter(BaseOptimizer):
    """Remove primers with poor foreground/background ratio.

    Iteratively drops the primer whose removal causes the least coverage
    loss relative to its background contribution, until the target size
    is reached or coverage would fall below a threshold.

    Designed as a first stage in serial cascades: it narrows a large
    candidate pool to a smaller, background-reduced set that downstream
    optimizers can process more efficiently.
    """

    def __init__(
        self,
        position_cache,
        fg_prefixes: List[str],
        fg_seq_lengths: List[int],
        bg_prefixes: Optional[List[str]] = None,
        bg_seq_lengths: Optional[List[int]] = None,
        config: Optional[OptimizerConfig] = None,
        background_weight: float = 2.0,
        min_coverage: float = 0.50,
        **kwargs,
    ):
        super().__init__(
            position_cache, fg_prefixes, fg_seq_lengths,
            bg_prefixes, bg_seq_lengths, config,
        )
        self.background_weight = background_weight
        self.min_coverage = min_coverage

    @property
    def name(self) -> str:
        return "bg-prefilter"

    @property
    def description(self) -> str:
        return "Background pruning pre-filter"

    @property
    def supports_background(self) -> bool:
        return True

    def optimize(
        self,
        candidates: List[str],
        target_size: Optional[int] = None,
        fixed_primers: Optional[List[str]] = None,
        **kwargs,
    ) -> OptimizationResult:
        candidates = self._validate_candidates(candidates)
        target = target_size or self.config.target_set_size

        if not self.bg_prefixes:
            # No background genome — pass through unchanged
            metrics = self.compute_metrics(candidates[:target])
            return OptimizationResult(
                primers=tuple(candidates[:target]),
                score=metrics.normalized_score(),
                status=OptimizationStatus.SUCCESS,
                metrics=metrics,
                iterations=0,
                optimizer_name=self.name,
                message="No background genome; skipped pruning",
            )

        # When the target is much smaller than the pool, greedy removal
        # is too slow (O(n^2) per removal). Use fast ranking instead.
        if len(candidates) > target * 5:
            pruned = self._rank_and_select(candidates, target)
            final_cov = self._calc_coverage(pruned)
            final_bg = self._count_bg(pruned)
            if self.config.verbose:
                logger.info(
                    f"  BG pre-filter (fast rank): {len(candidates)}->{len(pruned)} "
                    f"(coverage={final_cov:.1%}, bg={final_bg})"
                )
        else:
            pruned, final_cov, final_bg = self._prune_background(
                candidates, target, verbose=self.config.verbose,
            )

        metrics = self.compute_metrics(pruned)

        return OptimizationResult(
            primers=tuple(pruned),
            score=metrics.normalized_score(),
            status=OptimizationStatus.SUCCESS,
            metrics=metrics,
            iterations=len(candidates) - len(pruned),
            optimizer_name=self.name,
            message=(
                f"Pruned {len(candidates)}->{len(pruned)} primers "
                f"(coverage {final_cov:.1%}, bg sites {final_bg})"
            ),
        )

    # ------------------------------------------------------------------
    # Core pruning logic (extracted from BackgroundAwareOptimizer)
    # ------------------------------------------------------------------

    def _prune_background(
        self,
        primers: List[str],
        target_size: int,
        verbose: bool = False,
    ) -> Tuple[List[str], float, int]:
        """Greedy background pruning.

        Repeatedly removes the primer with the highest
        ``background_sites / coverage_loss`` ratio until the pool is
        reduced to *target_size* or coverage would drop below
        *min_coverage*.
        """
        current = list(primers)
        coverage = self._calc_coverage(current)

        if verbose:
            bg = self._count_bg(current)
            logger.info(
                f"  BG pre-filter: {len(current)} -> {target_size} "
                f"(coverage={coverage:.1%}, bg={bg})"
            )

        removed = 0

        while len(current) > target_size:
            best_primer = None
            best_score = -np.inf

            for primer in current:
                test_set = [p for p in current if p != primer]
                test_cov = self._calc_coverage(test_set)
                loss = coverage - test_cov

                if test_cov < self.min_coverage:
                    continue

                bg_sites = self._count_bg([primer])

                if loss > 0:
                    score = bg_sites / loss * self.background_weight
                else:
                    score = bg_sites * 1000.0

                if score > best_score:
                    best_score = score
                    best_primer = primer

            if best_primer is None:
                if verbose:
                    logger.info("  BG pre-filter: stopped (coverage threshold)")
                break

            current.remove(best_primer)
            coverage = self._calc_coverage(current)
            removed += 1

        final_bg = self._count_bg(current)
        if verbose:
            logger.info(
                f"  BG pre-filter done: removed {removed}, "
                f"coverage={coverage:.1%}, bg={final_bg}"
            )

        return current, coverage, final_bg

    def _rank_and_select(self, primers: List[str], target: int) -> List[str]:
        """Fast selection by ranking primers on fg/bg selectivity ratio.

        Used when the pool is much larger than target, where greedy removal
        would be too slow. Ranks by (fg_sites - bg_weight * bg_sites) and
        takes the top N.
        """
        scored = []
        for primer in primers:
            fg = sum(len(self.cache.get_positions(p, primer, 'both'))
                     for p in self.fg_prefixes)
            bg = sum(len(self.cache.get_positions(p, primer, 'both'))
                     for p in self.bg_prefixes)
            rank_score = fg - self.background_weight * bg
            scored.append((rank_score, primer))
        scored.sort(reverse=True)
        return [p for _, p in scored[:target]]

    def _calc_coverage(self, primers: List[str]) -> float:
        if not primers:
            return 0.0
        positions = set()
        total = sum(self.fg_seq_lengths)
        for primer in primers:
            for prefix in self.fg_prefixes:
                for pos in self.cache.get_positions(prefix, primer, 'both'):
                    positions.add(pos)
        return len(positions) / total if total > 0 else 0.0

    def _count_bg(self, primers: List[str]) -> int:
        total = 0
        for primer in primers:
            for prefix in self.bg_prefixes:
                total += len(self.cache.get_positions(prefix, primer, 'both'))
        return total
