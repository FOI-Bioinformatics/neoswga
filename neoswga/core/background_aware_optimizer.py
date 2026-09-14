#!/usr/bin/env python3
"""
Background-Aware Hybrid Optimizer for SWGA with Host Genome Suppression.

THE ULTIMATE SOLUTION for "longer oligos with fewer host genome hits":

Problem: Standard optimizers maximize coverage + network but DON'T minimize background.
Solution: Multi-objective optimization with EXPLICIT background minimization.

Three-stage optimization:
1. Coverage (dominating set) - maximize target genome coverage
2. Background Pruning - MINIMIZE background binding (NEW!)
3. Network Refinement - maximize amplification connectivity

Expected Impact:
- 10-20x reduction in background binding (vs standard hybrid optimizer)
- Holds coverage above `min_coverage_threshold` during pruning (default 0.95,
  and a setting the caller can lower)
- Maintains good network connectivity

Critical for:
- 16-18bp primers (more background binding opportunities)
- Human/mouse host backgrounds
- Challenging targets (endosymbionts, parasites, extreme-GC genomes)

Author: NeoSWGA Development Team
Date: November 2025
Version: 3.0 - Phase 2 Advanced (Background-Aware)
"""

import logging
import time
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from neoswga.core.dominating_set_optimizer import DominatingSetOptimizer
from neoswga.core.hybrid_optimizer import HybridOptimizer, HybridResult
from neoswga.core.network_optimizer import NetworkOptimizer

logger = logging.getLogger(__name__)


# The standalone three-stage `BackgroundAwareOptimizer`, its `BackgroundAwareResult`
# and the module-level `optimize()` / `compare_optimizers()` helpers were removed
# on 2026-09-10. Nothing dispatched to them: `OptimizerFactory.register` decorates
# `BackgroundAwareBaseOptimizer` below, which delegates to `HybridOptimizer`, and
# `unified_optimizer` imports this module only to trigger that registration.
#
# They were not merely unreachable. The removed `_prune_background` had DIVERGED
# from the one in `hybrid_optimizer.py`: it kept the absolute coverage floor that
# `hybrid_optimizer._prune_background` records as having made the whole stage a
# no-op, plus a different scoring expression and `strand="both"` counting where
# the live one counts per strand. Keeping a second three-stage optimizer that no
# dispatch path reaches, drifting away from the one that ships, is a maintenance
# cost with no user.
#
# `tests/test_background_aware_live_path.py` is the record of why this is gone.


# =============================================================================
# Factory Registration - BaseOptimizer Interface
# =============================================================================

from neoswga.core.base_optimizer import (
    BaseOptimizer,
    OptimizationResult,
    OptimizationStatus,
    OptimizerConfig,
    PrimerSetMetrics,
)
from neoswga.core.optimizer_factory import OptimizerFactory


@OptimizerFactory.register("background-aware", aliases=["clinical", "bg-aware"])
class BackgroundAwareBaseOptimizer(BaseOptimizer):
    """
    Background-aware optimizer implementing BaseOptimizer interface.

    Delegates to HybridOptimizer with background_pruning=True, providing
    three-stage optimization (coverage + background pruning + network
    refinement) for clinical applications requiring low background binding.
    The final network stage honours ReactionConditions additive corrections
    so this wrapper reports as additive-aware.
    """

    ADDITIVE_AWARE = True

    def __init__(
        self,
        position_cache,
        fg_prefixes: List[str],
        fg_seq_lengths: List[int],
        bg_prefixes: Optional[List[str]] = None,
        bg_seq_lengths: Optional[List[int]] = None,
        config: Optional[OptimizerConfig] = None,
        conditions=None,
        **kwargs,
    ):
        super().__init__(
            position_cache,
            fg_prefixes,
            fg_seq_lengths,
            bg_prefixes,
            bg_seq_lengths,
            config,
            conditions=conditions,
            # Forward the rest so background_profile / aggregate reach
            # BaseOptimizer. Dropping **kwargs here made a compositional
            # background silently inert on every optimizer.
            **kwargs,
        )

        # Delegate to HybridOptimizer with background pruning enabled
        from neoswga.core.hybrid_optimizer import HybridOptimizer

        self._hybrid = HybridOptimizer(
            position_cache=position_cache,
            fg_prefixes=fg_prefixes,
            fg_seq_lengths=fg_seq_lengths,
            bg_prefixes=bg_prefixes or [],
            bg_seq_lengths=bg_seq_lengths or [],
            background_pruning=True,
            # Realistic Stage-1 coverage reach (resolved by unified_optimizer),
            # matching the scored metric and this optimizer's own 3 kb docstring.
            coverage_reach=getattr(self.config, "extension_reach", None),
            background_weight=kwargs.get("background_weight", 2.0),
            min_coverage_threshold=kwargs.get("min_coverage_threshold", 0.95),
            min_coverage_drop=kwargs.get("min_coverage_drop", 0.05),
            absolute_coverage_floor=kwargs.get("absolute_coverage_floor", False),
            polymerase=kwargs.get("polymerase", "phi29"),
            # Forward ReactionConditions + mechanistic weight so the clinical
            # background-aware optimizer's stage-3 network refinement honours
            # DMSO / betaine / etc. Review I1: previously these dropped at
            # this seam even though self.conditions was set on the wrapper.
            conditions=conditions,
            mechanistic_weight=kwargs.get("mechanistic_weight", 0.0),
            # This wrapper builds the same HybridOptimizer class as the `hybrid`
            # wrapper does, and the two argument lists had drifted apart: six
            # values arrived at `hybrid` and were dropped here. The sharpest was
            # the Tm window, which the polymerase preset then replaced -- on an
            # equiphi29 pool configured for 20-70 C, `hybrid` screened at
            # 20.0-70.0 and passed 500/500 candidates while `background-aware`
            # screened at 37.0-62.0 and passed 167/500. Two thirds of a
            # user-approved pool, discarded before selection, by the method
            # documented as the clinical one.
            bin_size=kwargs.get("bin_size", 10000),
            max_extension=kwargs.get("max_extension"),
            min_tm=getattr(self.config, "min_tm", None),
            max_tm=getattr(self.config, "max_tm", None),
            genome_gc_content=kwargs.get("genome_gc_content"),
            uniformity_weight=kwargs.get("uniformity_weight", 0.0),
            reaction_temp=kwargs.get("reaction_temp"),
            tm_weight=kwargs.get("tm_weight", 0.0),
            dimer_penalty=kwargs.get("dimer_penalty", 0.0),
            max_dimer_bp=getattr(self.config, "max_dimer_bp", 4),
            allow_dimer_relaxation=self.config.allow_dimer_relaxation,
            refinement_method=self.config.refinement_method,
            swap_max_evaluations=self.config.swap_max_evaluations,
            swap_max_seconds=self.config.swap_max_seconds,
            template_gc=kwargs.get("template_gc", 0.5),
        )

        # A second BackgroundAwareOptimizer was built here "for backward
        # compat (standalone use)" and read by nothing. Constructing it built a
        # DominatingSetOptimizer and a NetworkOptimizer, so every run of this
        # method paid for two optimizers it never called. The standalone class
        # remains importable for callers that want it directly.

    @property
    def name(self) -> str:
        """Optimizer identifier for logging and factory registration."""
        return "background-aware"

    @property
    def description(self) -> str:
        """One-line summary of the optimization strategy."""
        return "Three-stage optimizer with background reduction for clinical use"

    @property
    def supports_background(self) -> bool:
        """Indicates this optimizer uses background genome data."""
        return True

    def optimize(
        self, candidates: List[str], target_size: Optional[int] = None, **kwargs
    ) -> OptimizationResult:
        """Run background-aware optimization via HybridOptimizer."""
        if not self.bg_prefixes:
            return OptimizationResult.failure(
                self.name, "Background-aware optimizer requires background genome data"
            )

        candidates = self._validate_candidates(candidates)
        target = target_size or self.config.target_set_size

        if self.config.verbose:
            logger.info(f"Running background-aware optimization: {len(candidates)} candidates")

        try:
            result = self._hybrid.optimize(
                candidates=candidates,
                final_count=target,
                verbose=self.config.verbose,
                apply_polymerase_multiplier=False,
            )

            primers = result.primers
            metrics = self.compute_metrics(primers)

            return OptimizationResult(
                primers=tuple(primers),
                score=result.final_predicted_amplification,
                status=(
                    OptimizationStatus.NO_CONVERGENCE
                    if not primers
                    else (
                        OptimizationStatus.PARTIAL
                        if len(primers) < target
                        else OptimizationStatus.SUCCESS
                    )
                ),
                metrics=metrics,
                iterations=3,  # Three stages
                optimizer_name=self.name,
                message=f"Coverage: {result.final_coverage:.1%}, "
                f"Connectivity: {result.final_connectivity:.2f}",
            )

        except Exception as e:
            logger.error(f"Background-aware optimization failed: {e}")
            return OptimizationResult.failure(self.name, str(e))
