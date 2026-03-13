"""
Serial cascade optimizer for multi-stage pipelines.

Chains two or more optimizers sequentially, passing each stage's output
primers as the candidate pool for the next stage.  This enables combining
the strengths of different algorithms — for example, fast coverage
selection followed by exhaustive dimer-free filtering.

Registered pipelines:
    coverage-then-dimerfree   DS -> Clique (coverage + dimer-free guarantee)
    dimerfree-scored          Clique -> Network (dimer-free + connectivity scoring)

Usage:
    neoswga optimize -j params.json --optimization-method coverage-then-dimerfree
"""

import logging
from typing import List, Optional

from .base_optimizer import (
    BaseOptimizer,
    OptimizerConfig,
    OptimizationResult,
    OptimizationStatus,
)
from .optimizer_factory import OptimizerFactory

logger = logging.getLogger(__name__)


class SerialCascadeOptimizer(BaseOptimizer):
    """
    Run optimizers sequentially, passing output as next input.

    Each stage narrows the candidate pool.  The final stage produces
    the returned result.  Intermediate target sizes are scaled so that
    each stage keeps enough candidates for the next.

    Args:
        stages: List of (method_name, target_size_multiplier) tuples.
            The multiplier is applied to the final target_size to compute
            how many primers each intermediate stage should keep.
            The last stage always uses the actual target_size.
    """

    def __init__(
        self,
        stages: List[tuple],
        position_cache,
        fg_prefixes: List[str],
        fg_seq_lengths: List[int],
        bg_prefixes: Optional[List[str]] = None,
        bg_seq_lengths: Optional[List[int]] = None,
        config: Optional[OptimizerConfig] = None,
        **kwargs,
    ):
        super().__init__(
            position_cache, fg_prefixes, fg_seq_lengths,
            bg_prefixes, bg_seq_lengths, config,
        )
        self.stages = stages
        self._extra_kwargs = kwargs

    @property
    def name(self) -> str:
        names = [s[0] for s in self.stages]
        return "serial:" + "->".join(names)

    @property
    def description(self) -> str:
        names = [s[0] for s in self.stages]
        return f"Serial cascade: {' -> '.join(names)}"

    def optimize(
        self,
        candidates: List[str],
        target_size: Optional[int] = None,
        fixed_primers: Optional[List[str]] = None,
        **kwargs,
    ) -> OptimizationResult:
        candidates = self._validate_candidates(candidates)
        target = target_size or self.config.target_set_size
        total_iterations = 0
        messages = []

        current_candidates = list(candidates)

        for stage_idx, (method, multiplier) in enumerate(self.stages):
            is_last = stage_idx == len(self.stages) - 1
            stage_target = target if is_last else max(target, int(target * multiplier))

            if self.config.verbose:
                logger.info(
                    f"  Stage {stage_idx + 1}/{len(self.stages)}: {method} "
                    f"({len(current_candidates)} candidates -> target {stage_target})"
                )

            try:
                sub_optimizer = OptimizerFactory.create(
                    name=method,
                    position_cache=self.cache,
                    fg_prefixes=self.fg_prefixes,
                    fg_seq_lengths=self.fg_seq_lengths,
                    bg_prefixes=self.bg_prefixes,
                    bg_seq_lengths=self.bg_seq_lengths,
                    config=OptimizerConfig(
                        target_set_size=stage_target,
                        max_iterations=self.config.max_iterations,
                        max_dimer_bp=self.config.max_dimer_bp,
                        min_tm=self.config.min_tm,
                        max_tm=self.config.max_tm,
                        verbose=self.config.verbose,
                    ),
                    **self._extra_kwargs,
                )
            except Exception as e:
                msg = f"Stage {stage_idx + 1} ({method}) failed to create: {e}"
                logger.error(msg)
                return OptimizationResult.failure(self.name, msg)

            result = sub_optimizer.optimize(
                current_candidates,
                target_size=stage_target,
                fixed_primers=fixed_primers,
                **kwargs,
            )
            total_iterations += result.iterations

            if not result.is_success:
                msg = f"Stage {stage_idx + 1} ({method}) failed: {result.message}"
                logger.warning(msg)
                return OptimizationResult.failure(self.name, msg)

            messages.append(
                f"{method}: {len(current_candidates)}->{result.num_primers}"
            )
            current_candidates = list(result.primers)

            if len(current_candidates) < target and not is_last:
                logger.warning(
                    f"Stage {stage_idx + 1} produced only {len(current_candidates)} "
                    f"primers (target {target}); skipping remaining stages"
                )
                break

        # Recompute metrics using base class for consistency
        final_primers = list(current_candidates[:target])
        metrics = self.compute_metrics(final_primers)

        return OptimizationResult(
            primers=tuple(final_primers),
            score=result.score,
            status=OptimizationStatus.SUCCESS,
            metrics=metrics,
            iterations=total_iterations,
            optimizer_name=self.name,
            message=f"Serial cascade: {' -> '.join(messages)}",
        )


# ---------------------------------------------------------------------------
# Pre-configured pipeline registrations
# ---------------------------------------------------------------------------

@OptimizerFactory.register(
    'coverage-then-dimerfree',
    aliases=['ds-clique', 'coverage-dimerfree'],
    description='Serial: dominating-set for coverage, then clique for dimer-free guarantee',
)
class CoverageThenDimerFreeOptimizer(SerialCascadeOptimizer):
    """Coverage-first, then dimer-free filtering via DS -> Clique cascade."""

    def __init__(self, position_cache, fg_prefixes, fg_seq_lengths,
                 bg_prefixes=None, bg_seq_lengths=None, config=None, **kwargs):
        stages = [
            ('dominating-set', 3),  # keep ~3x target for clique input
            ('clique', 1),          # final selection
        ]
        super().__init__(
            stages=stages,
            position_cache=position_cache,
            fg_prefixes=fg_prefixes,
            fg_seq_lengths=fg_seq_lengths,
            bg_prefixes=bg_prefixes,
            bg_seq_lengths=bg_seq_lengths,
            config=config,
            **kwargs,
        )


@OptimizerFactory.register(
    'dimerfree-scored',
    aliases=['clique-network'],
    description='Serial: clique for dimer-free sets, then network for connectivity scoring',
)
class DimerFreeScoredOptimizer(SerialCascadeOptimizer):
    """Dimer-free enumeration followed by network connectivity scoring."""

    def __init__(self, position_cache, fg_prefixes, fg_seq_lengths,
                 bg_prefixes=None, bg_seq_lengths=None, config=None, **kwargs):
        stages = [
            ('clique', 2),    # enumerate dimer-free, keep 2x target
            ('network', 1),   # score by connectivity
        ]
        super().__init__(
            stages=stages,
            position_cache=position_cache,
            fg_prefixes=fg_prefixes,
            fg_seq_lengths=fg_seq_lengths,
            bg_prefixes=bg_prefixes,
            bg_seq_lengths=bg_seq_lengths,
            config=config,
            **kwargs,
        )


@OptimizerFactory.register(
    'bg-prefilter-hybrid',
    aliases=['bgpf-hybrid'],
    description='Serial: background pruning pre-filter, then hybrid optimization',
)
class BgPrefilterHybridOptimizer(SerialCascadeOptimizer):
    """Background pruning followed by hybrid (DS + Network) optimization."""

    def __init__(self, position_cache, fg_prefixes, fg_seq_lengths,
                 bg_prefixes=None, bg_seq_lengths=None, config=None, **kwargs):
        stages = [
            ('bg-prefilter', 3),  # prune to ~3x target
            ('hybrid', 1),        # final hybrid optimization
        ]
        super().__init__(
            stages=stages,
            position_cache=position_cache,
            fg_prefixes=fg_prefixes,
            fg_seq_lengths=fg_seq_lengths,
            bg_prefixes=bg_prefixes,
            bg_seq_lengths=bg_seq_lengths,
            config=config,
            **kwargs,
        )
