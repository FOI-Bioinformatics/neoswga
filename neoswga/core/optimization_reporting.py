"""What a finished optimization tells the user at the console.

Extracted from `unified_optimizer.py` on 2026-09-22. That module sat at 1648 of
its 1650-line ceiling once the resolved request arrived, and `run_optimization`
sat one line over its own 470-line budget -- a budget met exactly is the shape
`tests/test_a_ceiling_needs_headroom.py` exists to refuse, because the next
correct small edit then fails on the merge rather than on the change.

The subject is its own. Reporting a delivered panel is not deciding anything
about it: this module reads the result and the validation record and writes
lines, and it must stay that way. `result_validation` owns the verdict,
`panel_evaluation` owns the measurements, and `design_result` owns the rule
that turns findings into a refusal.
"""

from __future__ import annotations

import logging

from .base_optimizer import OptimizationStatus

logger = logging.getLogger(__name__)

__all__ = ["log_optimization_outcome"]


def log_optimization_outcome(result, target_size, method, validation) -> None:
    """Report the delivered panel and every post-optimization finding.

    `validation` may be None: a caller driving the optimizer programmatically
    does not always produce one, and an absent record is reported as nothing
    rather than as no findings.
    """
    if result.is_success:
        logger.info(f"Selected {result.num_primers} primers")
        logger.info(f"Coverage: {result.metrics.fg_coverage:.1%} (measured)")
        logger.info(f"Score: {result.score:.4f}")
    elif result.status == OptimizationStatus.PARTIAL:
        if result.num_primers < target_size:
            logger.warning(
                f"PARTIAL result: found {result.num_primers} primers "
                f"but target was {target_size}"
            )
        else:
            logger.warning(
                f"PARTIAL result: found {result.num_primers} primers "
                f"but genome coverage is below threshold"
            )
        logger.warning("Suggestions to improve results:")
        logger.warning("  - Relax filtering thresholds (increase max_bg_freq or max_gini)")
        logger.warning("  - Widen the k-mer range (decrease min_k or increase max_k)")
        logger.warning("  - Increase the candidate pool (raise max_primer in filter step)")
        logger.warning(f"  - Try a different optimizer (current: {method})")
    else:
        logger.warning(f"Optimization failed: {result.message}")

    if validation is not None and validation.get("issues"):
        for issue in validation["issues"]:
            level = issue.get("level", "warning")
            message = f"Post-opt {level}: {issue.get('code')} - {issue.get('detail')}"
            (logger.error if level == "error" else logger.warning)(message)
