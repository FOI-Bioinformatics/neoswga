"""Deterministic, bounded one-for-one refinement of a primer panel."""

import math
import time
from collections import Counter
from dataclasses import dataclass


@dataclass(frozen=True)
class SwapResult:
    primers: tuple[str, ...]
    evaluations: int
    swaps: int
    stop_reason: str
    covered_bases: int
    background_sites: float


def refine_by_swaps(
    selected,
    candidates,
    bins_by_primer,
    bin_weights,
    dimers,
    *,
    fixed_primers=(),
    background_sites=None,
    max_evaluations=10000,
    max_seconds=10.0,
):
    """Improve covered bases, then background load, without increasing panel size.

    Each accepted swap strictly improves this lexicographic objective. New
    primers must be compatible with every retained primer. Coverage counts are
    updated only for bins touched by a swap. Time limits are cooperative,
    checked between evaluations; preprocessing belongs to the caller.
    """
    if max_evaluations < 0 or not math.isfinite(max_seconds) or max_seconds < 0:
        raise ValueError("Swap budgets must be finite and non-negative")
    current = list(dict.fromkeys(selected))
    pool = list(dict.fromkeys(candidates))
    fixed = set(fixed_primers)
    if not fixed.issubset(current):
        raise ValueError("Fixed primers must be present in the initial panel")
    bg = background_sites or {}
    counts = Counter(b for p in current for b in bins_by_primer.get(p, ()))
    covered = sum(bin_weights[b] for b in counts)
    background = sum(bg.get(p, 0.0) for p in current)
    evaluations = swaps = 0
    deadline = time.monotonic() + max_seconds
    reason = "local_optimum"
    while True:
        best = None
        best_gain = (0, 0.0)
        selected_set = set(current)
        exhausted = False
        for incoming in pool:
            if incoming in selected_set:
                continue
            new_bins = bins_by_primer.get(incoming, set())
            for outgoing in current:
                if outgoing in fixed:
                    continue
                if evaluations >= max_evaluations or time.monotonic() >= deadline:
                    reason = "evaluation_limit" if evaluations >= max_evaluations else "time_limit"
                    exhausted = True
                    break
                evaluations += 1
                retained = [p for p in current if p != outgoing]
                if dimers.dimerises(incoming, retained):
                    continue
                old_bins = bins_by_primer.get(outgoing, set())
                gain = sum(bin_weights[b] for b in new_bins if counts[b] == 0)
                loss = sum(bin_weights[b] for b in old_bins - new_bins if counts[b] == 1)
                improvement = (gain - loss, bg.get(outgoing, 0.0) - bg.get(incoming, 0.0))
                if improvement > best_gain:
                    best_gain = improvement
                    best = (outgoing, incoming)
            if exhausted:
                break
        if best is not None:
            outgoing, incoming = best
            for b in bins_by_primer.get(outgoing, ()):
                counts[b] -= 1
            counts.update(bins_by_primer.get(incoming, ()))
            current[current.index(outgoing)] = incoming
            covered += best_gain[0]
            background -= best_gain[1]
            swaps += 1
        if exhausted or best is None:
            break
    return SwapResult(tuple(current), evaluations, swaps, reason, covered, background)


def refine_hybrid_stage2(optimizer, primers, candidates, fixed_primers):
    """Stage-2 swap refinement for :class:`HybridOptimizer`.

    The glue that turns an optimizer's coverage bins into the plain
    ``(bins, weights, compatibility)`` inputs :func:`refine_by_swaps` takes.
    It lives beside the search rather than on the optimizer because it is
    specific to the swap method, and because ``hybrid_optimizer`` is the
    module the size ratchet is closest to.
    """
    import logging

    from .dimer_matrix import build

    logger = logging.getLogger(__name__)

    pool = list(dict.fromkeys(list(primers) + list(candidates)))
    regions = optimizer._coverage_bins_by_primer(pool)
    bins = {p: {optimizer._bin_key(r) for r in owned} for p, owned in regions.items()}
    weights = {optimizer._bin_key(r): r.end - r.start for owned in regions.values() for r in owned}
    background = (
        {p: optimizer._count_background_sites([p]) for p in pool}
        if optimizer.background_pruning and optimizer.bg_prefixes
        else {}
    )
    result = refine_by_swaps(
        primers,
        pool,
        bins,
        weights,
        build(pool, optimizer.max_dimer_bp),
        fixed_primers=fixed_primers,
        background_sites=background,
        max_evaluations=optimizer.swap_max_evaluations,
        max_seconds=optimizer.swap_max_seconds,
    )
    logger.info(
        "Swap refinement: %d swaps, %d evaluations, stopped at %s",
        result.swaps,
        result.evaluations,
        result.stop_reason,
    )
    return list(result.primers)
