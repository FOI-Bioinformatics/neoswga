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
    # `None` when the caller supplied no coverage bins, which is the
    # objective-only mode: the panel's coverage is the objective's to report
    # and inventing a base count here would look like a measurement.
    covered_bases: int | None
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
    objective=None,
    max_evaluations=10000,
    max_seconds=10.0,
):
    """Improve the panel without increasing its size.

    Without `objective`, the rule is the original lexicographic one: raw covered
    bases first, then background load. Existing callers are unchanged.

    With `objective`, acceptance uses the metric the design is judged on, and
    constraints come first:

        (fewer violations, higher coverage, lower background load)

    The ordering matters in both directions. A feasible panel is never swapped
    for an infeasible one however much coverage that would buy, because a
    constraint is not a scoring term to be outbid. And before a panel is
    feasible the useful direction is out of violation rather than up the
    coverage curve, so a swap that reduces violations is taken even when
    coverage falls.

    `bins_by_primer` and `bin_weights` may both be `None` when an `objective` is
    supplied. The objective answers every question acceptance asks, so a caller
    that has one need not also assemble a bin decomposition it will not read;
    `covered_bases` is then `None` rather than a zero that would read as a
    measured value.

    New primers must be compatible with every retained primer. Time limits are
    cooperative, checked between evaluations; preprocessing belongs to the
    caller.
    """
    if max_evaluations < 0 or not math.isfinite(max_seconds) or max_seconds < 0:
        raise ValueError("Swap budgets must be finite and non-negative")
    binless = bins_by_primer is None or bin_weights is None
    if binless and objective is None:
        raise ValueError("Swap refinement needs either coverage bins or an objective")
    bins_by_primer = {} if bins_by_primer is None else bins_by_primer
    bin_weights = {} if bin_weights is None else bin_weights
    current = list(dict.fromkeys(selected))
    pool = list(dict.fromkeys(candidates))
    fixed = set(fixed_primers)
    if not fixed.issubset(current):
        raise ValueError("Fixed primers must be present in the initial panel")
    bg = background_sites or {}
    counts = Counter(b for p in current for b in bins_by_primer.get(p, ()))
    covered = None if binless else sum(bin_weights[b] for b in counts)
    background = sum(bg.get(p, 0.0) for p in current)
    evaluations = swaps = 0
    deadline = time.monotonic() + max_seconds
    reason = "local_optimum"

    def _score(panel):
        """Lexicographic, constraints first. Higher is better throughout."""
        return (
            -len(objective.violations(panel)),
            objective.coverage(panel),
            -objective.metrics(panel).total_bg_sites,
        )

    while True:
        best = None
        best_score = None
        incumbent = _score(current) if objective is not None else None
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
                if objective is not None:
                    candidate_panel = [*retained, incoming]
                    candidate_score = _score(candidate_panel)
                    if candidate_score > incumbent and (
                        best is None or candidate_score > best_score
                    ):
                        best_score = candidate_score
                        best = (outgoing, incoming)
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
            if objective is None:
                covered += best_gain[0]
                background -= best_gain[1]
            else:
                covered = None if binless else sum(bin_weights[b] for b in counts if counts[b])
                background = objective.metrics(current).total_bg_sites
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
    # The objective the caller accepts the result on, when there is one. It
    # arrives as an attribute rather than a parameter because this function
    # already takes the optimizer, so nothing public changes and a plain
    # `optimize` run -- which declares no objective and no constraints -- keeps
    # the raw-bin rule it has always used. Only a caller that states what it is
    # judging on gets refined on it.
    #
    # Stage 2 chooses the delivered panel. Leaving it on covered bases while
    # `plan_pool` accepted on occupancy-weighted coverage meant the stage that
    # picked the panel used the rule the report does not quote, and the two
    # disagree in a known direction: a primer with many sites and a melting
    # temperature well below the reaction temperature touches many bins and
    # contributes little amplification.
    objective = getattr(optimizer, "pool_objective", None)
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
        objective=objective,
        max_evaluations=optimizer.swap_max_evaluations,
        max_seconds=optimizer.swap_max_seconds,
    )
    logger.info(
        "Swap refinement: %d swaps, %d evaluations, stopped at %s, scored on %s",
        result.swaps,
        result.evaluations,
        result.stop_reason,
        "the shared pool objective" if objective is not None else "covered bases",
    )
    return list(result.primers)
