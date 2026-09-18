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
    # Pairs the cheap pass looked at, and panels the objective actually scored.
    # `evaluations` counts pairs in unbounded mode and objective scorings in
    # bounded mode, which is ambiguous on its own, so both are reported.
    pairs_considered: int = 0
    objective_evaluations: int = 0
    # The scan width in force, or None when every surviving pair was scored.
    scan_width: int | None = None


def attach_search_config(optimizer, name, value):
    """Set an attribute on the optimizer AND on the delegate that searches.

    `OptimizerFactory` returns a wrapper -- `HybridBaseOptimizer` or
    `BackgroundAwareBaseOptimizer` -- and each delegates the search to an inner
    `HybridOptimizer`. `_swap_refine` is a method of the INNER one, so an
    attribute set on the wrapper is invisible to the stage that reads it.

    That is not hypothetical. `plan_pool` attached `pool_objective` to the
    wrapper and `refine_hybrid_stage2` read it off the delegate, so on every
    command-line design the refinement received None: it refined on raw covered
    bases while the row was accepted on occupancy-weighted coverage under a
    specificity floor, which is the two-rule split Phase 2 set out to remove.

    Two tests covered it and neither could see it. One asserted by AST that
    `plan_pool` assigns the attribute; the other asserted by source text that
    the refinement reads it. Both ends existed and the path did not.
    `tests/test_the_objective_reaches_the_stage_that_refines.py` asserts the
    path instead, through a real factory-built optimizer.
    """
    setattr(optimizer, name, value)
    inner = getattr(optimizer, "_hybrid", None)
    if inner is not None and inner is not optimizer:
        setattr(inner, name, value)


def coverage_bins(optimizer, pool):
    """The optimizer's coverage decomposition, as plain bins and weights.

    Returns ``(None, None)`` when the optimizer has no bin decomposition to
    give, which is the honest answer for the methods that do not build one; a
    caller then runs the unbounded scan and says so rather than inventing an
    empty decomposition that would read as "nothing is covered".

    The lookup walks through a wrapper because `OptimizerFactory` returns
    `HybridBaseOptimizer`, which delegates the coverage graph to an inner
    `HybridOptimizer`. Asking the wrapper directly returned nothing and was the
    reason `plan_pool` passed no bins at all.
    """
    provider = None
    for candidate in (optimizer, getattr(optimizer, "_hybrid", None)):
        if candidate is not None and hasattr(candidate, "_coverage_bins_by_primer"):
            provider = candidate
            break
    if provider is None:
        return None, None
    regions = provider._coverage_bins_by_primer(list(dict.fromkeys(pool)))
    bins = {p: {provider._bin_key(r) for r in owned} for p, owned in regions.items()}
    weights = {provider._bin_key(r): r.end - r.start for owned in regions.values() for r in owned}
    return bins, weights


def _validate_inputs(
    bins_by_primer, bin_weights, objective, objective_scan_width, max_evaluations, max_seconds
):
    """Check the budgets and the scan width. Returns whether bins are absent."""
    if max_evaluations < 0 or not math.isfinite(max_seconds) or max_seconds < 0:
        raise ValueError("Swap budgets must be finite and non-negative")
    binless = bins_by_primer is None or bin_weights is None
    if binless and objective is None:
        raise ValueError("Swap refinement needs either coverage bins or an objective")
    if objective_scan_width is None:
        return binless
    if objective is None:
        raise ValueError(
            "objective_scan_width bounds the objective-scored scan, so it needs an "
            "objective. Without one the whole scan is already the cheap bin rule."
        )
    if binless:
        raise ValueError(
            "objective_scan_width needs coverage bins: the prescreen that decides "
            "which pairs are worth scoring IS the bin rule. Pass bins_by_primer and "
            "bin_weights, or drop the width and accept an unbounded scan."
        )
    if not isinstance(objective_scan_width, int) or isinstance(objective_scan_width, bool):
        raise ValueError("objective_scan_width must be an integer")
    if objective_scan_width < 1:
        raise ValueError("objective_scan_width must be at least 1")
    return binless


def _bin_gain(incoming, outgoing, bins_by_primer, bin_weights, counts):
    """The cheap rule: bases this swap newly covers, less those it drops.

    The prescreen in front of the objective, and not a new criterion: it is
    exactly what `refine_by_swaps` has always used when given no objective.
    """
    new_bins = bins_by_primer.get(incoming, set())
    old_bins = bins_by_primer.get(outgoing, set())
    gain = sum(bin_weights[b] for b in new_bins if counts[b] == 0)
    loss = sum(bin_weights[b] for b in old_bins - new_bins if counts[b] == 1)
    return gain - loss


def _ranked_admissible_pairs(
    current, pool, fixed, selected_set, bins_by_primer, bin_weights, counts, dimers
):
    """Every swap that is allowed, best bin gain first, plus how many were seen.

    The whole frontier is walked. This is deliberately NOT bounded by
    `max_evaluations`: that budget was calibrated when every pair cost a full
    objective evaluation, and letting it cut this pass short would rank a
    prefix of the pool and reintroduce the blindness the scan width removes.

    Ties break on the pool order, which is the caller's ranking, so the
    traversal is deterministic and independent of dictionary layout.
    """
    ranked = []
    seen = 0
    for order, incoming in enumerate(pool):
        if incoming in selected_set:
            continue
        for outgoing in current:
            if outgoing in fixed:
                continue
            seen += 1
            retained = [p for p in current if p != outgoing]
            if dimers.dimerises(incoming, retained):
                continue
            gain = _bin_gain(incoming, outgoing, bins_by_primer, bin_weights, counts)
            ranked.append((-gain, order, incoming, outgoing))
    ranked.sort()
    return ranked, seen


def _refine_with_bounded_scan(
    current,
    pool,
    fixed,
    bins_by_primer,
    bin_weights,
    counts,
    dimers,
    score,
    objective,
    objective_scan_width,
    max_evaluations,
    deadline,
    covered,
    background,
):
    """Swap refinement with the objective scored only over the best pairs.

    One round is: rank every admissible pair by the cheap bin gain, score the
    best `objective_scan_width` of them with the objective, take the best
    improvement. Repeat until a round finds nothing or a budget stops it.

    Separate from the unbounded loop rather than folded into it, because the
    two differ in what they spend and in what bounds them: the unbounded one
    scores every surviving pair and is bounded by the pair count, this one
    scores a fixed number per round and is bounded by rounds.
    """
    evaluations = swaps = pairs_considered = objective_evaluations = 0
    reason = "local_optimum"
    while True:
        best = None
        best_score = None
        incumbent = score(current)
        exhausted = False
        ranked, seen = _ranked_admissible_pairs(
            current, pool, fixed, set(current), bins_by_primer, bin_weights, counts, dimers
        )
        pairs_considered += seen
        for _, _, incoming, outgoing in ranked[:objective_scan_width]:
            if evaluations >= max_evaluations or time.monotonic() >= deadline:
                reason = "evaluation_limit" if evaluations >= max_evaluations else "time_limit"
                exhausted = True
                break
            evaluations += 1
            objective_evaluations += 1
            candidate_score = score([p for p in current if p != outgoing] + [incoming])
            if candidate_score > incumbent and (best is None or candidate_score > best_score):
                best_score = candidate_score
                best = (outgoing, incoming)
        if best is not None:
            outgoing, incoming = best
            for b in bins_by_primer.get(outgoing, ()):
                counts[b] -= 1
            counts.update(bins_by_primer.get(incoming, ()))
            current[current.index(outgoing)] = incoming
            covered = sum(bin_weights[b] for b in counts if counts[b])
            background = objective.metrics(current).total_bg_sites
            swaps += 1
        if exhausted or best is None:
            return SwapResult(
                tuple(current),
                evaluations,
                swaps,
                reason,
                covered,
                background,
                pairs_considered=pairs_considered,
                objective_evaluations=objective_evaluations,
                scan_width=objective_scan_width,
            )


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
    objective_scan_width=None,
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

    `objective_scan_width` bounds the objective-scored part of each round. The
    scan is over `candidates x panel`, which at two thousand candidates and a
    twelve-primer panel is twenty-four thousand pairs, and scoring every pair
    that survived the dimer guard is what made a real repair run out of its
    deadline. With a width, the cheap bin gain ranks the surviving pairs and
    only the best `width` of them are scored.

    The prescreen is the other branch of this same function, which is why this
    needs no new criterion: without an objective the rule has always been bin
    gain minus bin loss. Measured on the real pool, its top ten pairs were
    exactly the objective's top ten, so it filters rather than decides.

    Two things about the bound are deliberate. The cheap pass always covers the
    whole frontier: `max_evaluations` was calibrated when every pair cost a full
    evaluation, and letting it cut the prescreen short would rank a prefix of
    the pool and reintroduce the blindness the bound removes. And a width
    without bins, or without an objective, is refused rather than quietly run
    unbounded, because a caller that asked to be bounded and silently was not is
    the shape of most of this audit.

    New primers must be compatible with every retained primer. Time limits are
    cooperative, checked between evaluations; preprocessing belongs to the
    caller.
    """
    binless = _validate_inputs(
        bins_by_primer, bin_weights, objective, objective_scan_width, max_evaluations, max_seconds
    )
    bounded = objective_scan_width is not None
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
    pairs_considered = objective_evaluations = 0
    deadline = time.monotonic() + max_seconds
    reason = "local_optimum"

    def _score(panel):
        """Lexicographic, constraints first. Higher is better throughout."""
        # `-shortfall` rather than `-len(violations)`: counting ties whenever
        # two panels fail the same single constraint, and coverage then decides,
        # which let more search move a panel further from the limit it was
        # chasing. Zero shortfall is exactly feasibility, so a feasible panel
        # still outranks every infeasible one however much coverage that would
        # buy.
        return (
            -objective.shortfall(panel),
            objective.coverage(panel),
            -objective.metrics(panel).total_bg_sites,
        )

    if bounded:
        return _refine_with_bounded_scan(
            current,
            pool,
            fixed,
            bins_by_primer,
            bin_weights,
            counts,
            dimers,
            _score,
            objective,
            objective_scan_width,
            max_evaluations,
            deadline,
            covered,
            background,
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
                pairs_considered += 1
                retained = [p for p in current if p != outgoing]
                if dimers.dimerises(incoming, retained):
                    continue
                if objective is not None:
                    objective_evaluations += 1
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
    return SwapResult(
        tuple(current),
        evaluations,
        swaps,
        reason,
        covered,
        background,
        pairs_considered=pairs_considered,
        objective_evaluations=objective_evaluations,
        scan_width=None,
    )


def refine_hybrid_stage2(optimizer, primers, candidates, fixed_primers):
    """Stage-2 swap refinement for :class:`HybridOptimizer`.

    The glue that turns an optimizer's coverage bins into the plain
    ``(bins, weights, compatibility)`` inputs :func:`refine_by_swaps` takes.
    It lives beside the search rather than on the optimizer because it is
    specific to the swap method, and because ``hybrid_optimizer`` is the
    module the size ratchet is closest to.
    """
    import logging

    from .lazy_dimer import dimer_screen

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
    bins, weights = coverage_bins(optimizer, pool)
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
        dimer_screen(
            pool,
            optimizer.max_dimer_bp,
            max_dimer_dg=getattr(optimizer, "max_dimer_dg", None),
            temp=float(getattr(optimizer, "reaction_temp", 37.0) or 37.0),
        ),
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
