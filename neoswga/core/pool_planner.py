"""Select the smallest evaluated oligo panel meeting explicit design targets."""

import math
import time

from .base_optimizer import OptimizationStatus
from .candidate_source import INVENTORY_EXHAUSTED, as_candidate_source
from .dimer_validator import DimerValidator
from .lazy_dimer import LazyDimerCompatibility
from .panel_beam import beam_search
from .pool_objective import PoolConstraints, PoolObjective
from .swap_refinement import attach_search_config, coverage_bins, refine_by_swaps

# Multiples of the configured reach to report coverage at. The window radius is
# a design-density convention rather than a measured extension distribution
# (see `coverage.polymerase_extension_reach`), and coverage is close to linear
# in it over this range, so a single figure without its reach says little. 1.0
# is included so the sweep contains the number the recommendation was made on.
REACH_SENSITIVITY_FACTORS = (1 / 3, 2 / 3, 1.0, 5 / 3, 10 / 3)

# How many partial panels the repair beam carries. Four rather than a
# larger number because the cost is linear in it and the beam runs only
# after the cheaper swap repair has already failed.
_BEAM_WIDTH = 4


def reach_sensitivity(cache, primers, prefixes, seq_lengths, reach, circular=False):
    """Coverage for one delivered panel across a range of extension reaches.

    Answers "how much of this number is the panel and how much is the window
    radius", which a single coverage figure cannot. Uses the same union-of-
    windows calculation the headline figure comes from, so the row at factor
    1.0 reproduces it rather than approximating it.

    Returns one dict per reach with ``reach``, ``factor`` and ``coverage``.
    """
    from .coverage import compute_per_prefix_coverage

    rows = []
    for factor in REACH_SENSITIVITY_FACTORS:
        scaled = max(1, int(round(reach * factor)))
        if any(row["reach"] == scaled for row in rows):
            continue
        aggregate, _ = compute_per_prefix_coverage(
            cache=cache,
            primers=list(primers),
            prefixes=list(prefixes),
            seq_lengths=list(seq_lengths),
            extension=scaled,
            circular=circular,
        )
        rows.append({"reach": scaled, "factor": factor, "coverage": aggregate})
    return sorted(rows, key=lambda row: row["reach"])


def _assess(objective, validator, primers):
    """Metrics, coverage and every constraint this panel fails.

    The dimer guard stays OUTSIDE the objective: it is a hard constraint on the
    delivered panel, not a scoring term. Folding it in among the others is how
    it became tradeable, and the relaxation that followed produced an 11 bp
    heterodimer against a configured 3.
    """
    metrics = objective.metrics(primers)
    pairs = validator.incompatible_pairs(primers)
    self_dimers = [p for p in primers if validator.has_self_dimer(p)]
    reasons = list(objective.violations(primers))
    if pairs or self_dimers:
        reasons.append("dimer constraint")
    return metrics, objective.coverage(primers), reasons, pairs, self_dimers


def _beam_candidates(pool, incumbent, size, budget):
    """The widest slice of the pool the remaining budget affords.

    The pool arrives in step-2 rank order, so a prefix of it is the best-ranked
    candidates rather than an arbitrary subset. The beam costs about
    `beam width * candidates * panel size` evaluations, so inverting that bound
    gives the number of candidates the budget will pay for.

    The incumbent is always included, wherever it ranks. Without that the beam
    could not reproduce the panel it was asked to improve on, and a repair that
    found nothing better would return something worse.
    """
    affordable = max(size, budget // (_BEAM_WIDTH * max(size, 1)))
    slice_pool = list(pool[:affordable])
    known = set(slice_pool)
    slice_pool.extend(primer for primer in incumbent if primer not in known)
    return slice_pool


def _repair(primers, pool, objective, reasons, config, target=None, bins=None, weights=None):
    """A bounded second attempt at a panel that missed a limit or a target.

    Returns the panel to use and a record of what was tried. The panel is only
    replaced when the search actually moved, and the caller re-evaluates
    whatever comes back rather than trusting this to have improved anything.

    Two different things bring a row here, and the difference stays visible in
    `reason`. A CONSTRAINT violation means the panel is not deliverable. A
    TARGET miss means it is deliverable and smaller than what was asked for.
    Both are worth a second attempt; only the first makes the row ineligible.
    The search order is the same either way -- fewer violations, then coverage,
    then background load -- so chasing a target can never buy coverage by giving
    up a constraint.

    A dimer violation is not repaired here. The objective does not see dimers,
    so the swap score cannot be steered by them, and a panel arriving with a
    dimerising pair means an upstream relaxation fired -- which is a thing to
    fix where it happens rather than to paper over at reporting time.
    """
    skipped = dict(
        attempted=False, method=None, reason=None, succeeded=False, swaps=0, evaluations=0
    )
    repairable = [r for r in reasons if r != "dimer constraint"]
    if not repairable and target is None:
        return list(primers), skipped
    reason = repairable[0] if repairable else "below coverage target"

    # The scan width needs the bin decomposition, because the prescreen that
    # decides which pairs are worth scoring IS the bin rule. This call passed
    # None for both bins for as long as it existed, so the width could not be
    # honoured here even once it existed; `plan_pool` now supplies them.
    width = getattr(config, "objective_scan_width", None) if bins is not None else None
    result = refine_by_swaps(
        primers,
        pool,
        bins,
        weights,
        LazyDimerCompatibility(config.max_dimer_bp),
        objective=objective,
        objective_scan_width=width,
        max_evaluations=config.swap_max_evaluations,
        max_seconds=config.swap_max_seconds,
    )
    repaired = list(result.primers)

    def _resolved(panel):
        """Whether this panel answers what brought the row here.

        A constraint repair succeeds by clearing the violations. A target repair
        succeeds by reaching the target, not merely by improving: the swap loop
        almost always improves something, and treating that as success returned
        before the beam had a turn, which is how the beam stayed unreached even
        once it had a budget of its own.
        """
        if objective.violations(panel):
            return False
        if repairable:
            return True
        covered = objective.coverage(panel)
        return covered is not None and covered >= target

    record = dict(
        attempted=True,
        method="swap",
        reason=reason,
        succeeded=_resolved(repaired),
        swaps=result.swaps,
        evaluations=result.evaluations,
        stop_reason=result.stop_reason,
        # What the scan actually spent. `evaluations` alone cannot separate a
        # narrow scan that finished from a wide one that ran out, and the
        # difference decides whether a disappointing row wants a wider scan or
        # a bigger budget.
        pairs_considered=result.pairs_considered,
        objective_evaluations=result.objective_evaluations,
        scan_width=result.scan_width,
    )
    if record["succeeded"]:
        return repaired, record

    # Swaps move one primer at a time from where the optimizer stopped, so they
    # cannot reach a panel that shares no primer with it. A beam rebuilds at the
    # same size and keeps several partial panels alive, which is what a
    # non-monotonic density floor needs.
    #
    # It used to be skipped whenever `beam width * pool size * panel size`
    # exceeded the budget left over from the swaps. On the 2,000-candidate
    # shortlist this tool actually produces, that permitted a panel of one
    # primer, so the beam never ran outside its own tests. Asking whether the
    # whole pool fits was the wrong question: the answer is to search the widest
    # slice the budget affords.
    size = len(repaired)
    budget = getattr(config, "beam_max_evaluations", config.swap_max_evaluations)
    if size < 1 or budget < _BEAM_WIDTH * size:
        record["beam"] = "not affordable within the remaining budget"
        # A repair that did not succeed returns the panel it was given, not the one
        # the search wandered to. The row is rejected either way, so a partial move
        # buys nothing, and it can cost: the search ranks on distance from
        # feasibility, so on a constraint no panel can meet it will trade real
        # coverage for a step toward a floor it will never reach. Keeping the
        # original makes that ordering safe -- progress toward feasibility is
        # pursued, and failing to get there costs nothing.
        return list(primers), record

    slice_pool = _beam_candidates(pool, repaired, size, budget)
    record["beam_candidates"] = len(slice_pool)
    beam = beam_search(
        slice_pool,
        objective,
        size,
        dimerises=LazyDimerCompatibility(config.max_dimer_bp).dimerises,
        beam_width=_BEAM_WIDTH,
        max_evaluations=budget,
        max_seconds=getattr(config, "beam_max_seconds", config.swap_max_seconds),
    )
    record["evaluations"] += beam.evaluations
    record["beam"] = beam.status
    # Only a qualifying panel of the SAME size is a repair, and only one that
    # answers what brought the row here. The beam also reports the best smaller
    # feasible panel it saw, which is a useful answer to a different question:
    # this row was asked for a panel of one size, and returning a shorter one
    # would show up as a different row's result. And a beam that rediscovers the
    # incumbent has repaired nothing, so it must not report success.
    if beam.violations or len(beam.primers) != size or not _resolved(list(beam.primers)):
        # Nothing qualified, so the panel this row reports is the one it arrived
        # with. See the note above.
        return list(primers), record
    record.update(method="beam", succeeded=True)
    return list(beam.primers), record


def _evaluate_size(
    requested,
    optimizer,
    pool,
    validator,
    objective,
    targets,
    repair,
    repair_bins,
    repair_weights,
    background_known,
    progress,
):
    """One size row: optimise, assess, repair, and report what happened.

    Extracted so the frontier refill loop can run a size more than once without
    duplicating any of this.
    """
    started = time.monotonic()
    if progress:
        progress(requested)
    result = optimizer.optimize(pool, target_size=requested)
    primers = list(dict.fromkeys(result.primers))
    if not primers or result.status not in {
        OptimizationStatus.SUCCESS,
        OptimizationStatus.PARTIAL,
    }:
        return dict(
            requested_size=requested,
            size=len(primers),
            primers=primers,
            status="no_panel",
            message=result.message,
            eligible=False,
            coverage=None,
            seconds=time.monotonic() - started,
        )
    if len(primers) > requested or not set(primers).issubset(pool):
        raise ValueError("Optimizer returned a panel outside the requested candidate/size bounds")

    metrics, coverage, reasons, violations, self_dimers = _assess(objective, validator, primers)
    repair_record = dict(
        attempted=False, method=None, reason=None, succeeded=False, swaps=0, evaluations=0
    )
    # The most demanding request drives the repair. Repairing to the lowest
    # would stop as soon as the easiest target was cleared, and every higher
    # row would report `not_found` beside a panel never asked to reach it.
    missed = max(targets) if coverage is not None and coverage < max(targets) else None
    if repair and (reasons or missed is not None):
        repaired, repair_record = _repair(
            primers,
            pool,
            objective,
            reasons,
            optimizer.config,
            target=missed,
            bins=repair_bins,
            weights=repair_weights,
        )
        if repaired != primers:
            primers = repaired
            metrics, coverage, reasons, violations, self_dimers = _assess(
                objective, validator, primers
            )
    for value in (metrics.fg_coverage, metrics.effective_fg_coverage):
        if value is not None and (not math.isfinite(value) or not 0 <= value <= 1):
            raise ValueError("Optimizer returned invalid coverage")
    if background_known and not math.isfinite(metrics.selectivity_density):
        raise ValueError("Optimizer returned non-finite specificity")
    return dict(
        requested_size=requested,
        size=len(primers),
        primers=primers,
        status="evaluated",
        eligible=not reasons,
        failed_constraints=reasons,
        coverage=coverage,
        raw_coverage=metrics.fg_coverage,
        effective_coverage=metrics.effective_fg_coverage,
        selectivity_density=metrics.selectivity_density if background_known else None,
        background_sites=metrics.total_bg_sites if background_known else None,
        violating_pairs=len(violations),
        self_dimers=len(self_dimers),
        repair=repair_record,
        max_gap=metrics.max_gap,
        seconds=time.monotonic() - started,
    )


def _prepare_candidate_pool(optimizer, source, primer_length):
    """The candidates a design may select, with the validator and the bins.

    Introduces the source to the cache its candidates will be scored on, which
    is what `ensure_positions` needed all along: it existed before this, was
    tested, and never ran, because no production caller attached a cache and
    the provider returned quietly when none was attached. The one configuration
    it was written to catch was the one it could not see. Attaching also gives
    the source what it needs to load positions for a candidate the frontier
    admits later.

    The position check is gated on the optimizer having brought a cache, which
    every command-line path does. A caller passing an optimizer-shaped object
    with a metrics table of its own has no index for this to check against.
    """
    position_cache = getattr(optimizer, "cache", None)
    if position_cache is not None and hasattr(source, "attach_positions"):
        source.attach_positions(position_cache)
    return _vet_frontier(optimizer, source, primer_length, source.initial())


def _vet_frontier(optimizer, source, primer_length, sequences):
    """Check and decompose one frontier, whether the first or a refilled one.

    Split from `_prepare_candidate_pool` because a refill must vet the widened
    window without re-drawing it: `source.initial()` resets the frontier, so
    calling it again after `advance()` would undo the refill.
    """
    position_cache = getattr(optimizer, "cache", None)
    pool = list(dict.fromkeys(p.upper() for p in sequences if len(p) == primer_length))
    if not pool:
        raise ValueError(
            f"No {primer_length}-mer candidates available; regenerate the candidate pool at that length"
        )
    if any(set(p) - set("ACGT") for p in pool):
        raise ValueError("Candidate oligos must contain only A, C, G and T")
    validator = DimerValidator(optimizer.config.max_dimer_bp, optimizer.config.max_self_dimer_bp)
    pool = validator.filter_self_dimers(pool)
    if not pool:
        raise ValueError("No candidates pass the configured self-dimer limit")

    # A candidate with no entry on one of the references is not a low-scoring
    # candidate, it is an unmeasured one, and the design would read its absent
    # host sites as perfect specificity.
    if position_cache is not None and hasattr(source, "ensure_positions"):
        source.ensure_positions(pool)

    # The coverage decomposition the repair's prescreen needs, built once per
    # run rather than once per size row: it is over the whole candidate pool
    # and does not depend on the panel. `(None, None)` for an optimizer that
    # builds no coverage graph, which then runs the unbounded scan and reports
    # `scan_width: null` rather than claiming a bound it did not have.
    bins, weights = coverage_bins(optimizer, pool)
    return pool, validator, bins, weights


def plan_pool(
    optimizer,
    candidates,
    sizes,
    coverage_targets,
    *,
    primer_length=12,
    min_selectivity_density=None,
    max_background_sites=None,
    coverage_metric="effective",
    repair=True,
    progress=None,
):
    """Re-optimize each size; retain panels satisfying all requested constraints.

    Recommendations are the smallest panels found, not proofs of minimum size.
    Binding density ratios are model metrics, not measured fold enrichment.
    """
    sizes = sorted(set(sizes))
    targets = sorted(set(coverage_targets))
    if not isinstance(primer_length, int) or primer_length < 1:
        raise ValueError("Primer length must be a positive integer")
    if not sizes or any(not isinstance(n, int) or n < 1 for n in sizes):
        raise ValueError("Panel sizes must be positive integers")
    if not targets or any(not math.isfinite(t) or not 0 < t <= 1 for t in targets):
        raise ValueError("Coverage targets must be fractions in (0, 1]")
    if coverage_metric not in {"raw", "effective"}:
        raise ValueError("Coverage metric must be raw or effective")
    for value in (min_selectivity_density, max_background_sites):
        if value is not None and (not math.isfinite(value) or value < 0):
            raise ValueError("Specificity limits must be finite and non-negative")
    background_known = bool(optimizer.bg_prefixes and sum(optimizer.bg_seq_lengths) > 0)
    # One contract for search and acceptance. `plan_pool` used to restate the
    # coverage choice and both specificity limits in its own words, and two
    # copies of a rule drift while each stays self-consistent.
    constraints = PoolConstraints(
        coverage_metric=coverage_metric,
        min_selectivity_density=min_selectivity_density,
        max_background_sites=max_background_sites,
    )
    constraints.require_background(available=background_known)
    # The focused evaluator when the optimizer has one, because the repair path
    # calls it thousands of times per design and reads five of its fields. A
    # caller passing its own optimizer-shaped object keeps the full one.
    evaluate = getattr(optimizer, "compute_pool_metrics", None) or optimizer.compute_metrics
    objective = PoolObjective(evaluate, constraints)
    if background_known and min_selectivity_density is None and max_background_sites is None:
        raise ValueError("Specify a minimum selectivity density or maximum background sites")
    if coverage_metric == "effective" and optimizer.conditions is None:
        raise ValueError("Effective coverage requires reaction conditions")
    # One contract for where candidates come from. A bare list still works and
    # becomes a `ListCandidateSource`, so every existing caller is unchanged;
    # a command that opened the inventory passes that instead, and the plan
    # then reports how much of the universe the run actually examined.
    source = as_candidate_source(candidates)
    pool, validator, repair_bins, repair_weights = _prepare_candidate_pool(
        optimizer, source, primer_length
    )
    # Stage-2 refinement inside the optimizer picks the delivered panel, so it
    # has to score on the same thing this function accepts on. See
    # `swap_refinement.refine_hybrid_stage2`, which reads this attribute.
    # Through `attach_search_config`, because the optimizer `plan_pool` is
    # handed is a wrapper and the stage that reads this is the delegate inside
    # it. Setting the attribute directly reached the wrapper only, so the
    # refinement received None on every command-line design.
    attach_search_config(optimizer, "pool_objective", objective)
    rows = []
    # A row that cannot satisfy its constraints asks the source for more
    # candidates before giving up. Without this the inventory is decorative:
    # 20,670 candidates are eligible on the Wolbachia design and 2,000 were
    # ever searched. A row that already qualifies never refills.
    max_refills = getattr(optimizer.config, "max_frontier_refills", 0)
    for requested in sizes:
        refills = 0
        while True:
            row = _evaluate_size(
                requested,
                optimizer,
                pool,
                validator,
                objective,
                targets,
                repair,
                repair_bins,
                repair_weights,
                background_known,
                progress,
            )
            qualified = row.get("status") == "evaluated" and not row.get("failed_constraints")
            if qualified or refills >= max_refills:
                break
            if not hasattr(source, "advance") or not source.advance(keep=row.get("primers") or ()):
                break
            refills += 1
            # The widened frontier is vetted, not assumed: it admits candidates
            # the position cache was not built over, which is the case
            # `PositionCache.load` and `ensure_positions` exist for.
            pool, validator, repair_bins, repair_weights = _vet_frontier(
                optimizer, source, primer_length, source.frontier()
            )
        row["frontier_refills"] = refills
        # Which kind of "nothing left" this row ended on. `frontier_exhausted`
        # means the universe still holds candidates and a bigger budget would
        # reach them; `inventory_exhausted` means it does not.
        row["candidates_exhausted"] = (
            source.exhaustion() if hasattr(source, "exhaustion") else INVENTORY_EXHAUSTED
        )
        rows.append(row)
    recommendations = []
    for target in targets:
        qualifying = [
            (i, r) for i, r in enumerate(rows) if r["eligible"] and r["coverage"] >= target
        ]
        chosen = (
            min(qualifying, key=lambda item: (item[1]["size"], -item[1]["coverage"], item[0]))
            if qualifying
            else None
        )
        recommendations.append(
            dict(
                target_coverage=target,
                row_index=chosen[0] if chosen else None,
                size=chosen[1]["size"] if chosen else None,
                status="smallest_found" if chosen else "not_found",
            )
        )
    # Reach sensitivity for the panel a reader is most likely to act on: the
    # largest qualifying recommendation. Computed once, saved with the plan, so
    # `report-pool` can render it without genome access.
    chosen = [r for r in recommendations if r["row_index"] is not None]
    sensitivity = []
    # Checked rather than caught: a caller may pass any optimizer-shaped object,
    # and one that cannot supply a position cache simply gets no sweep. A bare
    # `except` here would also hide a real cache failure.
    needed = ("cache", "fg_prefixes", "fg_seq_lengths", "config")
    if chosen and all(getattr(optimizer, name, None) is not None for name in needed):
        reach = getattr(optimizer.config, "extension_reach", 0)
        if reach > 0:
            panel = rows[max(chosen, key=lambda r: r["target_coverage"])["row_index"]]["primers"]
            sensitivity = reach_sensitivity(
                cache=optimizer.cache,
                primers=panel,
                prefixes=optimizer.fg_prefixes,
                seq_lengths=optimizer.fg_seq_lengths,
                reach=reach,
                circular=getattr(optimizer.config, "fg_circular", False),
            )

    return dict(
        primer_length=primer_length,
        candidate_count=len(pool),
        candidate_source=source.describe(),
        reach_sensitivity=sensitivity,
        coverage_metric=coverage_metric,
        extension_reach=optimizer.config.extension_reach,
        background_assessed=background_known,
        min_selectivity_density=min_selectivity_density,
        max_background_sites=max_background_sites,
        max_dimer_bp=optimizer.config.max_dimer_bp,
        max_self_dimer_bp=optimizer.config.max_self_dimer_bp,
        optimizer=optimizer.name,
        rows=rows,
        recommendations=recommendations,
        interpretation="Smallest qualifying panels found among evaluated sizes; not proven global minima. "
        "Coverage and specificity are computational estimates, not lab measurements.",
    )
