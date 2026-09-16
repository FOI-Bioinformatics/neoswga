"""Select the smallest evaluated oligo panel meeting explicit design targets."""

import math
import time

from .base_optimizer import OptimizationStatus
from .dimer_validator import DimerValidator
from .lazy_dimer import LazyDimerCompatibility
from .panel_beam import beam_search
from .pool_objective import PoolConstraints, PoolObjective
from .swap_refinement import refine_by_swaps

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


def _repair(primers, pool, objective, reasons, config):
    """A bounded second attempt at a panel that missed a repairable limit.

    Returns the panel to use and a record of what was tried. The panel is only
    replaced when the swap loop actually moved, and the caller re-evaluates
    whatever comes back rather than trusting this to have improved anything.

    A dimer violation is not repaired here. The objective does not see dimers,
    so the swap score cannot be steered by them, and a panel arriving with a
    dimerising pair means an upstream relaxation fired -- which is a thing to
    fix where it happens rather than to paper over at reporting time.
    """
    skipped = dict(attempted=False, method=None, succeeded=False, swaps=0, evaluations=0)
    repairable = [r for r in reasons if r != "dimer constraint"]
    if not repairable:
        return list(primers), skipped

    result = refine_by_swaps(
        primers,
        pool,
        None,
        None,
        LazyDimerCompatibility(config.max_dimer_bp),
        objective=objective,
        max_evaluations=config.swap_max_evaluations,
        max_seconds=config.swap_max_seconds,
    )
    repaired = list(result.primers)
    record = dict(
        attempted=True,
        method="swap",
        succeeded=not objective.violations(repaired),
        swaps=result.swaps,
        evaluations=result.evaluations,
        stop_reason=result.stop_reason,
    )
    if record["succeeded"]:
        return repaired, record

    # Swaps move one primer at a time from where the optimizer stopped, so they
    # cannot reach a panel that shares no primer with it. A beam rebuilds at the
    # same size and keeps several partial panels alive, which is what a
    # non-monotonic density floor needs. It is tried only when its bound fits
    # the same budget the swaps were given, because that bound is quadratic in
    # nothing but is still beam width times pool size times panel size.
    size = len(repaired)
    budget = config.swap_max_evaluations - result.evaluations
    if size < 1 or _BEAM_WIDTH * len(pool) * size > budget:
        record["beam"] = "not affordable within the remaining budget"
        return repaired, record

    beam = beam_search(
        pool,
        objective,
        size,
        dimerises=LazyDimerCompatibility(config.max_dimer_bp).dimerises,
        beam_width=_BEAM_WIDTH,
        max_evaluations=budget,
        max_seconds=config.swap_max_seconds,
    )
    record["evaluations"] += beam.evaluations
    record["beam"] = beam.status
    # Only a qualifying panel of the SAME size is a repair. The beam also
    # reports the best smaller feasible panel it saw, which is a useful answer
    # to a different question: this row was asked for a panel of one size, and
    # returning a shorter one would show up as a different row's result.
    if beam.violations or len(beam.primers) != size:
        return repaired, record
    record.update(method="beam", succeeded=True)
    return list(beam.primers), record


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
    objective = PoolObjective(optimizer.compute_metrics, constraints)
    if background_known and min_selectivity_density is None and max_background_sites is None:
        raise ValueError("Specify a minimum selectivity density or maximum background sites")
    if coverage_metric == "effective" and optimizer.conditions is None:
        raise ValueError("Effective coverage requires reaction conditions")
    pool = list(dict.fromkeys(p.upper() for p in candidates if len(p) == primer_length))
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
    rows = []
    for requested in sizes:
        started = time.monotonic()
        if progress:
            progress(requested)
        result = optimizer.optimize(pool, target_size=requested)
        primers = list(dict.fromkeys(result.primers))
        if not primers or result.status not in {
            OptimizationStatus.SUCCESS,
            OptimizationStatus.PARTIAL,
        }:
            rows.append(
                dict(
                    requested_size=requested,
                    size=len(primers),
                    primers=primers,
                    status="no_panel",
                    message=result.message,
                    eligible=False,
                    coverage=None,
                    seconds=time.monotonic() - started,
                )
            )
            continue
        if len(primers) > requested or not set(primers).issubset(pool):
            raise ValueError(
                "Optimizer returned a panel outside the requested candidate/size bounds"
            )
        metrics, coverage, reasons, violations, self_dimers = _assess(objective, validator, primers)
        repair_record = dict(attempted=False, method=None, succeeded=False, swaps=0, evaluations=0)
        if repair and reasons:
            repaired, repair_record = _repair(primers, pool, objective, reasons, optimizer.config)
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
        rows.append(
            dict(
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
        )
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
