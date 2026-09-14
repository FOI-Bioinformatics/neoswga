"""Select the smallest evaluated oligo panel meeting explicit design targets."""

import math
import time

from .base_optimizer import OptimizationStatus
from .dimer_validator import DimerValidator


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
    if not background_known and (
        min_selectivity_density is not None or max_background_sites is not None
    ):
        raise ValueError("Specificity limits require a background genome and index")
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
        metrics = optimizer.compute_metrics(primers)
        for value in (metrics.fg_coverage, metrics.effective_fg_coverage):
            if value is not None and (not math.isfinite(value) or not 0 <= value <= 1):
                raise ValueError("Optimizer returned invalid coverage")
        if background_known and not math.isfinite(metrics.selectivity_density):
            raise ValueError("Optimizer returned non-finite specificity")
        coverage = (
            metrics.effective_fg_coverage if coverage_metric == "effective" else metrics.fg_coverage
        )
        violations = validator.incompatible_pairs(primers)
        self_dimers = [p for p in primers if validator.has_self_dimer(p)]
        reasons = []
        if coverage is None:
            reasons.append("coverage unavailable")
        if violations or self_dimers:
            reasons.append("dimer constraint")
        if (
            min_selectivity_density is not None
            and metrics.selectivity_density < min_selectivity_density
        ):
            reasons.append("selectivity below minimum")
        if max_background_sites is not None and metrics.total_bg_sites > max_background_sites:
            reasons.append("background sites above maximum")
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
    return dict(
        primer_length=primer_length,
        candidate_count=len(pool),
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
