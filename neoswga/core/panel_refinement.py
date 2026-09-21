"""Shared objective for proposal refinement and pool-size reduction."""

from copy import copy
from dataclasses import is_dataclass, replace

from .pool_objective import PoolConstraints, PoolObjective
from .swap_refinement import attach_search_config, refine_by_swaps


def objective_for_optimizer(optimizer, constraints=None):
    """Reuse the design objective, or establish one before selection.

    Condition-free library evaluators retain geometric coverage. A configured
    design uses effective coverage; failure to measure it is not a raw fallback.
    """
    existing = getattr(optimizer, "pool_objective", None)
    if existing is not None:
        return existing
    if not callable(getattr(optimizer, "compute_metrics", None)):
        if constraints is not None:
            raise ValueError("Configured panel constraints require a panel evaluator")
        return None
    if constraints is None:
        constraints = PoolConstraints(
            coverage_metric=(
                "effective" if getattr(optimizer, "conditions", None) is not None else "raw"
            )
        )
    constraints.require_background(
        bool(
            getattr(optimizer, "bg_prefixes", None)
            and sum(getattr(optimizer, "bg_seq_lengths", []))
        )
    )
    # Reuse the planner's focused evaluator when it measures every requested
    # quantity. Gap-evenness and host-coverage limits require the full metrics.
    extended_limits = (
        constraints.max_mean_gap,
        constraints.max_evenness,
        constraints.max_host_coverage,
    )
    evaluate = optimizer.compute_metrics
    if all(value is None for value in extended_limits):
        evaluate = getattr(optimizer, "compute_pool_metrics", None) or evaluate
    objective = PoolObjective(evaluate, constraints)
    attach_search_config(optimizer, "pool_objective", objective)
    return objective


def refine_result(result, optimizer, candidates, fixed_primers=(), diagnostics=None):
    """Refine any method's proposal with the same panel acceptance contract.

    Preserve the incumbent if no feasible improvement is found. Dimer rules
    remain hard constraints, including the optional temperature-based floor.
    """
    from .base_optimizer import OptimizationStatus
    from .dimer_validator import DimerValidator
    from .lazy_dimer import dimer_screen

    if not result.primers or result.status not in {
        OptimizationStatus.SUCCESS,
        OptimizationStatus.PARTIAL,
    }:
        return result
    config = getattr(optimizer, "config", None)
    objective = objective_for_optimizer(optimizer)
    if config is None or objective is None:
        return result
    validator = DimerValidator(config.max_dimer_bp, config.max_self_dimer_bp)
    pool = validator.filter_self_dimers(list(dict.fromkeys([*result.primers, *candidates])))
    compatibility = dimer_screen(
        pool,
        config.max_dimer_bp,
        max_dimer_dg=config.max_dimer_dg,
        temp=float(getattr(getattr(optimizer, "conditions", None), "temp", 37.0)),
    )
    attempt = refine_by_swaps(
        result.primers,
        pool,
        None,
        None,
        compatibility,
        fixed_primers=fixed_primers,
        objective=objective,
        max_evaluations=config.swap_max_evaluations,
        max_seconds=config.swap_max_seconds,
    )
    if diagnostics is not None:
        diagnostics.update(
            evaluations=attempt.evaluations,
            stop_reason=attempt.stop_reason,
            swaps=attempt.swaps,
            max_evaluations=config.swap_max_evaluations,
            max_seconds=config.swap_max_seconds,
        )
    panel = list(attempt.primers)
    if tuple(panel) == tuple(result.primers) or objective.violations(panel):
        return result
    if any(validator.has_self_dimer(p) for p in panel):
        return result
    if any(compatibility.dimerises(p, panel[:i]) for i, p in enumerate(panel)):
        return result
    metrics = optimizer.compute_metrics(panel)
    if not is_dataclass(result):
        updated = copy(result)
        updated.primers = list(panel)
        return updated
    return replace(
        result,
        primers=tuple(panel),
        metrics=metrics,
        score=metrics.normalized_score(),
        message=f"{result.message} Shared-objective refinement: {attempt.swaps} swaps, "
        f"{attempt.evaluations} evaluations, {attempt.stop_reason}.",
    )
