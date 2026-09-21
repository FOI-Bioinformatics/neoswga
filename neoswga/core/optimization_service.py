"""Shared execution of panel proposals, refinement and constrained reduction.

Front ends resolve references and chemistry before constructing a request. The
optimizer owns that resolved evaluator; all stages reuse its objective.
"""

import math
import time
from copy import copy
from dataclasses import dataclass, is_dataclass, replace
from typing import Any

from .base_optimizer import OptimizationStatus
from .exceptions import NoCandidatesError
from .panel_refinement import objective_for_optimizer, refine_result


@dataclass(frozen=True)
class OptimizationRequest:
    """One panel search under a fixed evaluator and candidate frontier."""

    optimizer: Any
    candidates: tuple[str, ...]
    target_size: int
    constraints: Any = None
    fixed_primers: tuple[str, ...] = ()
    excluded_primers: tuple[str, ...] = ()
    refine: bool = True
    minimize: bool = False
    target_coverage: float = 0.7
    min_size: int = 1
    repair: bool = True
    repair_target: float | None = None
    repair_bins: Any = None
    repair_weights: Any = None
    budget: Any = None
    candidate_source: Any = None
    prepare_candidates: Any = None

    def __post_init__(self):
        for field in ("candidates", "fixed_primers", "excluded_primers"):
            object.__setattr__(
                self,
                field,
                tuple(dict.fromkeys(str(p).strip().upper() for p in getattr(self, field))),
            )
        if self.target_size < 1 or not 1 <= self.min_size <= self.target_size:
            raise ValueError("Invalid pool-size bounds")
        if not math.isfinite(self.target_coverage) or not 0 <= self.target_coverage <= 1:
            raise ValueError("Coverage target must be in [0, 1]")
        if len(set(self.fixed_primers)) > self.target_size:
            raise ValueError("Fixed primers exceed the requested pool size")
        if set(self.fixed_primers) & set(self.excluded_primers):
            raise ValueError("A fixed primer cannot also be excluded")


def panel_rank(objective, primers):
    """Constraint-first ordering, shared by proposal acceptance and search."""
    coverage = objective.coverage(primers)
    return (
        -objective.shortfall(primers),
        coverage if coverage is not None else -math.inf,
        -objective.metrics(primers).total_bg_sites,
    )


def panel_violations(optimizer, primers):
    """Panel limits plus the configured sequence and temperature dimer checks."""
    from .dimer_validator import DimerValidator
    from .lazy_dimer import LazyDimerCompatibility

    objective = objective_for_optimizer(optimizer)
    reasons = list(objective.violations(primers)) if objective else []
    config = getattr(optimizer, "config", None)
    if config is not None:
        validator = DimerValidator(config.max_dimer_bp, config.max_self_dimer_bp)
        compatibility = LazyDimerCompatibility(
            config.max_dimer_bp,
            max_dimer_dg=config.max_dimer_dg,
            temp=float(getattr(getattr(optimizer, "conditions", None), "temp", 37.0)),
        )
        if any(
            validator.has_self_dimer(p) or compatibility.dimerises(p, primers[:i])
            for i, p in enumerate(primers)
        ):
            reasons.append("dimer constraint")
    return tuple(reasons)


def reduce_result(
    result, optimizer, target_coverage, fixed_primers=(), min_size=1, diagnostics=None
):
    """Delete redundant primers without losing fixed oligos or panel limits."""
    objective = objective_for_optimizer(optimizer)
    if objective is None:
        return result
    current = list(result.primers)
    fixed = set(fixed_primers)
    config = getattr(optimizer, "config", None)
    budget = getattr(config, "swap_max_evaluations", 10000)
    seconds = getattr(config, "swap_max_seconds", 30.0)
    deadline = time.monotonic() + seconds
    evaluations = 0
    removals = []
    while len(current) > max(min_size, len(fixed)):
        proposals = []
        for primer in current:
            if primer in fixed:
                continue
            if evaluations >= budget or time.monotonic() >= deadline:
                break
            evaluations += 1
            panel = [p for p in current if p != primer]
            coverage = objective.coverage(panel)
            if (
                coverage is not None
                and coverage >= target_coverage
                and not panel_violations(optimizer, panel)
            ):
                proposals.append(panel)
        if not proposals:
            break
        following = max(proposals, key=lambda panel: panel_rank(objective, panel))
        removals.append(
            {
                "primer": next(p for p in current if p not in following),
                "resulting_coverage": objective.coverage(following),
            }
        )
        current = following
    if diagnostics is not None:
        diagnostics.update(
            evaluations=evaluations,
            removals=removals,
            max_evaluations=budget,
            max_seconds=seconds,
            stop_reason=(
                "evaluation_budget"
                if evaluations >= budget
                else "time_budget" if time.monotonic() >= deadline else "no_qualifying_deletion"
            ),
        )
    if tuple(current) == tuple(result.primers):
        return result
    metrics = optimizer.compute_metrics(current)
    return replace(
        result, primers=tuple(current), metrics=metrics, score=metrics.normalized_score()
    )


def repair_result(
    result, optimizer, candidates, target=None, fixed_primers=(), bins=None, weights=None
):
    """Apply the planner's bounded repair using this run's evaluator and chemistry."""
    from .dimer_validator import DimerValidator
    from .lazy_dimer import LazyDimerCompatibility
    from .pool_planner import repair_panel

    objective = objective_for_optimizer(optimizer)
    config = getattr(optimizer, "config", None)
    skipped = dict(
        attempted=False, method=None, reason=None, succeeded=False, swaps=0, evaluations=0
    )
    if objective is None or config is None or not result.primers:
        return result, skipped
    reasons = list(objective.violations(result.primers))
    coverage = objective.coverage(result.primers)
    missed = target if target is not None and coverage is not None and coverage < target else None
    if not reasons and missed is None:
        return result, skipped
    validator = DimerValidator(config.max_dimer_bp, config.max_self_dimer_bp)
    pool = validator.filter_self_dimers(list(candidates))
    temp = float(getattr(getattr(optimizer, "conditions", None), "temp", 37.0))
    panel, record = repair_panel(
        result.primers,
        pool,
        objective,
        reasons,
        config,
        target=missed,
        bins=bins,
        weights=weights,
        fixed_primers=fixed_primers,
        temp=temp,
    )
    compatibility = LazyDimerCompatibility(
        config.max_dimer_bp, max_dimer_dg=config.max_dimer_dg, temp=temp
    )
    if (
        not set(fixed_primers).issubset(panel)
        or not set(panel).issubset(candidates)
        or len(panel) > len(result.primers)
        or objective.violations(panel)
        or any(
            validator.has_self_dimer(p) or compatibility.dimerises(p, panel[:i])
            for i, p in enumerate(panel)
        )
    ):
        record["succeeded"] = False
        return result, record
    if tuple(panel) == tuple(result.primers):
        return result, record
    metrics = optimizer.compute_metrics(panel)
    if is_dataclass(result):
        return (
            replace(
                result, primers=tuple(panel), metrics=metrics, score=metrics.normalized_score()
            ),
            record,
        )
    updated = copy(result)
    updated.primers = list(panel)
    return updated, record


def _with_history(result, stages):
    if is_dataclass(result):
        return replace(result, stage_history=tuple(stages))
    # Planner library callers may provide a simple result-shaped object.
    updated = copy(result)
    updated.stage_history = tuple(stages)
    return updated


def _run_panel_stages(request, initial_result=None):
    """Select, refine, reduce and revisit refinement after each reduction.

    An existing proposal can enter at refinement. Reduction strictly decreases
    size, so the sequence terminates even when every smaller pool admits swaps.
    Each swap pass keeps the optimizer's configured evaluation/time budgets.
    """
    from .search_control import SearchBudgetExhausted, budgeted_objective

    optimizer = request.optimizer
    budget = request.budget
    objective = objective_for_optimizer(optimizer, request.constraints)
    pool = list(dict.fromkeys([*request.fixed_primers, *request.candidates]))
    pool = [p for p in pool if p not in set(request.excluded_primers)]
    fixed = set(request.fixed_primers)
    stages = list(getattr(initial_result, "stage_history", ()) or ())

    def check(result):
        panel = result.primers
        if (
            len(panel) != len(set(panel))
            or len(panel) > request.target_size
            or not set(panel).issubset(pool)
        ):
            raise ValueError(
                "Optimizer returned a panel outside the requested candidate/size bounds"
            )
        if panel and not fixed.issubset(panel):
            raise ValueError("Optimizer dropped fixed primers")

    def record(name, before, after, started):
        stages.append(
            dict(
                stage=name,
                before_size=len(before.primers) if before else 0,
                after_size=len(after.primers),
                coverage=objective.coverage(after.primers) if objective and after.primers else None,
                coverage_metric=(
                    getattr(objective, "metric_name", objective.constraints.coverage_metric)
                    if objective
                    else None
                ),
                failed_constraints=(
                    list(panel_violations(optimizer, after.primers)) if after.primers else []
                ),
                changed=before is None or before.primers != after.primers,
                conditions_fingerprint=(
                    optimizer.conditions.fingerprint()
                    if callable(
                        getattr(getattr(optimizer, "conditions", None), "fingerprint", None)
                    )
                    else None
                ),
                seconds=time.monotonic() - started,
                search_budget=budget.describe(),
            )
        )

    def execute(operation, incumbent):
        try:
            budget.check()
            with budgeted_objective(objective, budget):
                return operation()
        except SearchBudgetExhausted:
            return incumbent

    started = time.monotonic()
    result = initial_result
    if result is None:
        kwargs = {"fixed_primers": list(request.fixed_primers)} if request.fixed_primers else {}
        from .base_optimizer import OptimizationResult

        failure = OptimizationResult.failure(
            getattr(optimizer, "name", "search"),
            "Search budget exhausted before a panel was selected",
        )
        result = execute(
            lambda: optimizer.optimize(pool, target_size=request.target_size, **kwargs), failure
        )
        check(result)
        record("selection", None, result, started)
    else:
        check(result)
    if (
        result.status not in {OptimizationStatus.SUCCESS, OptimizationStatus.PARTIAL}
        or not result.primers
    ):
        return _with_history(result, stages)
    while True:
        try:
            budget.check()
        except SearchBudgetExhausted:
            break
        if request.repair:
            started = time.monotonic()
            previous = result
            result, details = execute(
                lambda: repair_result(
                    result,
                    optimizer,
                    pool,
                    request.repair_target,
                    request.fixed_primers,
                    request.repair_bins,
                    request.repair_weights,
                ),
                (previous, {"attempted": False, "stop_reason": budget.stop_reason}),
            )
            check(result)
            record("repair", previous, result, started)
            stages[-1].update(details)
        if request.refine:
            started = time.monotonic()
            previous = result
            diagnostics = {}
            result = execute(
                lambda: refine_result(
                    previous, optimizer, pool, request.fixed_primers, diagnostics
                ),
                previous,
            )
            check(result)
            record("refinement", previous, result, started)
            stages[-1].update(diagnostics)
        if not request.minimize:
            break
        started = time.monotonic()
        previous = result
        diagnostics = {}
        result = execute(
            lambda: reduce_result(
                result,
                optimizer,
                request.target_coverage,
                request.fixed_primers,
                request.min_size,
                diagnostics,
            ),
            previous,
        )
        check(result)
        record("reduction", previous, result, started)
        stages[-1].update(diagnostics)
        if len(result.primers) == len(previous.primers):
            break
    violations = panel_violations(optimizer, result.primers)
    if violations and is_dataclass(result):
        result = replace(
            result,
            status=OptimizationStatus.PARTIAL,
            message=f"{result.message} Unmet panel limits: {', '.join(violations)}.",
        )
    return _with_history(result, stages)


def run_panel_search(request, initial_result=None):
    """Run panel stages across one shared candidate source and search allowance."""
    from .base_optimizer import OptimizationResult
    from .search_control import SearchBudget, search_frontiers

    budget = request.budget or SearchBudget.from_config(getattr(request.optimizer, "config", None))
    request = replace(request, budget=budget)
    source = request.candidate_source
    if source is None:
        return _run_panel_stages(request, initial_result)
    optimizer = request.optimizer
    objective = objective_for_optimizer(optimizer, request.constraints)
    cache = getattr(optimizer, "cache", None)
    if cache is not None:
        source.attach_positions(cache)

    def prepare(sequences):
        from .dimer_validator import DimerValidator

        pool = list(dict.fromkeys(str(p).upper() for p in sequences))
        pool = [p for p in pool if p not in request.excluded_primers]
        if cache is not None:
            source.ensure_positions(pool)
        if request.prepare_candidates:
            pool = list(request.prepare_candidates(pool))
        config = optimizer.config
        kept = DimerValidator(config.max_dimer_bp, config.max_self_dimer_bp).filter_self_dimers(
            pool
        )
        if pool and not kept:
            # A screen that rejects every candidate must say so. Handing an
            # empty list onward produced "candidates list cannot be empty" from
            # deep inside an optimizer, which names neither the screen that
            # emptied the pool nor the threshold it applied, and reads like a
            # missing input rather than a rejected one.
            raise NoCandidatesError(
                len(pool),
                f"the self-dimer screen at max_self_dimer_bp={config.max_self_dimer_bp}",
            )
        return kept

    def assess(result):
        panel = result.primers
        valid = bool(panel) and not panel_violations(optimizer, panel)
        goal = request.repair_target
        if goal is None and request.minimize:
            goal = request.target_coverage
        coverage = objective.coverage(panel) if panel else None
        qualified = valid and len(panel) >= request.target_size
        if goal is not None:
            qualified = valid and coverage is not None and coverage >= goal
        rank = (
            (valid, *panel_rank(objective, panel))
            if panel
            else (False, -math.inf, -math.inf, -math.inf)
        )
        return qualified, rank, panel

    result, pool, history = search_frontiers(
        source,
        prepare(request.candidates),
        lambda pool: _run_panel_stages(
            replace(request, candidates=tuple(pool), candidate_source=None), initial_result
        ),
        assess,
        prepare,
        optimizer.config.max_frontier_refills,
        budget,
    )
    if result is None:
        result = initial_result or OptimizationResult.failure(
            getattr(optimizer, "name", "search"), "Search budget exhausted"
        )
    stage = dict(stage="frontier_search", **history)
    return _with_history(result, [*getattr(result, "stage_history", ()), stage])
