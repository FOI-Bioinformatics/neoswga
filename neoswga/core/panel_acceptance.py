"""Holding an `optimize` run to a limit, not only a `plan-pool` one.

`PoolConstraints` is enforced inside the search through `PoolObjective.shortfall`
and repaired by a bounded swap or beam, and for as long as it existed only
`plan-pool` constructed one. The audit of 2026-09-18 recorded that as the one
piece of architecture no published SWGA tool has, reaching one command out of
several: swga 1.0 has enforced two hard constraints inside its clique search
since 2017, and `optimize` enforced none.

This module is the bridge. It reads the limits a run configured, builds the same
`PoolConstraints` `plan-pool` builds, evaluates the delivered panel against them
and attempts the same bounded repair.

Every numeric limit is unset by default. The shared service still selects an
explicit coverage metric when no extra panel limits were requested. Configured
limits remain hard requirements rather than application-weighted score terms.

The dimer guarantee stays outside, as it does in `PoolObjective`: it is a hard
constraint on the delivered panel rather than a scoring term, and folding it in
among tradeable terms is how it came to be traded.
"""

from __future__ import annotations

import logging
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Any, List, Optional, Sequence, Tuple

from .pool_objective import PoolConstraints

logger = logging.getLogger(__name__)

# The configurable limits, which are also the `PoolConstraints` field names and
# the params.json keys. One tuple so a key cannot be declared in the schema and
# read here under another name, which is the transposition
# `test_every_key_reaches_the_field_of_the_same_name` exists to catch.
LIMIT_KEYS: Tuple[str, ...] = (
    "min_selectivity_density",
    "max_background_sites",
    "max_worst_hole",
    "max_mean_gap",
    "max_evenness",
    "max_host_coverage",
)

# Every limit `parameter._apply_params_only_keys` assigns from params.json with
# no fallback. `max_dimer_dg` is not a panel limit -- it is a term in the dimer
# screen, which stays outside the objective -- so it rides here for the
# configuration plumbing only and never reaches `PoolConstraints`.
CONFIGURED_LIMIT_KEYS: Tuple[str, ...] = LIMIT_KEYS + (
    "max_dimer_dg",
    "min_per_target_coverage",
)


@dataclass(frozen=True)
class AcceptanceReport:
    """What the limits said about the delivered panel."""

    primers: List[str]
    violations: Tuple[str, ...]
    shortfall: float
    values: Tuple[Tuple[str, Optional[float], float], ...]
    repaired: bool

    def lines(self) -> List[str]:
        """Lines for the CLI, naming each limit, its value and the panel's."""
        if not self.values:
            return []
        head = [
            "",
            "=" * 72,
            "Configured limits",
            "=" * 72,
            f"{'limit':<26} {'panel':>14} {'configured':>14}   met",
            "-" * 62,
        ]
        for name, value, limit in self.values:
            shown = "not measured" if value is None else f"{value:,.6g}"
            met = "yes" if name not in self._failing_fields() else "NO"
            head.append(f"{name:<26} {shown:>14} {limit:>14,.6g}   {met}")
        if self.violations:
            head.append("")
            head.append(f"Not met: {', '.join(self.violations)}.")
            head.append(
                "Repaired: "
                + ("yes" if self.repaired else "no, the panel returned is the one selected")
                + "."
            )
        head.append("=" * 72)
        return head

    def _failing_fields(self) -> set:
        """Field names behind the violation messages, for the `met` column."""
        from .pool_objective import _LIMITS

        by_message = {message: field for field, _m, _s, message, _b in _LIMITS}
        return {by_message[v] for v in self.violations if v in by_message}


def constraints_from_parameter(source: Any) -> Optional[PoolConstraints]:
    """The limits this run configured, or `None` when it configured none.

    `None` rather than an empty `PoolConstraints`, because the caller uses it to
    decide whether to build an objective at all. An empty one would be honoured
    identically and would still cost every panel an evaluation.
    """
    limits = {}
    for key in LIMIT_KEYS:
        value = getattr(source, key, None)
        if value is not None:
            limits[key] = value
    if not limits:
        return None
    return PoolConstraints(**limits)


def _configured_values(
    constraints: PoolConstraints, metrics: Any
) -> Tuple[Tuple[str, Optional[float], float], ...]:
    """Each configured limit beside the panel's own value for it."""
    from .pool_objective import _LIMITS

    rows = []
    for field, metric, _sense, _message, _needs_bg in _LIMITS:
        limit = getattr(constraints, field, None)
        if limit is None:
            continue
        value = getattr(metrics, metric, None)
        rows.append((field, None if value is None else float(value), float(limit)))
    return tuple(rows)


def enforce_constraints(
    primers: Sequence[str],
    optimizer: Any,
    *,
    candidates: Sequence[str],
    constraints: PoolConstraints,
    config: Any = None,
) -> AcceptanceReport:
    """Evaluate the delivered panel, and repair it once if it misses a limit.

    The repair is `pool_planner.repair_panel`, the same bounded swap
    `plan-pool` uses, so there is one repair in the codebase rather than two
    that can disagree. It is skipped when no `config` supplies a budget.

    A repair that does not resolve the violation returns its input. That rule
    is what makes the ordering safe: on a limit no panel can meet, chasing it
    would trade real coverage for a step toward a limit it never reaches, which
    was measured on the Wolbachia pool at an unreachable floor.
    """
    from .panel_refinement import objective_for_optimizer

    objective = objective_for_optimizer(optimizer, constraints)
    panel = list(primers)
    violations = objective.violations(panel)
    repaired = False

    if violations and config is not None:
        try:
            from .pool_planner import repair_panel

            candidate_pool = list(candidates) or panel
            attempt, _record = repair_panel(
                panel, candidate_pool, objective, list(violations), config
            )
            if not objective.violations(attempt):
                panel, repaired = list(attempt), True
        except Exception as exc:  # a repair must never lose the result
            logger.debug(f"Constraint repair skipped: {exc}")

    metrics = objective.metrics(panel)
    return AcceptanceReport(
        primers=panel,
        violations=objective.violations(panel),
        shortfall=objective.shortfall(panel),
        values=_configured_values(constraints, metrics),
        repaired=repaired,
    )


def report_acceptance(report: Optional[AcceptanceReport]) -> None:
    """Print the limits and whether they were met, or nothing at all."""
    if report is None:
        return
    for line in report.lines():
        logger.info(line)


@contextmanager
def _nothing():
    """A no-op scope, so the repair below reads the same with and without a
    budget rather than being written twice."""
    yield


def apply_configured_limits(
    result: Any,
    optimizer: Any,
    *,
    candidates: Sequence[str],
    config: Any,
    constraints: PoolConstraints,
    verbose: bool = False,
    background_available: bool = False,
    budget: Any = None,
) -> Tuple[Optional[Any], "AcceptanceReport"]:
    """Hold one delivered result to its limits, repairing once if it misses.

    Returns `(replacement, report)`. The replacement is a new
    `OptimizationResult` when the repair produced a different panel, and `None`
    when nothing changed, so the caller keeps the object it had rather than
    rebuilding an identical one.

    The report is always returned, and it is the half that says whether the
    delivered panel meets the limits. The replacement cannot: an unchanged
    result means either that nothing needed repairing or that the repair
    failed, and this function used to return the same `None` for both.

    Lives here rather than in `unified_optimizer` because it is acceptance
    logic, and because that module is at its size budget.
    """
    # A limit on a quantity nothing measured would pass every panel, which
    # reads as compliance rather than as an absent measurement. Refuse loudly
    # rather than report a limit that was never evaluated.
    constraints.require_background(background_available)

    from .optimization_service import repair_result
    from .panel_refinement import objective_for_optimizer

    objective = objective_for_optimizer(optimizer, constraints)

    # This repair runs AFTER `_run_panel_stages` returns, so it sat outside the
    # one seam that binds the shared ledger to an evaluator. It can spend up to
    # `swap_max_evaluations` -- 10,000 by default -- and none of it was counted,
    # so a run that declared `total_search_evaluations` got that allowance for
    # the search and a second, undeclared one afterwards.
    #
    # It was never unbounded, so this is the ledger telling the truth rather
    # than runaway cost. With no budget the binding is a no-op and the repair
    # runs exactly as before, which is every default run: the allowance is
    # None unless someone sets it.
    from .search_control import budgeted_objective

    with budgeted_objective(objective, budget) if budget is not None else _nothing():
        updated, details = repair_result(result, optimizer, list(candidates or result.primers))

    # Built on every run, not only a verbose one. This was constructed as an
    # argument to `report_acceptance` inside `if verbose:`, so on a programmatic
    # run the violations were never computed at all and there was nothing to
    # record even in principle. Printing stays verbose-only.
    report = AcceptanceReport(
        list(updated.primers),
        objective.violations(updated.primers),
        objective.shortfall(updated.primers),
        _configured_values(constraints, objective.metrics(updated.primers)),
        bool(details.get("succeeded")),
    )
    if verbose:
        report_acceptance(report)
    if updated is result:
        # NOT "nothing was wrong". The docstring above states the rule that
        # makes these two the same return: a repair that does not resolve the
        # violation returns its input. So an unchanged result arrives here both
        # when the panel met every limit and when it failed one the repair
        # could not fix, and returning a bare None for both is how a panel
        # missing the user's own limit reached the order form.
        return None, report
    from dataclasses import replace

    stage = dict(
        details,
        stage="repair",
        before_size=len(result.primers),
        after_size=len(updated.primers),
        coverage=objective.coverage(updated.primers),
        coverage_metric=objective.constraints.coverage_metric,
        failed_constraints=list(objective.violations(updated.primers)),
        changed=True,
    )
    return replace(updated, stage_history=(*updated.stage_history, stage)), report


@dataclass(frozen=True)
class PerTargetReport:
    """Whether every target cleared the floor, and which did not."""

    floor: float
    coverage: Tuple[Tuple[str, float], ...]
    below: Tuple[str, ...]

    @property
    def met(self) -> bool:
        return not self.below

    @property
    def worst_target(self) -> Optional[str]:
        return self.coverage[-1][0] if self.coverage else None

    @property
    def worst_coverage(self) -> Optional[float]:
        return self.coverage[-1][1] if self.coverage else None

    def lines(self) -> List[str]:
        out = [
            "",
            "=" * 72,
            f"Per-target coverage against a floor of {self.floor:.3g}",
            "=" * 72,
        ]
        for name, value in self.coverage:
            mark = "NO" if name in self.below else "yes"
            out.append(f"{name:<52} {value:>10.4f}   {mark}")
        if self.below:
            out.append("")
            out.append(
                f"Below the floor: {', '.join(self.below)}. Aggregate coverage "
                "can hide this: a panel covering one target well and another "
                "barely beats a balanced one on the mean."
            )
            out.append(
                "Not repaired, and not because nothing could be done: the "
                "repair scores candidate panels through `compute_metrics`, "
                "which does not populate per-target coverage, so chasing this "
                "floor would score every candidate against an empty dict. "
                "Raise the panel size or widen the candidate pool instead."
            )
        out.append("=" * 72)
        return out


def check_per_target_coverage(metrics: Any, floor: Optional[float]) -> Optional[PerTargetReport]:
    """Whether every target cleared `floor`. `None` when there is nothing to say.

    `None` rather than a passing report when no floor is set, when the floor is
    0.0 (how this option has always spelled "disabled"), or when
    `per_target_coverage` is empty -- which is the single-genome case, and a
    floor on an absent measurement must not read as satisfied.
    """
    if not floor:
        return None
    per_target = dict(getattr(metrics, "per_target_coverage", None) or {})
    if not per_target:
        return None
    ordered = tuple(sorted(per_target.items(), key=lambda kv: (-kv[1], kv[0])))
    below = tuple(name for name, value in ordered if value < floor)
    return PerTargetReport(floor=float(floor), coverage=ordered, below=below)


def per_target_floor(args: Any, params: Any) -> Optional[float]:
    """The per-target floor in force: the flag if given, else the config.

    `args` may be an argparse namespace or the kwargs dict `run_optimization`
    receives, which is the same value by another name.

    `--min-per-target-coverage` carried an argparse default of 0.0, which would
    have beaten any configured value on every run once the key existed. That is
    Known Issue 8's shape, so the flag now defaults to `None` and this resolver
    is what distinguishes "not asked" from "asked for zero".
    """
    if isinstance(args, dict):
        flag = args.get("min_per_target_coverage")
    else:
        flag = getattr(args, "min_per_target_coverage", None)
    if flag is not None:
        return float(flag)
    configured = getattr(params, "min_per_target_coverage", None)
    return None if configured is None else float(configured)


def report_per_target(report: Optional[PerTargetReport]) -> None:
    """Print the per-target table, or nothing at all."""
    if report is None:
        return
    for line in report.lines():
        logger.info(line)
