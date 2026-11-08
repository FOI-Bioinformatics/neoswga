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
from collections.abc import Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Any

from .pool_objective import PoolConstraints

logger = logging.getLogger(__name__)

# The configurable limits, which are also the `PoolConstraints` field names and
# the params.json keys. One tuple so a key cannot be declared in the schema and
# read here under another name, which is the transposition
# `test_every_key_reaches_the_field_of_the_same_name` exists to catch.
LIMIT_KEYS: tuple[str, ...] = (
    "min_selectivity_density",
    "max_background_sites",
    "max_worst_hole",
    "max_mean_gap",
    "max_evenness",
    "max_host_coverage",
    # Joined the enforced limits on 2026-10-02. It was in `CONFIGURED_LIMIT_KEYS`
    # alone -- configured, checked and reported, and deliberately not repaired,
    # because the repair scores candidate panels through `compute_metrics` and
    # `compute_metrics` did not populate per-target coverage. It does now, and
    # `worst_target_coverage` is the figure this limit is read against. See
    # `core/target_coverage.py` and Phase 6c of the 2026-10-02 diversity plan.
    "min_per_target_coverage",
)

# Every limit `parameter._apply_params_only_keys` assigns from params.json with
# no fallback. `max_dimer_dg` is not a panel limit -- it is a term in the dimer
# screen, which stays outside the objective -- so it rides here for the
# configuration plumbing only and never reaches `PoolConstraints`.
CONFIGURED_LIMIT_KEYS: tuple[str, ...] = LIMIT_KEYS + ("max_dimer_dg",)

# Keys where 0.0 means "disabled" rather than "a floor of zero". One key has
# ever been spelled that way: `--min-per-target-coverage 0` has meant off since
# before it was a limit, and reading it as a floor would build an objective for
# a run that asked for nothing -- which is the one thing that would move a
# delivered panel unasked.
DISABLED_WHEN_ZERO: frozenset[str] = frozenset({"min_per_target_coverage"})


@dataclass(frozen=True)
class AcceptanceReport:
    """What the limits said about the delivered panel."""

    primers: list[str]
    violations: tuple[str, ...]
    shortfall: float
    values: tuple[tuple[str, float | None, float], ...]
    repaired: bool
    #: Which reading of "the background" the three background-measured limits
    #: were judged under. Last, with a default, so neither existing positional
    #: construction moves. A table of limits that does not say this leaves a
    #: reader to assume a panel was held to each host when it was held to the
    #: pool, or the reverse.
    background_aggregate: str = "sum"

    def lines(self) -> list[str]:
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
        head.append("")
        head.append(f"background_aggregate: {self.background_aggregate}.")
        if self.background_aggregate != "sum":
            # The `panel` column above is then the WORST HOST's figure for each
            # background-measured limit, reduced per limit. Said here because
            # a number in a column reads as a property of the panel against the
            # whole background, and under this reading it is not.
            head.append(
                "  The panel column for min_selectivity_density, "
                "max_background_sites and max_host_coverage is the WORST HOST's "
                "figure, each reduced over hosts separately. The pooled figures "
                "are in the summary under their own names."
            )
        if any(name == "min_per_target_coverage" for name, _v, _l in self.values):
            # The one row whose panel figure is a reduction over several
            # measurements rather than one. Said here because "0.31" in a column
            # reads as a property of the panel, and this one is a minimum.
            head.append("")
            head.append(
                "min_per_target_coverage is read against worst_target_coverage: "
                "the MINIMUM over the foreground targets, shown as not measured "
                "when any target could not be measured, which does not pass."
            )
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


def _flag(flags: Any, key: str) -> Any:
    """One limit as this invocation gave it, or `None` for "not asked".

    `flags` is an argparse namespace or the kwargs dict `run_optimization`
    receives, which is the same value by another name.
    """
    if flags is None:
        return None
    if isinstance(flags, dict):
        return flags.get(key)
    return getattr(flags, key, None)


def constraints_from_parameter(source: Any, flags: Any = None) -> PoolConstraints | None:
    """The limits this run configured, or `None` when it configured none.

    `None` rather than an empty `PoolConstraints`, because the caller uses it to
    decide whether to build an objective at all. An empty one would be honoured
    identically and would still cost every panel an evaluation.

    `flags` is this invocation's CLI namespace or kwargs, consulted before
    `source`. `min_per_target_coverage` is the only panel limit that also has a
    CLI flag, and a flag the enforcing code never reads is Known Issue 8's
    shape: the floor would print a per-target table and not be repaired, while
    the same value in params.json would be. An absent flag is `None` and so does
    not beat the file.
    """
    limits = {}
    for key in LIMIT_KEYS:
        value = _flag(flags, key)
        if value is None:
            value = getattr(source, key, None)
        if value is None:
            continue
        if key in DISABLED_WHEN_ZERO and not float(value):
            continue
        limits[key] = value
    if not limits:
        return None
    return PoolConstraints(**limits)


def _configured_values(
    constraints: PoolConstraints, metrics: Any
) -> tuple[tuple[str, float | None, float], ...]:
    """Each configured limit beside the figure it was actually judged on.

    `per_host_metrics.limit_figure`, not a plain attribute lookup. The `met`
    column of the table these rows fill comes from `violations`, which asks
    `limit_figure`; reading the pooled attribute here put a different number in
    the `panel` column, and under `worst-case` the two disagreed. The printed
    symptom was `min_selectivity_density  51  10  NO` -- a value comfortably
    above its floor, marked not met -- with nothing in the table saying the
    verdict had been reached on a different figure.

    One dispatch for deciding and for displaying. A value and a verdict from
    different arithmetic is the defect `PoolObjective` exists to prevent, and
    this is the same defect one layer out.
    """
    from .per_host_metrics import limit_figure
    from .pool_objective import _LIMITS

    rows = []
    for field, metric, _sense, _message, _needs_bg in _LIMITS:
        limit = getattr(constraints, field, None)
        if limit is None:
            continue
        value = limit_figure(metrics, metric)
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
        background_aggregate=str(getattr(metrics, "background_aggregate", None) or "sum"),
    )


def report_acceptance(report: AcceptanceReport | None) -> None:
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
) -> tuple[Any | None, AcceptanceReport]:
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
    _final = objective.metrics(updated.primers)
    report = AcceptanceReport(
        list(updated.primers),
        objective.violations(updated.primers),
        objective.shortfall(updated.primers),
        _configured_values(constraints, _final),
        bool(details.get("succeeded")),
        str(getattr(_final, "background_aggregate", None) or "sum"),
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
    coverage: tuple[tuple[str, float], ...]
    below: tuple[str, ...]

    @property
    def met(self) -> bool:
        return not self.below

    @property
    def worst_target(self) -> str | None:
        return self.coverage[-1][0] if self.coverage else None

    @property
    def worst_coverage(self) -> float | None:
        return self.coverage[-1][1] if self.coverage else None

    def lines(self) -> list[str]:
        out = [
            "",
            "=" * 72,
            f"Per-target coverage against a floor of {self.floor:.3g}",
            "=" * 72,
        ]
        for name, value in self.coverage:
            mark = "NO" if name in self.below else "yes"
            out.append(f"{name:<52} {value:>10.4f}   {mark}")
        out.append("")
        out.append(
            "The floor is read against worst_target_coverage, the MINIMUM over "
            "the rows above. A target whose coverage could not be measured is "
            "not a row here and makes that minimum unknown, in which case the "
            "limit is refused rather than reported as met."
        )
        if self.below:
            out.append("")
            out.append(
                f"Below the floor: {', '.join(self.below)}. Aggregate coverage "
                "can hide this: a panel covering one target well and another "
                "barely beats a balanced one on the mean."
            )
            out.append(
                "This panel is what one bounded repair returned, or what was "
                "selected when the repair could not meet the floor. The "
                "'Configured limits' table above says which. Raising the panel "
                "size or widening the candidate pool is the other lever."
            )
        out.append("=" * 72)
        return out


def check_per_target_coverage(metrics: Any, floor: float | None) -> PerTargetReport | None:
    """Whether every target cleared `floor`. `None` when there is nothing to say.

    `None` rather than a passing report when no floor is set, when the floor is
    0.0 (how this option has always spelled "disabled"), or when
    `per_target_coverage` is empty, because a floor on an absent measurement
    must not read as satisfied.

    This is a REPORTER, and that last case is why it cannot be the only guard:
    "nothing to say" and "every target cleared the floor" print the same, which
    is nothing. `target_coverage.require_per_target` refuses a requested
    floor that nothing measured, and `pool_objective._unmet` refuses one whose
    `worst_target_coverage` is `None`; both run before this is called.
    """
    if not floor:
        return None
    per_target = dict(getattr(metrics, "per_target_coverage", None) or {})
    if not per_target:
        return None
    ordered = tuple(sorted(per_target.items(), key=lambda kv: (-kv[1], kv[0])))
    below = tuple(name for name, value in ordered if value < floor)
    return PerTargetReport(floor=float(floor), coverage=ordered, below=below)


def per_target_floor(args: Any, params: Any) -> float | None:
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


def report_per_target(report: PerTargetReport | None) -> None:
    """Print the per-target table, or nothing at all."""
    if report is None:
        return
    for line in report.lines():
        logger.info(line)
