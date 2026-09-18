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

**Every limit is unset by default and this module is inert until one is set.**
That is deliberate rather than cautious. No spacing threshold derived from the
polymerase reach separates the 18 published sets with wet-lab outcomes, the
winners included, and a fitted weight is wrong for one of the two benchmarks
either way; so NeoSWGA must not pick a limit, and a limit that changed a
delivered panel unasked would be exactly the scoring change that evidence
refuses. `constraints_from_parameter` returns `None` when nothing is configured,
and a run that asks for nothing acquires no objective at all. See
`docs/validation/getting_ahead_on_spacing_2026-09-18.md`.

The dimer guarantee stays outside, as it does in `PoolObjective`: it is a hard
constraint on the delivered panel rather than a scoring term, and folding it in
among tradeable terms is how it came to be traded.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, List, Optional, Sequence, Tuple

from .pool_objective import PoolConstraints, PoolObjective

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
CONFIGURED_LIMIT_KEYS: Tuple[str, ...] = LIMIT_KEYS + ("max_dimer_dg",)


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
    objective = PoolObjective(optimizer.compute_metrics, constraints)
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


def apply_configured_limits(
    result: Any,
    optimizer: Any,
    *,
    candidates: Sequence[str],
    config: Any,
    constraints: PoolConstraints,
    verbose: bool = False,
    background_available: bool = False,
) -> Optional[Any]:
    """Hold one delivered result to its limits, repairing once if it misses.

    Returns a replacement `OptimizationResult` when the repair produced a
    different panel, and `None` when nothing changed, so the caller keeps the
    object it had rather than rebuilding an identical one.

    Lives here rather than in `unified_optimizer` because it is acceptance
    logic, and because that module is at its size budget.
    """
    from dataclasses import replace as _dc_replace

    # A limit on a quantity nothing measured would pass every panel, which
    # reads as compliance rather than as an absent measurement. Refuse loudly
    # rather than report a limit that was never evaluated.
    constraints.require_background(background_available)

    report = enforce_constraints(
        list(result.primers),
        optimizer,
        candidates=list(candidates or []),
        constraints=constraints,
        config=config,
    )
    if verbose:
        report_acceptance(report)

    if list(report.primers) == list(result.primers):
        return None
    return _dc_replace(
        result,
        primers=tuple(report.primers),
        metrics=optimizer.compute_metrics(report.primers),
    )
