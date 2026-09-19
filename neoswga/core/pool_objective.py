"""One definition of what a panel is worth, shared by search and acceptance.

Task 5 of the condition-aware pool design plan.

The optimizers selected on one quantity and the report accepted on another. A
greedy maximising raw window coverage can prefer a panel the occupancy-weighted
metric scores lower, because occupancy weights each site by how much of the time
it is actually bound. When the number a design is accepted on is not the number
it was chosen on, improving the search does not reliably improve the result, and
a comparison between two designs measures the disagreement as much as the
designs.

`PoolObjective` is the contract. It wraps whatever evaluator computes a panel's
metrics -- `BaseOptimizer.compute_metrics` in the pipeline -- and answers three
questions from one set of definitions: what the metrics are, which coverage
figure counts, and which constraints a panel violates.

Two things are deliberately outside it.

The dimer guarantee is a hard constraint on the delivered panel, enforced
separately. Folding it in among scoring terms is how it came to be tradeable,
and the relaxation that followed produced an 11 bp heterodimer against a
configured 3.

Coverage here remains the occupancy-weighted proxy documented in
`base_optimizer._compute_effective_coverage`. Sharing a definition between search
and acceptance makes them consistent; it does not make either one a measurement
of genome recovery.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Callable, Dict, Optional, Sequence, Tuple

# Sense of each limit.
_AT_LEAST = "at_least"
_AT_MOST = "at_most"

# Every limit a panel can be held to, in one table, so `violations` and
# `shortfall` cannot drift apart. Each row is
# (constraint field, metrics field, sense, message, needs a background).
#
# The order is the order `violations` reports in, and the first two rows keep
# their original position because `pool_planner._repair` takes `reasons[0]` as
# the reason it records.
#
# Two quantities are deliberately absent, and the reason changed on
# 2026-09-19. `strand_coverage_ratio` and `strand_alternation_score` used to
# read 0.0 both when measured zero and when they could not be measured, so a
# limit on either would have rejected a panel for a MISSING MEASUREMENT while
# reporting a violated constraint -- the failure Known Issues 5, 6 and 13
# record. That ambiguity is now fixed at the source: alternation is None below
# two sites, where there is no adjacent pair to score, and the ratio is None
# only with no sites at all.
#
# They stay out because no threshold has a reference. That is the same reason
# `max_worst_hole` and `max_host_coverage` ship unset rather than defaulted,
# and adding a limit nobody can choose a value for would be worse than the
# reporting the regime table already gives. Adding them is now a decision
# about evidence rather than a blocked repair.
# The dimer limit is absent for a different reason: it is a hard constraint on
# the delivered panel rather than a scoring term, and folding it in among these
# is how it became tradeable.
_LIMITS = (
    (
        "min_selectivity_density",
        "selectivity_density",
        _AT_LEAST,
        "selectivity below minimum",
        True,
    ),
    ("max_background_sites", "total_bg_sites", _AT_MOST, "background sites above maximum", True),
    ("max_worst_hole", "max_gap", _AT_MOST, "worst hole above maximum", False),
    ("max_mean_gap", "mean_gap", _AT_MOST, "mean gap above maximum", False),
    ("max_evenness", "gap_gini", _AT_MOST, "evenness above maximum", False),
    ("max_host_coverage", "bg_coverage", _AT_MOST, "host coverage above maximum", True),
)


def _unmet(value, limit, sense) -> bool:
    """Whether this value misses its limit. An unmeasurable value misses it.

    `None` means the evaluator could not produce the quantity, and a limit on
    something nothing measured must not pass: every panel would clear it, which
    reads as compliance rather than as an absent measurement.
    """
    if value is None:
        return True
    if sense == _AT_LEAST:
        return value < limit
    return value > limit


def _distance(value, limit, sense) -> float:
    """How far this value is from its limit, relative to the limit itself.

    Scaled so a hole over its ceiling by half scores the same as evenness over
    its ceiling by half, and neither dominates merely by being measured in
    larger units. An unmeasurable value is infinitely far rather than large.
    """
    if value is None:
        return math.inf
    scale = abs(limit) or 1.0
    gap = (limit - value) if sense == _AT_LEAST else (value - limit)
    return gap / scale


@dataclass(frozen=True)
class PoolConstraints:
    """What a panel must satisfy to be accepted.

    Frozen: a design must not have its acceptance criteria changed underneath it
    while it runs, and a constraint set is part of the identity of a result.
    """

    coverage_metric: str = "effective"
    min_selectivity_density: Optional[float] = None
    max_background_sites: Optional[int] = None

    # The spacing quantities, added 2026-09-18. Every one is unset by default
    # and so inert: NeoSWGA must not pick a spacing threshold, because none
    # derived from the reach separates the 18 published sets with wet-lab
    # outcomes, and a fitted weight is wrong for one of the two benchmarks
    # either way. A user setting one is a different claim, and item 1's report
    # is what tells them these had no reference at all. See
    # `docs/validation/getting_ahead_on_spacing_2026-09-18.md`.
    max_worst_hole: Optional[float] = None
    max_mean_gap: Optional[float] = None
    max_evenness: Optional[float] = None
    max_host_coverage: Optional[float] = None

    def __post_init__(self) -> None:
        if self.coverage_metric not in {"effective", "raw"}:
            raise ValueError(
                f"coverage_metric must be 'effective' or 'raw', " f"not {self.coverage_metric!r}"
            )

    @property
    def needs_background(self) -> bool:
        """Whether any limit here is measured against the background genome."""
        return any(
            getattr(self, field) is not None
            for field, _metric, _sense, _message, needs_bg in _LIMITS
            if needs_bg
        )

    def require_background(self, available: bool) -> None:
        """Refuse a specificity limit when nothing measured the background.

        A limit on a quantity nothing measured can be neither satisfied nor
        refused; every panel would pass, which reads as specificity rather than
        as an absent measurement.
        """
        if self.needs_background and not available:
            raise ValueError(
                "A specificity limit was set but no background genome and index "
                "are available. Without them the limit cannot be evaluated, and "
                "every panel would pass as if it bound nothing off-target. "
                "Supply a background, or remove the limit: "
                + ", ".join(
                    field
                    for field, _m, _s, _msg, needs_bg in _LIMITS
                    if needs_bg and getattr(self, field) is not None
                )
                + "."
            )


class PoolObjective:
    """The metrics and constraints a design is both searched and accepted on."""

    def __init__(
        self,
        evaluate: Callable[[Sequence[str]], Any],
        constraints: PoolConstraints,
        cache_size: int = 4096,
    ):
        self._evaluate = evaluate
        self.constraints = constraints
        self._cache: Dict[Tuple[str, ...], Any] = {}
        self._cache_size = int(cache_size)

    def metrics(self, primers: Sequence[str]) -> Any:
        """The evaluator's metrics for one panel, computed once.

        Keyed on the sorted panel: search asks about the same set repeatedly and
        in different orders, and a panel is the same panel however it was built.

        Bounded rather than unbounded, because a long search over a large
        inventory would otherwise retain a metrics object for every panel it
        ever considered.
        """
        key = tuple(sorted(str(p).upper() for p in primers))
        if key not in self._cache:
            if len(self._cache) >= self._cache_size:
                self._cache.clear()
            self._cache[key] = self._evaluate(list(primers))
        return self._cache[key]

    def coverage(self, primers: Sequence[str]) -> Optional[float]:
        """The coverage figure this design is judged on.

        `None` when the evaluator could not compute it, which is not the same as
        zero: no reaction conditions means no temperature at which to evaluate
        occupancy, and a fabricated number there would be indistinguishable from
        a measured one.
        """
        metrics = self.metrics(primers)
        if self.constraints.coverage_metric == "raw":
            return metrics.fg_coverage
        return metrics.effective_fg_coverage

    def violations(self, primers: Sequence[str]) -> Tuple[str, ...]:
        """Every constraint this panel fails, named.

        A tuple rather than a boolean, so a caller can report which limit bound
        rather than only that the panel was rejected.
        """
        metrics = self.metrics(primers)
        reasons = []

        if self.coverage(primers) is None:
            reasons.append("coverage unavailable")

        for field, metric, sense, message, _needs_bg in _LIMITS:
            limit = getattr(self.constraints, field)
            if limit is None:
                continue
            if _unmet(getattr(metrics, metric, None), limit, sense):
                reasons.append(message)

        return tuple(reasons)

    def shortfall(self, primers: Sequence[str]) -> float:
        """How far this panel is from satisfying its constraints. 0 when it does.

        `violations` names what failed; this measures by how much, which is what
        a search needs when every candidate fails. The two agree exactly at the
        boundary: the shortfall is zero precisely when `violations` is empty,
        and `tests/test_the_objective_ranks_by_how_far_it_missed.py` pins that
        across the cases. Without that agreement a feasible panel could rank
        behind an infeasible one, which is what the ordering exists to prevent.

        Why it is needed. Both searches ranked on `len(violations)`, so two
        panels failing the SAME single constraint tied, coverage broke the tie,
        and the deciding metric was free to drift. Measured on the real
        Wolbachia pool, widening the candidate frontier then moved a panel away
        from the selectivity floor it was chasing: density 20.9 to 14.6 while
        coverage rose 0.740 to 0.765. More search made the answer worse. See
        `docs/validation/frontier_refill_2026-09-17.md`.

        Each term is RELATIVE to its own limit, so a density floor and a
        background-site ceiling are comparable and neither dominates merely by
        being measured on a larger scale. Missing a floor of 100 by half scores
        the same 0.5 as exceeding a ceiling of 10 by half. The terms are summed,
        so failing two limits is worse than failing one.

        An unmeasurable coverage is infinite rather than large. It is not a
        distance from feasibility: a panel that cannot be scored is not nearly
        acceptable, and a finite value would let coverage trade against it.
        """
        metrics = self.metrics(primers)
        if self.coverage(primers) is None:
            return math.inf

        # One loop over the same table `violations` reads, which is what keeps
        # the boundary agreement exact: a term is added exactly when `_unmet`
        # says the limit is missed, and `_distance` is positive there.
        total = 0.0
        for field, metric, sense, _message, _needs_bg in _LIMITS:
            limit = getattr(self.constraints, field)
            if limit is None:
                continue
            value = getattr(metrics, metric, None)
            if _unmet(value, limit, sense):
                total += _distance(value, limit, sense)

        return total
