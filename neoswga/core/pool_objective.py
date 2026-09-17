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


@dataclass(frozen=True)
class PoolConstraints:
    """What a panel must satisfy to be accepted.

    Frozen: a design must not have its acceptance criteria changed underneath it
    while it runs, and a constraint set is part of the identity of a result.
    """

    coverage_metric: str = "effective"
    min_selectivity_density: Optional[float] = None
    max_background_sites: Optional[int] = None

    def __post_init__(self) -> None:
        if self.coverage_metric not in {"effective", "raw"}:
            raise ValueError(
                f"coverage_metric must be 'effective' or 'raw', " f"not {self.coverage_metric!r}"
            )

    @property
    def needs_background(self) -> bool:
        return self.min_selectivity_density is not None or self.max_background_sites is not None

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
                "Supply a background, or remove min_selectivity_density and "
                "max_background_sites."
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

        floor = self.constraints.min_selectivity_density
        if floor is not None and metrics.selectivity_density < floor:
            reasons.append("selectivity below minimum")

        ceiling = self.constraints.max_background_sites
        if ceiling is not None and metrics.total_bg_sites > ceiling:
            reasons.append("background sites above maximum")

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

        total = 0.0
        floor = self.constraints.min_selectivity_density
        if floor is not None and metrics.selectivity_density < floor:
            # Scaled by the limit itself, and guarded for a zero limit so the
            # boundary keeps agreeing with `violations` rather than dividing.
            total += (floor - metrics.selectivity_density) / (abs(floor) or 1.0)

        ceiling = self.constraints.max_background_sites
        if ceiling is not None and metrics.total_bg_sites > ceiling:
            total += (metrics.total_bg_sites - ceiling) / max(abs(ceiling), 1)

        return total
