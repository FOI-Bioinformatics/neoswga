"""A bounded beam over partial panels.

Task 6 of the condition-aware pool design plan.

A greedy grows one panel and commits to every choice it makes. That is adequate
while the only thing being maximised is coverage, because coverage is monotonic
under additions and a prefix that looks good stays useful. It is not adequate
under the constraints this design carries.

Selectivity density is a ratio of target to background site density. Adding a
primer that binds the target well and the background little RAISES it, so a
partial panel below the floor may be above it two primers later, and the
locally best first primer can be exactly the one that makes the qualifying
panel unreachable. A greedy has no way back from that.

A beam keeps several partial panels alive instead of one. The cost is bounded
by construction: the beam width times the pool size times the panel size, in
objective evaluations, with an explicit evaluation and wall-clock budget on top.
That bound is the reason this is a beam and not a wider search.

Which violations may be abandoned is decided in :mod:`neoswga.core.partial_panel`
rather than here, so the monotonicity argument lives in one place and both the
greedy and the beam read the same rule.
"""

from __future__ import annotations

import math
import time
from collections.abc import Callable, Sequence
from dataclasses import dataclass

from neoswga.core.partial_panel import can_prune

# A feasible panel of the requested size was found.
TARGET_MET = "target_met"
# The evaluation or wall-clock budget ran out first.
BUDGET_EXHAUSTED = "budget_exhausted"
# No extension survived: the pool ran out, or every way forward was pruned.
# Not a proof that no qualifying panel exists.
INVENTORY_EXHAUSTED = "inventory_exhausted"
# Nothing to search over.
NO_CANDIDATES = "no_candidates"


@dataclass(frozen=True)
class BeamResult:
    """What the search delivered, and how it ended.

    `violations` describes the returned panel. An empty tuple means it satisfies
    every configured constraint; a non-empty one means the search returned the
    best it reached rather than a qualifying panel, and says which limits it
    fails. The two are deliberately separate from `status`: a search can end on
    an exhausted budget having already found a qualifying panel.
    """

    primers: tuple[str, ...]
    status: str
    evaluations: int
    pruned: int
    violations: tuple[str, ...]


def _bg_sites(objective, panel) -> float:
    """Background load, or zero when the evaluator does not report one."""
    return float(getattr(objective.metrics(panel), "total_bg_sites", 0) or 0)


def _rank(objective, panel):
    """Sort key, ascending, best first.

    Constraints come before coverage, which is the same ordering the swap
    refinement uses: a panel nearer to satisfying its limits is preferred
    however much coverage the alternative would buy, because a constraint is
    not a scoring term to be outbid. The panel itself is last, so a complete
    tie resolves to one answer rather than to whichever order the pool happened
    to arrive in.

    Ranked on `shortfall` rather than on the NUMBER of violated constraints.
    Counting ties whenever two panels fail the same single limit, and coverage
    then decides, which is how more search came to move a panel further from
    the floor it was chasing. Zero shortfall is exactly feasibility, so this
    keeps every feasible panel ahead of every infeasible one.
    """
    coverage = objective.coverage(panel)
    if coverage is None:
        coverage = -math.inf
    return (objective.shortfall(panel), -coverage, _bg_sites(objective, panel), panel)


def beam_search(
    candidates: Sequence[str],
    objective,
    size: int,
    *,
    dimerises: Callable[[str, Sequence[str]], bool] | None = None,
    beam_width: int = 4,
    fixed: Sequence[str] = (),
    max_evaluations: int = 10_000,
    max_seconds: float = 10.0,
) -> BeamResult:
    """Grow panels of `size` from `candidates`, keeping `beam_width` alive.

    `dimerises(candidate, selected)` excludes a candidate from a panel outright.
    It is a hard constraint on what may be built, not a term in `_rank`: folding
    it in among the others is how it became tradeable before, and the relaxation
    that followed produced an 11 bp heterodimer against a configured 3.

    `fixed` primers start every panel and are never removed.
    """
    if beam_width < 1:
        raise ValueError("Beam width must be at least 1")
    if not isinstance(size, int) or size < 1:
        raise ValueError("Panel size must be a positive integer")
    if max_evaluations < 0 or not math.isfinite(max_seconds) or max_seconds < 0:
        raise ValueError("Search budgets must be finite and non-negative")

    pool = list(dict.fromkeys(str(c).upper() for c in candidates))
    start = tuple(sorted(str(p).upper() for p in fixed))
    if not pool:
        return BeamResult(start, NO_CANDIDATES, 0, 0, ())

    constraints = objective.constraints
    deadline = time.monotonic() + max_seconds
    beam = [start]
    evaluations = pruned = 0
    best_feasible = None
    best_at_size = None
    status = INVENTORY_EXHAUSTED

    while beam and len(beam[0]) < size:
        seen = set()
        scored = []
        out_of_budget = False
        for prefix in beam:
            for candidate in pool:
                if candidate in prefix:
                    continue
                if dimerises is not None and dimerises(candidate, prefix):
                    continue
                panel = tuple(sorted(prefix + (candidate,)))
                if panel in seen:
                    continue
                seen.add(panel)
                if evaluations >= max_evaluations or time.monotonic() >= deadline:
                    out_of_budget = True
                    break
                evaluations += 1
                violations = objective.violations(panel)
                if can_prune(
                    violations,
                    constraints.max_background_sites,
                    constraints.min_selectivity_density,
                ):
                    pruned += 1
                    continue
                key = _rank(objective, panel)
                scored.append((key, panel))
                if not violations:
                    if best_feasible is None or key < best_feasible[0]:
                        best_feasible = (key, panel)
                    if len(panel) == size and (best_at_size is None or key < best_at_size[0]):
                        best_at_size = (key, panel)
            if out_of_budget:
                break
        if out_of_budget:
            status = BUDGET_EXHAUSTED
            break
        if not scored:
            break
        scored.sort()
        beam = [panel for _, panel in scored[:beam_width]]

    if best_at_size is not None:
        status = TARGET_MET if status != BUDGET_EXHAUSTED else status
        chosen = best_at_size[1]
    elif best_feasible is not None:
        chosen = best_feasible[1]
    elif beam:
        chosen = beam[0]
    else:
        chosen = start

    return BeamResult(
        primers=chosen,
        status=status,
        evaluations=evaluations,
        pruned=pruned,
        violations=tuple(objective.violations(chosen)) if chosen else (),
    )
