"""Missing the coverage target is a reason to try again.

Phase 2 of the plan for `docs/validation/pipeline_audit_2026-09-16/`, audit
finding F2.

`plan_pool` repaired a row only when `PoolObjective.violations` returned
something. That function knows the selectivity floor, the background cap and
unavailable coverage. It does not know `coverage_targets`, which are applied
later as a filter over finished rows. So the one failure a user is most likely
to care about -- "I asked for 90% coverage and did not get it" -- was the one
failure that never triggered the repair. The row was reported `not_found` while
the machinery that could have fixed it sat unused.

The distinction this draws, and keeps drawing, is between a constraint and a
target. A panel below the requested coverage is still deliverable; a panel over
the background cap is not. So a target miss is a reason to ATTEMPT a repair, and
never a violation to be traded against one. The repair objective is unchanged:
fewer violations first, then coverage, then background load. A feasible panel is
not swapped for an infeasible one to buy coverage, whatever the target says.
"""

from types import SimpleNamespace

import pytest

from neoswga.core.base_optimizer import OptimizationStatus, OptimizerConfig
from neoswga.core.pool_planner import plan_pool

GOOD = "ACGGACGGACGG"
POOR = "ACACACACACAC"
RICH = "AGGAGGAGGAGG"


class _Optimizer:
    """Metrics keyed on the panel's contents, so a swap changes the answer."""

    name = "test"
    bg_prefixes = ["bg"]
    bg_seq_lengths = [1000]
    conditions = object()

    def __init__(self, panel, table):
        self.panel = list(panel)
        self.table = table
        self.config = OptimizerConfig(max_dimer_bp=3, max_self_dimer_bp=4)

    def optimize(self, candidates, target_size):  # noqa: ARG002
        return SimpleNamespace(primers=self.panel, status=OptimizationStatus.PARTIAL, message="")

    def compute_metrics(self, primers):
        coverage, density, background = self.table[tuple(sorted(primers))]
        return SimpleNamespace(
            effective_fg_coverage=coverage,
            fg_coverage=coverage,
            selectivity_density=density,
            total_bg_sites=background,
            max_gap=100,
        )


def _table(rows):
    return {tuple(sorted(k)): v for k, v in rows.items()}


def test_a_feasible_panel_below_the_target_is_repaired():
    """The audit's reproduction.

    The optimizer returns a panel that satisfies every constraint and covers
    half the genome. Another panel of the same size covers 95% and satisfies the
    same constraints. Asked for 90%, the planner used to report `not_found`
    without ever attempting a repair.
    """
    table = _table(
        {
            (GOOD, POOR): (0.50, 90.0, 4),
            (GOOD, RICH): (0.95, 80.0, 6),
            (POOR, RICH): (0.40, 85.0, 5),
        }
    )
    opt = _Optimizer([GOOD, POOR], table)

    plan = plan_pool(opt, [GOOD, POOR, RICH], [2], [0.90], min_selectivity_density=10)
    row = plan["rows"][0]

    assert row["repair"]["attempted"] is True, "a target miss did not reach the repair"
    assert set(row["primers"]) == {GOOD, RICH}
    assert row["coverage"] == 0.95
    assert plan["recommendations"][0]["status"] == "smallest_found"


def test_the_row_says_the_target_was_why_it_tried():
    """A reader should be able to tell a constraint repair from a target repair."""
    table = _table(
        {
            (GOOD, POOR): (0.50, 90.0, 4),
            (GOOD, RICH): (0.95, 80.0, 6),
            (POOR, RICH): (0.40, 85.0, 5),
        }
    )
    opt = _Optimizer([GOOD, POOR], table)

    row = plan_pool(opt, [GOOD, POOR, RICH], [2], [0.90], min_selectivity_density=10)["rows"][0]

    assert row["repair"]["reason"] == "below coverage target"


def test_a_panel_already_over_the_target_is_left_alone():
    """No target miss, no constraint violation, no repair."""
    table = _table({(GOOD, RICH): (0.95, 80.0, 6)})
    opt = _Optimizer([GOOD, RICH], table)

    row = plan_pool(opt, [GOOD, RICH], [2], [0.90], min_selectivity_density=10)["rows"][0]

    assert row["repair"]["attempted"] is False


def test_a_target_miss_never_buys_coverage_with_a_constraint():
    """The ordering that must not move.

    The only higher-coverage panel available breaks the selectivity floor. A
    repair chasing the target must refuse it: a target is a request, a floor is
    a requirement, and the row is deliverable below the target but not below the
    floor.
    """
    # Singletons are listed because the swap cannot resolve this row, so the
    # beam runs and builds panels a prefix at a time. A table missing them would
    # fail on a KeyError rather than on the property under test.
    table = _table(
        {
            (GOOD,): (0.30, 95.0, 2),
            (POOR,): (0.20, 95.0, 2),
            (RICH,): (0.25, 95.0, 2),
            (GOOD, POOR): (0.50, 90.0, 4),
            (GOOD, RICH): (0.99, 2.0, 400),
            (POOR, RICH): (0.40, 85.0, 5),
        }
    )
    opt = _Optimizer([GOOD, POOR], table)

    row = plan_pool(opt, [GOOD, POOR, RICH], [2], [0.90], min_selectivity_density=10)["rows"][0]

    assert set(row["primers"]) == {GOOD, POOR}, "the repair traded the floor for coverage"
    assert row["eligible"] is True
    assert row["repair"]["succeeded"] is False


def test_the_highest_requested_target_is_the_one_that_counts():
    """Several targets are requested at once; the most demanding drives repair.

    Repairing to the lowest would stop improving a panel as soon as it cleared
    the easiest request, and the higher rows would report `not_found` beside a
    panel that was never asked to reach them.
    """
    table = _table(
        {
            (GOOD, POOR): (0.60, 90.0, 4),
            (GOOD, RICH): (0.95, 80.0, 6),
            (POOR, RICH): (0.40, 85.0, 5),
        }
    )
    opt = _Optimizer([GOOD, POOR], table)

    plan = plan_pool(opt, [GOOD, POOR, RICH], [2], [0.5, 0.9], min_selectivity_density=10)
    row = plan["rows"][0]

    assert row["repair"]["attempted"] is True, "0.60 cleared 0.5, so nothing was attempted"
    assert row["coverage"] == 0.95
    assert [r["status"] for r in plan["recommendations"]] == ["smallest_found", "smallest_found"]


def test_repair_off_still_reports_the_miss_without_acting():
    table = _table(
        {
            (GOOD, POOR): (0.50, 90.0, 4),
            (GOOD, RICH): (0.95, 80.0, 6),
            (POOR, RICH): (0.40, 85.0, 5),
        }
    )
    opt = _Optimizer([GOOD, POOR], table)

    plan = plan_pool(opt, [GOOD, POOR, RICH], [2], [0.90], min_selectivity_density=10, repair=False)

    assert plan["rows"][0]["repair"]["attempted"] is False
    assert plan["recommendations"][0]["status"] == "not_found"


def test_a_constraint_violation_still_reports_itself_as_the_reason():
    """The pre-existing trigger keeps its own label."""
    table = _table(
        {
            (GOOD, POOR): (0.95, 2.0, 400),
            (GOOD, RICH): (0.94, 80.0, 6),
            (POOR, RICH): (0.40, 85.0, 5),
        }
    )
    opt = _Optimizer([GOOD, POOR], table)

    row = plan_pool(opt, [GOOD, POOR, RICH], [2], [0.90], min_selectivity_density=10)["rows"][0]

    assert row["repair"]["reason"] == "selectivity below minimum"
    assert set(row["primers"]) == {GOOD, RICH}


@pytest.mark.parametrize("target", [0.30, 0.49])
def test_a_target_below_the_delivered_coverage_is_not_a_miss(target):
    table = _table({(GOOD, POOR): (0.50, 90.0, 4)})
    opt = _Optimizer([GOOD, POOR], table)

    row = plan_pool(opt, [GOOD, POOR], [2], [target], min_selectivity_density=10)["rows"][0]

    assert row["repair"]["attempted"] is False
