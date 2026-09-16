"""A panel that misses a repairable limit gets a bounded second attempt.

Task 6 of the condition-aware pool design plan.

`plan_pool` used to evaluate the optimizer's panel once and, if it missed a
constraint, record the row as ineligible and move on. The swap refinement that
prioritises reducing violations already existed and nothing in this path called
it, so a panel one swap away from qualifying was reported as not found.

Repair is bounded and conservative. It keeps the panel size, accepts only
strictly improving swaps under the same objective the row is judged on, and
re-evaluates the result through that objective rather than trusting the swap
loop's own bookkeeping. A panel that already qualifies is not touched.

Dimer violations are deliberately not repaired here: the objective does not see
them, so the swap score cannot be steered by them, and a panel that reaches this
point with a dimerising pair means an upstream relaxation fired. That is a
problem to fix where it happens, not to paper over at reporting time.
"""

from types import SimpleNamespace

from neoswga.core.base_optimizer import OptimizationStatus, OptimizerConfig
from neoswga.core.pool_planner import plan_pool

GOOD = "ACGGACGGACGG"
BAD = "ACACACACACAC"
SPARE = "AGGAGGAGGAGG"
# The only pair among these that the configured limit refuses.
HOOK = "AAGAAGAAGAAG"
BARB = "CTTCTTCTTCTT"


class _Optimizer:
    """Metrics depend on WHICH primers are in the panel, not only how many.

    The stub in `test_pool_planner.py` keys on panel size, which is enough for
    the questions asked there and useless here: a repair swap changes the
    composition and not the size, so a size-keyed stub would report the panel as
    unchanged whatever the swap did.
    """

    name = "test"
    bg_prefixes = ["bg"]
    bg_seq_lengths = [1000]
    conditions = object()

    def __init__(self, panel, table, max_dimer_bp=3):
        self.panel = list(panel)
        self.table = table
        self.config = OptimizerConfig(max_dimer_bp=max_dimer_bp, max_self_dimer_bp=4)
        self.evaluated = []

    def optimize(self, candidates, target_size):  # noqa: ARG002
        return SimpleNamespace(primers=self.panel, status=OptimizationStatus.PARTIAL, message="")

    def compute_metrics(self, primers):
        key = tuple(sorted(primers))
        self.evaluated.append(key)
        coverage, density, bg = self.table[key]
        return SimpleNamespace(
            effective_fg_coverage=coverage,
            fg_coverage=coverage,
            selectivity_density=density,
            total_bg_sites=bg,
            max_gap=100,
        )


def _table(rows, singletons=None):
    """Panels keyed on their sorted contents.

    Singletons must be listed too: the repair beam builds by addition, so it
    evaluates every prefix on the way to the requested size. A table that
    supplied a default for a missing panel would let a test pass because the
    default happened to agree, which is the fixture failure this project has
    met before.
    """
    table = {tuple(sorted(k)): v for k, v in rows.items()}
    for primer, value in (singletons or {}).items():
        table[(primer,)] = value
    return table


def test_a_panel_one_swap_short_of_qualifying_is_repaired():
    """The case the repair exists for."""
    table = _table(
        {
            (GOOD, BAD): (0.90, 2.0, 40),
            (GOOD, SPARE): (0.88, 30.0, 5),
            (BAD, SPARE): (0.50, 25.0, 6),
        }
    )
    opt = _Optimizer([GOOD, BAD], table)

    plan = plan_pool(opt, [GOOD, BAD, SPARE], [2], [0.8], min_selectivity_density=10)
    row = plan["rows"][0]

    assert set(row["primers"]) == {GOOD, SPARE}
    assert row["eligible"] is True
    assert row["repair"]["succeeded"] is True
    assert row["repair"]["swaps"] == 1
    assert row["coverage"] == 0.88, "the row must describe the panel it delivered"
    assert plan["recommendations"][0]["size"] == 2


def test_a_qualifying_panel_is_left_alone():
    """No repair, and no evaluations spent looking for one."""
    table = _table({(GOOD, SPARE): (0.88, 30.0, 5)})
    opt = _Optimizer([GOOD, SPARE], table)

    plan = plan_pool(opt, [GOOD, SPARE], [2], [0.8], min_selectivity_density=10)
    row = plan["rows"][0]

    assert row["repair"]["attempted"] is False
    assert set(row["primers"]) == {GOOD, SPARE}
    assert opt.evaluated.count((GOOD, SPARE)) >= 1


def test_an_unrepairable_panel_keeps_the_panel_it_had():
    """A failed repair must not deliver something worse than it started with."""
    table = _table(
        {
            (GOOD, BAD): (0.90, 2.0, 40),
            (GOOD, SPARE): (0.40, 3.0, 45),
            (BAD, SPARE): (0.30, 1.0, 50),
        },
        singletons={GOOD: (0.30, 2.0, 20), BAD: (0.20, 1.0, 25), SPARE: (0.10, 1.5, 30)},
    )
    opt = _Optimizer([GOOD, BAD], table)

    plan = plan_pool(opt, [GOOD, BAD, SPARE], [2], [0.8], min_selectivity_density=10)
    row = plan["rows"][0]

    assert set(row["primers"]) == {GOOD, BAD}
    assert row["eligible"] is False
    assert row["repair"]["attempted"] is True
    assert row["repair"]["succeeded"] is False
    assert row["repair"]["beam"] == "inventory_exhausted", "the beam was not reached"
    assert row["failed_constraints"] == ["selectivity below minimum"]


def test_a_smaller_qualifying_panel_is_not_reported_as_a_repair():
    """The row was asked for one size, and a shorter panel answers another row.

    Every panel of the requested size misses the floor here and every singleton
    clears it, so the beam has a qualifying answer to offer and it is the wrong
    shape for this row.
    """
    table = _table(
        {
            (GOOD, BAD): (0.90, 2.0, 40),
            (GOOD, SPARE): (0.85, 3.0, 45),
            (BAD, SPARE): (0.30, 1.0, 50),
        },
        singletons={GOOD: (0.30, 90.0, 2), BAD: (0.20, 90.0, 2), SPARE: (0.10, 90.0, 2)},
    )
    opt = _Optimizer([GOOD, BAD], table)

    row = plan_pool(opt, [GOOD, BAD, SPARE], [2], [0.8], min_selectivity_density=10)["rows"][0]

    assert row["size"] == 2
    assert set(row["primers"]) == {GOOD, BAD}
    assert row["repair"]["succeeded"] is False


def test_the_beam_repairs_what_the_swaps_cannot_reach():
    """A qualifying panel sharing no primer with the one the optimizer returned.

    One swap at a time cannot get there: every intermediate panel is worse than
    the incumbent under the same lexicographic rule, so the swap loop stops at a
    local optimum and the beam rebuilds instead.
    """
    table = _table(
        {
            (GOOD, BAD): (0.90, 2.0, 40),
            (GOOD, SPARE): (0.20, 1.0, 60),
            (BAD, SPARE): (0.25, 1.5, 55),
            (SPARE, HOOK): (0.85, 40.0, 3),
            (GOOD, HOOK): (0.30, 1.2, 50),
            (BAD, HOOK): (0.35, 1.1, 52),
        },
        singletons={
            GOOD: (0.10, 1.0, 20),
            BAD: (0.10, 1.0, 20),
            SPARE: (0.10, 1.0, 20),
            HOOK: (0.10, 1.0, 20),
        },
    )
    opt = _Optimizer([GOOD, BAD], table)

    row = plan_pool(opt, [GOOD, BAD, SPARE, HOOK], [2], [0.8], min_selectivity_density=10)["rows"][
        0
    ]

    assert set(row["primers"]) == {SPARE, HOOK}
    assert row["repair"]["method"] == "beam"
    assert row["eligible"] is True


def test_repair_does_not_change_the_panel_size():
    table = _table(
        {
            (GOOD, BAD): (0.90, 2.0, 40),
            (GOOD, SPARE): (0.88, 30.0, 5),
            (BAD, SPARE): (0.50, 25.0, 6),
        }
    )
    opt = _Optimizer([GOOD, BAD], table)

    row = plan_pool(opt, [GOOD, BAD, SPARE], [2], [0.8], min_selectivity_density=10)["rows"][0]

    assert row["size"] == 2


def test_repair_can_be_turned_off():
    """So the historical behaviour is still reachable for a comparison."""
    table = _table(
        {
            (GOOD, BAD): (0.90, 2.0, 40),
            (GOOD, SPARE): (0.88, 30.0, 5),
            (BAD, SPARE): (0.50, 25.0, 6),
        }
    )
    opt = _Optimizer([GOOD, BAD], table)

    row = plan_pool(
        opt,
        [GOOD, BAD, SPARE],
        [2],
        [0.8],
        min_selectivity_density=10,
        repair=False,
    )[
        "rows"
    ][0]

    assert set(row["primers"]) == {GOOD, BAD}
    assert row["eligible"] is False
    assert row["repair"]["attempted"] is False


def test_the_repaired_panel_is_re_evaluated_not_assumed():
    """Every returned pool goes back through the shared evaluator."""
    table = _table(
        {
            (GOOD, BAD): (0.90, 2.0, 40),
            (GOOD, SPARE): (0.88, 30.0, 5),
            (BAD, SPARE): (0.50, 25.0, 6),
        }
    )
    opt = _Optimizer([GOOD, BAD], table)

    row = plan_pool(opt, [GOOD, BAD, SPARE], [2], [0.8], min_selectivity_density=10)["rows"][0]

    assert (GOOD, SPARE) in [tuple(sorted(k)) for k in opt.evaluated]
    assert row["selectivity_density"] == 30.0
    assert row["background_sites"] == 5


def test_a_repair_never_introduces_a_dimerising_pair():
    """The swap-in guard is the dimer constraint, not a scoring term.

    The best-scoring repair here pairs HOOK with BARB, which the configured
    limit refuses. The second best qualifies and has lower coverage, so a
    repair that traded the constraint for coverage would be visible.
    """
    table = _table(
        {
            (GOOD, HOOK): (0.90, 2.0, 40),
            (BARB, HOOK): (0.99, 90.0, 1),
            (GOOD, BARB): (0.50, 25.0, 6),
        }
    )
    opt = _Optimizer([GOOD, HOOK], table)

    row = plan_pool(opt, [GOOD, HOOK, BARB], [2], [0.8], min_selectivity_density=10)["rows"][0]

    assert set(row["primers"]) != {BARB, HOOK}, "the dimerising pair was delivered"
    assert set(row["primers"]) == {GOOD, BARB}
    assert row["violating_pairs"] == 0


def test_the_report_names_the_rows_it_had_to_repair(tmp_path):
    """A repaired panel is not the one the optimizer first returned."""
    import csv

    from neoswga.core.pool_plan_report import write_pool_plan

    table = _table(
        {
            (GOOD, BAD): (0.90, 2.0, 40),
            (GOOD, SPARE): (0.88, 30.0, 5),
            (BAD, SPARE): (0.50, 25.0, 6),
        }
    )
    opt = _Optimizer([GOOD, BAD], table)
    plan = plan_pool(opt, [GOOD, BAD, SPARE], [2], [0.8], min_selectivity_density=10)

    write_pool_plan(plan, tmp_path / "out")

    rows = list(csv.DictReader((tmp_path / "out" / "pool_sizes.csv").open()))
    assert rows[0]["repaired_by"] == "swap"
    assert "repaired by swap" in (tmp_path / "out" / "pool_plan.html").read_text()


def test_repair_can_be_turned_off_from_the_command_line():
    """So the historical behaviour stays reachable without editing code."""
    import argparse

    from neoswga.cli.plan_pool import add_parsers

    parser = argparse.ArgumentParser()
    add_parsers(parser.add_subparsers(dest="command"))
    args = parser.parse_args(["plan-pool", "-j", "params.json", "--no-repair"])

    assert args.no_repair is True
