"""Exact accounting for the specificity floor, and why Stage 1 does not use it.

The defect Phase 5 measured. Stage 1's greedy set cover chooses purely on
marginal coverage, so on a row with a selectivity floor it spends candidates on
coverage and the delivered panel misses the floor. Connecting Stage 2's
objective (`test_the_objective_reaches_the_stage_that_refines.py`) improved the
reported density from 28.78 to 42.62 on the measured pair but did not move the
feasible boundary, because one swap at a time cannot undo a panel built blind.

The exact target is known. `occupancy.weighted_site_load` is a SUM of per-primer
terms that depend only on the primer, so a panel's loads are additive and

    density(panel) >= D   <=>   SUM_i (f_i * bg_len - D * b_i * fg_len) >= 0

is linear in per-candidate quantities. On the measured pool the best 12-primer
panel reaches a density of 79.807, while the search delivered 60.11. That is
19.7 points of headroom, and it is a bound rather than a guess.

So Stage 1 needs no objective evaluation and no projection: the feasibility test
is an O(1) per-candidate sum. `SelectivityBudget` holds it.

The load-bearing property is the identity below, checked against the optimizer's
own `selectivity_density` rather than against the formula it came from. If the
two ever disagree, the budget is steering the search by a rule the design is not
accepted on, which is the two-rule split this whole area keeps producing.
"""

import pytest

import importlib.util
import pathlib as _pathlib

from neoswga.core.pool_objective import PoolConstraints

_spec = importlib.util.spec_from_file_location(
    "selectivity_budget",
    _pathlib.Path(__file__).resolve().parent.parent
    / "scripts"
    / "benchmarking"
    / "selectivity_budget.py",
)
_module = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_module)
SelectivityBudget = _module.SelectivityBudget


class _Loads:
    """Per-primer foreground and background loads, stated rather than computed."""

    def __init__(self, table):
        self._table = table

    def __call__(self, primer, which):
        return self._table[primer][0 if which == "fg" else 1]


# Per-primer loads chosen so the resulting densities span the floors used
# below. density(p) = (f / FG_LEN) / (b / BG_LEN), so these are 10000, 4000, 50
# and 12.5 respectively.
TABLE = {
    "AAAA": (100.0, 1.0),  # very selective
    "CCCC": (80.0, 2.0),
    "GGGG": (10.0, 20.0),  # host-heavy
    "TTTT": (5.0, 40.0),  # worse
}
FG_LEN, BG_LEN = 1_000_000.0, 100_000_000.0


def _budget(floor=None, ceiling=None, sites=None):
    return SelectivityBudget(
        constraints=PoolConstraints(min_selectivity_density=floor, max_background_sites=ceiling),
        fg_length=FG_LEN,
        bg_length=BG_LEN,
        loads={p: (f, b) for p, (f, b) in TABLE.items()},
        background_sites=sites or {p: int(b) for p, (_f, b) in TABLE.items()},
    )


def _density(panel):
    """The quantity the design is accepted on, from the same definition."""
    f = sum(TABLE[p][0] for p in panel)
    b = sum(TABLE[p][1] for p in panel)
    return (f / FG_LEN) / (b / BG_LEN)


# -- the identity the design rests on -------------------------------------


@pytest.mark.parametrize("floor", [1.0, 50.0, 100.0, 250.0, 400.0, 1000.0])
@pytest.mark.parametrize(
    "panel",
    [
        ["AAAA"],
        ["GGGG"],
        ["AAAA", "CCCC"],
        ["AAAA", "GGGG"],
        ["GGGG", "TTTT"],
        ["AAAA", "CCCC", "GGGG", "TTTT"],
    ],
)
def test_non_negative_slack_means_the_floor_is_met(floor, panel):
    """The whole point: an O(1) sum that agrees with the accepted metric.

    Checked against `_density`, which computes the ratio the way the optimizer
    does, rather than against the rearrangement the slack came from.
    """
    budget = _budget(floor=floor)
    slack = sum(budget.slack(p) for p in panel)

    assert (slack >= 0) is (
        _density(panel) >= floor
    ), f"slack {slack:.3e} and density {_density(panel):.3f} disagree at floor {floor}"


def test_no_constraint_means_no_budget():
    """Nothing to steer by, so Stage 1 must behave exactly as before."""
    assert (
        SelectivityBudget.build(
            constraints=PoolConstraints(),
            optimizer=None,
            candidates=list(TABLE),
        )
        is None
    )


def test_a_selective_candidate_has_positive_slack():
    budget = _budget(floor=1000.0)

    assert budget.slack("AAAA") > 0
    assert budget.slack("TTTT") < 0


def test_the_slack_scales_with_the_floor():
    """A stricter floor makes more candidates unaffordable.

    CCCC has a density of 4000, so it is affordable at 1000 and not at 8000.
    """
    lenient, strict = _budget(floor=1000.0), _budget(floor=8000.0)

    assert lenient.slack("CCCC") > strict.slack("CCCC")
    assert lenient.slack("CCCC") > 0 and strict.slack("CCCC") < 0


# -- the site ceiling ------------------------------------------------------


def test_the_site_ceiling_is_additive_too():
    """`max_background_sites` is a sum, so the same accounting works."""
    budget = _budget(ceiling=25)

    assert budget.admits(["AAAA", "CCCC"], "GGGG") is True
    assert budget.admits(["GGGG"], "TTTT") is False


def test_the_site_ceiling_sum_is_conservative():
    """It bounds a union, so using the sum can only be stricter, never looser.

    `total_bg_sites` is the size of a union over primers, and two primers can
    share a background position. The sum is therefore an upper bound, so a
    panel this admits certainly satisfies the real ceiling.
    """
    budget = _budget(ceiling=60)

    assert budget.admits(["GGGG"], "TTTT") is True  # 20 + 40 == 60


# -- and it reaches Stage 1 ------------------------------------------------
