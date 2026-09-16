"""The beam has to run on the pools this tool actually has.

Phase 2 of the plan for `docs/validation/pipeline_audit_2026-09-16/`, the second
half of audit finding F2.

The beam was gated on `4 * len(pool) * size <= swap_max_evaluations - swaps_used`.
With the default budget of 10,000 that is `len(pool) * size <= 2500`:

    pool      largest panel the beam could run for
      30                                        83
     200                                        12
     500                                         5
    2000                                         1

The real shortlist is 2,000 candidates and panels are typically 8 to 24, so the
beam never ran. Its integration test used 30 candidates, where the bound allows
a panel of 83, which is why fourteen unit tests and a seven-test integration
suite all passed on a path production could not take.

The fix is not a larger budget. It is to stop asking whether the whole pool fits
and start searching the best-ranked slice the budget affords. The pool arrives in
step-2 rank order, so a prefix of it is the best candidates rather than an
arbitrary subset, and the incumbent panel is always included so the beam can
still return it.
"""

from types import SimpleNamespace

import pytest

from neoswga.core.base_optimizer import OptimizationStatus, OptimizerConfig
from neoswga.core.pool_planner import _BEAM_WIDTH, plan_pool


def _oligos(count):
    """Distinct 12-mers in base-4, so none collides and none self-dimerises."""
    out = []
    for index in range(count):
        digits, value = [], index
        for _ in range(12):
            digits.append("ACGT"[value % 4])
            value //= 4
        out.append("".join(reversed(digits)))
    return out


class _Optimizer:
    """Returns a deliberately poor panel, so the row needs repairing.

    Metrics are computed rather than tabulated: at these pool sizes a table
    would need more entries than the search touches, and computing them keeps
    the fixture honest about panels the search invents.
    """

    name = "test"
    bg_prefixes = ["bg"]
    bg_seq_lengths = [1000]
    conditions = object()

    def __init__(self, worst, evaluations=None):
        # The worst-ranked `worst` candidates of whatever pool it is handed. The
        # planner removes self-dimering candidates before optimizing, so a panel
        # picked from the raw list can fall outside the bounds it then checks.
        self.worst = int(worst)
        self.config = OptimizerConfig(max_dimer_bp=3, max_self_dimer_bp=4)
        self.evaluations = evaluations if evaluations is not None else []
        self.delivered = []

    def optimize(self, candidates, target_size):  # noqa: ARG002
        self.delivered = list(candidates)[-self.worst :]
        return SimpleNamespace(
            primers=self.delivered, status=OptimizationStatus.PARTIAL, message=""
        )

    def compute_metrics(self, primers):
        self.evaluations.append(tuple(sorted(primers)))
        # Coverage rises with how early the primers rank, so the best panel is
        # the top of the pool and a poor incumbent has somewhere to go.
        ranks = [sum(ord(c) for c in p) for p in primers]
        coverage = min(0.99, 0.30 + 0.0004 * len(primers) * (max(ranks) - min(ranks) + 1) / 100)
        return SimpleNamespace(
            effective_fg_coverage=coverage,
            fg_coverage=coverage,
            selectivity_density=90.0,
            total_bg_sites=len(primers),
            max_gap=100,
        )


@pytest.mark.parametrize("pool_size,panel", [(500, 8), (2000, 12), (2000, 24)])
def test_the_beam_runs_at_production_pool_sizes(pool_size, panel):
    """The defect, stated as the sizes that failed.

    Every one of these combinations reported "not affordable" before.
    """
    pool = _oligos(pool_size)
    opt = _Optimizer(panel)

    row = plan_pool(opt, pool, [panel], [0.99], min_selectivity_density=10)["rows"][0]

    assert row["repair"]["attempted"] is True
    assert row["repair"].get("beam") not in (
        None,
        "not affordable within the remaining budget",
    ), f"the beam still does not run at {pool_size} candidates and {panel} primers"


def test_the_beam_stays_inside_its_budget():
    """Reaching the beam must not be paid for with an unbounded search."""
    pool = _oligos(2000)
    seen = []
    opt = _Optimizer(12, evaluations=seen)

    row = plan_pool(opt, pool, [12], [0.99], min_selectivity_density=10)["rows"][0]

    ceiling = opt.config.swap_max_evaluations + opt.config.beam_max_evaluations
    assert row["repair"]["evaluations"] <= ceiling
    assert (
        len(seen) <= ceiling * 2
    ), "the evaluator was called far more often than the two budgets allow"


def test_the_row_reports_how_much_of_the_pool_the_beam_saw():
    """A panel chosen from a slice must not read as one chosen from the pool."""
    pool = _oligos(2000)
    opt = _Optimizer(12)

    row = plan_pool(opt, pool, [12], [0.99], min_selectivity_density=10)["rows"][0]

    searched = row["repair"]["beam_candidates"]
    assert 0 < searched <= len(pool)
    assert searched < len(pool), "at this size the beam cannot afford the whole pool"


def test_the_incumbent_is_always_inside_the_slice():
    """Otherwise the beam could not return the panel it started from.

    The incumbent here is the LAST twelve oligos in rank order, so a plain
    prefix of the pool excludes every one of them. If the slice were a bare
    prefix the beam could not reproduce the incumbent, and a repair that failed
    to improve would look like one that found something worse.
    """
    pool = _oligos(2000)
    seen = []
    opt = _Optimizer(12, evaluations=seen)

    plan_pool(opt, pool, [12], [0.99], min_selectivity_density=10)

    searched = {primer for panel in seen for primer in panel}
    missing = [p for p in opt.delivered if p not in searched]
    assert not missing, f"{len(missing)} incumbent primers were never evaluated"


def test_a_small_pool_still_searches_all_of_it():
    """The slice is a budget, not a cap. Nothing is excluded when it all fits."""
    pool = _oligos(24)
    opt = _Optimizer(4)

    row = plan_pool(opt, pool, [4], [0.99], min_selectivity_density=10)["rows"][0]

    assert row["repair"]["beam_candidates"] == len(pool)


def test_the_slice_scales_with_the_budget():
    """A bigger budget buys a wider search, which is the knob a user has."""
    pool = _oligos(2000)

    narrow = _Optimizer(12)
    narrow.config = OptimizerConfig(max_dimer_bp=3, max_self_dimer_bp=4, beam_max_evaluations=2000)
    wide = _Optimizer(12)
    wide.config = OptimizerConfig(max_dimer_bp=3, max_self_dimer_bp=4, beam_max_evaluations=60000)

    narrow_row = plan_pool(narrow, pool, [12], [0.99], min_selectivity_density=10)["rows"][0]
    wide_row = plan_pool(wide, pool, [12], [0.99], min_selectivity_density=10)["rows"][0]

    assert wide_row["repair"]["beam_candidates"] > narrow_row["repair"]["beam_candidates"]


def test_the_beam_width_is_still_what_the_bound_assumes():
    """Guard the guard: the slice arithmetic is derived from the width."""
    assert _BEAM_WIDTH == 4
