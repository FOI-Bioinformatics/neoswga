"""The candidate order handed to the optimizer must not be an accident.

Found while investigating F0 -- "the `score` step does not change the delivered
set". The conclusion about `amp_pred` is that it should NOT steer selection
(see below), but the investigation turned up something larger about the file it
writes.

`step3_df.csv` is written sorted by `gini`. On the plasmid pool 496 of 500 rows
share a gini value, so for almost every row the sort decides nothing and the
order is settled by tie-breaking. `DataFrame.sort_values` defaults to quicksort,
which is not stable, so the tie order is whatever the algorithm happens to
produce from the order the rows arrived in:

    same data, two different input orders, default sort
      identical output order? False
      overlap of the first 50: 7/50

That would be harmless if the optimizer ignored order. It does not. On the
E. coli pool at target size 24, dominating-set returns:

    candidate order          Jaccard against the as-written set
      as written (gini)        1.000
      amp_pred descending      1.000
      amp_pred ascending       0.600
      seeded shuffle           0.920

So up to 40% of the delivered oligos are decided by a tie-break nobody chose.

The fix is a TOTAL order: gini, then the primer sequence, which is unique. That
is deterministic, independent of how the rows arrived, and claims nothing about
primer quality.

WHY NOT ORDER BY `amp_pred`. It is the obvious candidate and the evidence is
against it. Under the default `fast_score=True` the delta-G features are zeroed,
so `amp_pred` is a pure function of the primer sequence -- blind to the genome
and to the reaction -- and 83-99% predictable from base composition alone. Asked
whether it identifies good primers, selecting the top half by `amp_pred` and
optimizing over it produced the WORST of five half-pools:

    1. random half 1           score=0.6580
    2. random half 2           score=0.6562
    3. random half 3           score=0.6538
    4. bottom half by amp_pred score=0.6530
    5. top half by amp_pred    score=0.6502

Using it as the tie-break would systematically prefer the primers that measured
worst. A neutral, stated tie-break is the honest choice.
"""

import pandas as pd
import pytest


def _pool(n=200):
    """A pool with heavy gini ties, which is the real shape."""
    rows = []
    for i in range(n):
        rows.append(
            {
                "primer": f"{i:04d}ACGTAC",
                # Four distinct values across 200 rows: ties everywhere.
                "gini": round(0.30 + 0.05 * (i % 4), 2),
                "on.target.pred": 10.0 + (i % 37) / 3.0,
            }
        )
    return pd.DataFrame(rows)


def test_the_pool_really_is_mostly_ties():
    """Guard the guard: with distinct gini values the sort is a total order
    already and this file would prove nothing."""
    df = _pool()
    assert df["gini"].duplicated(keep=False).sum() > 0.9 * len(df)


def test_the_order_does_not_depend_on_how_the_rows_arrived():
    """The property that was missing.

    Sorting the same data from two different starting orders must give the same
    result. With `sort_values(by="gini")` alone it does not.
    """
    from neoswga.core.pipeline import order_step3_rows

    df = _pool()
    a = order_step3_rows(df.sample(frac=1, random_state=0))["primer"].tolist()
    b = order_step3_rows(df.sample(frac=1, random_state=1))["primer"].tolist()

    assert a == b, (
        "the same rows in a different input order produced a different output "
        "order; the tie-break is an accident of the sort algorithm, and it "
        "changes which oligos the optimizer selects"
    )


def test_gini_still_leads_the_order():
    """The tie-break must not displace the actual sort key."""
    from neoswga.core.pipeline import order_step3_rows

    ordered = order_step3_rows(_pool())
    assert ordered["gini"].is_monotonic_increasing


def test_ties_are_broken_by_the_primer_sequence():
    """A stated, neutral rule rather than an emergent one."""
    from neoswga.core.pipeline import order_step3_rows

    ordered = order_step3_rows(_pool())
    for _gini, group in ordered.groupby("gini", sort=False):
        primers = group["primer"].tolist()
        assert primers == sorted(primers), f"ties within gini={_gini} are not ordered by primer"


def test_the_tie_break_is_not_amp_pred():
    """Ordering by `amp_pred` would prefer the primers that measured worst.

    Pinned so a future change does not reach for the obvious signal without
    re-running the half-pool comparison in this file's docstring.
    """
    from neoswga.core.pipeline import order_step3_rows

    ordered = order_step3_rows(_pool())
    first_group = ordered[ordered["gini"] == ordered["gini"].iloc[0]]
    assert not first_group["on.target.pred"].is_monotonic_decreasing, (
        "ties appear to be ordered by amp_pred descending; the evidence says "
        "that selects the worst-performing half of a pool"
    )


def test_an_empty_frame_is_not_an_error():
    from neoswga.core.pipeline import order_step3_rows

    assert len(order_step3_rows(_pool().iloc[:0])) == 0
