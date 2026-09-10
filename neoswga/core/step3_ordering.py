"""What `step3_df.csv` carries, and in what order.

Split out of ``core/pipeline.py``, which was over its module size budget.
The two functions here decide the columns the candidate pool carries into
step 4 and the row order it arrives in, and that order reaches an
order-sensitive optimizer. `pipeline` imports both under the same names.
"""

import logging

logger = logging.getLogger(__name__)


def order_step3_rows(df):
    """A deterministic total order for the scored candidate pool.

    `step3_df.csv` was written with `sort_values(by="gini")` alone. On a real
    pool almost every row ties: 496 of 500 on the plasmid example share a gini
    value. `sort_values` defaults to quicksort, which is not stable, so for
    those rows the order was whatever the algorithm produced from the order the
    rows happened to arrive in -- the same data from two different input orders
    gave two different files, sharing 7 of the first 50 primers.

    That would not matter if the optimizer ignored order. It does not. On the
    E. coli pool at target size 24, dominating-set returned a set with a Jaccard
    of 0.600 against the as-written order when the candidates were reversed, and
    0.920 when they were shuffled. Up to 40% of the delivered oligos were
    decided by a tie-break nobody chose.

    WHAT LEADS THE ORDER, changed by audit finding D1c on 2026-09-10. Step 2
    spends up to 50 s ranking candidates by occupancy-weighted background load
    (`_rank_by_occupancy`), a ranking its own docstring credits with moving
    coverage from 7.6% to 40.3% and selectivity from 0.69 to 2.54, and writes
    step2_df.csv in that order. Sorting by gini here discarded it. The two
    orders are Spearman -0.185 on the E. coli pool and share 2 of their top 24,
    and the head this handed the optimizer was the worse end: mean occupancy
    ratio 40.34 against 3.74. `dominating_set_optimizer` documents its tie-break
    as respecting "the caller's ranking", which had not arrived since the
    ordering work landed.

    NO CHANGE IN OUTPUT IS CLAIMED. Four orderings of the E. coli pool -- gini,
    step 2's rank, gini reversed and a seeded shuffle -- gave Jaccard 1.000 and
    identical coverage under both dominating-set and hybrid at target sizes 6,
    12 and 24. This removes an inversion that is indefensible on its face and
    costs nothing to remove; it is not a promise of a better panel.

    Gini is demoted to a tie-break and the primer sequence is last, which is
    unique and makes the order total. A frame carrying no `step2_rank` -- a
    library caller, or one that never passed through step 2's writer -- keeps
    the previous gini-then-primer order.

    WHY NOT ORDER BY `amp_pred`. It is the obvious alternative and the evidence
    is against it: selecting the top half of a pool by `amp_pred` and optimizing
    over it produced the WORST of five half-pools, behind all three random
    halves and behind the bottom half by the same measure.
    """
    if len(df) == 0:
        return df
    # `_candidate_carry_columns` sets `primer` as the index, so the sort key
    # has to be looked for on the index as well as on the columns. The previous
    # `sort_values(by=["gini", "primer"])` relied on pandas resolving `primer`
    # against the index; naming the keys explicitly must not lose that.
    available = set(df.columns)
    if df.index.name is not None:
        available.add(df.index.name)
    keys = [column for column in ("step2_rank", "gini", "primer") if column in available]
    if not keys:
        return df
    return df.sort_values(by=keys, kind="mergesort")


def _candidate_carry_columns(step2_df):
    """The step-2 measurements `step3_df.csv` carries forward.

    `step2_rank` is load-bearing: it leads the row order `order_step3_rows`
    establishes, and that order reaches the optimizer. `gini` follows it as a
    tie-break and is also read by the report.
    """
    df = step2_df.set_index("primer")
    keep = [c for c in ("step2_rank", "ratio", "gini", "fg_count", "bg_count") if c in df.columns]
    return df[keep]
