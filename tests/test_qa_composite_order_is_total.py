"""`score --enable-qa` must not reintroduce the unstable sort.

Audit finding D1d. `pipeline_qa_integration.py:718` sorted on composite_score
alone with the default quicksort and no secondary key. Measured on the
449-candidate E. coli pool, two input orders gave 65 differing positions, the
first at rank 14 where two primers swap, and 100 of the 449 rows share a
composite value with another row.

That is the same defect `order_step3_rows` exists to remove, on a different
path. The optimizer is order-sensitive: on the E. coli pool at target size 24,
dominating-set returned a set with a Jaccard of 0.600 against the as-written
order when the candidates were reversed.
"""

import pandas as pd
import pytest


def _tied_pool(n=200):
    """A pool with heavy composite ties, which is the real shape: 100 of the
    449 E. coli rows share a composite value with another row."""
    return pd.DataFrame(
        {
            "primer": [f"{i:04d}ACGTAC" for i in range(n)],
            # Five distinct values across 200 rows: ties everywhere.
            "composite_score": [round(0.90 + 0.01 * (i % 5), 2) for i in range(n)],
        }
    )


def test_the_pool_really_is_mostly_ties():
    """Guard the guard: with distinct values the sort is total already."""
    df = _tied_pool()
    assert df["composite_score"].duplicated(keep=False).sum() > 0.9 * len(df)


def test_the_order_does_not_depend_on_how_the_rows_arrived():
    """The regression."""
    from neoswga.core.pipeline_qa_integration import order_by_composite_score

    df = _tied_pool()
    a = order_by_composite_score(df.sample(frac=1, random_state=0))["primer"].tolist()
    b = order_by_composite_score(df.sample(frac=1, random_state=1))["primer"].tolist()
    assert a == b


def test_the_composite_still_leads():
    from neoswga.core.pipeline_qa_integration import order_by_composite_score

    ordered = order_by_composite_score(_tied_pool())
    assert ordered["composite_score"].is_monotonic_decreasing


def test_ties_are_broken_by_the_primer_sequence():
    from neoswga.core.pipeline_qa_integration import order_by_composite_score

    ordered = order_by_composite_score(_tied_pool())
    for _score, group in ordered.groupby("composite_score", sort=False):
        primers = group["primer"].tolist()
        assert primers == sorted(primers)


def test_the_qa_naming_of_the_primer_column_also_works():
    """`primer_column` prefers 'primer' and falls back to 'seq'."""
    from neoswga.core.pipeline_qa_integration import order_by_composite_score

    df = _tied_pool().rename(columns={"primer": "seq"})
    ordered = order_by_composite_score(df)
    for _score, group in ordered.groupby("composite_score", sort=False):
        assert group["seq"].tolist() == sorted(group["seq"].tolist())


def test_an_empty_frame_is_not_an_error():
    from neoswga.core.pipeline_qa_integration import order_by_composite_score

    assert len(order_by_composite_score(_tied_pool().iloc[:0])) == 0
