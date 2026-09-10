"""Step 3 must not discard step 2's ranking.

Audit finding D1c. `_rank_by_occupancy` (pipeline.py:837) spends up to 50 s
ranking candidates by occupancy-weighted background load per unit of foreground
binding, and its docstring credits that ranking with moving coverage from 7.6%
to 40.3% and selectivity from 0.69 to 2.54. `step2` writes step2_df.csv in that
order. `order_step3_rows` then re-sorted by Gini and threw it away.

The two orders are Spearman -0.185 on the E. coli pool and share 2 of their top
24. The head of the pool the optimizer received was measurably the worse end:

    top 24 by            mean fg sites   mean bg sites   mean occupancy ratio
    step 2's own rank    35.7            5.0             3.74
    Gini, as shipped     10.5            29.7            40.34

HONEST CAVEAT, and the reason every assertion here is about order rather than
about coverage. No effect on the delivered set could be detected. Four
orderings of the E. coli pool -- Gini, step 2's rank, Gini reversed, and a
seeded shuffle -- gave Jaccard 1.000 and identical coverage under both
dominating-set and hybrid at target sizes 6, 12 and 24. This is worth fixing
because the inversion is indefensible and the fix is free, not on a promise of
changed output.

The rank is read off the file's row order rather than recomputed, because step 2
uses occupancy_ratio when it can and falls back to ratio then fg_count when it
cannot. Row order captures whichever key actually ran.
"""

import pandas as pd
import pytest

# gini deliberately DISAGREES with the file order, so a test that passes here
# cannot be passing by coincidence.
_STEP2_ROWS = [
    # primer, ratio, gini, fg_count, bg_count, occupancy_ratio
    ("TTGGCCAATGCA", 4.0, 0.55, 12, 3, 0.8),
    ("CACCGACGACGA", 6.0, 0.45, 18, 3, 1.4),
    ("AAGGCCTTACGT", 7.0, 0.30, 21, 3, 2.9),
    ("ACCCGGGTTTAC", 9.5, 0.30, 33, 4, 5.1),
    ("AAACCCGGGTTT", 12.0, 0.20, 40, 3, 9.7),
]
_COLUMNS = ["primer", "ratio", "gini", "fg_count", "bg_count", "occupancy_ratio"]


@pytest.fixture
def primed(tmp_path, monkeypatch):
    """A data_dir holding only step2_df.csv, with the pipeline pointed at it."""
    import neoswga.core.pipeline as pipeline_mod
    from neoswga.core import parameter

    pd.DataFrame(_STEP2_ROWS, columns=_COLUMNS).to_csv(tmp_path / "step2_df.csv", index=False)

    monkeypatch.setattr(pipeline_mod, "_initialize", lambda: None)
    monkeypatch.setattr(parameter, "data_dir", str(tmp_path), raising=False)
    monkeypatch.setattr(parameter, "verbose", False, raising=False)
    monkeypatch.setattr(parameter, "min_amp_pred", None, raising=False)
    monkeypatch.setattr(parameter, "use_amp_model", False, raising=False)
    return tmp_path


def _run(primed):
    import neoswga.core.pipeline as pipeline_mod

    pipeline_mod.step3(validate_prerequisites=False)
    return pd.read_csv(primed / "step3_df.csv")


def test_the_fixture_really_does_invert_the_two_orders():
    """Guard the guard. If gini agreed with the file order this file would
    prove nothing."""
    df = pd.DataFrame(_STEP2_ROWS, columns=_COLUMNS)
    assert not df["gini"].is_monotonic_increasing


def test_step3_preserves_the_order_step2_wrote(primed):
    """The regression: step 3 re-sorted by gini and discarded step 2's rank."""
    out = _run(primed)
    assert out["primer"].tolist() == [row[0] for row in _STEP2_ROWS]


def test_the_rank_is_carried_as_a_column(primed):
    """So a reader of step3_df.csv can see the order is not an accident."""
    out = _run(primed)
    assert "step2_rank" in out.columns
    assert out["step2_rank"].tolist() == [0, 1, 2, 3, 4]


def test_gini_is_demoted_to_a_tie_break():
    """Rank leads; gini decides only where rank does not."""
    from neoswga.core.pipeline import order_step3_rows

    df = pd.DataFrame(
        {
            "primer": ["AAAA", "CCCC", "GGGG"],
            "step2_rank": [1, 1, 0],
            "gini": [0.9, 0.1, 0.5],
        }
    )
    assert order_step3_rows(df)["primer"].tolist() == ["GGGG", "CCCC", "AAAA"]


def test_the_primer_sequence_is_last_for_totality():
    """Rank and gini can both tie; the order must still be total."""
    from neoswga.core.pipeline import order_step3_rows

    df = pd.DataFrame(
        {
            "primer": ["CCCC", "AAAA"],
            "step2_rank": [3, 3],
            "gini": [0.4, 0.4],
        }
    )
    assert order_step3_rows(df)["primer"].tolist() == ["AAAA", "CCCC"]


def test_a_frame_without_a_rank_still_orders_by_gini_then_primer():
    """Library callers and the amp-model path can hand over a frame that never
    passed through step 2's writer."""
    from neoswga.core.pipeline import order_step3_rows

    df = pd.DataFrame(
        {
            "primer": ["CCCC", "AAAA", "GGGG"],
            "gini": [0.4, 0.4, 0.1],
        }
    )
    assert order_step3_rows(df)["primer"].tolist() == ["GGGG", "AAAA", "CCCC"]


def test_the_order_does_not_depend_on_how_the_rows_arrived(primed):
    """The property the ordering work established, which must survive this."""
    from neoswga.core.pipeline import order_step3_rows

    df = pd.DataFrame(_STEP2_ROWS, columns=_COLUMNS)
    df["step2_rank"] = range(len(df))
    a = order_step3_rows(df.sample(frac=1, random_state=0))["primer"].tolist()
    b = order_step3_rows(df.sample(frac=1, random_state=1))["primer"].tolist()
    assert a == b


def test_running_twice_gives_the_same_file(primed):
    assert _run(primed)["primer"].tolist() == _run(primed)["primer"].tolist()


def test_an_empty_frame_is_not_an_error():
    from neoswga.core.pipeline import order_step3_rows

    assert len(order_step3_rows(pd.DataFrame({"primer": [], "step2_rank": []}))) == 0
