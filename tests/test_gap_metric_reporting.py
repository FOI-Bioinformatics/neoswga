"""Guard: mean gap, max gap and Gini must all be surfaced, not collapsed.

Validation against two published primer-set datasets with wet-lab outcomes
(Clarke et al. 2017; Dwivedi-Yu et al. 2023) showed they disagree about which
gap statistic predicts success: density separated the winners in one, the worst
coverage hole and evenness in the other. The difference is that Clarke
pre-filtered their primer pool for even binding, leaving density as the limiting
property.

The conclusion was NOT to re-weight the composite score - six and twelve sets
are far too few to fit weights on, and a term tuned for one benchmark is wrong
for the other. It was to report all three with their benchmark context, so a
reader can see which regime they are in. These tests keep that reporting honest.

See docs/validation/published_primer_sets.md.
"""

import pytest

from neoswga.core.coverage import (
    BENCHMARK_GINI_GOOD,
    BENCHMARK_MEAN_GAP_BP,
    gap_regime_note,
    interpret_gap_metrics,
)


def test_all_three_statistics_are_reported():
    rows = interpret_gap_metrics(4100, 31500, 0.533)
    labels = [r["label"] for r in rows]
    assert any("Mean" in l for l in labels)
    assert any("Max" in l for l in labels)
    assert any("Gini" in l or "even" in l.lower() for l in labels)


def test_every_row_carries_a_verdict_and_a_reason():
    for row in interpret_gap_metrics(4100, 31500, 0.533):
        assert row["verdict"], row["label"]
        assert row["detail"], row["label"]
        # A verdict without a reference point is just an opinion.
        assert any(y in row["detail"] for y in ("2017", "2023", "reach"))


def test_mean_gap_is_judged_against_the_published_range():
    lo, hi = BENCHMARK_MEAN_GAP_BP
    inside = interpret_gap_metrics(hi - 500, 10_000, 0.5)[0]
    outside = interpret_gap_metrics(hi + 5_000, 10_000, 0.5)[0]
    assert "within" in inside["verdict"]
    assert "sparser" in outside["verdict"]


def test_max_gap_is_judged_against_the_reachable_span():
    """A gap up to twice the per-primer reach is still fully covered."""
    reach = 3000
    covered = interpret_gap_metrics(2000, 2 * reach - 1, 0.5, extension_reach=reach)[1]
    hole = interpret_gap_metrics(2000, 4 * reach, 0.5, extension_reach=reach)[1]
    assert covered["verdict"] == "fully reachable"
    assert "unreachable" in hole["verdict"]
    assert "6.0 kb" in hole["detail"], "should quantify the unreachable stretch"


def test_gini_is_judged_against_the_dwivedi_yu_winners():
    even = interpret_gap_metrics(3000, 10_000, BENCHMARK_GINI_GOOD - 0.02)[2]
    uneven = interpret_gap_metrics(3000, 10_000, BENCHMARK_GINI_GOOD + 0.05)[2]
    assert even["verdict"] == "even"
    assert uneven["verdict"] == "uneven"
    # The Clarke caveat must travel with the unfavourable verdict.
    assert "Clarke" in uneven["detail"]


def test_regime_note_points_at_the_limiting_property():
    uneven = gap_regime_note(BENCHMARK_GINI_GOOD + 0.1)
    even = gap_regime_note(BENCHMARK_GINI_GOOD - 0.1)
    assert "max gap" in uneven and "holes" in uneven
    assert "density" in even and "mean gap" in even


def test_zero_inputs_report_unknown_rather_than_good():
    """No measurement must not read as a clean bill of health."""
    for row in interpret_gap_metrics(0, 0, 0):
        assert row["verdict"] == "unknown"


# ----------------------------------------------------------------------
# Published-set sanity: the interpreter should agree with known outcomes
# ----------------------------------------------------------------------


@pytest.mark.parametrize(
    "name,mean_gap,max_gap,gini,expect_even",
    [
        ("Prev03 (96x, winner)", 4100, 31500, 0.533, True),
        ("Prev06 (120x, winner)", 4980, 31300, 0.552, True),
        ("Prev05 (2x, failure)", 2760, 37700, 0.638, False),
        ("Prev01 (2x, failure)", 7390, 65300, 0.600, False),
    ],
)
def test_interpreter_matches_known_prevotella_outcomes(name, mean_gap, max_gap, gini, expect_even):
    """The evenness verdict should track the measured outcome on this dataset."""
    rows = interpret_gap_metrics(mean_gap, max_gap, gini)
    is_even = rows[2]["verdict"] == "even"
    assert is_even == expect_even, name


def test_reporting_choice_is_documented_in_the_source():
    """The reason for not folding these into the score must stay visible."""
    import inspect

    from neoswga.core import coverage

    src = inspect.getsource(coverage)
    assert "Clarke" in src and "Dwivedi-Yu" in src
    assert "limiting" in src, (
        "the rationale - whichever property is limiting is the one that "
        "predicts - should be stated where the thresholds live"
    )


def test_gini_of_zero_with_gaps_present_means_even_not_unknown():
    """A single gap gives Gini 0 legitimately - that is even, not missing data.

    The first version reported "unknown / No gaps measured" whenever Gini was 0,
    which appeared in a real report alongside a measured 6.2 kb mean gap. Only
    the absence of any gap should read as unknown.
    """
    rows = interpret_gap_metrics(6200, 6200, 0.0)
    assert rows[0]["verdict"] != "unknown"
    assert rows[2]["verdict"] == "even"

    none_measured = interpret_gap_metrics(0, 0, 0.0)
    assert all(r["verdict"] == "unknown" for r in none_measured)
