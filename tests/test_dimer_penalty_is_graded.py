"""The dimer penalty must grade with duplex length, not saturate.

`calculate_dimer_score` used to be a step function: 1.0 above `max_dimer_bp`
and 0.0 below. At the settings used for the GC-tiered EquiPhi29 panels
(12-mers, max_dimer_bp=3) that put 328 of the 630 pairs at the ceiling, so a
pair sharing four complementary bases and a pair sharing eleven scored the
same. Selection had no gradient in the region that matters, and the network
optimizer put a near reverse-complement pair into the M. tuberculosis panel.

The pair below is the real one, taken from
runs/gc_tiers/high_mtb/step4_improved_df_summary.json.
"""

import pytest

from neoswga.core.dimer import max_complementary_run
from neoswga.core.network_optimizer import calculate_dimer_score

# Both pairs are from the shipped 36-primer M. tuberculosis panel.
WORST_PAIR = ("TCGTCAGACCCA", "GGGTCTGACGAC")  # 11 bp duplex
MILD_PAIR = ("CGGCACCGGCAA", "GCCGCCATCGAC")  # 4 bp duplex

PANEL_MAX_DIMER_BP = 3


def test_longer_duplex_is_penalised_more_than_shorter_one():
    """The regression test for the saturation defect.

    Before the fix both of these returned 1.0.
    """
    worst = calculate_dimer_score(*WORST_PAIR, PANEL_MAX_DIMER_BP)
    mild = calculate_dimer_score(*MILD_PAIR, PANEL_MAX_DIMER_BP)

    assert worst > mild, (
        f"11 bp duplex scored {worst}, 4 bp duplex scored {mild}: "
        "the penalty carries no gradient between them"
    )


def test_score_rises_monotonically_with_duplex_length():
    """A whole ladder of duplex lengths, not just two points."""
    # Each primer pairs with the first over a progressively longer stretch.
    anchor = "AAAACCCCGGGG"
    partners = [
        "CCCCGGGGTTTT",  # 12 bp, full-length duplex
        "TCCCGGGGTTTT",  # 11 bp
        "TTCCGGGGTTTT",  # 10 bp
    ]
    runs = [max_complementary_run(anchor, p) for p in partners]
    assert runs == [12, 11, 10], f"test fixture is wrong: got runs {runs}"

    scores = [calculate_dimer_score(anchor, p, PANEL_MAX_DIMER_BP) for p in partners]
    assert scores == sorted(scores, reverse=True)
    assert len(set(scores)) == 3, f"expected three distinct scores, got {scores}"


def test_no_penalty_at_or_below_the_threshold():
    """Preserved behaviour: at or under max_dimer_bp there is no risk."""
    pair = ("AAAAAAAAAAAA", "TTTAAAAAAAAA")
    assert max_complementary_run(*pair) == 3

    # Exactly at the threshold is still no risk, matching is_dimer_fast.
    assert calculate_dimer_score(*pair, 3) == 0.0
    # And comfortably under it.
    assert calculate_dimer_score(*pair, 4) == 0.0


def test_full_length_duplex_still_scores_one():
    """Preserved endpoint: a complete duplex is maximum risk."""
    score = calculate_dimer_score("AAAACCCCGGGG", "CCCCGGGGTTTT", PANEL_MAX_DIMER_BP)
    assert score == pytest.approx(1.0)


def test_score_stays_within_unit_interval():
    for a, b in (WORST_PAIR, MILD_PAIR, ("AAAACCCCGGGG", "CCCCGGGGTTTT")):
        for threshold in (2, 3, 4, 5):
            score = calculate_dimer_score(a, b, threshold)
            assert 0.0 <= score <= 1.0


class TestMaxComplementaryRun:
    def test_measures_the_longest_complementary_stretch(self):
        assert max_complementary_run(*WORST_PAIR) == 11
        assert max_complementary_run(*MILD_PAIR) == 4

    def test_a_primer_against_its_own_reverse_complement_runs_full_length(self):
        assert max_complementary_run("AAAACCCCGGGG", "CCCCGGGGTTTT") == 12

    def test_never_exceeds_the_shorter_primer(self):
        assert max_complementary_run("AAAACCCCGGGG", "CCCC") == 4

    def test_agrees_with_the_boolean_check(self):
        from neoswga.core.dimer import is_dimer_fast

        for a, b in (WORST_PAIR, MILD_PAIR):
            for threshold in (2, 3, 4, 5, 10, 11):
                assert is_dimer_fast(a, b, threshold) == (
                    max_complementary_run(a, b) > threshold
                ), f"disagreement on {a}/{b} at threshold {threshold}"
