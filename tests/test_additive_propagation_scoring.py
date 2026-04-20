"""Tests that additive concentrations propagate into the quality scorer.

Before this fix, IntegratedQualityScorer used Wallace's rule for Tm, ignoring
DMSO/betaine/etc. A primer ranked optimal at 30 C pure buffer might be +20 C
too high at 42 C + 1 M betaine, but the optimizer would not see the difference.
"""

import pytest

from neoswga.core.integrated_quality_scorer import IntegratedQualityScorer
from neoswga.core.reaction_conditions import ReactionConditions


def test_betaine_lowers_thermo_score_when_primer_is_on_high_edge():
    """A primer whose raw Tm sits at the upper bound of the acceptance window
    scores LOWER once a Tm-lowering additive (betaine 1 M) is present, because
    the corrected Tm drops below the window."""
    primer = "GCGCGCGCGCGC"  # 12-mer, 100% GC, high raw Tm

    cond_no_additive = ReactionConditions(temp=30.0, polymerase='phi29')
    cond_betaine = ReactionConditions(temp=30.0, polymerase='phi29', betaine_m=1.0)

    scorer_plain = IntegratedQualityScorer(conditions=cond_no_additive)
    scorer_betaine = IntegratedQualityScorer(conditions=cond_betaine)

    tm_plain = scorer_plain._primer_tm(primer, 1.0, 12)
    tm_betaine = scorer_betaine._primer_tm(primer, 1.0, 12)

    assert tm_betaine < tm_plain, (
        f"Betaine should lower Tm; got plain={tm_plain:.1f}, betaine={tm_betaine:.1f}"
    )
    # The 1 M betaine correction should be ~-1.2 C or more (vendor data)
    assert tm_plain - tm_betaine >= 1.0, (
        f"Betaine 1 M should shift Tm by at least 1 C; got {tm_plain - tm_betaine:.2f}"
    )


def test_dmso_and_betaine_together_shift_more():
    primer = "GCATCGATCGAT"

    cond_plain = ReactionConditions(temp=30.0, polymerase='phi29')
    cond_cocktail = ReactionConditions(
        temp=30.0, polymerase='phi29',
        dmso_percent=5.0, betaine_m=1.0,
    )

    scorer_plain = IntegratedQualityScorer(conditions=cond_plain)
    scorer_cocktail = IntegratedQualityScorer(conditions=cond_cocktail)

    tm_plain = scorer_plain._primer_tm(primer, 0.5, 12)
    tm_cocktail = scorer_cocktail._primer_tm(primer, 0.5, 12)

    # 5% DMSO ~ -2.75 C, 1 M betaine ~ -1.2 C, combined > 3 C
    assert tm_plain - tm_cocktail >= 3.0, (
        f"DMSO+betaine cocktail should shift Tm by >=3 C; got {tm_plain - tm_cocktail:.2f}"
    )


def test_scorer_without_conditions_still_works():
    """Legacy path: no conditions provided. Should fall back and not crash."""
    scorer = IntegratedQualityScorer()  # no conditions
    tm = scorer._primer_tm("ACGTACGTACGT", 0.5, 12)
    # Wallace fallback: 2*6 + 4*6 = 36
    assert 20 <= tm <= 60


def test_additives_change_thermo_score_for_borderline_primer():
    """Scoring path change: a primer sitting at the score cliff moves when
    additives shift its effective Tm."""
    # Balanced 12-mer
    primer = "ATCGATCGATCG"

    # Set target_temp so optimal range is tight around the raw Tm
    cond_plain = ReactionConditions(temp=30.0, polymerase='phi29')
    # Strong Tm-depressing cocktail
    cond_depressed = ReactionConditions(
        temp=30.0, polymerase='phi29',
        betaine_m=2.0, dmso_percent=8.0,
    )

    scorer_plain = IntegratedQualityScorer(conditions=cond_plain)
    scorer_depressed = IntegratedQualityScorer(conditions=cond_depressed)

    score_plain = scorer_plain._calculate_thermo_score(primer)
    score_depressed = scorer_depressed._calculate_thermo_score(primer)

    # The cocktail lowers Tm significantly; the score must differ.
    assert score_plain != score_depressed, (
        f"Thermo score should change with heavy additives; got {score_plain} vs {score_depressed}"
    )
