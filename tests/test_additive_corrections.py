"""
Consolidated additive correction tests for NeoSWGA.

Covers DMSO, betaine, formamide, trehalose, urea, TMAC, and ethanol
effects on Tm, GC-dependent corrections, Arrhenius temperature scaling,
combined additive effects, and validation against literature values.

Literature references:
    - Chester & Marshak (1993) Anal Biochem 209:284-290 (DMSO)
    - Varadaraj & Skinner (1994) Gene 140:1-5 (DMSO)
    - Rees et al. (1993) Biochemistry 32:137-144 (betaine)
    - Henke et al. (1997) NAR 25:3957-3958 (betaine)
    - Blake & Delcourt (1996) NAR 24:2095-2103 (formamide)
    - Spiess et al. (2004) Biotechniques 36:732-736 (trehalose)
    - Lesnick & Bhalla (1995) NAR 23:4665-4666 (urea)
    - Melchior & von Hippel (1973) PNAS 70:298-302 (TMAC)
    - Cheng et al. (1994) PNAS 91:5695-5699 (ethanol)
"""

import pytest

from neoswga.core.additives import (
    AdditiveConcentrations,
    ArrheniusTmCorrector,
)
from neoswga.core.mechanistic_params import (
    ADDITIVE_TM_PARAMS,
    get_additive_tm_params,
    list_additives,
)
from neoswga.core.reaction_conditions import (
    ReactionConditions,
    get_extreme_gc_conditions,
    get_high_gc_conditions,
    get_standard_conditions,
)
from neoswga.core.thermodynamics import calculate_tm_with_salt, gc_content


# =============================================================================
# DMSO Corrections
# =============================================================================


class TestDMSOCorrections:
    """Tests for DMSO Tm corrections."""

    def test_dmso_reduces_tm(self):
        """DMSO should reduce Tm by approximately 0.6C per percent."""
        seq = "ATCGATCGATCG"
        cond_0 = ReactionConditions(temp=30, dmso_percent=0)
        cond_5 = ReactionConditions(temp=30, dmso_percent=5)
        cond_10 = ReactionConditions(temp=30, dmso_percent=10)

        tm_0 = cond_0.calculate_effective_tm(seq)
        tm_5 = cond_5.calculate_effective_tm(seq)
        tm_10 = cond_10.calculate_effective_tm(seq)

        # 5% DMSO ~ -3C
        assert 2.5 < (tm_0 - tm_5) < 3.5
        # 10% DMSO ~ -6C
        assert 5.0 < (tm_0 - tm_10) < 6.5

    def test_dmso_effect_approximately_linear(self):
        """DMSO effect should be approximately linear with concentration."""
        base = ReactionConditions(polymerase="equiphi29", temp=42.0)
        dmso5 = ReactionConditions(polymerase="equiphi29", temp=42.0, dmso_percent=5.0)
        dmso10 = ReactionConditions(
            polymerase="equiphi29", temp=42.0, dmso_percent=10.0
        )

        c0 = base.calculate_tm_correction(0.5, 12)
        c5 = dmso5.calculate_tm_correction(0.5, 12)
        c10 = dmso10.calculate_tm_correction(0.5, 12)

        delta_5 = c5 - c0
        delta_10 = c10 - c0
        assert abs(delta_10 - 2 * delta_5) < 1.0

    def test_dmso_at_10_percent(self):
        """10% DMSO should reduce Tm by about 6C."""
        conditions = ReactionConditions(
            polymerase="equiphi29", temp=42.0, dmso_percent=10.0
        )
        correction = conditions.calculate_tm_correction(gc_content=0.5, primer_length=12)
        assert correction < -5

    def test_dmso_arrhenius_at_reference_temp(self):
        """DMSO Arrhenius correction at 37C should match ref_coef * concentration."""
        corrector = ArrheniusTmCorrector(37.0)
        correction = corrector.calculate_correction("dmso", 5.0)
        expected = -0.55 * 5.0
        assert abs(correction - expected) < 0.1

    def test_dmso_temperature_dependence(self):
        """DMSO correction should increase with temperature."""
        corr_30 = ArrheniusTmCorrector(30.0).calculate_correction("dmso", 5.0)
        corr_42 = ArrheniusTmCorrector(42.0).calculate_correction("dmso", 5.0)
        assert corr_42 < corr_30
        assert corr_30 < 0
        assert corr_42 < 0

    def test_dmso_literature_range(self):
        """
        DMSO correction should match literature.

        Chester & Marshak (1993): -0.5 to -0.6C per % at 37C.
        Varadaraj (1994): 10% DMSO lowers Tm by approximately 6C.
        """
        corrector = ArrheniusTmCorrector(37.0)
        correction = corrector.calculate_correction("dmso", 10.0)
        assert -7.0 < correction < -4.0


# =============================================================================
# Betaine Corrections
# =============================================================================


class TestBetaineCorrections:
    """Tests for betaine Tm corrections."""

    def test_betaine_reduces_tm_at_balanced_gc(self):
        """
        Betaine reduces Tm of 50% GC sequences by the uniform component
        plus residual GC-dependent correction.
        """
        seq = "ATCGATCGATCG"  # 50% GC
        cond_0 = ReactionConditions(temp=30, betaine_m=0)
        cond_1 = ReactionConditions(temp=30, betaine_m=1.0)
        cond_2 = ReactionConditions(temp=30, betaine_m=2.0)

        tm_0 = cond_0.calculate_effective_tm(seq)
        tm_1 = cond_1.calculate_effective_tm(seq)
        tm_2 = cond_2.calculate_effective_tm(seq)

        assert 0.3 < (tm_0 - tm_1) < 1.5
        assert 0.7 < (tm_0 - tm_2) < 3.0

    def test_betaine_gc_dependent_effect(self):
        """
        Betaine reduces Tm more for GC-rich sequences and may increase Tm
        for AT-rich sequences (GC normalization).
        """
        seq_gc = "GCGCGCGC"  # 100% GC
        seq_at = "AAAATTTT"  # 0% GC

        cond_0 = ReactionConditions(temp=30, betaine_m=0)
        cond_1 = ReactionConditions(temp=30, betaine_m=1.0)

        gc_reduction = cond_0.calculate_effective_tm(seq_gc) - cond_1.calculate_effective_tm(seq_gc)
        at_effect = cond_0.calculate_effective_tm(seq_at) - cond_1.calculate_effective_tm(seq_at)

        assert gc_reduction > 1.0
        assert at_effect < gc_reduction

    def test_betaine_at_2_5m(self):
        """Betaine near isostabilizing concentration (2.5M)."""
        conditions = ReactionConditions(
            polymerase="equiphi29", temp=42.0, betaine_m=2.5
        )
        corr_50 = conditions.calculate_tm_correction(gc_content=0.5, primer_length=12)
        corr_33 = conditions.calculate_tm_correction(gc_content=0.33, primer_length=12)
        assert isinstance(corr_50, float)
        assert isinstance(corr_33, float)

    def test_betaine_reduces_tm_range(self):
        """At high betaine, GC-rich and AT-rich sequences should converge in Tm."""
        seq_gc = "GCGCGCGCGC"
        seq_at = "AAAAAATTTT"

        cond_0 = ReactionConditions(temp=30, betaine_m=0)
        cond_2 = ReactionConditions(temp=30, betaine_m=2.0)

        range_no = (
            cond_0.calculate_effective_tm(seq_gc)
            - cond_0.calculate_effective_tm(seq_at)
        )
        range_2m = (
            cond_2.calculate_effective_tm(seq_gc)
            - cond_2.calculate_effective_tm(seq_at)
        )
        assert range_2m < range_no

    def test_betaine_arrhenius_correction(self):
        """Betaine Arrhenius correction at 37C, 50% GC, should be about -1.2C/M."""
        corrector = ArrheniusTmCorrector(37.0)
        correction = corrector.calculate_correction(
            "betaine", 1.0, gc_content=0.5
        )
        assert -1.5 < correction < -0.9

    def test_betaine_arrhenius_gc_ordering(self):
        """Betaine correction should be more negative for GC-rich sequences."""
        corrector = ArrheniusTmCorrector(37.0)
        corr_gc = corrector.calculate_correction(
            "betaine", 1.0, gc_content=0.6, primer_length=12
        )
        corr_50 = corrector.calculate_correction(
            "betaine", 1.0, gc_content=0.5, primer_length=12
        )
        corr_at = corrector.calculate_correction(
            "betaine", 1.0, gc_content=0.4, primer_length=12
        )
        assert corr_gc < corr_50 < corr_at

    def test_betaine_literature_range(self):
        """
        Betaine correction should match literature.

        Rees et al. (1993): -1.0 to -1.5C per M at 37C (GC-independent).
        """
        corrector = ArrheniusTmCorrector(37.0)
        correction = corrector.calculate_correction(
            "betaine", 1.0, gc_content=0.5, primer_length=12
        )
        assert -2.0 < correction < -0.8

    def test_betaine_isostabilizing_partial(self):
        """
        At 2.5M betaine (about half the isostabilizing 5.2M),
        Tm difference between GC-rich and balanced should decrease.
        """
        conditions = ReactionConditions(temp=30, betaine_m=2.5)
        cond_0 = ReactionConditions(temp=30, betaine_m=0)

        seq_gc = "GCGCGCGCGC"
        seq_50 = "GCGCATATAT"

        diff_no = abs(
            cond_0.calculate_effective_tm(seq_gc)
            - cond_0.calculate_effective_tm(seq_50)
        )
        diff_2m = abs(
            conditions.calculate_effective_tm(seq_gc)
            - conditions.calculate_effective_tm(seq_50)
        )
        assert diff_2m < diff_no
        reduction_pct = (diff_no - diff_2m) / diff_no * 100
        assert reduction_pct > 15

    def test_recalibrated_betaine_coefficient(self):
        """Betaine uses recalibrated coefficient (about -1.2C/M, not -0.5C/M)."""
        additives = AdditiveConcentrations(betaine_m=1.0)
        correction = additives.calculate_tm_correction(gc_content=0.5)
        assert -1.5 < correction < -0.8


# =============================================================================
# GC Normalization Mechanism
# =============================================================================


class TestGCNormalization:
    """Tests for the GC normalization internal mechanism."""

    def test_50_gc_no_normalization(self):
        """50% GC should have approximately zero GC normalization."""
        conditions = ReactionConditions(temp=30, betaine_m=1.0)
        correction = conditions._calculate_gc_normalization(0.5)
        assert abs(correction) < 0.01

    def test_100_gc_negative_correction(self):
        """100% GC should have negative correction."""
        conditions = ReactionConditions(temp=30, betaine_m=1.0)
        assert conditions._calculate_gc_normalization(1.0) < 0

    def test_0_gc_positive_correction(self):
        """0% GC should have positive correction."""
        conditions = ReactionConditions(temp=30, betaine_m=1.0)
        assert conditions._calculate_gc_normalization(0.0) > 0

    def test_symmetry_25_75(self):
        """25% and 75% GC corrections should be symmetric around zero."""
        conditions = ReactionConditions(temp=30, betaine_m=1.0)
        c25 = conditions._calculate_gc_normalization(0.25)
        c75 = conditions._calculate_gc_normalization(0.75)
        assert abs(c25 + c75) < 0.01

    def test_betaine_concentration_scaling(self):
        """Betaine effect should scale with concentration."""
        cond_0 = ReactionConditions(temp=30, betaine_m=0.0)
        cond_1 = ReactionConditions(temp=30, betaine_m=1.0)
        cond_2 = ReactionConditions(temp=30, betaine_m=2.0)

        gc = 0.8
        c0 = cond_0._calculate_gc_normalization(gc)
        c1 = cond_1._calculate_gc_normalization(gc)
        c2 = cond_2._calculate_gc_normalization(gc)

        assert abs(c0) < 0.01
        assert abs(c2) > abs(c1) > abs(c0)

    def test_tmac_concentration_scaling(self):
        """TMAC effect should scale with concentration."""
        cond_0 = ReactionConditions(temp=30, tmac_m=0.0)
        cond_001 = ReactionConditions(temp=30, tmac_m=0.01)
        cond_01 = ReactionConditions(temp=30, tmac_m=0.1)

        gc = 0.9
        c0 = cond_0._calculate_gc_normalization(gc)
        c001 = cond_001._calculate_gc_normalization(gc)
        c01 = cond_01._calculate_gc_normalization(gc)

        assert abs(c0) < 0.01
        assert abs(c01) > abs(c001) > abs(c0)

    def test_zero_additives_zero_correction(self):
        """Zero additive concentrations should give zero GC normalization."""
        conditions = ReactionConditions(temp=30, betaine_m=0, tmac_m=0)
        for gc in [0.0, 0.25, 0.5, 0.75, 1.0]:
            assert conditions._calculate_gc_normalization(gc) == 0

    def test_extreme_gc_symmetry(self):
        """0% and 100% GC corrections should be symmetric."""
        conditions = ReactionConditions(temp=30, betaine_m=2.0)
        c0 = conditions._calculate_gc_normalization(0.0)
        c100 = conditions._calculate_gc_normalization(1.0)
        assert c0 > 0
        assert c100 < 0
        assert abs(c0 + c100) < 0.01


# =============================================================================
# Length-Dependent GC Correction
# =============================================================================


class TestLengthDependentCorrection:
    """Tests for length-dependent GC correction scaling."""

    def test_correction_scales_with_length(self):
        """GC correction should scale linearly with primer length."""
        conditions = ReactionConditions(temp=30, betaine_m=2.0)
        gc = 0.8

        c6 = conditions._calculate_gc_normalization(gc, primer_length=6)
        c10 = conditions._calculate_gc_normalization(gc, primer_length=10)
        c18 = conditions._calculate_gc_normalization(gc, primer_length=18)

        assert all(c < 0 for c in [c6, c10, c18])
        assert abs(c18) > abs(c10) > abs(c6)

        # Check linear scaling: 2*length factor
        assert abs(c10 / c6 - 20.0 / 12.0) < 0.01
        assert abs(c18 / c10 - 36.0 / 20.0) < 0.01

    def test_effective_tm_uses_actual_length(self):
        """calculate_effective_tm should apply length-dependent correction."""
        conditions = ReactionConditions(temp=30, betaine_m=2.0)

        seq_8 = "GCGCGCGC"    # 8bp, 100% GC
        seq_16 = "GCGCGCGCGCGCGCGC"  # 16bp, 100% GC

        base_diff = (
            calculate_tm_with_salt(seq_16, 50, 0)
            - calculate_tm_with_salt(seq_8, 50, 0)
        )
        effective_diff = (
            conditions.calculate_effective_tm(seq_16)
            - conditions.calculate_effective_tm(seq_8)
        )
        assert effective_diff < base_diff

    def test_short_primer_correction(self):
        """6bp primer at 100% GC should get modest correction."""
        conditions = ReactionConditions(temp=30, betaine_m=1.0)
        corr = conditions._calculate_gc_normalization(1.0, primer_length=6)
        assert -2.0 < corr < -0.5

    def test_long_primer_correction(self):
        """18bp primer at 100% GC should get larger correction."""
        conditions = ReactionConditions(temp=30, betaine_m=1.0)
        corr = conditions._calculate_gc_normalization(1.0, primer_length=18)
        assert -5.0 < corr < -2.0


# =============================================================================
# Formamide Corrections
# =============================================================================


class TestFormamideCorrections:
    """Tests for formamide Tm corrections."""

    def test_formamide_reduces_tm(self):
        """Formamide should reduce Tm by approximately 0.65C per percent."""
        seq = "ATCGATCGATCG"
        cond_0 = ReactionConditions(temp=30, formamide_percent=0)
        cond_5 = ReactionConditions(temp=30, formamide_percent=5)

        reduction = cond_0.calculate_effective_tm(seq) - cond_5.calculate_effective_tm(seq)
        assert 2.8 < reduction < 3.7

    def test_formamide_arrhenius(self):
        """Formamide Arrhenius at 37C should be about -0.65C per %."""
        corrector = ArrheniusTmCorrector(37.0)
        correction = corrector.calculate_correction("formamide", 5.0)
        expected = -0.65 * 5.0
        assert abs(correction - expected) < 0.5

    def test_formamide_literature_range(self):
        """
        Formamide correction should match literature.

        Blake & Delcourt (1996): -0.6 to -0.72C per %.
        """
        corrector = ArrheniusTmCorrector(37.0)
        correction = corrector.calculate_correction("formamide", 10.0)
        assert -8.0 < correction < -5.0


# =============================================================================
# TMAC Corrections
# =============================================================================


class TestTMACCorrections:
    """Tests for TMAC Tm corrections."""

    def test_tmac_gc_normalization(self):
        """TMAC should affect high and low GC sequences differently."""
        tmac = ReactionConditions(polymerase="equiphi29", temp=42.0, tmac_m=0.05)
        corr_high = tmac.calculate_tm_correction(gc_content=0.7, primer_length=12)
        corr_low = tmac.calculate_tm_correction(gc_content=0.3, primer_length=12)
        assert corr_high != corr_low

    def test_tmac_gc_equalization_arrhenius(self):
        """TMAC primary effect is GC equalization via Arrhenius model."""
        corrector = ArrheniusTmCorrector(37.0)
        corr_50 = corrector.calculate_correction(
            "tmac", 0.05, gc_content=0.5, primer_length=12
        )
        corr_70 = corrector.calculate_correction(
            "tmac", 0.05, gc_content=0.7, primer_length=12
        )
        assert corr_70 < corr_50

    def test_tmac_at_practical_concentrations(self):
        """TMAC at practical concentrations should show measurable effect."""
        cond_0 = ReactionConditions(temp=30, tmac_m=0)
        cond_01 = ReactionConditions(temp=30, tmac_m=0.1)

        seq_gc = "GCGCGCGCGC"
        reduction = (
            cond_0.calculate_effective_tm(seq_gc)
            - cond_01.calculate_effective_tm(seq_gc)
        )
        assert reduction > 0.5

    def test_recalibrated_tmac_coefficient(self):
        """TMAC uses recalibrated coefficient."""
        additives = AdditiveConcentrations(tmac_m=0.1)
        correction = additives.calculate_tm_correction(gc_content=0.5)
        assert -0.5 < correction < 0.1


# =============================================================================
# Other Additives (Trehalose, Urea, Ethanol)
# =============================================================================


class TestOtherAdditives:
    """Tests for trehalose, urea, and ethanol corrections."""

    def test_trehalose_correction(self):
        """Trehalose uses recalibrated coefficient (about -3.0C/M)."""
        corrector = ArrheniusTmCorrector(37.0)
        correction = corrector.calculate_correction("trehalose", 0.5)
        expected = -3.0 * 0.5
        assert abs(correction - expected) < 0.5

    def test_trehalose_literature_range(self):
        """Trehalose should match literature: approximately 2-4C per M."""
        corrector = ArrheniusTmCorrector(37.0)
        correction = corrector.calculate_correction("trehalose", 0.5)
        assert -2.5 < correction < -0.5

    def test_urea_recalibrated_coefficient(self):
        """Urea uses recalibrated coefficient (-2.5C/M, not -5.0C/M)."""
        corrector = ArrheniusTmCorrector(37.0)
        correction = corrector.calculate_correction("urea", 1.0, gc_content=0.5)
        assert -3.5 < correction < -2.0

    def test_urea_literature_range(self):
        """Urea should match literature: -2.0 to -3.0C per M."""
        corrector = ArrheniusTmCorrector(37.0)
        correction = corrector.calculate_correction("urea", 1.0, gc_content=0.5)
        assert -4.0 < correction < -1.5

    def test_recalibrated_urea_via_additive_concentrations(self):
        """Urea via AdditiveConcentrations uses recalibrated coefficient."""
        additives = AdditiveConcentrations(urea_m=1.0)
        correction = additives.calculate_tm_correction(gc_content=0.5)
        assert -3.5 < correction < -2.0

    def test_ethanol_correction(self):
        """Ethanol should be about -0.4C per %."""
        corrector = ArrheniusTmCorrector(37.0)
        correction = corrector.calculate_correction("ethanol", 3.0)
        expected = -0.4 * 3.0
        assert abs(correction - expected) < 0.3


# =============================================================================
# Combined Additives
# =============================================================================


class TestCombinedAdditives:
    """Tests for combined additive effects."""

    def test_dmso_betaine_additive_effect(self):
        """Combined DMSO + betaine should be approximately additive."""
        seq = "ATCGATCGATCG"
        cond_none = ReactionConditions(temp=30)
        cond_dmso = ReactionConditions(temp=30, dmso_percent=5)
        cond_bet = ReactionConditions(temp=30, betaine_m=1.0)
        cond_both = ReactionConditions(temp=30, dmso_percent=5, betaine_m=1.0)

        tm_none = cond_none.calculate_effective_tm(seq)
        dmso_effect = tm_none - cond_dmso.calculate_effective_tm(seq)
        bet_effect = tm_none - cond_bet.calculate_effective_tm(seq)
        combined_effect = tm_none - cond_both.calculate_effective_tm(seq)

        expected = dmso_effect + bet_effect
        assert abs(combined_effect - expected) < 0.5

    def test_dmso_betaine_synergy(self):
        """Combined DMSO + betaine should exceed either alone."""
        ctrl = ReactionConditions(polymerase="equiphi29", temp=42.0)
        combined = ReactionConditions(
            polymerase="equiphi29", temp=42.0, dmso_percent=5.0, betaine_m=1.0
        )
        ctrl_corr = ctrl.calculate_tm_correction(gc_content=0.5, primer_length=12)
        comb_corr = combined.calculate_tm_correction(gc_content=0.5, primer_length=12)
        assert comb_corr < ctrl_corr

    def test_all_additives_moderate(self):
        """Multiple additives at moderate concentrations."""
        conditions = ReactionConditions(
            polymerase="equiphi29",
            temp=45.0,
            dmso_percent=5.0,
            betaine_m=1.0,
            trehalose_m=0.3,
        )
        correction = conditions.calculate_tm_correction(gc_content=0.5, primer_length=15)
        assert correction <= -4.9

    def test_combined_tmac_betaine(self):
        """TMAC + betaine together should give at least as much effect as either alone."""
        seq_gc = "GCGCGCGCGC"
        cond_none = ReactionConditions(temp=30, betaine_m=0, tmac_m=0)
        cond_bet = ReactionConditions(temp=30, betaine_m=1.0, tmac_m=0)
        cond_tmac = ReactionConditions(temp=30, betaine_m=0, tmac_m=0.1)
        cond_both = ReactionConditions(temp=30, betaine_m=1.0, tmac_m=0.1)

        tm_none = cond_none.calculate_effective_tm(seq_gc)
        red_bet = tm_none - cond_bet.calculate_effective_tm(seq_gc)
        red_tmac = tm_none - cond_tmac.calculate_effective_tm(seq_gc)
        red_both = tm_none - cond_both.calculate_effective_tm(seq_gc)

        assert red_both >= max(red_bet, red_tmac) - 0.5

    def test_arrhenius_total_correction(self):
        """Arrhenius total correction should be sum of individual corrections."""
        corrector = ArrheniusTmCorrector(42.0)
        additives = {"dmso": 5.0, "betaine": 1.0}
        total = corrector.calculate_total_correction(additives, gc_content=0.5)

        dmso_corr = corrector.calculate_correction("dmso", 5.0, gc_content=0.5)
        bet_corr = corrector.calculate_correction("betaine", 1.0, gc_content=0.5)
        assert abs(total - (dmso_corr + bet_corr)) < 0.01

    def test_combined_extremes_reasonable(self):
        """Combined additives at moderate concentrations should stay reasonable."""
        additives = AdditiveConcentrations(
            dmso_percent=5.0, betaine_m=1.5, formamide_percent=2.0
        )
        correction = additives.calculate_tm_correction(
            gc_content=0.5, reaction_temp_celsius=42.0
        )
        assert -15.0 < correction < 0.0


# =============================================================================
# GC-Rich / AT-Rich Sequence-Level Tests
# =============================================================================


class TestGCRichSequences:
    """Tests for additive corrections on GC-rich sequences."""

    def test_betaine_reduces_tm_of_gc_rich(self):
        """Betaine should substantially reduce Tm of GC-rich sequences."""
        seq = "GCGCGCGCGC"
        assert gc_content(seq) == 1.0

        cond_0 = ReactionConditions(temp=30, betaine_m=0)
        cond_1 = ReactionConditions(temp=30, betaine_m=1.0)
        cond_2 = ReactionConditions(temp=30, betaine_m=2.0)

        tm_0 = cond_0.calculate_effective_tm(seq)
        tm_1 = cond_1.calculate_effective_tm(seq)
        tm_2 = cond_2.calculate_effective_tm(seq)

        assert tm_2 < tm_1 < tm_0
        reduction = tm_0 - tm_1
        assert 1.5 < reduction < 4.0

    def test_high_gc_preset_effectiveness(self):
        """High GC preset should reduce Tm of GC-rich sequences."""
        seq = "GCGCGCGCGCGC"
        standard = get_standard_conditions()
        high_gc = get_high_gc_conditions()

        assert high_gc.calculate_effective_tm(seq) < standard.calculate_effective_tm(seq)
        assert high_gc.betaine_m >= 2.0
        assert high_gc.dmso_percent >= 5.0


class TestATRichSequences:
    """Tests for additive corrections on AT-rich sequences."""

    def test_betaine_modest_effect_on_at_rich(self):
        """Betaine effect on AT-rich should be modest (may even increase Tm)."""
        seq = "AAAAAATTTT"
        assert gc_content(seq) == 0.0

        cond_0 = ReactionConditions(temp=30, betaine_m=0)
        cond_1 = ReactionConditions(temp=30, betaine_m=1.0)

        effect = (
            cond_0.calculate_effective_tm(seq) - cond_1.calculate_effective_tm(seq)
        )
        assert effect < 2.0


class TestBalancedSequences:
    """Tests for additive corrections on balanced (50% GC) sequences."""

    def test_balanced_modest_betaine_effect(self):
        """50% GC sequences should show modest effect for betaine."""
        seq = "ATCGATCGATCG"
        assert abs(gc_content(seq) - 0.5) < 0.01

        cond_0 = ReactionConditions(temp=30, betaine_m=0)
        cond_2 = ReactionConditions(temp=30, betaine_m=2.0)

        reduction = (
            cond_0.calculate_effective_tm(seq) - cond_2.calculate_effective_tm(seq)
        )
        assert 0.5 < reduction < 3.0


# =============================================================================
# Extreme GC Conditions Preset
# =============================================================================


class TestExtremeGCPreset:
    """Tests for the extreme GC conditions preset."""

    def test_has_tmac_betaine_urea(self):
        """Extreme GC preset should have TMAC, betaine, and urea."""
        extreme = get_extreme_gc_conditions()
        assert extreme.tmac_m > 0
        assert extreme.betaine_m >= 2.0
        assert extreme.urea_m > 0

    def test_effectiveness_on_extreme_gc(self):
        """Extreme GC preset should substantially reduce Tm of GC-rich sequences."""
        seq = "GCGCGCGCGCGCGCGC"
        standard = get_standard_conditions()
        extreme = get_extreme_gc_conditions()

        tm_std = standard.calculate_effective_tm(seq)
        tm_ext = extreme.calculate_effective_tm(seq)
        assert tm_ext < tm_std - 5


# =============================================================================
# Arrhenius Temperature Behavior
# =============================================================================


class TestArrheniusTemperatureBehavior:
    """Tests for Arrhenius model temperature dependence."""

    def test_init(self):
        """Test initialization."""
        corrector = ArrheniusTmCorrector(30.0)
        assert corrector.reaction_temp_celsius == 30.0
        assert abs(corrector.reaction_temp_k - 303.15) < 0.01

    def test_higher_temp_larger_effect(self):
        """Higher temperatures should produce larger destabilizing effects."""
        temps = [25, 30, 37, 42, 45]
        corrections = [
            ArrheniusTmCorrector(float(t)).calculate_correction("dmso", 5.0)
            for t in temps
        ]
        for c in corrections:
            assert c < 0
        for i in range(len(corrections) - 1):
            assert corrections[i + 1] <= corrections[i] * 0.98

    def test_reasonable_temperature_sensitivity(self):
        """Temperature sensitivity ratio between 30C and 42C should be moderate."""
        _, _, ratio = ArrheniusTmCorrector(37.0).compare_temperatures(
            "dmso", 5.0, 30.0, 42.0
        )
        assert 1.0 < ratio < 1.5

    def test_compare_temperatures(self):
        """compare_temperatures should return coherent values."""
        corr1, corr2, ratio = ArrheniusTmCorrector(37.0).compare_temperatures(
            "dmso", 5.0, 30.0, 42.0
        )
        assert corr1 < 0
        assert corr2 < 0
        assert abs(corr2) > abs(corr1)
        assert ratio > 1

    def test_unknown_additive_raises(self):
        """Unknown additive should raise ValueError."""
        corrector = ArrheniusTmCorrector(37.0)
        with pytest.raises(ValueError, match="Unknown additive"):
            corrector.calculate_correction("invalid_additive", 1.0)


# =============================================================================
# Additive Parameters Structure
# =============================================================================


class TestAdditiveParameters:
    """Tests for ADDITIVE_TM_PARAMS structure and helpers."""

    def test_all_required_fields(self):
        """All additives should have required parameter fields."""
        required = [
            "ref_coef", "ref_temp", "activation_energy",
            "max_concentration", "gc_dependent",
        ]
        for additive, params in ADDITIVE_TM_PARAMS.items():
            for field in required:
                assert field in params, f"{additive} missing {field}"

    def test_activation_energies_positive(self):
        """All activation energies should be positive."""
        for additive, params in ADDITIVE_TM_PARAMS.items():
            assert params["activation_energy"] > 0

    def test_ref_temps_reasonable(self):
        """Reference temperatures should be between 25C and 40C."""
        for additive, params in ADDITIVE_TM_PARAMS.items():
            ref_temp_c = params["ref_temp"] - 273.15
            assert 25 <= ref_temp_c <= 40

    def test_get_additive_tm_params(self):
        """get_additive_tm_params should return correct values."""
        params = get_additive_tm_params("dmso")
        assert params["ref_coef"] == -0.55
        assert params["activation_energy"] == 2500.0

    def test_get_additive_tm_params_case_insensitive(self):
        """Helper should be case-insensitive."""
        assert get_additive_tm_params("DMSO") == get_additive_tm_params("dmso")

    def test_get_additive_tm_params_unknown_raises(self):
        """Unknown additive should raise ValueError."""
        with pytest.raises(ValueError, match="Unknown additive"):
            get_additive_tm_params("unknown")

    def test_list_additives(self):
        """list_additives should return all known additives."""
        additives = list_additives()
        assert "dmso" in additives
        assert "betaine" in additives
        assert len(additives) == len(ADDITIVE_TM_PARAMS)

    def test_zero_concentration_gives_zero_correction(self):
        """Zero concentration should give zero correction for all additives."""
        corrector = ArrheniusTmCorrector(37.0)
        for additive in ADDITIVE_TM_PARAMS:
            assert corrector.calculate_correction(additive, 0.0) == 0.0

    def test_max_concentration_reasonable(self):
        """Max concentrations should produce reasonable corrections (< 20C)."""
        corrector = ArrheniusTmCorrector(37.0)
        for additive, params in ADDITIVE_TM_PARAMS.items():
            correction = corrector.calculate_correction(
                additive, params["max_concentration"], gc_content=0.5
            )
            assert correction < 0
            assert correction > -20.0


# =============================================================================
# ReactionConditions with Arrhenius
# =============================================================================


class TestReactionConditionsArrhenius:
    """Tests for ReactionConditions with Arrhenius-based corrections."""

    def test_default_uses_arrhenius(self):
        """Default behavior should use Arrhenius corrections."""
        conditions = ReactionConditions(
            temp=42.0, dmso_percent=5.0, polymerase="equiphi29"
        )
        correction = conditions.calculate_tm_correction(gc_content=0.5)
        arrhenius_corr = conditions.additives.calculate_tm_correction(
            gc_content=0.5, reaction_temp_celsius=42.0
        )
        assert abs(correction - arrhenius_corr) < 0.01

    def test_can_disable_arrhenius(self):
        """Arrhenius can be disabled for backward compatibility."""
        conditions = ReactionConditions(
            temp=30.0, dmso_percent=5.0, polymerase="phi29"
        )
        arrhenius = conditions.calculate_tm_correction(
            gc_content=0.5, use_arrhenius=True
        )
        fixed = conditions.calculate_tm_correction(
            gc_content=0.5, use_arrhenius=False
        )
        assert abs(arrhenius - fixed) < 1.0

    def test_phi29_vs_equiphi29_temperature_effect(self):
        """Additive effects should differ between phi29 and equiphi29."""
        phi29 = ReactionConditions(
            temp=30.0, dmso_percent=5.0, betaine_m=1.0, polymerase="phi29"
        )
        equiphi29 = ReactionConditions(
            temp=42.0, dmso_percent=5.0, betaine_m=1.0, polymerase="equiphi29"
        )
        assert abs(
            phi29.calculate_tm_correction(gc_content=0.5)
            - equiphi29.calculate_tm_correction(gc_content=0.5)
        ) > 0.1

    def test_fixed_vs_arrhenius_similar_at_reference(self):
        """At reference temp (37C), fixed and Arrhenius should be similar."""
        additives = AdditiveConcentrations(dmso_percent=5.0)
        fixed = additives.calculate_tm_correction(gc_content=0.5)
        arrhenius = additives.calculate_tm_correction(
            gc_content=0.5, reaction_temp_celsius=37.0
        )
        assert abs(fixed - arrhenius) < abs(fixed) * 0.2

    def test_arrhenius_varies_with_temperature(self):
        """Arrhenius correction should vary with temperature."""
        additives = AdditiveConcentrations(dmso_percent=5.0, betaine_m=1.0)
        corr_30 = additives.calculate_tm_correction(
            gc_content=0.5, reaction_temp_celsius=30.0
        )
        corr_42 = additives.calculate_tm_correction(
            gc_content=0.5, reaction_temp_celsius=42.0
        )
        assert corr_42 != corr_30


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
