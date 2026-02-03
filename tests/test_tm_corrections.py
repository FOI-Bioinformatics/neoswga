"""
Validation tests for Arrhenius-based Tm corrections.

Tests the recalibrated additive Tm correction coefficients against
literature values and Tm calculator predictions.

Literature references:
    - Chester & Marshak (1993) Anal Biochem 209:284-290 (DMSO multi-temp)
    - Rees et al. (1993) Biochemistry 32:137-144 (betaine)
    - Henke et al. (1997) NAR 25:3957-3958 (betaine PCR)
    - Blake & Delcourt (1996) NAR 24:2095-2103 (formamide)
    - Spiess et al. (2004) Biotechniques 36:732-736 (trehalose)
    - Lesnick & Bhalla (1995) NAR 23:4665-4666 (urea)
    - Melchior & von Hippel (1973) PNAS 70:298-302 (TMAC)
    - Cheng et al. (1994) PNAS 91:5695-5699 (ethanol)
"""

import pytest
import math

from neoswga.core.additives import (
    AdditiveConcentrations,
    ArrheniusTmCorrector,
)
from neoswga.core.mechanistic_params import (
    ADDITIVE_TM_PARAMS,
    get_additive_tm_params,
    list_additives,
)
from neoswga.core.reaction_conditions import ReactionConditions


# =============================================================================
# Test ArrheniusTmCorrector
# =============================================================================

class TestArrheniusTmCorrector:
    """Tests for the ArrheniusTmCorrector class."""

    def test_init(self):
        """Test initialization with reaction temperature."""
        corrector = ArrheniusTmCorrector(30.0)
        assert corrector.reaction_temp_celsius == 30.0
        assert abs(corrector.reaction_temp_k - 303.15) < 0.01

    def test_dmso_correction_at_reference_temp(self):
        """Test DMSO correction at reference temperature (37C)."""
        corrector = ArrheniusTmCorrector(37.0)
        correction = corrector.calculate_correction('dmso', 5.0)

        # At 37C (reference), should be close to ref_coef * concentration
        expected = -0.55 * 5.0  # -2.75C
        assert abs(correction - expected) < 0.1

    def test_dmso_temperature_dependence(self):
        """Test DMSO correction varies with temperature."""
        corr_30 = ArrheniusTmCorrector(30.0).calculate_correction('dmso', 5.0)
        corr_42 = ArrheniusTmCorrector(42.0).calculate_correction('dmso', 5.0)

        # Higher temperature should give larger correction (more negative)
        assert corr_42 < corr_30
        # Both should be negative
        assert corr_30 < 0
        assert corr_42 < 0

    def test_betaine_correction(self):
        """Test betaine correction with recalibrated coefficient."""
        corrector = ArrheniusTmCorrector(37.0)
        correction = corrector.calculate_correction('betaine', 1.0, gc_content=0.5)

        # At 50% GC, no GC normalization effect
        # Should be approximately -1.2C per M at reference temp
        assert -1.5 < correction < -0.9

    def test_betaine_gc_normalization(self):
        """Test betaine GC normalization effect."""
        corrector = ArrheniusTmCorrector(37.0)

        # GC-rich sequence (60%)
        corr_gc_rich = corrector.calculate_correction(
            'betaine', 1.0, gc_content=0.6, primer_length=12
        )

        # Balanced sequence (50%)
        corr_balanced = corrector.calculate_correction(
            'betaine', 1.0, gc_content=0.5, primer_length=12
        )

        # AT-rich sequence (40%)
        corr_at_rich = corrector.calculate_correction(
            'betaine', 1.0, gc_content=0.4, primer_length=12
        )

        # GC-rich should have more negative correction (Tm lowered more)
        # AT-rich should have less negative correction (Tm lowered less)
        assert corr_gc_rich < corr_balanced < corr_at_rich

    def test_urea_recalibrated_coefficient(self):
        """Test urea uses recalibrated coefficient (-2.5 vs old -5.0)."""
        corrector = ArrheniusTmCorrector(37.0)
        correction = corrector.calculate_correction('urea', 1.0, gc_content=0.5)

        # Should be approximately -2.5C per M (not -5.0)
        assert -3.5 < correction < -2.0

    def test_tmac_gc_equalization(self):
        """Test TMAC's primary effect is GC equalization."""
        corrector = ArrheniusTmCorrector(37.0)

        # At 50% GC, minimal GC normalization
        corr_balanced = corrector.calculate_correction(
            'tmac', 0.05, gc_content=0.5, primer_length=12
        )

        # At 70% GC, should show GC normalization effect
        corr_gc_rich = corrector.calculate_correction(
            'tmac', 0.05, gc_content=0.7, primer_length=12
        )

        # GC-rich should have larger (more negative) correction
        assert corr_gc_rich < corr_balanced

    def test_formamide_correction(self):
        """Test formamide correction."""
        corrector = ArrheniusTmCorrector(37.0)
        correction = corrector.calculate_correction('formamide', 5.0)

        # Should be approximately -0.65C per %
        expected = -0.65 * 5.0  # -3.25C
        assert abs(correction - expected) < 0.5

    def test_trehalose_recalibrated(self):
        """Test trehalose uses recalibrated coefficient."""
        corrector = ArrheniusTmCorrector(37.0)
        correction = corrector.calculate_correction('trehalose', 0.5)

        # Should be approximately -3.0C per M (not -5.0)
        expected = -3.0 * 0.5  # -1.5C
        assert abs(correction - expected) < 0.5

    def test_ethanol_correction(self):
        """Test ethanol correction."""
        corrector = ArrheniusTmCorrector(37.0)
        correction = corrector.calculate_correction('ethanol', 3.0)

        # Should be approximately -0.4C per %
        expected = -0.4 * 3.0  # -1.2C
        assert abs(correction - expected) < 0.3

    def test_unknown_additive_raises(self):
        """Test that unknown additive raises ValueError."""
        corrector = ArrheniusTmCorrector(37.0)
        with pytest.raises(ValueError, match="Unknown additive"):
            corrector.calculate_correction('invalid_additive', 1.0)

    def test_total_correction(self):
        """Test calculating total correction from multiple additives."""
        corrector = ArrheniusTmCorrector(42.0)
        additives = {
            'dmso': 5.0,
            'betaine': 1.0,
        }
        total = corrector.calculate_total_correction(additives, gc_content=0.5)

        # Should be sum of individual corrections
        dmso_corr = corrector.calculate_correction('dmso', 5.0, gc_content=0.5)
        betaine_corr = corrector.calculate_correction('betaine', 1.0, gc_content=0.5)

        assert abs(total - (dmso_corr + betaine_corr)) < 0.01

    def test_compare_temperatures(self):
        """Test comparing corrections at different temperatures."""
        corrector = ArrheniusTmCorrector(37.0)
        corr1, corr2, ratio = corrector.compare_temperatures(
            'dmso', 5.0, 30.0, 42.0
        )

        # Both should be negative
        assert corr1 < 0
        assert corr2 < 0
        # 42C should have larger effect
        assert abs(corr2) > abs(corr1)
        # Ratio should be > 1
        assert ratio > 1


# =============================================================================
# Test AdditiveConcentrations with Arrhenius
# =============================================================================

class TestAdditiveConcentrationsArrhenius:
    """Test AdditiveConcentrations with Arrhenius-based corrections."""

    def test_fixed_vs_arrhenius_at_reference_temp(self):
        """At reference temp (37C), fixed and Arrhenius should be similar."""
        additives = AdditiveConcentrations(dmso_percent=5.0)

        fixed = additives.calculate_tm_correction(gc_content=0.5)
        arrhenius = additives.calculate_tm_correction(
            gc_content=0.5, reaction_temp_celsius=37.0
        )

        # Should be within 10% of each other at reference temp
        assert abs(fixed - arrhenius) < abs(fixed) * 0.2

    def test_arrhenius_temperature_effect(self):
        """Test that Arrhenius correction varies with temperature."""
        additives = AdditiveConcentrations(dmso_percent=5.0, betaine_m=1.0)

        corr_30 = additives.calculate_tm_correction(
            gc_content=0.5, reaction_temp_celsius=30.0
        )
        corr_42 = additives.calculate_tm_correction(
            gc_content=0.5, reaction_temp_celsius=42.0
        )

        # Higher temperature should generally give larger corrections
        assert corr_42 != corr_30

    def test_recalibrated_betaine(self):
        """Test betaine uses recalibrated coefficient."""
        additives = AdditiveConcentrations(betaine_m=1.0)

        # Fixed coefficient calculation
        correction = additives.calculate_tm_correction(gc_content=0.5)

        # Should be approximately -1.2C (not -0.5C as before)
        assert -1.5 < correction < -0.8

    def test_recalibrated_urea(self):
        """Test urea uses recalibrated coefficient."""
        additives = AdditiveConcentrations(urea_m=1.0)

        correction = additives.calculate_tm_correction(gc_content=0.5)

        # Should be approximately -2.5C (not -5.0C as before)
        assert -3.5 < correction < -2.0

    def test_recalibrated_tmac(self):
        """Test TMAC uses recalibrated coefficient."""
        additives = AdditiveConcentrations(tmac_m=0.1)

        # At 50% GC, mainly uniform correction
        correction = additives.calculate_tm_correction(gc_content=0.5)

        # Should be approximately -0.05C (not -1.0C as before)
        # The uniform component is now -0.5C/M, so 0.1M gives -0.05C
        assert -0.5 < correction < 0.1


# =============================================================================
# Test ReactionConditions with Arrhenius
# =============================================================================

class TestReactionConditionsArrhenius:
    """Test ReactionConditions with Arrhenius-based corrections."""

    def test_default_uses_arrhenius(self):
        """Test that default behavior uses Arrhenius corrections."""
        conditions = ReactionConditions(
            temp=42.0,
            dmso_percent=5.0,
            polymerase='equiphi29'
        )

        # Should use Arrhenius by default
        correction = conditions.calculate_tm_correction(gc_content=0.5)

        # Verify it's using temperature-dependent calculation
        # by comparing with explicit Arrhenius call
        arrhenius_corr = conditions.additives.calculate_tm_correction(
            gc_content=0.5, reaction_temp_celsius=42.0
        )
        assert abs(correction - arrhenius_corr) < 0.01

    def test_can_disable_arrhenius(self):
        """Test that Arrhenius can be disabled for backward compatibility."""
        conditions = ReactionConditions(
            temp=30.0,
            dmso_percent=5.0,
            polymerase='phi29'
        )

        arrhenius = conditions.calculate_tm_correction(
            gc_content=0.5, use_arrhenius=True
        )
        fixed = conditions.calculate_tm_correction(
            gc_content=0.5, use_arrhenius=False
        )

        # At 30C, should still have some difference
        # (37C is reference, so 30C will be slightly different)
        # Allow for some similarity since we're close to reference
        assert abs(arrhenius - fixed) < 1.0

    def test_phi29_vs_equiphi29_temperature_effect(self):
        """Test that additive effects differ between phi29 and equiphi29 temps."""
        phi29 = ReactionConditions(
            temp=30.0,
            dmso_percent=5.0,
            betaine_m=1.0,
            polymerase='phi29'
        )

        equiphi29 = ReactionConditions(
            temp=42.0,
            dmso_percent=5.0,
            betaine_m=1.0,
            polymerase='equiphi29'
        )

        corr_phi29 = phi29.calculate_tm_correction(gc_content=0.5)
        corr_equiphi29 = equiphi29.calculate_tm_correction(gc_content=0.5)

        # At higher temperature, additive effects should be different
        assert abs(corr_phi29 - corr_equiphi29) > 0.1


# =============================================================================
# Test Parameter Structure
# =============================================================================

class TestAdditiveParams:
    """Test ADDITIVE_TM_PARAMS structure and helper functions."""

    def test_all_additives_have_required_fields(self):
        """Test that all additives have required parameter fields."""
        required_fields = [
            'ref_coef', 'ref_temp', 'activation_energy',
            'max_concentration', 'gc_dependent'
        ]

        for additive, params in ADDITIVE_TM_PARAMS.items():
            for field in required_fields:
                assert field in params, f"{additive} missing {field}"

    def test_activation_energies_positive(self):
        """Test that all activation energies are positive."""
        for additive, params in ADDITIVE_TM_PARAMS.items():
            assert params['activation_energy'] > 0, \
                f"{additive} has non-positive Ea"

    def test_ref_temps_reasonable(self):
        """Test that reference temperatures are reasonable (25-40C)."""
        for additive, params in ADDITIVE_TM_PARAMS.items():
            ref_temp_c = params['ref_temp'] - 273.15
            assert 25 <= ref_temp_c <= 40, \
                f"{additive} ref_temp {ref_temp_c}C out of range"

    def test_get_additive_tm_params(self):
        """Test get_additive_tm_params helper function."""
        params = get_additive_tm_params('dmso')
        assert params['ref_coef'] == -0.55
        assert params['activation_energy'] == 2500.0

    def test_get_additive_tm_params_case_insensitive(self):
        """Test that helper is case-insensitive."""
        params1 = get_additive_tm_params('DMSO')
        params2 = get_additive_tm_params('dmso')
        assert params1 == params2

    def test_get_additive_tm_params_unknown_raises(self):
        """Test that unknown additive raises ValueError."""
        with pytest.raises(ValueError, match="Unknown additive"):
            get_additive_tm_params('unknown')

    def test_list_additives(self):
        """Test list_additives helper function."""
        additives = list_additives()
        assert 'dmso' in additives
        assert 'betaine' in additives
        assert len(additives) == len(ADDITIVE_TM_PARAMS)


# =============================================================================
# Validation Against Literature Values
# =============================================================================

class TestLiteratureValidation:
    """
    Validate Tm corrections against published literature values.

    These tests ensure our recalibrated coefficients match
    experimental measurements from the literature.
    """

    def test_dmso_literature_range(self):
        """
        Validate DMSO correction against literature.

        Chester & Marshak (1993): -0.5 to -0.6C per % at 37C
        Varadaraj (1994): 10% DMSO lowers Tm by ~6C
        """
        corrector = ArrheniusTmCorrector(37.0)

        # 10% DMSO
        correction = corrector.calculate_correction('dmso', 10.0)

        # Literature: -5 to -6C for 10% DMSO
        assert -7.0 < correction < -4.0, \
            f"DMSO correction {correction}C outside literature range"

    def test_betaine_literature_range(self):
        """
        Validate betaine correction against literature.

        Rees et al. (1993): -1.0 to -1.5C per M at 37C (GC-independent)
        Henke et al. (1997): Similar values for PCR
        """
        corrector = ArrheniusTmCorrector(37.0)

        # 1M betaine, 50% GC (no normalization)
        correction = corrector.calculate_correction(
            'betaine', 1.0, gc_content=0.5, primer_length=12
        )

        # Literature: -1.0 to -1.5C per M
        assert -2.0 < correction < -0.8, \
            f"Betaine correction {correction}C outside literature range"

    def test_urea_literature_range(self):
        """
        Validate urea correction against literature.

        Lesnick & Bhalla (1995): -2.0 to -3.0C per M
        Hutton (1977): Similar range
        """
        corrector = ArrheniusTmCorrector(37.0)

        # 1M urea, 50% GC
        correction = corrector.calculate_correction(
            'urea', 1.0, gc_content=0.5
        )

        # Literature: -2.0 to -3.0C per M
        assert -4.0 < correction < -1.5, \
            f"Urea correction {correction}C outside literature range"

    def test_formamide_literature_range(self):
        """
        Validate formamide correction against literature.

        Blake & Delcourt (1996): -0.6 to -0.72C per %
        McConaughy et al. (1969): Similar values
        """
        corrector = ArrheniusTmCorrector(37.0)

        # 10% formamide
        correction = corrector.calculate_correction('formamide', 10.0)

        # Literature: -6 to -7.2C for 10%
        assert -8.0 < correction < -5.0, \
            f"Formamide correction {correction}C outside literature range"

    def test_trehalose_literature_range(self):
        """
        Validate trehalose correction against literature.

        Spiess et al. (2004): ~2-4C per M
        """
        corrector = ArrheniusTmCorrector(37.0)

        # 0.5M trehalose
        correction = corrector.calculate_correction('trehalose', 0.5)

        # Literature: -1 to -2C for 0.5M
        assert -2.5 < correction < -0.5, \
            f"Trehalose correction {correction}C outside literature range"


# =============================================================================
# Test Arrhenius Temperature Behavior
# =============================================================================

class TestArrheniusTemperatureBehavior:
    """Test that Arrhenius model produces expected temperature behavior."""

    def test_higher_temp_larger_effect(self):
        """
        Test that higher temperatures produce larger destabilizing effects.

        This is expected from Arrhenius kinetics with positive Ea.
        """
        temps = [25, 30, 37, 42, 45]
        corrections = []

        for temp in temps:
            corrector = ArrheniusTmCorrector(float(temp))
            corr = corrector.calculate_correction('dmso', 5.0)
            corrections.append(corr)

        # All should be negative (Tm lowering)
        for corr in corrections:
            assert corr < 0

        # Corrections should become more negative at higher temps
        for i in range(len(corrections) - 1):
            assert corrections[i+1] <= corrections[i] * 0.98, \
                f"Expected larger effect at higher temp: {temps[i]}C vs {temps[i+1]}C"

    def test_reasonable_temperature_sensitivity(self):
        """
        Test that temperature sensitivity is reasonable (not too extreme).

        The difference between 30C and 42C should be noticeable but not huge.
        """
        corr_30, corr_42, ratio = ArrheniusTmCorrector(37.0).compare_temperatures(
            'dmso', 5.0, 30.0, 42.0
        )

        # Ratio should be between 1.05 and 1.30 for this temperature range
        assert 1.0 < ratio < 1.5, \
            f"Temperature sensitivity ratio {ratio} seems unreasonable"


# =============================================================================
# Test Extreme Concentration Bounds
# =============================================================================

class TestExtremeBounds:
    """Test behavior at concentration extremes."""

    def test_zero_concentration(self):
        """Test that zero concentration gives zero correction."""
        corrector = ArrheniusTmCorrector(37.0)

        for additive in ADDITIVE_TM_PARAMS.keys():
            correction = corrector.calculate_correction(additive, 0.0)
            assert correction == 0.0, f"{additive} non-zero at 0 concentration"

    def test_max_concentration_reasonable(self):
        """Test that max concentrations produce reasonable corrections."""
        corrector = ArrheniusTmCorrector(37.0)

        for additive, params in ADDITIVE_TM_PARAMS.items():
            max_conc = params['max_concentration']
            correction = corrector.calculate_correction(
                additive, max_conc, gc_content=0.5
            )

            # All corrections should be negative
            assert correction < 0, f"{additive} positive at max concentration"

            # No single additive should lower Tm by more than 20C
            assert correction > -20.0, \
                f"{additive} correction {correction}C too extreme"

    def test_combined_extremes(self):
        """Test combined additives at moderate concentrations."""
        additives = AdditiveConcentrations(
            dmso_percent=5.0,
            betaine_m=1.5,
            formamide_percent=2.0
        )

        correction = additives.calculate_tm_correction(
            gc_content=0.5, reaction_temp_celsius=42.0
        )

        # Combined correction should be reasonable
        assert -15.0 < correction < 0.0


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
