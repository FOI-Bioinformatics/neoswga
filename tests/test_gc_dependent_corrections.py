"""
Integration tests for GC-dependent Tm corrections (TMAC and betaine).

Tests the scientifically accurate implementation of:
1. TMAC AT/GC equalization effect (Melchior & von Hippel 1973)
2. Betaine GC normalization effect (Rees et al. 1993)
3. Combined effects under realistic SWGA conditions

Scientific basis:
- TMAC at 3M: complete equalization of AT and GC contributions to Tm
- Betaine at 5.2M: complete equalization ("isostabilizing concentration")
- At lower concentrations: proportional partial equalization
- Effect: GC-rich sequences get reduced Tm, AT-rich sequences get increased Tm

References:
- Melchior & von Hippel (1973) PNAS 70:298-302
- Rees et al. (1993) Biochemistry 32:137-144
- Henke et al. (1997) Nucleic Acids Res 25:3957-3958
"""

import pytest
import numpy as np
from neoswga.core.reaction_conditions import (
    ReactionConditions,
    get_standard_conditions,
    get_enhanced_conditions,
    get_high_gc_conditions,
    get_extreme_gc_conditions,
)
from neoswga.core.thermodynamics import calculate_tm_with_salt, gc_content


class TestGCNormalizationMechanism:
    """Test the fundamental GC normalization mechanism."""

    def test_gc_normalization_calculation(self):
        """Test the _calculate_gc_normalization internal method."""
        conditions = ReactionConditions(temp=30, betaine_m=1.0)

        # 50% GC: no normalization effect
        correction_50 = conditions._calculate_gc_normalization(0.5)
        assert abs(correction_50) < 0.01, \
            f"50% GC should have ~0 correction, got {correction_50:.3f}"

        # 100% GC: maximum negative correction
        correction_100 = conditions._calculate_gc_normalization(1.0)
        assert correction_100 < 0, \
            f"100% GC should have negative correction, got {correction_100:.3f}"

        # 0% GC: maximum positive correction
        correction_0 = conditions._calculate_gc_normalization(0.0)
        assert correction_0 > 0, \
            f"0% GC should have positive correction, got {correction_0:.3f}"

        # Symmetry: 25% and 75% GC should have opposite corrections
        correction_25 = conditions._calculate_gc_normalization(0.25)
        correction_75 = conditions._calculate_gc_normalization(0.75)
        assert abs(correction_25 + correction_75) < 0.01, \
            f"25% ({correction_25:.2f}) and 75% ({correction_75:.2f}) should be symmetric"

    def test_betaine_concentration_scaling(self):
        """Betaine effect scales linearly up to isostabilizing concentration."""
        conditions_0m = ReactionConditions(temp=30, betaine_m=0.0)
        conditions_1m = ReactionConditions(temp=30, betaine_m=1.0)
        conditions_2m = ReactionConditions(temp=30, betaine_m=2.0)

        # GC-rich sequence to see clear effect
        gc_content = 0.8

        corr_0m = conditions_0m._calculate_gc_normalization(gc_content)
        corr_1m = conditions_1m._calculate_gc_normalization(gc_content)
        corr_2m = conditions_2m._calculate_gc_normalization(gc_content)

        # 0M betaine: no GC normalization
        assert abs(corr_0m) < 0.01, f"0M betaine should have no effect"

        # Higher concentration = stronger effect
        assert abs(corr_2m) > abs(corr_1m) > abs(corr_0m), \
            f"Effect should increase with concentration: {corr_0m:.2f} < {corr_1m:.2f} < {corr_2m:.2f}"

    def test_tmac_concentration_scaling(self):
        """TMAC effect scales linearly up to full equalization at 3M."""
        conditions_0m = ReactionConditions(temp=30, tmac_m=0.0)
        conditions_001m = ReactionConditions(temp=30, tmac_m=0.01)
        conditions_01m = ReactionConditions(temp=30, tmac_m=0.1)

        # GC-rich sequence
        gc_content = 0.9

        corr_0m = conditions_0m._calculate_gc_normalization(gc_content)
        corr_001m = conditions_001m._calculate_gc_normalization(gc_content)
        corr_01m = conditions_01m._calculate_gc_normalization(gc_content)

        # 0M TMAC: no GC normalization
        assert abs(corr_0m) < 0.01, f"0M TMAC should have no effect"

        # Higher concentration = stronger effect
        assert abs(corr_01m) > abs(corr_001m) > abs(corr_0m), \
            "Effect should increase with TMAC concentration"


class TestGCRichSequences:
    """Test GC-dependent corrections on GC-rich sequences."""

    def test_betaine_reduces_tm_of_gc_rich(self):
        """Betaine should significantly reduce Tm of GC-rich sequences."""
        seq_gc_rich = 'GCGCGCGCGC'  # 100% GC, 10bp
        assert gc_content(seq_gc_rich) == 1.0

        conditions_no = ReactionConditions(temp=30, betaine_m=0)
        conditions_1m = ReactionConditions(temp=30, betaine_m=1.0)
        conditions_2m = ReactionConditions(temp=30, betaine_m=2.0)

        tm_no = conditions_no.calculate_effective_tm(seq_gc_rich)
        tm_1m = conditions_1m.calculate_effective_tm(seq_gc_rich)
        tm_2m = conditions_2m.calculate_effective_tm(seq_gc_rich)

        # GC-rich: substantial reduction
        assert tm_1m < tm_no, "1M betaine should reduce Tm of GC-rich"
        assert tm_2m < tm_1m < tm_no, "Higher betaine = lower Tm for GC-rich"

        # Magnitude check: 100% GC deviation of 0.5, so at 1M/5.2M equalization
        # expect roughly 0.5 * 20 * (1/5.2) = ~1.9C plus uniform 0.5C
        reduction_1m = tm_no - tm_1m
        assert 1.5 < reduction_1m < 4.0, \
            f"1M betaine on 100% GC: expected 1.5-4C reduction, got {reduction_1m:.1f}C"

    def test_high_gc_conditions_target_gc_rich(self):
        """High GC condition presets should help with GC-rich genomes."""
        seq_gc_rich = 'GCGCGCGCGCGC'  # 100% GC, 12bp

        standard = get_standard_conditions()
        high_gc = get_high_gc_conditions()

        tm_standard = standard.calculate_effective_tm(seq_gc_rich)
        tm_high_gc = high_gc.calculate_effective_tm(seq_gc_rich)

        # High GC conditions should reduce Tm of GC-rich sequences
        assert tm_high_gc < tm_standard, \
            f"High GC conditions should reduce Tm: standard={tm_standard:.1f}, high_gc={tm_high_gc:.1f}"


class TestATRichSequences:
    """Test GC-dependent corrections on AT-rich sequences."""

    def test_betaine_effect_on_at_rich(self):
        """Betaine should increase (or minimally reduce) Tm of AT-rich sequences."""
        seq_at_rich = 'AAAAAATTTT'  # 0% GC, 10bp
        assert gc_content(seq_at_rich) == 0.0

        conditions_no = ReactionConditions(temp=30, betaine_m=0)
        conditions_1m = ReactionConditions(temp=30, betaine_m=1.0)

        tm_no = conditions_no.calculate_effective_tm(seq_at_rich)
        tm_1m = conditions_1m.calculate_effective_tm(seq_at_rich)

        # AT-rich: GC normalization should increase Tm (counteracting uniform reduction)
        # Net effect depends on relative magnitudes
        # At 1M betaine: uniform -0.5C vs GC correction +1.9C = net +1.4C
        effect = tm_no - tm_1m
        # Effect should be less than for GC-rich (could be positive or negative)
        assert effect < 2.0, f"AT-rich effect should be modest, got {effect:.1f}C"


class TestBalancedSequences:
    """Test GC-dependent corrections on balanced (50% GC) sequences."""

    def test_balanced_minimal_gc_correction(self):
        """50% GC sequences should have minimal GC-dependent correction."""
        seq_balanced = 'ATCGATCGATCG'  # 50% GC, 12bp
        assert abs(gc_content(seq_balanced) - 0.5) < 0.01

        conditions_no = ReactionConditions(temp=30, betaine_m=0)
        conditions_2m = ReactionConditions(temp=30, betaine_m=2.0)

        tm_no = conditions_no.calculate_effective_tm(seq_balanced)
        tm_2m = conditions_2m.calculate_effective_tm(seq_balanced)

        # For 50% GC: uniform effect plus residual GC-dependent correction (~1-2.5C for 2M)
        reduction = tm_no - tm_2m
        assert 0.5 < reduction < 3.0, \
            f"50% GC should show modest effect for 2M betaine, got {reduction:.1f}C"


class TestTmConvergence:
    """Test that high additive concentrations cause Tm convergence."""

    def test_betaine_reduces_tm_range(self):
        """At high betaine, GC-rich and AT-rich sequences should have closer Tm."""
        seq_gc = 'GCGCGCGCGC'  # 100% GC
        seq_at = 'AAAAAATTTT'  # 0% GC
        seq_50 = 'GCGCATATAT'  # 50% GC

        conditions_no = ReactionConditions(temp=30, betaine_m=0)
        conditions_2m = ReactionConditions(temp=30, betaine_m=2.0)

        # Without betaine
        tm_gc_no = conditions_no.calculate_effective_tm(seq_gc)
        tm_at_no = conditions_no.calculate_effective_tm(seq_at)
        tm_50_no = conditions_no.calculate_effective_tm(seq_50)

        # With 2M betaine
        tm_gc_2m = conditions_2m.calculate_effective_tm(seq_gc)
        tm_at_2m = conditions_2m.calculate_effective_tm(seq_at)
        tm_50_2m = conditions_2m.calculate_effective_tm(seq_50)

        # Range without betaine
        range_no = tm_gc_no - tm_at_no

        # Range with betaine (should be smaller)
        range_2m = tm_gc_2m - tm_at_2m

        assert range_2m < range_no, \
            f"Betaine should reduce Tm range: no betaine={range_no:.1f}C, 2M betaine={range_2m:.1f}C"


class TestExtremeGCConditionsPreset:
    """Test the extreme GC conditions preset."""

    def test_extreme_gc_preset_has_tmac_and_betaine(self):
        """Extreme GC preset should have both TMAC and betaine."""
        extreme_gc = get_extreme_gc_conditions()

        assert extreme_gc.tmac_m > 0, "Extreme GC should have TMAC"
        assert extreme_gc.betaine_m >= 2.0, "Extreme GC should have high betaine"
        assert extreme_gc.urea_m > 0, "Extreme GC should have urea"

    def test_extreme_gc_preset_effectiveness(self):
        """Extreme GC preset should effectively handle high-GC sequences."""
        seq_extreme_gc = 'GCGCGCGCGCGCGCGC'  # 100% GC, 16bp

        standard = get_standard_conditions()
        extreme = get_extreme_gc_conditions()

        tm_standard = standard.calculate_effective_tm(seq_extreme_gc)
        tm_extreme = extreme.calculate_effective_tm(seq_extreme_gc)

        # Extreme conditions should substantially reduce Tm
        assert tm_extreme < tm_standard - 5, \
            f"Extreme GC conditions should reduce Tm: standard={tm_standard:.1f}, extreme={tm_extreme:.1f}"


class TestScientificValidity:
    """Test scientific validity against literature values."""

    def test_betaine_isostabilizing_concentration(self):
        """At high betaine, Tm difference between GC-rich and balanced reduces (Rees 1993)."""
        # Note: We can't actually use 5.2M because validation limits to 2.5M
        # Test the concept at maximum allowed concentration
        conditions = ReactionConditions(temp=30, betaine_m=2.5)

        # Compare GC-rich and balanced
        seq_gc = 'GCGCGCGCGC'
        seq_50 = 'GCGCATATAT'

        tm_gc = conditions.calculate_effective_tm(seq_gc)
        tm_50 = conditions.calculate_effective_tm(seq_50)

        # At 2.5M (about half of isostabilizing), expect partial equalization
        conditions_no = ReactionConditions(temp=30, betaine_m=0)
        tm_gc_no = conditions_no.calculate_effective_tm(seq_gc)
        tm_50_no = conditions_no.calculate_effective_tm(seq_50)

        diff_no = abs(tm_gc_no - tm_50_no)
        diff_2m = abs(tm_gc - tm_50)

        # At 2.5M betaine, expect measurable reduction in Tm difference
        # The Wallace rule approximation predicts ~48% reduction at 2.5/5.2 = 48% equalization
        # The actual NN model shows about 20-25% reduction due to different GC contributions
        # Both demonstrate the GC-normalization effect direction
        assert diff_2m < diff_no, \
            f"2.5M betaine should reduce Tm difference: no betaine={diff_no:.1f}C, 2.5M={diff_2m:.1f}C"

        # Expect at least 15% reduction
        reduction_percent = (diff_no - diff_2m) / diff_no * 100
        assert reduction_percent > 15, \
            f"Expected >15% reduction in Tm difference, got {reduction_percent:.1f}%"

    def test_tmac_at_practical_concentrations(self):
        """TMAC effect at practical SWGA concentrations."""
        # TMAC validation range is 0-0.1M (3M is full equalization)
        conditions = ReactionConditions(temp=30, tmac_m=0.1)

        # At 0.1M = 3.3% of full equalization
        # Should see modest but measurable effect
        seq_gc = 'GCGCGCGCGC'

        conditions_no = ReactionConditions(temp=30, tmac_m=0)
        tm_no = conditions_no.calculate_effective_tm(seq_gc)
        tm_01 = conditions.calculate_effective_tm(seq_gc)

        # Expect small reduction (TMAC has uniform effect + small GC normalization)
        reduction = tm_no - tm_01
        assert reduction > 0.5, f"0.1M TMAC should reduce Tm by at least 0.5C, got {reduction:.2f}C"


class TestCombinedEffects:
    """Test combined TMAC and betaine effects."""

    def test_combined_tmac_betaine(self):
        """TMAC and betaine together: use maximum effect (not additive)."""
        seq_gc = 'GCGCGCGCGC'

        # Just betaine
        cond_bet = ReactionConditions(temp=30, betaine_m=1.0, tmac_m=0)

        # Just TMAC (at max practical concentration)
        cond_tmac = ReactionConditions(temp=30, betaine_m=0, tmac_m=0.1)

        # Both
        cond_both = ReactionConditions(temp=30, betaine_m=1.0, tmac_m=0.1)

        # Control
        cond_none = ReactionConditions(temp=30, betaine_m=0, tmac_m=0)

        tm_none = cond_none.calculate_effective_tm(seq_gc)
        tm_bet = cond_bet.calculate_effective_tm(seq_gc)
        tm_tmac = cond_tmac.calculate_effective_tm(seq_gc)
        tm_both = cond_both.calculate_effective_tm(seq_gc)

        # Both should give at least as much reduction as either alone
        reduction_bet = tm_none - tm_bet
        reduction_tmac = tm_none - tm_tmac
        reduction_both = tm_none - tm_both

        # Combined effect uses max of individual GC normalizations (not sum)
        # Plus both uniform effects
        assert reduction_both >= max(reduction_bet, reduction_tmac) - 0.5, \
            f"Combined should not be less than individual: bet={reduction_bet:.1f}, tmac={reduction_tmac:.1f}, both={reduction_both:.1f}"


class TestLengthDependentCorrection:
    """Test length-dependent GC correction scaling."""

    def test_correction_scales_with_length(self):
        """GC correction should scale linearly with primer length."""
        conditions = ReactionConditions(temp=30, betaine_m=2.0)

        # At 80% GC, the GC correction should be larger for longer primers
        gc_content = 0.8

        corr_6bp = conditions._calculate_gc_normalization(gc_content, primer_length=6)
        corr_10bp = conditions._calculate_gc_normalization(gc_content, primer_length=10)
        corr_18bp = conditions._calculate_gc_normalization(gc_content, primer_length=18)

        # Corrections should be negative (GC-rich gets reduced Tm)
        assert corr_6bp < 0
        assert corr_10bp < 0
        assert corr_18bp < 0

        # Longer primers should have larger (more negative) corrections
        assert abs(corr_18bp) > abs(corr_10bp) > abs(corr_6bp)

        # Check scaling: correction should scale with 2*length
        # 6bp: scale = 12, 10bp: scale = 20, 18bp: scale = 36
        ratio_10_6 = corr_10bp / corr_6bp
        ratio_18_10 = corr_18bp / corr_10bp

        expected_ratio_10_6 = 20.0 / 12.0  # 1.67
        expected_ratio_18_10 = 36.0 / 20.0  # 1.80

        assert abs(ratio_10_6 - expected_ratio_10_6) < 0.01
        assert abs(ratio_18_10 - expected_ratio_18_10) < 0.01

    def test_effective_tm_uses_actual_length(self):
        """calculate_effective_tm should use actual primer length."""
        conditions = ReactionConditions(temp=30, betaine_m=2.0)

        # Same GC content, different lengths
        seq_8bp = 'GCGCGCGC'   # 8bp, 100% GC
        seq_16bp = 'GCGCGCGCGCGCGCGC'  # 16bp, 100% GC

        # Both have same GC content, but 16bp should get larger correction
        tm_8bp = conditions.calculate_effective_tm(seq_8bp)
        tm_16bp = conditions.calculate_effective_tm(seq_16bp)

        # Base Tm without correction: longer sequence has higher Tm
        base_tm_8bp = calculate_tm_with_salt(seq_8bp, 50, 0)
        base_tm_16bp = calculate_tm_with_salt(seq_16bp, 50, 0)

        # The Tm difference should be reduced by length-dependent correction
        # because longer primer gets larger GC correction
        base_diff = base_tm_16bp - base_tm_8bp
        effective_diff = tm_16bp - tm_8bp

        # With length-dependent correction, the effective difference should be less
        # because 16bp gets ~2x the GC correction of 8bp
        assert effective_diff < base_diff

    def test_short_primer_correction(self):
        """Short primers (6bp) should get smaller GC correction."""
        conditions = ReactionConditions(temp=30, betaine_m=1.0)

        # 6bp primer at 100% GC
        seq_6bp = 'GCGCGC'
        gc = gc_content(seq_6bp)
        assert gc == 1.0

        # Calculate correction for 6bp (scale = 12)
        # gc_deviation = 0.5, betaine_factor = 1/5.2 ≈ 0.192
        # correction = -12 * 0.5 * 0.192 ≈ -1.15
        corr = conditions._calculate_gc_normalization(gc, primer_length=6)
        assert -2.0 < corr < -0.5

    def test_long_primer_correction(self):
        """Long primers (18bp) should get larger GC correction."""
        conditions = ReactionConditions(temp=30, betaine_m=1.0)

        # 18bp primer at 100% GC
        gc = 1.0

        # Calculate correction for 18bp (scale = 36)
        # gc_deviation = 0.5, betaine_factor = 1/5.2 ≈ 0.192
        # correction = -36 * 0.5 * 0.192 ≈ -3.46
        corr = conditions._calculate_gc_normalization(gc, primer_length=18)
        assert -5.0 < corr < -2.0


class TestEdgeCases:
    """Test edge cases for GC-dependent corrections."""

    def test_zero_concentrations(self):
        """Zero additive concentrations should give no GC normalization."""
        conditions = ReactionConditions(temp=30, betaine_m=0, tmac_m=0)

        for gc in [0.0, 0.25, 0.5, 0.75, 1.0]:
            correction = conditions._calculate_gc_normalization(gc)
            assert correction == 0, f"Zero additives should give zero correction at {gc} GC"

    def test_extreme_gc_values(self):
        """Handle extreme GC values (0% and 100%)."""
        conditions = ReactionConditions(temp=30, betaine_m=2.0)

        # 0% GC
        corr_0 = conditions._calculate_gc_normalization(0.0)
        assert corr_0 > 0, "0% GC should get positive correction"

        # 100% GC
        corr_100 = conditions._calculate_gc_normalization(1.0)
        assert corr_100 < 0, "100% GC should get negative correction"

        # Symmetry
        assert abs(corr_0 + corr_100) < 0.01, "0% and 100% GC should be symmetric"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
