"""
Unit tests for filter module.

Tests:
- filter_extra function (primer design rules)
- GC content filtering
- Homopolymer detection
- Dinucleotide repeat detection
- GC clamp rules
"""

import pytest
from unittest.mock import patch, MagicMock

from neoswga.core.filter import filter_extra


# =============================================================================
# Test Fixtures
# =============================================================================

@pytest.fixture(autouse=True)
def reset_reaction_conditions():
    """Reset cached reaction conditions between tests."""
    from neoswga.core.filter import reset_reaction_conditions
    reset_reaction_conditions()
    yield
    reset_reaction_conditions()


@pytest.fixture
def mock_parameter():
    """Mock the parameter module with standard settings."""
    with patch('neoswga.core.filter.parameter') as mock_param:
        mock_param.min_tm = 15
        mock_param.max_tm = 55
        mock_param.gc_min = 0.375
        mock_param.gc_max = 0.625
        mock_param.verbose = False
        mock_param.max_self_dimer_bp = 4
        mock_param.genome_gc = None
        mock_param.polymerase = 'phi29'
        mock_param.reaction_temp = 30.0
        mock_param.na_conc = 50.0
        mock_param.mg_conc = 0.0
        mock_param.dmso_percent = 0.0
        mock_param.betaine_m = 0.0
        mock_param.trehalose_m = 0.0
        mock_param.formamide_percent = 0.0
        mock_param.ethanol_percent = 0.0
        mock_param.urea_m = 0.0
        mock_param.tmac_m = 0.0
        yield mock_param


@pytest.fixture
def mock_dimer():
    """Mock the dimer module."""
    with patch('neoswga.core.filter.dimer') as mock_d:
        mock_d.is_dimer.return_value = False
        mock_d.is_dimer_fast.return_value = False
        yield mock_d


# =============================================================================
# Homopolymer (Rule 5) Tests
# =============================================================================

class TestHomopolymerFilter:
    """Tests for homopolymer detection (Rule 5)."""

    def test_rejects_five_g(self, mock_parameter, mock_dimer):
        """Test rejection of 5 consecutive Gs."""
        result = filter_extra("ATCGGGGGAT")

        assert result is False

    def test_rejects_five_c(self, mock_parameter, mock_dimer):
        """Test rejection of 5 consecutive Cs."""
        result = filter_extra("ATCCCCCGAT")

        assert result is False

    def test_rejects_five_a(self, mock_parameter, mock_dimer):
        """Test rejection of 5 consecutive As."""
        result = filter_extra("GCAAAAAATC")

        assert result is False

    def test_rejects_five_t(self, mock_parameter, mock_dimer):
        """Test rejection of 5 consecutive Ts."""
        result = filter_extra("GCTTTTTGAT")

        assert result is False

    def test_accepts_four_consecutive(self, mock_parameter, mock_dimer):
        """Test acceptance of 4 consecutive same bases.

        ATCGGGGAT has 4 Gs (below the 5-base homopolymer threshold),
        55.6% GC (within 0.375-0.625), valid GC clamp, and no other
        rule violations.
        """
        result = filter_extra("ATCGGGGAT")

        assert result is True


# =============================================================================
# GC Content (Rule 2) Tests
# =============================================================================

class TestGCContentFilter:
    """Tests for GC content filtering (Rule 2)."""

    def test_accepts_50_percent_gc(self, mock_parameter, mock_dimer):
        """Test acceptance of 50% GC content.

        ATCGATCG is 50% GC, within the default 0.375-0.625 range,
        and passes all other filtering rules. Tm range is relaxed
        because this short primer has an effective Tm below the
        default min_tm of 15.
        """
        mock_parameter.min_tm = 0
        mock_parameter.max_tm = 100

        result = filter_extra("ATCGATCG")

        assert result is True

    def test_rejects_all_a(self, mock_parameter, mock_dimer):
        """Test rejection of 0% GC content."""
        result = filter_extra("ATATATAT")

        assert result is False

    def test_rejects_all_gc(self, mock_parameter, mock_dimer):
        """Test rejection of 100% GC content."""
        result = filter_extra("GCGCGCGC")

        assert result is False

    def test_gc_boundaries(self, mock_parameter, mock_dimer):
        """Test GC content at boundaries.

        ATCGATCGAT is 40% GC (4 of 10 bases), within the 0.35-0.65
        range. Last 5 bases (TCGAT) have 2 GC, last 3 (GAT) have 1 GC,
        and no dinucleotide repeats are present.
        """
        mock_parameter.gc_min = 0.35
        mock_parameter.gc_max = 0.65

        result = filter_extra("ATCGATCGAT")

        assert result is True


# =============================================================================
# GC Clamp (Rule 3) Tests
# =============================================================================

class TestGCClampFilter:
    """Tests for GC clamp rules (Rule 3)."""

    def test_rejects_no_gc_in_last_five(self, mock_parameter, mock_dimer):
        """Test rejection when no G/C in last 5 bases."""
        # Last 5 bases are AAAAT - no G or C
        result = filter_extra("GCGCAAAAT")

        assert result is False

    def test_rejects_too_many_gc_in_last_five(self, mock_parameter, mock_dimer):
        """Test rejection when >3 G/C in last 5 bases."""
        # Last 5 bases are GCGCG - 5 G/C
        result = filter_extra("ATATGCGCG")

        assert result is False

    def test_rejects_three_gc_in_last_three(self, mock_parameter, mock_dimer):
        """Test rejection when last 3 bases are all G/C.

        ATCGATGCG passes the GC clamp rule (last 5 = ATGCG, 3 GC),
        but fails Rule 1 because all 3 bases at the 3' end (GCG) are
        G or C, exceeding the MAX_GC_AT_3PRIME_END limit of 2.
        """
        mock_parameter.gc_min = 0.3
        mock_parameter.gc_max = 0.7

        primer = "ATCGATGCG"
        result = filter_extra(primer)

        assert result is False


# =============================================================================
# 3' End GC (Rule 1) Tests
# =============================================================================

class TestThreePrimeGCFilter:
    """Tests for 3' end GC rules (Rule 1)."""

    def test_rejects_three_gc_at_end(self, mock_parameter, mock_dimer):
        """Test rejection when last 3 bases are all G/C."""
        # Last 3 bases are GCG - all G/C
        result = filter_extra("ATATATGCG")

        assert result is False

    def test_rejects_three_c_at_end(self, mock_parameter, mock_dimer):
        """Test rejection when last 3 bases are CCC."""
        result = filter_extra("ATATATCCC")

        assert result is False

    def test_accepts_two_gc_at_end(self, mock_parameter, mock_dimer):
        """Test acceptance of 2 G/C in last 3 bases.

        ATCGATAGC has last 3 = AGC (2 GC, within limit), last 5 = ATAGC
        (2 GC, within 1-3 range), and 44.4% GC overall.
        """
        mock_parameter.gc_min = 0.3
        mock_parameter.gc_max = 0.7

        primer = "ATCGATAGC"
        result = filter_extra(primer)

        assert result is True


# =============================================================================
# Dinucleotide Repeat (Rule 4) Tests
# =============================================================================

class TestDinucleotideRepeatFilter:
    """Tests for dinucleotide repeat detection (Rule 4)."""

    def test_rejects_at_repeat(self, mock_parameter, mock_dimer):
        """Test rejection of 5 consecutive AT repeats."""
        result = filter_extra("ATATATATAT")

        assert result is False

    def test_rejects_ta_repeat(self, mock_parameter, mock_dimer):
        """Test rejection of 5 consecutive TA repeats."""
        result = filter_extra("TATATATATA")

        assert result is False

    def test_rejects_gc_repeat(self, mock_parameter, mock_dimer):
        """Test rejection of 5 consecutive GC repeats."""
        result = filter_extra("GCGCGCGCGC")

        assert result is False

    def test_rejects_ag_repeat(self, mock_parameter, mock_dimer):
        """Test rejection of 5 consecutive AG repeats."""
        result = filter_extra("AGAGAGAGAG")

        assert result is False

    def test_rejects_low_gc_with_short_repeat(self, mock_parameter, mock_dimer):
        """Test that a short repeat primer is rejected for low GC content.

        GCATATAT has only 2/8 = 25% GC, which falls below gc_min=0.3.
        Although the dinucleotide repeat (ATATAT) is short enough to pass
        Rule 4 (length < 10), the primer fails the GC content filter.
        """
        mock_parameter.gc_min = 0.3
        mock_parameter.gc_max = 0.7

        primer = "GCATATAT"
        result = filter_extra(primer)

        assert result is False


# =============================================================================
# Self-Dimer Tests
# =============================================================================

class TestSelfDimerFilter:
    """Tests for self-dimer detection."""

    def test_rejects_self_dimer(self, mock_parameter):
        """Test rejection of primer that forms self-dimer."""
        with patch('neoswga.core.filter.dimer') as mock_dimer:
            mock_dimer.is_dimer_fast.return_value = True

            result = filter_extra("ATCGATCG")

            assert result is False

    def test_accepts_no_self_dimer(self, mock_parameter):
        """Test acceptance when no self-dimer.

        ATCGATCG (50% GC) passes all filtering rules when the dimer
        check returns False, and GC and Tm ranges are relaxed.
        """
        with patch('neoswga.core.filter.dimer') as mock_dimer:
            mock_dimer.is_dimer_fast.return_value = False
            mock_parameter.gc_min = 0.3
            mock_parameter.gc_max = 0.7
            mock_parameter.min_tm = 0
            mock_parameter.max_tm = 100

            result = filter_extra("ATCGATCG")

            assert result is True


# =============================================================================
# Tm Filter Tests
# =============================================================================

class TestTmFilter:
    """Tests for melting temperature filtering."""

    def test_rejects_low_tm(self, mock_parameter, mock_dimer):
        """Test rejection of primer with low Tm."""
        mock_parameter.min_tm = 50
        mock_parameter.max_tm = 70

        # Short primer will have low Tm
        result = filter_extra("ATCGATCG")

        assert result is False

    def test_rejects_high_tm(self, mock_parameter, mock_dimer):
        """Test rejection of primer with high Tm.

        GCGCGCGCGCGCGCGC (100% GC, 16 bp) has a Tm far above the
        max_tm=10 threshold and is rejected. It would also fail the
        GC content filter, but the Tm check runs first.
        """
        mock_parameter.min_tm = 0
        mock_parameter.max_tm = 10

        result = filter_extra("GCGCGCGCGCGCGCGC")

        assert result is False


# =============================================================================
# Integration Tests
# =============================================================================

class TestFilterIntegration:
    """Integration tests for filter_extra."""

    def test_good_primer_passes(self, mock_parameter, mock_dimer):
        """Test that a well-designed primer passes all filters."""
        mock_parameter.gc_min = 0.3
        mock_parameter.gc_max = 0.7
        mock_parameter.min_tm = 10
        mock_parameter.max_tm = 60

        # Design a primer that should pass all rules:
        # - No homopolymers (5+ same base)
        # - 40-60% GC
        # - 1-3 G/C in last 5 bases
        # - Not all G/C in last 3 bases
        # - No dinucleotide repeats
        good_primer = "ATCGTAGC"

        result = filter_extra(good_primer)

        # Should pass all rules
        assert result is True

    def test_multiple_violations(self, mock_parameter, mock_dimer):
        """Test primer with multiple rule violations."""
        # AAAAAAAAAAGGG violates:
        # - Homopolymer rule (6 As)
        # - GC content rule (too low)
        # - Last 3 bases all G/C
        result = filter_extra("AAAAAAAAGGG")

        assert result is False


# =============================================================================
# Edge Cases
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases."""

    def test_minimum_length_primer(self, mock_parameter, mock_dimer):
        """Test with minimum length primer.

        ATCG (4 bp, 50% GC) passes all rules when GC and Tm ranges
        are fully relaxed. Last 3 = TCG (2 GC), last 5 = ATCG (2 GC).
        """
        mock_parameter.gc_min = 0.0
        mock_parameter.gc_max = 1.0
        mock_parameter.min_tm = -100
        mock_parameter.max_tm = 100

        result = filter_extra("ATCG")

        assert result is True

    def test_exact_boundary_gc(self, mock_parameter, mock_dimer):
        """Test primer at exact GC boundary.

        Boundaries are inclusive: primers at exactly gc_min or gc_max are accepted.
        """
        # Widen Tm range so short primer Tm doesn't interfere with this GC test
        mock_parameter.min_tm = 0
        mock_parameter.max_tm = 100
        mock_parameter.gc_min = 0.5
        mock_parameter.gc_max = 0.5

        # Exactly 50% GC should be accepted (inclusive boundaries)
        result = filter_extra("ATCGATCG")

        assert result is True  # Inclusive boundary

        # But just below gc_min should be rejected
        mock_parameter.gc_min = 0.6
        result_below = filter_extra("ATCGATCG")  # 50% GC
        assert result_below is False  # Below gc_min

        # And just above gc_max should be rejected
        mock_parameter.gc_min = 0.4
        mock_parameter.gc_max = 0.49
        result_above = filter_extra("ATCGATCG")  # 50% GC
        assert result_above is False  # Above gc_max

    def test_short_primer_skips_dinucleotide_check(self, mock_parameter, mock_dimer):
        """Test that short primers skip dinucleotide repeat check.

        ATATATAT (8 bp) would fail the dinucleotide repeat check if it
        were applied, but Rule 4 only checks primers with length >= 10.
        However, this primer has 0% GC (0/8), which falls below gc_min=0.3,
        so it is rejected by the GC content filter instead.
        """
        mock_parameter.gc_min = 0.3
        mock_parameter.gc_max = 0.7
        mock_parameter.min_tm = 0
        mock_parameter.max_tm = 100

        primer = "ATATATAT"

        result = filter_extra(primer)

        assert result is False


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
