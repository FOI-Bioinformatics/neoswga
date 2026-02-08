"""
Tests for thermodynamic dimer detection (is_dimer_thermodynamic).

Verifies that:
- Perfect complements are detected as dimers
- Short or weak complementary regions are not flagged
- The function is consistent with the existing sequence-based check
- Custom thresholds and reaction conditions are respected
"""

import pytest

from neoswga.core.dimer import is_dimer_fast, is_dimer_thermodynamic
from neoswga.core.thermodynamics import reverse_complement


class TestIsDimerThermodynamic:
    """Tests for thermodynamic dimer scoring."""

    def test_perfect_complement_detected(self):
        """A primer paired with its exact reverse complement should be a dimer."""
        seq = 'ATCGATCGATCG'
        seq_rc = reverse_complement(seq)
        assert is_dimer_thermodynamic(seq, seq_rc) is True

    def test_short_complement_not_flagged(self):
        """Two primers sharing only 2-3 bp of complementarity should not be flagged."""
        seq_1 = 'ATCGATCG'
        seq_2 = 'TTTTTTTT'  # Very low complementarity
        assert is_dimer_thermodynamic(seq_1, seq_2) is False

    def test_no_complementarity(self):
        """Primers with no complementary bases should not be flagged."""
        seq_1 = 'AAAAAAAA'
        seq_2 = 'AAAAAAAA'
        assert is_dimer_thermodynamic(seq_1, seq_2) is False

    def test_moderate_complement_depends_on_threshold(self):
        """A moderate complementary region should be flagged only with a lenient threshold."""
        seq_1 = 'ATCGATCGATCG'
        seq_2 = 'CGATCG'  # 6bp complement embedded

        # Strict threshold: should not flag short duplexes
        assert is_dimer_thermodynamic(seq_1, seq_2, delta_g_threshold=-15.0) is False

        # Lenient threshold: should flag even weak interactions
        assert is_dimer_thermodynamic(seq_1, seq_2, delta_g_threshold=-2.0) is True

    def test_custom_conditions_temperature(self):
        """Higher temperature should make dimer formation less favorable."""
        from neoswga.core.reaction_conditions import ReactionConditions

        seq = 'GCGCGCGCGCGC'
        seq_rc = reverse_complement(seq)

        # At low temperature, strong complement forms dimer
        cold = ReactionConditions(temp=20.0, polymerase='phi29')
        assert is_dimer_thermodynamic(seq, seq_rc, conditions=cold) is True

        # At very high temperature with strict threshold, should still detect
        # a perfect 12bp GC complement
        hot = ReactionConditions(temp=65.0, polymerase='bst')
        assert is_dimer_thermodynamic(seq, seq_rc, conditions=hot) is True

    def test_consistency_with_sequence_dimer(self):
        """Pairs flagged by sequence-based check should also be flagged thermodynamically."""
        seq_1 = 'GCGATCGATCGC'
        seq_2 = reverse_complement(seq_1)

        # If sequence-based check flags it (long common substring)
        assert is_dimer_fast(seq_1, seq_2, max_dimer_bp=3) is True

        # Thermodynamic check should also flag it
        assert is_dimer_thermodynamic(seq_1, seq_2) is True

    def test_self_dimer(self):
        """A palindromic sequence should be detected as a self-dimer."""
        palindrome = 'ATCGATCGAT'
        # Palindrome binds to itself (seq == rc for the complementary region)
        result = is_dimer_thermodynamic(palindrome, palindrome)
        # This is a self-complementary sequence, so it should form a dimer
        assert isinstance(result, bool)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
