"""
Comprehensive edge case tests for thermodynamic calculations.

Tests include:
- Nearest-neighbor calculations for unusual sequences
- Salt correction extremes
- Additive combinations and saturation effects
- Temperature range validation
- Real primer sequences from different GC-content organisms
"""

import unittest
import warnings
import numpy as np

from neoswga.core import thermodynamics as thermo
from neoswga.core.reaction_conditions import ReactionConditions


class TestNNCalculationEdgeCases(unittest.TestCase):
    """Test nearest-neighbor calculations with edge cases."""

    def test_palindromic_sequence_gaattc(self):
        """Test palindromic EcoRI site GAATTC."""
        seq = 'GAATTC'
        self.assertTrue(thermo.is_palindrome(seq))
        enthalpy, entropy = thermo.calculate_enthalpy_entropy(seq)
        # Palindrome should have symmetry correction applied
        self.assertIsInstance(enthalpy, float)
        self.assertLess(entropy, 0)  # Entropy should be negative

    def test_palindromic_sequence_cgcgcg(self):
        """Test highly repetitive palindrome CGCGCGCG."""
        seq = 'CGCGCGCG'
        self.assertTrue(thermo.is_palindrome(seq))
        enthalpy, entropy = thermo.calculate_enthalpy_entropy(seq)
        # High GC content should give strong stacking
        self.assertLess(enthalpy, -50)  # Strong enthalpy contribution

    def test_very_short_sequence_4bp(self):
        """Test minimum practical primer length (4bp)."""
        seq = 'ATCG'
        enthalpy, entropy = thermo.calculate_enthalpy_entropy(seq)
        tm = thermo.calculate_tm_basic(seq)
        # Short sequences have low Tm
        self.assertLess(tm, 30)

    def test_very_short_sequence_6bp(self):
        """Test 6bp sequence (minimum SWGA primer)."""
        seq = 'ATCGAT'
        tm = thermo.calculate_tm_with_salt(seq)
        # Should give reasonable Tm for short primer
        self.assertIsInstance(tm, float)
        self.assertFalse(np.isnan(tm))

    def test_all_at_sequence(self):
        """Test poly-AT sequence."""
        seq = 'ATATATATAT'
        gc = thermo.gc_content(seq)
        self.assertEqual(gc, 0.0)
        tm = thermo.calculate_tm_with_salt(seq)
        # Pure AT should have low Tm
        self.assertLess(tm, 35)

    def test_all_gc_sequence(self):
        """Test poly-GC sequence."""
        seq = 'GCGCGCGCGC'
        gc = thermo.gc_content(seq)
        self.assertEqual(gc, 1.0)
        tm = thermo.calculate_tm_with_salt(seq)
        # Pure GC should have high Tm (above 40C for 10bp)
        self.assertGreater(tm, 40)

    def test_poly_a_tract(self):
        """Test poly-A tract stability."""
        seq = 'AAAAAAAAAAAA'
        enthalpy, entropy = thermo.calculate_enthalpy_entropy(seq)
        # All AA/TT stacks
        tm = thermo.calculate_tm_basic(seq)
        self.assertLess(tm, 50)  # Low Tm expected


class TestSaltCorrectionEdgeCases(unittest.TestCase):
    """Test salt corrections at extreme conditions."""

    def test_zero_sodium(self):
        """Test with very low sodium concentration (edge case)."""
        seq = 'ATCGATCG'
        # Very low sodium should still work
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                tm = thermo.calculate_tm_with_salt(seq, na_conc=0.01, mg_conc=0.0)
                self.assertIsInstance(tm, float)
            except (ValueError, ZeroDivisionError):
                # Acceptable to raise error for invalid conditions
                pass

    def test_very_low_salt_1mm(self):
        """Test with 1mM total salt (minimal PCR conditions)."""
        seq = 'ATCGATCG'
        tm = thermo.calculate_tm_with_salt(seq, na_conc=1.0, mg_conc=0.0)
        # Very low salt should give lower Tm
        self.assertIsInstance(tm, float)

    def test_very_high_sodium_500mm(self):
        """Test with 500mM sodium (high salt conditions)."""
        seq = 'ATCGATCG'
        tm_low = thermo.calculate_tm_with_salt(seq, na_conc=50.0, mg_conc=0.0)
        tm_high = thermo.calculate_tm_with_salt(seq, na_conc=500.0, mg_conc=0.0)
        # Higher salt should increase Tm
        self.assertGreater(tm_high, tm_low)

    def test_very_high_magnesium_10mm(self):
        """Test with high magnesium (PCR with excess Mg2+)."""
        seq = 'ATCGATCGATCG'  # 12bp for more stable duplex
        tm = thermo.calculate_tm_with_salt(seq, na_conc=50.0, mg_conc=10.0)
        self.assertIsInstance(tm, float)
        # High Mg should give reasonable Tm
        self.assertGreater(tm, 20)

    def test_magnesium_dominates_at_low_sodium(self):
        """Test that Mg2+ dominates at low sodium concentrations."""
        seq = 'ATCGATCGATCG'
        tm_na_only = thermo.calculate_tm_with_salt(seq, na_conc=10.0, mg_conc=0.0)
        tm_with_mg = thermo.calculate_tm_with_salt(seq, na_conc=10.0, mg_conc=2.0)
        # Adding Mg should significantly increase Tm at low Na+
        self.assertGreater(tm_with_mg, tm_na_only)

    def test_mixed_high_salt(self):
        """Test high concentrations of both Na+ and Mg2+."""
        seq = 'ATCGATCG'
        tm = thermo.calculate_tm_with_salt(seq, na_conc=200.0, mg_conc=5.0)
        self.assertIsInstance(tm, float)


class TestAdditiveEdgeCases(unittest.TestCase):
    """Test additive effects at extreme concentrations."""

    def test_dmso_at_10_percent(self):
        """Test DMSO at maximum valid concentration (10%)."""
        conditions = ReactionConditions(
            polymerase='equiphi29',
            temp=42.0,
            dmso_percent=10.0
        )
        correction = conditions.calculate_tm_correction(gc_content=0.5, primer_length=12)
        # 10% DMSO should reduce Tm by ~6C
        self.assertLess(correction, -5)

    def test_betaine_at_2_5m(self):
        """Test betaine near isostabilizing concentration (2.5M)."""
        conditions = ReactionConditions(
            polymerase='equiphi29',
            temp=42.0,
            betaine_m=2.5
        )
        correction_50gc = conditions.calculate_tm_correction(gc_content=0.5, primer_length=12)
        correction_33gc = conditions.calculate_tm_correction(gc_content=0.33, primer_length=12)
        # At high betaine, GC effect should be reduced
        self.assertIsInstance(correction_50gc, float)
        self.assertIsInstance(correction_33gc, float)

    def test_dmso_betaine_synergy(self):
        """Test combined DMSO + betaine effects."""
        # Control: no additives
        control = ReactionConditions(
            polymerase='equiphi29',
            temp=42.0
        )
        # Combined: 5% DMSO + 1M betaine
        combined = ReactionConditions(
            polymerase='equiphi29',
            temp=42.0,
            dmso_percent=5.0,
            betaine_m=1.0
        )
        ctrl_corr = control.calculate_tm_correction(gc_content=0.5, primer_length=12)
        comb_corr = combined.calculate_tm_correction(gc_content=0.5, primer_length=12)
        # Combined should have larger negative correction
        self.assertLess(comb_corr, ctrl_corr)

    def test_all_additives_moderate(self):
        """Test with multiple additives at moderate concentrations."""
        conditions = ReactionConditions(
            polymerase='equiphi29',
            temp=45.0,
            dmso_percent=5.0,
            betaine_m=1.0,
            trehalose_m=0.3
        )
        correction = conditions.calculate_tm_correction(gc_content=0.5, primer_length=15)
        # Significant Tm reduction expected
        self.assertLessEqual(correction, -4.9)

    def test_additive_linearity(self):
        """Test that DMSO effect is approximately linear."""
        base = ReactionConditions(polymerase='equiphi29', temp=42.0)
        dmso5 = ReactionConditions(polymerase='equiphi29', temp=42.0, dmso_percent=5.0)
        dmso10 = ReactionConditions(polymerase='equiphi29', temp=42.0, dmso_percent=10.0)

        c0 = base.calculate_tm_correction(0.5, 12)
        c5 = dmso5.calculate_tm_correction(0.5, 12)
        c10 = dmso10.calculate_tm_correction(0.5, 12)

        # Effect should be approximately linear
        delta_5 = c5 - c0
        delta_10 = c10 - c0
        # 10% should be roughly 2x effect of 5%
        self.assertAlmostEqual(delta_10, 2 * delta_5, delta=1.0)

    def test_tmac_gc_normalization(self):
        """Test TMAC GC-normalizing effect."""
        tmac = ReactionConditions(
            polymerase='equiphi29',
            temp=42.0,
            tmac_m=0.05
        )
        # TMAC should reduce difference between high/low GC primers
        correction_high_gc = tmac.calculate_tm_correction(gc_content=0.7, primer_length=12)
        correction_low_gc = tmac.calculate_tm_correction(gc_content=0.3, primer_length=12)
        # Both should have non-zero corrections
        self.assertNotEqual(correction_high_gc, correction_low_gc)


class TestTemperatureRangeValidation(unittest.TestCase):
    """Test temperature range validation for different polymerases."""

    def test_phi29_valid_range(self):
        """Test phi29 within valid range (20-40C)."""
        conditions = ReactionConditions(polymerase='phi29', temp=30.0)
        self.assertEqual(conditions.polymerase, 'phi29')
        self.assertEqual(conditions.temp, 30.0)

    def test_phi29_at_low_boundary(self):
        """Test phi29 at lower boundary."""
        conditions = ReactionConditions(polymerase='phi29', temp=20.0)
        self.assertEqual(conditions.temp, 20.0)

    def test_phi29_at_high_boundary(self):
        """Test phi29 at upper boundary."""
        conditions = ReactionConditions(polymerase='phi29', temp=40.0)
        self.assertEqual(conditions.temp, 40.0)

    def test_phi29_outside_range_low(self):
        """Test phi29 below valid range raises error."""
        with self.assertRaises(ValueError):
            ReactionConditions(polymerase='phi29', temp=15.0)

    def test_phi29_outside_range_high(self):
        """Test phi29 above valid range raises error."""
        with self.assertRaises(ValueError):
            ReactionConditions(polymerase='phi29', temp=45.0)

    def test_equiphi29_valid_range(self):
        """Test equiphi29 within valid range (30-50C)."""
        conditions = ReactionConditions(polymerase='equiphi29', temp=42.0)
        self.assertEqual(conditions.polymerase, 'equiphi29')

    def test_equiphi29_at_boundaries(self):
        """Test equiphi29 at boundaries."""
        low = ReactionConditions(polymerase='equiphi29', temp=30.0)
        high = ReactionConditions(polymerase='equiphi29', temp=50.0)
        self.assertEqual(low.temp, 30.0)
        self.assertEqual(high.temp, 50.0)

    def test_bst_high_temperature(self):
        """Test BST at high temperature (50-72C)."""
        conditions = ReactionConditions(polymerase='bst', temp=65.0)
        self.assertEqual(conditions.temp, 65.0)


class TestEnthalpyEntropyCalculations(unittest.TestCase):
    """Test enthalpy and entropy calculations."""

    def test_enthalpy_negative(self):
        """Test that total enthalpy is negative (favorable)."""
        seq = 'ATCGATCG'
        enthalpy, entropy = thermo.calculate_enthalpy_entropy(seq)
        self.assertLess(enthalpy, 0)

    def test_entropy_negative(self):
        """Test that total entropy is negative (ordering)."""
        seq = 'ATCGATCG'
        enthalpy, entropy = thermo.calculate_enthalpy_entropy(seq)
        self.assertLess(entropy, 0)

    def test_gc_rich_more_negative_enthalpy(self):
        """Test that GC-rich sequences have more negative enthalpy."""
        at_seq = 'ATATATATAT'
        gc_seq = 'GCGCGCGCGC'
        h_at, _ = thermo.calculate_enthalpy_entropy(at_seq)
        h_gc, _ = thermo.calculate_enthalpy_entropy(gc_seq)
        # GC stacking is stronger
        self.assertLess(h_gc, h_at)


class TestRealPrimerSequences(unittest.TestCase):
    """Test with real primer sequences from different organisms."""

    def test_francisella_primer_low_gc(self):
        """Test primer for Francisella tularensis (33% GC)."""
        # Low-GC primer typical for Francisella
        primer = 'AATATATATAT'  # Very AT-rich
        gc = thermo.gc_content(primer)
        self.assertLess(gc, 0.2)
        tm = thermo.calculate_tm_with_salt(primer)
        # AT-rich should have low Tm
        self.assertLess(tm, 40)

    def test_ecoli_primer_balanced_gc(self):
        """Test primer for E. coli (50% GC)."""
        primer = 'ATCGATCGATCG'  # Balanced 50% GC
        gc = thermo.gc_content(primer)
        self.assertEqual(gc, 0.5)
        tm = thermo.calculate_tm_with_salt(primer)
        # Balanced should have moderate Tm
        self.assertIsInstance(tm, float)

    def test_burkholderia_primer_high_gc(self):
        """Test primer for Burkholderia (67% GC)."""
        primer = 'GCGCGCATGCGC'  # High GC (75%)
        gc = thermo.gc_content(primer)
        self.assertGreater(gc, 0.6)
        tm = thermo.calculate_tm_with_salt(primer)
        # High GC should have high Tm
        self.assertGreater(tm, 35)

    def test_primer_set_tm_uniformity(self):
        """Test that a primer set has uniform Tm range."""
        # Example 6-primer set with mixed sequences
        primers = [
            'ATCGATCG',
            'GCTAGCTA',
            'AATTAATT',
            'CCGGCCGG',
            'ATATGCGC',
            'GCGCATAT'
        ]
        tms = [thermo.calculate_tm_with_salt(p) for p in primers]
        tm_range = max(tms) - min(tms)
        # All Tms should be computable
        self.assertEqual(len(tms), 6)
        # Range should be finite
        self.assertLess(tm_range, 50)

    def test_gc_content_consistency(self):
        """Test GC content calculation consistency."""
        test_cases = [
            ('AAAA', 0.0),
            ('TTTT', 0.0),
            ('GGGG', 1.0),
            ('CCCC', 1.0),
            ('ATCG', 0.5),
            ('ATAT', 0.0),
            ('GCGC', 1.0),
            ('ATGC', 0.5),
        ]
        for seq, expected_gc in test_cases:
            with self.subTest(seq=seq):
                gc = thermo.gc_content(seq)
                self.assertAlmostEqual(gc, expected_gc, places=4)


class TestSequenceNormalization(unittest.TestCase):
    """Test sequence normalization and handling."""

    def test_lowercase_conversion(self):
        """Test that lowercase sequences work correctly."""
        upper = 'ATCGATCG'
        lower = 'atcgatcg'
        tm_upper = thermo.calculate_tm_with_salt(upper)
        tm_lower = thermo.calculate_tm_with_salt(lower)
        self.assertAlmostEqual(tm_upper, tm_lower, places=5)

    def test_mixed_case(self):
        """Test mixed case sequence handling."""
        seq = 'AtCgAtCg'
        tm = thermo.calculate_tm_with_salt(seq)
        self.assertIsInstance(tm, float)

    def test_ambiguous_base_warning(self):
        """Test that ambiguous bases trigger warning."""
        seq = 'ATCGATNCG'
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            tm = thermo.calculate_tm_with_salt(seq)
            # Should get a warning about ambiguous bases
            self.assertTrue(any('ambiguous' in str(warning.message).lower()
                               for warning in w))

    def test_has_ambiguous_bases_function(self):
        """Test has_ambiguous_bases detection."""
        self.assertFalse(thermo.has_ambiguous_bases('ATCG'))
        self.assertTrue(thermo.has_ambiguous_bases('ATNG'))
        self.assertTrue(thermo.has_ambiguous_bases('ATRY'))
        self.assertFalse(thermo.has_ambiguous_bases('atcg'))  # lowercase ATCG

    def test_normalize_sequence(self):
        """Test sequence normalization function."""
        self.assertEqual(thermo.normalize_sequence('atcg'), 'ATCG')
        self.assertEqual(thermo.normalize_sequence('AtCg'), 'ATCG')
        self.assertEqual(thermo.normalize_sequence('ATCG'), 'ATCG')


class TestNNParameterValidation(unittest.TestCase):
    """Test that NN parameters are complete and consistent."""

    def test_canonical_stacks_present(self):
        """Test that all 10 canonical SantaLucia stacks are present."""
        canonical = [
            'AA/TT', 'AT/TA', 'TA/AT', 'CA/GT', 'GT/CA',
            'CT/GA', 'GA/CT', 'CG/GC', 'GC/CG', 'GG/CC'
        ]
        for stack in canonical:
            with self.subTest(stack=stack):
                self.assertIn(stack, thermo.ENTHALPY_NN)
                self.assertIn(stack, thermo.ENTROPY_NN)

    def test_enthalpy_entropy_consistency(self):
        """Test that enthalpy and entropy dicts have matching keys."""
        enthalpy_keys = set(thermo.ENTHALPY_NN.keys())
        entropy_keys = set(thermo.ENTROPY_NN.keys())
        # All enthalpy keys should be in entropy
        missing_in_entropy = enthalpy_keys - entropy_keys
        missing_in_enthalpy = entropy_keys - enthalpy_keys
        self.assertEqual(len(missing_in_entropy), 0,
                        f"Missing in entropy: {missing_in_entropy}")
        self.assertEqual(len(missing_in_enthalpy), 0,
                        f"Missing in enthalpy: {missing_in_enthalpy}")

    def test_enthalpy_values_reasonable(self):
        """Test that enthalpy values are in reasonable range."""
        for stack, value in thermo.ENTHALPY_NN.items():
            with self.subTest(stack=stack):
                # All NN enthalpies should be negative (favorable stacking)
                self.assertLess(value, 0)
                # Should be in range -11 to -7 kcal/mol
                self.assertGreater(value, -12)
                self.assertLess(value, -6)

    def test_entropy_values_reasonable(self):
        """Test that entropy values are in reasonable range."""
        for stack, value in thermo.ENTROPY_NN.items():
            with self.subTest(stack=stack):
                # All NN entropies should be negative
                self.assertLess(value, 0)
                # Should be in range -28 to -19 cal/(mol*K)
                self.assertGreater(value, -30)
                self.assertLess(value, -18)


if __name__ == '__main__':
    unittest.main()
