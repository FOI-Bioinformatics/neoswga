"""
Comprehensive unit tests for reaction_conditions module extensions.

Tests new functionality including:
- Additional PCR additives (glycerol, BSA, PEG)
- Mg2+ optimization based on genome GC
- Commercial mix presets
- Parameter validation
"""

import unittest

from neoswga.core.reaction_conditions import (
    ReactionConditions,
    get_q_solution_equivalent,
    get_gc_melt_conditions,
    get_crude_sample_conditions
)


class TestNewAdditives(unittest.TestCase):
    """Test new additive parameters."""

    def test_glycerol_default(self):
        """Test glycerol defaults to 0."""
        conditions = ReactionConditions()
        self.assertEqual(conditions.glycerol_percent, 0.0)

    def test_glycerol_valid_range(self):
        """Test glycerol accepts valid values."""
        conditions = ReactionConditions(glycerol_percent=10.0)
        self.assertEqual(conditions.glycerol_percent, 10.0)

    def test_glycerol_max(self):
        """Test glycerol maximum value."""
        conditions = ReactionConditions(glycerol_percent=15.0)
        self.assertEqual(conditions.glycerol_percent, 15.0)

    def test_glycerol_too_high(self):
        """Test glycerol above maximum raises error."""
        with self.assertRaises(ValueError) as context:
            ReactionConditions(glycerol_percent=20.0)
        self.assertIn('Glycerol', str(context.exception))
        self.assertIn('15', str(context.exception))

    def test_glycerol_negative(self):
        """Test negative glycerol raises error."""
        with self.assertRaises(ValueError):
            ReactionConditions(glycerol_percent=-5.0)

    def test_bsa_default(self):
        """Test BSA defaults to 0."""
        conditions = ReactionConditions()
        self.assertEqual(conditions.bsa_ug_ml, 0.0)

    def test_bsa_valid_range(self):
        """Test BSA accepts valid values."""
        conditions = ReactionConditions(bsa_ug_ml=200.0)
        self.assertEqual(conditions.bsa_ug_ml, 200.0)

    def test_bsa_max(self):
        """Test BSA maximum value."""
        conditions = ReactionConditions(bsa_ug_ml=400.0)
        self.assertEqual(conditions.bsa_ug_ml, 400.0)

    def test_bsa_too_high(self):
        """Test BSA above maximum raises error."""
        with self.assertRaises(ValueError) as context:
            ReactionConditions(bsa_ug_ml=500.0)
        self.assertIn('BSA', str(context.exception))
        self.assertIn('400', str(context.exception))

    def test_bsa_negative(self):
        """Test negative BSA raises error."""
        with self.assertRaises(ValueError):
            ReactionConditions(bsa_ug_ml=-50.0)

    def test_peg_default(self):
        """Test PEG defaults to 0."""
        conditions = ReactionConditions()
        self.assertEqual(conditions.peg_percent, 0.0)

    def test_peg_valid_range(self):
        """Test PEG accepts valid values."""
        conditions = ReactionConditions(peg_percent=5.0)
        self.assertEqual(conditions.peg_percent, 5.0)

    def test_peg_max(self):
        """Test PEG maximum value."""
        conditions = ReactionConditions(peg_percent=15.0)
        self.assertEqual(conditions.peg_percent, 15.0)

    def test_peg_too_high(self):
        """Test PEG above maximum raises error."""
        with self.assertRaises(ValueError) as context:
            ReactionConditions(peg_percent=20.0)
        self.assertIn('PEG', str(context.exception))
        self.assertIn('15', str(context.exception))

    def test_peg_negative(self):
        """Test negative PEG raises error."""
        with self.assertRaises(ValueError):
            ReactionConditions(peg_percent=-5.0)

    # Ethanol tests
    def test_ethanol_default(self):
        """Test ethanol defaults to 0."""
        conditions = ReactionConditions()
        self.assertEqual(conditions.ethanol_percent, 0.0)

    def test_ethanol_valid_range(self):
        """Test ethanol accepts valid values."""
        conditions = ReactionConditions(ethanol_percent=2.5)
        self.assertEqual(conditions.ethanol_percent, 2.5)

    def test_ethanol_max(self):
        """Test ethanol maximum value."""
        conditions = ReactionConditions(ethanol_percent=5.0)
        self.assertEqual(conditions.ethanol_percent, 5.0)

    def test_ethanol_too_high(self):
        """Test ethanol above maximum raises error."""
        with self.assertRaises(ValueError) as context:
            ReactionConditions(ethanol_percent=6.0)
        self.assertIn('Ethanol', str(context.exception))
        self.assertIn('5', str(context.exception))

    def test_ethanol_negative(self):
        """Test negative ethanol raises error."""
        with self.assertRaises(ValueError):
            ReactionConditions(ethanol_percent=-1.0)

    # Urea tests
    def test_urea_default(self):
        """Test urea defaults to 0."""
        conditions = ReactionConditions()
        self.assertEqual(conditions.urea_m, 0.0)

    def test_urea_valid_range(self):
        """Test urea accepts valid values."""
        conditions = ReactionConditions(urea_m=1.0)
        self.assertEqual(conditions.urea_m, 1.0)

    def test_urea_max(self):
        """Test urea maximum value."""
        conditions = ReactionConditions(urea_m=2.0)
        self.assertEqual(conditions.urea_m, 2.0)

    def test_urea_too_high(self):
        """Test urea above maximum raises error."""
        with self.assertRaises(ValueError) as context:
            ReactionConditions(urea_m=2.5)
        self.assertIn('Urea', str(context.exception))
        self.assertIn('2.0', str(context.exception))

    def test_urea_negative(self):
        """Test negative urea raises error."""
        with self.assertRaises(ValueError):
            ReactionConditions(urea_m=-0.5)

    # TMAC tests
    def test_tmac_default(self):
        """Test TMAC defaults to 0."""
        conditions = ReactionConditions()
        self.assertEqual(conditions.tmac_m, 0.0)

    def test_tmac_valid_range(self):
        """Test TMAC accepts valid values."""
        conditions = ReactionConditions(tmac_m=0.05)
        self.assertEqual(conditions.tmac_m, 0.05)

    def test_tmac_max(self):
        """Test TMAC maximum value."""
        conditions = ReactionConditions(tmac_m=0.1)
        self.assertEqual(conditions.tmac_m, 0.1)

    def test_tmac_too_high(self):
        """Test TMAC above maximum raises error."""
        with self.assertRaises(ValueError) as context:
            ReactionConditions(tmac_m=0.2)
        self.assertIn('TMAC', str(context.exception))
        self.assertIn('0.1', str(context.exception))

    def test_tmac_negative(self):
        """Test negative TMAC raises error."""
        with self.assertRaises(ValueError):
            ReactionConditions(tmac_m=-0.01)

    # Formamide tests
    def test_formamide_default(self):
        """Test formamide defaults to 0."""
        conditions = ReactionConditions()
        self.assertEqual(conditions.formamide_percent, 0.0)

    def test_formamide_valid_range(self):
        """Test formamide accepts valid values."""
        conditions = ReactionConditions(formamide_percent=5.0)
        self.assertEqual(conditions.formamide_percent, 5.0)

    def test_formamide_max(self):
        """Test formamide maximum value."""
        conditions = ReactionConditions(formamide_percent=10.0)
        self.assertEqual(conditions.formamide_percent, 10.0)

    def test_formamide_too_high(self):
        """Test formamide above maximum raises error."""
        with self.assertRaises(ValueError) as context:
            ReactionConditions(formamide_percent=15.0)
        self.assertIn('Formamide', str(context.exception))
        self.assertIn('10', str(context.exception))

    def test_formamide_negative(self):
        """Test negative formamide raises error."""
        with self.assertRaises(ValueError):
            ReactionConditions(formamide_percent=-2.0)


class TestMg2Optimization(unittest.TestCase):
    """Test Mg2+ optimization method."""

    def test_at_rich_genome(self):
        """Test Mg2+ for AT-rich genome (<35% GC)."""
        conditions = ReactionConditions()
        mg_conc = conditions.optimize_mg_concentration(0.25)  # 25% GC

        self.assertEqual(mg_conc, 2.5)

    def test_at_rich_boundary(self):
        """Test Mg2+ at AT-rich boundary."""
        conditions = ReactionConditions()
        mg_conc = conditions.optimize_mg_concentration(0.34)  # 34% GC

        self.assertEqual(mg_conc, 2.5)

    def test_balanced_genome(self):
        """Test Mg2+ for balanced genome."""
        conditions = ReactionConditions()
        mg_conc = conditions.optimize_mg_concentration(0.50)  # 50% GC

        self.assertEqual(mg_conc, 2.0)

    def test_gc_rich_genome(self):
        """Test Mg2+ for GC-rich genome (>65% GC)."""
        conditions = ReactionConditions()
        mg_conc = conditions.optimize_mg_concentration(0.70)  # 70% GC

        self.assertEqual(mg_conc, 1.5)

    def test_gc_rich_boundary(self):
        """Test Mg2+ at GC-rich boundary."""
        conditions = ReactionConditions()
        mg_conc = conditions.optimize_mg_concentration(0.66)  # 66% GC

        self.assertEqual(mg_conc, 1.5)

    def test_francisella_genome(self):
        """Test Mg2+ for Francisella (33% GC)."""
        conditions = ReactionConditions()
        mg_conc = conditions.optimize_mg_concentration(0.33)

        self.assertEqual(mg_conc, 2.5)

    def test_burkholderia_genome(self):
        """Test Mg2+ for Burkholderia (67% GC)."""
        conditions = ReactionConditions()
        mg_conc = conditions.optimize_mg_concentration(0.67)

        self.assertEqual(mg_conc, 1.5)

    def test_ecoli_genome(self):
        """Test Mg2+ for E. coli (50% GC)."""
        conditions = ReactionConditions()
        mg_conc = conditions.optimize_mg_concentration(0.50)

        self.assertEqual(mg_conc, 2.0)

    def test_invalid_gc_too_low(self):
        """Test invalid GC content (<0)."""
        conditions = ReactionConditions()

        with self.assertRaises(ValueError):
            conditions.optimize_mg_concentration(-0.1)

    def test_invalid_gc_too_high(self):
        """Test invalid GC content (>1)."""
        conditions = ReactionConditions()

        with self.assertRaises(ValueError):
            conditions.optimize_mg_concentration(1.5)

    def test_boundary_gc_zero(self):
        """Test GC content of 0 (all AT)."""
        conditions = ReactionConditions()
        mg_conc = conditions.optimize_mg_concentration(0.0)

        self.assertEqual(mg_conc, 2.5)  # AT-rich

    def test_boundary_gc_one(self):
        """Test GC content of 1.0 (all GC)."""
        conditions = ReactionConditions()
        mg_conc = conditions.optimize_mg_concentration(1.0)

        self.assertEqual(mg_conc, 1.5)  # GC-rich


class TestCommercialPresets(unittest.TestCase):
    """Test commercial mix preset functions."""

    def test_q_solution_preset(self):
        """Test Q-Solution equivalent preset."""
        conditions = get_q_solution_equivalent()

        self.assertEqual(conditions.temp, 42.0)
        self.assertEqual(conditions.betaine_m, 1.5)
        self.assertEqual(conditions.glycerol_percent, 10.0)
        self.assertEqual(conditions.bsa_ug_ml, 200.0)
        self.assertEqual(conditions.polymerase, 'equiphi29')

    def test_gc_melt_preset(self):
        """Test GC-Melt preset."""
        conditions = get_gc_melt_conditions()

        self.assertEqual(conditions.temp, 45.0)
        self.assertEqual(conditions.dmso_percent, 5.0)
        self.assertEqual(conditions.betaine_m, 2.0)
        self.assertEqual(conditions.trehalose_m, 0.5)
        self.assertEqual(conditions.polymerase, 'equiphi29')

    def test_crude_sample_preset(self):
        """Test crude sample preset."""
        conditions = get_crude_sample_conditions()

        self.assertEqual(conditions.temp, 30.0)
        self.assertEqual(conditions.betaine_m, 1.0)
        self.assertEqual(conditions.glycerol_percent, 10.0)
        self.assertEqual(conditions.bsa_ug_ml, 400.0)
        self.assertEqual(conditions.peg_percent, 5.0)
        self.assertEqual(conditions.polymerase, 'phi29')

    def test_presets_are_valid(self):
        """Test that all presets pass validation."""
        # Should not raise errors
        conditions1 = get_q_solution_equivalent()
        conditions2 = get_gc_melt_conditions()
        conditions3 = get_crude_sample_conditions()

        # All should be valid
        self.assertIsInstance(conditions1, ReactionConditions)
        self.assertIsInstance(conditions2, ReactionConditions)
        self.assertIsInstance(conditions3, ReactionConditions)


class TestCombinedAdditives(unittest.TestCase):
    """Test using multiple new additives together."""

    def test_all_new_additives(self):
        """Test using glycerol, BSA, and PEG together."""
        conditions = ReactionConditions(
            glycerol_percent=10.0,
            bsa_ug_ml=200.0,
            peg_percent=5.0
        )

        self.assertEqual(conditions.glycerol_percent, 10.0)
        self.assertEqual(conditions.bsa_ug_ml, 200.0)
        self.assertEqual(conditions.peg_percent, 5.0)

    def test_new_with_existing_additives(self):
        """Test new additives with existing ones."""
        conditions = ReactionConditions(
            temp=42.0,
            polymerase='equiphi29',  # equiphi29 supports 40-50C
            betaine_m=1.5,
            dmso_percent=5.0,
            glycerol_percent=10.0,
            bsa_ug_ml=300.0,
            peg_percent=5.0
        )

        # Old additives
        self.assertEqual(conditions.betaine_m, 1.5)
        self.assertEqual(conditions.dmso_percent, 5.0)

        # New additives
        self.assertEqual(conditions.glycerol_percent, 10.0)
        self.assertEqual(conditions.bsa_ug_ml, 300.0)
        self.assertEqual(conditions.peg_percent, 5.0)

    def test_max_all_additives(self):
        """Test using maximum of all additives."""
        conditions = ReactionConditions(
            betaine_m=2.5,
            dmso_percent=10.0,
            trehalose_m=0.5,
            formamide_percent=5.0,
            glycerol_percent=15.0,
            bsa_ug_ml=400.0,
            peg_percent=15.0
        )

        # Should all be at max
        self.assertEqual(conditions.betaine_m, 2.5)
        self.assertEqual(conditions.dmso_percent, 10.0)
        self.assertEqual(conditions.glycerol_percent, 15.0)
        self.assertEqual(conditions.bsa_ug_ml, 400.0)
        self.assertEqual(conditions.peg_percent, 15.0)


class TestBackwardCompatibility(unittest.TestCase):
    """Test that new features don't break existing code."""

    def test_default_conditions_unchanged(self):
        """Test that default conditions still work."""
        conditions = ReactionConditions()

        # New additives should default to 0
        self.assertEqual(conditions.glycerol_percent, 0.0)
        self.assertEqual(conditions.bsa_ug_ml, 0.0)
        self.assertEqual(conditions.peg_percent, 0.0)

        # Existing defaults should be unchanged
        self.assertEqual(conditions.temp, 30.0)
        self.assertEqual(conditions.betaine_m, 0.0)
        self.assertEqual(conditions.dmso_percent, 0.0)

    def test_existing_conditions_work(self):
        """Test that existing condition creation still works."""
        # This should work exactly as before
        conditions = ReactionConditions(
            temp=42.0,
            betaine_m=1.5,
            dmso_percent=5.0,
            polymerase='equiphi29'
        )

        self.assertEqual(conditions.temp, 42.0)
        self.assertEqual(conditions.betaine_m, 1.5)
        self.assertEqual(conditions.dmso_percent, 5.0)
        self.assertEqual(conditions.polymerase, 'equiphi29')


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and boundary conditions."""

    def test_zero_additives(self):
        """Test all additives set to zero."""
        conditions = ReactionConditions(
            glycerol_percent=0.0,
            bsa_ug_ml=0.0,
            peg_percent=0.0
        )

        self.assertEqual(conditions.glycerol_percent, 0.0)
        self.assertEqual(conditions.bsa_ug_ml, 0.0)
        self.assertEqual(conditions.peg_percent, 0.0)

    def test_fractional_values(self):
        """Test fractional additive values."""
        conditions = ReactionConditions(
            glycerol_percent=7.5,
            bsa_ug_ml=150.5,
            peg_percent=2.5
        )

        self.assertEqual(conditions.glycerol_percent, 7.5)
        self.assertEqual(conditions.bsa_ug_ml, 150.5)
        self.assertEqual(conditions.peg_percent, 2.5)

    def test_extreme_gc_values(self):
        """Test Mg2+ optimization with extreme but valid GC values."""
        conditions = ReactionConditions()

        # Very AT-rich
        mg1 = conditions.optimize_mg_concentration(0.01)
        self.assertEqual(mg1, 2.5)

        # Very GC-rich
        mg2 = conditions.optimize_mg_concentration(0.99)
        self.assertEqual(mg2, 1.5)


class TestParameterConsistency(unittest.TestCase):
    """Test parameter consistency and relationships."""

    def test_mg_optimization_returns_float(self):
        """Test that Mg2+ optimization returns float."""
        conditions = ReactionConditions()
        mg_conc = conditions.optimize_mg_concentration(0.5)

        self.assertIsInstance(mg_conc, (int, float))
        self.assertGreater(mg_conc, 0)

    def test_preset_temperatures(self):
        """Test that presets use appropriate temperatures."""
        q_sol = get_q_solution_equivalent()
        gc_melt = get_gc_melt_conditions()
        crude = get_crude_sample_conditions()

        # All should have valid temperatures
        self.assertGreaterEqual(q_sol.temp, 20.0)
        self.assertLessEqual(q_sol.temp, 50.0)

        self.assertGreaterEqual(gc_melt.temp, 20.0)
        self.assertLessEqual(gc_melt.temp, 50.0)

        self.assertGreaterEqual(crude.temp, 20.0)
        self.assertLessEqual(crude.temp, 50.0)

    def test_preset_polymerase_choices(self):
        """Test that presets use valid polymerases."""
        q_sol = get_q_solution_equivalent()
        gc_melt = get_gc_melt_conditions()
        crude = get_crude_sample_conditions()

        valid_polymerases = ['phi29', 'equiphi29', 'bst', 'klenow']

        self.assertIn(q_sol.polymerase, valid_polymerases)
        self.assertIn(gc_melt.polymerase, valid_polymerases)
        self.assertIn(crude.polymerase, valid_polymerases)


class TestDocumentation(unittest.TestCase):
    """Test that functions have proper documentation."""

    def test_optimize_mg_has_docstring(self):
        """Test that optimize_mg_concentration has docstring."""
        conditions = ReactionConditions()
        self.assertIsNotNone(conditions.optimize_mg_concentration.__doc__)
        self.assertIn('genome', conditions.optimize_mg_concentration.__doc__.lower())

    def test_preset_functions_have_docstrings(self):
        """Test that preset functions have docstrings."""
        self.assertIsNotNone(get_q_solution_equivalent.__doc__)
        self.assertIsNotNone(get_gc_melt_conditions.__doc__)
        self.assertIsNotNone(get_crude_sample_conditions.__doc__)


if __name__ == '__main__':
    unittest.main()
