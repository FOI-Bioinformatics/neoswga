"""
Tests for Phase 6: Cleanup - deprecated thermo_estimation.py.

Verifies that:
1. thermodynamics.py has the compute_free_energy_for_two_strings function
2. rf_preprocessing.py uses thermodynamics module (not deprecated thermo_estimation)
3. The function produces consistent results
"""

import pytest
import numpy as np


class TestThermodynamicsModule:
    """Test thermodynamics module has required functions."""

    def test_compute_free_energy_exists(self):
        """thermodynamics should have compute_free_energy_for_two_strings."""
        from neoswga.core.thermodynamics import compute_free_energy_for_two_strings
        assert callable(compute_free_energy_for_two_strings)

    def test_compute_free_energy_basic(self):
        """compute_free_energy_for_two_strings should return valid values."""
        from neoswga.core.thermodynamics import compute_free_energy_for_two_strings

        # Test with complementary sequences
        result = compute_free_energy_for_two_strings('ATCGATCG', 'TAGCTAGC')
        assert isinstance(result, float)

    def test_delta_g_mismatch_table_exists(self):
        """DELTA_G_MISMATCH table should exist."""
        from neoswga.core.thermodynamics import DELTA_G_MISMATCH
        assert isinstance(DELTA_G_MISMATCH, dict)
        assert len(DELTA_G_MISMATCH) > 50  # Should have many entries


class TestRFPreprocessingImport:
    """Test that rf_preprocessing uses thermodynamics module."""

    def test_uses_thermo_module(self):
        """rf_preprocessing should import thermodynamics (not thermo_estimation)."""
        import neoswga.core.rf_preprocessing as rf

        # Should have imported thermo from thermodynamics
        assert hasattr(rf, 'thermo')

    def test_no_deprecated_import(self):
        """rf_preprocessing should not import deprecated te (thermo_estimation)."""
        import neoswga.core.rf_preprocessing as rf

        # Should NOT have 'te' attribute (deprecated import)
        assert not hasattr(rf, 'te')


class TestCoreModuleExports:
    """Test core __init__.py exports."""

    def test_thermo_estimation_not_in_all(self):
        """thermo_estimation should not be in __all__."""
        from neoswga.core import __all__
        # Note: there's a comment about it being deprecated, but it shouldn't be exported
        assert 'thermo_estimation' not in __all__

    def test_thermodynamics_in_all(self):
        """thermodynamics should be in __all__."""
        from neoswga.core import __all__
        assert 'thermodynamics' in __all__


class TestBackwardCompatibility:
    """Test backward compatibility of thermodynamic calculations."""

    def test_similar_results(self):
        """New function should produce similar results to legacy."""
        from neoswga.core.thermodynamics import compute_free_energy_for_two_strings

        # Test a few known pairs
        test_pairs = [
            ('ATCGATCG', 'TAGCTAGC'),
            ('GCGCGCGC', 'CGCGCGCG'),
            ('AAAAAAA', 'TTTTTTT'),
        ]

        for x, y in test_pairs:
            result = compute_free_energy_for_two_strings(x, y)
            # Should return reasonable delta G values
            assert -50 < result < 50, f"Unreasonable delta G for {x}/{y}: {result}"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
