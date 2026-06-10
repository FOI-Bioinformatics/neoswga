"""
Tests for GC-adaptive parameter application in pipeline.py.

Verifies that _apply_gc_adaptive_defaults correctly reads
betaine_concentration and dmso_concentration from GCAdaptiveParameters.
"""

import pytest
from types import SimpleNamespace
from unittest.mock import patch, MagicMock

from neoswga.core.gc_adaptive_strategy import GCAdaptiveParameters, GenomeClass


def _make_adaptive_params(betaine=1.5, dmso=5.0):
    """Create a GCAdaptiveParameters with specified additive concentrations."""
    return GCAdaptiveParameters(
        genome_class=GenomeClass.GC_RICH,
        gc_content=0.70,
        recommended_polymerase='equiphi29',
        kmer_range=(12, 18),
        optimal_kmer=15,
        reaction_temp=42.0,
        min_primer_tm=50.0,
        max_primer_tm=70.0,
        min_gc=0.40,
        max_gc=0.80,
        betaine_concentration=betaine,
        dmso_concentration=dmso,
        max_homodimer_dg=-5.0,
        max_heterodimer_dg=-5.0,
        max_hairpin_dg=-3.0,
        primer_count_multiplier=1.5,
        max_extension=50000,
        extension_rate=1000.0,
        confidence=0.85,
    )


class TestGCAdaptiveDefaults:
    """Tests for _apply_gc_adaptive_defaults in pipeline.py."""

    @patch('neoswga.core.gc_adaptive_strategy.GCAdaptiveStrategy')
    @patch('neoswga.core.pipeline.parameter')
    def test_betaine_and_dmso_applied(self, mock_param, MockStrategy):
        """Verify betaine_concentration and dmso_concentration are read correctly."""
        from neoswga.core.pipeline import _apply_gc_adaptive_defaults

        # Set up mock parameter module
        mock_param.genome_gc = 0.70
        mock_param.polymerase = 'phi29'
        mock_param.reaction_temp = None
        mock_param.betaine_m = 0.0
        mock_param.dmso_percent = 0.0
        mock_param._json_data = {}
        mock_param.min_k = 6
        mock_param.max_k = 12

        adaptive = _make_adaptive_params(betaine=1.5, dmso=5.0)
        MockStrategy.return_value.get_parameters.return_value = adaptive

        _apply_gc_adaptive_defaults()

        assert mock_param.betaine_m == 1.5
        assert mock_param.dmso_percent == 5.0

    @patch('neoswga.core.gc_adaptive_strategy.GCAdaptiveStrategy')
    @patch('neoswga.core.pipeline.parameter')
    def test_zero_concentrations_not_applied(self, mock_param, MockStrategy):
        """When adaptive params have zero concentrations, values stay at zero."""
        from neoswga.core.pipeline import _apply_gc_adaptive_defaults

        mock_param.genome_gc = 0.70
        mock_param.polymerase = 'phi29'
        mock_param.reaction_temp = None
        mock_param.betaine_m = 0.0
        mock_param.dmso_percent = 0.0
        mock_param._json_data = {}
        mock_param.min_k = 6
        mock_param.max_k = 12

        adaptive = _make_adaptive_params(betaine=0.0, dmso=0.0)
        MockStrategy.return_value.get_parameters.return_value = adaptive

        _apply_gc_adaptive_defaults()

        assert mock_param.betaine_m == 0.0
        assert mock_param.dmso_percent == 0.0

    def test_adaptive_params_has_correct_attribute_names(self):
        """GCAdaptiveParameters uses betaine_concentration and dmso_concentration."""
        params = _make_adaptive_params(betaine=1.0, dmso=3.0)
        # These must not raise AttributeError (the original bug)
        assert params.betaine_concentration == 1.0
        assert params.dmso_concentration == 3.0
        # The old incorrect names must not exist
        assert not hasattr(params, 'betaine_m')
        assert not hasattr(params, 'dmso_percent')
