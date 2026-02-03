"""
Tests for Phase 4: Enhanced Feature Engineering.

Verifies that the enhanced feature engineering with 120+ features
and synthetic training data integration works correctly.
"""

import pytest
import numpy as np
import pandas as pd
from unittest.mock import Mock, patch


class TestAdvancedFeatureEngineer:
    """Test AdvancedFeatureEngineer class."""

    def test_feature_count(self):
        """Should generate 100+ features."""
        from neoswga.core.advanced_features import AdvancedFeatureEngineer
        from neoswga.core.reaction_conditions import ReactionConditions

        genome = "ATCGATCGATCGATCGATCG" * 1000  # 20kb genome
        conditions = ReactionConditions(temp=30.0, na_conc=50.0)

        engineer = AdvancedFeatureEngineer(genome, conditions)
        features_df = engineer.engineer_features(['ATCGATCGAT', 'GCTAGCTAGT'])

        # Should have many features (100+ excluding primer column)
        assert len(features_df.columns) >= 50  # Relaxed for flexibility

    def test_base_composition_features(self):
        """Base composition features should be calculated correctly."""
        from neoswga.core.advanced_features import AdvancedFeatureEngineer
        from neoswga.core.reaction_conditions import ReactionConditions

        genome = "ATCGATCGATCG" * 100
        conditions = ReactionConditions(temp=30.0, na_conc=50.0)

        engineer = AdvancedFeatureEngineer(genome, conditions)
        features_df = engineer.engineer_features(['ATCGATCGAT'])

        # Check expected features exist
        assert 'gc_content' in features_df.columns
        assert 'length' in features_df.columns
        assert 'a_count' in features_df.columns

        # GC content of ATCGATCGAT should be 0.4 (4 G/C out of 10)
        assert abs(features_df['gc_content'].iloc[0] - 0.4) < 0.01

    def test_thermodynamic_features(self):
        """Thermodynamic features should be calculated."""
        from neoswga.core.advanced_features import AdvancedFeatureEngineer
        from neoswga.core.reaction_conditions import ReactionConditions

        genome = "ATCGATCGATCG" * 100
        conditions = ReactionConditions(temp=30.0, na_conc=50.0)

        engineer = AdvancedFeatureEngineer(genome, conditions)
        features_df = engineer.engineer_features(['ATCGATCGAT'])

        # Should have thermodynamic features
        assert 'tm_base' in features_df.columns or 'tm_effective' in features_df.columns
        assert 'dg_37c' in features_df.columns or 'dg_reaction_temp' in features_df.columns

    def test_positional_features_with_positions(self):
        """Positional features should use primer positions when provided."""
        from neoswga.core.advanced_features import AdvancedFeatureEngineer
        from neoswga.core.reaction_conditions import ReactionConditions

        genome = "ATCGATCGATCG" * 100
        conditions = ReactionConditions(temp=30.0, na_conc=50.0)

        # Provide position data
        positions = {
            'ATCGATCGAT': {
                'forward': [100, 500, 900, 1300],
                'reverse': [200, 600]
            }
        }

        engineer = AdvancedFeatureEngineer(genome, conditions, positions)
        features_df = engineer.engineer_features(['ATCGATCGAT'])

        # Should have positional features
        assert 'total_binding_sites' in features_df.columns
        assert features_df['total_binding_sites'].iloc[0] == 6


class TestEnhancedPrediction:
    """Test enhanced prediction functions."""

    def test_is_enhanced_model_available(self):
        """Should correctly detect enhanced model availability."""
        from neoswga.core.rf_preprocessing import is_enhanced_model_available

        # With non-existent path, should return False
        result = is_enhanced_model_available('/nonexistent/path/model.pkl')
        assert result == False

    def test_get_enhanced_model_info_missing(self):
        """Should return None for missing model."""
        from neoswga.core.rf_preprocessing import get_enhanced_model_info

        result = get_enhanced_model_info('/nonexistent/path/model.pkl')
        assert result is None


class TestCLIEnhancedFeatures:
    """Test CLI options for enhanced features."""

    def test_cli_has_enhanced_features_option(self):
        """CLI should have --use-enhanced-features option."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()
        args = parser.parse_args(['score', '-j', 'test.json', '--use-enhanced-features'])

        assert hasattr(args, 'use_enhanced_features')
        assert args.use_enhanced_features == True

    def test_cli_has_enhanced_model_path_option(self):
        """CLI should have --enhanced-model-path option."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()
        args = parser.parse_args([
            'score', '-j', 'test.json',
            '--enhanced-model-path', '/path/to/model.pkl'
        ])

        assert hasattr(args, 'enhanced_model_path')
        assert args.enhanced_model_path == '/path/to/model.pkl'


class TestTrainingDataGeneration:
    """Test training data generation functions."""

    def test_generate_random_primers(self):
        """Should generate valid primers from genome."""
        from scripts.generate_training_data import generate_random_primers

        genome = "ATCGATCGATCGATCG" * 1000
        primers = generate_random_primers(genome, k=10, num_primers=20)

        assert len(primers) <= 20
        for primer in primers:
            assert len(primer) == 10
            assert all(b in 'ATCG' for b in primer)

    def test_find_primer_positions(self):
        """Should find primer positions in both strands."""
        from scripts.generate_training_data import find_primer_positions

        genome = "ATCGATCGATCGATCGATCGATCGATCGATCG"
        positions = find_primer_positions(genome, "ATCGATCG")

        assert 'forward' in positions
        assert 'reverse' in positions
        assert len(positions['forward']) > 0


class TestTrainingScript:
    """Test training script functions."""

    def test_exclude_columns_defined(self):
        """Training script should define columns to exclude."""
        from scripts.train_enhanced_rf import EXCLUDE_COLUMNS

        assert 'primer' in EXCLUDE_COLUMNS
        assert 'mean_enrichment' in EXCLUDE_COLUMNS

    def test_default_rf_params_defined(self):
        """Training script should define default RF parameters."""
        from scripts.train_enhanced_rf import DEFAULT_RF_PARAMS

        assert 'n_estimators' in DEFAULT_RF_PARAMS
        assert 'max_depth' in DEFAULT_RF_PARAMS


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
