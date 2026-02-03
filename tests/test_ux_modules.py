"""
Unit tests for UX improvement modules.

Tests:
- param_validator.py
- condition_suggester.py
- results_interpreter.py
"""

import pytest
import json
import tempfile
from pathlib import Path

from neoswga.core.param_validator import (
    ParamValidator,
    ValidationLevel,
    validate_params_file
)
from neoswga.core.condition_suggester import (
    ConditionSuggester,
    suggest_conditions,
    classify_gc,
    classify_primer_length
)
from neoswga.core.results_interpreter import (
    ResultsInterpreter,
    rate_metric,
    QualityRating,
    COVERAGE_THRESHOLDS
)


class TestParamValidator:
    """Test parameter validation."""

    def test_validate_required_params(self):
        """Test detection of missing required parameters."""
        validator = ParamValidator()

        # Missing fg_genomes
        params = {'fg_prefixes': ['data/test'], 'data_dir': 'results'}
        messages = validator.validate_params(params)

        errors = [m for m in messages if m.level == ValidationLevel.ERROR]
        assert any('fg_genomes' in m.parameter for m in errors)

    def test_validate_polymerase(self):
        """Test polymerase validation."""
        validator = ParamValidator()

        # Invalid polymerase
        params = {
            'fg_genomes': ['test.fasta'],
            'fg_prefixes': ['data/test'],
            'data_dir': 'results',
            'polymerase': 'invalid_poly'
        }
        messages = validator.validate_params(params)

        errors = [m for m in messages if m.level == ValidationLevel.ERROR]
        assert any('polymerase' in m.parameter for m in errors)

    def test_validate_temp_polymerase_compatibility(self):
        """Test temperature vs polymerase compatibility warning."""
        validator = ParamValidator()

        # phi29 at high temperature
        params = {
            'fg_genomes': ['test.fasta'],
            'fg_prefixes': ['data/test'],
            'data_dir': 'results',
            'polymerase': 'phi29',
            'reaction_temp': 50.0  # Too high for phi29
        }
        messages = validator.validate_params(params)

        warnings = [m for m in messages if m.level == ValidationLevel.WARNING]
        assert any('reaction_temp' in m.parameter for m in warnings)

    def test_validate_primer_additive_warning(self):
        """Test warning for long primers without additives."""
        validator = ParamValidator()

        params = {
            'fg_genomes': ['test.fasta'],
            'fg_prefixes': ['data/test'],
            'data_dir': 'results',
            'max_k': 16,
            'betaine_m': 0.0,
            'dmso_percent': 0.0
        }
        messages = validator.validate_params(params)

        warnings = [m for m in messages if m.level == ValidationLevel.WARNING]
        assert any('max_k' in m.parameter for m in warnings)

    def test_validate_range(self):
        """Test parameter range validation."""
        validator = ParamValidator()

        params = {
            'fg_genomes': ['test.fasta'],
            'fg_prefixes': ['data/test'],
            'data_dir': 'results',
            'dmso_percent': 15.0  # Too high (max 10)
        }
        messages = validator.validate_params(params)

        errors = [m for m in messages if m.level == ValidationLevel.ERROR]
        assert any('dmso_percent' in m.parameter for m in errors)

    def test_validate_file(self):
        """Test validation from file."""
        # Create temp params file
        params = {
            'fg_genomes': ['nonexistent.fasta'],
            'fg_prefixes': ['data/test'],
            'data_dir': 'results'
        }

        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump(params, f)
            temp_path = f.name

        try:
            success, messages = validate_params_file(temp_path, verbose=False)
            # Should have error for missing genome file
            assert not success
        finally:
            Path(temp_path).unlink()


class TestConditionSuggester:
    """Test condition suggestion logic."""

    def test_classify_gc(self):
        """Test GC classification."""
        assert classify_gc(0.19) == 'extreme_at'
        assert classify_gc(0.32) == 'at_rich'
        assert classify_gc(0.50) == 'balanced'
        assert classify_gc(0.67) == 'gc_rich'
        assert classify_gc(0.75) == 'extreme_gc'

    def test_classify_primer_length(self):
        """Test primer length classification."""
        assert classify_primer_length(8) == 'short'
        assert classify_primer_length(12) == 'medium'
        assert classify_primer_length(18) == 'long'

    def test_suggest_at_rich(self):
        """Test suggestions for AT-rich genome."""
        suggester = ConditionSuggester(genome_gc=0.25, primer_length=10)
        rec = suggester.suggest()

        # AT-rich should recommend phi29 for short primers
        assert rec.polymerase == 'phi29'
        assert rec.reaction_temp == 30.0

    def test_suggest_gc_rich(self):
        """Test suggestions for GC-rich genome."""
        suggester = ConditionSuggester(genome_gc=0.70, primer_length=15)
        rec = suggester.suggest()

        # GC-rich should recommend equiphi29 and high additives
        assert rec.polymerase == 'equiphi29'
        assert rec.betaine_m >= 1.5
        assert rec.dmso_percent >= 5.0

    def test_suggest_extreme_gc_warnings(self):
        """Test warnings for extreme GC genomes."""
        suggester = ConditionSuggester(genome_gc=0.80, primer_length=16)
        rec = suggester.suggest()

        assert len(rec.warnings) > 0
        assert any('extreme' in w.lower() or 'gc' in w.lower() for w in rec.warnings)

    def test_suggest_long_primers(self):
        """Test suggestions for long primers."""
        suggester = ConditionSuggester(genome_gc=0.50, primer_length=18)
        rec = suggester.suggest()

        # Long primers should require additives
        assert rec.betaine_m >= 1.0 or rec.dmso_percent >= 3.0

    def test_confidence_varies(self):
        """Test that confidence varies by scenario."""
        # Balanced genome should have higher confidence
        balanced = ConditionSuggester(genome_gc=0.50, primer_length=12).suggest()

        # Extreme genome should have lower confidence
        extreme = ConditionSuggester(genome_gc=0.80, primer_length=12).suggest()

        assert balanced.confidence > extreme.confidence


class TestResultsInterpreter:
    """Test results interpretation."""

    def test_rate_metric_coverage(self):
        """Test coverage rating."""
        assert rate_metric(0.98, COVERAGE_THRESHOLDS) == QualityRating.EXCELLENT
        assert rate_metric(0.90, COVERAGE_THRESHOLDS) == QualityRating.GOOD
        assert rate_metric(0.75, COVERAGE_THRESHOLDS) == QualityRating.ACCEPTABLE
        assert rate_metric(0.55, COVERAGE_THRESHOLDS) == QualityRating.POOR
        assert rate_metric(0.30, COVERAGE_THRESHOLDS) == QualityRating.CRITICAL

    def test_rate_metric_lower_is_better(self):
        """Test rating where lower is better (like Gini index)."""
        from neoswga.core.results_interpreter import UNIFORMITY_THRESHOLDS

        # Lower Gini is better
        assert rate_metric(0.2, UNIFORMITY_THRESHOLDS, lower_is_better=True) == QualityRating.EXCELLENT
        assert rate_metric(0.5, UNIFORMITY_THRESHOLDS, lower_is_better=True) == QualityRating.ACCEPTABLE
        assert rate_metric(0.9, UNIFORMITY_THRESHOLDS, lower_is_better=True) == QualityRating.CRITICAL

    def test_interpreter_missing_file(self):
        """Test interpreter with missing results."""
        interpreter = ResultsInterpreter('/nonexistent/path')

        with pytest.raises(FileNotFoundError):
            interpreter.analyze()


class TestIntegration:
    """Integration tests combining multiple modules."""

    def test_full_workflow_params(self):
        """Test creating and validating params."""
        # Create valid params
        params = {
            'fg_genomes': ['test.fasta'],
            'fg_prefixes': ['data/test'],
            'data_dir': 'results',
            'polymerase': 'equiphi29',
            'reaction_temp': 42.0,
            'min_k': 12,
            'max_k': 15,
            'betaine_m': 1.0,
            'dmso_percent': 5.0,
            'genome_gc': 0.65
        }

        # Validate (should have warnings but no errors about params themselves)
        validator = ParamValidator()
        messages = validator.validate_params(params)

        # Only file-related errors expected (test.fasta doesn't exist)
        param_errors = [m for m in messages
                       if m.level == ValidationLevel.ERROR
                       and 'fg_genomes' not in m.parameter]
        assert len(param_errors) == 0

    def test_suggestion_to_validation(self):
        """Test that suggested params pass validation."""
        # Get suggestion
        rec = suggest_conditions(genome_gc=0.50, primer_length=12, verbose=False)

        # Create params from suggestion
        params = {
            'fg_genomes': ['test.fasta'],
            'fg_prefixes': ['data/test'],
            'data_dir': 'results',
            'polymerase': rec.polymerase,
            'reaction_temp': rec.reaction_temp,
            'min_k': rec.primer_length_range[0],
            'max_k': rec.primer_length_range[1],
            'betaine_m': rec.betaine_m,
            'dmso_percent': rec.dmso_percent,
            'genome_gc': 0.50
        }

        # Validate
        validator = ParamValidator()
        messages = validator.validate_params(params)

        # Should have no param-related errors or warnings
        issues = [m for m in messages
                 if m.level in (ValidationLevel.ERROR, ValidationLevel.WARNING)
                 and 'fg_genomes' not in m.parameter]
        assert len(issues) == 0, f"Unexpected issues: {issues}"


class TestProgress:
    """Test progress indicators."""

    def test_progress_bar_basic(self):
        """Test basic progress bar usage."""
        from neoswga.core.progress import ProgressBar

        # Test with known total
        pbar = ProgressBar(total=10, desc="Testing", disable=True)
        pbar.n = 0
        pbar.update(5)
        assert pbar.n == 5

    def test_progress_bar_wrap(self):
        """Test progress bar wrap function."""
        from neoswga.core.progress import ProgressBar

        items = [1, 2, 3, 4, 5]
        result = []
        for item in ProgressBar.wrap(items, desc="Testing", disable=True):
            result.append(item)
        assert result == items

    def test_step_progress(self):
        """Test step progress indicator."""
        from neoswga.core.progress import StepProgress

        progress = StepProgress(
            steps=["Step 1", "Step 2", "Step 3"],
            disable=True
        )
        assert len(progress.steps) == 3
        assert progress.step_status == ['pending', 'pending', 'pending']

        progress.start_step(0)
        assert progress.step_status[0] == 'running'

        progress.complete_step(0)
        assert progress.step_status[0] == 'done'

    def test_progress_context(self):
        """Test progress context manager."""
        from neoswga.core.progress import progress_context

        # Should not raise
        with progress_context("Testing", disable=True):
            pass

    def test_progress_context_exception(self):
        """Test progress context with exception."""
        from neoswga.core.progress import progress_context

        with pytest.raises(ValueError):
            with progress_context("Testing", disable=True):
                raise ValueError("test error")

    def test_step_progress_fail(self):
        """Test step progress failure tracking."""
        from neoswga.core.progress import StepProgress

        progress = StepProgress(steps=["Step 1"], disable=True)
        progress.start_step(0)
        assert progress.step_status[0] == 'running'
        progress.fail_step(0, "Error message")
        assert progress.step_status[0] == 'failed'

    def test_progress_bar_context_manager(self):
        """Test progress bar as context manager."""
        from neoswga.core.progress import ProgressBar

        with ProgressBar(total=10, disable=True) as pbar:
            pbar.update(5)
            assert pbar.n == 5

    def test_format_time(self):
        """Test time formatting."""
        from neoswga.core.progress import _format_time

        assert _format_time(30) == "30s"
        assert _format_time(90) == "1m 30s"
        assert _format_time(3700) == "1h 1m"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
