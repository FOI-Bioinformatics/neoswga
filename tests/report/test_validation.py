"""
Tests for neoswga.core.report.validation module.

Tests validation of results directories and metrics completeness.
"""

import json
import pytest

from neoswga.core.report.validation import (
    validate_results_directory,
    validate_metrics,
    ValidationLevel,
)
from neoswga.core.report.metrics import (
    PipelineMetrics,
    CoverageMetrics,
    PrimerMetrics,
)


class TestValidateResultsDirectory:
    """Tests for validate_results_directory()."""

    def test_summary_json_missing_warns(self, tmp_path):
        """Missing step4_improved_df_summary.json produces a warning."""
        (tmp_path / 'step4_improved_df.csv').write_text('sequence\nATCG\n')
        (tmp_path / 'params.json').write_text('{}')
        result = validate_results_directory(str(tmp_path))
        warning_msgs = [i.message for i in result.warnings]
        assert any('step4_improved_df_summary.json' in m for m in warning_msgs)

    def test_summary_json_present_no_warn(self, tmp_path):
        """Present step4_improved_df_summary.json produces no warning about it."""
        (tmp_path / 'step4_improved_df.csv').write_text('sequence\nATCG\n')
        (tmp_path / 'params.json').write_text('{}')
        (tmp_path / 'step4_improved_df_summary.json').write_text('{}')
        result = validate_results_directory(str(tmp_path))
        warning_msgs = [i.message for i in result.warnings]
        assert not any('step4_improved_df_summary.json' in m for m in warning_msgs)


class TestValidateMetricsCoverageSource:
    """Tests for coverage source warnings in validate_metrics()."""

    def test_estimated_coverage_warns(self):
        """Estimated coverage (from_optimizer=False) produces a warning."""
        metrics = PipelineMetrics(
            results_dir='/tmp',
            generated_at='now',
            primers=[PrimerMetrics(
                sequence='ATCG', length=4, gc_content=0.5,
                tm=30.0, fg_freq=0.001, bg_freq=0.0001,
                fg_sites=10, bg_sites=1, gini=0.3, specificity=10.0,
            )],
            primer_count=1,
            coverage=CoverageMetrics(overall_coverage=0.5, from_optimizer=False),
        )
        result = validate_metrics(metrics)
        warning_msgs = [i.message for i in result.warnings]
        assert any('estimated' in m.lower() for m in warning_msgs)

    def test_measured_coverage_info(self):
        """Measured coverage (from_optimizer=True) produces an info message."""
        metrics = PipelineMetrics(
            results_dir='/tmp',
            generated_at='now',
            primers=[PrimerMetrics(
                sequence='ATCG', length=4, gc_content=0.5,
                tm=30.0, fg_freq=0.001, bg_freq=0.0001,
                fg_sites=10, bg_sites=1, gini=0.3, specificity=10.0,
            )],
            primer_count=1,
            coverage=CoverageMetrics(overall_coverage=0.85, from_optimizer=True),
        )
        result = validate_metrics(metrics)
        info_msgs = [i.message for i in result.infos]
        assert any('measured' in m.lower() for m in info_msgs)

    def test_no_coverage_no_crash(self):
        """No coverage object does not crash."""
        metrics = PipelineMetrics(
            results_dir='/tmp',
            generated_at='now',
            primers=[],
            primer_count=0,
            coverage=None,
        )
        result = validate_metrics(metrics)
        # Should not crash, just report primer error
        assert not result.is_valid  # no primers is an error
