"""
End-to-end tests for report generation.

Tests the full pipeline from results directory to HTML output
using realistic Bacillus genome test data.
"""

import pytest
from pathlib import Path

# Use the same import trick as conftest.py
import sys
project_root = Path(__file__).parent.parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

sys.modules.setdefault('neoswga', type(sys)('neoswga'))
sys.modules.setdefault('neoswga.core', type(sys)('neoswga.core'))
sys.modules['neoswga'].__path__ = [str(project_root / 'neoswga')]
sys.modules['neoswga.core'].__path__ = [str(project_root / 'neoswga' / 'core')]

from neoswga.core.report.metrics import collect_pipeline_metrics
from neoswga.core.report.quality import calculate_quality_grade, QualityGrade
from neoswga.core.report.executive_summary import (
    generate_executive_summary,
    create_executive_summary,
    render_executive_summary,
)
from neoswga.core.report.validation import (
    validate_results_directory,
    validate_metrics,
    ValidationLevel,
)


class TestE2EReportGeneration:
    """End-to-end tests with Bacillus genome fixtures."""

    def test_collect_metrics_from_bacillus(self, bacillus_results_dir):
        """Collect metrics from Bacillus test data."""
        metrics = collect_pipeline_metrics(str(bacillus_results_dir))

        assert len(metrics.primers) == 6
        assert metrics.primer_count == 6
        assert metrics.parameters.get("polymerase") == "phi29"
        assert metrics.parameters.get("fg_size") == 5220000

        # Check primer data
        first_primer = metrics.primers[0]
        assert first_primer.sequence == "TTTAAATTTTTA"
        assert first_primer.tm > 30.0
        assert first_primer.gc_content > 0

    def test_calculate_quality_grade_bacillus(self, bacillus_results_dir):
        """Calculate quality grade for Bacillus data."""
        metrics = collect_pipeline_metrics(str(bacillus_results_dir))
        quality = calculate_quality_grade(metrics)

        # Test data has realistic values - grade depends on coverage calculation
        # which is estimated from primer count
        assert quality.grade in [QualityGrade.A, QualityGrade.B, QualityGrade.C, QualityGrade.D]
        assert 0.0 <= quality.composite_score <= 1.0
        assert len(quality.components) == 5

        # Check all components present
        component_names = {c.name for c in quality.components}
        assert "Coverage" in component_names
        assert "Specificity" in component_names
        assert "Uniformity" in component_names
        assert "Thermodynamics" in component_names
        assert "Dimer Risk" in component_names

    def test_generate_executive_summary_bacillus(self, bacillus_results_dir, tmp_path):
        """Generate executive summary from Bacillus data."""
        output_file = tmp_path / "bacillus_summary.html"

        summary = generate_executive_summary(
            str(bacillus_results_dir),
            str(output_file)
        )

        # Check summary data
        assert summary.metrics.primer_count == 6
        assert summary.quality.grade is not None
        assert summary.generated_at != ""

        # Check output file created
        assert output_file.exists()
        html_content = output_file.read_text()
        assert "<!DOCTYPE html>" in html_content
        assert "NeoSWGA" in html_content

    def test_render_executive_summary_contains_key_sections(self, bacillus_results_dir):
        """Verify HTML output contains all required sections."""
        metrics = collect_pipeline_metrics(str(bacillus_results_dir))
        quality = calculate_quality_grade(metrics)
        summary = create_executive_summary(metrics, quality)
        html = render_executive_summary(summary)

        # Check structure
        assert "<!DOCTYPE html>" in html
        assert '<html lang="en">' in html

        # Check header
        assert "NeoSWGA Primer Set Report" in html

        # Check grade section
        assert "grade-letter" in html
        assert quality.grade.value in html  # A, B, C, D, or F

        # Check metrics grid
        assert "Coverage" in html
        assert "Enrichment" in html
        assert "Uniformity" in html
        assert "Dimer Risk" in html

        # Check primer table
        assert "Primer Set" in html
        assert "TTTAAATTTTTA" in html  # First primer sequence

        # Check recommendation
        assert quality.recommendation in html

    def test_html_output_valid_structure(self, bacillus_results_dir):
        """Verify HTML output has valid structure."""
        metrics = collect_pipeline_metrics(str(bacillus_results_dir))
        quality = calculate_quality_grade(metrics)
        summary = create_executive_summary(metrics, quality)
        html = render_executive_summary(summary)

        # Basic HTML validation
        assert html.count('<html') == 1
        assert html.count('</html>') == 1
        assert html.count('<head>') == 1
        assert html.count('</head>') == 1
        assert html.count('<body>') == 1
        assert html.count('</body>') == 1

        # Check for unclosed tags (basic check)
        assert html.count('<div') == html.count('</div>')
        assert html.count('<table>') == html.count('</table>')


class TestValidation:
    """Tests for validation utilities."""

    def test_validate_valid_directory(self, bacillus_results_dir):
        """Validate complete results directory."""
        result = validate_results_directory(str(bacillus_results_dir))

        assert result.is_valid
        assert len(result.errors) == 0

    def test_validate_empty_directory(self, tmp_path):
        """Validate empty directory."""
        result = validate_results_directory(str(tmp_path))

        assert not result.is_valid
        assert len(result.errors) > 0
        assert any("step4_improved_df.csv" in e.message for e in result.errors)

    def test_validate_nonexistent_directory(self, tmp_path):
        """Validate nonexistent directory."""
        result = validate_results_directory(str(tmp_path / "nonexistent"))

        assert not result.is_valid
        assert len(result.errors) == 1
        assert "not found" in result.errors[0].message

    def test_validate_metrics_complete(self, bacillus_results_dir):
        """Validate complete metrics."""
        metrics = collect_pipeline_metrics(str(bacillus_results_dir))
        result = validate_metrics(metrics)

        # May have warnings but no errors
        assert result.is_valid

    def test_validate_metrics_empty_primers(self, minimal_pipeline_metrics):
        """Validate metrics with no primers."""
        result = validate_metrics(minimal_pipeline_metrics)

        assert not result.is_valid
        assert any("No primers" in e.message for e in result.errors)


class TestEdgeCases:
    """Edge case tests for E2E functionality."""

    def test_minimal_step4_only(self, minimal_step4_csv):
        """Generate report with only step4 file."""
        metrics = collect_pipeline_metrics(str(minimal_step4_csv))
        quality = calculate_quality_grade(metrics)
        summary = create_executive_summary(metrics, quality)
        html = render_executive_summary(summary)

        # Should generate valid HTML
        assert "<!DOCTYPE html>" in html
        assert len(html) > 1000

    def test_report_with_missing_optional_files(self, tmp_path):
        """Generate report when optional files are missing."""
        # Create minimal step4
        step4 = tmp_path / "step4_improved_df.csv"
        step4.write_text(
            "sequence,score,fg_freq,bg_freq,tm,gini\n"
            "ATCGATCG,0.8,0.001,0.00001,35.0,0.2\n"
        )

        metrics = collect_pipeline_metrics(str(tmp_path))
        quality = calculate_quality_grade(metrics)
        summary = create_executive_summary(metrics, quality)
        html = render_executive_summary(summary)

        # Should generate valid HTML even with missing optional files
        assert "<!DOCTYPE html>" in html
        assert "ATCGATCG" in html

    def test_report_with_single_primer(self, tmp_path):
        """Generate report with single primer."""
        step4 = tmp_path / "step4_improved_df.csv"
        step4.write_text(
            "sequence,score,fg_freq,bg_freq,tm,gini\n"
            "ATCGATCG,0.8,0.001,0.00001,35.0,0.2\n"
        )

        metrics = collect_pipeline_metrics(str(tmp_path))
        quality = calculate_quality_grade(metrics)

        # Should handle single primer without error
        assert quality.grade is not None
        assert len(metrics.primers) == 1

    def test_report_with_edge_values(self, tmp_path):
        """Generate report with edge case values."""
        step4 = tmp_path / "step4_improved_df.csv"
        step4.write_text(
            "sequence,score,fg_freq,bg_freq,tm,gini,dimer_score\n"
            "GGGGGGGG,0.0,0.0,0.0,0.0,0.0,0.0\n"
            "CCCCCCCC,1.0,1.0,1.0,100.0,1.0,-20.0\n"
        )

        metrics = collect_pipeline_metrics(str(tmp_path))
        quality = calculate_quality_grade(metrics)

        # Should handle extreme values without error
        assert quality.grade is not None
        assert 0.0 <= quality.composite_score <= 1.0


class TestSecurityEscaping:
    """Tests for XSS and format string injection prevention."""

    def test_xss_prevention_in_primer_sequence(self, tmp_path):
        """Verify XSS is prevented in primer sequences."""
        step4 = tmp_path / "step4_improved_df.csv"
        step4.write_text(
            'sequence,score,fg_freq,bg_freq\n'
            '<script>alert("xss")</script>,0.8,0.001,0.00001\n'
        )

        metrics = collect_pipeline_metrics(str(tmp_path))
        quality = calculate_quality_grade(metrics)
        summary = create_executive_summary(metrics, quality)
        html = render_executive_summary(summary)

        # Script tag should be escaped
        assert '<script>' not in html
        assert '&lt;script&gt;' in html

    def test_format_string_injection_prevention(self, tmp_path):
        """Verify format string injection is prevented."""
        # Create params with braces in genome name
        params_file = tmp_path / "params.json"
        params_file.write_text('{"fg": "genome_{test}.fna"}')

        step4 = tmp_path / "step4_improved_df.csv"
        step4.write_text("sequence,score\nATCG,0.8\n")

        metrics = collect_pipeline_metrics(str(tmp_path))
        quality = calculate_quality_grade(metrics)
        summary = create_executive_summary(metrics, quality)

        # Should not raise KeyError from format string injection
        html = render_executive_summary(summary)
        assert "<!DOCTYPE html>" in html


class TestOutputQuality:
    """Tests for output quality and usability."""

    def test_primer_table_ordered(self, bacillus_results_dir):
        """Verify primers appear in expected order."""
        metrics = collect_pipeline_metrics(str(bacillus_results_dir))
        quality = calculate_quality_grade(metrics)
        summary = create_executive_summary(metrics, quality)
        html = render_executive_summary(summary)

        # First primer should appear before second
        first_pos = html.find("TTTAAATTTTTA")
        second_pos = html.find("AAAATTTAAAAT")

        assert first_pos < second_pos

    def test_grade_color_appropriate(self, bacillus_results_dir):
        """Verify grade colors are appropriate."""
        metrics = collect_pipeline_metrics(str(bacillus_results_dir))
        quality = calculate_quality_grade(metrics)
        summary = create_executive_summary(metrics, quality)
        html = render_executive_summary(summary)

        # Should have color styling
        assert 'style' in html or 'class=' in html

        # Grade letter should be present
        assert quality.grade.value in html

    def test_responsive_design_classes(self, bacillus_results_dir):
        """Verify responsive design is included."""
        metrics = collect_pipeline_metrics(str(bacillus_results_dir))
        quality = calculate_quality_grade(metrics)
        summary = create_executive_summary(metrics, quality)
        html = render_executive_summary(summary)

        # Should have media queries for responsive design
        assert "@media" in html
        assert "max-width" in html
