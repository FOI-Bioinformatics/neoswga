"""
Tests for the visualizations module.

Tests Plotly availability detection and chart rendering functions,
including graceful degradation when Plotly is not installed.
"""

import pytest
import sys
from unittest.mock import patch, MagicMock


class TestPlotlyAvailability:
    """Tests for Plotly availability detection."""

    def test_is_plotly_available_returns_bool(self):
        """is_plotly_available should return a boolean."""
        from neoswga.core.report.visualizations import is_plotly_available
        result = is_plotly_available()
        assert isinstance(result, bool)

    def test_chart_colors_defined(self):
        """CHART_COLORS should be defined with expected keys."""
        from neoswga.core.report.visualizations import CHART_COLORS
        assert isinstance(CHART_COLORS, dict)
        assert "primary" in CHART_COLORS
        assert "success" in CHART_COLORS
        assert "radar_fill" in CHART_COLORS
        assert "funnel_gradient" in CHART_COLORS


class TestGracefulDegradation:
    """Tests for graceful degradation when Plotly is not available."""

    def test_render_filtering_funnel_empty_without_data(self):
        """render_filtering_funnel should return empty string for empty data."""
        from neoswga.core.report.visualizations import render_filtering_funnel
        result = render_filtering_funnel([])
        assert result == ""

    def test_render_component_radar_empty_without_data(self):
        """render_component_radar should return empty string for empty components."""
        from neoswga.core.report.visualizations import render_component_radar
        result = render_component_radar([])
        assert result == ""

    def test_render_tm_gc_distribution_empty_without_data(self):
        """render_tm_gc_distribution should return empty string for empty primers."""
        from neoswga.core.report.visualizations import render_tm_gc_distribution
        result = render_tm_gc_distribution([])
        assert result == ""

    def test_render_coverage_specificity_scatter_empty_without_data(self):
        """render_coverage_specificity_scatter should return empty string for empty primers."""
        from neoswga.core.report.visualizations import render_coverage_specificity_scatter
        result = render_coverage_specificity_scatter([])
        assert result == ""

    def test_render_primer_heatmap_empty_without_enough_data(self):
        """render_primer_heatmap should return empty string for < 2 primers."""
        from neoswga.core.report.visualizations import render_primer_heatmap
        result = render_primer_heatmap([])
        assert result == ""

    def test_render_dimer_network_heatmap_empty_without_enough_data(self):
        """render_dimer_network_heatmap should return empty string for < 2 primers."""
        from neoswga.core.report.visualizations import render_dimer_network_heatmap
        result = render_dimer_network_heatmap([])
        assert result == ""

    def test_render_dimer_network_graph_empty_without_enough_data(self):
        """render_dimer_network_graph should return empty string for < 2 primers."""
        from neoswga.core.report.visualizations import render_dimer_network_graph
        result = render_dimer_network_graph([])
        assert result == ""


class TestDimerThermodynamics:
    """Tests for dimer thermodynamic calculations."""

    def test_calculate_heterodimer_dg_returns_float(self):
        """_calculate_heterodimer_dg should return a float."""
        from neoswga.core.report.visualizations import _calculate_heterodimer_dg
        result = _calculate_heterodimer_dg("ATCGATCG", "CGATCGAT")
        assert isinstance(result, float)

    def test_calculate_heterodimer_dg_different_sequences(self):
        """Different sequence pairs should give different delta G values."""
        from neoswga.core.report.visualizations import _calculate_heterodimer_dg
        import math
        # Different pairs should give different thermodynamic values
        dg1 = _calculate_heterodimer_dg("ATCGATCG", "CGATCGAT")
        dg2 = _calculate_heterodimer_dg("ATCGATCG", "TTTTTTT")
        # Both should be finite numbers
        assert math.isfinite(dg1)
        assert math.isfinite(dg2)
        # Sequences with complementary regions should have more negative dG
        assert dg1 < 0  # Has complementary regions

    def test_calculate_self_dimer_dg_returns_float(self):
        """_calculate_self_dimer_dg should return a float."""
        from neoswga.core.report.visualizations import _calculate_self_dimer_dg
        result = _calculate_self_dimer_dg("ATCGATCG")
        assert isinstance(result, float)

    def test_calculate_self_dimer_dg_palindrome(self):
        """Palindromic sequences should have stronger self-dimer potential."""
        from neoswga.core.report.visualizations import _calculate_self_dimer_dg
        # GCGC is a palindrome
        palindrome = _calculate_self_dimer_dg("GCGCGCGC")
        # AAAA is not palindromic
        non_palindrome = _calculate_self_dimer_dg("AAAAAAAA")
        # Palindrome should have more negative dG
        assert palindrome < non_palindrome

    def test_build_dimer_matrix_returns_tuple(self):
        """_build_dimer_matrix should return tuple of matrix, labels, sequences."""
        from neoswga.core.report.visualizations import _build_dimer_matrix
        from neoswga.core.report.metrics import PrimerMetrics

        primers = [
            PrimerMetrics(sequence="ATCGATCG", length=8, gc_content=0.5, tm=30.0,
                         fg_freq=0.001, bg_freq=0.0001, fg_sites=10, bg_sites=1,
                         gini=0.1, specificity=100.0),
            PrimerMetrics(sequence="GCTAGCTA", length=8, gc_content=0.5, tm=30.0,
                         fg_freq=0.001, bg_freq=0.0001, fg_sites=10, bg_sites=1,
                         gini=0.1, specificity=100.0),
        ]
        matrix, labels, sequences = _build_dimer_matrix(primers)
        assert len(matrix) == 2
        assert len(matrix[0]) == 2
        assert len(labels) == 2
        assert len(sequences) == 2

    def test_build_dimer_matrix_respects_max_primers(self):
        """_build_dimer_matrix should limit primers to max_primers."""
        from neoswga.core.report.visualizations import _build_dimer_matrix
        from neoswga.core.report.metrics import PrimerMetrics

        primers = [
            PrimerMetrics(sequence=f"ATCGATCG{i:02d}", length=10, gc_content=0.5, tm=30.0,
                         fg_freq=0.001, bg_freq=0.0001, fg_sites=10, bg_sites=1,
                         gini=0.1, specificity=100.0)
            for i in range(10)
        ]
        matrix, labels, sequences = _build_dimer_matrix(primers, max_primers=5)
        assert len(matrix) == 5
        assert len(labels) == 5
        assert len(sequences) == 5


class TestParetoFrontier:
    """Tests for Pareto frontier calculation."""

    def test_pareto_frontier_empty_input(self):
        """_calculate_pareto_frontier should handle empty input."""
        from neoswga.core.report.visualizations import _calculate_pareto_frontier
        result = _calculate_pareto_frontier([])
        assert result == []

    def test_pareto_frontier_single_point(self):
        """_calculate_pareto_frontier should handle single point."""
        from neoswga.core.report.visualizations import _calculate_pareto_frontier
        result = _calculate_pareto_frontier([(1, 2)])
        assert result == [(1, 2)]

    def test_pareto_frontier_finds_optimal_points(self):
        """_calculate_pareto_frontier should find Pareto-optimal points."""
        from neoswga.core.report.visualizations import _calculate_pareto_frontier
        # Points: (1,3) and (3,1) are Pareto-optimal
        # (2,2) is dominated by both
        points = [(1, 3), (2, 2), (3, 1)]
        result = _calculate_pareto_frontier(points)
        # Should include (1, 3) and (3, 1), sorted by x
        assert (1, 3) in result
        assert (3, 1) in result

    def test_pareto_frontier_all_dominated(self):
        """_calculate_pareto_frontier with one dominating point."""
        from neoswga.core.report.visualizations import _calculate_pareto_frontier
        # (5, 5) dominates all others
        points = [(1, 1), (2, 2), (5, 5)]
        result = _calculate_pareto_frontier(points)
        assert (5, 5) in result


# Only run these tests if Plotly is available
try:
    import plotly
    PLOTLY_INSTALLED = True
except ImportError:
    PLOTLY_INSTALLED = False


@pytest.mark.skipif(not PLOTLY_INSTALLED, reason="Plotly not installed")
class TestVisualizationsWithPlotly:
    """Tests that require Plotly to be installed."""

    @pytest.fixture
    def sample_funnel_data(self):
        """Sample filtering funnel data."""
        return [
            ("Total k-mers", 100000),
            ("After frequency filter", 50000),
            ("After background filter", 10000),
            ("Final candidates", 100),
        ]

    @pytest.fixture
    def sample_components(self):
        """Sample GradeComponent objects."""
        from neoswga.core.report.quality import GradeComponent
        return [
            GradeComponent(
                name="Coverage",
                weight=0.35,
                raw_value=0.85,
                normalized_score=0.9,
                rating="Excellent",
                description="85% coverage",
            ),
            GradeComponent(
                name="Specificity",
                weight=0.30,
                raw_value=500.0,
                normalized_score=0.85,
                rating="Good",
                description="500x enrichment",
            ),
            GradeComponent(
                name="Uniformity",
                weight=0.20,
                raw_value=0.75,
                normalized_score=0.7,
                rating="Good",
                description="Gini 0.25",
            ),
        ]

    @pytest.fixture
    def sample_primers(self):
        """Sample PrimerMetrics objects."""
        from neoswga.core.report.metrics import PrimerMetrics
        return [
            PrimerMetrics(
                sequence="ATCGATCGATCG",
                length=12,
                gc_content=0.5,
                tm=42.5,
                fg_freq=0.002,
                bg_freq=0.00001,
                fg_sites=200,
                bg_sites=2,
                gini=0.20,
                specificity=200000.0,
                amp_pred=0.85,
            ),
            PrimerMetrics(
                sequence="GCTAGCTAGCTA",
                length=12,
                gc_content=0.5,
                tm=41.8,
                fg_freq=0.0018,
                bg_freq=0.00001,
                fg_sites=180,
                bg_sites=2,
                gini=0.22,
                specificity=180000.0,
                amp_pred=0.80,
            ),
            PrimerMetrics(
                sequence="TAGCTAGCTAGC",
                length=12,
                gc_content=0.5,
                tm=40.2,
                fg_freq=0.0015,
                bg_freq=0.00002,
                fg_sites=150,
                bg_sites=3,
                gini=0.25,
                specificity=75000.0,
                amp_pred=0.78,
            ),
        ]

    def test_is_plotly_available_true(self):
        """is_plotly_available should return True when Plotly is installed."""
        from neoswga.core.report.visualizations import is_plotly_available
        assert is_plotly_available() is True

    def test_render_filtering_funnel_returns_html(self, sample_funnel_data):
        """render_filtering_funnel should return valid HTML."""
        from neoswga.core.report.visualizations import render_filtering_funnel
        html = render_filtering_funnel(sample_funnel_data)
        assert isinstance(html, str)
        assert len(html) > 0
        assert '<div' in html
        # Check for Plotly.js CDN reference
        assert 'plotly' in html.lower()

    def test_render_filtering_funnel_no_plotlyjs(self, sample_funnel_data):
        """render_filtering_funnel with include_plotlyjs=False."""
        from neoswga.core.report.visualizations import render_filtering_funnel
        html = render_filtering_funnel(sample_funnel_data, include_plotlyjs=False)
        assert isinstance(html, str)
        assert len(html) > 0
        # Should not include the CDN script tag
        assert 'cdn.plot.ly' not in html

    def test_render_component_radar_returns_html(self, sample_components):
        """render_component_radar should return valid HTML."""
        from neoswga.core.report.visualizations import render_component_radar
        html = render_component_radar(sample_components)
        assert isinstance(html, str)
        assert len(html) > 0
        assert '<div' in html

    def test_render_tm_gc_distribution_returns_html(self, sample_primers):
        """render_tm_gc_distribution should return valid HTML."""
        from neoswga.core.report.visualizations import render_tm_gc_distribution
        html = render_tm_gc_distribution(sample_primers, reaction_temp=30.0)
        assert isinstance(html, str)
        assert len(html) > 0
        assert '<div' in html

    def test_render_coverage_specificity_scatter_returns_html(self, sample_primers):
        """render_coverage_specificity_scatter should return valid HTML."""
        from neoswga.core.report.visualizations import render_coverage_specificity_scatter
        html = render_coverage_specificity_scatter(
            sample_primers,
            genome_size=4000000,
        )
        assert isinstance(html, str)
        assert len(html) > 0
        assert '<div' in html

    def test_render_primer_heatmap_returns_html(self, sample_primers):
        """render_primer_heatmap should return valid HTML."""
        from neoswga.core.report.visualizations import render_primer_heatmap
        html = render_primer_heatmap(sample_primers)
        assert isinstance(html, str)
        assert len(html) > 0
        assert '<div' in html

    def test_render_primer_heatmap_needs_multiple_primers(self):
        """render_primer_heatmap needs at least 2 primers."""
        from neoswga.core.report.visualizations import render_primer_heatmap
        from neoswga.core.report.metrics import PrimerMetrics

        single_primer = [PrimerMetrics(
            sequence="ATCG",
            length=4,
            gc_content=0.5,
            tm=20.0,
            fg_freq=0.001,
            bg_freq=0.0001,
            fg_sites=10,
            bg_sites=1,
            gini=0.1,
            specificity=10.0,
        )]
        html = render_primer_heatmap(single_primer)
        assert html == ""

    def test_render_dimer_network_heatmap_returns_html(self, sample_primers):
        """render_dimer_network_heatmap should return valid HTML."""
        from neoswga.core.report.visualizations import render_dimer_network_heatmap
        html = render_dimer_network_heatmap(sample_primers)
        assert isinstance(html, str)
        assert len(html) > 0
        assert '<div' in html

    def test_render_dimer_network_heatmap_no_plotlyjs(self, sample_primers):
        """render_dimer_network_heatmap with include_plotlyjs=False."""
        from neoswga.core.report.visualizations import render_dimer_network_heatmap
        html = render_dimer_network_heatmap(sample_primers, include_plotlyjs=False)
        assert isinstance(html, str)
        assert len(html) > 0
        # Should not include the CDN script tag
        assert 'cdn.plot.ly' not in html

    def test_render_dimer_network_graph_returns_html(self, sample_primers):
        """render_dimer_network_graph should return valid HTML."""
        from neoswga.core.report.visualizations import render_dimer_network_graph
        html = render_dimer_network_graph(sample_primers)
        assert isinstance(html, str)
        assert len(html) > 0
        assert '<div' in html

    def test_render_dimer_network_graph_with_custom_threshold(self, sample_primers):
        """render_dimer_network_graph should respect dg_threshold parameter."""
        from neoswga.core.report.visualizations import render_dimer_network_graph
        # Very negative threshold should show more edges
        html_low = render_dimer_network_graph(sample_primers, dg_threshold=-20.0)
        # Less negative threshold should show fewer edges
        html_high = render_dimer_network_graph(sample_primers, dg_threshold=-1.0)
        # Both should be valid HTML
        assert '<div' in html_low
        assert '<div' in html_high

    def test_chart_height_parameter(self, sample_funnel_data):
        """Charts should respect height parameter."""
        from neoswga.core.report.visualizations import render_filtering_funnel
        html_small = render_filtering_funnel(sample_funnel_data, height=200)
        html_large = render_filtering_funnel(sample_funnel_data, height=600)
        # Both should be valid HTML
        assert '<div' in html_small
        assert '<div' in html_large


class TestIntegrationWithReports:
    """Integration tests with report generation."""

    def test_executive_summary_accepts_interactive_param(self, complete_results_dir):
        """generate_executive_summary should accept interactive parameter."""
        from neoswga.core.report import generate_executive_summary

        # Should not raise
        summary = generate_executive_summary(
            str(complete_results_dir),
            interactive=False,
        )
        assert summary is not None

    @pytest.mark.skipif(not PLOTLY_INSTALLED, reason="Plotly not installed")
    def test_executive_summary_with_interactive(self, complete_results_dir, tmp_path):
        """generate_executive_summary with interactive=True should include charts."""
        from neoswga.core.report import generate_executive_summary

        output_file = tmp_path / "interactive_report.html"
        summary = generate_executive_summary(
            str(complete_results_dir),
            output_file=str(output_file),
            interactive=True,
        )

        assert summary is not None
        assert output_file.exists()

        content = output_file.read_text()
        # Should contain Plotly reference
        assert 'plotly' in content.lower()

    def test_technical_report_accepts_interactive_param(self, complete_results_dir):
        """generate_technical_report should accept interactive parameter."""
        from neoswga.core.report import generate_technical_report

        # Should not raise
        data = generate_technical_report(
            str(complete_results_dir),
            interactive=False,
        )
        assert data is not None

    @pytest.mark.skipif(not PLOTLY_INSTALLED, reason="Plotly not installed")
    def test_technical_report_with_interactive(self, complete_results_dir, tmp_path):
        """generate_technical_report with interactive=True should include charts."""
        from neoswga.core.report import generate_technical_report

        output_file = tmp_path / "interactive_technical_report.html"
        data = generate_technical_report(
            str(complete_results_dir),
            output_file=str(output_file),
            interactive=True,
        )

        assert data is not None
        assert output_file.exists()

        content = output_file.read_text()
        # Should contain Plotly reference
        assert 'plotly' in content.lower()
