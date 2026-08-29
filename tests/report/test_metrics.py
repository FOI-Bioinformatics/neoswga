"""
Unit tests for neoswga.core.report.metrics module.

Tests safe type conversion, PrimerMetrics parsing, and metrics collection.
"""

import json
import math
import pytest
from pathlib import Path

from neoswga.core.report.metrics import (
    _safe_float,
    _safe_int,
    _load_csv,
    _load_optimizer_summary,
    PrimerMetrics,
    PipelineMetrics,
    CoverageMetrics,
    FilteringStats,
    collect_pipeline_metrics,
)


class TestSafeFloat:
    """Tests for _safe_float() function."""

    def test_valid_float(self):
        """Normal float values."""
        assert _safe_float(3.14) == 3.14
        assert _safe_float(0.0) == 0.0
        assert _safe_float(-2.5) == -2.5

    def test_valid_int(self):
        """Integer values converted to float."""
        assert _safe_float(42) == 42.0
        assert _safe_float(0) == 0.0
        assert _safe_float(-10) == -10.0

    def test_valid_string(self):
        """String representations of floats."""
        assert _safe_float("3.14") == 3.14
        assert _safe_float("0.0") == 0.0
        assert _safe_float("-2.5") == -2.5
        assert _safe_float("1e-5") == 1e-5

    def test_nan_returns_default(self):
        """NaN values return default."""
        assert _safe_float(float("nan")) == 0.0
        assert _safe_float(float("nan"), default=1.0) == 1.0
        assert _safe_float("nan") == 0.0
        assert _safe_float("NaN") == 0.0
        assert _safe_float("NAN") == 0.0

    def test_inf_returns_default(self):
        """Infinity values return default."""
        assert _safe_float(float("inf")) == 0.0
        assert _safe_float(float("-inf")) == 0.0
        assert _safe_float("inf") == 0.0
        assert _safe_float("Inf") == 0.0
        assert _safe_float("-inf") == 0.0
        assert _safe_float("infinity") == 0.0
        assert _safe_float("-infinity") == 0.0

    def test_none_returns_default(self):
        """None returns default."""
        assert _safe_float(None) == 0.0
        assert _safe_float(None, default=5.0) == 5.0

    def test_empty_string_returns_default(self):
        """Empty string returns default."""
        assert _safe_float("") == 0.0
        assert _safe_float("", default=2.0) == 2.0

    def test_invalid_string_returns_default(self):
        """Invalid strings return default."""
        assert _safe_float("abc") == 0.0
        assert _safe_float("not a number") == 0.0
        assert _safe_float("none") == 0.0
        assert _safe_float("null") == 0.0

    def test_custom_default(self):
        """Custom default values work."""
        assert _safe_float(None, default=-1.0) == -1.0
        assert _safe_float("invalid", default=999.0) == 999.0


class TestSafeInt:
    """Tests for _safe_int() function."""

    def test_valid_int(self):
        """Normal integer values."""
        assert _safe_int(42) == 42
        assert _safe_int(0) == 0
        assert _safe_int(-10) == -10

    def test_valid_string(self):
        """String representations of integers."""
        assert _safe_int("42") == 42
        assert _safe_int("0") == 0
        assert _safe_int("-10") == -10

    def test_float_string_truncated(self):
        """Float strings are truncated to int."""
        assert _safe_int("12.5") == 12
        assert _safe_int("12.9") == 12
        assert _safe_int("0.999") == 0

    def test_float_value_truncated(self):
        """Float values are truncated to int."""
        assert _safe_int(12.5) == 12
        assert _safe_int(12.9) == 12

    def test_none_returns_default(self):
        """None returns default."""
        assert _safe_int(None) == 0
        assert _safe_int(None, default=5) == 5

    def test_empty_string_returns_default(self):
        """Empty string returns default."""
        assert _safe_int("") == 0
        assert _safe_int("", default=10) == 10

    def test_invalid_string_returns_default(self):
        """Invalid strings return default."""
        assert _safe_int("abc") == 0
        assert _safe_int("not a number") == 0


class TestLoadCSV:
    """Tests for _load_csv() function."""

    def test_valid_csv(self, tmp_path):
        """Load valid CSV file."""
        csv_file = tmp_path / "test.csv"
        csv_file.write_text("name,value\nfoo,1\nbar,2\n")

        rows = _load_csv(csv_file)
        assert len(rows) == 2
        assert rows[0]["name"] == "foo"
        assert rows[0]["value"] == "1"
        assert rows[1]["name"] == "bar"

    def test_nonexistent_file(self, tmp_path):
        """Nonexistent file returns empty list."""
        csv_file = tmp_path / "does_not_exist.csv"
        rows = _load_csv(csv_file)
        assert rows == []

    def test_empty_file(self, tmp_path):
        """Empty file returns empty list."""
        csv_file = tmp_path / "empty.csv"
        csv_file.write_text("")
        rows = _load_csv(csv_file)
        assert rows == []

    def test_header_only(self, tmp_path):
        """File with only header returns empty list."""
        csv_file = tmp_path / "header_only.csv"
        csv_file.write_text("col1,col2,col3\n")
        rows = _load_csv(csv_file)
        assert rows == []


class TestPrimerMetricsFromRow:
    """Tests for PrimerMetrics.from_row() class method."""

    def test_complete_row(self):
        """Parse row with all fields."""
        row = {
            "sequence": "ATCGATCG",
            "fg_freq": "0.001",
            "bg_freq": "0.00001",
            "tm": "35.5",
            "gini": "0.25",
            "gc": "0.5",
            "fg_count": "100",
            "bg_count": "1",
            "amp_pred": "17.0",
            "dimer_score": "-3.5",
            "hairpin_dg": "-1.2",
            "self_dimer_dg": "-2.8",
            "strand_ratio": "1.2",
        }

        primer = PrimerMetrics.from_row(row)

        assert primer.sequence == "ATCGATCG"
        assert primer.length == 8
        assert primer.fg_freq == 0.001
        assert primer.bg_freq == 0.00001
        assert primer.tm == 35.5
        assert primer.gini == 0.25
        assert primer.gc_content == 0.5
        assert primer.fg_sites == 100
        assert primer.bg_sites == 1
        assert primer.amp_pred == 0.85  # 17.0 / 20.0 = 0.85 (normalized)
        assert primer.dimer_score == -3.5
        assert primer.strand_ratio == 1.2

    def test_alternative_column_names(self):
        """Parse row with alternative column names."""
        row = {
            "primer": "ATCGATCG",
            "Tm": "35.5",
            "gini_index": "0.25",
            "gc_content": "0.5",
            "fg_sites": "100",
            "bg_sites": "1",
            "on.target.pred": "17.0",
        }

        primer = PrimerMetrics.from_row(row)

        assert primer.sequence == "ATCGATCG"
        assert primer.tm == 35.5
        assert primer.gini == 0.25
        assert primer.gc_content == 0.5
        assert primer.fg_sites == 100
        assert primer.amp_pred == 0.85  # 17.0 / 20.0 = 0.85 (normalized)

    def test_gc_calculated_from_sequence(self):
        """GC content calculated when not provided."""
        row = {
            "sequence": "GGCC",  # 100% GC
        }

        primer = PrimerMetrics.from_row(row)
        assert primer.gc_content == 1.0

        row2 = {"sequence": "ATAT"}  # 0% GC
        primer2 = PrimerMetrics.from_row(row2)
        assert primer2.gc_content == 0.0

        row3 = {"sequence": "ATGC"}  # 50% GC
        primer3 = PrimerMetrics.from_row(row3)
        assert primer3.gc_content == 0.5

    def test_gc_case_insensitive(self):
        """GC calculation is case-insensitive."""
        row = {"sequence": "atgc"}  # lowercase
        primer = PrimerMetrics.from_row(row)
        assert primer.gc_content == 0.5

        row2 = {"sequence": "AtGc"}  # mixed case
        primer2 = PrimerMetrics.from_row(row2)
        assert primer2.gc_content == 0.5

    def test_missing_columns_use_defaults(self):
        """Missing columns use safe defaults."""
        row = {"sequence": "ATCG"}

        primer = PrimerMetrics.from_row(row)

        assert primer.sequence == "ATCG"
        assert primer.fg_freq == 0.0
        assert primer.bg_freq == 0.0
        assert primer.tm == 0.0
        assert primer.gini == 0.0
        assert primer.strand_ratio == 1.0  # default for strand_ratio

    def test_invalid_values_use_defaults(self, edge_case_primer_row):
        """Invalid values use safe defaults."""
        primer = PrimerMetrics.from_row(edge_case_primer_row)

        assert primer.sequence == "atcgatcg"
        assert primer.fg_freq == 0.0  # nan -> default
        assert primer.bg_freq == 0.0  # empty -> default
        assert primer.tm == 0.0  # inf -> default
        assert primer.fg_sites == 0  # "abc" -> default
        assert primer.bg_sites == 12  # "12.5" truncated
        assert primer.dimer_score == 0.0  # "None" -> default
        assert primer.strand_ratio == 1.0  # "null" -> default

    def test_specificity_calculated(self):
        """Specificity is calculated as fg_freq/bg_freq."""
        row = {
            "sequence": "ATCG",
            "fg_freq": "0.001",
            "bg_freq": "0.00001",
        }

        primer = PrimerMetrics.from_row(row)

        assert primer.specificity == pytest.approx(100.0, rel=0.01)

    def test_specificity_zero_background(self):
        """Specificity handles zero background safely."""
        row = {
            "sequence": "ATCG",
            "fg_freq": "0.001",
            "bg_freq": "0",
        }

        primer = PrimerMetrics.from_row(row)

        # Should not raise, uses max(bg_freq, 1e-10)
        assert primer.specificity > 0
        assert math.isfinite(primer.specificity)

    def test_dimer_risk_score_fallback(self):
        """dimer_risk_score is used when dimer_score is absent."""
        row = {
            "sequence": "ATCGATCG",
            "dimer_risk_score": "-6.2",
        }

        primer = PrimerMetrics.from_row(row)
        assert primer.dimer_score == -6.2

    def test_dimer_score_preferred_over_risk_score(self):
        """dimer_score takes precedence over dimer_risk_score."""
        row = {
            "sequence": "ATCGATCG",
            "dimer_score": "-3.5",
            "dimer_risk_score": "-6.2",
        }

        primer = PrimerMetrics.from_row(row)
        assert primer.dimer_score == -3.5


class TestCollectPipelineMetrics:
    """Tests for collect_pipeline_metrics() function."""

    def test_valid_results_dir(self, complete_results_dir):
        """Collect metrics from valid results directory."""
        metrics = collect_pipeline_metrics(str(complete_results_dir))

        assert len(metrics.primers) == 6
        assert metrics.primer_count == 6
        assert metrics.parameters.get("polymerase") == "phi29"
        assert metrics.parameters.get("fg_size") == 4000000

    def test_minimal_results_dir(self, minimal_step4_csv):
        """Collect metrics from minimal results."""
        metrics = collect_pipeline_metrics(str(minimal_step4_csv))

        assert len(metrics.primers) == 2
        assert metrics.primers[0].sequence == "ATCGATCG"
        assert metrics.primers[0].tm == 35.5

    def test_nonexistent_directory_raises(self, tmp_path):
        """Nonexistent directory raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            collect_pipeline_metrics(str(tmp_path / "nonexistent"))

    def test_empty_directory(self, tmp_path):
        """Empty directory returns metrics with empty primers."""
        metrics = collect_pipeline_metrics(str(tmp_path))
        assert len(metrics.primers) == 0
        assert metrics.primer_count == 0

    def test_filtering_stats_loaded(self, complete_results_dir):
        """Filtering stats loaded from filter_stats.json."""
        metrics = collect_pipeline_metrics(str(complete_results_dir))

        assert metrics.filtering is not None
        assert metrics.filtering.total_kmers == 100000
        assert metrics.filtering.final_candidates == 6

    def test_filtering_stats_invalid_fields_filtered(self, tmp_path):
        """Invalid fields in filter_stats.json are ignored."""
        # Create minimal step4
        step4 = tmp_path / "step4_improved_df.csv"
        step4.write_text("sequence\nATCG\n")

        # Create filter_stats with extra invalid fields
        filter_stats = tmp_path / "filter_stats.json"
        filter_stats.write_text(
            json.dumps(
                {
                    "total_kmers": 1000,
                    "final_candidates": 10,
                    "invalid_field": "should be ignored",
                    "another_invalid": 123,
                }
            )
        )

        metrics = collect_pipeline_metrics(str(tmp_path))

        assert metrics.filtering is not None
        assert metrics.filtering.total_kmers == 1000
        # Should not have raised an error

    def test_coverage_metrics_calculated(self, complete_results_dir):
        """Coverage metrics are calculated."""
        metrics = collect_pipeline_metrics(str(complete_results_dir))

        assert metrics.coverage is not None
        assert metrics.coverage.total_bases == 4000000

    def test_specificity_metrics_calculated(self, complete_results_dir):
        """Specificity metrics are calculated."""
        metrics = collect_pipeline_metrics(str(complete_results_dir))

        assert metrics.specificity is not None
        assert metrics.specificity.target_sites > 0

    def test_thermodynamic_metrics_calculated(self, complete_results_dir):
        """Thermodynamic metrics are calculated."""
        metrics = collect_pipeline_metrics(str(complete_results_dir))

        assert metrics.thermodynamics is not None
        assert metrics.thermodynamics.mean_tm > 0
        assert metrics.thermodynamics.tm_range >= 0

    def test_uniformity_metrics_calculated(self, complete_results_dir):
        """Uniformity metrics are calculated."""
        metrics = collect_pipeline_metrics(str(complete_results_dir))

        assert metrics.uniformity is not None
        assert metrics.uniformity.mean_gini >= 0
        assert metrics.uniformity.strand_ratio > 0


class TestFilteringStats:
    """Tests for FilteringStats dataclass."""

    def test_as_funnel(self):
        """Test funnel representation."""
        stats = FilteringStats(
            total_kmers=1000,
            after_frequency=500,
            after_background=200,
            after_gini=100,
            after_thermodynamic=50,
            after_complexity=25,
            final_candidates=10,
        )

        funnel = stats.as_funnel()

        assert len(funnel) == 7
        assert funnel[0] == ("Total k-mers", 1000)
        assert funnel[-1] == ("Final candidates", 10)


class TestEdgeCases:
    """Edge case tests for metrics collection."""

    def test_all_nan_tm_values(self, tmp_path):
        """Handle all NaN Tm values."""
        csv_path = tmp_path / "step4_improved_df.csv"
        csv_path.write_text("sequence,tm\n" "ATCG,nan\n" "GCTA,NaN\n")

        metrics = collect_pipeline_metrics(str(tmp_path))

        assert len(metrics.primers) == 2
        assert metrics.primers[0].tm == 0.0
        assert metrics.thermodynamics is not None

    def test_unicode_in_csv(self, tmp_path):
        """Handle unicode characters in CSV (should not crash)."""
        csv_path = tmp_path / "step4_improved_df.csv"
        csv_path.write_text("sequence,comment\n" "ATCG,test primer\n")

        metrics = collect_pipeline_metrics(str(tmp_path))
        assert len(metrics.primers) == 1

    def test_very_large_numbers(self, tmp_path):
        """Handle very large numbers."""
        csv_path = tmp_path / "step4_improved_df.csv"
        csv_path.write_text("sequence,fg_freq,bg_freq\n" "ATCG,1e10,1e-15\n")

        metrics = collect_pipeline_metrics(str(tmp_path))
        assert metrics.primers[0].fg_freq == 1e10
        assert metrics.primers[0].bg_freq == 1e-15

    def test_strand_ratio_edge_cases(self, tmp_path):
        """Test strand ratio calculation edge cases."""
        # All primers favor forward strand
        csv_path = tmp_path / "step4_improved_df.csv"
        csv_path.write_text("sequence,strand_ratio\n" "ATCG,1.5\n" "GCTA,2.0\n" "TAGC,1.1\n")

        metrics = collect_pipeline_metrics(str(tmp_path))
        assert metrics.uniformity is not None
        # All strand_ratio > 1.0, so all favor forward
        # strand_ratio should be 3/0 but capped
        assert metrics.uniformity.strand_ratio <= 10.0

    def test_params_in_parent_directory(self, tmp_path):
        """Load params.json from parent directory."""
        results_dir = tmp_path / "results"
        results_dir.mkdir()

        # Create step4 in results dir
        step4 = results_dir / "step4_improved_df.csv"
        step4.write_text("sequence\nATCG\n")

        # Create params in parent
        params = tmp_path / "params.json"
        params.write_text(json.dumps({"polymerase": "equiphi29"}))

        metrics = collect_pipeline_metrics(str(results_dir))

        assert metrics.parameters.get("polymerase") == "equiphi29"


class TestOptimizerSummaryIngestion:
    """Test loading and merging of step4_summary.json."""

    def test_load_optimizer_summary_missing_file(self, tmp_path):
        assert _load_optimizer_summary(tmp_path) is None

    def test_load_optimizer_summary_valid(self, tmp_path):
        summary = {"metrics": {"fg_coverage": 0.85, "mean_gap": 5000.0}}
        (tmp_path / "step4_improved_df_summary.json").write_text(json.dumps(summary))
        result = _load_optimizer_summary(tmp_path)
        assert result["metrics"]["fg_coverage"] == 0.85

    def test_load_optimizer_summary_corrupt_json(self, tmp_path):
        (tmp_path / "step4_improved_df_summary.json").write_text("not json{")
        assert _load_optimizer_summary(tmp_path) is None

    def test_coverage_uses_optimizer_data(self, tmp_path):
        """When summary JSON exists, coverage should come from optimizer."""
        csv_content = (
            "sequence,score,fg_freq,bg_freq,tm,gini\nATCGATCG,1.0,0.001,0.00001,35.0,0.2\n"
        )
        (tmp_path / "step4_improved_df.csv").write_text(csv_content)
        params = {"fg_genome": "target.fna", "fg_size": 100000}
        (tmp_path / "params.json").write_text(json.dumps(params))
        summary = {
            "metrics": {
                "fg_coverage": 0.72,
                "mean_gap": 3500.0,
                "max_gap": 12000.0,
                "gap_gini": 0.35,
                "gap_entropy": 2.1,
            }
        }
        (tmp_path / "step4_improved_df_summary.json").write_text(json.dumps(summary))

        metrics = collect_pipeline_metrics(str(tmp_path))
        assert metrics.coverage.overall_coverage == 0.72
        assert metrics.coverage.from_optimizer is True
        assert metrics.coverage.mean_gap == 3500.0
        assert metrics.coverage.max_gap == 12000.0
        assert metrics.coverage.gap_gini == 0.35
        assert metrics.coverage.gap_entropy == 2.1
        assert metrics.coverage.covered_bases == int(0.72 * 100000)

    def test_coverage_without_summary_is_estimated(self, tmp_path):
        """Without summary JSON, coverage is estimated and from_optimizer is False."""
        csv_content = "sequence,fg_freq,bg_freq,tm,gini\nATCGATCG,0.001,0.00001,35.0,0.2\n"
        (tmp_path / "step4_improved_df.csv").write_text(csv_content)
        params = {"fg_size": 100000}
        (tmp_path / "params.json").write_text(json.dumps(params))

        metrics = collect_pipeline_metrics(str(tmp_path))
        assert metrics.coverage.from_optimizer is False
        assert metrics.coverage.overall_coverage > 0  # estimated

    def test_coverage_metrics_new_fields_default(self):
        """New CoverageMetrics fields have correct defaults."""
        cm = CoverageMetrics()
        assert cm.mean_gap == 0.0
        assert cm.max_gap == 0.0
        assert cm.gap_gini == 0.0
        assert cm.gap_entropy == 0.0
        assert cm.from_optimizer is False


# ---------------------------------------------------------------------------
# The report is fed step4_improved_df.csv, whose real schema is the one the
# optimizer writes. It has no `tm`, `fg_count` or `bg_count` column, so the
# parser's `_safe_*(..., 0)` defaults turned every absent column into a
# confident zero: every primer showed Tm 0.0 and specificity 0x, and enrichment
# came out 0x / "Critical" for a set with no background sites at all -- the
# best possible specificity result, graded worst-possible.
# ---------------------------------------------------------------------------

REAL_STEP4_HEADER = (
    "primer,set_index,score,normalized_score,coverage,bg_coverage,selectivity,"
    "mean_gap,max_gap,coverage_uniformity,gap_gini,gap_entropy,dimer_risk_score,"
    "strand_alternation,strand_coverage_ratio,mean_tm,optimizer\n"
)
REAL_STEP4_ROWS = (
    "CGGCGACGATGC,0,1048576.0,0.617,0.5048,0.0,73.55,10331.4,61444,"
    "0.5419,0.5419,4.377,0.6515,0.4741,0.9768,53.68,hybrid\n"
    "CGACGACATCGA,0,1048576.0,0.617,0.5048,0.0,73.55,10331.4,61444,"
    "0.5419,0.5419,4.377,0.6515,0.4741,0.9768,53.68,hybrid\n"
)


@pytest.fixture
def real_step4_results(tmp_path):
    """A results directory in the schema the optimizer actually writes."""
    (tmp_path / "step4_improved_df.csv").write_text(REAL_STEP4_HEADER + REAL_STEP4_ROWS)
    (tmp_path / "params.json").write_text(
        json.dumps(
            {
                "polymerase": "equiphi29",
                "reaction_temp": 42.0,
                "betaine_m": 2.0,
                "fg_genomes": ["mtb.fna"],
                "bg_genomes": ["human_chr21.fna"],
                "fg_seq_lengths": [4411532],
                "bg_seq_lengths": [46709983],
            }
        )
    )
    (tmp_path / "step4_improved_df_summary.json").write_text(
        json.dumps(
            {
                "primers": ["CGGCGACGATGC", "CGACGACATCGA"],
                "score": 1048576.0,
                "metrics": {
                    "fg_coverage": 0.5048,
                    "bg_coverage": 0.0,
                    "total_fg_sites": 427,
                    "total_bg_sites": 0,
                    "selectivity_ratio": 73.55,
                    "selectivity_density": 778.78,
                    "mean_tm": 53.68,
                },
            }
        )
    )
    return tmp_path


class TestRealStep4SchemaIsNotSilentlyZero:
    def test_per_primer_tm_is_not_zero(self, real_step4_results):
        """step4 has no `tm` column, so Tm is derived from the sequence."""
        metrics = collect_pipeline_metrics(str(real_step4_results))

        assert metrics.primers, "no primers parsed"
        assert all(p.tm > 0 for p in metrics.primers)

    def test_per_primer_tm_matches_the_reaction_conditions(self, real_step4_results):
        from neoswga.core.reaction_conditions import ReactionConditions

        conditions = ReactionConditions(temp=42.0, polymerase="equiphi29", betaine_m=2.0)
        metrics = collect_pipeline_metrics(str(real_step4_results))

        for primer in metrics.primers:
            expected = conditions.calculate_effective_tm(primer.sequence)
            assert primer.tm == pytest.approx(expected, abs=0.05)

    def test_zero_background_is_not_reported_as_zero_enrichment(self, real_step4_results):
        """A set with no background sites is the best case, not the worst."""
        metrics = collect_pipeline_metrics(str(real_step4_results))

        assert metrics.specificity.enrichment_ratio > 0

    def test_enrichment_uses_the_measured_density_ratio(self, real_step4_results):
        """`selectivity_density` is the same quantity the report grades:
        target sites per Mb over background sites per Mb."""
        metrics = collect_pipeline_metrics(str(real_step4_results))

        assert metrics.specificity.enrichment_ratio == pytest.approx(778.78, rel=1e-6)


@pytest.fixture
def real_step4_with_step2(real_step4_results):
    """The same run, with the per-primer measurements step2 recorded.

    step4 carries only the set-level result, but the filter step wrote the
    per-primer foreground/background counts for these same primers.
    """
    (real_step4_results / "step2_df.csv").write_text(
        "primer,fg_count,bg_count,gini,gini_bool,ratio,occupancy_ratio\n"
        "CGGCGACGATGC,90,0,0.41,True,90.0,0.21\n"
        "CGACGACATCGA,32,2,0.48,True,16.0,0.19\n"
    )
    return real_step4_results


class TestPerPrimerCountsComeFromTheFilterStep:
    def test_site_counts_are_backfilled_from_step2(self, real_step4_with_step2):
        metrics = collect_pipeline_metrics(str(real_step4_with_step2))

        by_seq = {p.sequence: p for p in metrics.primers}
        assert by_seq["CGGCGACGATGC"].fg_sites == 90
        assert by_seq["CGACGACATCGA"].fg_sites == 32
        assert by_seq["CGACGACATCGA"].bg_sites == 2

    def test_gini_is_backfilled_from_step2(self, real_step4_with_step2):
        metrics = collect_pipeline_metrics(str(real_step4_with_step2))

        by_seq = {p.sequence: p for p in metrics.primers}
        assert by_seq["CGGCGACGATGC"].gini == pytest.approx(0.41)

    def test_per_primer_specificity_is_not_zero(self, real_step4_with_step2):
        """It read 0x for every primer, including ones with no background."""
        metrics = collect_pipeline_metrics(str(real_step4_with_step2))

        assert all(p.specificity > 0 for p in metrics.primers)

    def test_specificity_is_the_density_ratio(self, real_step4_with_step2):
        """Same quantity as enrichment: sites per bp of target over background."""
        metrics = collect_pipeline_metrics(str(real_step4_with_step2))

        by_seq = {p.sequence: p for p in metrics.primers}
        expected = (32 / 4411532) / (2 / 46709983)
        assert by_seq["CGACGACATCGA"].specificity == pytest.approx(expected, rel=1e-6)

    def test_absent_step2_leaves_counts_alone(self, real_step4_results):
        """Nothing to backfill from is not an error."""
        metrics = collect_pipeline_metrics(str(real_step4_results))

        assert all(p.fg_sites == 0 for p in metrics.primers)


class TestReportUsesTheRecordedReactionConditions:
    """The report must correct Tm for the same reaction the optimizer used."""

    @pytest.fixture
    def results_without_betaine(self, real_step4_results):
        """params.json as the user wrote it: no betaine anywhere."""
        params = json.loads((real_step4_results / "params.json").read_text())
        params.pop("betaine_m", None)
        (real_step4_results / "params.json").write_text(json.dumps(params))
        return real_step4_results

    def _manifest(self, path):
        (path / "run_manifest.json").write_text(
            json.dumps(
                {
                    "steps": [
                        {
                            "step": "optimize",
                            "effective_conditions": {
                                "polymerase": "equiphi29",
                                "reaction_temp": 42.0,
                                "mg_conc": 10.0,
                                "betaine_m": 2.0,
                            },
                        }
                    ]
                }
            )
        )

    def test_recorded_betaine_lowers_the_reported_tm(self, results_without_betaine):
        """params.json here carries no betaine; the run used 2 M."""
        without = collect_pipeline_metrics(str(results_without_betaine))
        self._manifest(results_without_betaine)
        with_manifest = collect_pipeline_metrics(str(results_without_betaine))

        assert with_manifest.primers[0].tm < without.primers[0].tm

    def test_reported_tm_matches_the_recorded_reaction(self, results_without_betaine):
        from neoswga.core.reaction_conditions import ReactionConditions

        self._manifest(results_without_betaine)
        metrics = collect_pipeline_metrics(str(results_without_betaine))

        conditions = ReactionConditions(
            temp=42.0, polymerase="equiphi29", mg_conc=10.0, betaine_m=2.0
        )
        for primer in metrics.primers:
            expected = conditions.calculate_effective_tm(primer.sequence)
            assert primer.tm == pytest.approx(expected, abs=0.05)


class TestGenomeInfoReadsTheRealParamsSchema:
    """`fg_genomes` and `fg_seq_lengths` are what params.json carries.

    `_extract_genome_info` looked for `fg_genome` and `fg_size`, neither of
    which the pipeline writes, so the report had no genome at all: no name, no
    size, and `covered_bases` stayed 0 because it is guarded on total_bases.
    """

    def test_target_genome_is_found(self, real_step4_results):
        metrics = collect_pipeline_metrics(str(real_step4_results))

        assert metrics.target_genome is not None
        assert metrics.target_genome.size == 4411532
        assert metrics.target_genome.name == "mtb"

    def test_background_genome_is_found(self, real_step4_results):
        metrics = collect_pipeline_metrics(str(real_step4_results))

        assert metrics.background_genome is not None
        assert metrics.background_genome.size == 46709983

    def test_covered_bases_follows_from_a_known_genome_size(self, real_step4_results):
        metrics = collect_pipeline_metrics(str(real_step4_results))

        assert metrics.coverage.total_bases == 4411532
        assert metrics.coverage.covered_bases == pytest.approx(0.5048 * 4411532, rel=1e-3)

    def test_missing_genome_params_stay_none(self, tmp_path):
        (tmp_path / "step4_improved_df.csv").write_text(REAL_STEP4_HEADER + REAL_STEP4_ROWS)
        (tmp_path / "params.json").write_text(json.dumps({"polymerase": "phi29"}))

        metrics = collect_pipeline_metrics(str(tmp_path))

        assert metrics.target_genome is None
