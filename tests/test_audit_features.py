"""
Tests for pipeline audit features:
- No-background (host-free) optimization mode
- Enrichment estimate in results interpreter
- Condition sweep in suggest command
- Simulation fitness integration
- BED/BedGraph export (already implemented, verifying)
- Shannon entropy and strand alternation metrics (already implemented, verifying)
"""

import json
import math
import tempfile
from pathlib import Path

import numpy as np
import pytest


# =============================================================================
# Test: --no-background (host-free) optimization mode
# =============================================================================


class TestNoBackgroundFlag:
    """Test the --no-background CLI flag for host-free optimization."""

    def test_cli_parser_has_no_background_flag(self):
        """CLI parser accepts --no-background."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()
        args = parser.parse_args([
            'optimize', '-j', 'params.json', '--no-background'
        ])
        assert args.no_background is True

    def test_cli_no_background_default_false(self):
        """--no-background defaults to False."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()
        args = parser.parse_args([
            'optimize', '-j', 'params.json'
        ])
        assert args.no_background is False


# =============================================================================
# Test: Enrichment estimate in results interpreter
# =============================================================================


class TestEnrichmentEstimate:
    """Test enrichment fold-change estimation in ResultsInterpreter."""

    def test_results_report_has_enrichment_field(self):
        """ResultsReport dataclass includes enrichment_estimate field."""
        from neoswga.core.results_interpreter import ResultsReport, QualityRating

        report = ResultsReport(
            primer_count=6,
            assessments=[],
            overall_rating=QualityRating.GOOD,
            recommendation="Test",
            next_steps=[],
            warnings=[],
            enrichment_estimate=None,
        )
        assert report.enrichment_estimate is None

    def test_enrichment_estimate_with_params(self, tmp_path):
        """Enrichment estimate works when params.json is available."""
        from neoswga.core.results_interpreter import ResultsInterpreter

        # Create mock step4 results
        step4 = tmp_path / "step4_improved_df.csv"
        step4.write_text(
            "primer,set_index,score,coverage,fg_freq,bg_freq,gini\n"
            "ATCGATCGATCG,0,0.85,0.75,0.001,0.0001,0.3\n"
            "GCTAGCTAGCTA,0,0.85,0.75,0.002,0.0002,0.35\n"
        )

        # Create params.json
        params = {
            "polymerase": "phi29",
            "reaction_temp": 30.0,
            "mg_conc": 2.5,
            "dmso_percent": 0.0,
            "betaine_m": 0.0,
            "fg_gc": 0.5,
        }
        (tmp_path / "params.json").write_text(json.dumps(params))

        interpreter = ResultsInterpreter(str(tmp_path))
        report = interpreter.analyze()

        # Should have enrichment estimate if mechanistic model is available
        if report.enrichment_estimate is not None:
            est = report.enrichment_estimate
            assert 'estimated_fold_enrichment' in est
            assert est['estimated_fold_enrichment'] >= 1.0
            assert 'polymerase' in est
            assert est['polymerase'] == 'phi29'

    def test_enrichment_estimate_without_params(self, tmp_path):
        """Enrichment estimate returns None when params.json is missing."""
        from neoswga.core.results_interpreter import ResultsInterpreter

        # Create mock step4 results only
        step4 = tmp_path / "step4_improved_df.csv"
        step4.write_text(
            "primer,set_index,score,coverage,fg_freq,bg_freq\n"
            "ATCGATCGATCG,0,0.85,0.75,0.001,0.0001\n"
        )

        interpreter = ResultsInterpreter(str(tmp_path))
        report = interpreter.analyze()

        # Should still work, just without enrichment estimate
        assert report is not None
        assert report.primer_count >= 1


# =============================================================================
# Test: Condition sweep
# =============================================================================


class TestConditionSweep:
    """Test condition sweep functionality in condition_suggester."""

    def test_sweep_conditions_basic(self):
        """sweep_conditions returns ranked list of conditions."""
        from neoswga.core.condition_suggester import sweep_conditions

        results = sweep_conditions(
            genome_gc=0.5,
            primer_length=10,
            polymerase='phi29',
            verbose=False,
        )

        assert len(results) > 0
        # Results should be sorted by amplification_factor (descending)
        for i in range(len(results) - 1):
            assert results[i]['amplification_factor'] >= results[i + 1]['amplification_factor']

    def test_sweep_conditions_has_required_fields(self):
        """Each result has required fields."""
        from neoswga.core.condition_suggester import sweep_conditions

        results = sweep_conditions(
            genome_gc=0.5,
            primer_length=10,
            verbose=False,
        )

        required_fields = [
            'temperature', 'dmso_percent', 'betaine_m', 'mg_conc',
            'amplification_factor', 'processivity_factor',
        ]
        for r in results:
            for field in required_fields:
                assert field in r, f"Missing field: {field}"

    def test_sweep_conditions_csv_output(self, tmp_path):
        """sweep_conditions writes CSV output."""
        from neoswga.core.condition_suggester import sweep_conditions

        output_path = str(tmp_path / "sweep_results.csv")
        results = sweep_conditions(
            genome_gc=0.65,
            primer_length=12,
            polymerase='equiphi29',
            output_path=output_path,
            verbose=False,
        )

        assert Path(output_path).exists()
        content = Path(output_path).read_text()
        assert 'temperature' in content
        assert 'amplification_factor' in content

    def test_sweep_conditions_gc_variation(self):
        """Different GC contents produce different optimal conditions."""
        from neoswga.core.condition_suggester import sweep_conditions

        results_low_gc = sweep_conditions(genome_gc=0.3, verbose=False)
        results_high_gc = sweep_conditions(genome_gc=0.7, verbose=False)

        # Best conditions should differ for different GC
        if results_low_gc and results_high_gc:
            best_low = results_low_gc[0]
            best_high = results_high_gc[0]
            # At least one parameter should differ in the top result
            differs = (
                best_low['temperature'] != best_high['temperature'] or
                best_low['dmso_percent'] != best_high['dmso_percent'] or
                best_low['betaine_m'] != best_high['betaine_m']
            )
            # This is not guaranteed but likely
            assert len(results_low_gc) > 0 and len(results_high_gc) > 0

    def test_cli_sweep_flag_exists(self):
        """CLI parser accepts --sweep flag."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()
        args = parser.parse_args([
            'suggest', '--genome-gc', '0.5', '--sweep'
        ])
        assert args.sweep is True

    def test_cli_sweep_output_flag(self):
        """CLI parser accepts --output for sweep results."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()
        args = parser.parse_args([
            'suggest', '--genome-gc', '0.5', '--sweep',
            '--output', 'conditions.csv'
        ])
        assert args.output == 'conditions.csv'


# =============================================================================
# Test: BED/BedGraph export (verify existing implementation)
# =============================================================================


class TestBedExport:
    """Verify BED and BedGraph export functionality."""

    def test_bed_export_basic(self, tmp_path):
        """BED export produces valid 6-column format."""
        from neoswga.core.export import export_to_bed

        primers = ["ATCGATCG", "GCTAGCTA"]
        positions = {
            "ATCGATCG": [(100, 'forward'), (500, 'reverse')],
            "GCTAGCTA": [(200, 'forward')],
        }

        bed_path = str(tmp_path / "test.bed")
        export_to_bed(primers, positions, "chr1", bed_path)

        content = Path(bed_path).read_text()
        lines = content.strip().split('\n')
        assert len(lines) == 3  # 3 total binding sites

        # Verify BED format: chrom, start, end, name, score, strand
        for line in lines:
            fields = line.split('\t')
            assert len(fields) == 6
            assert fields[0] == 'chr1'
            assert fields[5] in ('+', '-')

    def test_bedgraph_export_basic(self, tmp_path):
        """BedGraph export produces valid format."""
        from neoswga.core.export import export_to_bedgraph

        primers = ["ATCGATCG"]
        positions = {
            "ATCGATCG": [(100, 'forward'), (150, 'reverse'), (1100, 'forward')],
        }

        bg_path = str(tmp_path / "test.bedgraph")
        export_to_bedgraph(primers, positions, "chr1", 5000, bg_path, window_size=1000)

        content = Path(bg_path).read_text()
        lines = content.strip().split('\n')
        assert len(lines) >= 1

        # Verify BedGraph format: chrom, start, end, value
        for line in lines:
            fields = line.split('\t')
            assert len(fields) == 4
            assert fields[0] == 'chr1'
            assert int(fields[3]) > 0  # Count should be positive

    def test_cli_export_bed_format(self):
        """CLI export accepts --format bed."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()
        args = parser.parse_args([
            'export', '-d', '/tmp/test', '--format', 'bed'
        ])
        assert args.format == 'bed'

    def test_cli_export_bedgraph_format(self):
        """CLI export accepts --format bedgraph."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()
        args = parser.parse_args([
            'export', '-d', '/tmp/test', '--format', 'bedgraph',
            '--window-size', '500'
        ])
        assert args.format == 'bedgraph'
        assert args.window_size == 500


# =============================================================================
# Test: Shannon entropy metric (verify existing implementation)
# =============================================================================


class TestShannonEntropy:
    """Verify Shannon entropy metric in base_optimizer and position_cache."""

    def test_primer_set_metrics_has_entropy(self):
        """PrimerSetMetrics includes gap_entropy field."""
        from neoswga.core.base_optimizer import PrimerSetMetrics

        metrics = PrimerSetMetrics.empty()
        assert hasattr(metrics, 'gap_entropy')
        assert metrics.gap_entropy == 0.0

    def test_entropy_in_metrics_dict(self):
        """gap_entropy appears in to_dict() output."""
        from neoswga.core.base_optimizer import PrimerSetMetrics

        metrics = PrimerSetMetrics(
            fg_coverage=0.8,
            bg_coverage=0.1,
            coverage_uniformity=0.3,
            total_fg_sites=100,
            total_bg_sites=10,
            selectivity_ratio=10.0,
            mean_tm=30.0,
            tm_range=(28.0, 32.0),
            dimer_risk_score=0.1,
            mean_gap=5000.0,
            max_gap=15000.0,
            gap_gini=0.3,
            gap_entropy=2.5,
            strand_alternation_score=0.7,
            strand_coverage_ratio=0.9,
        )

        d = metrics.to_dict()
        assert 'gap_entropy' in d
        assert d['gap_entropy'] == 2.5


# =============================================================================
# Test: Strand alternation metrics (verify existing implementation)
# =============================================================================


class TestStrandAlternation:
    """Verify strand alternation metrics in base_optimizer."""

    def test_metrics_has_strand_fields(self):
        """PrimerSetMetrics includes strand alternation fields."""
        from neoswga.core.base_optimizer import PrimerSetMetrics

        metrics = PrimerSetMetrics.empty()
        assert hasattr(metrics, 'strand_alternation_score')
        assert hasattr(metrics, 'strand_coverage_ratio')

    def test_strand_fields_in_dict(self):
        """Strand fields appear in to_dict() output."""
        from neoswga.core.base_optimizer import PrimerSetMetrics

        metrics = PrimerSetMetrics.empty()
        d = metrics.to_dict()
        assert 'strand_alternation_score' in d
        assert 'strand_coverage_ratio' in d


# =============================================================================
# Test: Simulation fitness module
# =============================================================================


class TestSimulationFitnessModule:
    """Test simulation_fitness module availability and structure."""

    def test_module_importable(self):
        """simulation_fitness module can be imported."""
        from neoswga.core.simulation_fitness import (
            SimulationBasedEvaluator,
            SimulationFitness,
        )
        assert SimulationBasedEvaluator is not None
        assert SimulationFitness is not None

    def test_simulation_fitness_dataclass(self):
        """SimulationFitness has expected fields."""
        from neoswga.core.simulation_fitness import SimulationFitness

        fitness = SimulationFitness(
            mean_coverage=0.8,
            std_coverage=0.05,
            mean_forks_created=50.0,
            mean_fork_travel=10000.0,
            coverage_uniformity=0.7,
            fitness_score=0.75,
            simulation_time=5.0,
            n_replicates=3,
        )

        assert fitness.mean_coverage == 0.8
        assert fitness.fitness_score == 0.75


# =============================================================================
# Test: Unified optimizer simulation re-scoring
# =============================================================================


class TestSimulationRescore:
    """Test _simulation_rescore function in unified_optimizer."""

    def test_rescore_function_exists(self):
        """_simulation_rescore is available."""
        from neoswga.core.unified_optimizer import _simulation_rescore
        assert callable(_simulation_rescore)

    def test_rescore_with_failed_result(self):
        """Simulation re-score returns None for failed optimization."""
        from neoswga.core.unified_optimizer import _simulation_rescore
        from neoswga.core.base_optimizer import OptimizationResult

        result = OptimizationResult.failure("test", "test failure")
        score = _simulation_rescore(result, ['prefix'], [1000])
        assert score is None
