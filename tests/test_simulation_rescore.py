"""Smoke tests for _simulation_rescore in unified_optimizer."""

from unittest.mock import patch, MagicMock
import pytest

from neoswga.core.base_optimizer import (
    OptimizationResult,
    OptimizationStatus,
    PrimerSetMetrics,
)
from neoswga.core.unified_optimizer import _simulation_rescore


def _make_success_result(primers=("ATCGATCG", "GCTAGCTA")):
    """Helper to build a successful OptimizationResult."""
    return OptimizationResult(
        primers=primers,
        score=1.0,
        status=OptimizationStatus.SUCCESS,
        metrics=PrimerSetMetrics.empty(),
        iterations=1,
        optimizer_name="test",
    )


class TestSimulationRescoreFailurePaths:
    """Tests for early-return None cases."""

    def test_returns_none_for_failed_result(self):
        result = OptimizationResult.failure("test", "failed")
        out = _simulation_rescore(result, ["prefix"], [1000])
        assert out is None

    def test_returns_none_for_empty_primers(self):
        result = _make_success_result(primers=tuple())
        out = _simulation_rescore(result, ["prefix"], [1000])
        assert out is None

    def test_returns_none_when_import_fails(self):
        """Simulate ImportError by making the module entry None in sys.modules."""
        result = _make_success_result()
        with patch.dict(
            "sys.modules",
            {"neoswga.core.simulation_fitness": None},
        ):
            out = _simulation_rescore(result, ["prefix"], [1000], verbose=False)
            assert out is None


class TestSimulationRescoreHappyPath:
    """Test the happy path with mocked dependencies."""

    def test_returns_expected_keys(self):
        """Mock GenomeLoader, PositionCache, and SimulationBasedEvaluator."""
        # Build a fake fitness result returned by the evaluator
        fake_fitness = MagicMock()
        fake_fitness.mean_coverage = 0.85
        fake_fitness.coverage_uniformity = 0.72
        fake_fitness.fitness_score = 0.80
        fake_fitness.mean_forks_created = 42.0
        fake_fitness.simulation_time = 1800.0

        mock_evaluator_cls = MagicMock()
        mock_evaluator_cls.return_value.evaluate.return_value = fake_fitness

        mock_loader_cls = MagicMock()
        mock_loader_cls.return_value.load_genome.return_value = "ATCG" * 250

        mock_cache_cls = MagicMock()

        result = _make_success_result()

        # Patch the classes in the source modules so that the local imports
        # inside _simulation_rescore pick up the mocks.
        with patch(
            "neoswga.core.simulation_fitness.SimulationBasedEvaluator",
            mock_evaluator_cls,
        ), patch(
            "neoswga.core.position_cache.PositionCache",
            mock_cache_cls,
        ), patch(
            "neoswga.core.genome_io.GenomeLoader",
            mock_loader_cls,
        ), patch(
            "neoswga.core.unified_optimizer.parameter",
        ) as mock_parameter:
            mock_parameter.fg_genomes = ["/fake/genome.fna"]

            out = _simulation_rescore(
                result,
                ["prefix"],
                [1000],
                simulation_time=1800.0,
                verbose=False,
            )

        assert out is not None, "Expected a dict but got None"
        expected_keys = {
            "simulation_coverage",
            "simulation_uniformity",
            "simulation_fitness",
            "simulation_forks",
            "simulation_time",
        }
        assert set(out.keys()) == expected_keys

        assert out["simulation_coverage"] == 0.85
        assert out["simulation_uniformity"] == 0.72
        assert out["simulation_fitness"] == 0.80
        assert out["simulation_forks"] == 42.0
        assert out["simulation_time"] == 1800.0

    def test_returns_none_when_genome_not_available(self):
        """When fg_genomes is empty, genome_seq stays None and function returns None."""
        result = _make_success_result()

        with patch(
            "neoswga.core.simulation_fitness.SimulationBasedEvaluator",
            MagicMock(),
        ), patch(
            "neoswga.core.position_cache.PositionCache",
            MagicMock(),
        ), patch(
            "neoswga.core.genome_io.GenomeLoader",
            MagicMock(),
        ), patch(
            "neoswga.core.unified_optimizer.parameter",
        ) as mock_parameter:
            mock_parameter.fg_genomes = []

            out = _simulation_rescore(
                result,
                ["prefix"],
                [1000],
                verbose=False,
            )

        assert out is None
