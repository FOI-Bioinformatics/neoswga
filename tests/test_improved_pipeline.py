"""Tests for neoswga.core.improved_pipeline module.

Covers PipelineConfig defaults, optimizer routing in design_primers,
genome length parsing, and the parameter migration utility.
"""

import json

import pytest
from unittest.mock import patch, MagicMock, PropertyMock

from neoswga.core.improved_pipeline import (
    ImprovedPipeline,
    PipelineConfig,
    migrate_from_old_pipeline,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_fasta(path, sequences):
    """Write a multi-sequence FASTA file.

    Args:
        sequences: list of (name, seq) tuples.
    """
    with open(path, "w") as f:
        for name, seq in sequences:
            f.write(f">{name}\n{seq}\n")


def _mock_position_cache():
    """Return a mock PositionCache that returns empty positions."""
    cache = MagicMock()
    cache.get_positions.return_value = []
    return cache


def _mock_network_optimizer():
    """Return a mock NetworkOptimizer with optimize_greedy and score_primer_set."""
    opt = MagicMock()
    opt.optimize_greedy.return_value = ["ATCG", "GCTA"]
    opt.score_primer_set.return_value = {
        "target_largest_component": 2,
        "target_connectivity": 0.5,
        "target_amplification": 10.0,
        "background_avg_component": 1.0,
        "background_amplification": 1.0,
        "enrichment": 10.0,
    }
    return opt


# ---------------------------------------------------------------------------
# PipelineConfig
# ---------------------------------------------------------------------------

class TestPipelineConfig:

    def test_defaults(self):
        cfg = PipelineConfig()
        assert cfg.optimization_method == "hybrid"
        assert cfg.num_primers == 10
        assert cfg.reaction_temp == 30.0
        assert cfg.use_position_cache is True
        assert cfg.skip_adaptive_filter is False
        assert cfg.max_extension == 70000

    def test_custom_values(self):
        cfg = PipelineConfig(
            optimization_method="genetic",
            num_primers=6,
            reaction_temp=42.0,
            tm_weight=0.5,
            dimer_penalty=0.3,
        )
        assert cfg.optimization_method == "genetic"
        assert cfg.num_primers == 6
        assert cfg.tm_weight == 0.5

    def test_background_filter_defaults(self):
        cfg = PipelineConfig()
        assert cfg.use_background_filter is True
        assert cfg.background_bloom_path is None
        assert cfg.max_bg_exact_matches == 10
        assert cfg.max_bg_1mm_matches == 100


# ---------------------------------------------------------------------------
# ImprovedPipeline init
# ---------------------------------------------------------------------------

class TestImprovedPipelineInit:

    def test_default_config(self):
        pipeline = ImprovedPipeline()
        assert pipeline.config.optimization_method == "hybrid"
        assert pipeline.stats == {}

    def test_custom_config(self):
        cfg = PipelineConfig(num_primers=5)
        pipeline = ImprovedPipeline(cfg)
        assert pipeline.config.num_primers == 5


# ---------------------------------------------------------------------------
# _get_genome_lengths
# ---------------------------------------------------------------------------

class TestGetGenomeLengths:

    def test_single_sequence(self, tmp_path):
        fasta = tmp_path / "single.fna"
        _write_fasta(fasta, [("seq1", "ATCG" * 100)])
        pipeline = ImprovedPipeline()
        lengths = pipeline._get_genome_lengths(str(fasta))
        assert lengths == [400]

    def test_multiple_sequences(self, tmp_path):
        fasta = tmp_path / "multi.fna"
        _write_fasta(fasta, [("chr1", "A" * 500), ("chr2", "T" * 300)])
        pipeline = ImprovedPipeline()
        lengths = pipeline._get_genome_lengths(str(fasta))
        assert sorted(lengths) == [300, 500]


# ---------------------------------------------------------------------------
# design_primers — optimizer routing
# ---------------------------------------------------------------------------

class TestDesignPrimersRouting:
    """Verify that optimization_method dispatches to the correct optimizer."""

    @pytest.fixture()
    def genome_fasta(self, tmp_path):
        fasta = tmp_path / "genome.fna"
        _write_fasta(fasta, [("seq1", "ATCGATCG" * 50)])
        return str(fasta)

    def _run_with_method(self, method, genome_fasta):
        """Run design_primers with heavy deps mocked, using given method."""
        cfg = PipelineConfig(
            optimization_method=method,
            skip_adaptive_filter=True,
            use_background_filter=False,
            num_primers=2,
        )
        pipeline = ImprovedPipeline(cfg)

        mock_cache = _mock_position_cache()
        mock_net_opt = _mock_network_optimizer()

        patches = {
            "cache": patch(
                "neoswga.core.improved_pipeline.PositionCache",
                return_value=mock_cache,
            ),
            "net_opt": patch(
                "neoswga.core.improved_pipeline.NetworkOptimizer",
                return_value=mock_net_opt,
            ),
            "milp_fb": patch(
                "neoswga.core.improved_pipeline.MILPFallbackOptimizer",
            ),
        }

        with patches["cache"], patches["net_opt"] as net_cls, \
             patches["milp_fb"] as milp_fb_cls:
            # Configure MILP fallback mock
            milp_fb_instance = MagicMock()
            milp_fb_instance.optimize.return_value = ["ATCG", "GCTA"]
            milp_fb_cls.return_value = milp_fb_instance

            result = pipeline.design_primers(
                fg_genome_path=genome_fasta,
                fg_prefixes=["fg_prefix"],
                bg_genome_path=None,
                bg_prefixes=None,
                candidates=["ATCGATCG", "GCTAGCTA", "TTTTAAAA"],
            )

        return result, milp_fb_cls, net_cls

    def test_hybrid_uses_milp_fallback(self, genome_fasta):
        result, milp_fb_cls, net_cls = self._run_with_method("hybrid", genome_fasta)
        milp_fb_cls.assert_called_once()
        assert "primers" in result
        assert "timing" in result

    def test_greedy_uses_network_optimizer(self, genome_fasta):
        result, milp_fb_cls, net_cls = self._run_with_method("greedy", genome_fasta)
        milp_fb_cls.assert_not_called()
        # NetworkOptimizer should be instantiated for greedy/network path
        assert net_cls.called

    def test_skip_adaptive_filter(self, genome_fasta):
        """With skip_adaptive_filter=True, phase 2 should be skipped."""
        cfg = PipelineConfig(
            skip_adaptive_filter=True,
            use_background_filter=False,
            num_primers=2,
        )
        pipeline = ImprovedPipeline(cfg)

        with patch("neoswga.core.improved_pipeline.PositionCache"), \
             patch("neoswga.core.improved_pipeline.NetworkOptimizer",
                   return_value=_mock_network_optimizer()), \
             patch("neoswga.core.improved_pipeline.MILPFallbackOptimizer") as milp_cls, \
             patch("neoswga.core.improved_pipeline.AdaptiveFilterPipeline") as af_cls:

            milp_instance = MagicMock()
            milp_instance.optimize.return_value = ["ATCG"]
            milp_cls.return_value = milp_instance

            result = pipeline.design_primers(
                fg_genome_path=genome_fasta,
                fg_prefixes=["fg"],
                bg_genome_path=None,
                bg_prefixes=None,
                candidates=["ATCG", "GCTA"],
            )

        af_cls.assert_not_called()
        assert result["timing"]["adaptive_filter"] == 0.0

    def test_result_structure(self, genome_fasta):
        """design_primers should return dict with expected keys."""
        result, _, _ = self._run_with_method("hybrid", genome_fasta)
        assert "primers" in result
        assert "evaluation" in result
        assert "timing" in result
        assert "config" in result
        assert "cache_build" in result["timing"]
        assert "optimization" in result["timing"]
        assert "total" in result["timing"]


# ---------------------------------------------------------------------------
# design_primers — genetic algorithm routing
# ---------------------------------------------------------------------------

class TestGeneticAlgorithmRouting:

    def test_genetic_dispatches_to_ga(self, tmp_path):
        fasta = tmp_path / "genome.fna"
        _write_fasta(fasta, [("seq1", "ATCGATCG" * 50)])

        cfg = PipelineConfig(
            optimization_method="genetic",
            skip_adaptive_filter=True,
            use_background_filter=False,
            num_primers=2,
        )
        pipeline = ImprovedPipeline(cfg)

        mock_individual = MagicMock()
        mock_individual.primers = ["ATCG", "GCTA"]

        with patch("neoswga.core.improved_pipeline.PositionCache"), \
             patch("neoswga.core.improved_pipeline.NetworkOptimizer",
                   return_value=_mock_network_optimizer()), \
             patch("neoswga.core.improved_pipeline.PrimerSetGA",
                   create=True) as ga_cls:

            ga_instance = MagicMock()
            ga_instance.evolve.return_value = mock_individual
            ga_cls.return_value = ga_instance

            # Patch the lazy import inside design_primers
            with patch.dict("sys.modules", {
                "neoswga.core.genetic_algorithm": MagicMock(
                    PrimerSetGA=ga_cls,
                    GAConfig=MagicMock(),
                ),
            }):
                result = pipeline.design_primers(
                    fg_genome_path=str(fasta),
                    fg_prefixes=["fg"],
                    bg_genome_path=None,
                    bg_prefixes=None,
                    candidates=["ATCG", "GCTA", "AAAA"],
                )

        assert result["primers"] == ["ATCG", "GCTA"]


# ---------------------------------------------------------------------------
# migrate_from_old_pipeline
# ---------------------------------------------------------------------------

class TestMigrateFromOldPipeline:

    def test_creates_output_file(self, tmp_path):
        old_json = tmp_path / "old_params.json"
        old_json.write_text(json.dumps({
            "fg_genomes": ["target.fna"],
            "bg_genomes": ["host.fna"],
            "fg_prefixes": ["data/target"],
            "bg_prefixes": ["data/host"],
            "data_dir": "data/",
            "max_sets": 8,
        }))

        new_json = tmp_path / "new_params.json"
        migrate_from_old_pipeline(str(old_json), str(new_json))

        assert new_json.exists()
        data = json.loads(new_json.read_text())
        assert data["fg_genomes"] == ["target.fna"]
        assert data["bg_genomes"] == ["host.fna"]
        assert data["config"]["num_primers"] == 8
        assert data["config"]["optimization_method"] == "hybrid"

    def test_preserves_prefixes(self, tmp_path):
        old_json = tmp_path / "old.json"
        old_json.write_text(json.dumps({
            "fg_genomes": ["t.fna"],
            "fg_prefixes": ["data/t"],
            "data_dir": "data/",
        }))

        new_json = tmp_path / "new.json"
        migrate_from_old_pipeline(str(old_json), str(new_json))

        data = json.loads(new_json.read_text())
        assert data["fg_prefixes"] == ["data/t"]
        assert data["data_dir"] == "data/"
