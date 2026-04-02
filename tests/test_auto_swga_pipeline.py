"""Tests for neoswga.core.auto_swga_pipeline module.

Exercises the automated pipeline orchestration: genome analysis,
GC-adaptive strategy selection, and result serialization.

Note: The auto-pipeline's run() method raises NotImplementedError because
candidate generation requires Jellyfish. These tests verify that the error
is raised with a helpful message directing users to the standard pipeline.
"""

import json
import random

import pytest

from neoswga.core.auto_swga_pipeline import (
    AutoSWGAPipeline,
    PipelineConfig,
    PipelineResult,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_fasta(path, name, sequence):
    """Write a single-sequence FASTA file."""
    with open(path, "w") as f:
        f.write(f">{name}\n{sequence}\n")


def _make_genome(gc_fraction, length=2000):
    """Generate a deterministic genome with approximate GC content."""
    rng = random.Random(42)
    bases_gc = ["G", "C"]
    bases_at = ["A", "T"]
    seq = []
    for _ in range(length):
        if rng.random() < gc_fraction:
            seq.append(rng.choice(bases_gc))
        else:
            seq.append(rng.choice(bases_at))
    return "".join(seq)


# ---------------------------------------------------------------------------
# PipelineConfig
# ---------------------------------------------------------------------------

class TestPipelineConfig:

    def test_defaults(self):
        cfg = PipelineConfig(target_genome="t.fna")
        assert cfg.kmer_min == 8
        assert cfg.kmer_max == 15
        assert cfg.thermodynamic_filter is True
        assert cfg.target_primer_count is None
        assert cfg.force_polymerase is None

    def test_custom_values(self):
        cfg = PipelineConfig(
            target_genome="t.fna",
            kmer_min=12,
            kmer_max=18,
            target_primer_count=6,
            force_polymerase="equiphi29",
        )
        assert cfg.kmer_min == 12
        assert cfg.force_polymerase == "equiphi29"


# ---------------------------------------------------------------------------
# PipelineResult
# ---------------------------------------------------------------------------

class TestPipelineResult:

    @pytest.fixture()
    def result(self):
        return PipelineResult(
            primers=["ATCGATCG", "GCTAGCTA"],
            primer_count=2,
            target_genome_name="test",
            target_genome_size=5000,
            target_gc_content=0.45,
            genome_classification="balanced",
            polymerase="phi29",
            reaction_temp=30.0,
            betaine_concentration=0.0,
            kmer_range=(8, 12),
            coverage=0.55,
            connectivity=0.8,
            predicted_amplification=40.0,
            total_candidates=100,
            thermodynamic_filtered=50,
            optimization_method="hybrid",
            total_runtime=1.5,
            protocol="mock protocol",
        )

    def test_to_dict(self, result):
        d = result.to_dict()
        assert d["primers"] == ["ATCGATCG", "GCTAGCTA"]
        assert d["primer_count"] == 2
        assert d["target_gc_content"] == 0.45

    def test_str_contains_key_sections(self, result):
        s = str(result)
        assert "Genome Information" in s
        assert "Strategy" in s
        assert "Selected Primers" in s
        assert "ATCGATCG" in s

    def test_str_with_simulation(self, result):
        result.simulation_coverage = 0.60
        result.simulation_fitness = 0.85
        s = str(result)
        assert "Simulation Validation" in s

    def test_to_dict_roundtrip(self, result):
        d = result.to_dict()
        serialized = json.dumps(d)
        restored = json.loads(serialized)
        assert restored["primers"] == result.primers
        assert restored["polymerase"] == result.polymerase


# ---------------------------------------------------------------------------
# AutoSWGAPipeline.__init__ validation
# ---------------------------------------------------------------------------

class TestAutoPipelineInit:

    def test_missing_target_raises(self):
        with pytest.raises(FileNotFoundError):
            AutoSWGAPipeline("/nonexistent/genome.fna")

    def test_missing_background_raises(self, tmp_path):
        fg = tmp_path / "target.fna"
        _write_fasta(fg, "seq1", "ATCG" * 50)
        with pytest.raises(FileNotFoundError):
            AutoSWGAPipeline(str(fg), background_genome="/nonexistent/bg.fna")

    def test_valid_init(self, tmp_path):
        fg = tmp_path / "target.fna"
        _write_fasta(fg, "seq1", "ATCG" * 50)
        pipeline = AutoSWGAPipeline(str(fg), output_dir=str(tmp_path / "out"))
        assert pipeline.target_genome == fg


# ---------------------------------------------------------------------------
# Full pipeline run raises NotImplementedError
# ---------------------------------------------------------------------------

class TestAutoPipelineRun:

    def _make_pipeline(self, tmp_path, gc_fraction):
        fg = tmp_path / "target.fna"
        _write_fasta(fg, "test_genome", _make_genome(gc_fraction, length=3000))
        return AutoSWGAPipeline(str(fg), output_dir=str(tmp_path / "out"))

    def test_run_raises_not_implemented(self, tmp_path):
        """run() raises NotImplementedError directing users to standard pipeline."""
        pipeline = self._make_pipeline(tmp_path, 0.50)
        with pytest.raises(NotImplementedError, match="standard pipeline"):
            pipeline.run(target_primer_count=5, verbose=False)

    def test_error_message_includes_commands(self, tmp_path):
        """Error message includes the four standard pipeline commands."""
        pipeline = self._make_pipeline(tmp_path, 0.50)
        with pytest.raises(NotImplementedError) as exc_info:
            pipeline.run(target_primer_count=5, verbose=False)
        msg = str(exc_info.value)
        assert "count-kmers" in msg
        assert "filter" in msg
        assert "score" in msg
        assert "optimize" in msg

    def test_error_message_includes_gc_info(self, tmp_path):
        """Error message includes genome GC content and classification."""
        pipeline = self._make_pipeline(tmp_path, 0.70)
        with pytest.raises(NotImplementedError) as exc_info:
            pipeline.run(target_primer_count=5, verbose=False)
        msg = str(exc_info.value)
        assert "GC content" in msg or "%" in msg
