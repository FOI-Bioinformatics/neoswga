"""Tests for neoswga.core.auto_swga_pipeline module.

Exercises the automated pipeline orchestration: genome analysis,
GC-adaptive strategy selection, thermodynamic filtering, protocol
generation, and result serialization.
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
# _generate_mock_candidates
# ---------------------------------------------------------------------------

class TestGenerateMockCandidates:

    @pytest.fixture()
    def pipeline(self, tmp_path):
        fg = tmp_path / "target.fna"
        _write_fasta(fg, "seq1", _make_genome(0.50))
        return AutoSWGAPipeline(str(fg), output_dir=str(tmp_path / "out"))

    def test_correct_count(self, pipeline):
        candidates = pipeline._generate_mock_candidates(0.50, (8, 12), count=50)
        assert len(candidates) == 50

    def test_kmer_lengths_in_range(self, pipeline):
        candidates = pipeline._generate_mock_candidates(0.50, (10, 15), count=200)
        for c in candidates:
            assert 10 <= len(c) <= 15

    def test_only_valid_bases(self, pipeline):
        candidates = pipeline._generate_mock_candidates(0.50, (8, 12), count=100)
        valid = set("ATCG")
        for c in candidates:
            assert set(c).issubset(valid), f"Invalid bases in {c}"


# ---------------------------------------------------------------------------
# Full pipeline run
# ---------------------------------------------------------------------------

class TestAutoPipelineRun:

    def _make_pipeline(self, tmp_path, gc_fraction):
        fg = tmp_path / "target.fna"
        _write_fasta(fg, "test_genome", _make_genome(gc_fraction, length=3000))
        return AutoSWGAPipeline(str(fg), output_dir=str(tmp_path / "out"))

    def test_balanced_genome(self, tmp_path):
        pipeline = self._make_pipeline(tmp_path, 0.50)
        result = pipeline.run(target_primer_count=5, verbose=False)
        assert isinstance(result, PipelineResult)
        assert result.primer_count <= 5
        assert len(result.primers) == result.primer_count
        assert result.target_genome_name == "target"
        assert result.total_runtime > 0

    def test_gc_rich_genome(self, tmp_path):
        pipeline = self._make_pipeline(tmp_path, 0.70)
        result = pipeline.run(target_primer_count=5, verbose=False)
        assert result.target_gc_content > 0.60
        assert result.genome_classification == "gc_rich"

    def test_at_rich_genome(self, tmp_path):
        pipeline = self._make_pipeline(tmp_path, 0.25)
        result = pipeline.run(target_primer_count=5, verbose=False)
        assert result.target_gc_content < 0.35
        assert result.genome_classification == "at_rich"

    def test_force_polymerase(self, tmp_path):
        pipeline = self._make_pipeline(tmp_path, 0.50)
        result = pipeline.run(
            target_primer_count=5,
            force_polymerase="equiphi29",
            verbose=False,
        )
        assert result.polymerase == "equiphi29"

    def test_auto_primer_count_small_genome(self, tmp_path):
        """Genome < 2 Mbp should get base count of 8."""
        pipeline = self._make_pipeline(tmp_path, 0.50)
        result = pipeline.run(verbose=False)
        # 3000 bp genome < 2 Mbp → base_count=8, multiplied by primer_count_multiplier
        assert result.primer_count >= 4  # at least reasonable


# ---------------------------------------------------------------------------
# Protocol generation
# ---------------------------------------------------------------------------

class TestGenerateProtocol:

    def _run_pipeline(self, tmp_path, gc_fraction, **kwargs):
        fg = tmp_path / "target.fna"
        _write_fasta(fg, "test_genome", _make_genome(gc_fraction, length=2000))
        pipeline = AutoSWGAPipeline(str(fg), output_dir=str(tmp_path / "out"))
        return pipeline.run(target_primer_count=3, verbose=False, **kwargs)

    def test_protocol_contains_primers(self, tmp_path):
        result = self._run_pipeline(tmp_path, 0.50)
        for primer in result.primers:
            assert primer in result.protocol

    def test_protocol_contains_temperature(self, tmp_path):
        result = self._run_pipeline(tmp_path, 0.50)
        assert f"{result.reaction_temp}" in result.protocol

    def test_gc_rich_note(self, tmp_path):
        result = self._run_pipeline(tmp_path, 0.72)
        assert "GC-rich" in result.protocol

    def test_at_rich_note(self, tmp_path):
        result = self._run_pipeline(tmp_path, 0.22)
        assert "AT-rich" in result.protocol


# ---------------------------------------------------------------------------
# save_results
# ---------------------------------------------------------------------------

class TestSaveResults:

    def test_saves_all_files(self, tmp_path):
        fg = tmp_path / "target.fna"
        _write_fasta(fg, "test_genome", _make_genome(0.50, length=2000))
        pipeline = AutoSWGAPipeline(str(fg), output_dir=str(tmp_path / "out"))
        result = pipeline.run(target_primer_count=3, verbose=False)

        out_dir = tmp_path / "save_output"
        pipeline.save_results(result, output_dir=str(out_dir))

        assert (out_dir / "primers.fasta").exists()
        assert (out_dir / "results.json").exists()
        assert (out_dir / "protocol.txt").exists()
        assert (out_dir / "summary.txt").exists()

    def test_primers_fasta_format(self, tmp_path):
        fg = tmp_path / "target.fna"
        _write_fasta(fg, "test_genome", _make_genome(0.50, length=2000))
        pipeline = AutoSWGAPipeline(str(fg), output_dir=str(tmp_path / "out"))
        result = pipeline.run(target_primer_count=3, verbose=False)

        out_dir = tmp_path / "save_output"
        pipeline.save_results(result, output_dir=str(out_dir))

        content = (out_dir / "primers.fasta").read_text()
        assert content.startswith(">primer_1\n")
        assert content.count(">") == result.primer_count

    def test_results_json_valid(self, tmp_path):
        fg = tmp_path / "target.fna"
        _write_fasta(fg, "test_genome", _make_genome(0.50, length=2000))
        pipeline = AutoSWGAPipeline(str(fg), output_dir=str(tmp_path / "out"))
        result = pipeline.run(target_primer_count=3, verbose=False)

        out_dir = tmp_path / "save_output"
        pipeline.save_results(result, output_dir=str(out_dir))

        data = json.loads((out_dir / "results.json").read_text())
        assert data["primer_count"] == result.primer_count
        assert data["polymerase"] == result.polymerase
