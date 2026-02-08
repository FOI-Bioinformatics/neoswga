"""
Unit tests for pipeline modules.

Tests:
- pipeline.py: Step validation, prerequisite checking
- improved_pipeline.py: ImprovedPipeline class, PipelineConfig
"""

import pytest
import tempfile
import os
from pathlib import Path

from neoswga.core.pipeline import (
    StepValidationResult,
    StepPrerequisiteError,
    validate_step1_prerequisites,
    validate_step2_prerequisites,
    validate_step3_prerequisites,
    validate_step4_prerequisites,
)
from neoswga.core.improved_pipeline import (
    PipelineConfig,
    ImprovedPipeline,
)


# =============================================================================
# StepValidationResult Tests
# =============================================================================

class TestStepValidationResult:
    """Tests for StepValidationResult dataclass."""

    def test_valid_result(self):
        """Test creating a valid result."""
        result = StepValidationResult(
            valid=True,
            missing_files=[],
            error_message="",
            remediation=""
        )
        assert result.valid is True
        assert result.missing_files == []

    def test_invalid_result(self):
        """Test creating an invalid result with missing files."""
        result = StepValidationResult(
            valid=False,
            missing_files=["file1.fasta", "file2.fasta"],
            error_message="Files not found",
            remediation="Check paths"
        )
        assert result.valid is False
        assert len(result.missing_files) == 2
        assert "file1.fasta" in result.missing_files


class TestStepPrerequisiteError:
    """Tests for StepPrerequisiteError exception."""

    def test_exception_message(self):
        """Test that exception message contains step info."""
        validation = StepValidationResult(
            valid=False,
            missing_files=["missing.fasta"],
            error_message="Test error",
            remediation="Fix it"
        )
        error = StepPrerequisiteError(step=2, validation=validation)

        assert "STEP 2" in str(error)
        assert "Test error" in str(error)
        assert "missing.fasta" in str(error)
        assert "Fix it" in str(error)

    def test_exception_with_many_files(self):
        """Test that exception truncates long file lists."""
        validation = StepValidationResult(
            valid=False,
            missing_files=[f"file{i}.txt" for i in range(10)],
            error_message="Many files missing",
            remediation="Check paths"
        )
        error = StepPrerequisiteError(step=1, validation=validation)

        # Should show first 5 and indicate more
        message = str(error)
        assert "file0.txt" in message
        assert "file4.txt" in message
        assert "5 more" in message


# =============================================================================
# Step Validation Tests
# =============================================================================

class TestValidateStep1Prerequisites:
    """Tests for step 1 prerequisite validation."""

    def test_valid_with_existing_files(self, tmp_path):
        """Test validation passes with existing files."""
        # Create temporary genome files
        fg_genome = tmp_path / "foreground.fasta"
        fg_genome.write_text(">seq1\nACGT\n")

        bg_genome = tmp_path / "background.fasta"
        bg_genome.write_text(">seq1\nACGT\n")

        result = validate_step1_prerequisites(
            data_dir=str(tmp_path),
            fg_genomes=[str(fg_genome)],
            bg_genomes=[str(bg_genome)]
        )

        assert result.valid is True
        assert result.missing_files == []

    def test_invalid_with_missing_fg(self, tmp_path):
        """Test validation fails with missing foreground genome."""
        result = validate_step1_prerequisites(
            data_dir=str(tmp_path),
            fg_genomes=["/nonexistent/genome.fasta"],
            bg_genomes=[]
        )

        assert result.valid is False
        assert "/nonexistent/genome.fasta" in result.missing_files

    def test_valid_without_background(self, tmp_path):
        """Test validation passes without background genomes."""
        fg_genome = tmp_path / "foreground.fasta"
        fg_genome.write_text(">seq1\nACGT\n")

        result = validate_step1_prerequisites(
            data_dir=str(tmp_path),
            fg_genomes=[str(fg_genome)],
            bg_genomes=[]
        )

        assert result.valid is True

    def test_creates_data_dir_if_missing(self, tmp_path):
        """Test that data directory is created if it doesn't exist."""
        new_dir = tmp_path / "new_data_dir"
        fg_genome = tmp_path / "foreground.fasta"
        fg_genome.write_text(">seq1\nACGT\n")

        result = validate_step1_prerequisites(
            data_dir=str(new_dir),
            fg_genomes=[str(fg_genome)],
            bg_genomes=[]
        )

        assert result.valid is True
        assert new_dir.exists()


class TestValidateStep2Prerequisites:
    """Tests for step 2 prerequisite validation."""

    def test_valid_with_kmer_files(self, tmp_path):
        """Test validation passes with k-mer files."""
        # Create mock k-mer file with full path prefix
        kmer_file = tmp_path / "target_6mer_all.txt"
        kmer_file.write_text("ACGTAC\t100\n")

        result = validate_step2_prerequisites(
            data_dir=str(tmp_path),
            fg_prefixes=[str(tmp_path / "target")],
            bg_prefixes=[],
            min_k=6,
            max_k=6
        )

        assert result.valid is True

    def test_invalid_without_kmer_files(self, tmp_path):
        """Test validation fails without k-mer files."""
        result = validate_step2_prerequisites(
            data_dir=str(tmp_path),
            fg_prefixes=[str(tmp_path / "target")],
            bg_prefixes=[],
            min_k=6,
            max_k=8
        )

        assert result.valid is False
        assert "Step 1" in result.remediation


class TestValidateStep3Prerequisites:
    """Tests for step 3 prerequisite validation."""

    def test_valid_with_step2_output(self, tmp_path):
        """Test validation passes with step2 output."""
        step2_file = tmp_path / "step2_df.csv"
        step2_file.write_text("primer,fg_freq,bg_freq\nACGT,100,10\n")

        result = validate_step3_prerequisites(data_dir=str(tmp_path))

        assert result.valid is True

    def test_invalid_without_step2_output(self, tmp_path):
        """Test validation fails without step2 output."""
        result = validate_step3_prerequisites(data_dir=str(tmp_path))

        assert result.valid is False
        assert "step2_df.csv" in str(result.missing_files)


class TestValidateStep4Prerequisites:
    """Tests for step 4 prerequisite validation."""

    def test_valid_with_step3_output(self, tmp_path):
        """Test validation passes with step3 output and positions."""
        step3_file = tmp_path / "step3_df.csv"
        step3_file.write_text("primer,fg_freq,bg_freq,amp_pred\nACGT,100,10,5\n")

        # Position file format is {prefix}_{k}mer_positions.h5
        positions_file = tmp_path / "target_4mer_positions.h5"
        positions_file.write_bytes(b"dummy")  # Just need file to exist

        result = validate_step4_prerequisites(
            data_dir=str(tmp_path),
            fg_prefixes=[str(tmp_path / "target")]
        )

        assert result.valid is True

    def test_invalid_without_step3_output(self, tmp_path):
        """Test validation fails without step3 output."""
        result = validate_step4_prerequisites(
            data_dir=str(tmp_path),
            fg_prefixes=[str(tmp_path / "target")]
        )

        assert result.valid is False
        assert "step3_df.csv" in str(result.missing_files)

    def test_detects_missing_position_files_for_specific_k(self, tmp_path):
        """Test that missing position files for specific primer lengths are detected."""
        # Create step3_df.csv with 8bp and 10bp primers
        step3_file = tmp_path / "step3_df.csv"
        step3_file.write_text(
            "primer,fg_freq,bg_freq,amp_pred\n"
            "ATCGATCG,100,10,5\n"       # 8bp primer
            "GCTAGCTAGC,80,8,4\n"        # 10bp primer
        )

        # Only create position file for k=8, not for k=10
        pos_8 = tmp_path / "target_8mer_positions.h5"
        pos_8.write_bytes(b"dummy")

        result = validate_step4_prerequisites(
            data_dir=str(tmp_path),
            fg_prefixes=[str(tmp_path / "target")]
        )

        assert result.valid is False
        assert any("10mer_positions.h5" in f for f in result.missing_files)

    def test_valid_when_all_required_position_files_exist(self, tmp_path):
        """Test validation passes when position files match primer lengths."""
        step3_file = tmp_path / "step3_df.csv"
        step3_file.write_text(
            "primer,fg_freq,bg_freq,amp_pred\n"
            "ATCGATCG,100,10,5\n"       # 8bp
            "GCTAGCTAGC,80,8,4\n"        # 10bp
        )

        # Create position files for both k=8 and k=10
        (tmp_path / "target_8mer_positions.h5").write_bytes(b"dummy")
        (tmp_path / "target_10mer_positions.h5").write_bytes(b"dummy")

        result = validate_step4_prerequisites(
            data_dir=str(tmp_path),
            fg_prefixes=[str(tmp_path / "target")]
        )

        assert result.valid is True


# =============================================================================
# PipelineConfig Tests
# =============================================================================

class TestPipelineConfig:
    """Tests for PipelineConfig dataclass."""

    def test_default_values(self):
        """Test default configuration values."""
        config = PipelineConfig()

        assert config.gc_tolerance == 0.15
        assert config.reaction_temp == 30.0
        assert config.num_primers == 10
        assert config.optimization_method == 'hybrid'
        assert config.verbose is True

    def test_custom_values(self):
        """Test custom configuration values."""
        config = PipelineConfig(
            gc_tolerance=0.2,
            reaction_temp=42.0,
            num_primers=15,
            optimization_method='greedy'
        )

        assert config.gc_tolerance == 0.2
        assert config.reaction_temp == 42.0
        assert config.num_primers == 15
        assert config.optimization_method == 'greedy'

    def test_thermodynamic_options(self):
        """Test thermodynamic optimization options."""
        config = PipelineConfig(
            tm_weight=0.5,
            dimer_penalty=0.3,
            max_dimer_bp=5
        )

        assert config.tm_weight == 0.5
        assert config.dimer_penalty == 0.3
        assert config.max_dimer_bp == 5


# =============================================================================
# ImprovedPipeline Tests
# =============================================================================

class TestImprovedPipeline:
    """Tests for ImprovedPipeline class."""

    def test_initialization_default_config(self):
        """Test pipeline initialization with default config."""
        pipeline = ImprovedPipeline()

        assert pipeline.config is not None
        assert pipeline.config.num_primers == 10
        assert pipeline.stats == {}

    def test_initialization_custom_config(self):
        """Test pipeline initialization with custom config."""
        config = PipelineConfig(num_primers=8, reaction_temp=37.0)
        pipeline = ImprovedPipeline(config=config)

        assert pipeline.config.num_primers == 8
        assert pipeline.config.reaction_temp == 37.0

    def test_get_genome_lengths(self, tmp_path):
        """Test genome length calculation."""
        # Create a simple FASTA file
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">seq1\nACGTACGT\n>seq2\nGCTA\n")

        pipeline = ImprovedPipeline()
        lengths = pipeline._get_genome_lengths(str(fasta_file))

        # Should return list of sequence lengths
        assert isinstance(lengths, list)
        assert len(lengths) > 0

    def test_get_genome_lengths_missing_file(self, tmp_path):
        """Test genome length calculation with missing file raises error."""
        pipeline = ImprovedPipeline()

        # Should raise FileNotFoundError for missing file
        with pytest.raises(FileNotFoundError):
            pipeline._get_genome_lengths("/nonexistent/file.fasta")


# =============================================================================
# Integration Tests
# =============================================================================

class TestPipelineIntegration:
    """Integration tests for pipeline components."""

    def test_config_affects_pipeline_behavior(self):
        """Test that config changes affect pipeline behavior."""
        # Default config
        default_pipeline = ImprovedPipeline()

        # Verbose disabled config
        quiet_config = PipelineConfig(verbose=False)
        quiet_pipeline = ImprovedPipeline(config=quiet_config)

        assert default_pipeline.config.verbose is True
        assert quiet_pipeline.config.verbose is False

    def test_skip_adaptive_filter_option(self):
        """Test skip_adaptive_filter configuration option."""
        config = PipelineConfig(skip_adaptive_filter=True)
        pipeline = ImprovedPipeline(config=config)

        assert pipeline.config.skip_adaptive_filter is True


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
