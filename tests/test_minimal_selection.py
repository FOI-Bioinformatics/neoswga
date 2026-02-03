"""
Tests for Phase 2: Minimal Primer Selection.

Verifies that the --minimize-primers option properly reduces primer count
while maintaining target coverage.
"""

import pytest
import numpy as np
from unittest.mock import Mock


class TestMinimalPrimerSelector:
    """Test MinimalPrimerSelector functionality."""

    def test_select_minimal_set_reduces_count(self):
        """Minimal selection should reduce primer count."""
        from neoswga.core.minimal_primer_selector import MinimalPrimerSelector, Primer

        genome_length = 100000
        selector = MinimalPrimerSelector(genome_length=genome_length)

        # Create primers with overlapping coverage
        primers = []
        for i in range(10):
            # Each primer covers 10k positions starting at i*8k
            # This creates overlapping coverage
            start = i * 8000
            positions = list(range(start, min(start + 10000, genome_length)))
            primers.append(Primer(
                sequence=f"PRIMER{i:02d}",
                binding_sites=positions,
                quality_score=0.8,
                gc_content=0.5
            ))

        # Select minimal set for 50% coverage
        result = selector.select_minimal_set(
            primers=primers,
            target_coverage=0.50,
            max_primers=10
        )

        # Should need fewer than all 10 primers
        assert len(result.selected_primers) < 10
        # Should achieve target coverage
        assert result.coverage >= 0.50

    def test_select_minimal_set_respects_target_coverage(self):
        """Selection should stop once target coverage is reached."""
        from neoswga.core.minimal_primer_selector import MinimalPrimerSelector, Primer

        genome_length = 100000
        selector = MinimalPrimerSelector(genome_length=genome_length)

        # Create non-overlapping primers (each covers unique region)
        primers = []
        for i in range(10):
            start = i * 10000
            positions = list(range(start, start + 10000))
            primers.append(Primer(
                sequence=f"PRIMER{i:02d}",
                binding_sites=positions,
                quality_score=0.8,
                gc_content=0.5
            ))

        # Request 30% coverage - should need 3 primers
        result = selector.select_minimal_set(
            primers=primers,
            target_coverage=0.30,
            max_primers=10
        )

        # Should select approximately 3 primers for 30% coverage
        assert len(result.selected_primers) <= 5  # Allow some margin
        assert result.coverage >= 0.30

    def test_coverage_uniformity_calculation(self):
        """Coverage uniformity should be calculated correctly."""
        from neoswga.core.minimal_primer_selector import MinimalPrimerSelector, Primer

        genome_length = 100000
        selector = MinimalPrimerSelector(genome_length=genome_length)

        # Create primers with even coverage
        primers = []
        for i in range(5):
            start = i * 20000
            positions = list(range(start, start + 15000))
            primers.append(Primer(
                sequence=f"PRIMER{i:02d}",
                binding_sites=positions,
                quality_score=0.8,
                gc_content=0.5
            ))

        result = selector.select_minimal_set(
            primers=primers,
            target_coverage=0.70,
            max_primers=5
        )

        # Uniformity should be calculated
        assert result.coverage_uniformity >= 0
        assert result.coverage_uniformity <= 1

    def test_gap_analysis(self):
        """Gap analysis should identify uncovered regions."""
        from neoswga.core.minimal_primer_selector import MinimalPrimerSelector, Primer

        genome_length = 100000
        selector = MinimalPrimerSelector(genome_length=genome_length)

        # Create primers that leave a gap in the middle
        primers = [
            Primer(
                sequence="PRIMER01",
                binding_sites=list(range(0, 30000)),
                quality_score=0.8,
                gc_content=0.5
            ),
            Primer(
                sequence="PRIMER02",
                binding_sites=list(range(70000, 100000)),
                quality_score=0.8,
                gc_content=0.5
            )
        ]

        result = selector.select_minimal_set(
            primers=primers,
            target_coverage=0.50,
            max_primers=2
        )

        # Should have identified a gap
        assert len(result.gaps) > 0


class TestMinimalSelectionCLI:
    """Test CLI options for minimal selection."""

    def test_cli_has_minimize_primers_option(self):
        """CLI should have --minimize-primers option."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()
        args = parser.parse_args(['optimize', '-j', 'test.json', '--minimize-primers'])

        assert hasattr(args, 'minimize_primers')
        assert args.minimize_primers == True

    def test_cli_has_target_coverage_option(self):
        """CLI should have --target-coverage option."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()
        args = parser.parse_args(['optimize', '-j', 'test.json', '--target-coverage', '0.80'])

        assert hasattr(args, 'target_coverage')
        assert args.target_coverage == 0.80

    def test_cli_target_coverage_default(self):
        """Default target-coverage should be 0.70."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()
        args = parser.parse_args(['optimize', '-j', 'test.json'])

        assert args.target_coverage == 0.70


class TestPipelineIntegrationMinimal:
    """Test pipeline integration with minimal selection."""

    def test_optimize_step4_accepts_minimize_params(self):
        """optimize_step4 should accept minimize_primers and target_coverage."""
        import inspect
        from neoswga.core.unified_optimizer import optimize_step4

        sig = inspect.signature(optimize_step4)
        params = sig.parameters

        assert 'minimize_primers' in params
        assert 'target_coverage' in params

        # Check defaults
        assert params['minimize_primers'].default == False
        assert params['target_coverage'].default == 0.70


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
