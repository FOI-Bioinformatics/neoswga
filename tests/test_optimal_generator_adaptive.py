#!/usr/bin/env python3
"""
Test adaptive QA integration in OptimalOligoGenerator.

Verifies that:
1. Genome analysis detects GC composition correctly
2. Adaptive QA recommendations are generated
3. Reaction conditions include genome_gc parameter
4. AT-rich and GC-rich genomes trigger adaptive QA
5. Balanced genomes do not trigger adaptive QA
"""

import unittest
import tempfile
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from neoswga.core.optimal_oligo_generator import OptimalOligoGenerator


class TestOptimalGeneratorAdaptiveQA(unittest.TestCase):
    """Test adaptive QA integration in OptimalOligoGenerator."""

    def setUp(self):
        """Create temporary test genomes."""
        self.temp_dir = tempfile.mkdtemp()

    def create_test_genome(self, name: str, length: int, gc_content: float) -> Path:
        """
        Create a test genome FASTA file with specified GC content.

        Args:
            name: Genome name
            length: Sequence length (bp)
            gc_content: Target GC content (0-1)

        Returns:
            Path to created FASTA file
        """
        # Calculate base counts
        gc_count = int(length * gc_content)
        at_count = length - gc_count
        g_count = gc_count // 2
        c_count = gc_count - g_count
        a_count = at_count // 2
        t_count = at_count - a_count

        # Create sequence
        sequence = 'G' * g_count + 'C' * c_count + 'A' * a_count + 'T' * t_count

        # Shuffle to avoid patterns
        import random
        sequence_list = list(sequence)
        random.shuffle(sequence_list)
        sequence = ''.join(sequence_list)

        # Create SeqRecord
        record = SeqRecord(
            Seq(sequence),
            id=name,
            description=f"Test genome {gc_content:.0%} GC"
        )

        # Write to file
        filepath = Path(self.temp_dir) / f"{name}.fasta"
        SeqIO.write([record], filepath, 'fasta')

        return filepath

    def test_at_rich_genome_triggers_adaptive_qa(self):
        """AT-rich genome (32% GC) should trigger adaptive QA."""
        # Create Francisella-like genome
        francisella = self.create_test_genome('francisella', 50000, 0.32)

        # Create generator
        generator = OptimalOligoGenerator(
            target_genome=str(francisella),
            output_dir=Path(self.temp_dir) / 'output'
        )

        # Check genome characteristics
        chars = generator.target_chars
        self.assertAlmostEqual(chars.gc_content, 0.32, places=2)
        self.assertEqual(chars.gc_class, 'at_rich')
        self.assertTrue(chars.use_adaptive_qa)
        self.assertIn('AT-rich', chars.adaptive_qa_reason)
        self.assertIn('coverage improvement', chars.expected_improvement)

        # Generate conditions and verify genome_gc is included
        polymerase = generator.select_polymerase()
        conditions = generator.generate_reaction_conditions(polymerase)

        self.assertIn('genome_gc', conditions)
        self.assertAlmostEqual(conditions['genome_gc'], 0.32, places=2)
        self.assertIn('use_adaptive_qa', conditions)
        self.assertTrue(conditions['use_adaptive_qa'])

    def test_gc_rich_genome_triggers_adaptive_qa(self):
        """GC-rich genome (67% GC) should trigger adaptive QA."""
        # Create Burkholderia-like genome
        burkholderia = self.create_test_genome('burkholderia', 50000, 0.67)

        # Create generator
        generator = OptimalOligoGenerator(
            target_genome=str(burkholderia),
            output_dir=Path(self.temp_dir) / 'output'
        )

        # Check genome characteristics
        chars = generator.target_chars
        self.assertAlmostEqual(chars.gc_content, 0.67, places=2)
        self.assertEqual(chars.gc_class, 'gc_rich')
        self.assertTrue(chars.use_adaptive_qa)
        self.assertIn('GC-rich', chars.adaptive_qa_reason)

        # Verify conditions include genome_gc
        polymerase = generator.select_polymerase()
        conditions = generator.generate_reaction_conditions(polymerase)

        self.assertIn('genome_gc', conditions)
        self.assertAlmostEqual(conditions['genome_gc'], 0.67, places=2)
        self.assertTrue(conditions['use_adaptive_qa'])

    def test_balanced_genome_no_adaptive_qa(self):
        """Balanced genome (50% GC) should NOT trigger adaptive QA."""
        # Create E. coli-like genome
        ecoli = self.create_test_genome('ecoli', 50000, 0.50)

        # Create generator
        generator = OptimalOligoGenerator(
            target_genome=str(ecoli),
            output_dir=Path(self.temp_dir) / 'output'
        )

        # Check genome characteristics
        chars = generator.target_chars
        self.assertAlmostEqual(chars.gc_content, 0.50, places=2)
        self.assertEqual(chars.gc_class, 'balanced')
        self.assertFalse(chars.use_adaptive_qa)
        self.assertIn('Balanced', chars.adaptive_qa_reason)

        # Verify conditions include genome_gc (even when not using adaptive QA)
        polymerase = generator.select_polymerase()
        conditions = generator.generate_reaction_conditions(polymerase)

        self.assertIn('genome_gc', conditions)
        self.assertAlmostEqual(conditions['genome_gc'], 0.50, places=2)
        self.assertFalse(conditions['use_adaptive_qa'])

    def test_extreme_at_genome_triggers_critical_adaptive_qa(self):
        """Extreme AT-rich genome (19% GC) should trigger CRITICAL adaptive QA."""
        # Create Plasmodium-like genome
        plasmodium = self.create_test_genome('plasmodium', 50000, 0.19)

        # Create generator
        generator = OptimalOligoGenerator(
            target_genome=str(plasmodium),
            output_dir=Path(self.temp_dir) / 'output'
        )

        # Check genome characteristics
        chars = generator.target_chars
        self.assertAlmostEqual(chars.gc_content, 0.19, places=2)
        self.assertEqual(chars.gc_class, 'extreme_at')
        self.assertTrue(chars.use_adaptive_qa)
        self.assertIn('Extreme AT-rich', chars.adaptive_qa_reason)
        self.assertIn('CRITICAL', chars.adaptive_qa_reason)

        # Verify conditions
        polymerase = generator.select_polymerase()
        conditions = generator.generate_reaction_conditions(polymerase)

        self.assertAlmostEqual(conditions['genome_gc'], 0.19, places=2)
        self.assertTrue(conditions['use_adaptive_qa'])

    def test_extreme_gc_genome_triggers_critical_adaptive_qa(self):
        """Extreme GC-rich genome (72% GC) should trigger CRITICAL adaptive QA."""
        # Create Streptomyces-like genome
        streptomyces = self.create_test_genome('streptomyces', 50000, 0.72)

        # Create generator
        generator = OptimalOligoGenerator(
            target_genome=str(streptomyces),
            output_dir=Path(self.temp_dir) / 'output'
        )

        # Check genome characteristics
        chars = generator.target_chars
        self.assertAlmostEqual(chars.gc_content, 0.72, places=2)
        self.assertEqual(chars.gc_class, 'extreme_gc')
        self.assertTrue(chars.use_adaptive_qa)
        self.assertIn('Extreme GC-rich', chars.adaptive_qa_reason)
        self.assertIn('CRITICAL', chars.adaptive_qa_reason)

        # Verify conditions
        polymerase = generator.select_polymerase()
        conditions = generator.generate_reaction_conditions(polymerase)

        self.assertAlmostEqual(conditions['genome_gc'], 0.72, places=2)
        self.assertTrue(conditions['use_adaptive_qa'])

    def test_polymerase_selection_for_at_rich(self):
        """AT-rich genomes should prefer Phi29."""
        francisella = self.create_test_genome('francisella', 50000, 0.32)

        generator = OptimalOligoGenerator(
            target_genome=str(francisella),
            output_dir=Path(self.temp_dir) / 'output'
        )

        polymerase = generator.select_polymerase()
        self.assertEqual(polymerase, 'phi29')

    def test_polymerase_selection_for_gc_rich(self):
        """GC-rich genomes should prefer EquiPhi29."""
        burkholderia = self.create_test_genome('burkholderia', 50000, 0.67)

        generator = OptimalOligoGenerator(
            target_genome=str(burkholderia),
            output_dir=Path(self.temp_dir) / 'output'
        )

        polymerase = generator.select_polymerase()
        self.assertEqual(polymerase, 'equiphi29')

    def test_conditions_dict_completeness(self):
        """Verify all necessary keys are in conditions dict."""
        francisella = self.create_test_genome('francisella', 50000, 0.32)

        generator = OptimalOligoGenerator(
            target_genome=str(francisella),
            output_dir=Path(self.temp_dir) / 'output'
        )

        polymerase = generator.select_polymerase()
        conditions = generator.generate_reaction_conditions(polymerase)

        # Required keys for adaptive QA integration
        required_keys = [
            'genome_gc',
            'use_adaptive_qa',
            'polymerase',
            'temperature',
            'kmer_range',
            'optimal_kmer'
        ]

        for key in required_keys:
            self.assertIn(key, conditions, f"Missing required key: {key}")


if __name__ == '__main__':
    # Run tests
    unittest.main(verbosity=2)
