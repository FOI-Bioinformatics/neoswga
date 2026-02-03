#!/usr/bin/env python3
"""
Create test dataset from real Francisella and Burkholderia genomes.

Extracts representative genome regions and generates candidate primers
for systematic benchmarking of different argument combinations.

The dataset is realistic (uses actual sequences) but sized for practical
benchmarking (~500-1000 primers, completes in minutes not hours).
"""

import sys
import os
import h5py
import json
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Set
from collections import defaultdict
from Bio import SeqIO
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


class TestDatasetGenerator:
    """Generate test datasets from real genomes"""

    def __init__(self,
                 output_dir: str = "./test_data",
                 primer_length: int = 12,
                 target_primers: int = 1000):
        """
        Initialize generator.

        Args:
            output_dir: Directory for output files
            primer_length: Length of primers to generate
            target_primers: Target number of candidate primers
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)

        self.primer_length = primer_length
        self.target_primers = target_primers

        # Paths to real genomes
        self.genome_dir = Path("/Users/andreassjodin/Code/swga-dev/test")
        self.genome_files = {
            'francisella': self.genome_dir / "francisella_GCF_000008985.1_ASM898v1_genomic.fna",
            'burkholderia': self.genome_dir / "burkholderia_GCF_030297255.1_ASM3029725v1_genomic.fna",
            'bacillus': self.genome_dir / "bacillus_GCF_000008445.1_ASM844v1_genomic.fna",
            'ecoli': self.genome_dir / "ecoli_GCF_000005845.2_ASM584v2_genomic.fna",
            'yersinia': self.genome_dir / "yersinia_GCF_000222975.1_ASM22297v1_genomic.fna"
        }

    def generate_dataset(self, genome_file: Path, genome_name: str,
                        gc_tolerance: float = 0.15) -> Dict:
        """
        Generate test dataset for one genome.

        Args:
            genome_file: Path to genome FASTA
            genome_name: Name for output files
            gc_tolerance: GC filter tolerance

        Returns:
            Dataset metadata dictionary
        """
        logger.info(f"\nGenerating dataset for {genome_name}...")

        # Load genome
        records = list(SeqIO.parse(genome_file, 'fasta'))
        genome_seq = ''.join(str(rec.seq) for rec in records).upper()
        genome_length = len(genome_seq)

        # Calculate GC content
        gc_count = genome_seq.count('G') + genome_seq.count('C')
        genome_gc = gc_count / genome_length

        logger.info(f"  Length: {genome_length:,} bp")
        logger.info(f"  GC content: {genome_gc:.1%}")
        logger.info(f"  Chromosomes: {len(records)}")

        # Adaptive GC range
        gc_min = max(0.20, genome_gc - gc_tolerance)
        gc_max = min(0.80, genome_gc + gc_tolerance)

        logger.info(f"  GC range: {gc_min:.1%} - {gc_max:.1%}")

        # Extract candidate primers
        logger.info(f"  Extracting {self.primer_length}-mers...")

        candidates = self._extract_primers(
            genome_seq,
            gc_min,
            gc_max,
            target_count=self.target_primers
        )

        logger.info(f"  Extracted {len(candidates)} candidate primers")

        # Compute binding positions for each primer
        logger.info("  Computing binding positions...")
        primer_positions = self._compute_positions(genome_seq, candidates)

        # Save to HDF5
        h5_file = self.output_dir / f"{genome_name}_test.h5"
        logger.info(f"  Saving to {h5_file}...")
        self._save_hdf5(h5_file, genome_name, candidates, primer_positions, genome_length)

        # Save primer list
        primer_file = self.output_dir / f"{genome_name}_primers.txt"
        with open(primer_file, 'w') as f:
            for primer in candidates:
                f.write(f"{primer}\n")

        # Metadata
        metadata = {
            'genome_name': genome_name,
            'genome_length': genome_length,
            'gc_content': genome_gc,
            'gc_min': gc_min,
            'gc_max': gc_max,
            'gc_tolerance': gc_tolerance,
            'n_primers': len(candidates),
            'primer_length': self.primer_length,
            'h5_file': str(h5_file),
            'primer_file': str(primer_file),
            'total_positions': sum(len(pos) for pos in primer_positions.values())
        }

        logger.info(f"  Total binding sites: {metadata['total_positions']:,}")
        logger.info(f"  Avg sites per primer: {metadata['total_positions']/len(candidates):.1f}")

        return metadata

    def _extract_primers(self, genome_seq: str, gc_min: float, gc_max: float,
                        target_count: int) -> List[str]:
        """
        Extract candidate primers from genome.

        Args:
            genome_seq: Genome sequence
            gc_min: Minimum GC content
            gc_max: Maximum GC content
            target_count: Target number of primers

        Returns:
            List of candidate primers
        """
        candidates = set()

        # Calculate sampling step to get roughly target_count primers
        genome_length = len(genome_seq)
        step = max(1, (genome_length - self.primer_length) // (target_count * 2))

        logger.info(f"    Sampling every {step} bp...")

        for i in range(0, genome_length - self.primer_length, step):
            primer = genome_seq[i:i+self.primer_length]

            # Filter: no ambiguous bases
            if 'N' in primer:
                continue

            # Filter: GC content
            gc = (primer.count('G') + primer.count('C')) / self.primer_length
            if not (gc_min <= gc <= gc_max):
                continue

            # Filter: no homopolymers >4
            if any(base*5 in primer for base in 'ACGT'):
                continue

            candidates.add(primer)

            if len(candidates) >= target_count * 1.5:  # Get a bit extra
                break

        # Convert to sorted list for reproducibility
        candidates = sorted(list(candidates))

        # If we have too many, sample evenly
        if len(candidates) > target_count:
            step_sample = len(candidates) // target_count
            candidates = [candidates[i] for i in range(0, len(candidates), step_sample)][:target_count]

        return candidates

    def _compute_positions(self, genome_seq: str, primers: List[str]) -> Dict[str, List[int]]:
        """
        Compute all binding positions for primers.

        Args:
            genome_seq: Genome sequence
            primers: List of primers

        Returns:
            Dictionary: primer -> list of positions
        """
        positions = {}

        for i, primer in enumerate(primers):
            if (i + 1) % 100 == 0:
                logger.info(f"    Processed {i+1}/{len(primers)} primers...")

            # Find all occurrences (forward strand only for simplicity)
            pos_list = []
            start = 0
            while True:
                pos = genome_seq.find(primer, start)
                if pos == -1:
                    break
                pos_list.append(pos)
                start = pos + 1

            positions[primer] = pos_list

        return positions

    def _save_hdf5(self, filename: Path, genome_name: str, primers: List[str],
                   positions: Dict[str, List[int]], genome_length: int):
        """
        Save dataset to HDF5 file.

        Args:
            filename: Output HDF5 file
            genome_name: Genome identifier
            primers: List of primers
            positions: Primer positions dictionary
            genome_length: Total genome length
        """
        with h5py.File(filename, 'w') as f:
            # Metadata
            f.attrs['genome_name'] = genome_name
            f.attrs['genome_length'] = genome_length
            f.attrs['n_primers'] = len(primers)
            f.attrs['primer_length'] = self.primer_length

            # Create group for this genome
            genome_group = f.create_group(genome_name)

            # Save each primer's positions
            for primer in primers:
                pos_array = np.array(positions[primer], dtype=np.int32)

                # Forward strand
                genome_group.create_dataset(
                    f"{primer}/+",
                    data=pos_array,
                    compression='gzip'
                )

                # Reverse strand (empty for now, can be populated if needed)
                genome_group.create_dataset(
                    f"{primer}/-",
                    data=np.array([], dtype=np.int32),
                    compression='gzip'
                )

    def generate_combined_dataset(self) -> Dict:
        """
        Generate datasets for all genomes.

        Returns:
            Combined metadata dictionary
        """
        logger.info("="*80)
        logger.info("GENERATING TEST DATASETS FROM REAL GENOMES")
        logger.info("="*80)

        # Generate dataset for each genome
        all_metadata = {}
        for genome_name, genome_file in self.genome_files.items():
            metadata = self.generate_dataset(
                genome_file,
                genome_name,
                gc_tolerance=0.15
            )
            all_metadata[genome_name] = metadata

        # Combined metadata
        combined_metadata = {
            'created': str(Path(__file__).name),
            'primer_length': self.primer_length,
            'target_primers_per_genome': self.target_primers,
            'genomes': all_metadata
        }

        # Save metadata
        metadata_file = self.output_dir / "test_metadata.json"
        with open(metadata_file, 'w') as f:
            json.dump(combined_metadata, f, indent=2)

        logger.info(f"\nMetadata saved to: {metadata_file}")

        # Summary
        logger.info("\n" + "="*80)
        logger.info("DATASET GENERATION COMPLETE")
        logger.info("="*80)
        logger.info("\nSummary:")
        for genome_name, metadata in all_metadata.items():
            logger.info(f"  {genome_name.capitalize()}: {metadata['n_primers']} primers, "
                       f"{metadata['total_positions']:,} binding sites, "
                       f"{metadata['gc_content']:.1%} GC")
        logger.info(f"\nOutput directory: {self.output_dir.absolute()}")
        logger.info(f"Files created:")
        for genome_name in all_metadata.keys():
            logger.info(f"  - {genome_name}_test.h5")
            logger.info(f"  - {genome_name}_primers.txt")
        logger.info(f"  - test_metadata.json")

        return combined_metadata


def main():
    """Main entry point"""
    import argparse

    parser = argparse.ArgumentParser(
        description='Generate test datasets from real genomes'
    )
    parser.add_argument(
        '--output-dir',
        default='./test_data',
        help='Output directory for test datasets'
    )
    parser.add_argument(
        '--target-primers',
        type=int,
        default=1000,
        help='Target number of candidate primers per genome (default: 1000)'
    )
    parser.add_argument(
        '--primer-length',
        type=int,
        default=12,
        help='Primer length (default: 12)'
    )

    args = parser.parse_args()

    # Generate datasets
    generator = TestDatasetGenerator(
        output_dir=args.output_dir,
        primer_length=args.primer_length,
        target_primers=args.target_primers
    )

    metadata = generator.generate_combined_dataset()

    logger.info("\nDatasets ready for benchmarking!")
    logger.info("Next: Run benchmark_suite.py to test all argument combinations")


if __name__ == '__main__':
    main()
