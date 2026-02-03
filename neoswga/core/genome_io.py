#!/usr/bin/env python3
"""
Genome I/O utilities with support for multiple formats.

Supports:
- Plain text: .fasta, .fa, .fna
- Gzip compressed: .fasta.gz, .fa.gz, .fna.gz, .gz
- Zip compressed: .fasta.zip, .fa.zip, .fna.zip, .zip
- Auto-detection based on file extension and magic bytes

Features:
- Automatic format detection
- Memory-efficient streaming for large files
- Multi-sequence FASTA handling
- Quality validation
- GC content calculation
- Statistics generation

Author: NeoSWGA Development Team
Date: November 2025
Version: 3.0 - Genome I/O Enhancement
"""

import gzip
import zipfile
import logging
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Union, Iterator
from dataclasses import dataclass
from io import TextIOWrapper
from Bio import SeqIO

logger = logging.getLogger(__name__)


@dataclass
class GenomeStats:
    """Statistics for a loaded genome"""
    file_path: Path
    file_format: str  # "fasta", "fasta.gz", "fasta.zip"
    total_sequences: int
    total_length: int
    gc_content: float
    n_count: int
    n_fraction: float
    sequence_lengths: List[int]
    mean_length: float
    max_length: int
    min_length: int

    def __str__(self):
        return f"""Genome Statistics:
  File: {self.file_path.name}
  Format: {self.file_format}
  Sequences: {self.total_sequences}
  Total length: {self.total_length:,} bp
  GC content: {self.gc_content:.1%}
  N content: {self.n_fraction:.2%} ({self.n_count:,} bp)
  Sequence lengths: {self.min_length:,} - {self.max_length:,} bp (mean: {self.mean_length:,.0f})"""


class GenomeLoader:
    """
    Unified genome loader supporting multiple compression formats.

    Automatically detects and handles:
    - Plain FASTA/FA/FNA
    - Gzip compressed (.gz)
    - Zip compressed (.zip)

    Usage:
        loader = GenomeLoader()
        sequence = loader.load_genome("genome.fasta.gz")
        stats = loader.get_stats()
    """

    def __init__(self):
        """Initialize genome loader"""
        self.last_stats: Optional[GenomeStats] = None

    def detect_format(self, file_path: Path) -> str:
        """
        Detect file format from extension and content.

        Args:
            file_path: Path to genome file

        Returns:
            Format string: "fasta", "fasta.gz", "fasta.zip"
        """
        if not file_path.exists():
            raise FileNotFoundError(f"Genome file not found: {file_path}")

        # Check extension
        suffixes = file_path.suffixes
        name_lower = file_path.name.lower()

        # Check for compression
        if name_lower.endswith('.gz'):
            return "fasta.gz"
        elif name_lower.endswith('.zip'):
            return "fasta.zip"
        elif any(ext in name_lower for ext in ['.fasta', '.fa', '.fna']):
            return "fasta"

        # Check magic bytes if extension unclear
        with open(file_path, 'rb') as f:
            magic = f.read(2)

        if magic == b'\x1f\x8b':  # Gzip magic bytes
            return "fasta.gz"
        elif magic == b'PK':  # Zip magic bytes
            return "fasta.zip"

        # Default to plain FASTA
        return "fasta"

    def _open_plain(self, file_path: Path) -> TextIOWrapper:
        """Open plain text file"""
        return open(file_path, 'r')

    def _open_gzip(self, file_path: Path) -> TextIOWrapper:
        """Open gzip compressed file"""
        return gzip.open(file_path, 'rt')

    def _open_zip(self, file_path: Path) -> TextIOWrapper:
        """
        Open zip compressed file.

        For zip files, extracts the first .fasta/.fa/.fna file found.
        """
        zf = zipfile.ZipFile(file_path, 'r')

        # Find FASTA file in zip
        fasta_files = [
            name for name in zf.namelist()
            if name.lower().endswith(('.fasta', '.fa', '.fna'))
        ]

        if not fasta_files:
            raise ValueError(f"No FASTA file found in zip: {file_path}")

        if len(fasta_files) > 1:
            logger.warning(f"Multiple FASTA files in zip, using: {fasta_files[0]}")

        # Return text wrapper for the FASTA file
        return TextIOWrapper(zf.open(fasta_files[0], 'r'), encoding='utf-8')

    def load_genome(self, file_path: Union[str, Path],
                   return_stats: bool = True) -> str:
        """
        Load genome sequence from file (auto-detects compression).

        Args:
            file_path: Path to genome file (plain, .gz, or .zip)
            return_stats: Calculate and store statistics

        Returns:
            Complete genome sequence as string (concatenated if multi-sequence)
        """
        file_path = Path(file_path)

        logger.info(f"Loading genome: {file_path}")

        # Detect format
        format_type = self.detect_format(file_path)
        logger.info(f"  Detected format: {format_type}")

        # Open file with appropriate method
        if format_type == "fasta.gz":
            file_handle = self._open_gzip(file_path)
        elif format_type == "fasta.zip":
            file_handle = self._open_zip(file_path)
        else:
            file_handle = self._open_plain(file_path)

        # Parse FASTA
        sequences = []
        sequence_lengths = []
        total_n = 0

        try:
            for record in SeqIO.parse(file_handle, "fasta"):
                seq_str = str(record.seq).upper()
                sequences.append(seq_str)
                sequence_lengths.append(len(seq_str))
                total_n += seq_str.count('N')

        finally:
            file_handle.close()

        if not sequences:
            raise ValueError(f"No sequences found in: {file_path}")

        # Concatenate sequences
        full_sequence = "".join(sequences)
        total_length = len(full_sequence)

        logger.info(f"  Loaded {len(sequences)} sequence(s), {total_length:,} bp total")

        # Calculate statistics if requested
        if return_stats:
            gc_count = full_sequence.count('G') + full_sequence.count('C')
            gc_content = gc_count / total_length if total_length > 0 else 0.0

            self.last_stats = GenomeStats(
                file_path=file_path,
                file_format=format_type,
                total_sequences=len(sequences),
                total_length=total_length,
                gc_content=gc_content,
                n_count=total_n,
                n_fraction=total_n / total_length if total_length > 0 else 0.0,
                sequence_lengths=sequence_lengths,
                mean_length=sum(sequence_lengths) / len(sequence_lengths),
                max_length=max(sequence_lengths),
                min_length=min(sequence_lengths)
            )

            logger.info(f"  GC content: {gc_content:.1%}")
            if total_n > 0:
                logger.info(f"  N content: {total_n:,} bp ({total_n/total_length:.2%})")

        return full_sequence

    def get_stats(self) -> Optional[GenomeStats]:
        """
        Get statistics from last loaded genome.

        Returns:
            GenomeStats object or None if no genome loaded yet
        """
        return self.last_stats

    def load_genome_streaming(self, file_path: Union[str, Path]) -> Iterator[str]:
        """
        Load genome sequence by sequence (memory-efficient for large files).

        Args:
            file_path: Path to genome file

        Yields:
            Individual sequences as strings
        """
        file_path = Path(file_path)
        format_type = self.detect_format(file_path)

        if format_type == "fasta.gz":
            file_handle = self._open_gzip(file_path)
        elif format_type == "fasta.zip":
            file_handle = self._open_zip(file_path)
        else:
            file_handle = self._open_plain(file_path)

        try:
            for record in SeqIO.parse(file_handle, "fasta"):
                yield str(record.seq).upper()
        finally:
            file_handle.close()

    def validate_genome(self, file_path: Union[str, Path],
                       max_n_fraction: float = 0.10,
                       min_length: int = 100000) -> Tuple[bool, List[str]]:
        """
        Validate genome file quality.

        Args:
            file_path: Path to genome file
            max_n_fraction: Maximum acceptable fraction of N bases
            min_length: Minimum acceptable genome length (bp)

        Returns:
            Tuple of (is_valid, list_of_issues)
        """
        issues = []

        try:
            sequence = self.load_genome(file_path, return_stats=True)
            stats = self.last_stats

            # Check length
            if stats.total_length < min_length:
                issues.append(f"Genome too short: {stats.total_length:,} bp "
                            f"(minimum: {min_length:,} bp)")

            # Check N content
            if stats.n_fraction > max_n_fraction:
                issues.append(f"Too many N bases: {stats.n_fraction:.1%} "
                            f"(maximum: {max_n_fraction:.1%})")

            # Check for empty sequences
            if any(length == 0 for length in stats.sequence_lengths):
                issues.append("Contains empty sequences")

            # Warn about unusual GC
            if stats.gc_content < 0.15 or stats.gc_content > 0.85:
                issues.append(f"Unusual GC content: {stats.gc_content:.1%} "
                            f"(expected 15-85%)")

        except Exception as e:
            issues.append(f"Failed to load genome: {str(e)}")

        is_valid = len(issues) == 0

        return is_valid, issues


def load_genome(file_path: Union[str, Path]) -> str:
    """
    Convenience function to load a genome file (auto-detects compression).

    Args:
        file_path: Path to genome file (.fasta, .fasta.gz, .fasta.zip)

    Returns:
        Complete genome sequence as string
    """
    loader = GenomeLoader()
    return loader.load_genome(file_path)


def get_genome_stats(file_path: Union[str, Path]) -> GenomeStats:
    """
    Convenience function to get genome statistics.

    Args:
        file_path: Path to genome file

    Returns:
        GenomeStats object with complete statistics
    """
    loader = GenomeLoader()
    loader.load_genome(file_path, return_stats=True)
    return loader.get_stats()


def validate_genome_file(file_path: Union[str, Path]) -> bool:
    """
    Convenience function to validate genome file.

    Args:
        file_path: Path to genome file

    Returns:
        True if valid, raises exception with details if invalid
    """
    loader = GenomeLoader()
    is_valid, issues = loader.validate_genome(file_path)

    if not is_valid:
        raise ValueError(f"Invalid genome file:\n" + "\n".join(f"  - {issue}" for issue in issues))

    return True


class GenomeCache:
    """
    Cache for loaded genomes to avoid re-reading.

    Useful when same genome is used multiple times.
    """

    def __init__(self, max_cache_size: int = 5):
        """
        Initialize genome cache.

        Args:
            max_cache_size: Maximum number of genomes to cache
        """
        self.cache: Dict[str, Tuple[str, GenomeStats]] = {}
        self.max_cache_size = max_cache_size
        self.access_order: List[str] = []
        self.loader = GenomeLoader()

    def get(self, file_path: Union[str, Path]) -> Tuple[str, GenomeStats]:
        """
        Get genome from cache or load if not cached.

        Args:
            file_path: Path to genome file

        Returns:
            Tuple of (sequence, stats)
        """
        key = str(Path(file_path).resolve())

        # Check cache
        if key in self.cache:
            logger.debug(f"Using cached genome: {file_path}")
            # Update access order
            if key in self.access_order:
                self.access_order.remove(key)
            self.access_order.append(key)
            return self.cache[key]

        # Load genome
        logger.debug(f"Loading genome (not cached): {file_path}")
        sequence = self.loader.load_genome(file_path, return_stats=True)
        stats = self.loader.get_stats()

        # Add to cache
        self.cache[key] = (sequence, stats)
        self.access_order.append(key)

        # Evict oldest if cache full
        if len(self.cache) > self.max_cache_size:
            oldest_key = self.access_order.pop(0)
            del self.cache[oldest_key]
            logger.debug(f"Evicted from cache: {oldest_key}")

        return sequence, stats

    def clear(self):
        """Clear the cache"""
        self.cache.clear()
        self.access_order.clear()
        logger.debug("Cache cleared")


if __name__ == '__main__':
    # Example usage
    logging.basicConfig(level=logging.INFO)

    print("="*80)
    print("Genome I/O - Example Usage")
    print("="*80)

    # Create test files to demonstrate (in practice, use real genome files)
    from pathlib import Path

    test_dir = Path("/Users/andreassjodin/Code/swga-dev/test")

    if test_dir.exists():
        genomes = list(test_dir.glob("*.fna"))[:3]

        print(f"\nFound {len(genomes)} test genomes\n")

        for genome_file in genomes:
            print("-"*80)

            # Load genome
            loader = GenomeLoader()
            sequence = loader.load_genome(genome_file)
            stats = loader.get_stats()

            # Print statistics
            print(stats)

            # Validate
            is_valid, issues = loader.validate_genome(genome_file)
            if is_valid:
                print("\n✓ Genome validation passed")
            else:
                print("\n✗ Genome validation failed:")
                for issue in issues:
                    print(f"  - {issue}")

        # Demonstrate cache
        print("\n" + "="*80)
        print("Testing genome cache...")
        print("="*80)

        cache = GenomeCache(max_cache_size=2)

        for i in range(3):
            print(f"\nIteration {i+1}:")
            for genome_file in genomes[:2]:
                seq, stats = cache.get(genome_file)
                print(f"  Loaded {genome_file.name}: {len(seq):,} bp")

    else:
        print(f"Test directory not found: {test_dir}")
        print("\nDemonstrating API usage:")

        print("\n# Simple loading")
        print("sequence = load_genome('genome.fasta.gz')")

        print("\n# Get statistics")
        print("stats = get_genome_stats('genome.fasta.gz')")
        print("print(f'GC content: {stats.gc_content:.1%}')")

        print("\n# Validation")
        print("try:")
        print("    validate_genome_file('genome.fasta')")
        print("    print('Valid genome')")
        print("except ValueError as e:")
        print("    print(f'Invalid: {e}')")

        print("\n# Using cache for repeated access")
        print("cache = GenomeCache()")
        print("seq, stats = cache.get('genome.fasta.gz')")
