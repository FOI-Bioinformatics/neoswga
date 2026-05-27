"""
Genome library for reusable pre-calculated k-mer data.

Stores pre-computed Jellyfish k-mer counts and optional Bloom filters
in a central directory (~/.neoswga/genomes/) so they can be reused
across projects. This avoids recomputing expensive k-mer data for
large genomes such as human.

Usage:
    from neoswga.core.genome_library import GenomeLibrary

    library = GenomeLibrary()

    # Pre-calculate k-mers for human genome
    library.add('human-grch38', '/path/to/human.fna', role='bg')

    # List available genomes
    for entry in library.list():
        print(f"{entry.name}: {entry.genome_size} bp")

    # Use in pipeline (resolves name to prefix)
    prefix = library.get_kmer_prefix('human-grch38')
"""

import json
import logging
import os
import shutil
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Tuple

logger = logging.getLogger(__name__)


@dataclass
class GenomeLibraryEntry:
    """Metadata for a pre-calculated genome in the library."""

    name: str
    species: str
    fasta_path: str
    genome_size: int
    gc_content: float
    role: str = "bg"  # "bg" or "bl"
    kmer_dir: str = ""
    kmer_prefix: str = ""
    computed_k_ranges: List[Tuple[int, int]] = field(default_factory=list)
    bloom_path: Optional[str] = None
    created_date: str = ""

    def __post_init__(self) -> None:
        if not self.created_date:
            self.created_date = datetime.now().isoformat()

    def to_dict(self) -> dict:
        """Serialize to dictionary."""
        return {
            'name': self.name,
            'species': self.species,
            'fasta_path': self.fasta_path,
            'genome_size': self.genome_size,
            'gc_content': self.gc_content,
            'role': self.role,
            'kmer_dir': self.kmer_dir,
            'kmer_prefix': self.kmer_prefix,
            'computed_k_ranges': [list(r) for r in self.computed_k_ranges],
            'bloom_path': self.bloom_path,
            'created_date': self.created_date,
        }

    @classmethod
    def from_dict(cls, data: dict) -> 'GenomeLibraryEntry':
        """Deserialize from dictionary."""
        ranges = [tuple(r) for r in data.get('computed_k_ranges', [])]
        return cls(
            name=data['name'],
            species=data.get('species', ''),
            fasta_path=data['fasta_path'],
            genome_size=data.get('genome_size', 0),
            gc_content=data.get('gc_content', 0.0),
            role=data.get('role', 'bg'),
            kmer_dir=data.get('kmer_dir', ''),
            kmer_prefix=data.get('kmer_prefix', ''),
            computed_k_ranges=ranges,
            bloom_path=data.get('bloom_path'),
            created_date=data.get('created_date', ''),
        )


class GenomeLibrary:
    """Central library for pre-calculated genome k-mer data.

    Stores Jellyfish k-mer counts and optional Bloom filters in
    ~/.neoswga/genomes/ for reuse across projects.
    """

    def __init__(self, library_dir: Optional[str] = None) -> None:
        """Initialize genome library.

        Args:
            library_dir: Library directory path. Defaults to
                ~/.neoswga/genomes/.
        """
        if library_dir is None:
            library_dir = str(Path.home() / '.neoswga' / 'genomes')
        self.library_dir = library_dir
        os.makedirs(self.library_dir, exist_ok=True)

    def _metadata_path(self, name: str) -> str:
        """Path to metadata JSON for a genome entry."""
        return os.path.join(self.library_dir, name, 'metadata.json')

    def _load_entry(self, name: str) -> Optional[GenomeLibraryEntry]:
        """Load a single entry from its metadata file."""
        meta_path = self._metadata_path(name)
        if not os.path.exists(meta_path):
            return None
        try:
            with open(meta_path) as f:
                data = json.load(f)
            return GenomeLibraryEntry.from_dict(data)
        except Exception as e:
            logger.warning(f"Failed to load metadata for '{name}': {e}")
            return None

    def _save_entry(self, entry: GenomeLibraryEntry) -> None:
        """Save entry metadata to disk."""
        entry_dir = os.path.join(self.library_dir, entry.name)
        os.makedirs(entry_dir, exist_ok=True)
        meta_path = self._metadata_path(entry.name)
        with open(meta_path, 'w') as f:
            json.dump(entry.to_dict(), f, indent=2)

    def add(
        self,
        name: str,
        fasta_path: str,
        role: str = 'bg',
        species: str = '',
        k_ranges: Optional[List[Tuple[int, int]]] = None,
        build_bloom: str = 'auto',
    ) -> GenomeLibraryEntry:
        """Add a genome to the library and pre-calculate k-mer data.

        Args:
            name: Human-readable identifier (e.g., 'human-grch38').
            fasta_path: Path to genome FASTA file.
            role: 'bg' for background, 'bl' for blacklist.
            species: Species name (optional).
            k_ranges: List of (min_k, max_k) tuples to compute.
                Defaults to [(6, 12), (12, 18)].
            build_bloom: 'auto' builds Bloom filter for genomes > 50 Mbp,
                'yes' always builds, 'no' never builds.

        Returns:
            The created GenomeLibraryEntry.
        """
        if k_ranges is None:
            k_ranges = [(6, 12), (12, 18)]

        fasta_path = os.path.abspath(fasta_path)
        if not os.path.exists(fasta_path):
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

        entry_dir = os.path.join(self.library_dir, name)
        os.makedirs(entry_dir, exist_ok=True)

        # Symlink FASTA into library directory
        link_path = os.path.join(entry_dir, 'genome.fna')
        if os.path.exists(link_path) or os.path.islink(link_path):
            os.remove(link_path)
        os.symlink(fasta_path, link_path)

        # Calculate genome stats
        genome_size, gc_content = self._calculate_genome_stats(fasta_path)
        logger.info(
            f"Genome '{name}': {genome_size:,} bp, GC={gc_content:.1%}"
        )

        # K-mer prefix
        kmer_prefix = os.path.join(entry_dir, name)

        # Run Jellyfish for each k-range
        computed_ranges: List[Tuple[int, int]] = []
        for min_k, max_k in k_ranges:
            logger.info(
                f"Computing k-mers for {name} (k={min_k}-{max_k})..."
            )
            try:
                from neoswga.core.kmer_counter import run_jellyfish

                run_jellyfish(fasta_path, kmer_prefix, min_k, max_k)
                computed_ranges.append((min_k, max_k))
            except Exception as e:
                logger.error(f"Jellyfish failed for k={min_k}-{max_k}: {e}")

        # Optional Bloom filter
        bloom_path: Optional[str] = None
        should_bloom = (
            build_bloom == 'yes'
            or (build_bloom == 'auto' and genome_size > 50_000_000)
        )
        if should_bloom:
            bloom_path = os.path.join(entry_dir, 'bloom.pkl')
            logger.info(f"Building Bloom filter for {name}...")
            try:
                from neoswga.core.background_filter import (
                    BackgroundBloomFilter,
                )

                bloom = BackgroundBloomFilter()
                bloom.add_genome(fasta_path)
                bloom.save(bloom_path)
            except Exception as e:
                logger.warning(f"Bloom filter build failed: {e}")
                bloom_path = None

        entry = GenomeLibraryEntry(
            name=name,
            species=species,
            fasta_path=fasta_path,
            genome_size=genome_size,
            gc_content=gc_content,
            role=role,
            kmer_dir=entry_dir,
            kmer_prefix=kmer_prefix,
            computed_k_ranges=computed_ranges,
            bloom_path=bloom_path,
        )
        self._save_entry(entry)
        logger.info(f"Added genome '{name}' to library")
        return entry

    def list(self) -> List[GenomeLibraryEntry]:
        """List all genomes in the library.

        Returns:
            List of GenomeLibraryEntry objects.
        """
        entries: List[GenomeLibraryEntry] = []
        if not os.path.exists(self.library_dir):
            return entries
        for item in sorted(os.listdir(self.library_dir)):
            entry = self._load_entry(item)
            if entry is not None:
                entries.append(entry)
        return entries

    def get(self, name: str) -> Optional[GenomeLibraryEntry]:
        """Get a genome entry by name.

        Args:
            name: Genome identifier.

        Returns:
            GenomeLibraryEntry or None if not found.
        """
        return self._load_entry(name)

    def remove(self, name: str) -> bool:
        """Remove a genome and its data from the library.

        Args:
            name: Genome identifier.

        Returns:
            True if removed, False if not found.
        """
        entry_dir = os.path.join(self.library_dir, name)
        if not os.path.exists(entry_dir):
            return False
        shutil.rmtree(entry_dir)
        logger.info(f"Removed genome '{name}' from library")
        return True

    def get_kmer_prefix(self, name: str) -> Optional[str]:
        """Get the k-mer file prefix for a genome.

        Args:
            name: Genome identifier.

        Returns:
            Prefix path for k-mer files, or None if not found.
        """
        entry = self.get(name)
        if entry is None:
            return None
        return entry.kmer_prefix

    def has_kmers_for_range(self, name: str, min_k: int, max_k: int) -> bool:
        """Check whether pre-calculated k-mers cover a given range.

        Args:
            name: Genome identifier.
            min_k: Minimum k-mer length needed.
            max_k: Maximum k-mer length needed.

        Returns:
            True if all k-mer files exist for the range.
        """
        entry = self.get(name)
        if entry is None:
            return False
        # Verify actual files exist for every k in range
        for k in range(min_k, max_k + 1):
            kmer_file = f"{entry.kmer_prefix}_{k}mer_all.txt"
            if not os.path.exists(kmer_file):
                return False
        return True

    def _calculate_genome_stats(
        self, fasta_path: str
    ) -> Tuple[int, float]:
        """Calculate genome size and GC content from FASTA.

        Args:
            fasta_path: Path to FASTA file.

        Returns:
            Tuple of (genome_size, gc_content).
        """
        try:
            from neoswga.core.genome_io import GenomeLoader

            loader = GenomeLoader()
            sequence = loader.load_genome(fasta_path, return_stats=False)
            seq_upper = sequence.upper()
            gc_count = seq_upper.count('G') + seq_upper.count('C')
            total = len(seq_upper)
            gc = gc_count / total if total > 0 else 0.0
            return total, gc
        except Exception as e:
            logger.warning(f"Could not calculate genome stats: {e}")
            # Fallback: estimate from file size
            file_size = os.path.getsize(fasta_path)
            return file_size, 0.5
