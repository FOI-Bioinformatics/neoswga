"""
Background genome registry for SWGA.

Provides a database of pre-computed background genome data (Bloom filters,
k-mer files) that can be reused across projects. This avoids recomputing
expensive background data for common hosts like human.

Usage:
    from neoswga.core.background_registry import BackgroundRegistry

    registry = BackgroundRegistry()

    # List available backgrounds
    for bg in registry.list_all():
        print(f"{bg.name}: {bg.species}")

    # Add a new background
    registry.add(
        name="Human GRCh38",
        species="Homo sapiens",
        bloom_path="/path/to/human_bloom.pkl",
        genome_size=3_000_000_000,
    )

    # Get background by name
    human = registry.get("Human GRCh38")
    print(f"Bloom filter: {human.bloom_path}")
"""

import json
import logging
import os
from dataclasses import dataclass, asdict, field
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from glob import glob

logger = logging.getLogger(__name__)


@dataclass
class BackgroundEntry:
    """
    Registry entry for a pre-computed background genome.

    Attributes:
        name: Human-readable name (e.g., "Human GRCh38")
        species: Species name
        genome_size: Genome size in bp
        bloom_path: Path to Bloom filter pickle file (optional)
        kmer_prefix: Prefix for k-mer files (optional)
        k_range: Range of k-mer lengths available
        created_date: Date entry was created
        description: Additional description
        source: Source/reference for the genome
    """
    name: str
    species: str
    genome_size: int
    bloom_path: Optional[str] = None
    kmer_prefix: Optional[str] = None
    k_range: Tuple[int, int] = (6, 12)
    created_date: str = ""
    description: str = ""
    source: str = ""

    def __post_init__(self):
        if not self.created_date:
            self.created_date = datetime.now().isoformat()

    @property
    def has_bloom(self) -> bool:
        """Whether Bloom filter is available."""
        if not self.bloom_path:
            return False
        return os.path.exists(self.bloom_path)

    @property
    def has_kmers(self) -> bool:
        """Whether k-mer files are available."""
        if not self.kmer_prefix:
            return False
        # Check if at least one k-mer file exists
        for k in range(self.k_range[0], self.k_range[1] + 1):
            kmer_file = f"{self.kmer_prefix}_{k}mer_all.txt"
            if os.path.exists(kmer_file):
                return True
        return False

    def to_dict(self) -> Dict:
        """Convert to dictionary for serialization."""
        return {
            'name': self.name,
            'species': self.species,
            'genome_size': self.genome_size,
            'bloom_path': self.bloom_path,
            'kmer_prefix': self.kmer_prefix,
            'k_range': list(self.k_range),
            'created_date': self.created_date,
            'description': self.description,
            'source': self.source,
        }

    @classmethod
    def from_dict(cls, data: Dict) -> 'BackgroundEntry':
        """Create from dictionary."""
        k_range = tuple(data.get('k_range', [6, 12]))
        return cls(
            name=data['name'],
            species=data['species'],
            genome_size=data['genome_size'],
            bloom_path=data.get('bloom_path'),
            kmer_prefix=data.get('kmer_prefix'),
            k_range=k_range,
            created_date=data.get('created_date', ''),
            description=data.get('description', ''),
            source=data.get('source', ''),
        )

    def __str__(self) -> str:
        size_gb = self.genome_size / 1e9
        lines = [
            f"{self.name} ({self.species})",
            f"  Size: {size_gb:.2f} Gbp",
        ]
        if self.has_bloom:
            lines.append(f"  Bloom filter: {self.bloom_path}")
        if self.has_kmers:
            lines.append(f"  K-mer prefix: {self.kmer_prefix}")
            lines.append(f"  K range: {self.k_range[0]}-{self.k_range[1]}")
        if self.description:
            lines.append(f"  Description: {self.description}")
        return "\n".join(lines)


class BackgroundRegistry:
    """
    Registry for pre-computed background genome data.

    Maintains a database of available background genomes with their
    Bloom filters and k-mer file locations.
    """

    # Default search directories
    DEFAULT_DIRS = [
        '~/.neoswga/backgrounds',
        './backgrounds',
        './filters',
    ]

    def __init__(
        self,
        registry_path: Optional[str] = None,
        auto_discover: bool = True,
    ):
        """
        Initialize background registry.

        Args:
            registry_path: Path to registry JSON file.
                          If None, uses ~/.neoswga/background_registry.json
            auto_discover: Whether to auto-discover backgrounds on init
        """
        if registry_path is None:
            registry_dir = Path.home() / '.neoswga'
            registry_dir.mkdir(exist_ok=True)
            registry_path = str(registry_dir / 'background_registry.json')

        self.registry_path = registry_path
        self.entries: Dict[str, BackgroundEntry] = {}
        self._load()

        if auto_discover:
            self.discover()

    def _load(self) -> None:
        """Load registry from file."""
        if os.path.exists(self.registry_path):
            try:
                with open(self.registry_path) as f:
                    data = json.load(f)
                    for entry_data in data.get('entries', []):
                        entry = BackgroundEntry.from_dict(entry_data)
                        self.entries[entry.name] = entry
                logger.debug(f"Loaded {len(self.entries)} background entries")
            except Exception as e:
                logger.warning(f"Failed to load registry: {e}")
                self.entries = {}
        else:
            self.entries = {}

    def _save(self) -> None:
        """Save registry to file."""
        try:
            data = {
                'version': '1.0',
                'updated': datetime.now().isoformat(),
                'entries': [entry.to_dict() for entry in self.entries.values()],
            }
            with open(self.registry_path, 'w') as f:
                json.dump(data, f, indent=2)
        except Exception as e:
            logger.error(f"Failed to save registry: {e}")

    def add(
        self,
        name: str,
        species: str,
        genome_size: int,
        bloom_path: Optional[str] = None,
        kmer_prefix: Optional[str] = None,
        k_range: Tuple[int, int] = (6, 12),
        description: str = "",
        source: str = "",
        overwrite: bool = False,
    ) -> BackgroundEntry:
        """
        Add a background to the registry.

        Args:
            name: Human-readable name
            species: Species name
            genome_size: Genome size in bp
            bloom_path: Path to Bloom filter
            kmer_prefix: Prefix for k-mer files
            k_range: Range of k-mer lengths
            description: Additional description
            source: Source/reference
            overwrite: Whether to overwrite existing entry

        Returns:
            Created BackgroundEntry
        """
        if name in self.entries and not overwrite:
            logger.warning(f"Entry '{name}' already exists. Use overwrite=True to replace.")
            return self.entries[name]

        entry = BackgroundEntry(
            name=name,
            species=species,
            genome_size=genome_size,
            bloom_path=bloom_path,
            kmer_prefix=kmer_prefix,
            k_range=k_range,
            description=description,
            source=source,
        )

        self.entries[name] = entry
        self._save()

        logger.info(f"Added background: {name}")
        return entry

    def remove(self, name: str) -> bool:
        """
        Remove a background from the registry.

        Does not delete the actual files, only the registry entry.

        Args:
            name: Name of background to remove

        Returns:
            True if removed, False if not found
        """
        if name in self.entries:
            del self.entries[name]
            self._save()
            logger.info(f"Removed background: {name}")
            return True
        return False

    def get(self, name: str) -> Optional[BackgroundEntry]:
        """
        Get a background by name.

        Args:
            name: Name of background

        Returns:
            BackgroundEntry or None if not found
        """
        return self.entries.get(name)

    def list_all(self) -> List[BackgroundEntry]:
        """
        List all backgrounds in registry.

        Returns:
            List of BackgroundEntry objects
        """
        return list(self.entries.values())

    def discover(
        self,
        directories: Optional[List[str]] = None,
        verbose: bool = False,
    ) -> int:
        """
        Auto-discover background data in directories.

        Searches for:
        - Bloom filter files (*_bloom.pkl, *_bloom_filter.pkl)
        - K-mer file sets (*_Xmer_all.txt)

        Args:
            directories: Directories to search (default: DEFAULT_DIRS)
            verbose: Print discovery progress

        Returns:
            Number of new backgrounds discovered
        """
        if directories is None:
            directories = self.DEFAULT_DIRS

        discovered = 0

        for dir_pattern in directories:
            dir_path = os.path.expanduser(dir_pattern)

            if not os.path.isdir(dir_path):
                continue

            if verbose:
                logger.info(f"Searching {dir_path}...")

            # Look for Bloom filters
            bloom_files = glob(os.path.join(dir_path, '*_bloom*.pkl'))
            for bloom_file in bloom_files:
                name = self._infer_name_from_bloom(bloom_file)
                if name and name not in self.entries:
                    # Try to infer genome info
                    self.add(
                        name=name,
                        species=name,  # Use name as species if unknown
                        genome_size=0,  # Unknown
                        bloom_path=bloom_file,
                        description="Auto-discovered from Bloom filter",
                    )
                    discovered += 1
                    if verbose:
                        logger.info(f"  Found: {name}")

            # Look for k-mer file sets
            kmer_files = glob(os.path.join(dir_path, '*_6mer_all.txt'))
            for kmer_file in kmer_files:
                prefix = kmer_file.replace('_6mer_all.txt', '')
                name = self._infer_name_from_prefix(prefix)
                if name and name not in self.entries:
                    # Determine k-mer range
                    k_min, k_max = self._detect_k_range(prefix)
                    self.add(
                        name=name,
                        species=name,
                        genome_size=0,
                        kmer_prefix=prefix,
                        k_range=(k_min, k_max),
                        description="Auto-discovered from k-mer files",
                    )
                    discovered += 1
                    if verbose:
                        logger.info(f"  Found: {name} (k={k_min}-{k_max})")

        if discovered > 0:
            self._save()

        return discovered

    def _infer_name_from_bloom(self, bloom_path: str) -> Optional[str]:
        """Infer background name from Bloom filter path."""
        basename = os.path.basename(bloom_path)
        # Remove common suffixes
        name = basename
        for suffix in ['_bloom_filter.pkl', '_bloom.pkl', '.pkl']:
            if name.endswith(suffix):
                name = name[:-len(suffix)]
                break
        return name if name else None

    def _infer_name_from_prefix(self, prefix: str) -> Optional[str]:
        """Infer background name from k-mer file prefix."""
        basename = os.path.basename(prefix)
        return basename if basename else None

    def _detect_k_range(self, prefix: str) -> Tuple[int, int]:
        """Detect available k-mer range from files."""
        k_min = 30
        k_max = 1

        for k in range(4, 31):
            kmer_file = f"{prefix}_{k}mer_all.txt"
            if os.path.exists(kmer_file):
                k_min = min(k_min, k)
                k_max = max(k_max, k)

        if k_min > k_max:
            return (6, 12)  # Default

        return (k_min, k_max)

    def search(self, query: str) -> List[BackgroundEntry]:
        """
        Search for backgrounds by name or species.

        Args:
            query: Search string (case-insensitive)

        Returns:
            List of matching BackgroundEntry objects
        """
        query_lower = query.lower()
        results = []

        for entry in self.entries.values():
            if (query_lower in entry.name.lower() or
                query_lower in entry.species.lower()):
                results.append(entry)

        return results


# Common pre-defined backgrounds
COMMON_BACKGROUNDS = {
    'human': {
        'name': 'Human GRCh38',
        'species': 'Homo sapiens',
        'genome_size': 3_088_286_401,
        'source': 'NCBI GRCh38',
    },
    'mouse': {
        'name': 'Mouse GRCm39',
        'species': 'Mus musculus',
        'genome_size': 2_728_222_451,
        'source': 'NCBI GRCm39',
    },
    'ecoli': {
        'name': 'E. coli K-12',
        'species': 'Escherichia coli',
        'genome_size': 4_641_652,
        'source': 'NCBI NC_000913',
    },
}


def get_common_background(name: str) -> Optional[Dict]:
    """
    Get pre-defined information for common backgrounds.

    Args:
        name: Background name (e.g., 'human', 'mouse')

    Returns:
        Dictionary with background info or None
    """
    return COMMON_BACKGROUNDS.get(name.lower())
