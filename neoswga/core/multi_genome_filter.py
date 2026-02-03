#!/usr/bin/env python3
"""
Multi-genome filtering for SWGA with differential penalties.

Supports three genome categories:
1. Target genome(s): Maximize primer binding (high frequency desired)
2. Background genome(s): Minimize but tolerate (e.g., human, mouse, tick host)
3. Blacklist genome(s): Strongly avoid (e.g., co-existing bacteria we don't want)

Key innovation: Differential penalty weights
- Background genomes: penalty_weight = 1.0 (standard avoidance)
- Blacklist genomes: penalty_weight = 5.0-10.0 (strong avoidance)

This allows primers that occasionally bind to host DNA (acceptable) while
completely avoiding primers that bind to blacklisted organisms (unacceptable).

Example use cases:
- Detect Borrelia (Lyme) in tick, avoid amplifying tick DNA and Rickettsia
- Detect Plasmodium in human blood, avoid human DNA and other Plasmodium species
- Detect Mycobacterium in sputum, avoid human DNA and other mycobacteria

Author: NeoSWGA Development Team
Date: November 2025
Version: 3.0 - Multi-Genome Support
"""

import logging
import numpy as np
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Set
from dataclasses import dataclass, field
from enum import Enum

logger = logging.getLogger(__name__)


class GenomeRole(Enum):
    """Role of a genome in the SWGA experiment"""
    TARGET = "target"           # Maximize binding
    BACKGROUND = "background"   # Minimize but tolerate (host)
    BLACKLIST = "blacklist"     # Strongly avoid (contaminants)


@dataclass
class GenomeEntry:
    """
    Single genome in a multi-genome SWGA experiment.

    Attributes:
        name: Descriptive name (e.g., "Human", "Borrelia_burgdorferi")
        fasta_path: Path to genome FASTA file
        role: TARGET, BACKGROUND, or BLACKLIST
        penalty_weight: Multiplier for off-target binding penalty
        gc_content: Pre-calculated GC content (optional)
        size: Genome size in bp (optional)
    """
    name: str
    fasta_path: Path
    role: GenomeRole
    penalty_weight: float = 1.0
    gc_content: Optional[float] = None
    size: Optional[int] = None

    def __post_init__(self):
        """Validate and set defaults"""
        if isinstance(self.fasta_path, str):
            self.fasta_path = Path(self.fasta_path)

        if not self.fasta_path.exists():
            raise FileNotFoundError(f"Genome file not found: {self.fasta_path}")

        # Set default penalty weights based on role
        if self.role == GenomeRole.TARGET:
            self.penalty_weight = 0.0  # No penalty for target binding
        elif self.role == GenomeRole.BACKGROUND:
            if self.penalty_weight == 1.0:  # Use default
                self.penalty_weight = 1.0
        elif self.role == GenomeRole.BLACKLIST:
            if self.penalty_weight == 1.0:  # Use default
                self.penalty_weight = 5.0  # Much stronger avoidance


@dataclass
class GenomeSet:
    """
    Collection of genomes for multi-genome SWGA.

    Manages target, background, and blacklist genomes with validation.
    """
    targets: List[GenomeEntry] = field(default_factory=list)
    backgrounds: List[GenomeEntry] = field(default_factory=list)
    blacklists: List[GenomeEntry] = field(default_factory=list)

    def add_genome(self, name: str, fasta_path: str, role: str,
                   penalty_weight: Optional[float] = None):
        """
        Add a genome to the set.

        Args:
            name: Genome name
            fasta_path: Path to FASTA file
            role: "target", "background", or "blacklist"
            penalty_weight: Optional custom penalty weight
        """
        role_enum = GenomeRole(role.lower())

        entry = GenomeEntry(
            name=name,
            fasta_path=Path(fasta_path),
            role=role_enum,
            penalty_weight=penalty_weight if penalty_weight is not None else 1.0
        )

        if role_enum == GenomeRole.TARGET:
            self.targets.append(entry)
        elif role_enum == GenomeRole.BACKGROUND:
            self.backgrounds.append(entry)
        else:  # BLACKLIST
            self.blacklists.append(entry)

    def get_all_genomes(self) -> List[GenomeEntry]:
        """Get all genomes in a single list"""
        return self.targets + self.backgrounds + self.blacklists

    def get_non_targets(self) -> List[GenomeEntry]:
        """Get all background and blacklist genomes"""
        return self.backgrounds + self.blacklists

    def validate(self) -> bool:
        """
        Validate genome set.

        Returns:
            True if valid, raises ValueError otherwise
        """
        if len(self.targets) == 0:
            raise ValueError("At least one target genome required")

        # Check for duplicate names
        all_names = [g.name for g in self.get_all_genomes()]
        if len(all_names) != len(set(all_names)):
            raise ValueError("Duplicate genome names detected")

        # Check all files exist
        for genome in self.get_all_genomes():
            if not genome.fasta_path.exists():
                raise FileNotFoundError(f"Genome not found: {genome.fasta_path}")

        return True

    def summary(self) -> str:
        """Get summary of genome set"""
        lines = []
        lines.append("="*80)
        lines.append("MULTI-GENOME SWGA CONFIGURATION")
        lines.append("="*80)

        lines.append(f"\nTarget Genomes ({len(self.targets)}):")
        for g in self.targets:
            lines.append(f"  • {g.name}")
            lines.append(f"    Path: {g.fasta_path}")
            if g.gc_content:
                lines.append(f"    GC: {g.gc_content:.1%}")
            if g.size:
                lines.append(f"    Size: {g.size:,} bp")

        if self.backgrounds:
            lines.append(f"\nBackground Genomes ({len(self.backgrounds)}):")
            for g in self.backgrounds:
                lines.append(f"  • {g.name} (penalty: {g.penalty_weight}x)")
                lines.append(f"    Path: {g.fasta_path}")

        if self.blacklists:
            lines.append(f"\nBlacklist Genomes ({len(self.blacklists)}):")
            for g in self.blacklists:
                lines.append(f"  • {g.name} (penalty: {g.penalty_weight}x)")
                lines.append(f"    Path: {g.fasta_path}")

        lines.append("="*80)
        return "\n".join(lines)


@dataclass
class MultiGenomeScore:
    """
    Score for a primer across multiple genomes.

    Attributes:
        primer: Primer sequence
        target_frequency: Binding frequency in target genome(s)
        background_frequency: Binding frequency in background genome(s)
        blacklist_frequency: Binding frequency in blacklist genome(s)
        enrichment_score: Target/background ratio
        penalty_score: Weighted penalty from non-target binding
        passes: Whether primer meets all criteria
        details: Per-genome binding frequencies
    """
    primer: str
    target_frequency: float
    background_frequency: float
    blacklist_frequency: float
    enrichment_score: float
    penalty_score: float
    passes: bool
    details: Dict[str, float] = field(default_factory=dict)

    def __str__(self):
        return f"""Primer: {self.primer}
  Target frequency: {self.target_frequency:.2e}
  Background frequency: {self.background_frequency:.2e}
  Blacklist frequency: {self.blacklist_frequency:.2e}
  Enrichment: {self.enrichment_score:.1f}x
  Penalty: {self.penalty_score:.3f}
  Passes: {self.passes}"""


class MultiGenomeFilter:
    """
    Filter primers based on binding to multiple genomes with differential penalties.

    Strategy:
    1. Calculate binding frequency in each genome
    2. Apply role-based penalty weights
    3. Compute enrichment score (target/background ratio)
    4. Filter based on thresholds with genome-specific criteria

    Key parameters:
    - min_target_freq: Minimum binding in target (e.g., 1e-5)
    - max_background_freq: Maximum binding in background (e.g., 1e-4)
    - max_blacklist_freq: Maximum binding in blacklist (e.g., 1e-6, much stricter)
    - min_enrichment: Minimum target/background ratio (e.g., 10x)
    """

    def __init__(self,
                 genome_set: GenomeSet,
                 min_target_freq: float = 1e-5,
                 max_background_freq: float = 1e-4,
                 max_blacklist_freq: float = 1e-6,
                 min_enrichment: float = 10.0,
                 use_gini: bool = True,
                 max_gini: float = 0.6):
        """
        Initialize multi-genome filter.

        Args:
            genome_set: Collection of target/background/blacklist genomes
            min_target_freq: Minimum frequency in target genome(s)
            max_background_freq: Maximum frequency in background genome(s)
            max_blacklist_freq: Maximum frequency in blacklist genome(s) (stricter)
            min_enrichment: Minimum enrichment ratio (target/background)
            use_gini: Use Gini index for binding uniformity
            max_gini: Maximum Gini index (0 = perfect uniformity, 1 = all in one site)
        """
        self.genome_set = genome_set
        self.genome_set.validate()

        self.min_target_freq = min_target_freq
        self.max_background_freq = max_background_freq
        self.max_blacklist_freq = max_blacklist_freq
        self.min_enrichment = min_enrichment
        self.use_gini = use_gini
        self.max_gini = max_gini

        # Cache for k-mer counts and positions
        self.kmer_counts: Dict[str, Dict[str, int]] = {}  # {genome_name: {kmer: count}}
        self.genome_sizes: Dict[str, int] = {}

        logger.info("Initialized multi-genome filter")
        logger.info(f"  Targets: {len(genome_set.targets)}")
        logger.info(f"  Backgrounds: {len(genome_set.backgrounds)}")
        logger.info(f"  Blacklists: {len(genome_set.blacklists)}")
        logger.info(f"  Min target freq: {min_target_freq:.2e}")
        logger.info(f"  Max background freq: {max_background_freq:.2e}")
        logger.info(f"  Max blacklist freq: {max_blacklist_freq:.2e}")
        logger.info(f"  Min enrichment: {min_enrichment}x")

    def load_genome_counts(self, genome_name: str, kmer_counts: Dict[str, int],
                          genome_size: int):
        """
        Load pre-computed k-mer counts for a genome.

        Args:
            genome_name: Name of genome
            kmer_counts: Dictionary mapping k-mer to count
            genome_size: Total genome size in bp
        """
        self.kmer_counts[genome_name] = kmer_counts
        self.genome_sizes[genome_name] = genome_size

        logger.info(f"Loaded counts for {genome_name}: {len(kmer_counts)} k-mers, "
                   f"{genome_size:,} bp")

    def _calculate_frequency(self, primer: str, genome_name: str) -> float:
        """
        Calculate binding frequency for a primer in a genome.

        Args:
            primer: Primer sequence
            genome_name: Name of genome

        Returns:
            Binding frequency (count / genome_size)
        """
        if genome_name not in self.kmer_counts:
            logger.warning(f"No counts loaded for {genome_name}")
            return 0.0

        count = self.kmer_counts[genome_name].get(primer, 0)
        size = self.genome_sizes[genome_name]

        return count / size if size > 0 else 0.0

    def score_primer(self, primer: str) -> MultiGenomeScore:
        """
        Score a primer across all genomes.

        Args:
            primer: Primer sequence

        Returns:
            MultiGenomeScore with frequencies and pass/fail
        """
        details = {}

        # Calculate frequencies for each genome
        target_freqs = []
        for genome in self.genome_set.targets:
            freq = self._calculate_frequency(primer, genome.name)
            target_freqs.append(freq)
            details[genome.name] = freq

        background_freqs = []
        for genome in self.genome_set.backgrounds:
            freq = self._calculate_frequency(primer, genome.name)
            background_freqs.append(freq)
            details[genome.name] = freq

        blacklist_freqs = []
        for genome in self.genome_set.blacklists:
            freq = self._calculate_frequency(primer, genome.name)
            blacklist_freqs.append(freq)
            details[genome.name] = freq

        # Aggregate frequencies
        target_freq = np.mean(target_freqs) if target_freqs else 0.0
        background_freq = np.max(background_freqs) if background_freqs else 0.0
        blacklist_freq = np.max(blacklist_freqs) if blacklist_freqs else 0.0

        # Calculate enrichment
        if background_freq > 0:
            enrichment = target_freq / background_freq
        else:
            enrichment = float('inf') if target_freq > 0 else 0.0

        # Calculate weighted penalty score
        penalty = 0.0

        # Background penalty
        for genome, freq in zip(self.genome_set.backgrounds, background_freqs):
            penalty += freq * genome.penalty_weight

        # Blacklist penalty (much higher weight)
        for genome, freq in zip(self.genome_set.blacklists, blacklist_freqs):
            penalty += freq * genome.penalty_weight

        # Apply filters
        passes = (
            target_freq >= self.min_target_freq and
            background_freq <= self.max_background_freq and
            blacklist_freq <= self.max_blacklist_freq and
            enrichment >= self.min_enrichment
        )

        return MultiGenomeScore(
            primer=primer,
            target_frequency=target_freq,
            background_frequency=background_freq,
            blacklist_frequency=blacklist_freq,
            enrichment_score=enrichment,
            penalty_score=penalty,
            passes=passes,
            details=details
        )

    def filter_primers(self, candidates: List[str],
                      verbose: bool = False) -> Tuple[List[str], List[MultiGenomeScore]]:
        """
        Filter primer candidates across all genomes.

        Args:
            candidates: List of primer sequences
            verbose: Print detailed progress

        Returns:
            Tuple of (passing_primers, all_scores)
        """
        scores = []
        passing = []

        for i, primer in enumerate(candidates):
            score = self.score_primer(primer)
            scores.append(score)

            if score.passes:
                passing.append(primer)

            if verbose and (i + 1) % 1000 == 0:
                logger.info(f"Processed {i+1}/{len(candidates)} primers, "
                           f"{len(passing)} passing ({len(passing)/(i+1)*100:.1f}%)")

        logger.info(f"Multi-genome filtering complete:")
        logger.info(f"  Input: {len(candidates)} primers")
        logger.info(f"  Passing: {len(passing)} primers ({len(passing)/len(candidates)*100:.1f}%)")

        # Summary statistics
        if scores:
            target_freqs = [s.target_frequency for s in scores if s.passes]
            background_freqs = [s.background_frequency for s in scores if s.passes]
            blacklist_freqs = [s.blacklist_frequency for s in scores if s.passes]

            if target_freqs:
                logger.info(f"  Target frequency: {np.mean(target_freqs):.2e} ± "
                           f"{np.std(target_freqs):.2e}")
            if background_freqs:
                logger.info(f"  Background frequency: {np.mean(background_freqs):.2e} ± "
                           f"{np.std(background_freqs):.2e}")
            if blacklist_freqs and any(f > 0 for f in blacklist_freqs):
                logger.info(f"  Blacklist frequency: {np.mean(blacklist_freqs):.2e} ± "
                           f"{np.std(blacklist_freqs):.2e}")

        return passing, scores

    def rank_primers(self, primers: List[str],
                    top_n: Optional[int] = None) -> List[Tuple[str, MultiGenomeScore]]:
        """
        Rank primers by composite score.

        Composite score = target_freq * enrichment / penalty

        Args:
            primers: List of primer sequences
            top_n: Return only top N primers

        Returns:
            List of (primer, score) tuples sorted by composite score
        """
        scored = []

        for primer in primers:
            score = self.score_primer(primer)

            # Composite score: maximize target, maximize enrichment, minimize penalty
            if score.penalty_score > 0:
                composite = (score.target_frequency * score.enrichment_score) / score.penalty_score
            else:
                composite = score.target_frequency * score.enrichment_score * 1e6

            scored.append((primer, score, composite))

        # Sort by composite score descending
        scored.sort(key=lambda x: x[2], reverse=True)

        if top_n:
            scored = scored[:top_n]

        return [(p, s) for p, s, c in scored]

    def generate_report(self, primers: List[str], output_file: Optional[str] = None) -> str:
        """
        Generate detailed multi-genome report.

        Args:
            primers: List of primers to report on
            output_file: Optional file to write report to

        Returns:
            Report string
        """
        lines = []
        lines.append("="*80)
        lines.append("MULTI-GENOME SWGA REPORT")
        lines.append("="*80)

        # Configuration
        lines.append("\nGenome Configuration:")
        lines.append(f"  Target genomes: {len(self.genome_set.targets)}")
        for g in self.genome_set.targets:
            lines.append(f"    • {g.name}")

        if self.genome_set.backgrounds:
            lines.append(f"  Background genomes: {len(self.genome_set.backgrounds)}")
            for g in self.genome_set.backgrounds:
                lines.append(f"    • {g.name} (penalty: {g.penalty_weight}x)")

        if self.genome_set.blacklists:
            lines.append(f"  Blacklist genomes: {len(self.genome_set.blacklists)}")
            for g in self.genome_set.blacklists:
                lines.append(f"    • {g.name} (penalty: {g.penalty_weight}x)")

        # Filtering criteria
        lines.append("\nFiltering Criteria:")
        lines.append(f"  Min target frequency: {self.min_target_freq:.2e}")
        lines.append(f"  Max background frequency: {self.max_background_freq:.2e}")
        lines.append(f"  Max blacklist frequency: {self.max_blacklist_freq:.2e}")
        lines.append(f"  Min enrichment: {self.min_enrichment}x")

        # Primer details
        lines.append(f"\nSelected Primers ({len(primers)}):")
        lines.append("="*80)

        for i, primer in enumerate(primers, 1):
            score = self.score_primer(primer)
            lines.append(f"\nPrimer {i}: {primer}")
            lines.append(f"  Target frequency: {score.target_frequency:.2e}")
            lines.append(f"  Background frequency: {score.background_frequency:.2e}")
            lines.append(f"  Blacklist frequency: {score.blacklist_frequency:.2e}")
            lines.append(f"  Enrichment: {score.enrichment_score:.1f}x")
            lines.append(f"  Penalty score: {score.penalty_score:.3e}")

            # Per-genome breakdown
            lines.append(f"  Binding frequencies:")
            for genome_name, freq in score.details.items():
                lines.append(f"    {genome_name:30s}: {freq:.2e}")

        lines.append("\n" + "="*80)

        report = "\n".join(lines)

        if output_file:
            with open(output_file, 'w') as f:
                f.write(report)
            logger.info(f"Report written to {output_file}")

        return report


if __name__ == '__main__':
    # Example usage
    logging.basicConfig(level=logging.INFO)

    print("Multi-Genome Filter - Example Usage\n")
    print("Demonstrates filtering primers for Borrelia detection in tick")
    print("while avoiding tick DNA and other tick-borne bacteria\n")

    # Create genome set
    genome_set = GenomeSet()

    # Add target: Borrelia burgdorferi (Lyme disease)
    genome_set.add_genome(
        name="Borrelia_burgdorferi",
        fasta_path="/path/to/borrelia.fasta",
        role="target"
    )

    # Add background: Ixodes tick (host)
    genome_set.add_genome(
        name="Ixodes_scapularis",
        fasta_path="/path/to/tick.fasta",
        role="background",
        penalty_weight=1.0  # Standard avoidance
    )

    # Add blacklist: Rickettsia (co-existing pathogen)
    genome_set.add_genome(
        name="Rickettsia_rickettsii",
        fasta_path="/path/to/rickettsia.fasta",
        role="blacklist",
        penalty_weight=5.0  # Strong avoidance
    )

    print(genome_set.summary())

    # Create filter
    mfilter = MultiGenomeFilter(
        genome_set=genome_set,
        min_target_freq=1e-5,      # Must bind Borrelia
        max_background_freq=1e-4,   # Tolerate some tick binding
        max_blacklist_freq=1e-6,    # Minimal Rickettsia binding
        min_enrichment=10.0         # 10x more Borrelia than tick
    )

    print("\nFilter configured successfully")
    print("Ready to process primers with multi-genome filtering")
