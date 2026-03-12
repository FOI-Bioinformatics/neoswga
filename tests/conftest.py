"""Shared fixtures for the neoswga test suite.

Centralizes commonly duplicated test data: primer sequences, genome strings,
reaction conditions, mock caches, and temporary FASTA files.
"""

import numpy as np
import pytest
from pathlib import Path

from neoswga.core.reaction_conditions import ReactionConditions


# ---------------------------------------------------------------------------
# Primer sequences
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_primers():
    """List of 8-mer primer sequences spanning different GC contents."""
    return [
        "ATCGATCG",   # 50% GC, balanced
        "GCTAGCTA",   # 50% GC, balanced
        "AATTCCGG",   # 50% GC, clustered
        "TTAACCGG",   # 50% GC, clustered
        "GGCCAATT",   # 50% GC, clustered
    ]


@pytest.fixture
def diverse_primers():
    """Primers spanning a wider range of GC content and lengths."""
    return [
        "ATCGATCG",       # 50% GC, 8-mer
        "GCGCGCGC",       # 100% GC, 8-mer
        "AAAATTTT",       # 0% GC, 8-mer
        "ATCGATCGATCG",   # 50% GC, 12-mer
        "GCTAGCTAGC",     # 50% GC, 10-mer
    ]


# ---------------------------------------------------------------------------
# Genome sequences
# ---------------------------------------------------------------------------

@pytest.fixture
def simple_genome():
    """A short repeating genome string (12 bp)."""
    return "ATCGATCGATCG"


@pytest.fixture
def repeat_genome():
    """A homopolymer genome string (10 bp)."""
    return "AAAAAAAAAA"


# ---------------------------------------------------------------------------
# Reaction conditions
# ---------------------------------------------------------------------------

@pytest.fixture
def phi29_conditions():
    """Standard phi29 reaction conditions (30 C)."""
    return ReactionConditions(temp=30.0, polymerase="phi29", mg_conc=2.5)


@pytest.fixture
def equiphi29_conditions():
    """EquiPhi29 reaction conditions (42 C) with DMSO and betaine."""
    return ReactionConditions(
        temp=42.0,
        polymerase="equiphi29",
        mg_conc=2.5,
        dmso_percent=5.0,
        betaine_m=1.0,
    )


# ---------------------------------------------------------------------------
# Mock position cache
# ---------------------------------------------------------------------------

class MockPositionCache:
    """Lightweight mock for PositionCache used across optimizer tests.

    Returns deterministic positions based on an optional positions dict,
    falling back to evenly spaced sites derived from primer length.
    """

    def __init__(self, positions_dict=None, genome_length=100_000):
        self._positions = positions_dict or {}
        self._genome_length = genome_length

    def get_positions(self, prefix, primer, strand="both"):
        key = (prefix, primer, strand)
        if key in self._positions:
            return np.asarray(self._positions[key])
        # Default: shorter primers bind more frequently
        n_sites = max(5, 20 - len(primer))
        spacing = self._genome_length // (n_sites + 1)
        return np.arange(spacing, spacing * (n_sites + 1), spacing)


@pytest.fixture
def mock_position_cache():
    """A MockPositionCache with default settings."""
    return MockPositionCache()


# ---------------------------------------------------------------------------
# Temporary FASTA files
# ---------------------------------------------------------------------------

@pytest.fixture
def tmp_fasta_file(tmp_path):
    """Create a minimal single-sequence FASTA file. Returns the Path."""
    fasta = tmp_path / "test_genome.fasta"
    fasta.write_text(">seq1\nATCGATCGATCGATCGATCG\nGCTAGCTAGCTAGCTAGCTA\n")
    return fasta


@pytest.fixture
def tmp_fasta_pair(tmp_path):
    """Create a foreground and background FASTA pair. Returns (fg, bg) Paths."""
    fg = tmp_path / "target.fasta"
    fg.write_text(">target_seq1\nATCGATCGATCGATCGATCG\nGCTAGCTAGCTAGCTAGCTA\n")
    bg = tmp_path / "background.fasta"
    bg.write_text(">bg_seq1\nAAAACCCCGGGGTTTTAAAA\nCCCCGGGGTTTTAAAACCCC\n")
    return fg, bg
