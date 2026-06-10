"""Shared fixtures for the neoswga test suite.

Centralizes commonly duplicated test data: primer sequences, genome strings,
reaction conditions, mock caches, and temporary FASTA files.
"""

import glob
import os
from pathlib import Path

import numpy as np
import pytest

from neoswga.core.reaction_conditions import ReactionConditions

_EXAMPLE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "examples", "plasmid_example"
)


@pytest.fixture(scope="session", autouse=True)
def _prime_plasmid_example():
    """Prime the plasmid example with pipeline outputs once per session.

    Several tests depend on generated files in ``examples/plasmid_example/`` that
    are NOT committed (they are gitignored build artifacts):

    - Integration tests copy the directory (incl. the ``*_Nmer_all.txt`` k-mer
      files) into a tmpdir.
    - The swap-primer / contract-set / rescore-set CLI smoke tests run the CLI
      with ``cwd=examples/plasmid_example`` and need ``step2_df.csv`` /
      ``step3_df.csv`` and the ``*_positions.h5`` files.

    On a clean checkout / CI these do not exist, so those tests fail. This
    session-scoped autouse fixture runs ``count-kmers`` -> ``filter`` ->
    ``score`` once in the example directory to generate them (the default 6-12
    k range covers every range the tests use). It is a no-op when the outputs
    already exist (a dev's primed directory) or when jellyfish is unavailable
    (the pipeline's hard dependency); dependent tests then skip or surface the
    missing dependency themselves.
    """
    if not os.path.isdir(_EXAMPLE_DIR):
        return

    have_kmers = bool(glob.glob(os.path.join(_EXAMPLE_DIR, "*mer_all.txt")))
    have_step3 = os.path.exists(os.path.join(_EXAMPLE_DIR, "step3_df.csv"))
    if have_kmers and have_step3:
        return  # already primed

    try:
        from neoswga.core.kmer_counter import check_jellyfish_available
    except Exception:
        return
    if not check_jellyfish_available():
        return

    import neoswga.core.pipeline as pipeline_mod
    from neoswga.core import parameter

    def _reset():
        pipeline_mod._initialized = False
        for attr in ("fg_prefixes", "bg_prefixes", "fg_genomes", "bg_genomes",
                     "fg_seq_lengths", "bg_seq_lengths", "fg_circular", "bg_circular"):
            setattr(pipeline_mod, attr, None)
        parameter.json_file = "params.json"

    cwd = os.getcwd()
    try:
        os.chdir(_EXAMPLE_DIR)
        _reset(); pipeline_mod._initialize(); pipeline_mod.step1()
        _reset(); pipeline_mod.step2()
        _reset(); pipeline_mod.step3()
    except Exception:
        # Best-effort priming; dependent tests will report any real problem.
        pass
    finally:
        _reset()
        os.chdir(cwd)


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
