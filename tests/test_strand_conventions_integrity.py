"""Strand handling across the codebase.

Three strand vocabularies coexist here, which `strand_conventions` documents:

    PositionCache          'forward' / 'reverse' / 'both'
    strand_bias_analyzer   '+' / '-'
    amplicon_network       sign-encoded positions (negative = reverse)

Underneath all of them is one storage convention: an HDF5 position file holds
forward hits under the primer, and reverse-strand hits under its REVERSE
COMPLEMENT -- because that is where the revcomp matches the forward strand.

Getting this wrong has been the single most productive source of bugs in this
audit. Four modules independently read only the forward strand and reported the
result as a measurement: PositionCache's scan fallback, the Bloom background
filter, SwgaSimulator's loader, and AmpliconNetwork. A primer binds
double-stranded DNA, so forward-only means roughly half of every count.

These tests fix the storage convention in place and check each consumer against
it.
"""

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

from neoswga.core.strand_conventions import from_sign, to_bed_strand, to_cache_strand
from neoswga.core.thermodynamics import reverse_complement

PRIMER = "TTGACCATGA"
FORWARD_SITES = [1_000, 4_000]
REVERSE_SITES = [2_500, 7_000]
GENOME_LENGTH = 10_000


@pytest.fixture
def position_file(tmp_path):
    """An HDF5 in the layout the pipeline writes."""
    prefix = str(tmp_path / "g")
    with h5py.File(f"{prefix}_{len(PRIMER)}mer_positions.h5", "w") as f:
        f.create_dataset(PRIMER, data=np.array(FORWARD_SITES, dtype=np.int32))
        f.create_dataset(reverse_complement(PRIMER), data=np.array(REVERSE_SITES, dtype=np.int32))
    return prefix


# ----------------------------------------------------------------------
# The adapters
# ----------------------------------------------------------------------


@pytest.mark.parametrize(
    "label,expected",
    [
        ("forward", "forward"),
        ("reverse", "reverse"),
        ("both", "both"),
        ("+", "forward"),
        ("-", "reverse"),
        ("F", "forward"),
        ("REV", "reverse"),
    ],
)
def test_cache_strand_accepts_every_vocabulary(label, expected):
    assert to_cache_strand(label) == expected


def test_cache_strand_rejects_nonsense():
    """Silently defaulting is what makes these bugs invisible."""
    with pytest.raises(ValueError):
        to_cache_strand("sideways")


@pytest.mark.parametrize(
    "label,expected", [("forward", "+"), ("reverse", "-"), ("both", "."), ("+", "+")]
)
def test_bed_strand_mapping(label, expected):
    """'both' is unstranded in BED, not reverse."""
    assert to_bed_strand(label) == expected


def test_sign_encoding_cannot_represent_a_reverse_hit_at_zero():
    """A real limitation of the sign convention, worth stating.

    `-0 == 0` in Python, so a reverse-strand site at position 0 decodes as
    forward. It is an edge case rather than a live bug -- position 0 requires a
    match at the very first base -- but it is a reason to prefer the
    forward/reverse vocabulary in new code.
    """
    assert from_sign(-5) == "reverse"
    assert from_sign(5) == "forward"
    assert from_sign(-0) == "forward"  # not "reverse": -0 is 0


# ----------------------------------------------------------------------
# The storage convention every consumer must agree on
# ----------------------------------------------------------------------


def test_position_cache_reads_both_strands(position_file):
    from neoswga.core.position_cache import PositionCache

    cache = PositionCache([position_file], [PRIMER])

    assert sorted(cache.get_positions(position_file, PRIMER, "forward")) == FORWARD_SITES
    assert sorted(cache.get_positions(position_file, PRIMER, "reverse")) == REVERSE_SITES
    assert sorted(cache.get_positions(position_file, PRIMER, "both")) == sorted(
        FORWARD_SITES + REVERSE_SITES
    )


def test_position_cache_rejects_the_wrong_vocabulary(position_file):
    """Passing '+' returned an empty array -- a confident zero.

    That is indistinguishable from "this primer binds nowhere", and it is the
    exact trap CLAUDE.md warns about in prose. Prose is not a guard.
    """
    from neoswga.core.position_cache import PositionCache

    cache = PositionCache([position_file], [PRIMER])

    with pytest.raises(ValueError, match="forward"):
        cache.get_positions(position_file, PRIMER, "+")


def test_primer_attributes_counts_both_strands(position_file):
    from neoswga.core import primer_attributes as pa

    assert pa.get_rate_from_h5py(PRIMER, [position_file]) == len(FORWARD_SITES) + len(REVERSE_SITES)


def test_amplicon_network_sees_both_strands(position_file):
    """It read only the primer key and split it by sign.

    Real position files hold no negative values, so `positions_reverse` was
    always empty -- and an amplicon edge is "forward primer -> downstream
    reverse primer", so the network could never form a single edge. Every metric
    derived from it was computed on an edgeless graph.
    """
    from neoswga.core.amplicon_network import AmpliconNetwork

    net = AmpliconNetwork(primers=[PRIMER], genome_length=GENOME_LENGTH)
    net.load_positions_from_hdf5(position_file)

    assert sorted(net.positions_forward[PRIMER]) == FORWARD_SITES
    assert sorted(net.positions_reverse[PRIMER]) == REVERSE_SITES


def test_amplicon_network_forms_edges(position_file):
    """The whole point of the network. It had zero edges for every input."""
    from neoswga.core.amplicon_network import AmpliconNetwork

    net = AmpliconNetwork(primers=[PRIMER], genome_length=GENOME_LENGTH)
    net.load_positions_from_hdf5(position_file)
    net.build_network()

    assert net.G.number_of_nodes() == len(FORWARD_SITES) + len(REVERSE_SITES)
    assert net.G.number_of_edges() > 0, "no amplicon edges formed"


def test_amplicon_network_still_reads_sign_encoded_files(tmp_path):
    """The documented sign convention is honoured where a file uses it."""
    from neoswga.core.amplicon_network import AmpliconNetwork

    prefix = str(tmp_path / "signed")
    with h5py.File(f"{prefix}_{len(PRIMER)}mer_positions.h5", "w") as f:
        signed = FORWARD_SITES + [-p for p in REVERSE_SITES]
        f.create_dataset(PRIMER, data=np.array(signed, dtype=np.int32))

    net = AmpliconNetwork(primers=[PRIMER], genome_length=GENOME_LENGTH)
    net.load_positions_from_hdf5(prefix)

    assert sorted(net.positions_forward[PRIMER]) == FORWARD_SITES
    assert sorted(net.positions_reverse[PRIMER]) == REVERSE_SITES


def test_swga_simulator_sees_both_strands(position_file, tmp_path):
    """Fixed earlier in this audit; kept here so all four sit together."""
    from neoswga.core.swga_simulator import SwgaSimulator

    fasta = tmp_path / "g.fasta"
    fasta.write_text(">g\n" + "A" * GENOME_LENGTH + "\n")

    sim = SwgaSimulator(
        primers=[PRIMER],
        fg_genome=str(fasta),
        bg_genome=str(fasta),
        fg_positions_h5=f"{position_file}_{len(PRIMER)}mer_positions.h5",
        bg_positions_h5=f"{position_file}_{len(PRIMER)}mer_positions.h5",
        bin_size=1_000,
    )

    assert sorted(sim.fg_positions[PRIMER]["+"]) == FORWARD_SITES
    assert sorted(sim.fg_positions[PRIMER]["-"]) == REVERSE_SITES


def test_background_bloom_filter_matches_either_strand(tmp_path):
    """Also fixed earlier: a primer whose revcomp is in the background binds it."""
    pytest.importorskip("pybloom_live")
    import random

    from neoswga.core.background_filter import BackgroundBloomFilter

    rng = random.Random(11)
    sequence = "".join(rng.choice("ACGT") for _ in range(4_000))
    fasta = tmp_path / "bg.fasta"
    fasta.write_text(">bg\n" + sequence + "\n")

    bloom = BackgroundBloomFilter(capacity=100_000, error_rate=0.01)
    bloom.add_genome(str(fasta), min_k=10, max_k=10)

    probe = sequence[500:510]
    assert bloom.contains(probe)
    assert bloom.contains(reverse_complement(probe)), "reverse strand not matched"


# ----------------------------------------------------------------------
# Export
# ----------------------------------------------------------------------


def test_bed_export_uses_the_adapter(tmp_path):
    """`"+" if strand == "forward" else "-"` mapped everything else to reverse,
    including "both", so a BED file could carry confidently wrong orientations.
    """
    from neoswga.core.export import export_to_bed

    out = tmp_path / "sites.bed"
    positions = {PRIMER: [(1_000, "forward"), (2_500, "reverse"), (4_000, "both")]}
    export_to_bed([PRIMER], positions, "chr1", str(out))

    strands = [line.split("\t")[5].strip() for line in out.read_text().splitlines() if line.strip()]
    assert strands == ["+", "-", "."]
