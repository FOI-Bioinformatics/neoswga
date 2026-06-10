"""Phase 2 hardening tests: determinism, strand adapters, column resolver."""

import pandas as pd
import pytest


# ---------------------------------------------------------------------------
# RF k-mer sampling determinism
# ---------------------------------------------------------------------------

def test_kmer_sampling_seed_is_reproducible():
    """Same seed -> identical sampling decisions; different seed differs."""
    from neoswga.core import rf_preprocessing as rf

    def draws(seed, n=200):
        rf.set_kmer_sampling_seed(seed)
        return [rf._RNG.random() for _ in range(n)]

    a = draws(42)
    b = draws(42)
    c = draws(43)
    assert a == b           # reproducible
    assert a != c           # seed actually matters


def test_set_kmer_sampling_seed_none_reseeds_from_entropy():
    """Passing None must not raise and yields a usable RNG."""
    from neoswga.core import rf_preprocessing as rf

    rf.set_kmer_sampling_seed(None)
    val = rf._RNG.random()
    assert 0.0 <= val <= 1.0


# ---------------------------------------------------------------------------
# Strand convention adapters
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("raw,expected", [
    ("+", "forward"), ("-", "reverse"),
    ("forward", "forward"), ("reverse", "reverse"), ("both", "both"),
    ("F", "forward"), ("r", "reverse"), ("FWD", "forward"),
])
def test_to_cache_strand(raw, expected):
    from neoswga.core.strand_conventions import to_cache_strand
    assert to_cache_strand(raw) == expected


def test_to_cache_strand_rejects_unknown():
    from neoswga.core.strand_conventions import to_cache_strand
    with pytest.raises(ValueError):
        to_cache_strand("sideways")


@pytest.mark.parametrize("raw,expected", [
    ("forward", "+"), ("+", "+"),
    ("reverse", "-"), ("-", "-"),
    ("both", "."),
])
def test_to_bed_strand(raw, expected):
    from neoswga.core.strand_conventions import to_bed_strand
    assert to_bed_strand(raw) == expected


@pytest.mark.parametrize("pos,expected", [
    (10, "forward"), (0, "forward"), (-1, "reverse"), (-9999, "reverse"),
])
def test_from_sign(pos, expected):
    from neoswga.core.strand_conventions import from_sign
    assert from_sign(pos) == expected


# ---------------------------------------------------------------------------
# primer_column resolver
# ---------------------------------------------------------------------------

def test_primer_column_prefers_primer():
    from neoswga.core.io_utils import primer_column
    df = pd.DataFrame({"primer": ["ATCG"], "seq": ["GGGG"]})
    assert primer_column(df) == "primer"


def test_primer_column_falls_back_to_seq():
    from neoswga.core.io_utils import primer_column
    df = pd.DataFrame({"seq": ["ATCG"]})
    assert primer_column(df) == "seq"


def test_primer_column_raises_when_neither():
    from neoswga.core.io_utils import primer_column
    df = pd.DataFrame({"other": [1]})
    with pytest.raises(KeyError):
        primer_column(df)
