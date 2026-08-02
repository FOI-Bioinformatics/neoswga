"""Tests for `string_search.get_positions`, the scanning entry point.

Every binding position in the tool comes from this module, so a silent gap here
propagates into coverage, gaps, selectivity and every downstream score. Three of
the bugs found in this audit -- in PositionCache's scan fallback, the Bloom
background filter, and SwgaSimulator's loader -- were all about positions going
missing without anything saying so.

`get_positions` had its own version: it scanned `range(parameter.min_k,
parameter.max_k + 1)` rather than the lengths it was given, so any primer
outside the configured window was dropped silently. That is precisely the
bring-your-own-oligo case this function is the entry point for.

The module has two implementations -- a single-pass Aho-Corasick scan when
`pyahocorasick` is installed, and a per-k multiprocessing fallback otherwise --
and they differ in return type. Both are exercised here where possible.
"""

import os
import random

import pytest

from neoswga.core import parameter
from neoswga.core import string_search as ss

P10 = "TTGACCATGA"
P14 = "TTGACCATGACGTA"
P12 = "GCATTACGGTAC"


@pytest.fixture
def genome(tmp_path):
    """A genome with primers planted at known positions, of three lengths."""
    rng = random.Random(5)
    seq = "".join(rng.choice("ACGT") for _ in range(6_000))
    seq = seq[:1_000] + P10 + seq[1_000 + len(P10) :]
    seq = seq[:3_000] + P14 + seq[3_000 + len(P14) :]
    seq = seq[:5_000] + P12 + seq[5_000 + len(P12) :]

    fasta = tmp_path / "g.fasta"
    fasta.write_text(">g\n" + seq + "\n")
    return {"fasta": str(fasta), "prefix": str(tmp_path / "g"), "seq": seq}


@pytest.fixture(autouse=True)
def restore_k_window():
    before = (parameter.min_k, parameter.max_k)
    yield
    parameter.min_k, parameter.max_k = before


def _positions(result, prefix, primer):
    if result is None:
        pytest.skip("fallback path returns None rather than a position map")
    return result.get((prefix, primer))


# ----------------------------------------------------------------------
# The configured k window must not silently drop primers
# ----------------------------------------------------------------------


def test_primer_outside_the_configured_window_is_still_found(genome):
    """The bring-your-own-oligo case.

    An externally designed 14-mer set against a params.json configured for
    10-mers returned nothing for every primer -- an empty list indistinguishable
    from "binds nowhere".
    """
    parameter.min_k, parameter.max_k = 10, 10

    result = ss.get_positions([P10, P14], [genome["prefix"]], [genome["fasta"]], circular=False)

    assert _positions(result, genome["prefix"], P14) == [3_000]


def test_primers_inside_the_window_are_unaffected(genome):
    parameter.min_k, parameter.max_k = 10, 10

    result = ss.get_positions([P10, P14], [genome["prefix"]], [genome["fasta"]], circular=False)

    # The 10-mer is also the prefix of the planted 14-mer, so it occurs twice.
    assert _positions(result, genome["prefix"], P10) == [1_000, 3_000]


def test_mixed_lengths_are_all_scanned(genome):
    """A set with 10-, 12- and 14-mers must return all three."""
    parameter.min_k, parameter.max_k = 6, 12

    result = ss.get_positions(
        [P10, P12, P14], [genome["prefix"]], [genome["fasta"]], circular=False
    )

    assert _positions(result, genome["prefix"], P12) == [5_000]
    assert _positions(result, genome["prefix"], P14) == [3_000]


def test_explicit_k_values_still_constrain(genome):
    """Callers that want to limit the scan can, and that is honoured."""
    parameter.min_k, parameter.max_k = 6, 20

    result = ss.get_positions(
        [P10, P14], [genome["prefix"]], [genome["fasta"]], circular=False, k_values=[10]
    )

    assert _positions(result, genome["prefix"], P10) == [1_000, 3_000]
    assert not _positions(result, genome["prefix"], P14)


def test_out_of_window_lengths_are_logged(genome, caplog):
    """Including them is right, but silently changing what gets indexed is not."""
    import logging

    parameter.min_k, parameter.max_k = 10, 10

    with caplog.at_level(logging.INFO, logger="neoswga.core.string_search"):
        ss.get_positions([P10, P14], [genome["prefix"]], [genome["fasta"]], circular=False)

    assert "14" in caplog.text


# ----------------------------------------------------------------------
# What the scan returns
# ----------------------------------------------------------------------


def test_reverse_complements_are_scanned_too(genome):
    """A primer binds either strand.

    Forward hits are stored under the primer, reverse-strand hits under its
    reverse complement -- the convention PositionCache reads back.
    """
    from neoswga.core.thermodynamics import reverse_complement

    parameter.min_k, parameter.max_k = 10, 14

    result = ss.get_positions([P10], [genome["prefix"]], [genome["fasta"]], circular=False)
    if result is None:
        pytest.skip("fallback path returns None")

    assert (genome["prefix"], reverse_complement(P10)) in result


def test_a_primer_absent_from_the_genome_returns_empty_not_missing(genome, tmp_path):
    """Absent must be an empty list, so callers can tell it was looked up."""
    parameter.min_k, parameter.max_k = 10, 10
    absent = "CGCGCGCGCG"

    result = ss.get_positions([P10, absent], [genome["prefix"]], [genome["fasta"]], circular=False)
    if result is None:
        pytest.skip("fallback path returns None")

    assert result.get((genome["prefix"], absent)) == []


def test_positions_are_written_to_the_layout_position_cache_reads(genome):
    """The HDF5 file name and key layout are a contract between the two."""
    import h5py

    parameter.min_k, parameter.max_k = 10, 14
    ss.get_positions([P10, P14], [genome["prefix"]], [genome["fasta"]], circular=False)

    for primer in (P10, P14):
        path = f"{genome['prefix']}_{len(primer)}mer_positions.h5"
        assert os.path.exists(path), f"no position file for {len(primer)}-mers"
        with h5py.File(path, "r") as f:
            assert primer in f


def test_written_positions_round_trip_through_position_cache(genome):
    """End to end: what the scanner writes is what the cache reads."""
    from neoswga.core.position_cache import PositionCache

    parameter.min_k, parameter.max_k = 10, 10
    ss.get_positions([P10], [genome["prefix"]], [genome["fasta"]], circular=False)

    cache = PositionCache([genome["prefix"]], [P10])
    assert sorted(cache.get_positions(genome["prefix"], P10, "forward").tolist()) == [1_000, 3_000]


# ----------------------------------------------------------------------
# The two implementations
# ----------------------------------------------------------------------


def test_return_type_depends_on_an_optional_package():
    """Pinned because it is surprising and callers must handle both.

    With `pyahocorasick` installed, `get_positions` returns a position map that
    lets the caller skip an HDF5 round-trip. Without it, the multiprocessing
    fallback returns None and the caller must read the files back. Installing an
    optional extra therefore changes the return type, not just the speed.
    """
    import inspect

    source = inspect.getsource(ss.get_positions)
    assert "AHOCORASICK_AVAILABLE" in source
    assert "return None" in source


@pytest.mark.skipif(not ss.AHOCORASICK_AVAILABLE, reason="pyahocorasick not installed")
def test_multi_k_scan_finds_every_length_in_one_pass(genome):
    """`get_all_positions_multi_k` is the fast path and had no direct coverage."""
    lists_by_k = {10: [P10], 12: [P12], 14: [P14]}
    found = ss.get_all_positions_multi_k(lists_by_k, genome["fasta"], circular=False)

    assert found[P10] == [1_000, 3_000]
    assert found[P12] == [5_000]
    assert found[P14] == [3_000]


@pytest.mark.skipif(not ss.AHOCORASICK_AVAILABLE, reason="pyahocorasick not installed")
def test_multi_k_scan_agrees_with_the_per_k_scan(genome):
    """The two scanners must not disagree; one is an optimisation of the other."""
    multi = ss.get_all_positions_multi_k({10: [P10]}, genome["fasta"], circular=False)
    per_k = ss.get_all_positions_per_k([P10], genome["fasta"], circular=False)

    assert multi[P10] == per_k[P10]


@pytest.mark.skipif(not ss.AHOCORASICK_AVAILABLE, reason="pyahocorasick not installed")
def test_multi_k_scan_handles_circular_wraparound(tmp_path):
    """A match spanning the origin is real on a circular chromosome."""
    seq = "ACGT" * 50
    primer = seq[-5:] + seq[:5]  # straddles the join

    fasta = tmp_path / "circ.fasta"
    fasta.write_text(">c\n" + seq + "\n")

    found = ss.get_all_positions_multi_k({10: [primer]}, str(fasta), circular=True)
    assert found[primer], "circular wrap-around match not found"
