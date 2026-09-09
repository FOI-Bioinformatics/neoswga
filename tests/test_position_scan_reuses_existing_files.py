"""Position files that already exist must not be scanned again.

`pipeline.step2` computed `position_files_exist` and fed it to nothing but a log
line reading "Reusing existing position files (incremental update only)". On the
Aho-Corasick branch the automaton was rebuilt from the whole primer list every
time, so re-running `filter` after a threshold tweak re-paid the entire scan
while announcing that it had not.
`check_which_primers_absent_in_h5py` was reached only from the fallback path.

The reuse has to keep the returned position map complete.
`primer_attributes.get_gini_from_txt_for_one_k` reads it with
`position_cache.get((prefix, primer), [])`, so a primer that is skipped and left
out of the map comes back with no positions, which the Gini gate turns into NaN
and drops. On a second run that would empty the entire candidate pool.
"""

import math
import random

import pytest

from neoswga.core import parameter
from neoswga.core import primer_attributes
from neoswga.core import string_search as ss

A12 = "GCATTACGGTAC"
B12 = "TTGACCATGACG"


@pytest.fixture
def genome(tmp_path):
    """Two 12-mers planted at known offsets in a fixed random genome.

    Seed 11 was checked: each primer occurs exactly once and neither reverse
    complement occurs at all, so the expected position lists are exact.
    """
    rng = random.Random(11)
    seq = "".join(rng.choice("ACGT") for _ in range(6_000))
    seq = seq[:1_000] + A12 + seq[1_000 + len(A12) :]
    seq = seq[:5_000] + B12 + seq[5_000 + len(B12) :]

    fasta = tmp_path / "g.fasta"
    fasta.write_text(">g\n" + seq + "\n")
    ss.clear_genome_cache()
    yield {"fasta": str(fasta), "prefix": str(tmp_path / "g")}
    ss.clear_genome_cache()


@pytest.fixture(autouse=True)
def k_window():
    before = (parameter.min_k, parameter.max_k)
    parameter.min_k, parameter.max_k = 6, 12
    yield
    parameter.min_k, parameter.max_k = before


def _record_scans(monkeypatch):
    """Record how many patterns each automaton build was handed."""
    calls = []
    real = ss.get_all_positions_multi_k

    def spy(primer_lists_by_k, seq_fname, circular, chunk_size=None):
        calls.append(sum(len(v) for v in primer_lists_by_k.values()))
        return real(primer_lists_by_k, seq_fname, circular, chunk_size)

    monkeypatch.setattr(ss, "get_all_positions_multi_k", spy)
    return calls


def test_a_repeated_scan_builds_no_automaton(genome, monkeypatch):
    calls = _record_scans(monkeypatch)

    ss.get_positions([A12, B12], [genome["prefix"]], [genome["fasta"]], circular=False)
    ss.get_positions([A12, B12], [genome["prefix"]], [genome["fasta"]], circular=False)

    assert calls == [4], (
        "the second call rebuilt the automaton from the full primer list; "
        f"pattern counts per scan were {calls} (two primers plus two reverse "
        "complements is one scan of 4)"
    )


def test_the_reused_result_is_identical_to_the_scanned_one(genome, monkeypatch):
    """The returned map must stay complete, or the Gini gate empties the pool."""
    first = ss.get_positions(
        [A12, B12], [genome["prefix"]], [genome["fasta"]], circular=False
    )
    second = ss.get_positions(
        [A12, B12], [genome["prefix"]], [genome["fasta"]], circular=False
    )

    assert second == first, (
        "the reused map differs from the scanned one; a primer missing from it "
        "reads downstream as a primer that binds nowhere"
    )
    assert first[(genome["prefix"], A12)] == [1000]
    assert first[(genome["prefix"], B12)] == [5000]


def test_only_the_primers_not_already_scanned_are_handed_to_the_automaton(
    genome, monkeypatch
):
    """A pool that grows by one primer costs one primer's worth of scanning."""
    ss.get_positions([A12], [genome["prefix"]], [genome["fasta"]], circular=False)

    calls = _record_scans(monkeypatch)
    result = ss.get_positions(
        [A12, B12], [genome["prefix"]], [genome["fasta"]], circular=False
    )

    assert calls == [2], (
        "expected one scan of B12 and its reverse complement only; pattern "
        f"counts per scan were {calls}"
    )
    assert result[(genome["prefix"], A12)] == [1000], "the cached primer was dropped"
    assert result[(genome["prefix"], B12)] == [5000]


def test_a_second_run_does_not_empty_the_candidate_pool(genome):
    """The failure this reuse could introduce, measured where it would be felt.

    A primer skipped by the reuse and left out of the returned map reaches
    `get_gini_from_txt_for_one_k` as a primer with no binding sites, which
    scores NaN, and `filter.get_gini` drops every NaN. On a second run every
    foreground primer is a skipped primer, so the whole pool would be dropped
    and `filter` would write an empty result while logging success. This test
    fails if the map comes back short, whatever the reason.
    """
    primers = [A12, B12]
    ss.get_positions(primers, [genome["prefix"]], [genome["fasta"]], circular=False)
    reused = ss.get_positions(primers, [genome["prefix"]], [genome["fasta"]], circular=False)

    ginis = primer_attributes.get_gini_from_txt_for_one_k(
        primers,
        genome["prefix"],
        genome["fasta"],
        seq_length=6_000,
        circular=False,
        position_cache=reused,
    )

    survivors = [
        primer
        for primer, (forward, reverse) in ginis.items()
        if not (math.isnan(forward) and math.isnan(reverse))
    ]
    assert survivors == primers, (
        "the second run dropped primers from the Gini gate; the reused position "
        f"map was short and the surviving pool was {survivors}"
    )


def test_overwrite_still_rescans_everything(genome, monkeypatch):
    """`overwrite=True` takes the fallback path and must not be short-circuited."""
    ss.get_positions([A12], [genome["prefix"]], [genome["fasta"]], circular=False)

    calls = _record_scans(monkeypatch)
    result = ss.get_positions(
        [A12], [genome["prefix"]], [genome["fasta"]], circular=False, overwrite=True
    )

    assert calls == [], "overwrite=True must not reach the Aho-Corasick branch"
    assert result is None, "the fallback path returns None by contract"
