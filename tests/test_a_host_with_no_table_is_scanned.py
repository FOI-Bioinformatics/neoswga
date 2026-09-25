"""A background with no k-mer table is measured, not refused and not zeroed.

Counting a host builds a table bounded by the reference rather than by the
k-mer space above about k=15, so hg38 at k=18 is tens of gigabytes of text.
The counts a design actually needs are those of its own candidate list, which
one pass over the reference answers with nothing written to disk.

Two things are asserted here and neither is visible from a unit test of
`query_scan` itself. That the path is REACHED from `filter`, which is the
defect class this repository calls Known Issue 8: both ends existing while
nothing walks between them. And that the counts it produces equal the counts
the table path produces, because a background measured one way and gated
another is worse than either.
"""

import os

import pytest

from neoswga.core import filter as filter_module
from neoswga.core import kmer_tables, query_scan
from neoswga.core.kmer_backend import JellyfishBackend, KmcBackend


def _a_counter_is_available():
    return JellyfishBackend().available() or KmcBackend().available()


PRIMERS = ["ACGTTGCAAC", "GGGGCCCCAA", "TTTTTTTTTT"]


@pytest.fixture
def reference(tmp_path):
    """One host FASTA, two records, holding two of the three primers."""
    import random

    rng = random.Random(3)
    filler = "".join(rng.choice("ACGT") for _ in range(4000))
    records = [
        filler[:2000] + "ACGTTGCAAC" + filler[2000:3000],
        "GGGGCCCCAA" + filler[3000:] + "GGGGCCCCAA",
    ]
    path = tmp_path / "host.fasta"
    path.write_text("".join(f">rec{i}\n{seq}\n" for i, seq in enumerate(records)))
    return str(path)


def test_the_scan_is_reached_from_the_rate_lookup(tmp_path, reference):
    """`get_rates_for_one_species` answers for a prefix that was never counted."""
    prefix = str(tmp_path / "host")
    assert not kmer_tables.table_exists(prefix, 10), "the fixture must have no table"

    counts = filter_module.get_rates_for_one_species(PRIMERS, [prefix], [reference])

    assert set(counts) == set(PRIMERS), "every requested primer must be answered for"
    assert counts["ACGTTGCAAC"] == 1
    assert counts["GGGGCCCCAA"] == 2
    assert counts["TTTTTTTTTT"] == 0


def test_without_a_genome_an_uncounted_prefix_still_refuses(tmp_path):
    """The silent zero this fallback must not introduce.

    Reporting zeros for an uncounted background would say every candidate is
    absent from the host, which is the most permissive answer available.
    """
    prefix = str(tmp_path / "never_counted")
    with pytest.raises(FileNotFoundError, match="No 10-mer table"):
        filter_module.get_rates_for_one_species(PRIMERS, [prefix])


def test_mismatched_prefix_and_genome_lists_refuse_rather_than_pair_by_position(
    tmp_path, reference, caplog
):
    """Two prefixes and one genome cannot be paired, so neither is scanned.

    Pairing what happens to line up would score one reference's counts under
    another's name, which this repository has already done once.
    """
    prefixes = [str(tmp_path / "a"), str(tmp_path / "b")]
    with pytest.raises(FileNotFoundError):
        filter_module.get_rates_for_one_species(PRIMERS, prefixes, [reference])


@pytest.mark.skipif(not _a_counter_is_available(), reason="no k-mer counter installed")
def test_the_scanned_counts_equal_the_counted_ones(tmp_path, reference):
    """Same host, same primers: the two routes must not disagree."""
    from neoswga.core import kmer_counter

    prefix = str(tmp_path / "host")
    scanned = filter_module.get_rates_for_one_species(PRIMERS, [prefix], [reference])

    kmer_counter.run_jellyfish(reference, prefix, min_k=10, max_k=10, cpus=1)
    assert kmer_tables.table_exists(prefix, 10)
    counted = filter_module.get_rates_for_one_species(PRIMERS, [prefix], [reference])

    assert scanned == counted


@pytest.mark.skipif(not _a_counter_is_available(), reason="no k-mer counter installed")
def test_an_existing_table_is_preferred_over_scanning(tmp_path, reference, monkeypatch):
    """A counted prefix must not pay for a genome pass it does not need."""
    from neoswga.core import kmer_counter

    prefix = str(tmp_path / "host")
    kmer_counter.run_jellyfish(reference, prefix, min_k=10, max_k=10, cpus=1)

    def fail(*args, **kwargs):
        raise AssertionError("the table exists, so nothing should be scanned")

    monkeypatch.setattr(query_scan, "count_kmers", fail)
    counts = filter_module.get_rates_for_one_species(PRIMERS, [prefix], [reference])
    assert counts["GGGGCCCCAA"] == 2


def test_a_named_genome_that_is_absent_is_an_error_not_a_zero(tmp_path):
    prefix = str(tmp_path / "host")
    with pytest.raises((FileNotFoundError, OSError)):
        filter_module.get_rates_for_one_species(
            PRIMERS, [prefix], [str(tmp_path / "not_here.fasta")]
        )


def test_the_cost_line_names_the_reference_and_the_query_count(reference):
    message = query_scan.describe_cost(reference, 18, 2000)
    assert os.path.basename(reference) in message
    assert "2,000" in message
