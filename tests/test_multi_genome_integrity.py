"""Regression tests that ensure multi-genome inputs are fully honored.

The audit surfaced two [0]-indexed sites that silently ignored additional
target or background genomes:

- the position-file cache check, which is now taken per foreground prefix in
  string_search rather than once for all of them in pipeline.py
- background_aware_optimizer.compare_optimizers (now aggregates all bg_prefixes)

These tests lock in the fix.
"""


def test_position_file_cache_is_decided_per_fg_prefix(tmp_path):
    """A prefix without a cached HDF5 is scanned even when another prefix has one.

    `pipeline.step2` used to answer this with one `all(...)` over the prefixes at
    `parameter.min_k`, which was read by nothing but a log line. The decision now
    lives in `string_search._reusable_positions` and is taken per prefix and per
    k, so this asserts it there.
    """
    import h5py

    from neoswga.core import string_search

    k = 8
    primers = ["ACGTACGT", "TTTTGGGG"]
    cached = str(tmp_path / "fg_a")
    uncached = str(tmp_path / "fg_b")
    genome = tmp_path / "g.fasta"
    genome.write_text(">g\n" + "ACGT" * 100 + "\n")

    with h5py.File(string_search.position_file_path(cached, k), "w") as handle:
        for primer in primers:
            handle.create_dataset(primer, data=[10])
    string_search.write_position_provenance(cached, str(genome), k, False)

    reused, to_scan, replace = string_search._reusable_positions(
        primers, cached, str(genome), k, False
    )
    assert to_scan == [] and sorted(reused) == sorted(primers) and replace is False

    reused, to_scan, replace = string_search._reusable_positions(
        primers, uncached, str(genome), k, False
    )
    assert (
        to_scan == primers and reused == {}
    ), "a prefix with no cached position file must still be scanned in full"


# `test_background_aware_optimizer_aggregates_all_bg_prefixes` was removed on
# 2026-09-10 with `background_aware_optimizer.compare_optimizers`, which nothing
# called. It read the function's SOURCE TEXT with `inspect.getsource` and
# asserted on the strings in it -- which is a tell in itself: a function that
# can only be checked by reading it is a function nothing runs.
