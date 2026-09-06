"""Removing a primer from a coverage set does not rebuild the whole set.

Finding A3. _prune_background evaluated every candidate removal by rebuilding a
BipartiteGraph over the whole set, so pruning 160 primers to 24 was on the order
of a million primer-binning operations to make 136 decisions. The counter form is
arithmetic-identical: the coverage lost by removing a primer is the number of
bins it is the only coverer of.
"""

import json

import pytest

from neoswga.core.coverage_counter import CoverageCounter


def test_counter_matches_a_recount_after_each_removal():
    coverage = {
        "P1": {1, 2, 3},
        "P2": {3, 4},
        "P3": {4, 5},
        "P4": {5, 6, 7},
    }
    counter = CoverageCounter(total_bins=10)
    for primer, bins in coverage.items():
        counter.add(primer, bins)

    assert counter.covered_fraction() == pytest.approx(7 / 10)

    for primer in ["P2", "P1", "P3"]:
        expected_after = len(set().union(*(b for p, b in coverage.items() if p != primer)))
        assert counter.loss_if_removed(primer) == counter.covered_count() - expected_after
        counter.remove(primer)
        coverage.pop(primer)
        assert counter.covered_count() == expected_after


def test_loss_is_zero_for_a_primer_whose_bins_are_all_shared():
    counter = CoverageCounter(total_bins=4)
    counter.add("P1", {1, 2})
    counter.add("P2", {1, 2})
    assert counter.loss_if_removed("P1") == 0
    assert counter.loss_if_removed("P2") == 0


def test_removing_an_unknown_primer_is_a_no_op():
    counter = CoverageCounter(total_bins=4)
    counter.add("P1", {1})
    counter.remove("absent")
    assert counter.covered_count() == 1


def test_empty_counter_reports_zero_not_one():
    """An empty region set is not full coverage; see finding A7."""
    counter = CoverageCounter(total_bins=10)
    assert counter.covered_fraction() == 0.0


def test_counter_fraction_matches_calculate_coverage_on_a_real_pool():
    """The counter must compute the same number the rebuild computed.

    Uses the plasmid example, which conftest primes, so the position files
    exist.
    """
    pytest.importorskip("h5py")
    import os

    example = os.path.join(os.path.dirname(__file__), "..", "examples", "plasmid_example")
    if not os.path.exists(os.path.join(example, "step3_df.csv")):
        pytest.skip("plasmid example not primed")

    import pandas as pd

    from neoswga.core.hybrid_optimizer import HybridOptimizer

    params = json.load(open(os.path.join(example, "params.json")))
    primers = pd.read_csv(os.path.join(example, "step3_df.csv"))["primer"].tolist()[:12]

    from neoswga.core.position_cache import PositionCache

    # fg_prefixes in the plasmid params is the relative "pcDNA"; the HDF5 files
    # sit beside params.json, so resolve against the example directory.
    prefixes = [os.path.join(example, p) for p in params["fg_prefixes"]]

    cache = PositionCache(prefixes, primers)
    optimizer = HybridOptimizer(
        cache,
        fg_prefixes=prefixes,
        fg_seq_lengths=params["fg_seq_lengths"],
        polymerase=params.get("polymerase", "phi29"),
    )
    counter = optimizer._build_coverage_counter(primers)
    assert counter.covered_fraction() == pytest.approx(
        optimizer._calculate_coverage(primers), abs=1e-9
    )
