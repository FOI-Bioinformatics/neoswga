"""Removing a primer from a coverage set does not rebuild the whole set.

Finding A3. _prune_background evaluated every candidate removal by rebuilding a
BipartiteGraph over the whole set, so pruning 160 primers to 24 was on the order
of a million primer-binning operations to make 136 decisions. The counter form is
arithmetic-identical: the coverage lost by removing a primer is the number of
bins it is the only coverer of.
"""

import random

import pytest

from neoswga.core.coverage_counter import CoverageCounter

h5py = pytest.importorskip("h5py")

import numpy as np

from neoswga.core.hybrid_optimizer import HybridOptimizer
from neoswga.core.position_cache import PositionCache
from neoswga.core.thermodynamics import reverse_complement


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


def _build_optimizer_and_primers(tmp_path):
    """A real `HybridOptimizer` over synthetic HDF5 position data.

    Mirrors the genome/cache pattern in tests/test_hybrid_optimizer_run.py so
    the tests below need no jellyfish. The previous version of this test read
    `examples/plasmid_example`, which `tests/conftest.py` only primes when
    jellyfish is on PATH -- so on a runner without it, this was the only test
    pinning the `CoverageCounter` denominator against `_calculate_coverage`,
    and it disappeared silently. h5py is a core dependency (unlike jellyfish,
    an external binary), so `pytest.importorskip` on it above is not expected
    to skip in practice; the fixture itself needs no external tool at all.
    """
    genome_length = 200_000
    rng = random.Random(4471)
    seq = list("".join(rng.choice("ACGT") for _ in range(genome_length)))

    primers = []
    for i in range(14):
        primer = "".join(rng.choice("ACGT") for _ in range(10))
        if primer in primers:
            continue
        primers.append(primer)
        # Each primer gets 2-4 sites, clustered differently per primer.
        for j in range(2 + (i % 3)):
            pos = (i * 14_000 + j * 700) % (genome_length - 20)
            seq[pos : pos + 10] = list(primer)
    seq = "".join(seq)

    prefix = str(tmp_path / "target")
    with h5py.File(f"{prefix}_10mer_positions.h5", "w") as f:
        for primer in primers:
            for key in {primer, reverse_complement(primer)}:
                positions, i = [], seq.find(key)
                while i != -1:
                    positions.append(i)
                    i = seq.find(key, i + 1)
                if positions:
                    f.create_dataset(key, data=np.array(positions, dtype=np.int32))

    cache = PositionCache([prefix], primers)
    optimizer = HybridOptimizer(
        position_cache=cache,
        fg_prefixes=[prefix],
        fg_seq_lengths=[genome_length],
        bin_size=1_000,
        coverage_reach=3_000,
    )
    return optimizer, primers


def test_counter_fraction_matches_calculate_coverage_through_several_removals(tmp_path):
    """The counter must track `_calculate_coverage` through a sequence of
    removals, not just at the entry state.

    Finding 2. The committed version of this test compared a single number at
    the entry-state counter and never called `remove()`, so the incremental
    property this whole module exists to make cheap -- exactness after a
    removal, not just at the start -- was pinned only abstractly, by the
    hand-built dict in `test_counter_matches_a_recount_after_each_removal`.
    This walks every primer out one at a time and checks agreement with a real
    `_calculate_coverage` rebuild of the true remaining set at every step,
    ending at the empty set.

    Finding 3. Rebuilt on synthetic HDF5 data (see `_build_optimizer_and_primers`)
    rather than the plasmid example, which needs jellyfish to have been primed.
    """
    optimizer, primers = _build_optimizer_and_primers(tmp_path)

    counter = optimizer._build_coverage_counter(primers)
    remaining = list(primers)
    assert counter.covered_fraction() == pytest.approx(
        optimizer._calculate_coverage(remaining), abs=1e-9
    )

    for primer in list(primers):
        counter.remove(primer)
        remaining.remove(primer)
        assert counter.covered_fraction() == pytest.approx(
            optimizer._calculate_coverage(remaining), abs=1e-9
        )


def test_duplicate_primers_do_not_corrupt_the_counter(tmp_path):
    """A duplicated primer in the input must not desynchronise the counter
    from the true remaining set.

    Finding 4. `_prune_background` used to rebuild `_calculate_coverage` from
    `current_primers` on every step, so a duplicate entry was harmless -- the
    rebuild always saw the true remaining set. The counter is keyed by primer
    sequence: `current_primers.remove(p)` drops one occurrence of a duplicate
    while `counter.remove(p)` drops that primer's bins entirely, so the two
    disagreed once a duplicate reached the loop. `_prune_background` now
    dedupes its input at entry so the list and the counter start, and stay, in
    agreement.

    `target_size=6` is load-bearing, not arbitrary. At `target_size=4` on this
    7-primer, one-duplicate input, both copies of the duplicate are gone
    before the loop ends, so the list and the counter never get a chance to
    disagree -- the dedupe fix could be deleted and this would still pass.
    `target_size=6` (one removal from a 7-entry, 6-distinct-primer list) is
    the smallest size that discriminates: with the dedupe removed, one
    duplicate copy of the surviving primer is still present after that single
    removal, so `current_primers.remove` and `counter.remove` disagree and the
    reported coverage understates the true rebuild by about 7.6%. Verified by
    temporarily removing the dedupe line and confirming this assertion fails.
    """
    optimizer, primers = _build_optimizer_and_primers(tmp_path)
    duplicated = primers[:6] + [primers[0]]

    kept, coverage, _bg = optimizer._prune_background(duplicated, target_size=6, verbose=False)

    assert len(kept) == len(set(kept)), "a duplicate survived pruning"
    assert coverage == pytest.approx(optimizer._calculate_coverage(kept), abs=1e-9)
