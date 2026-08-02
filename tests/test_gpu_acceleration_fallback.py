"""Tests for the optional GPU-acceleration module.

This module is documented as experimental, and the audit plan put it out of
scope on the grounds that its current Python-loop form gives no real speedup.
That is a fair reason not to optimise it. It is not a reason to leave it
unexamined, because it contains a `PositionDatabase` that reads the same HDF5
position files everything else reads -- and reading those wrongly has been the
most productive bug class in this codebase.

Two things are worth asserting about an "acceleration" module regardless of
whether the acceleration is real:

1. It must agree with the canonical implementation. A faster path that returns
   different numbers is worse than no fast path, because the difference is
   invisible -- nothing fails, the scores are just quietly different depending
   on whether an optional package happens to be installed.

2. It must degrade cleanly when the optional dependency is absent, which is the
   normal case.

WHAT THESE TESTS DO NOT COVER: the CuPy kernels themselves. `cupy` is not
installed here, so every path below exercises the NumPy fallback. The GPU
branches remain unverified, and a test suite cannot pretend otherwise.
"""

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

from neoswga.core import reaction_conditions as rc
from neoswga.core.gpu_acceleration import (
    GPUThermodynamics,
    PositionDatabase,
    get_gpu_info,
    is_gpu_available,
)
from neoswga.core.thermodynamics import gc_content, reverse_complement

PRIMER = "TTGACCATGA"
PRIMERS = ["TTGACCATGA", "GCATTACGGT", "AACCGGTTAC", "ACGTACGTACGT"]
FORWARD_SITES = [1_000, 4_000]
REVERSE_SITES = [2_500, 7_000]


@pytest.fixture
def conditions():
    return rc.ReactionConditions(temp=30.0, polymerase="phi29")


@pytest.fixture
def position_file(tmp_path):
    """A position file in the layout the pipeline writes."""
    prefix = str(tmp_path / "g")
    with h5py.File(f"{prefix}_{len(PRIMER)}mer_positions.h5", "w") as f:
        f.create_dataset(PRIMER, data=np.array(FORWARD_SITES, dtype=np.int32))
        f.create_dataset(reverse_complement(PRIMER), data=np.array(REVERSE_SITES, dtype=np.int32))
    return prefix


# ----------------------------------------------------------------------
# PositionDatabase: the same strand bug, a third time
# ----------------------------------------------------------------------


def test_position_database_reads_both_strands(position_file):
    """It split the primer's array by sign and never looked at the revcomp.

    `strand_conventions` documents this module as using sign-encoded positions
    (negative means reverse), but the project writes no negative values --
    forward hits live under the primer, reverse hits under its reverse
    complement. So the reverse array was always empty and every reverse-strand
    site was invisible. The identical mistake was found in `amplicon_network`.
    """
    db = PositionDatabase([position_file])
    forward, reverse = db.get_positions(PRIMER, position_file)

    assert sorted(forward) == FORWARD_SITES
    assert sorted(reverse) == REVERSE_SITES


def test_batch_lookup_agrees_with_the_single_lookup(position_file):
    """Two code paths did the same split; both had the bug, so both are pinned.

    They must also agree with each other, or a caller gets different positions
    depending on whether it batched.
    """
    db = PositionDatabase([position_file])

    single = db.get_positions(PRIMER, position_file)
    db_fresh = PositionDatabase([position_file])
    batched = db_fresh.batch_get_positions([PRIMER], position_file)[PRIMER]

    assert sorted(single[0]) == sorted(batched[0])
    assert sorted(single[1]) == sorted(batched[1])


def test_sign_encoded_files_are_still_read(tmp_path):
    """The documented convention is honoured where a file actually uses it."""
    prefix = str(tmp_path / "signed")
    with h5py.File(f"{prefix}_{len(PRIMER)}mer_positions.h5", "w") as f:
        signed = FORWARD_SITES + [-p for p in REVERSE_SITES]
        f.create_dataset(PRIMER, data=np.array(signed, dtype=np.int32))

    db = PositionDatabase([prefix])
    forward, reverse = db.get_positions(PRIMER, prefix)

    assert sorted(forward) == FORWARD_SITES
    assert sorted(reverse) == REVERSE_SITES


def test_absent_primer_returns_empty_arrays(position_file):
    """Absent must be empty, not missing -- callers index the result."""
    db = PositionDatabase([position_file])
    forward, reverse = db.get_positions("CGCGCGCGCG", position_file)

    assert len(forward) == 0 and len(reverse) == 0


def test_repeated_lookups_are_served_from_the_cache(position_file):
    """The cache is this class's actual contribution over reading HDF5 directly.

    There is no hit/miss instrumentation -- `get_cache_stats` reports only size,
    capacity and utilisation -- so effectiveness is checked by identity: a
    cached lookup returns the very same tuple object rather than a fresh read.
    """
    db = PositionDatabase([position_file])

    first = db.get_positions(PRIMER, position_file)
    assert db.get_cache_stats()["size"] == 1

    second = db.get_positions(PRIMER, position_file)
    assert second is first, "second lookup re-read instead of using the cache"


def test_cached_result_matches_the_first_read(position_file):
    """A cache that returns something different from the source is worse than
    none -- and the strand fix has to survive the cache round trip."""
    db = PositionDatabase([position_file])

    first = db.get_positions(PRIMER, position_file)
    second = db.get_positions(PRIMER, position_file)

    assert sorted(first[0]) == sorted(second[0])
    assert sorted(first[1]) == sorted(second[1])


# ----------------------------------------------------------------------
# Agreement with the canonical implementation
# ----------------------------------------------------------------------


def test_batch_tm_matches_the_canonical_calculation(conditions):
    """A faster path returning different numbers is worse than no fast path.

    Nothing fails when they diverge; the scores just quietly differ depending on
    whether an optional package is installed.
    """
    calculator = GPUThermodynamics(conditions)
    batched = calculator.batch_calculate_tm(PRIMERS)

    for primer, value in zip(PRIMERS, batched):
        assert value == pytest.approx(conditions.calculate_effective_tm(primer), abs=1e-6)


def test_batch_gc_matches_the_canonical_calculation(conditions):
    calculator = GPUThermodynamics(conditions)
    batched = calculator.calculate_gc_content_batch(PRIMERS)

    for primer, value in zip(PRIMERS, batched):
        assert float(value) == pytest.approx(gc_content(primer))


def test_batch_results_are_ordered_like_the_input(conditions):
    """Results are matched back to primers positionally, so a reordering would
    attach every value to the wrong sequence."""
    calculator = GPUThermodynamics(conditions)

    forward = calculator.batch_calculate_tm(PRIMERS)
    reversed_order = calculator.batch_calculate_tm(list(reversed(PRIMERS)))

    assert list(forward) == pytest.approx(list(reversed(reversed_order)))


def test_batch_of_one_matches_a_scalar_call(conditions):
    calculator = GPUThermodynamics(conditions)
    assert calculator.batch_calculate_tm([PRIMER])[0] == pytest.approx(
        conditions.calculate_effective_tm(PRIMER), abs=1e-6
    )


def test_empty_batch_returns_empty(conditions):
    calculator = GPUThermodynamics(conditions)
    assert len(calculator.batch_calculate_tm([])) == 0


def test_batch_results_are_finite(conditions):
    """A NaN here propagates into scoring without failing."""
    import math

    calculator = GPUThermodynamics(conditions)
    for value in calculator.batch_calculate_tm(PRIMERS):
        assert math.isfinite(float(value))


# ----------------------------------------------------------------------
# Degrading without the optional dependency
# ----------------------------------------------------------------------


def test_gpu_availability_is_reported_honestly():
    """`is_gpu_available` gates every fast path; a wrong answer either skips
    acceleration silently or tries to use a device that is not there."""
    available = is_gpu_available()
    assert isinstance(available, bool)

    try:
        import cupy  # noqa: F401

        has_cupy = True
    except ImportError:
        has_cupy = False

    if not has_cupy:
        assert available is False


def test_gpu_info_is_returned_even_without_a_device():
    """Reported in diagnostics; it must not raise when there is no GPU."""
    info = get_gpu_info()
    assert isinstance(info, dict)


def test_calculations_work_with_no_gpu_present(conditions):
    """The normal case. Everything above already runs on the NumPy fallback --
    this states that plainly rather than leaving it implied.
    """
    assert is_gpu_available() is False

    calculator = GPUThermodynamics(conditions)
    assert len(calculator.batch_calculate_tm(PRIMERS)) == len(PRIMERS)


@pytest.mark.skipif(is_gpu_available(), reason="a GPU is present; fallback is not in use")
def test_the_gpu_kernels_themselves_are_not_covered_here():
    """Recorded so the coverage figure is not mistaken for verification.

    Every test in this file exercises the NumPy fallback, because cupy is not
    installed. The CuPy branches are unverified, and would need a machine with a
    device to test properly.
    """
    import inspect

    from neoswga.core import gpu_acceleration

    source = inspect.getsource(gpu_acceleration)
    assert "cupy" in source, "the module no longer has a GPU path to leave untested"
