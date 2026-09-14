"""The two position caches must answer an unindexed prefix the same way.

Known Issue 13 records five call sites that read a background genome as empty
because `PositionCache.get_positions` answered a prefix it was never built over
with an empty array. The fix made that raise. It went one class deep:
`StreamingPositionCache` is the sibling `unified_optimizer` selects on a flag,
and it still returned `np.array([])` for a prefix it holds no file for.

An empty array is indistinguishable from a primer that genuinely binds nowhere,
which is what made the original defect invisible: a panel scored against a
background nobody had indexed looked perfectly specific.

The AST guard in `tests/test_expansion_counts_background.py` cannot catch this.
It matches construction of the literal name `PositionCache`, so every call site
that builds the streaming class is outside it.
"""

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

from neoswga.core.position_cache import (
    POSITION_DTYPE,
    MissingPositionsError,
    PositionCache,
    StreamingPositionCache,
)

PRIMER = "ACGTACGTACGT"


@pytest.fixture
def indexed(tmp_path):
    """One indexed prefix, and the name of one that was never built."""
    prefix = str(tmp_path / "present")
    with h5py.File(f"{prefix}_12mer_positions.h5", "w") as db:
        db.create_dataset(PRIMER, data=np.array([10, 200], dtype=np.int64))
    return prefix, str(tmp_path / "never_indexed")


def test_the_main_cache_refuses_an_unindexed_prefix(indexed):
    """The behaviour Known Issue 13 established, as the reference."""
    present, absent = indexed
    cache = PositionCache([present], [PRIMER])

    with pytest.raises(MissingPositionsError):
        cache.get_positions(absent, PRIMER, "both")


def test_the_streaming_cache_refuses_it_too(indexed):
    present, absent = indexed
    cache = StreamingPositionCache([present], [PRIMER])

    with pytest.raises(MissingPositionsError):
        cache.get_positions(absent, PRIMER, "both")


@pytest.mark.parametrize("strand", ["both", "forward", "reverse"])
def test_both_caches_agree_on_an_indexed_prefix(indexed, strand):
    """Parity is the point; refusing everything would also pass the tests above."""
    present, _ = indexed
    main = PositionCache([present], [PRIMER]).get_positions(present, PRIMER, strand)
    streaming = StreamingPositionCache([present], [PRIMER]).get_positions(present, PRIMER, strand)

    assert list(main) == list(streaming)


def test_a_primer_that_binds_nowhere_is_still_an_empty_answer(indexed):
    """An indexed prefix with no hits is a measurement, and must not raise."""
    present, _ = indexed
    other = "TTTTTTTTTTTT"

    assert len(StreamingPositionCache([present], [other]).get_positions(present, other)) == 0


def test_the_empty_answer_carries_the_genome_coordinate_dtype(indexed):
    """int64, so a caller concatenating results cannot be silently downcast.

    Known Issue 7 records genome coordinates saturating at the int32 ceiling;
    an empty float64 array from this path is a smaller version of the same
    class of problem.
    """
    present, _ = indexed
    other = "TTTTTTTTTTTT"
    cache = StreamingPositionCache([present], [other])

    for strand in ("both", "forward", "reverse"):
        assert cache.get_positions(present, other, strand).dtype == POSITION_DTYPE
