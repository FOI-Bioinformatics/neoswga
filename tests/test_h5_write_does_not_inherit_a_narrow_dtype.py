"""Rewriting an index must not truncate genome coordinates, or grow the file.

Known Issue 7 records genome coordinates being stored as int32, where every
offset past 2,147,483,647 saturated at the ceiling. Human and mouse are past it.

The per-dataset writer overwrote a dataset IN PLACE when the new site count
equalled the old, to avoid HDF5 fragmentation, and in place meant into the
existing dtype: an int32 index rescanned to the same site count could not
hold what it scanned. The sorted-blocks writer (`core/position_index.py`)
rewrites the whole index into a new file, so it inherits no dtype from the old
one and leaves no freed space behind. These tests pin both properties that the
in-place branch existed to protect, now that the branch is gone.
"""

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

from neoswga.core.position_index import open_index
from neoswga.core.string_search import write_to_h5py

PRIMER = "ACGTACGTACGT"
# Past the int32 ceiling of 2,147,483,647; reachable on human and mouse.
FAR = [2_147_483_600, 3_000_000_000, 3_100_000_000]


def test_an_equal_length_rewrite_keeps_far_coordinates(tmp_path):
    prefix = str(tmp_path / "idx")
    path = f"{prefix}_12mer_positions.h5"

    # An index written before coordinates were int64, with the same site count.
    with h5py.File(path, "w") as db:
        db.create_dataset(PRIMER, data=np.array([1, 2, 3], dtype=np.int32))

    write_to_h5py({PRIMER: FAR}, prefix)

    with open_index(path) as db:
        assert db[PRIMER].tolist() == FAR, "coordinates were truncated by an inherited dtype"


def test_an_int32_entry_that_is_not_rewritten_is_widened_on_conversion(tmp_path):
    prefix = str(tmp_path / "idx")
    path = f"{prefix}_12mer_positions.h5"
    with h5py.File(path, "w") as db:
        db.create_dataset(PRIMER, data=np.array([1, 2, 3], dtype=np.int32))

    write_to_h5py({"CCCCCCCCCCCC": FAR}, prefix)

    with open_index(path) as db:
        assert db[PRIMER].tolist() == [1, 2, 3]
        assert db[PRIMER].dtype == np.int64
        assert db["CCCCCCCCCCCC"].tolist() == FAR


def test_repeated_rewrites_do_not_grow_the_file(tmp_path):
    """What the in-place branch was for: HDF5 does not reclaim deleted space."""
    import os

    prefix = str(tmp_path / "idx")
    path = f"{prefix}_12mer_positions.h5"
    entries = {f"{i:012b}".replace("0", "A").replace("1", "C"): [i, i + 1] for i in range(200)}
    write_to_h5py(entries, prefix)
    first = os.path.getsize(path)
    for round_ in range(10):
        write_to_h5py({key: [round_, round_ + 1] for key in entries}, prefix)
    assert os.path.getsize(path) == first
