"""Rewriting a dataset in place must not truncate genome coordinates.

Known Issue 7 records genome coordinates being stored as int32, where every
offset past 2,147,483,647 saturated at the ceiling. Human and mouse are past it.
The stored dtype is int64 now, but one write path can still inherit a narrow one
from a file written earlier.

`write_to_h5py` overwrites a dataset in place when the new site count equals the
old, to avoid HDF5 fragmentation. In place means into the existing dataset, with
the existing dtype. So an older int32 index plus a rescan that happens to find
the same number of sites cannot store what it scanned. Observed on h5py 3.x it
raises OverflowError, so the run dies rather than lying; on a build that casts
instead, the coordinates would be truncated into the file and the truncation
would be indistinguishable from a measurement. Neither is acceptable, and the
same fix covers both: a dataset that cannot hold the values is recreated.

The round-trip test for far coordinates does not reach this branch: it creates
the dataset, and the other write test changes the site count, so both take the
create path instead.
"""

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

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

    with h5py.File(path, "r") as db:
        assert (
            list(db[PRIMER]) == FAR
        ), "coordinates were truncated by the dataset's inherited dtype"


def test_the_branch_under_test_is_the_in_place_one(tmp_path):
    """Guard the guard: a different site count would take the create path."""
    prefix = str(tmp_path / "idx")
    path = f"{prefix}_12mer_positions.h5"
    with h5py.File(path, "w") as db:
        db.create_dataset(PRIMER, data=np.array([1, 2, 3], dtype=np.int32))

    with h5py.File(path, "r") as db:
        assert len(db[PRIMER]) == len(FAR)


def test_an_equal_length_rewrite_still_avoids_fragmentation(tmp_path):
    """The optimisation is kept where the dtype already fits."""
    prefix = str(tmp_path / "idx")
    path = f"{prefix}_12mer_positions.h5"
    with h5py.File(path, "w") as db:
        db.create_dataset(PRIMER, data=np.array([1, 2, 3], dtype=np.int64))
        original = db[PRIMER].id.get_offset()

    write_to_h5py({PRIMER: [10, 20, 30]}, prefix)

    with h5py.File(path, "r") as db:
        assert list(db[PRIMER]) == [10, 20, 30]
        assert db[PRIMER].id.get_offset() == original, "dataset was needlessly recreated"
