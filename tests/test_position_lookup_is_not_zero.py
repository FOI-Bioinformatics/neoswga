"""A position lookup that cannot be answered is not "binds nowhere".

`BaseOptimizer.get_primer_positions` caught every exception and returned an
empty array. `PositionCache` raises `MissingPositionsError` for a prefix that
was never indexed (Known Issue 13), and this wrapper turned that back into zero
sites for every optimizer: an unindexed background read as perfectly specific
(2026-11-23 sweep, tier 1 item 5).
"""

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

from neoswga.core.exceptions import ReferenceDataError
from neoswga.core.hybrid_optimizer import HybridBaseOptimizer
from neoswga.core.position_cache import MissingPositionsError, PositionCache

PRIMER = "ACGTACGTAC"


@pytest.fixture
def optimizer(tmp_path):
    prefix = str(tmp_path / "target")
    with h5py.File(f"{prefix}_10mer_positions.h5", "w") as f:
        f.create_dataset(PRIMER, data=np.array([10, 500], dtype=np.int64))
    cache = PositionCache([prefix], [PRIMER])
    return HybridBaseOptimizer(
        position_cache=cache, fg_prefixes=[prefix], fg_seq_lengths=[5_000]
    ), prefix


def test_missing_index_raises_reference_data_error(optimizer, monkeypatch):
    opt, prefix = optimizer

    def _missing(*_args, **_kwargs):
        raise MissingPositionsError("never indexed", missing=[PRIMER])

    monkeypatch.setattr(opt.cache, "get_positions", _missing)
    with pytest.raises(ReferenceDataError) as info:
        opt.get_primer_positions(PRIMER, "host", "both")
    assert info.value.artifact == "host"


def test_an_unreadable_index_raises(optimizer, monkeypatch):
    opt, _prefix = optimizer

    def _broken(*_args, **_kwargs):
        raise OSError("truncated HDF5")

    monkeypatch.setattr(opt.cache, "get_positions", _broken)
    with pytest.raises(ReferenceDataError):
        opt.get_primer_positions(PRIMER, "host", "both")


def test_absent_cache_raises(optimizer):
    opt, prefix = optimizer
    opt.cache = None
    with pytest.raises(ReferenceDataError):
        opt.get_primer_positions(PRIMER, prefix, "both")


def test_an_indexed_primer_still_answers(optimizer):
    opt, prefix = optimizer
    positions = opt.get_primer_positions(PRIMER, prefix, "both")
    assert sorted(positions.tolist()) == [10, 500]
