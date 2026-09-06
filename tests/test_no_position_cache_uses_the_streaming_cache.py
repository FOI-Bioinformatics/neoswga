"""`--no-position-cache` has to select a different cache, and that cache has to work.

The flag advertises "Disable position cache (slower)". `optimize_step4` accepted
the value as `use_cache` and never referenced it -- its own docstring said
"always True with new system" -- while `run_step4` logged "Position cache: False"
first, which reads as confirmation the request was honoured. Nothing changed.

`StreamingPositionCache` is the memory-mapped alternative the flag should
select, and it had its own defect: it discovered HDF5 files with
`for k in range(6, 13)`, so a 12-mer design worked and any longer primer -- the
equiphi29 (12-18) and bst (15-25) ranges the tool ships presets for -- opened no
files at all and reported every position as absent. That is indistinguishable
downstream from a primer that genuinely does not bind.
"""

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

from neoswga.core.position_cache import StreamingPositionCache

SITES = [100, 5_000, 12_000]


def _write(prefix, k, primer):
    with h5py.File(f"{prefix}_{k}mer_positions.h5", "w") as handle:
        handle.create_dataset(primer, data=np.array(SITES, dtype=np.int64))


@pytest.mark.parametrize(
    "k,primer",
    [
        (10, "ACCACAGATA"),
        (12, "ACCACAGATAGC"),
        # equiphi29 runs at 12-18 and bst at 15-25; both were invisible.
        (16, "ACCACAGATAGCTTGA"),
        (20, "ACCACAGATAGCTTGACCTA"),
    ],
)
def test_the_streaming_cache_finds_primers_of_any_shipped_length(tmp_path, k, primer):
    prefix = str(tmp_path / "target")
    _write(prefix, k, primer)

    cache = StreamingPositionCache([prefix])
    positions = cache.get_positions(prefix, primer, "forward")

    assert list(positions) == SITES, (
        f"a {k}-mer's positions came back as {list(positions)}; the cache is "
        "not discovering the HDF5 file for that primer length"
    )


def test_the_optimizer_selects_the_streaming_cache_when_asked(monkeypatch, tmp_path):
    """The flag must change which class is built."""
    from neoswga.core import unified_optimizer as uo

    built = []

    class _Spy:
        def __init__(self, *args, **kwargs):
            built.append(type(self).__name__)

        def get_positions(self, *a, **k):
            return np.array([], dtype=np.int64)

    class _Plain(_Spy):
        pass

    class _Streaming(_Spy):
        pass

    monkeypatch.setattr(uo, "PositionCache", _Plain)
    monkeypatch.setattr(uo, "StreamingPositionCache", _Streaming, raising=False)

    assert uo._make_position_cache(["fg"], ["ACGT"], use_cache=True) is not None
    assert built == ["_Plain"]

    built.clear()
    uo._make_position_cache(["fg"], ["ACGT"], use_cache=False)
    assert built == ["_Streaming"], (
        "--no-position-cache did not select the streaming cache; the flag is "
        f"still inert. Built: {built}"
    )
