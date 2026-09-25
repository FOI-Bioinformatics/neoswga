"""One place every consumer reads or writes a binding-position index.

A position index maps each k-mer AS WRITTEN to the forward-strand offsets at
which it occurs in one reference. Reverse-strand sites are read from the
reverse complement's entry, so keys are never canonicalised: storing only the
canonical form would lose one strand of every non-palindromic primer.

Two layouts are read, and the newer one is written.

  per_dataset    one HDF5 dataset per k-mer. What every index written before
                 2026-09-25 holds. HDF5 spends several hundred bytes of object
                 header and B-tree node per dataset, so on the Wolbachia
                 design 376 MB of a 381 MB file was metadata around 5 MB of
                 positions, and loading it took 27.6 s.
  sorted_blocks  three datasets: the k-mers sorted as bytes, an offsets array
                 of length n + 1, and every position concatenated in key
                 order. Entry i is positions[offsets[i]:offsets[i + 1]]. The
                 same index is 24.7 MB and loads in 0.36 s, entry for entry
                 identical
                 (docs/validation/position_index_layout_2026-09-25.md).

Absence and emptiness stay distinct in both. A k-mer whose entry is empty was
scanned and occurs nowhere, which is a measurement. A k-mer with no entry was
never scanned, and reporting zero sites for it would read as perfect
specificity, which is the silent-zero shape of Known Issues 5, 6, 13 and 15.
`get` returns an empty array for the first and None for the second.

`index_format_version` is NOT the layout. It records whether the SITES are
sound (version 2 is the first whose scan respects record joins), and a layout
change alters no site, so it is left alone and `index_layout` names the layout.
"""

from __future__ import annotations

import contextlib
import os
from collections.abc import Iterable, Iterator, Mapping, Sequence

import h5py
import numpy as np

# Datasets whose names begin with '#' are never k-mers, so they cannot collide
# with an entry in the per-dataset layout.
RECORD_STARTS_KEY = "#record_starts"
KEYS_KEY = "#keys"
OFFSETS_KEY = "#offsets"
POSITIONS_KEY = "#positions"

LAYOUT_ATTR = "index_layout"
SORTED_BLOCKS = "sorted_blocks"
PER_DATASET = "per_dataset"

POSITION_DTYPE = np.int64

# Reads of the concatenated positions are coalesced: two requested entries
# separated by at most this many unrequested positions are read in one slice.
# Each HDF5 read costs tens of microseconds regardless of size, so a request
# for half a million entries read one at a time would spend most of its time
# in call overhead. The waste this admits is bounded by the gap per entry.
_COALESCE_GAP = 4096
# And no single read materialises more than this many positions (128 MB).
_MAX_RUN = 16 * 1024 * 1024


class MixedLayoutError(OSError):
    """An index holding sorted blocks AND per-k-mer datasets.

    Only a release predating the sorted-blocks layout produces one, by
    writing its datasets into a file this release converted. The reader would
    see the blocks and never the datasets, so the file is refused. It is an
    OSError so that every caller treating an unreadable index as one to
    rescan does so here too.
    """


def index_path(prefix: str, k: int) -> str:
    """Where the binding positions for one prefix and one k live."""
    return f"{prefix}_{k}mer_positions.h5"


def _encode(key: str) -> bytes:
    return key.encode("ascii")


class PositionIndex:
    """Read access to one position index file, in either layout.

    Open it as a context manager. The file handle is held for the object's
    lifetime; nothing is cached beyond the sorted keys and offsets of the
    sorted-blocks layout, which are what a lookup needs.
    """

    def __init__(self, path: str, handle: h5py.File | None = None):
        self.path = path
        self._owns_handle = handle is None
        self._handle = h5py.File(path, "r") if handle is None else handle
        self.layout = SORTED_BLOCKS if KEYS_KEY in self._handle else PER_DATASET
        if self.layout == SORTED_BLOCKS:
            # Root names iterate in byte order and '#' sorts before every
            # base, so this reads past at most the four bookkeeping names.
            stray = next((name for name in self._handle if not name.startswith("#")), None)
            if stray is not None:
                if self._owns_handle:
                    self._handle.close()
                raise MixedLayoutError(
                    f"{path} holds sorted blocks and also a per-k-mer dataset "
                    f"({stray!r}), which only an older NeoSWGA release writes. "
                    f"Its entries cannot be trusted as a whole. Re-run 'neoswga "
                    f"filter -j params.json' with this release to rebuild it."
                )
        self._keys: np.ndarray | None = None
        self._offsets: np.ndarray | None = None
        if self.layout == SORTED_BLOCKS:
            self._keys = self._handle[KEYS_KEY][()]
            self._offsets = np.asarray(self._handle[OFFSETS_KEY][()], dtype=np.int64)

    # -- lifecycle ------------------------------------------------------------

    def close(self) -> None:
        if self._owns_handle:
            self._handle.close()

    def __enter__(self) -> PositionIndex:
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    # -- metadata -------------------------------------------------------------

    @property
    def attrs(self) -> dict:
        return dict(self._handle.attrs)

    def record_starts(self) -> list[int] | None:
        """Where each FASTA record begins, or None for an index without them."""
        if RECORD_STARTS_KEY not in self._handle:
            return None
        return [int(v) for v in np.asarray(self._handle[RECORD_STARTS_KEY][()])]

    # -- entries --------------------------------------------------------------

    def _slot(self, key: str) -> int | None:
        assert self._keys is not None
        try:
            encoded = _encode(key)
        except UnicodeEncodeError:
            return None
        i = int(np.searchsorted(self._keys, encoded))
        if i < len(self._keys) and self._keys[i] == encoded:
            return i
        return None

    def __contains__(self, key: object) -> bool:
        if not isinstance(key, str) or key.startswith("#"):
            return False
        if self.layout == PER_DATASET:
            return key in self._handle
        return self._slot(key) is not None

    def __len__(self) -> int:
        if self.layout == SORTED_BLOCKS:
            assert self._keys is not None
            return len(self._keys)
        return sum(1 for name in self._handle if not name.startswith("#"))

    def keys(self) -> list[str]:
        if self.layout == SORTED_BLOCKS:
            assert self._keys is not None
            return [bytes(k).decode("ascii") for k in self._keys]
        return [name for name in self._handle if not name.startswith("#")]

    def get(self, key: str, default=None):
        """The positions of one k-mer.

        Empty when it was scanned and occurs nowhere; `default` (None unless
        given) when it was never scanned.
        """
        if key not in self:
            return default
        if self.layout == PER_DATASET:
            return np.asarray(self._handle[key][()], dtype=POSITION_DTYPE)
        i = self._slot(key)
        assert i is not None and self._offsets is not None
        start, stop = int(self._offsets[i]), int(self._offsets[i + 1])
        if start == stop:
            return np.empty(0, dtype=POSITION_DTYPE)
        return np.asarray(self._handle[POSITIONS_KEY][start:stop], dtype=POSITION_DTYPE)

    def __getitem__(self, key: str) -> np.ndarray:
        """Mapping access: the entry, or KeyError when it was never scanned."""
        found = self.get(key)
        if found is None:
            raise KeyError(key)
        return found

    def get_many(self, keys: Iterable[str]) -> dict[str, np.ndarray]:
        """Positions for every requested k-mer that has an entry.

        A k-mer with no entry is left out of the result rather than given an
        empty array, so a caller can still tell never-scanned from absent.
        """
        wanted = list(dict.fromkeys(keys))
        if self.layout == PER_DATASET:
            present = {}
            for key in wanted:
                value = self.get(key)
                if value is not None:
                    present[key] = value
            return present

        assert self._keys is not None and self._offsets is not None
        slots = []
        for key in wanted:
            i = self._slot(key)
            if i is not None:
                slots.append((key, i))
        if not slots:
            return {}

        starts = np.array([self._offsets[i] for _, i in slots], dtype=np.int64)
        stops = np.array([self._offsets[i + 1] for _, i in slots], dtype=np.int64)
        order = np.argsort(starts, kind="stable")

        found: dict[str, np.ndarray] = {}
        dataset = self._handle[POSITIONS_KEY]
        run: list[int] = []
        run_start = run_stop = 0

        def flush() -> None:
            if not run:
                return
            block = np.asarray(dataset[run_start:run_stop], dtype=POSITION_DTYPE)
            for j in run:
                key = slots[j][0]
                found[key] = block[starts[j] - run_start : stops[j] - run_start].copy()

        for j in order.tolist():
            s, e = int(starts[j]), int(stops[j])
            if run and s - run_stop <= _COALESCE_GAP and e - run_start <= _MAX_RUN:
                run.append(j)
                run_stop = max(run_stop, e)
                continue
            flush()
            run, run_start, run_stop = [j], s, e
        flush()
        return found

    def items(self) -> Iterator[tuple[str, np.ndarray]]:
        """Every entry, in key order for sorted blocks."""
        if self.layout == PER_DATASET:
            for name in self.keys():
                yield name, np.asarray(self._handle[name][()], dtype=POSITION_DTYPE)
            return
        assert self._keys is not None and self._offsets is not None
        positions = np.asarray(self._handle[POSITIONS_KEY][()], dtype=POSITION_DTYPE)
        for i, raw in enumerate(self._keys):
            yield (
                bytes(raw).decode("ascii"),
                positions[self._offsets[i] : self._offsets[i + 1]].copy(),
            )

    def _as_arrays(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """(sorted keys as bytes, offsets, positions) for the whole index."""
        if self.layout == SORTED_BLOCKS:
            assert self._keys is not None and self._offsets is not None
            positions = np.asarray(self._handle[POSITIONS_KEY][()], dtype=POSITION_DTYPE)
            return self._keys, self._offsets, positions
        return _pack(dict(self.items()))


def open_index(path: str) -> PositionIndex:
    return PositionIndex(path)


def _pack(entries: Mapping[str, Sequence[int] | np.ndarray]):
    """Sorted keys, offsets and concatenated positions for a mapping."""
    names = sorted(entries)
    arrays = [np.asarray(entries[name], dtype=POSITION_DTYPE).ravel() for name in names]
    lengths = np.array([len(a) for a in arrays], dtype=np.int64)
    offsets = np.zeros(len(names) + 1, dtype=np.int64)
    np.cumsum(lengths, out=offsets[1:])
    positions = np.concatenate(arrays) if arrays else np.empty(0, dtype=POSITION_DTYPE)
    width = max((len(n) for n in names), default=1)
    keys = np.array([_encode(n) for n in names], dtype=f"S{width}")
    return keys, offsets, positions.astype(POSITION_DTYPE, copy=False)


def _merge(old, new):
    """Entries of `old` not replaced by `new`, plus all of `new`, in key order."""
    old_keys, old_offsets, old_positions = old
    new_keys, new_offsets, new_positions = new
    if len(old_keys) == 0:
        return new
    keep = ~np.isin(old_keys, new_keys)
    old_lengths = np.diff(old_offsets)
    kept_positions = old_positions[np.repeat(keep, old_lengths)]

    keys = np.concatenate([old_keys[keep], new_keys])
    lengths = np.concatenate([old_lengths[keep], np.diff(new_offsets)])
    positions = np.concatenate([kept_positions, new_positions])

    starts = np.zeros(len(lengths), dtype=np.int64)
    if len(lengths) > 1:
        np.cumsum(lengths[:-1], out=starts[1:])
    order = np.argsort(keys, kind="stable")
    keys, lengths, starts = keys[order], lengths[order], starts[order]
    offsets = np.zeros(len(keys) + 1, dtype=np.int64)
    np.cumsum(lengths, out=offsets[1:])
    gather = np.repeat(starts - offsets[:-1], lengths) + np.arange(offsets[-1], dtype=np.int64)
    return keys, offsets, positions[gather]


def _hold(path: str, replace: bool) -> h5py.File | None:
    """Open the existing index read-write, taking HDF5's exclusive lock.

    A file that cannot be opened for another reason is kept as an error when
    its entries are to be merged, since merging would otherwise silently drop
    them. Under `replace` they are being discarded anyway, and a corrupt file
    is simply written over, as opening it with "w" always did. A LOCK failure
    is raised in both cases: it means another process has the file open.
    """
    if not os.path.exists(path):
        return None
    try:
        return h5py.File(path, "r+")
    except BlockingIOError:
        raise
    except OSError as exc:
        if not replace or "lock" in str(exc).lower():
            raise
        return None


def write_entries(
    path: str,
    entries: Mapping[str, Sequence[int] | np.ndarray],
    *,
    replace: bool = False,
    record_starts: Sequence[int] | None = None,
    attrs: Mapping[str, object] | None = None,
) -> None:
    """Write entries into the index at `path`, in the sorted-blocks layout.

    Existing entries are kept unless `replace`, and an entry named in
    `entries` overwrites the stored one. An index in the per-dataset layout is
    converted on its first write. Existing root attributes and record starts
    survive a merge; `replace` discards them with the entries, as opening the
    file with "w" used to.

    The new file is written beside the old one and moved into place, so an
    interrupted write leaves the previous index intact rather than half
    rewritten. The old file is held open read-write for the duration, which
    takes HDF5's exclusive lock: a second run using the same data directory
    still fails with the lock error `core/concurrent_runs.py` translates
    (Known Issue 21), rather than both writers succeeding and the later one
    silently discarding the earlier one's entries.
    """
    new = _pack(entries)
    held = _hold(path, replace)
    tmp = f"{path}.{os.getpid()}.tmp"
    try:
        old_attrs: dict = {}
        old_starts = None
        merged = new
        if held is not None and not replace:
            existing = PositionIndex(path, handle=held)
            old_attrs = existing.attrs
            old_starts = existing.record_starts()
            merged = _merge(existing._as_arrays(), new)

        keys, offsets, positions = merged
        starts = record_starts if record_starts is not None else old_starts
        with h5py.File(tmp, "w") as out:
            for name, value in old_attrs.items():
                out.attrs[name] = value
            for name, value in (attrs or {}).items():
                out.attrs[name] = value
            out.attrs[LAYOUT_ATTR] = SORTED_BLOCKS
            if starts is not None:
                out.create_dataset(RECORD_STARTS_KEY, data=np.asarray(starts, dtype=np.int64))
            out.create_dataset(KEYS_KEY, data=keys)
            out.create_dataset(OFFSETS_KEY, data=offsets)
            out.create_dataset(POSITIONS_KEY, data=positions, dtype=POSITION_DTYPE)
        os.replace(tmp, path)
    finally:
        if held is not None:
            held.close()
        with contextlib.suppress(FileNotFoundError):
            os.remove(tmp)
