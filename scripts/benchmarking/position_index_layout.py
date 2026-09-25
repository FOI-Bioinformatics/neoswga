"""Per-dataset against sorted-blocks position indexes, on a real index.

Copies a per-dataset index, converts the copy, checks every entry, attribute
and record start agrees, and then times `PositionCache` loading the same
primers from each copy in a fresh process, comparing a digest of everything
it cached.

Usage:
    python scripts/benchmarking/position_index_layout.py \
        <prefix_of_a_per_dataset_index> <k> <step3_df.csv> <workdir> [--all]

`--all` loads every entry in the index instead of the step 3 primers.
Measurements: docs/validation/position_index_layout_2026-09-25.md
"""

import os
import shutil
import subprocess
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.getcwd())
from neoswga.core import position_index as pi  # noqa: E402

LOAD = """
import hashlib, os, resource, sys, time
sys.path.insert(0, os.getcwd())
import numpy as np
from neoswga.core.position_cache import PositionCache
primers = open(sys.argv[2]).read().split()
start = time.perf_counter()
cache = PositionCache([sys.argv[1]], primers)
elapsed = time.perf_counter() - start
digest = hashlib.sha256()
for key in sorted(cache.cache, key=lambda k: k[1:]):
    arr = np.asarray(cache.cache[key])
    digest.update(repr(key[1:]).encode()); digest.update(str(arr.dtype).encode())
    digest.update(arr.tobytes())
scale = 1 if sys.platform == "darwin" else 1024
rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale / 1e6
print(f"load={elapsed:.2f}s arrays={len(cache.cache):,} peakRSS={rss:.0f}MB "
      f"digest={digest.hexdigest()[:16]}")
"""


def main():
    prefix, k, step3, workdir = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
    everything = "--all" in sys.argv
    os.makedirs(workdir, exist_ok=True)
    src = pi.index_path(prefix, k)
    old_prefix = os.path.join(workdir, "old")
    new_prefix = os.path.join(workdir, "new")
    old, new = pi.index_path(old_prefix, k), pi.index_path(new_prefix, k)
    shutil.copy(src, old)
    shutil.copy(src, new)

    start = time.perf_counter()
    pi.write_entries(new, {})
    convert = time.perf_counter() - start

    with pi.open_index(old) as a, pi.open_index(new) as b:
        assert a.layout == pi.PER_DATASET, "the source must be a per-dataset index"
        assert sorted(a.keys()) == b.keys()
        assert a.record_starts() == b.record_starts()
        assert all(b.attrs[name] == value for name, value in a.attrs.items())
        entries = dict(a.items())
        for key, arr in b.items():
            assert np.array_equal(arr, entries[key]), key
        empty = sum(1 for v in entries.values() if len(v) == 0)
        print(
            f"entries={len(entries):,} empty={empty:,} "
            f"per_dataset={os.path.getsize(old) / 1e6:.1f}MB "
            f"sorted_blocks={os.path.getsize(new) / 1e6:.1f}MB convert={convert:.1f}s "
            f"(identical entry for entry)"
        )
        primers = b.keys() if everything else pd.read_csv(step3)["primer"].tolist()

    listing = os.path.join(workdir, "primers.txt")
    with open(listing, "w") as fh:
        fh.write("\n".join(p for p in primers if len(p) == k))
    for label, layout_prefix in (("per_dataset", old_prefix), ("sorted_blocks", new_prefix)):
        result = subprocess.run(
            [sys.executable, "-c", LOAD, layout_prefix, listing],
            check=True,
            capture_output=True,
            text=True,
        )
        print(f"{label:14s} {result.stdout.strip()}")


if __name__ == "__main__":
    main()
