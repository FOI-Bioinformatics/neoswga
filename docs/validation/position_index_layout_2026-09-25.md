# Position index layout: one dataset per k-mer against sorted blocks

2026-09-25. Measured on the shipped Wolbachia design
(`examples/wolbachia_pool_design/work`) and the bundled plasmid example, on
macOS, Python 3.13, h5py 3.x.

## What changed

A position index (`{prefix}_{k}mer_positions.h5`) used to store one HDF5
dataset per k-mer. It now stores three: the k-mers sorted as bytes, an offsets
array of length n + 1, and every position concatenated in key order. Entry i
is `positions[offsets[i]:offsets[i + 1]]`. Root attributes and
`#record_starts` are unchanged. `core/position_index.py` reads both layouts
and writes only the new one.

The per-dataset layout spent most of its bytes on HDF5 bookkeeping. `h5stat`
on the wMel index reported 376.2 MB of metadata around 5.02 MB of raw
positions: object headers 267.6 MB (145.9 MB of it unused), B-tree nodes
85.6 MB and the local heap 23.1 MB.

## Results

### The wMel index, 983,672 entries of which 446,071 are empty

From `scripts/benchmarking/position_index_layout.py`, each load in a fresh
process. Two runs are shown where they were made, to give the spread.

| | per dataset | sorted blocks |
|---|---|---|
| file size | 381.3 MB | 24.7 MB |
| `PositionCache` load, 2,000-primer shortlist (4,000 arrays) | 0.15 s | 0.02 s |
| peak RSS for that load | 164 MB | 112 MB |
| `PositionCache` load, every entry (1,967,344 arrays) | 37.0 s, 38.5 s | 5.4 s, 5.5 s |
| peak RSS for that load | 968 MB, 1,047 MB | 927 MB, 986 MB |

A SHA-256 over every cached key, dtype and array is identical for the two
layouts in both loads. The script also checks every entry, every root
attribute and the record starts one by one against the original.

The full load barely moves peak memory because `PositionCache` holds about two
million separate numpy arrays whatever the file layout. Most of the remaining
5.4 s is building those objects, not reading the file.

Converting the wMel index took 28.9 s, most of it reading the old layout.
This happens once, on the first write to an index. Reading does not convert,
so an index is converted when `filter` next writes to it and is read in place
until then.

### The plasmid example, run end to end

All four steps were run twice in one directory, on the code before and after
this change. `step2_df.csv`, `step3_df.csv` and `step4_improved_df.csv` were
byte-identical in both rounds. The second round reused the indexes the first
had written, which is the merge path. The summary JSON differs only in its
timing fields. The twelve indexes total 7.9 MB before and 592 KB after.

## Decisions and their reasons

- **Keys are stored as written, never canonicalised.** The reverse strand of
  a primer is read from its reverse complement's entry. An index holding only
  canonical k-mers would lose one strand of every non-palindromic primer.
- **An empty entry and an absent one stay distinct.** Empty means scanned and
  found nowhere. Absent means never scanned. `get` returns an empty array for
  the first and None for the second, and `get_many` leaves absent keys out.
- **`index_format_version` stays 2.** It records whether the sites are
  join-safe (version 2 is the first whose scan respects record joins). A
  layout change alters no site. The layout is named by `index_layout`, and
  readers detect it from the presence of `#keys`.
- **The whole index is rewritten on every write.** The old writer changed
  datasets in place, and HDF5 does not reclaim the space of a deleted
  dataset. The new writer merges in memory, writes a new file beside the old
  one and renames it into place. An interrupted write therefore leaves the
  previous index intact, and repeated writes do not grow the file. Writes
  happen once per prefix and k per `filter` run, so the cost of rewriting is
  bounded by the index size.
- **The old file is held open read-write during a write.** This takes HDF5's
  exclusive lock. Without it, two runs sharing a data directory would both
  succeed and the later rename would silently drop the earlier run's entries.
  Holding it keeps the lock error that `core/concurrent_runs.py` translates
  (Known Issue 21).
- **Reads of the concatenated positions are coalesced.** Requested entries
  closer than 4,096 positions apart are read in one slice, and no slice
  exceeds 16 million positions (128 MB). Each HDF5 read has a fixed overhead
  of tens of microseconds, so reading half a million entries one at a time
  would be dominated by it.
- **No compression.** The layout change alone removed 94% of the file. With
  2-bit key codes the wMel index would be about 17 MB rather than 24.7 MB.
  That was not judged worth a key encoding every reader must share.

## Compatibility with older releases

An older release finds no entry in a sorted-blocks index. Measured: the
previous `optimize` on the converted plasmid directory exits 1 with its step 4
prerequisite error ("37 of 37 candidate primers have no binding positions in
the foreground index"). It does not report a result.

An older `filter` would rescan and write per-k-mer datasets into the
converted file. The new reader refuses such a mixed file (`MixedLayoutError`,
an `OSError`), and the scan path then rebuilds it rather than trusting either
half. Both cases are pinned by tests.

## Not measured

- No index against a whole host genome has been converted. The largest
  measured is wMel at 983,672 entries.
- Only macOS was measured. Linux behaviour of the lock and the rename is
  assumed to match, since both are POSIX. CI covers the tests on Linux but
  not the timings.

## Reproducing

```bash
python scripts/benchmarking/position_index_layout.py \
    examples/wolbachia_pool_design/work/wmel 12 \
    examples/wolbachia_pool_design/work/step3_df.csv /tmp/pidx
python scripts/benchmarking/position_index_layout.py \
    examples/wolbachia_pool_design/work/wmel 12 unused /tmp/pidx --all
```

The source must still be a per-dataset index. The script converts a copy and
never writes to the source.
