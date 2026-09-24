"""Resident memory of SampledGenomeIndex, per stored entry.

The Bloom filter exists so a host-sized background need not be held as an exact
index. The sampled index beside it is a plain Python dict keyed by k-mer
string, and at host scale it is the larger of the two, so a user who reached
for Bloom to save memory should know where the memory went.

`ru_maxrss` is a high-water mark that never falls, so one size per process:
run this script once per size and read the per-entry figure off each run.

Usage:
    python scripts/benchmarking/sampled_index_rss.py <n_entries> [k]
"""

import resource
import sys

from neoswga.core.background_filter import SampledGenomeIndex

# macOS reports ru_maxrss in bytes, Linux in kilobytes.
_RSS_SCALE = 1 if sys.platform == "darwin" else 1024


def _rss_bytes():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * _RSS_SCALE


def main():
    n = int(sys.argv[1])
    k = int(sys.argv[2]) if len(sys.argv) > 2 else 12

    baseline = _rss_bytes()
    index = SampledGenomeIndex(sample_rate=100)

    bases = "ACGT"
    for i in range(n):
        # Distinct k-mers without holding a list of them: base-4 of the index.
        value = i
        chars = []
        for _ in range(k):
            chars.append(bases[value & 3])
            value >>= 2
        index.kmers["".join(chars)] += 1

    peak = _rss_bytes()
    delta = peak - baseline
    print(
        f"entries={n:,} k={k} distinct={len(index.kmers):,} "
        f"delta={delta / 1e6:.1f} MB  {delta / max(1, n):.1f} B/entry"
    )


if __name__ == "__main__":
    main()
