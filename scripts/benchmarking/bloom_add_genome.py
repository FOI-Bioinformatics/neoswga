"""Cost of BackgroundBloomFilter.add_genome, before and after the run rewrite.

The old implementation walked every position, validated the k-mer base by base
with `all(b in valid_bases for b in kmer)`, appended it to a list, and then
looped over that list calling `bloom.add` one at a time. pybloom has no bulk
insert, so the batching bought nothing; the validity test is O(k) per position
per k value.

The new one locates maximal ACGT runs once per record and slides inside them,
which makes the validity test free.

Both are timed here on the same sequence, in one process, so the comparison is
not across machine states. Quote the measured figure and the genome it was
measured on. A host-genome figure from these numbers is an extrapolation.

Usage:
    python scripts/benchmarking/bloom_add_genome.py <genome.fna> [k] [max_bp]
"""

import re
import sys
import time

from pybloom_live import BloomFilter

_ACGT_RUN = re.compile(r"[ACGT]+")


def _load(path, max_bp):
    from Bio import SeqIO

    chunks = []
    total = 0
    for record in SeqIO.parse(path, "fasta"):
        seq = str(record.seq).upper()
        chunks.append(seq)
        total += len(seq)
        if total >= max_bp:
            break
    return "".join(chunks)[:max_bp]


def old_scan(seq, k, bloom):
    valid_bases = set("ATCG")
    n_positions = len(seq) - k + 1
    batch = []
    count = 0
    for i in range(n_positions):
        kmer = seq[i : i + k]
        if all(b in valid_bases for b in kmer):
            batch.append(kmer)
            if len(batch) >= 10000:
                for km in batch:
                    bloom.add(km)
                count += len(batch)
                batch = []
    if batch:
        for km in batch:
            bloom.add(km)
        count += len(batch)
    return count


def new_scan(seq, k, bloom):
    add = bloom.add
    count = 0
    for start, end in (m.span() for m in _ACGT_RUN.finditer(seq)):
        stop = end - k + 1
        if stop <= start:
            continue
        for i in range(start, stop):
            add(seq[i : i + k])
        count += stop - start
    return count


def main():
    path = sys.argv[1]
    k = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    max_bp = int(sys.argv[3]) if len(sys.argv) > 3 else 5_000_000

    seq = _load(path, max_bp)
    print(f"sequence: {len(seq):,} bp from {path}, k={k}")

    results = {}
    for name, scan in (("old", old_scan), ("new", new_scan)):
        # A fresh filter each time: an already-populated one short-circuits
        # differently and would not be the same measurement.
        bloom = BloomFilter(capacity=max(1, len(seq) * 2), error_rate=0.01)
        start = time.perf_counter()
        inserted = scan(seq, k, bloom)
        elapsed = time.perf_counter() - start
        results[name] = (elapsed, inserted)
        print(
            f"  {name}: {elapsed:8.3f} s for {inserted:,} positions "
            f"({1e6 * elapsed / max(1, inserted):.3f} us/position)"
        )

    old_t, old_n = results["old"]
    new_t, new_n = results["new"]
    assert old_n == new_n, f"the two scans disagree: {old_n} against {new_n}"
    print(f"  ratio: {old_t / new_t:.2f}x on this sequence at k={k}")


if __name__ == "__main__":
    main()
