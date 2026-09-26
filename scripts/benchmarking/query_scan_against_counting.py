"""Scanning a reference for a known k-mer set, against counting it first.

The two answer the same question and pay for it differently. Counting builds a
table bounded by the reference rather than by the k-mer space above about
k=15, then answers every later batch cheaply. Scanning reads the reference
once per batch, writes nothing, and holds memory proportional to the query set.

This measures both on one reference and CHECKS THE COUNTS AGREE. A speed
comparison between two routes that disagree is worthless, and the ways they
could disagree are real: canonical form, windows spanning a record join, and
windows holding an ambiguous base are each a decision both routes have to make
the same way.

Usage:
    python scripts/benchmarking/query_scan_against_counting.py \
        <reference.fna> <k> <queries.txt> [--skip-count]

`queries.txt` is one k-mer per line, or any file whose first whitespace field
is a k-mer, so a `*_all.txt` k-mer table works as input. `--skip-count` scans
only, for a reference whose table would not fit.

Measurements: docs/validation/query_scan_2026-09-25.md
"""

import os
import resource
import sys
import tempfile
import time

sys.path.insert(0, os.getcwd())

from neoswga.core import kmer_tables, query_scan  # noqa: E402

_RSS_SCALE = 1 if sys.platform == "darwin" else 1024


def _peak_rss_bytes() -> int:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * _RSS_SCALE


_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _is_canonical(kmer: str) -> bool:
    return kmer <= kmer.translate(_COMPLEMENT)[::-1]


def _read_queries(path: str, k: int) -> list[str]:
    queries = []
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if fields and len(fields[0]) == k:
                queries.append(fields[0].upper())
    return queries


def main() -> None:
    reference, k, query_file = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    skip_count = "--skip-count" in sys.argv

    queries = _read_queries(query_file, k)
    if not queries:
        raise SystemExit(f"{query_file} holds no {k}-mers")
    size = os.path.getsize(reference)
    print(
        f"reference={os.path.basename(reference)} ({size / 1e9:.2f} Gb) "
        f"k={k} queries={len(queries):,}"
    )

    before = _peak_rss_bytes()
    start = time.perf_counter()
    scanned = query_scan.count_kmers(reference, k, queries)
    scan_seconds = time.perf_counter() - start
    scan_rss = max(_peak_rss_bytes(), before)
    present = sum(1 for count in scanned.values() if count)
    print(
        f"  scan      {scan_seconds:8.1f}s  peakRSS={scan_rss / 1e6:8.1f}MB  "
        f"disk=0MB  present={present:,}  hits={sum(scanned.values()):,}"
    )

    if skip_count:
        return

    with tempfile.TemporaryDirectory(prefix="query_scan_bench_") as workdir:
        from neoswga.core import kmer_counter

        prefix = os.path.join(workdir, "ref")
        start = time.perf_counter()
        kmer_counter.run_jellyfish(reference, prefix, min_k=k, max_k=k, cpus=4)
        count_seconds = time.perf_counter() - start

        start = time.perf_counter()
        counted = kmer_tables.counts_for(prefix, k, queries)
        lookup_seconds = time.perf_counter() - start
        table_bytes = sum(
            os.path.getsize(os.path.join(workdir, name)) for name in os.listdir(workdir)
        )
        print(
            f"  count     {count_seconds:8.1f}s  table={table_bytes / 1e6:8.1f}MB\n"
            f"  lookup    {lookup_seconds:8.1f}s  (per batch, once the table exists)"
        )

    # Split the comparison, because one difference between the two routes is
    # known and the other would be a defect.
    #
    # A canonical table stores one spelling of each reverse-complement pair, so
    # `counts_for` answers 0 for the other spelling. The scan answers with the
    # pair's count. Every caller in this package reads its k-mers from such a
    # table, so every query the pipeline makes is canonical and the two agree;
    # a disagreement on a CANONICAL query would mean they are counting
    # different things, and that is what fails here.
    canonical = [q for q in queries if _is_canonical(q)]
    disagreements = [q for q in canonical if scanned[q] != counted[q]]
    if disagreements:
        example = disagreements[0]
        raise SystemExit(
            f"the two routes disagree on {len(disagreements):,} of "
            f"{len(canonical):,} canonical k-mers, so the timings describe "
            f"different quantities; {example}: scan={scanned[example]} "
            f"table={counted[example]}"
        )
    print(f"  AGREE on all {len(canonical):,} canonical counts")

    other = [q for q in queries if not _is_canonical(q) and scanned[q] != counted[q]]
    if other:
        print(
            f"  note: {len(other):,} of {len(queries) - len(canonical):,} "
            f"non-canonical spellings differ, where the table answers 0 and "
            f"the scan answers the pair's count. The pipeline asks only "
            f"canonical k-mers, so it never sees this."
        )


if __name__ == "__main__":
    main()
