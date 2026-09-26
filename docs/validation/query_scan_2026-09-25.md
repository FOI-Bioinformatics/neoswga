# Counting a known k-mer set without counting the reference

2026-09-25. Measured on the references in
`examples/wolbachia_pool_design/input/`, macOS arm64, Python 3.13,
jellyfish 2.3.1 and KMC 3.2.4.

## What this is

`filter` asks one question of a background: the counts of its own candidate
list. That is a set question, and it was answered by first counting every
k-mer in the reference. `core/query_scan.py` answers it by reading the
reference once instead, with no table built and nothing written to disk.

It is reached from `filter` only where a prefix has **no** table. A counted
prefix takes the table path exactly as before, so no run that works today
changes. The case it serves is a host that cannot practically be counted,
where the alternative was a hard refusal.

## Results

Reproduced with
`scripts/benchmarking/query_scan_against_counting.py <reference> <k> <queries>`.

### Counting is faster wherever it is possible

Drosophila, 144 Mb, 2,000 candidate 12-mers from the shipped design:

| route | time | kept on disk |
|---|---|---|
| scan | 4.0 s | none |
| count with jellyfish | 0.9 s | 33.7 MB table |
| look the batch up in that table | 0.2 s | |

The same at k=18, 1,272 queries:

| route | time | kept on disk |
|---|---|---|
| scan | 5.0 s | none |
| count with jellyfish | 3.0 s | 818 MB table |
| look the batch up in that table | 2.4 s | |

So the scan does not win on time. It wins on what it does not need: no
counter, no table, no disk. At k=12 that is 33.7 MB and the scan is not worth
it; at k=18 the same 144 Mb reference costs 818 MB, and the gap grows with
both the reference and k. hg38 at k=18 is roughly 78 to 84 GB as text, which
is the point at which "count it first" stops being an option.

Peak memory for the scan was 332 MB on Drosophila and 46.6 MB on wMel
(1.27 Mb). It is bounded by the largest RECORD, held once as bytes and once as
codes, plus the chunk. It is **not** proportional to the query set; an earlier
draft of the module said it was, and the first measurement refuted it.

### The chunk size was measured, not chosen

Drosophila at k=12. A larger chunk is both slower and heavier, because the
working set stops fitting in cache:

| positions per chunk | time | peak RSS |
|---|---|---|
| 30,000 | 4.0 s | 312 MB |
| 60,000 | 4.0 s | 329 MB |
| 125,000 | 3.9 s | 351 MB |
| 250,000 | 4.4 s | 400 MB |
| 1,000,000 | 4.6 s | 560 MB |
| 4,000,000 | 6.5 s | 683 MB |
| 8,000,000 | 6.6 s | 1,194 MB |

The default is 125,000, in the flat part of the range rather than at its edge.
The prototype used 8,000,000, which was the worst value measured on both axes.

## Agreement

A fallback that disagrees with the path it replaces is worse than no fallback,
so agreement is checked rather than assumed, three ways.

- Against a **brute-force oracle** written inside
  `tests/test_query_scan_agrees_with_counting.py`, which calls nothing from
  the module: every window built by Python string slicing, canonical form by
  string comparison. It covers overlapping occurrences, a palindrome, a
  reverse complement, record joins, ambiguous bases, lower case, and chunk
  edges at five sizes.
- Against **jellyfish**, and separately against **KMC**, on generated
  references.
- Against the production `kmer_tables.counts_for` on real references: wMel at
  k=12 agreed on all 2,000 counts, Drosophila at k=18 on all 608 canonical
  counts.

### One known difference, and it is not the scan's

Of 1,272 queries at k=18, 13 disagreed, and **all 13 were non-canonical
spellings**. A canonical table stores one spelling of each reverse-complement
pair, so `counts_for` answers 0 for the other one. The scan answers with the
pair's count.

The scan is right and the 0 is not a measurement. It is left alone here
because every caller in this package reads its k-mers from a canonical table,
so no query the pipeline makes is affected -- verified on the bundled plasmid
design, where 0 of 37 step-2 primers carry a non-canonical spelling. Both
behaviours are pinned by a test, so the gap cannot widen unnoticed and a
deliberate fix to `counts_for` would show up as a change there.

## Not measured

- No scan has been run against hg38 through the package. The prototype this
  module was written from scanned hg38 at k=12 in 707 s and at k=18 in 797 s,
  both at about 2.3 GB peak, and agreed with KMC exactly at k=12 on 921,263
  present k-mers and 774,238,630 occurrences. Those figures predate the chunk
  change, which lowered memory substantially on every reference measured here,
  so treat them as an upper bound rather than as this module's numbers.
- The scan is single-threaded. Both counters use several cores, which is part
  of why counting wins on time.
