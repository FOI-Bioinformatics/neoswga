# A design reads 0.2% of the background table, and Known Issue 3 names the wrong path

Measured 2 November 2026, on `main` at `ecbc67c`. Phase 1.1 of
[the plan](../reports/PLAN_2026-11-02_known_issues.md), against Known Issue 3.
Reproduce with:

```
PYTHONPATH=. python docs/validation/background_table_reads_2026-11-02/measure.py \
    --params <work>/params.json --out results.json --weigh
```

Results and the pool's configuration in
[`background_table_reads_2026-11-02/`](background_table_reads_2026-11-02/).
**No code changed.** Phase 1.1 was a measurement to decide whether to write
Phase 1.2, and it decided against writing it as specified.

## The entry names a path that was fixed a month ago

Known Issue 3 reads, in full: *"Memory usage: The filter command loads all
background k-mers into memory. Use Bloom filter for large backgrounds."*

`filter`'s background-frequency gate does not load the table. It calls
`kmer_tables.counts_for`, which asks KMC to intersect the table with the
candidate list, and `filter._get_rate_for_one_file` carries the measurement
that motivated the change — on *Drosophila* at k=18 the old line-by-line scan
read 116,702,442 lines in 8.01 s where the intersection costs 2.5 s
([`kmer_counter_comparison_2026-09-25.md`](kmer_counter_comparison_2026-09-25.md)).
So the named path has been batched since 2026-09-25 and the entry was not
updated.

Six places still read a whole table. Four of them stream it one pair at a time
and hold nothing: the two Bloom-filter builds in `background_filter`, which
must see every k-mer by construction, one counting internal, and the exclusion
gate in `pipeline`, which breaks on the first match. **Two materialise a
dict:**

| site | purpose | on the design path? |
|---|---|---|
| `mismatch_counts.load_kmer_counts` | the occupancy-weighted site load | **yes**, and hot |
| `kmer_counter.MultiGenomeKmerCounter.count_kmers` | a per-(genome, k) cache | only `multi_genome_pipeline` |

So the real consumer is the occupancy path, and `mismatch_counts`' own
docstring already says so: *"This still materialises the whole table as a dict,
which is its own ceiling on a host genome."*

## How much of the table is read: 0.2%

The occupancy path asks for the exact match plus the 3k one-mismatch
neighbours of every shortlisted candidate, deduplicated by canonical form.
300 candidates, 4 Mb host, `max_mismatches` 1:

| k | candidates | table holds | looked up | fraction |
|---|---|---|---|---|
| 11 | 2 | 1,770,838 | 68 | 0.004% |
| 12 | 298 | 3,158,687 | 10,676 | 0.34% |
| **total** | 300 | **4,929,525** | **10,744** | **0.22%** |

The dict holds about **460x** what the run reads.

## What it costs, and what a batched lookup costs instead

Per k, same table, same forms:

| k | stream whole table | one batched lookup | speedup | dict holds |
|---|---|---|---|---|
| 11 | 2.56 s | 0.60 s | 4.3x | 174 MB |
| 12 | 4.52 s | 1.59 s | 2.8x | 233 MB |

**The counts agree exactly** on every one of the 10,744 forms, which is the
check that matters: a batched lookup returning 0 for a form the table holds
would be a silent wrong answer. `counts_for` guarantees every REQUESTED k-mer
appears in its result, so the hazard is asking outside the batch, not the
batch itself.

Measured memory, holding both tables: 4,929,525 entries for **515 MB**, about
109 bytes an entry. A large host's k=12 table approaches 4^12/2 = 8,388,608
canonical forms — CLAUDE.md records that 99.7% of all canonical 12-mers occur
in the human genome — which projects to about **876 MB for k=12 alone**, and
roughly 1.2 GB across k 10-12. To serve eleven thousand lookups.

## Why Phase 1.2 is not written: the plan mis-specified it

The plan said to "move the callers to `kmer_tables.counts_for`, batched per k".
That cannot be done as stated, and the reason is the call frequency. Counted
over one 12-primer `hybrid` run on this pool:

| | calls |
|---|---|
| `occupancy.weighted_site_load` | 3,616 |
| `mismatch_counts.mismatch_class_counts` | 43,392 |
| `load_kmer_counts` | 43,392 |
| **whole-table reads actually performed** | **4** |

`load_kmer_counts` is `lru_cache`d, so 43,392 calls become four table reads and
43,388 cache hits. **That cache is what makes the per-evaluation path
affordable at all.** A per-call `counts_for` would be 43,392 KMC intersections
at roughly 1.5 s each — about eighteen hours for a run that currently takes
minutes.

So the change is not "move the callers". It is an accumulating per-(prefix, k)
cache that fetches only the forms it has not already seen, in batches — which
would converge on the pool's 10,744 forms after the first few evaluations and
hit thereafter. That is a different and larger piece of work than the plan
assumed, and `mismatch_counts`' docstring under-describes it too when it calls
the fix "a restructuring of the callers rather than of this function".

It also carries a correctness hazard the current code does not have. An
accumulating cache answers `.get(form, 0)` for a form it never fetched, and 0
is a measured background count — the absent-as-zero defect this module's
canonicalisation notes are entirely about. An honest version has to refuse an
unasked form rather than return zero.

**And its only demonstrable benefit cannot be measured here.** The time saved
end to end is four table reads, about 12 s of a run dominated by 43,392
evaluations. The prize is the ~0.9-1.2 GB at host scale, and no host genome is
reachable from this container (`tests/validation/genomes/` is gitignored and
absent, NCBI is refused by the egress proxy). Writing an optimisation against a
43,000-call hot path whose benefit this container cannot show, and whose
failure mode is a silently wrong background count, is the trade this project's
conventions argue against.

## What this changes

- **Known Issue 3 should name the occupancy path, not `filter`'s frequency
  gate**, and can now carry figures: 0.22% of the table read, 515 MB measured
  for 4.9M entries, ~0.9-1.2 GB projected at host scale.
- **The Bloom filter remains the answer at host scale**, which is what the
  entry already says — and Known Issue 22 records that the Bloom path has not
  been built against a host genome here either.
- **If the accumulating cache is written**, the measurement to make first is
  how many batches a real search needs before it stops missing, because that
  decides whether the KMC call overhead stays amortised. The ratio above says
  the memory is worth chasing; it does not say the design is free.

## Scope of the claim

One synthetic pool, one host size, one candidate count, `max_mismatches` 1,
one method, one panel size for the call counts. The host is synthetic, so the
number of DISTINCT k-mers in its table is near the maximum for its length and
the 0.22% fraction is therefore close to the best case for the dict — a real
genome of the same length holds fewer distinct forms, which would raise the
fraction, not lower it. At `max_mismatches` 2 the numerator rises by roughly an
order of magnitude and the fraction with it, which is the regime where this
answer could change. Nothing here is a wet-lab outcome and no delivered panel
moved: no code was modified.
