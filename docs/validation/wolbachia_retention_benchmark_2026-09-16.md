# Candidate retention on the Wolbachia/Drosophila pair

Measured 2026-09-16. Plan steps 249 and 250 of
`docs/superpowers/plans/2026-09-15-condition-aware-pool-design.md`.

Reproduce with `scripts/benchmarking/wolbachia_retention_benchmark.py` from the
repository root. Raw numbers are in
`examples/wolbachia_pool_design/benchmark_2026-09-16/results.json`.

**`legacy` was removed after this measurement.** It is kept here because it is
the measurement that justified removing it, and because the figures below are
the only published record of what the historical behaviour cost. The script now
compares `post_gini` against `all_qc`; re-running it will not reproduce the
`legacy` column.

## What was run

`count-kmers` and `filter` for both `candidate_retention` modes, into separate
output directories, from the same `examples/wolbachia_pool_design/params.json`:
wMel `NC_002978.6` as the target, Drosophila `GCF_000001215.4` as the
background, k = 12, phi29 at 30 C, `max_primer` 2000, `max_gini` 0.7. Each mode
counted its own k-mer tables rather than reusing the other's, so the counting
time is a real measurement in both rows.

Runtime is wall clock around the subprocess. Peak memory is the maximum
resident set size reported by `/usr/bin/time -l`.

The run stops after `filter`. It does not exercise `score`, `optimize` or the
report, so it measures the cost of retention rather than any change in a
delivered panel.

## Stage counts

The three counts the plan asked to re-check recur exactly after the chemistry
corrections, and they are identical in both modes:

| Stage | Count |
|---|---|
| Distinct canonical 12-mers counted | 874,596 |
| After the foreground frequency gate | 874,596 |
| After the background frequency gate | 706,756 |
| After thermodynamic QC | 491,836 |
| After the Gini gate | 20,670 |
| After the `max_primer` cut | 2,000 |

The foreground gate removes nothing here because `min_fg_freq` is 1e-6 and the
target is 1.27 Mb, so a single occurrence already clears it.

## What the retention mode changes

| Quantity | `legacy` | `all_qc` |
|---|---|---|
| `filter` runtime, seconds | 85.0 | 160.4 |
| `filter` peak resident set, MB | 2650 | 2761 |
| Candidates indexed | 2,000 | 491,836 |
| Background sequences indexed | 2,000 | 491,836 |
| Background index, MB | 1.7 | 443.1 |
| Foreground index, MB | 381.3 | 381.3 |
| Inventory, MB | 190.2 | 189.4 |
| `step2_df.csv` rows | 2,000 | 2,000 |

Two of these rows carry most of the information.

The **foreground index is the same size in both modes**. It already held every
candidate that cleared hard QC: 983,672 datasets, one per sequence and one per
reverse complement, for 491,836 candidates. Only the background index was cut
to the shortlist. That asymmetry is what the retention flag now addresses.

The **delivered shortlist is unchanged**. Step 2 writes the same 2,000 rows in
both modes, so `all_qc` is not a different design, it is the same design with
the rest of the inventory still addressable.

## The measurement the legacy index cannot make

Under `legacy` the background index holds the 2,000 shortlisted candidates.
Every other candidate is absent from it. `PositionCache.get_positions` answers
an absent sequence on an indexed prefix with an empty array, which is a
plausible measurement of zero rather than an error, so a candidate reached by
expansion reads as having no host binding at all.

Counting against the full index shows how large that is:

| Set | Candidates | Host sites the legacy index reports as zero |
|---|---|---|
| Cleared hard QC, outside the shortlist | 489,836 | 8,334,025 |
| Cleared the Gini gate, outside the shortlist | 18,670 | 390,482 |

963,931 of the 979,672 hard-QC datasets outside the shortlist have at least one
real host site. The second row is the conservative reading: even an expansion
policy that never looks past the Gini gate reaches 18,670 candidates whose host
load the legacy index cannot see.

This is the same failure shape as Known Issues 5, 6 and 13: a quantity that
reads zero where the true value is not zero, with nothing in the output to say
so. It is the reason the retention default changed rather than a performance
argument.

## The intermediate retention point

The plan asks for three modes: the historical cap, an expanded post-Gini
inventory, and all of hard QC. When this was run only two existed in the code,
so only two were run.

Its index cost is estimated here rather than measured. Fitting bytes against
dataset count and site count across the two measured background indexes gives
about 403 bytes per indexed sequence and 5.6 bytes per site, so per-sequence
HDF5 overhead dominates. At the 20,670 post-Gini candidates, 41,340 datasets
and 414,096 sites, that predicts about 19 MB. This is an interpolation between
two points, not a measurement, and it says nothing about runtime.

`post_gini` was added on 2026-09-16 as a real mode. The distinction it draws is
between a gate and a ranking: the evenness gate is a declared requirement a
candidate either meets or does not, while `max_primer` is a cut through a
ranking chosen for the size of the working set. Keeping the first and dropping
the second retains what was set aside arbitrarily and not what was set aside on
a stated rule. It does not reopen the silent zero, because the candidate
provider does not expand past a hard gate.

## Why `legacy` was removed

Its only distinguishing property is the measurement fault in the section above.
It existed to reproduce the historical truncation for comparison against an
existing result, that comparison is the table in this document, and a mode kept
only so a known fault stays reachable is how the fault comes back. A
configuration still naming it is refused with a message saying what replaced it
and why, rather than with a bare schema error.

## Limitations

- Two modes, one organism pair, one length (k = 12), one chemistry, one
  machine. Single runs, not repeated timings, so the runtime figures carry no
  spread.
- The benchmark stops at `filter`. Whether the retained inventory changes a
  delivered panel is a separate question this run does not answer.
- The 443 MB background index is the cost at one length on a 144 Mb background.
  It grows with both, and the plan's per-length grid multiplies it.
- Peak memory is the peak of the whole `filter` process, so it includes the
  k-mer table load. The 4 percent difference between modes is not an isolated
  measurement of the indexing step.
- Host site counts above are index contents, not an independent recount from
  the FASTA. They agree with the jellyfish tables by construction.

## Artifacts

The output directories held 1.7 GB, most of it k-mer tables and HDF5 indexes.
Those are regenerable and were deleted. What is committed is `results.json`,
which carries every number quoted above, and the two run logs per mode.

The per-run `step2_df.csv`, `filter_stats.json`, `params.json` and
`run_manifest.json` are covered by existing `.gitignore` rules for pipeline
artifacts and were left uncommitted rather than forced in. `results.json`
already holds the funnel counts and the step-2 row count, so nothing quoted
here depends on them.
