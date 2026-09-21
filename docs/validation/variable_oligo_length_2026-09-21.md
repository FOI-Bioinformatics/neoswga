# Mixing oligo lengths: it works, and on this pair it bought nothing

Measured 21 September 2026. Prompted by an external review claiming neoswga is
"restricted to fixed k-mer lengths (6-12 bp)".

## The claim is wrong twice over

`params.schema.json` admits k of 4 to 30, and 6-12 is phi29's default window
rather than a limit: bst is 15-25 and equiphi29 10-18. Nothing in the pipeline
requires one k. The scan writes `{prefix}_{k}mer_positions.h5` per length,
`PositionCache.load` groups its primers by length and opens the matching file,
and the dimer screen codes t-mers so unequal lengths compare with no special
handling.

What was missing was a demonstration. Every saved run in this repository used
one k, so "mixed length works" was an argument, and an argument is what a
reviewer is entitled to disbelieve.

`tests/integration/test_variable_oligo_length.py` now settles it through a real
four-step run on the packaged plasmid pair at k 7-11. The delivered panel holds
three 7-mers, four 8-mers and one 10-mer, carries no pair above the configured
`max_dimer_bp` of 3, and survives export, interpret and report.

The reach is shortened to 300 bp in that test on purpose. At phi29's realistic
3 kb one primer covers the 6 kb plasmid completely and Stage 1 returns after a
single pick, which is why this was never visible on the bundled example. It is
the same limitation Known Issue 14 records for `expand-primers`.

## Does mixing help? Not here

Prevotella melaninogenica (3,168,282 bp) against human chr21 (46,709,983 bp),
equiphi29 at 42 C, Tm window 42-54, `max_primer` 5000, seed unset. Four full
pipeline runs, identical but for `min_k`/`max_k`.

Request 8 primers:

| design | panel lengths | coverage | effective | density | host sites |
|---|---|---|---|---|---|
| mixed k 10-12 | 6x11, 2x12 | 0.5309 | 0.4832 | 20.55 | 17 |
| k=10 only | 8x10 | 0.2847 | 0.1667 | 4.10 | 80 |
| k=11 only | 8x11 | **0.5409** | 0.4508 | 16.12 | 26 |
| k=12 only | 8x12 | 0.5083 | 0.4814 | **29.87** | **1** |

Request 20 primers:

| design | panel lengths | coverage | effective | density | host sites |
|---|---|---|---|---|---|
| mixed k 10-12 | 12x11, 8x12 | 0.6374 | 0.5693 | 19.69 | 24 |
| k=10 only | 20x10 | 0.4654 | 0.2841 | 3.89 | 221 |
| k=11 only | 20x11 | **0.6688** | 0.5591 | 13.05 | 56 |
| k=12 only | 20x12 | 0.6025 | 0.5637 | **24.83** | **5** |

The mixed design lands BETWEEN the single-length designs on both axes, at both
sizes. It is beaten on coverage by the all-11-mer panel and beaten on
specificity by the all-12-mer panel, which binds the host once at n=8 against
the mixed panel's seventeen.

That is what a union pool should be expected to do, and it means **choosing k
matters far more than choosing whether to mix**. The gap between k=10 and k=12
is larger than anything mixing recovers: at n=20 the all-10-mer panel reaches
0.4654 coverage at density 3.89 against 0.6025 at 24.83.

The one axis where the mixed panel leads is occupancy-weighted effective
coverage, at both sizes: 0.4832 against 0.4814 and 0.5693 against 0.5637. That
is about 1% relative, from one run per cell with no seed sweep, and should not
be read as a result.

The optimizer does select across lengths when given the chance -- the mixed
panels are genuinely mixed rather than collapsing to one k -- so this is a
finding about the chemistry and the pool, not about selection refusing to mix.

## The occupancy spread, and the control that changed the conclusion

Occupancy is the fraction of the time a site is bound, and in the raw k-mer
space it rises very steeply with length. Median over 3,000 random k-mers per
length:

| k | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|
| phi29 30 C | 0.017 | 0.145 | 0.661 | 0.960 | 0.996 | 1.000 |
| equiphi29 42 C | 0.001 | 0.005 | 0.032 | 0.205 | 0.611 | 0.917 |

The first draft of this document stopped there and concluded that a
mixed-length panel is inherently unbalanced. **A control refuted that.** How
much of the spread survives into a real candidate pool is set by the Tm window,
not by the length range:

| pool | Tm window | median Tm spread across lengths | occupancy range |
|---|---|---|---|
| plasmid, phi29 30 C | 15-45, so 30 C | 23.5 C | 0.042 - 0.994 |
| Prevotella, equiphi29 42 C | 42-54, so 12 C | 2.0 C | 0.580 - 0.794 |

A Tm gate selects on the very quantity occupancy is computed from. The wide
window admits a 7-mer bound 4% of the time beside an 11-mer bound 99%; the
narrow one admits oligos comparable at every length. So the lever is the Tm
window at least as much as the reaction temperature, and a mixed-length design
under a narrow window is not unbalanced at all.

Length is not the only driver either. Two 8-mers in the plasmid panel sit at
0.066 and 0.693, which is composition.

What the per-length view adds over the pool-wide discrimination mean `filter`
already printed: on the plasmid pool that mean is 2.405 while the per-length
values run 3.010 down to 1.015, so the single number describes no length in the
pool.

`core/length_occupancy.py` reports this from `filter` and from `optimize`, and
warns when a weak length sits beside a saturated one. It is reported and never
enforced, the resolution Known Issue 17 reached for the same quantity after
measuring that an occupancy gate makes delivered panels worse on both axes. A
test walks the package and fails if any caller reads the threshold. It is
silent on a single-length design, so no existing run gains output.

## The locking failure was a concurrent run, not mixed length

`tests/validation/genomes/f_mixed.log` records a `BlockingIOError: [Errno 35]
unable to lock file` from a mixed-length run. Mixed length was the suspect
because it is what the config changed.

| runs | attempts | reproduced the error |
|---|---|---|
| two processes, one data directory | 11 | 10 |
| one process, including multi-k | 11 | 0 |

Two independent sets of trials, agreeing. A mixed k 10-12 `filter` over the
same Prevotella and chr21 pair, in a clean directory, completed in 146 s with
exit status 0, and a Prevotella-against-Wolbachia run at k 10-12 with
`max_primer` 5000 completed in 253 s. The recorded traceback differs from the
reproduction only in `h5f.open` against `h5f.create`, which is whether the
target file already existed, and its frames are the sequential Aho-Corasick
branch, so the per-k multiprocessing fallback that was first suspected did not
run.

**The errno is cross-process by construction**, which is what makes the
conclusion firm rather than statistical. A foreign process holding the file,
even read-only, produces exactly the recorded message. A second handle inside
ONE process produces something else, `OSError: ... file is already open for
read-only`, with no errno 35. So a leaked handle within a single run cannot
produce this message at all.

**A reader is enough to stop a writer**, and that is the shape a user will
actually hit: a command that only reads the index can stop a `filter`. The
error message says so, because "another writer" would send someone looking for
a second filter.

How wide the window is depends on which cache is in use, and a first version of
this section overstated it by naming `optimize` flatly. The default
`PositionCache` opens each file inside a `with` block and closes it promptly,
so the collision is a race rather than a certainty. `StreamingPositionCache`
keeps its handles until `close()` and holds them for a whole run, but it is
selected only when the in-memory cache is disabled, which is not the default.

**The recorded directory shows the collision directly.** Its
`run_manifest.json` records `score` on that same config completing 2.2 s before
the failing filter's last log write and `optimize` completing 8.7 s after, with
no `filter` entry at all, because the manifest is written on step completion.
Beside them sits a `prevotella_13mer_positions.h5` of 800 bytes holding zero
datasets, with the same mtime, for a k that run never requested -- a third
process still running under a config edited minutes earlier.

**Which process held the lock is not determined.** The manifest proves that
runs overlapped; it does not say which handle blocked the filter. The `score`
that completed 2.2 s earlier, the `optimize` 8.7 s later and the third process
that left the stray 13-mer index are all candidates, and nothing in the
artifacts separates them. The correction above makes this MORE open rather than
less: with the default cache holding each file only briefly, a collision is a
race, so the holder is whichever process happened to have the file open in that
instant rather than whichever one ran longest.

That limit is worth stating because the cause is settled and the mechanism is
not fully reconstructed, and those are different claims. Nothing downstream
depends on the answer: the remedy is the same whichever process it was.

The remedy is a message, not a change to the scan. `core/concurrent_runs.py`
translates the error at the step boundary; steps 2, 3 and 4 each consult it.
Verified against a real collision: four concurrent filters over one directory,
three lost the lock, and all three printed the message naming the directory and
the cause. Nothing takes a lock or retries, because a retry loop would hide a
genuine second process.

## Caveats

One organism pair, one reaction, one Tm window, two panel sizes, one run per
cell. The mixed pool's shortlist was dominated by 11- and 12-mers (2,056 and
2,896 against 48 10-mers), so "mixed" here is largely an 11-and-12 comparison
rather than a wide one. A target whose GC content makes one length clearly
unsuitable would be the interesting next case, as would a wide Tm window, where
the occupancy spread above says the lengths are least comparable.

Nothing here says mixing cannot help. It says it did not help on the one pair
measured, and that the choice of k dominated it.
