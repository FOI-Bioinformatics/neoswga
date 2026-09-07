# Audit: the four pipeline steps, for speed, robustness and pool quality

**Date:** 2026-09-06
**Branch:** `audit-alternatives-and-scaling`
**Scope:** `count-kmers`, `filter`, `score`, `optimize`, and what falls between them.

This audit asks a different question from
[AUDIT_optimize_step_2026-09-03.md](AUDIT_optimize_step_2026-09-03.md), which
hunted defects. Here the question is where the pipeline spends time it does not
need to spend, where it can fail without saying so, and where the delivered
oligo pool is worse than the candidate pool allows. Findings the earlier audit
already recorded are not repeated; the four per-step reports dated 2026-03-31
are superseded by this one.

## Method

Everything quantified below was measured on this machine, not estimated.

A fresh four-step run: *Wolbachia* (1.27 Mb, circular) as foreground against
human chr21 (46.7 Mb) as background, k=12, `max_primer` 3000, `num_primers` 12,
`max_sets` 3, EquiPhi29 at 42 C, `cpus` 8.

Three existing runs under `runs/gc_tiers/` supplied the funnel counts, the
delivered panels and a set-size sweep: *S. aureus* (32.9% GC), *E. coli*
(50.8%), *M. tuberculosis* (65.6%), each against whole hg38.

Per-pair and per-primer costs were timed directly against the installed
package.

Not covered here, and worth a later pass: the `-C` canonical convention in
`count-kmers` against the strand handling in the position scanner, where a
factor-of-two disagreement would mis-set every frequency threshold; the network
optimizer's two stable attractors at 0.778 and 0.756; and which of the 22
findings in the September 3 audit remain open.

Two smaller items surfaced and are recorded here rather than given sections:
`--seed` and `--full-score` log work they do not do (`cli/pipeline.py:384` and
`:410`), and `validate_step3_prerequisites` does not check that `step2_df.csv`
carries a `primer` or a `gini` column.

## Where the time actually goes

| Step | Wall | Step's own figure | Gap |
|---|---|---|---|
| `count-kmers` | 5.09 s | 3.3 s | process start-up |
| `filter` | 30.37 s | 28.9 s | |
| `score` | 1.71 s | 0.2 s | 1.5 s of start-up around 0.2 s of work |
| `optimize` (hybrid) | 30 min, abandoned unfinished | | see A1 |
| `optimize` (dominating-set, same pool) | 36.61 s | 27.7 s | |

The hybrid run was stopped at thirty minutes on its third of three
heterodimer-screening passes, having selected nothing. `dominating-set` on the
identical 3000-candidate pool returned a 12-primer panel at 0.816 foreground
coverage in 36.6 s.

`neoswga --help` alone costs 2.25 s, so the four-invocation workflow pays that
import four times.

The headline is the last two rows. On an identical candidate pool the default
method had not finished in thirty minutes what the fast method finished in
twenty-eight seconds, and the reason is a single screen, described in A1.

---

## A. `optimize`

### A1 (CRITICAL, speed) An all-pairs thermodynamic screen runs before selection, and again for every alternative set

`hybrid_optimizer.py:485` runs a pre-stage over the whole candidate pool, and
`:1179` asks it for heterodimers unconditionally:

```python
filtered, stats = thermo_filter.filter_candidates(
    candidates,
    check_heterodimers=True,
    max_heterodimer_fraction=0.3,
)
```

`thermodynamic_filter.py:314-319` then materialises every pair at once:

```python
pairs = [
    (passing[i].sequence, passing[j].sequence, i, j, conditions_dict)
    for i in range(len(passing))
    for j in range(i + 1, len(passing))
]
```

One pair costs 627 us through `secondary_structure.check_heterodimer`, measured
over 800 pairs of the real *E. coli* candidate pool on an otherwise idle
machine. That gives, serially:

| pool | pairs | serial cost |
|---|---|---|
| 449 | 100,576 | 1.1 min |
| 1,222 | 746,031 | 7.8 min |
| 2,789 | 3,887,866 | 40.6 min |

**These figures were revised down on 2026-09-06.** The first measurement of this
audit reported 1372 us a pair and 88.9 min at n=2789. It was taken while five
audit subagents were running concurrently on the same machine, and it is
inflated by roughly a factor of two. The re-measurement above was taken on a
quiet machine against the same real pool. The conclusion does not change: the
screen still costs tens of minutes per pass where the method without it finishes
in seconds, and the run that failed to complete in thirty minutes was observed
under the same load as the inflated timing, so that observation stands as
recorded.

`dimer.is_dimer_fast`, the cheap substring test, was likewise revised from 105 us
to 12.9 us a pair. The ratio between the two tests is what matters for the
proposed pre-screen, and it is 49x.

**The attribution in this finding was wrong, and the measurement that settled it
is in [docs/validation/optimizer_cost_2026-09.md](validation/optimizer_cost_2026-09.md)
(2026-09-07).** This finding treated the thermodynamic screen as the reason
`hybrid` fails to finish where `dominating-set` finishes in seconds. It is not.
After the screen was made exact and vectorised, stopped from repeating once per
alternative set, and given a pre-screen that discards 60.1% of pairs, `hybrid`
on the 1215-candidate *S. aureus* pool went from 2499.91 s to 2392.31 s. That is
4.3%, inside run-to-run variation on this machine, against 13.05 s for
`dominating-set` on the same input: a factor of 183 rather than the comparable
time this finding implied was within reach.

Summing the per-set `Total runtime` lines shows where the time actually goes:
Stage-2 network refinement, which none of the work above touched. The screen
fixes are real and verifiable in the logs -- five screens became one, and the
pair count is now reported as `188170 of 471906 pairs among 972 primers` -- but
they address a minor term.

The error was one of attribution, not of measurement. The screen does cost what
this finding says it costs. What was never established, and what the thirty
minute observation cannot establish, is that the screen was the dominant term.
A cost model built from one component's timing needs the other components timed
before it can name a bottleneck.

Resident memory during the screen was 2.0 GB for the eager pair list alone.

The screen is not run once. `unified_optimizer.py:254-258` re-enters
`optimizer.optimize(remaining, target_size)` for each alternative set, up to
`iterations` times (default 8), and each re-entry repeats the whole screen on a
pool one set smaller. The observed run, configured for `max_sets` 3, logged
three passes at 2789, 2777 and 2765 primers. A run at the default `max_sets` 5
would pay up to nine.

Scope: the pre-stage is gated on `poly_config.thermo_filter`, which
`registry/polymerases.py` sets True for EquiPhi29, Bst 2.0 and Bst 3.0 and
False for Phi29, Bsu and Klenow. So Phi29 runs escape it and EquiPhi29 runs do
not, and the GC-adaptive strategy selects EquiPhi29 on its own for a balanced
genome, which is how a user who never chose it lands in the expensive path.
`background-aware` wraps the same `HybridOptimizer` and pays the same cost;
`dominating-set` never constructs it.

This is the mechanism behind the 260x hybrid-versus-dominating-set gap the
project already measured, and behind the Jaccard 1.000 that accompanied it: the
screen removed 16 hairpin primers out of 2805 and did not change the answer.

**Proposed change (M).** Compute the pairwise verdict once for the pool and
cache it. Every later pool is a subset of the first, so cached verdicts are
reusable and nine passes become one. Then put a cheap test in front of the
expensive one: `dimer.is_dimer_fast` costs 12.9 us against 627 us, a 49x
reduction, and only pairs it flags need the thermodynamic calculation. Even the
cheap test should be vectorised, since a longest-common-substring test on
12-mers has no business costing 12.9 us at scale.

(This paragraph originally quoted 105 us, 1372 us and 13x. Those are the
concurrent-load figures this section revises above; the ratio, which is what the
proposal turns on, is 49x rather than 13x.)

Two caveats on the caching half, both established after this audit was written.
The reuse is not exact: `filter_candidates` removes heterodimer hubs at a
threshold proportional to pool size, so reusing a full-pool verdict on a subset
under-removes. Set 0 is computed from the first pool and is unaffected;
alternative sets 1 to 4 can differ. See
`hybrid_thermo_screen.ThermoScreenMixin._thermo_filter_with_cache`.

**Verify.** Time `optimize -m hybrid` on the *S. aureus* pool (1222 candidates)
before and after, and assert the delivered set is unchanged.

That last criterion is settled and it did not hold: the delivered `hybrid` set
changed (Jaccard 0.352 on *S. aureus*), by design, because the dimer rejection
guard added in the same branch rejects primers the unguarded greedy accepted.
The unchanged-set check applies to the caching and pre-screen work alone; see
[docs/validation/optimizer_cost_2026-09.md](validation/optimizer_cost_2026-09.md).

### A2 (HIGH, pool quality) The expensive screen tests a criterion the user never configured, and the criterion the user did configure reaches no optimizer

`hybrid_optimizer.py:1147` hardcodes `max_heterodimer_dg=-10.0`;
`thermodynamic_filter.py:93` defaults it to -9.0. Neither is derived from
`max_dimer_bp`, the parameter a user sets, which is consumed only by
`dimer.is_dimer` as a longest-common-substring length.

`dominating_set_optimizer.py` and `dominating_set_adapter.py` contain no
occurrence of the string "dimer" at all. That is the method every panel in
`runs/gc_tiers/` was built with.

The result reaches the delivered product. The *E. coli* panel's summary records

```json
"worst_heterodimer": {"length_bp": 11, "pair": ["ATCAGCGCCAGA", "TCTGGCGCTGAA"]}
```

against a configured `max_dimer_bp` of 3. The fresh Wolbachia panel from
`dominating-set` records 10 bp against the same configured 3. Both shipped with
a warning.

So the pipeline can spend an hour computing a dimer criterion, deliver a panel
selected without reference to any dimer criterion, and report it against a
third criterion that it then violates.

**Proposed change (M).** Make `max_dimer_bp` the single dimer criterion, thread
it into the pre-stage, and give the greedy in `dominating_set_optimizer` a
rejection test against the already-selected set. That test is cheap because it
runs against the growing set, not the pool: selecting 160 primers costs at most
12,720 comparisons, against 4.5 million for the pool.

### A3 (HIGH, speed) Background pruning recomputes whole-set coverage once per candidate per removal

`hybrid_optimizer.py:1263-1266`, inside the `while len(current_primers) >
target_size` loop:

```python
for primer in current_primers:
    test_set = [p for p in current_primers if p != primer]
    test_coverage = self._calculate_coverage(test_set)
```

`_calculate_coverage` at `:1050` rebuilds a fresh `BipartiteGraph` over every
primer in the set on every call, re-reading positions and re-binning them.
Pruning 160 primers to 24 is 136 removals, each evaluating roughly 90
candidates, each rebuilding over roughly 90 primers: on the order of a million
primer-binning operations to make 136 decisions.

**Proposed change (M).** Hold a per-bin coverage count. Removing a primer
decrements only the bins it occupies, so evaluating a candidate costs the
number of that primer's own sites rather than the whole set's. This is
arithmetic-identical, not an approximation.

Scope: `background_pruning` is enabled for `background-aware`, the method
documented as slow and intended for clinical work, which is where a user is
least willing to wait.

### A4 (MEDIUM, pool quality) The score the tool reports is not the rule by which the set was chosen

`base_optimizer.py:410-427` composes `normalized_score` from coverage,
selectivity, dimer safety, gap evenness and Tm tightness, weighted by
`--application`. That is a defensible objective, and it is what the ensemble
compares methods on and what the summary reports.

It is not what any optimizer maximises. `dominating_set_optimizer.py:603-611`
states the program it solves:

```
maximise    sum_b w_b * y_b
subject to  y_b <= sum_{i covers b} x_i     for every bin b
            sum_i x_i <= S
```

Weighted foreground bin coverage, subject to a set-size budget. Nothing else.
The greedy path is the same rule in cheaper form: it ranks candidates on
`len(new_regions)`, the count of newly covered bins, and stops when no candidate
adds one.

So the delivered set is the argmax of foreground coverage alone, and
`normalized_score` is a description applied to it afterwards. Selectivity,
dimer safety, evenness and Tm tightness appear in the number the user reads and
in no decision the tool makes. A2's 10 and 11 bp heterodimers are the visible
consequence.

**Proposed change (L).** Either say plainly in the output that
`normalized_score` describes the result rather than selecting it, or move the
composite into the greedy's ranking so the two agree. The first is a
documentation change and honest. The second is the real fix, and it is where
the pool quality is.

### A5 (MEDIUM, pool quality) Set size is the user's most consequential choice and the tool does not help with it

Existing sweep, `runs/gc_tiers/*/sweep/`, `dominating-set`:

*M. tuberculosis*

| n | fg coverage | selectivity density |
|---|---|---|
| 8 | 0.463 | 417 |
| 16 | 0.692 | 413 |
| 32 | 0.881 | 421 |
| 40 | 0.923 | 407 |
| 48 | 0.949 | 392 |
| 64 | 0.975 | 364 |

*E. coli*

| n | fg coverage | selectivity density | synthesis USD |
|---|---|---|---|
| 12 | 0.407 | 54.8 | 132 |
| 32 | 0.627 | 37.8 | 352 |
| 96 | 0.890 | 34.2 | 1056 |
| 160 | 0.943 | 33.1 | 1760 |

Coverage rises monotonically while selectivity density peaks near n=32 and
then decays, so there is a real knee and it is not where coverage alone would
put it. Two flags exist that sound like they address this, `--auto-size` at
`cli/pipeline.py:717` and `--show-frontier` at `:982`. Neither appears in
CLAUDE.md, the README quick start or the parameter reference, and this
project's own GC-tier sweep was performed by hand with repeated `-n` values.

Neither would have helped if it had been found.

`--auto-size` cannot see the knee. `recommend_set_size`
(`set_size_optimizer.py:883`) calls `estimate_optimal_set_size` at `:788`,
which inverts a single saturation curve for the application profile's coverage
target, applies a fixed multiplier of 0.85, 1.1 or 1.0 for priority, and clamps
to the profile's `typical_size`. It reads no candidate pool, no binding
positions and no background counts, and `min_fg_bg_ratio` is returned in its
result dict without entering the calculation. The clamp tops out at 20 primers,
so it cannot even name the n=32 where selectivity density peaks for
*M. tuberculosis*.

`--show-frontier` does weigh both axes and does read the real pool and the
position cache, so it is the right tool. But `cli/pipeline.py:1036` passes
`max_size=min(20, len(primer_pool))`. Between them the two flags cannot look
past 20 primers, which is the concrete reason a sweep to n=160 was done by hand.

`--auto-size` is also unreachable from the one-shot entry point. `run_design`
fixes `"auto_size": False` in the defaults namespace it builds
(`cli/commands.py:590`), and `auto_size` is not a params.json key.

**Proposed change (S to M).** Let `design` pass both flags through, lift the
20-primer ceiling on `--show-frontier`, document both where a user will meet
them, and describe `--auto-size` honestly as an estimate of how many primers
reach a coverage target rather than as a recommendation that weighs
specificity. Separately, have `optimize` print the marginal foreground coverage
per primer at the delivered size, which is the only one of the three that has
no size limit.

### A8 (HIGH, pool quality) The frontier tool reads a genome length that does not exist

`cli/pipeline.py:1010-1011`, inside the `--show-frontier` block:

```python
                fg_lengths = getattr(parameter, "fg_lengths", [1_000_000])
                bg_lengths = getattr(parameter, "bg_lengths", [])
```

`parameter.fg_lengths` does not exist. `fg_lengths` is a local inside
`get_params` at `parameter.py:1070`; the module global is `fg_seq_lengths`.
Confirmed after a real initialization against the *E. coli* configuration:

```
fg_lengths      False  None
fg_seq_lengths  True   [4641652]
bg_lengths      False  None
bg_seq_lengths  True   [3298430636]
```

So every `--show-frontier` run has computed coverage against a hardcoded
1,000,000 bp genome, whatever the target actually is, and specificity against no
background at all. On the *E. coli* design the real target is 4.6 times that
length, and the background it is being made specific against is hg38.

This is the class CLAUDE.md's Known Issue 8 called the last known instance of
"a config key or flag that is documented, accepted, and read by nothing". It is
a variant rather than a repeat: the flag is read, but a `getattr` default
silently substitutes for a global whose name is wrong. `additionalProperties`
is not the mechanism here and no schema check would have caught it.

**Proposed change (S).** Read `fg_seq_lengths` and `bg_seq_lengths`. Add a test
that asserts `parameter.fg_lengths` is absent, so the wrong name cannot come
back.

### A6 (LOW, pool quality) Near-duplicate primers reach the delivered panel

Delivered *E. coli* set 0, 160 primers drawn from a 449-candidate pool:

- 9 pairs at Hamming distance 1 or less
- 5 primers share the 3' hexamer GCGAAA, 4 share GCCAGA, 4 share GGATAA
- no reverse-complement pairs, which is the case already handled

The pool itself holds 65 pairs at Hamming distance 1 or less. Primers differing
by one base bind largely the same sites, so the second one buys little coverage
while adding synthesis cost and dimer surface. Shared 3' ends correlate
mispriming behaviour across the panel.

The two offered alternatives are also not comparable panels: set 0 has 160
primers and set 1 has 112, because alternatives are found by removing the
chosen primers and selecting again from what remains.

**Proposed change (S).** Add a redundancy penalty to the greedy: reject a
candidate whose site set is largely contained in the already-selected union.
The positions are already in the cache the optimizer holds.

---

### A7 (LOW, robustness) An empty coverage graph is announced as full coverage

Failure on an inadequate pool is otherwise sound. Given a five-primer pool and
a request for fifty, `optimize` exits 1, writes no summary, warns
`PARTIAL result: found 0 primers but target was 50` and offers four specific
remedies naming the parameters to change. That is the behaviour you want.

One line in the middle of it is not:

```
INFO: Graph: 0 primers, 0 regions of 1284 genome bins
INFO: Full coverage achieved
INFO: Selected 0 primers
INFO: Coverage: 0.0% (0/0 regions)
```

"Full coverage achieved" is vacuously true of zero regions and is printed
immediately before a zero. In this reproduction the empty graph came from
absent position files, which is E3's failure mode, and the run was saved only
by selecting nothing at all. A pool where some primers are indexed and some are
not would pass this point quietly.

## B. `filter`

### B1 (HIGH, speed) The background is scanned for candidates the step is about to discard

`pipeline.py:1103-1107` scans background positions for every candidate that
survived the frequency and thermodynamic gates. The call's return value is
discarded; it exists for the HDF5 side effect that step 4 reads.
`pipeline.py:1129` cuts the pool to `max_primer` only afterwards.

In the fresh run that meant scanning for 20,301 primers and keeping 3,000.

The cut does not need scanned positions. `_rank_by_occupancy` at
`pipeline.py:885` reaches `_occupancy_ratio_column` at `:937`, then
`occupancy.weighted_site_load` at `:147`, then
`mismatch_counts.mismatch_class_counts`, which reads the jellyfish count
tables. Neither `occupancy.py` nor `mismatch_counts.py` imports h5py. The Gini
gate is the only step-2 consumer of scanned positions and it uses the
foreground only.

So the ranking and the cut can move ahead of the background scan with nothing
lost and no approximation anywhere. The correct order is: frequency and quality
gates, foreground scan, Gini, occupancy ranking and cut, then background scan
on the cut pool.

Measured on chr21 with the genome already parsed:

| primers scanned | scan time | HDF5 written |
|---|---|---|
| 3,000 | 2.5 s | 2.2 MB |
| 20,301 | 7.9 s | 16.0 MB |

Aho-Corasick scaling measured on the same background, extrapolated to hg38:

| patterns | chr21 scan | chr21 matches | projected hg38 |
|---|---|---|---|
| 2,000 | 0.86 s | 4,975 | 1.0 min |
| 39,993 | 1.86 s | 113,184 | 2.2 min |
| 799,897 | 12.75 s | 2,139,894 | 15.0 min |

The shipped `tests/validation/genomes/filter_stats.json` records 369,459
primers reaching the scan and 2,000 surviving. Counting reverse complements
that is the bottom row against the top one, and it accounts for very nearly all
of the 924 s the project has recorded for a filter step against a large host.

**Proposed change (S).** Move the background `get_positions` call to after
`_rank_and_cut_candidates`. The foreground scan stays where it is.

**Not recommended: a count-based pre-cut ahead of the foreground scan.** It
looks like the same idea and it is lossy. The exact-count `ratio` key correlates
only weakly with the occupancy key it would stand in for, Spearman 0.721 on
*E. coli*, 0.340 on *S. aureus*, 0.703 on *M. tuberculosis*, and a pre-cut to
the top 25% discards between 5% and 39% of what the occupancy ranking rates
best. The key also degenerates exactly where this tool is used most: on
*Prevotella* against human chr21 all 369,431 survivors tie at `ratio` 0.0, so
the pre-cut silently becomes "the most abundant N".

### B2 (HIGH, robustness) One join holds the whole genome twice, and that is the 8.5 GB ceiling

`runs/gc_tiers/mid_ecoli/step2.log`:

```
Loading genome sequence from .../human_full.fna...
Cached genome sequence: 3,298,430,636 bp
done (52.1s)
...
Step 2 complete in 97.7s
```

The 52 s is not parsing. A Biopython pass over chr21 takes 0.19 s, which
projects to about 13 s for hg38. The cost is at `string_search.py:68`:

```python
sequence = "".join(loader.load_genome_streaming(seq_fname))
```

`str.join` calls `PySequence_Fast` on its argument, so a generator is
materialised into a list of every record before any concatenation begins. Every
chromosome string is therefore alive at the moment the finished sequence is
allocated. Measured with incompressible input, resident memory grows linearly
with the parts and does not fall as they are consumed. For hg38 that is roughly
3.3 GB of records plus a 3.3 GB result, which is where the recorded 8.5 GB peak
comes from and why hg38 filter runs have to be serialised.

The comment directly above that line states the opposite, that "the peak here
is the finished sequence plus the largest chromosome". That is what a streaming
join would do and is not what this line does. The comment is correct about the
defect it replaced, a character-at-a-time generator costing 26 GB, and wrong
about the one it introduced.

`string_search.py:30` also holds `_genome_cache` as a module-level dict, so the
cache is per-process and every CLI invocation repeats the whole thing. There is
no on-disk parsed form, and `clear_genome_cache` at `:91` has no caller outside
the tests, so the genome stays resident alongside the 1.24 GB k-mer count table
that `mismatch_counts.py:95` caches for the occupancy ranking.

**Proposed change (M).** Write the uppercased, header-stripped sequence once to
a uint8 sidecar beside the k-mer tables and memory-map it thereafter, decoding
one window at a time inside the loop already chunked at `MAX_SCAN_CHUNK`. A
hard constraint applies: `Automaton.iter()` raises `TypeError: string required`
for both `bytes` and `bytearray`, so the mmap cannot be handed to the scanner
directly, and decoding a window costs about 0.7 s per 3.3 GB. Two-bit packing
is not worth pursuing, since it cannot represent N and the scanner needs a
`str` regardless. Release the genome cache once the scans are done.

### B3 (MEDIUM, reporting) The filtering funnel names a stage that cannot filter and hides the stage that does most of the cutting

`pipeline.py:1044` sets `after_thermodynamic` and `pipeline.py:1111` sets
`after_background` from the same dataframe. The only things between them are the
exclusion and blacklist blocks at `:1055` and `:1071`, both guarded on
configuration that is absent in every run examined, so the two counts are equal
whenever no blacklist is configured, which is the normal case. Confirmed on four
independent runs: 458/458 (*E. coli*), 1231/1231 (*S. aureus*), 337/337
(*M. tuberculosis*), 20329/20329 (Wolbachia). The real background gate is the
`bg_bool` term at `pipeline.py:1029`, computed at `filter.py:648` and folded
into `after_frequency`. The report labels the inert stage "After
background/blacklist" at `report/metrics.py:202`, which is what makes it read
as a background result.

What the funnel should record instead, with *E. coli* values:

| stage | count |
|---|---|
| total_kmers | 2,542,354 |
| after_fg_frequency | not currently measured |
| after_bg_frequency | 1,639 |
| after_thermodynamic | 458 |
| after_exclusion_blacklist, omitted when unconfigured | 458 |
| after_gini | 449 |
| after_max_primer_cut | 449 |

Both booleans exist as columns before they are dropped at `:1030`, so splitting
the frequency row is two calls to `len`.

One further number belongs in the report. At k=12 the length scaling at
`filter.py:145` divides `max_bg_freq` by 16, so a configured 2.425e-07 is
applied as 1.5156e-08, which against hg38 permits 50 sites rather than the 800
the configured value implies. `describe_freq_gate` exists precisely to make this
visible; the point is that the funnel should carry it.

Meanwhile the `max_primer` cut appears only as `final_candidates`. In the fresh
run it dropped 20,301 to 3,000, 85% of survivors, under no stage name at all.
The Gini gate, which does have a name, removed 9, 9, 18 and 28.

A reader of the quality report is therefore told that background filtering did
nothing, which is false, and is not told that a size cap removed most of the
pool, which is true and is the stage most worth tuning.

**Proposed change (S).** Record the background contribution separately at
`:1029`, and add an explicit `after_max_primer_cut` stage. Consider whether
`max_gini` earns its default place given what it removes.

### B4 (HIGH, pool quality) The Gini gate rewards the primers it exists to exclude

`primer_attributes.py:175` returns NaN when a primer has no sites at all, but a
primer with exactly one site passes that guard and scores a Gini of 0.0, the
best value available. A measure of how evenly a primer's sites are spread is
meaningless for a primer with one site, and as written it ranks that primer
first.

How much this bites depends on the pool, and the two regimes are far apart. On
the shipped *Prevotella*-against-chr21 pool, 86% of the 10,000 kept primers
score exactly 0.0 and every one of those has two or fewer foreground sites,
against a pool mean of 1.29. On the plasmid example the figure is 96.2% of 500
rows, again all at two sites or fewer. On the three whole-genome GC-tier pools,
which reach step 3 already cut to a few hundred abundant primers, no primer
scores 0.0 at all and mean foreground counts run 5.7 to 26.4.

So this is a defect of small targets and of large sparse pools, which is the
regime a distant background and a generous `max_primer` produce, and it is
invisible in the whole-genome runs this project usually inspects. Any change to
the gate has to hold in both regimes.

The gate also does not agree with the optimizer about what is good. Median Gini
of the delivered set against the rest of the pool: 0.429 against 0.435 on
*E. coli*, 0.310 against 0.384 on *S. aureus*, 0.460 against 0.485 on
*M. tuberculosis*. And the threshold is either inert or harmful:

| `max_gini` | *E. coli* removed from pool / from delivered set | *S. aureus* | *M. tuberculosis* |
|---|---|---|---|
| 0.7, the configured value | 0 of 449 / 0 of 160 | 0 of 1222 / 0 of 200 | 0 of 319 / 0 of 36 |
| 0.6, the documented default | 29 / 9 of 160 | 6 / 0 of 200 | 29 / 5 of 36 |
| 0.5 | 110 / 41 of 160 | 306 / 19 of 200 | 132 / 13 of 36 |

At the value these runs used it removes nothing that matters. At the documented
default it deletes primers the optimizer went on to choose.

**Proposed change (S for the code).** Return NaN below a minimum site count,
which the existing `.notna()` guard at `filter.py:753` already rejects. Then
re-derive the threshold against delivered coverage rather than against a
distribution the single-site primers dominate.

### B5 (HIGH, pool quality) The candidate loader filters on a melting temperature it knows to be wrong

`melting_temp.py:153`:

```python
fgc = 1.0 / len(s)  # replicates original melt bug for compatibility
```

This shim is used at `kmer_counter.py:565` with a fixed 15 C margin, to decide
which k-mers become candidates at all. The compatibility it preserves is with
the random forest model retired from the default path on 2026-09-05, so the
reason for the bug outlived the bug's purpose.

The systematic offset is 10.0 C at k=12 and 20.6 C at k=6 with no additives.
Candidates that pass the real Tm window and are discarded before ever reaching
it:

| conditions | k=8 | k=12 |
|---|---|---|
| plain phi29 | 0.0% | 0.0% |
| DMSO 10% | 0.3% | 1.8% |
| DMSO 10% and betaine 1.5 M | 0.2% | 9.6% |

So the loss is nil on the simplest chemistry and grows with exactly the
additives the tool exists to model. It is silent in every case.

**Proposed change (S).** Pre-filter with the same
`ReactionConditions.calculate_effective_tm` that `filter_extra` uses.

The same call is also the loader's dominant cost: 10.5 us per k-mer, which is
21 s over a 2.03 M-kmer *Prevotella* pool and 27 s over a 2.54 M-kmer *E. coli*
one. `melting_temp._overcount` at `:40` runs sixteen separate `str.find` scans
per k-mer where one pass over dinucleotides would do.

### B6 (MEDIUM, pool quality) Cross-primer dimers are not screened at all by default, and the screen that exists rejects most of the pool

`filter.py:562` screens self-dimers only, each primer against itself. No default
path screens a primer against another primer. That is the gap the clique and
network optimizers exist to work around, and it is why A2's delivered panels
carry 10 and 11 bp heterodimers.

The `--enable-qa` screen that would fill the gap is not usable as it stands.
Measured rejection on the three shipped pools: 297 of 449 on *E. coli* (66.1%),
897 of 1222 on *S. aureus* (73.4%), 188 of 319 on *M. tuberculosis* (58.9%),
almost entirely on 3' stability. `pipeline_qa_integration.py:128` calls
`create_three_prime_analyzer`, while `three_prime_stability.py:621` offers
`create_three_prime_analyzer_adaptive` and describes it as the one to use below
35% or above 65% GC. *S. aureus* at 33% and *M. tuberculosis* at 65% are
precisely those cases, and *S. aureus* has the highest rejection rate.

**Proposed change (S for the factory, L for calibration).** Switch to the
adaptive factory and re-measure. Do not make this default-on yet: a screen that
rejects two thirds of a delivered pool is either finding something systemic or
is miscalibrated, and the factory mismatch makes the second more likely.

### B7 (MEDIUM, speed) "Reusing existing position files" is printed and the scan re-runs in full

`pipeline.py:1094` computes `position_files_exist` and feeds it to nothing but a
log line. On the Aho-Corasick branch the automaton is rebuilt from the whole
primer list every time; `check_which_primers_absent_in_h5py` is reached only
from the fallback path at `string_search.py:445`. Re-running `filter` after a
threshold tweak re-pays B1 in full while announcing that it did not.

**Proposed change (S).** Intersect the primer list with the existing HDF5 keys
before building the automaton.

### B8 (LOW, speed) Melting temperature is computed about twice per candidate

`runs/gc_tiers/mid_ecoli/step2.log`:

```
[Step 2] Thermo cache 'enthalpy_entropy': 52.3% hit rate
(1,796 hits, 1,639 misses, 1,639/1,000,000 entries)
```

and in the fresh run, 81,204 hits against 31,578 misses. The miss count equals
the surviving-candidate count exactly in both, so the cache is working as
designed and the primers are simply evaluated more than once per step. A cache
sized at a million entries holding 1,639 also suggests the sizing was never
matched to the workload.

---

### B9 (LOW, robustness) A circular genome is held twice during the scan

`string_search.py:157`:

```python
    if circular:
        search_seq = sequence + sequence[: max_k - 1]
```

The concatenation allocates a full second copy of the genome, and the original
stays alive in `_genome_cache`, so a circular scan peaks at twice the target's
size. This is not on the hg38 path, because a host background is linear and
`bg_circular` is False in every shipped configuration. It is on the foreground
path for every bacterial design, where `fg_circular` is normally True, and a
foreground is megabases rather than gigabytes, so the cost is real and small.

The wrap only needs `max_k - 1` bases. Scanning the body and then a short
join-region window separately would remove the copy. Recorded here rather than
planned: it is the same shape as B2 but two orders of magnitude smaller, and it
should be fixed alongside B2 if that work touches this function anyway.

## C. `count-kmers`

### C1 (HIGH, robustness) Jellyfish output is reused on filename alone

`kmer_counter.py:336` and `:189` both guard the call with
`if not os.path.exists(txt_file):`. Nothing compares the FASTA against what
produced the existing table. The manifest writes `input_checksums` but no step
reads them back.

The failure is silent and complete: point `fg_genomes` at a new assembly while
leaving `fg_prefixes` unchanged, and the whole design is built from the
previous organism's k-mer counts with no warning. Step 1 reporting "complete in
2.0s" with hg38 as the background is this path firing.

**Proposed change (S).** Store the input checksum beside the table and compare
on reuse. Skip when it matches, recount and say so when it does not.

---

## D. `score`

### D0 (HIGH, regression) Retiring the amplification model gave every primer the same quality rating again

This is the defect commit `a2a58c0` fixed on 2026-09-04, reintroduced by the
retirement commit `c92888d` the following day, from the other end.

`report/metrics.py:761` backfills `amp_pred` from `step3_df.csv`. The default
path no longer writes that column, so the backfill finds nothing and every
primer keeps 0.0. The rendering layer then substitutes a constant:

```python
# executive_summary.py:503
quality_score = primer.amp_pred if primer.amp_pred > 0 else 0.5
stars = min(5, max(1, int(quality_score * 5 + 0.5)))
```

`int(0.5 * 5 + 0.5)` is 3, so every row renders three stars. Measured on the
shipped post-retirement run at `runs/gc_tiers/mid_ecoli`: 272 primers, one
distinct `amp_pred` value, `{0.0}`, and an identical rating on every line.
`technical_report.py:186` loses its 0.2 amplification term uniformly, and
`visualizations.py:382` and `:537` colour every point 0.5.

`tests/test_score_stage_is_retired.py` asserts `amp_pred == 0` under a docstring
promising that "the report must omit the rating rather than reinstate that
constant". The metrics layer complies. The rendering layer does not, and no test
looks at it.

**Proposed change (S).** Drop the star column and the amplification term when no
primer carries a score, rather than falling back to a midpoint.

### D1 (MEDIUM, speed) The step's own work is a tenth of the cost of invoking it

Measured on the fresh run: 1.71 s wall for a step that reports 0.2 s. Broken
down on a separate invocation:

| phase | cost |
|---|---|
| `import parameter` | 0.38 s |
| `import pipeline` | 1.89 s, of which sklearn is 1.40 s |
| `_initialize()` | 0.30 s |
| the step body | 0.01 s |

So of a 2.6 s invocation, 0.01 s is the work. The 1.40 s is scikit-learn
imported for a model this path never loads, which is E6. The 0.30 s is
`_apply_gc_adaptive_defaults` needing the foreground GC, which makes
`parameter.get_params` read the foreground FASTA at `parameter.py:1354`. That
read scales with target size, so it is 0.30 s on a 4.6 Mb genome and
proportionally worse on a large one. The 3.3 Gb background is not read; only
`bg_seq_lengths` is used, which is what raises the misleading "Background genome
is large" warning in a step that never touches it.

This is not an argument for deleting the step. Twenty or more documentation
files and the `neoswga-cli` skill teach it as step 3, six modules read
`step3_df.csv`, `neoswga start --start-step/--stop-step` addresses it by number,
`tests/conftest.py:44` primes the example directory by running it, it writes its
own manifest entry with input checksums, and it is the interposition point where
a user can hand-edit `step2_df.csv` and re-derive the pool without repeating a
long filter. Folding it in saves 0.3 s and breaks a published interface.

**Proposed change (S).** Make `rf_preprocessing` a lazy import inside
`_score_with_amp_model`, and let the default path skip the GC-adaptive
derivation it never applies. One caveat: the manifest records
`effective_conditions` per step, and those conditions come from that derivation,
so skipping it makes the `score` entry differ from the others. Recording
conditions a step never applied is arguably worse than omitting them, but that
is a judgement to make deliberately.

### D1b (MEDIUM, speed) `--full-score` costs 126 times as much for no change

Measured on the 449-candidate *E. coli* pool: 767.6 s against 6.1 s for
`--amp-model` alone, for a mean absolute change of 0.0016 on an 11.1 to 19.1
scale, Pearson 1.0000, Spearman 0.9999, and an identical delivered order.
Nothing in `tests/` exercises `create_augmented_df(skip_delta_g=False)`. The
flag's own help text claims the delta-G features are "under 2% of model accuracy
but over 99% of scoring compute time"; the measurement puts the accuracy figure
under 0.02%.

**Proposed change (S).** Delete `--full-score`. Keep `--amp-model`, which works
correctly and is currently the only thing that makes D0's quality column vary.

### D1c (MEDIUM, pool quality) Step 3 discards step 2's ranking and hands the optimizer the worse end of the pool

Step 2 spends up to 50 s ranking candidates by occupancy-weighted background
load (`_rank_by_occupancy`, `pipeline.py:837`), a ranking its own docstring
credits with moving coverage from 7.6% to 40.3% and selectivity from 0.69 to
2.54. It writes `step2_df.csv` in that order. `order_step3_rows` at
`pipeline.py:607` then re-sorts by Gini and throws it away.

The two orders are Spearman -0.185 on the *E. coli* pool and share 2 of their
top 24. The head of the pool the optimizer actually receives is measurably the
worse end:

| top 24 by | mean fg sites | mean bg sites | mean occupancy ratio |
|---|---|---|---|
| step 2's own rank | 35.7 | 5.0 | 3.74 |
| Gini, as shipped | 10.5 | 29.7 | 40.34 |

This matters because `dominating_set_optimizer.py:49` documents its tie-break as
deliberately respecting the caller's ranking, on the grounds that "the caller's
list arrives ranked by the scoring step". That ranking no longer arrives.

**Honest caveat.** No effect on the delivered set could be detected. Four
orderings of the *E. coli* pool, Gini, step 2's rank, Gini reversed and a seeded
shuffle, gave Jaccard 1.000 and identical coverage under both `dominating-set`
and `hybrid` at target sizes 6, 12 and 24. The evidence pinned in
`tests/test_step3_order_is_deterministic.py` does not reproduce on the run its
docstring names. So this is worth fixing because it is free and because the
inversion is indefensible on its face, not on a promise of changed output.

**Proposed change (S).** Sort by step 2's rank, with Gini demoted to a
tie-break and the primer sequence last for totality.

### D1d (MEDIUM, robustness) `score --enable-qa` reintroduces the unstable sort the ordering work removed

`pipeline_qa_integration.py:718` sorts on `composite_score` alone with the
default quicksort and no secondary key. Measured on the 449-candidate *E. coli*
pool: two input orders give 65 differing positions, the first at rank 14, where
two primers swap. 100 of the 449 rows share a composite value with another row.

Underneath it, 55% of the composite's declared weight is a constant.
`integrated_quality_scorer.py:211` declares dimer at 0.35 and strand bias at
0.20, but `score_primer` sets both to 1.0 for every primer at `:269` and `:251`,
because `rank_by_quality` calls it with no binding sites. Measured over the 449
candidates, `dimer` and `strand` each take exactly one distinct value; only 3'
stability, thermodynamics and complexity vary, and the composite spans 0.8505 to
0.9644.

So the QA ordering is an 11% band, 45% of whose declared weight cannot move,
sorted unstably. It is not a better default than Gini and it is not a better
default than step 2's rank.

### D2 The ordering hypothesis, tested and largely refuted

Since `order_step3_rows` sorts by Gini and then by primer sequence, and since
near-identical sequences have near-identical Gini, an order-sensitive greedy
might be expected to take runs of near-duplicates. Measured on the three
delivered pools, the effect exists and is small.

| pool | Hamming-1 pairs | median gap | null median | pairs within gap 5 | expected | adjacent | expected |
|---|---|---|---|---|---|---|---|
| *E. coli*, n=449 | 65 | 84 | 130 | 6 | 1.5 | 2 | 0.3 |
| *S. aureus*, n=1215 | 156 | 274 | 352 | 2 | 0.9 | 1 | 0.2 |
| *M. tuberculosis*, n=319 | 185 | 81 | 92 | 3 | 6.3 | 0 | 1.8 |

Near-duplicates sit closer together than chance would put them, by 12% to 35%
on the median gap. But the near-adjacency excess is tiny in absolute terms, and
on *M. tuberculosis* it reverses: three pairs within gap 5 against 6.3 expected,
and no adjacent pairs at all against 1.8 expected. Nothing here is enough to let
a greedy take several near-duplicates in a row.

The related claim that Gini is so often exactly 0.0 that the order is
effectively alphabetical holds only on the plasmid example, where 96.2% of the
500-row pool sits at 0.0 with nine distinct values and 99.2% of adjacent rows
tied. All three whole-genome pools have no rows at 0.0 at all, and 384, 910 and
319 distinct Gini values. So a rule justified on the plasmid would be tuned on
the atypical case.

The ordering is therefore not the cause of A6. The absence of any similarity or
dimer rejection test in the selection rule is.

---

## E. Between the steps

### E1 (HIGH, robustness) The manifest cannot answer "what produced this file"

Each step entry carries exactly `cli_invocation`, `effective_conditions`,
`git_sha`, `input_checksums`, `jellyfish_version`, `neoswga_version`,
`params_path`, `platform`, `python_version`, `resolved_params`, `seed`, `step`
and `timestamp_utc`. There is no timing field of any kind, so the only record
of where a run spent its time is a line in a step log.

`runs/gc_tiers/mid_ecoli/run_manifest.json` holds 19 appended entries spanning
three git SHAs, and nothing marks which of the several `optimize` entries
produced the `step4_improved_df.csv` sitting beside it.

The last of those entries has a `cli_invocation` ending `-n 160` and produced a
160-primer panel, while its `resolved_params` still records `num_primers: 24`.
Command-line overrides reach the invocation string and not the resolved
parameters, which is the same class of drift the file already documents for
`effective_conditions` against `resolved_params`.

**Proposed change (S).** Record wall time per step. Record the effective set
size the way `effective_conditions` records the effective chemistry. Record the
output checksum so an entry can be matched to the file it wrote.

Note that `run_manifest.py:87` already accepts an `extra` field, merged into the
entry at `:151`, and every step handler already measures its own elapsed time
into a local that is then discarded (`cli/pipeline.py:141`, `:359`, `:454`,
`:1169`). The timing finding is a one-line change at four call sites.

`read_effective_conditions` (`run_manifest.py:198`) returns the last entry
carrying conditions whatever step it belongs to. In `mid_ecoli` two `score`
entries sit after an `optimize` entry, so a report generated at that point would
describe the optimize result under the score step's reaction. They agree today
only because both derive from params.json, which stops being true under E4.

### E2 (HIGH, pool quality) The GC-adaptive strategy overrides an explicitly configured additive concentration

`pipeline.py:559`:

```python
current_betaine = getattr(parameter, "betaine_m", 0.0)
if current_betaine == 0.0 and adaptive_params.betaine_concentration > 0:
    parameter.betaine_m = adaptive_params.betaine_concentration
```

The guard cannot distinguish a user's explicit `"betaine_m": 0.0` from the
unset default. The DMSO branch below it is identical. The k-mer branch twenty
lines above was fixed for exactly this and consults the raw JSON at `:534`,
which is why one run logs "Preserving user-specified k-mer range" and the next
line logs "Setting betaine to 1.0M".

This is visible in the shipped runs. `runs/gc_tiers/mid_ecoli/params.json` sets
`"betaine_m": 0.0`, and all eighteen manifest entries record
`effective_conditions.betaine_m = 1.0`. Every melting temperature in that design
was computed against a buffer additive the user had excluded.

It is not confined to GC-extreme genomes. `get_params` computes `genome_gc` from
the FASTA when params.json omits it (`parameter.py:1341-1360`), so the early
return at `pipeline.py:471` is never taken and the branch runs on every design.

**Proposed change (S).** Test `"betaine_m" in parameter._json_data`, as the
k-mer branch does, for both additives.

### E3 (HIGH, robustness) Two of the four prerequisite validators are never called, and step 4 is one of them

`validate_step1_prerequisites` at `pipeline.py:147` and
`validate_step4_prerequisites` at `:260` are defined and have no call site.
Only steps 2 and 3 check, at `:983` and `:1330`.

Step 4 therefore builds its position cache with the default `on_missing="warn"`,
and `position_cache.py:163` warns that the affected primers "will score as zero
coverage" and that "the coverage number is meaningless, not low". The run then
selects a set and reports that coverage. The validator written to catch this,
with the per-primer-length file list already assembled, sits a hundred lines
above the one that is used.

**Proposed change (S).** Call it and raise `StepPrerequisiteError(4, ...)`.

**Verify.** Delete the `_positions.h5` files and re-run `optimize`; expect a
named failure rather than a panel scored against nothing.

### E4 (MEDIUM, pool quality) Chemistry flags exist on `filter` and on no other step

The `filter` subparser registers sixteen chemistry flags plus `--preset`
(`cli/pipeline.py:1364-1456`). `_optimize_parser.py` registers none; `optimize`
inherits only `--polymerase`. The optimizer builds its conditions from the
parameter globals at `unified_optimizer.py:916`, so
`filter --preset high_gc_genome` followed by a plain `optimize` filters under
one reaction and scores under another.

This is a flag asymmetry, not an adaptation problem. The GC-adaptive strategy
itself is not a source of step-2-against-step-4 drift, because
`unified_optimizer.py:891` calls `_initialize()` before building conditions, so
the same adaptation runs in both steps.

**Proposed change (S).** Compare the `effective_conditions` recorded for the
filter step against the one the optimize step is about to use, and warn when
they differ.

### E5 (MEDIUM, robustness) There is no memory guard anywhere on the pipeline path

A search across the package for `virtual_memory`, `psutil`, `MemoryError`,
`getrlimit` and `estimated_memory` returns two lines, both about GPU memory.
`MemoryLimitError` is defined at `exceptions.py:333` and raised nowhere. The
only large-input signal is advisory: `pipeline.py:412` sets
`_BLOOM_THRESHOLD = 50_000_000`, logs a suggestion and proceeds. A filter run
that exceeds available memory gets an OOM kill, which is how B2's ceiling
presents to a user.

More broadly, `exceptions.py` defines 29 classes over 413 lines and only two of
them are ever raised. `PositionFileNotFoundError`, `NoCandidatesError`,
`JellyfishError`, `MemoryLimitError` and `PipelineStateError` all describe
failures this pipeline really has; the real surface is `RuntimeError` and
`ValueError` caught by broad handlers at `cli_unified.py:391-423`.

### E6 (MEDIUM, speed) The CLI imports scikit-learn on every invocation to serve a retired model

`cli_unified.py:60` imports `core.pipeline` at module scope solely to make one
exception class catchable. `core/pipeline.py:16` imports `rf_preprocessing`,
which imports `sklearn` and `sklearn.ensemble` at module scope.

Measured with `-X importtime` on a cold run: 4.67 s total, of which
`rf_preprocessing` is 3.16 s and `sklearn` 2.47 s, against a bare interpreter of
0.03 s. Warm, the figure is the 2.25 s that `neoswga --help` costs.

The amplification model was retired from the default path on 2026-09-05, so the
four-step workflow pays this four times for scoring it does not perform, and
`--help` and `show-presets` pay it too.

**Proposed change (S).** Define `StepPrerequisiteError` in `core/exceptions.py`
and import it from there; move `import sklearn` inside the functions that use
it.

Relatedly, `design` already calls all four handlers in one process
(`cli/commands.py:605`) and so pays the import once. The four-step split is the
right unit of caching, because `filter` is expensive and a user wants to re-run
`optimize` without repeating it, and the wrong unit of invocation. `design`
lacks `--optimization-method` and the chemistry flags, which is the gap to close
before recommending it as the default entry point.

### E7 (LOW, usability) The tool is verbose about absent parameters and silent about misspelled ones

`parameter.py:820` logs a debug line for each of the twelve names in
`OPTIONAL_PARAMS` that a params.json omits. Seven are absent from the shipped
`mid_ecoli` file, which is the seven lines at the head of every step log, and
one of them announces a default for `min_amp_pred`, a gate that no longer runs
without `--amp-model`.

The other direction has no code at all. The schema declares 79 properties with
`"additionalProperties": true` (`schema/params.schema.json:7`) and
`param_validator.py` contains no unknown-key or close-match logic, so
`max_bg_freqency` is accepted in silence and the default applies, changing the
design.

**Proposed change (S).** Drop the debug line where the default is documented and
static, and add an unknown-key pass using `difflib.get_close_matches`. Leave
`additionalProperties: true` so it stays a warning.

### E8 (LOW, robustness) Two parallel implementations, one unreachable and one divergent

`improved_pipeline.py` is unreachable from any pipeline command; its only
non-test importer is `core/validation.py:441`, which constructs it with
`optimization_method="greedy"`, a retired method. Three CLI flags are documented
as unimplemented precisely because they reach only that module.

`multi_genome_pipeline.py` is reachable, from `cli/simulate.py:23` and
`cli/commands.py:549`, and uses none of `pipeline.step1` through `step4`. It
counts k-mers in a triple-nested Python loop at `:305` without canonicalisation,
where the main path passes `-C`. It has no Gini filter, no occupancy ranking, no
`max_primer` cut, no HDF5 positions, no step CSVs and no manifest. Only
`HybridOptimizer` is shared, so the two paths will keep diverging.

### E9 (LOW, robustness) The manifest read-modify-write is unlocked

`run_manifest.py:153-167` reads, appends and rewrites with no lock, so two runs
sharing a `data_dir` lose entries. `fcntl.flock` is already used correctly at
`experimental_tracker.py:243`.

### E10 (LOW, speed) There is no fast way to check a config before committing to a long run

`validate --quick` (`core/validation.py:523`) runs three synthetic tests and
opens no genome, no params.json and no pipeline step. `validate params -j` runs
the schema validator, which by E7 cannot catch an unknown key. So there is no
seconds-long check of a configuration before a thirty-second filter, or a
sixteen-minute one against hg38, although `tests/integration/` already has the
fixtures for it.

---

## What was checked and found sound

These were investigated and are not defects. They are recorded so the next
audit does not spend time on them.

- **`--enable-qa` does not carry between steps.** It is assigned, not merely
  set, at `cli/pipeline.py:75`, `:173`, `:387` and `:540`, each with a comment
  naming the in-process leak it prevents.
- **Two params files in one process do not contaminate each other.**
  `_initialize` resets on a changed JSON path (`pipeline.py:377-379`), which
  also drops the reaction-conditions singleton, and every reaction global is
  reassigned rather than merged (`parameter.py:1191-1193`). The narrow
  remaining gap is an in-process caller that edits params.json without changing
  its path.
- **Progress reporting is honest.** `progress_context` prints a real measured
  elapsed time (`progress.py:227-239`). The `ProgressBar` class with the
  estimating ETA has no caller on the pipeline path.
- **Step 2 and step 3 fail with actionable messages.**
  `StepPrerequisiteError` formats the missing files plus a remediation command
  (`pipeline.py:126-144`), and `check_genome_inputs` pre-flights genomes so a
  truncated download fails at step 1 rather than two steps later.
- **Exception swallowing is not a systemic problem.** Five `except: pass` sites
  across the package, none on the four-step main path in a position that could
  hide a step failure.
- **HDF5 writing is fine.** The file is opened once per prefix and k, not once
  per primer (`string_search.py:246`). Absence of chunking and compression is
  correct for datasets this small: 54 us and 400 bytes per dataset measured.
- **Reverse-complement duplicates do not reach the delivered panel.** Zero in
  the *E. coli* set.
- **The default k range is not a blanket 6 to 12.** `parameter.py:173` and
  `:1036` derive it from the polymerase preset, so 6-12 applies to Phi29, 10-18
  to EquiPhi29 and 15-25 to Bst. The seven-table sweep a blanket default would
  imply does not happen. Even for Phi29 the short end contributes little,
  since a balanced 6-mer melts at 4.0 C against a 20 C floor, but a k=6 table
  has 4096 rows and costs nothing to hold.
- **The thermodynamic cache is behaving as designed.** Its miss count equals
  the candidate count exactly, and the hits are accounted for by four
  evaluations per Gini survivor in the occupancy ranking. The only thing worth
  changing is the one-million-entry provisioning at `thermodynamics.py:182`.

---

## Three things worth saying plainly

**The delivered panel is chosen on foreground coverage and nothing else.** Every
other property in the summary, selectivity, dimer safety, evenness, Tm spread,
is measured after the fact and reported as though it had been optimised. That
is why panels ship with heterodimers three to four times the configured limit,
and why near-duplicate primers occupy 6% of a delivered set.

**The most expensive computation in the pipeline does not affect the answer.**
The all-pairs thermodynamic heterodimer screen costs more than everything else
combined, tests a criterion the user did not set, and in the measured cases
removed sixteen hairpin primers and left the delivered set identical to what the
method without it produces.

**Three silent-wrong-answer paths exist and none of them fails the run.** Stale
k-mer tables reused on filename alone (C1). Missing position files scored as
zero coverage (E3). An explicitly configured additive silently overridden (E2).
Each produces a complete, plausible, wrong result.

**One finding is a regression from this week.** Retiring the amplification model
on 2026-09-05 reintroduced, from the other end, the constant per-primer quality
rating that had been fixed the day before. Every primer in the shipped report
now shows three stars (D0).

## Recommendations, in order of payoff

0. **Stop the report giving every primer three stars** (D0). This is a
   regression from `c92888d` on 2026-09-05 and it is in front of users now. The
   test that was meant to prevent it checks the metrics layer and not the
   rendering layer.
1. **Cache the pairwise dimer verdict and put a cheap test in front of it**
   (A1). One change removes the largest cost in the pipeline, turns up to nine
   quadratic passes into one, and drops 2 GB of transient allocation. Nothing
   about the delivered set needs to change.
2. **Make `max_dimer_bp` the one dimer criterion and let it reach selection**
   (A2, B6). Two of the three delivered panels examined carry heterodimers of
   10 and 11 bp against a configured 3, because no default path screens a
   primer against another primer and the method that selects them contains no
   dimer test at all.
3. **Move the background scan behind the `max_primer` cut** (B1). A few lines.
   On the shipped validation configuration it is the difference between 15
   minutes and 1 minute, and it accounts for very nearly all of the 924 s
   recorded for a filter against a large host.
4. **Close the three silent-wrong-answer paths** (C1, E3, E2). Compare the
   recorded input checksum before reusing a k-mer table. Call
   `validate_step4_prerequisites`, which is already written. Test key presence
   rather than value equality before overriding an additive. All three are
   small and each one can otherwise deliver a confident wrong design.
5. **Fix the genome join** (B2). One line holds hg38 twice, which is the 8.5 GB
   peak that forces filter runs to be serialised. The sidecar and mmap are the
   larger follow-on.
6. **Make coverage incremental in background pruning** (A3). Arithmetic
   identical, and it is the difference between a clinical run finishing and a
   clinical run being abandoned.
7. **Move `sklearn` off the import path** (E6). Half the CLI start-up, paid four
   times per workflow, to serve a model the default path retired.
8. **Fix the Gini gate's degenerate case, then re-derive its threshold** (B4).
   As configured it removes nothing; at its documented default it removes
   primers the optimizer chose.
9. **Give the candidate loader the real melting temperature** (B5). The
   deliberate bug it preserves was for a model that no longer runs, and the loss
   grows with exactly the additives this tool exists to model.
10. **Make the funnel and the manifest tell the truth** (B3, E1). Neither
    changes a result. Both change whether anyone can tell what happened, and
    the manifest already has the field and the measurement.
11. **Fix the frontier's genome length, then surface the sizing tools** (A8,
    A5). `--show-frontier` has been scoring every design against a hardcoded
    1 Mb genome and no background, so fix that before documenting it. Then lift
    its 20-primer ceiling, let `design` reach it, and print marginal coverage at
    the delivered size. The knee is real and this project's own sweep was done
    by hand.
12. **Decide about `normalized_score`** (A4). Either say it describes the result
    or make it select it. The present arrangement reports an objective the tool
    does not pursue.
13. **Stop step 3 inverting step 2's ranking, and delete `--full-score`**
    (D1c, D1b). The first hands the optimizer the worse end of the pool by an
    eleven-fold margin on occupancy ratio, though no change in the delivered set
    could be measured. The second costs 767 s against 6 s to move a score by
    0.0016.

Nothing in this list has been implemented. The audit changed no code.
