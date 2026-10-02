# A single-reference design loses most of its coverage on another supergroup, and the host it is given decides its specificity against a host 1,000 times larger

Phase 3 of `docs/superpowers/plans/2026-10-02-genomic-diversity-and-multi-host.md`.
This record reports measurements, answers Q1 to Q4, and reads decision gate A
against them: the target-diversity half is decided, the host half is not. Q5
and the full D3 are not answered, for the reason given below.

## The instance

One chemistry, one optimizer seed, these five strains, two design hosts and one
evaluation-only host. Every claim below is bounded by that.

| | |
|---|---|
| Targets | *Wolbachia* wMel (`wolbachia`, the reference), wRi, wPip, wAlbB, wBm |
| Supergroups | A: wMel, wRi. B: wPip, wAlbB. D: wBm |
| Design hosts | *L. plantarum* WCFS1 (3,348,624 bp, 4 records), *Drosophila* (143,726,002 bp, 1,870 records) |
| Evaluation-only host | hg38 (3,298,430,636 bp, 705 records), on the counts route. **No design has had it as a background** |
| Host size span | 3.3 Mb to 3.30 Gb, a factor of 985 |
| k | 12 (`min_k` = `max_k` = 12) |
| Chemistry | phi29, 30.0 C, Na 50 mM, Mg 10 mM, no additives |
| Panel sizes requested | 6, 12, 24 |
| `optimize` seed | 0, one run per design and size |
| C2 control seeds | 20 per design and size |
| Optimizer | `hybrid`, `max_dimer_bp` 3, no panel limit set |
| Coverage reach | 3,000 bp (phi29, realistic) |
| Designs | D1 (reference only, one host) x2; D2 (five strains pooled, one host) x2; C1 (leave one strain out, vs *Drosophila*) x5; D3 preliminary (reference, both small hosts) |
| git sha of the run | `66c4c19c92d2d3c3d170eac4baf8caf9fff6756c` |
| Script | `scripts/benchmarking/diversity_baseline.py` |
| Results | `tests/validation/genomes/diversity_baseline/results.json` (gitignored) |

Ten designs were built from `count-kmers` onward, each delivered a set at all
three requested sizes, and `delivered_size` equalled `requested_size`
everywhere. No design failed.

## Two densities, and which is which

Two different quantities are called a selectivity density in this repository
and they are not interchangeable. Every density column below says which one it
is.

- **exact-site density**: the `target_host_pairs` figure of
  `reference_panel_evaluation`, formed from exact k-mer site counts, each
  divided by its own reference length.
- **occupancy-weighted density**: derived by the script from the two
  `weighted_site_load` figures that the same function returned, through
  `selectivity.selectivity_density_from_loads`. This is the quantity the
  optimizer reports as `selectivity_density` in occupancy mode.

They differ by more than a factor of ten on these panels and must never be
read as one series.

**The exact-site ceiling.** A control panel that binds a host nowhere has no
measured ratio: the package substitutes `MAX_SELECTIVITY` (1e6). A median over
a set of seeds containing one of those is not a median of measured ratios, so
where that happened the table reports `ceiling in N/20` and no spread. The
count of such seeds is the figure, not a number formed partly from the ceiling.

## Host coverage is not tabulated across routes

A host read from a position index and a host scanned from its FASTA do not
give comparable coverage: in a scanned multi-record reference the reach windows
cross record joins, which an index built per record does not do. The previous
measurement on one panel was 0.010683 against 0.010931 for the same quantity by
the two routes. Host coverage is therefore not tabulated here, and host
specificity is reported as sites, site density and the two pair densities,
which are route-independent for exact counts.

## Q1. How far does coverage fall from the reference to the other strains under D1?

D1, reference wMel only in the foreground, *Drosophila* as the host, set 0.
Coverage per strain, with the C2 median and range over 20 random size-matched
panels from the same candidate pool.

| n | wMel (ref) | wRi (A) | wPip (B) | wAlbB (B) | wBm (D) |
|---|---|---|---|---|---|
| 6 | 0.6569 | 0.5872 | 0.2690 | 0.2832 | 0.2120 |
| 12 | 0.7381 | 0.6678 | 0.3120 | 0.3131 | 0.2515 |
| 24 | 0.7792 | 0.7107 | 0.3238 | 0.3227 | 0.2931 |
| C2 median at n=12 | 0.330 | 0.298 | 0.054 | 0.055 | 0.058 |
| C2 range at n=12 | [0.233, 0.474] | [0.207, 0.425] | [0.020, 0.154] | [0.020, 0.149] | [0.017, 0.186] |

Read as a fraction of the reference's own coverage at n=12: wRi 90 percent,
wPip 42 percent, wAlbB 42 percent, wBm 34 percent. The ordering follows the
supergroup structure recorded in Phase 1 (A, then B, then D) at all three
sizes, and the two supergroup B strains are indistinguishable from each other
to within 0.001.

The control matters here. The designed panel is above the C2 median on every
strain, including the ones it was not designed against (0.3131 against a C2
median of 0.055 on wAlbB), and outside the C2 range on all five. So the
retained off-strain coverage is a property of the selected panel and not of any
panel drawn from this pool. What falls with divergence is the amount retained,
not whether anything is retained.

Size does not recover it: going from 6 to 24 primers adds 0.12 on the reference
and 0.08 on wBm, leaving the gap roughly where it was.

## Q2. Under D2, what is the worst strain, and is it the smallest or the most divergent?

D2, all five strains pooled in the foreground, *Drosophila* as the host.

| n | wMel | wRi | wPip | wAlbB | wBm | worst | best |
|---|---|---|---|---|---|---|---|
| 6 | 0.5514 | 0.4697 | 0.3057 | 0.6434 | 0.2375 | wBm 0.2375 | wAlbB 0.6434 |
| 12 | 0.6257 | 0.5352 | 0.4256 | 0.6964 | 0.3445 | wBm 0.3445 | wAlbB 0.6964 |
| 24 | 0.6900 | 0.6244 | 0.5493 | 0.7728 | 0.4485 | wBm 0.4485 | wAlbB 0.7728 |

The worst strain is wBm at every size. The spread between worst and best is
0.41 at n=6, 0.35 at n=12 and 0.32 at n=24: pooling the five foregrounds does
not equalise them, and the gap narrows only slowly with size.

**This instance cannot say whether the cause is size or divergence.** wBm is
both the smallest target (1,080,084 bp) and the only supergroup D member. The
two explanations are confounded in this panel and no measurement here separates
them. Deciding it needs a strain that is small and close, or large and distant.

Against *L. plantarum* instead of *Drosophila* the same ordering holds (wBm
worst at 0.3211, 0.4642, 0.6088), so the finding is not specific to one host.

## Q3. Under C1, how does the held-out strain compare?

Each C1 design pools four strains and leaves one out; the fifth is then
evaluated. The comparison is the held-out strain's coverage in C1 against its
coverage in D2, where the same strain was in the foreground. Host *Drosophila*,
set 0.

| Held out | n=6 C1 / D2 | n=12 C1 / D2 | n=24 C1 / D2 | n=12 change |
|---|---|---|---|---|
| wMel (A) | 0.5514 / 0.5514 | 0.6219 / 0.6257 | 0.6924 / 0.6900 | -0.004 |
| wRi (A) | 0.4697 / 0.4697 | 0.5394 / 0.5352 | 0.5932 / 0.6244 | +0.004 |
| wPip (B) | 0.2853 / 0.3057 | 0.3903 / 0.4256 | 0.4445 / 0.5493 | -0.035 |
| wBm (D) | 0.1937 / 0.2375 | 0.2690 / 0.3445 | 0.3624 / 0.4485 | -0.076 |
| wAlbB (B) | 0.3087 / 0.6434 | 0.3737 / 0.6964 | 0.4637 / 0.7728 | -0.323 |

Two regimes, not one.

Leaving out wMel or wRi costs nothing measurable: at n=6 the delivered set is
byte-identical to the pooled D2 set (see "Identical sets"), and at n=12 and
n=24 the change is within 0.004 in either direction. Each is covered by the
other, so supergroup A is represented whichever of the two is dropped.

Leaving out wAlbB costs 0.32 to 0.35 at every size, which is the largest effect
in this record. wPip remained in the foreground and is the same supergroup, so
within-supergroup representation did not substitute for wAlbB here. That is the
case the plan named as the one that matters for a strain not yet sequenced, and
on this panel it is a large loss rather than a small one.

The size of the claim: five leave-one-out designs, one seed each, one host, one
chemistry. The spread across the five held-out strains (0.004 to 0.323) is
itself the result, and it is not predicted by supergroup membership alone.

## Q4. Under D1, how does density vary across hosts that were not in the design?

Reference strain wMel against each of the three hosts, n=12, set 0. The host
that was in the design is marked. hg38 is on the counts route and was in no
design.

| Design | Host | Host bp | In design | Exact sites | exact-site density | C2 exact-site median [min, max] | occupancy-weighted density | C2 weighted median [min, max] |
|---|---|---|---|---|---|---|---|---|
| D1 vs *L. plantarum* | *L. plantarum* | 3.3 M | yes | 3 | 413.81 | ceiling in 3/20 | 21.03 | 14.33 [11.42, 16.53] |
| D1 vs *L. plantarum* | *Drosophila* | 144 M | no | 717 | 74.31 | 25.95 [8.68, 91.14] | 7.87 | 4.26 [2.86, 8.12] |
| D1 vs *L. plantarum* | hg38 | 3,298 M | no | 22,993 | 53.18 | 17.17 [6.44, 48.36] | 5.28 | 3.07 [1.58, 4.80] |
| D1 vs *Drosophila* | *L. plantarum* | 3.3 M | no | 5 | 272.06 | 61.41 [14.34, 212.63] | 20.45 | 6.43 [3.76, 16.59] |
| D1 vs *Drosophila* | *Drosophila* | 144 M | yes | 141 | 414.07 | 88.64 [39.59, 149.57] | 26.51 | 8.10 [6.32, 20.47] |
| D1 vs *Drosophila* | hg38 | 3,298 M | no | 6,096 | 219.80 | 47.41 [22.64, 120.24] | 11.33 | 5.57 [3.91, 11.20] |

Three readings, and the third one revises what two hosts suggested.

**Which host was in the design matters more than any single density.** On hg38,
which neither design saw, the panel designed against the 144 Mb host reaches
exact-site density 219.80 and the one designed against the 3.3 Mb host only
53.18: a factor of 4.1 on the same host, from the same candidate pool, at the
same size. The exact site counts behind that are 6,096 against 22,993. The
weighted densities differ by a factor of 2.1 in the same direction. So on this
instance the host a design is given should be the largest one the design will
meet, and giving it the small host costs most of the specificity against the
large one.

**The transfer is still directional.** Designing against *Drosophila* puts the
panel outside the control spread on both unseen hosts (272.06 against a C2
maximum of 212.63 on *L. plantarum*; 219.80 against 120.24 on hg38). Designing
against *L. plantarum* leaves the panel INSIDE the control spread on
*Drosophila* (74.31 within [8.68, 91.14], and 7.87 within [2.86, 8.12]) -
indistinguishable there from a random panel from the same pool.

**But "indistinguishable from random" does not extend to hg38, which is what
the two-host table wrongly suggested.** On hg38 the *L. plantarum*-designed
panel is 53.18 against a C2 maximum of 48.36, and 5.28 weighted against a
maximum of 4.80: outside the control spread on both densities, though by only
1.10 times the maximum against 1.83 times for the *Drosophila*-designed panel.
So every D1 panel here is measurably better than every control panel on hg38;
what differs between them is the margin, not the sign. The earlier two-host
reading of this record overstated the case and is corrected here.

A likely reason the margin survives on hg38 and not on *Drosophila*: hg38 is
the only host on which no control panel binds nowhere (zero ceiling seeds in
every cell), so its control spread is a complete set of measured ratios and is
narrow, whereas the *Drosophila* spread is wide enough to contain the panel.
That is an observation about the controls, not a mechanism.

Host site densities behind these figures, exact counts, n=12: D1 vs
*L. plantarum* gives 0.896 sites/Mb on *L. plantarum*, 4.989 on *Drosophila*
and 6.971 on hg38; D1 vs *Drosophila* gives 1.493, 0.981 and 1.848.

### Which host sets the worst case

With all three hosts measured, the worst target-against-host density of a panel
is almost always set by hg38, the host no design was given:

| Density | hg38 is the worst host | *L. plantarum* is | *Drosophila* is |
|---|---|---|---|
| occupancy-weighted | 30 of 30 panels | 0 | 0 |
| exact-site | 23 of 30 panels | 7 of 30 | 0 |

The seven exceptions are all at n=6 or n=12 on pooled or leave-one-out designs
against *Drosophila*, where *L. plantarum* carries so few sites that its exact
ratio is noisy; the weighted density, which is the comparable figure, puts hg38
worst everywhere.

This is NOT an answer to Q5. Q5 asks which host determines the figure for a
design whose background POOLS hosts of very different size, and no design here
had hg38 as a background. What this table says is narrower and still useful: a
panel designed against a bacterial or insect host, and then measured against a
mammalian one, has its worst case set by the mammalian host in every case
measured here.

## Q5, and the full D3: still not answered

Every panel is now measured against hg38, but **no design in this record has had
hg38 as a background.** Q5 asks which host determines the pooled figure for a
design whose background pools hosts of very different size. The D3 here pools
only the two small hosts (3.3 Mb and 144 Mb), a factor of 43, and its delivered
set is identical to the D1-vs-*Drosophila* set at 6 and 12 primers, so it is not
a second measurement at those sizes. The D3 rows are labelled preliminary
throughout `results.json`.

Why no design was given hg38 as a background: one `filter` run against hg38
peaks at about 8.5 GB (Known Issue 1 and the Phase 3 cost note), which is
outside the memory plan of this stage. It was never attempted, and hg38 was
never given to `filter`, `optimize` or a position scan.

What the hg38 evaluation does and does not support: the worst-host table under
Q4 says hg38 sets the worst case in 30 of 30 panels by the comparable density.
That is an after-the-fact measurement of panels selected without it. It says
nothing about what a design would select if hg38 were in its background, which
is what Phase 7 would build and what gate A has to weigh.

## The hg38 k-mer table

Built through the pipeline's own counting step, `neoswga count-kmers`, with
hg38 named as the foreground and no background, into the script's evaluation
directory. KMC is not installed here, so jellyfish counted it; `kmer_counter`
was not set.

**Two attempts, and the first is the more informative measurement.**

| | First attempt | Second attempt |
|---|---|---|
| Child ceiling | 6.00 GiB | 7.50 GiB |
| Outcome | stopped by the watchdog at 46.7 s | **ok**, 116.25 s |
| Peak, whole step | 6.18 GiB | 6.21 GiB (6,667,534,336 B) |
| Peak while the `neoswga` process was alone | 6.18 GiB | **6.210 GiB**, at t = 23.6 s |
| Peak once jellyfish was running | not reached | **0.168 GiB** (180,027,392 B), at t = 109.7 s |
| jellyfish first seen | never | t = 24.6 s |
| System free, before / lowest / after | 66 / 36 / 68 percent | 75 / 49 / 73 percent |
| Table | none | 132.0 MiB (138,391,683 B), 8,368,476 distinct canonical 12-mers |

The phase split is the result. **The entire memory cost of counting hg38 is the
Python-side genome load, not the counter.** jellyfish at k = 12 peaked at 180 MB,
about 2.7 percent of the step's peak, and the counting itself took about 92 of
the 116 seconds.

The 6.21 GiB is accounted for arithmetically: two Python string copies of a
3,298,430,636 bp reference are 6.144 GiB, against 6.210 GiB observed, a 1.1
percent difference. `GenomeLoader.load_genome` parses the 705 records into a
list and then forms `"".join(sequences)`, so the per-record list and the
concatenation are both live at that moment; the later `.upper()` copies are a
second instance of the same shape. `count-kmers` reaches this twice, once to
resolve `genome_gc` (`neoswga/core/parameter.py`) and once in
`check_genome_inputs` (`neoswga/core/pipeline.py`). **Setting `genome_gc` in
params.json removes one of the two call sites but not the peak**, because the
peak is inside a single load.

Consequences worth separating. The table size matches the 138 MB recorded in
Known Issue 1 almost exactly, so that figure transfers from KMC to jellyfish;
the 44 s at 10 threads recorded beside it does not, and 116 s at 4 threads with
jellyfish is the figure measured here. The memory of this step was never
recorded before; it is 6.21 GiB, and it is a property of the loader rather than
of the counter or of k.

No product code was changed. Recorded at
`tests/validation/genomes/diversity_baseline/eval_counts/human_full/`, with the
per-sample trace on the step record.

## Reading hg38: what is measured and what is not

hg38 is admitted on the counts route only, and the record of every reference
says which route answered for it. For hg38 that is
`site_source: counts`, `route: "k-mer counts; no positions"`.

Measured: exact site counts, sites per Mb, both pair densities and the
occupancy-weighted load (298,959.88 for the D1-vs-*Drosophila* 12-primer set,
available because the table exists).

Not available, with the reason carried on each figure rather than a zero:
coverage, mean gap, maximum gap and gap Gini, all reading "sites were counted
rather than located, so no positional figure exists for this reference; scan it
to measure coverage". hg38 host coverage is therefore absent from this record,
which is independently the right outcome: see "Host coverage is not tabulated
across routes".

## The counts route was checked before being relied on

The counts route and the position route answer the same question about exact
binding sites at very different cost, and only the position route can also
answer a positional question. That was checked on this data rather than assumed,
because every hg38 figure would have depended on it.

One panel, the delivered D1 wMel-vs-*Drosophila* 12-primer set 0, measured
against *Drosophila* both ways:

| Route | `site_source` | Sites | Coverage |
|---|---|---|---|
| Position index of the design | `positions` | 141 | 0.005752 |
| k-mer counts, `eval/` prefix | `counts` | 141 | unavailable, with the reason |

Exact agreement, difference 0.0. Coverage is `unavailable` on the counts route
with "sites were counted rather than located" as its reason, not 0.0.

The counts route is the one that asks both spellings of each k-mer and keeps
the larger, because a canonical table stores one of each reverse-complement
pair and answers 0 for the other
(`reference_panel_evaluation._measure_counts`). The record confirms it: the
basis string reads "k-mer counts, canonical, both strands", and the equality
with the position route is the measurement that the canonical 0 is not being
taken as an answer.

**One thing this check established that was not expected.** `scan=False` alone
does not put a reference on the counts route: `evaluate_reference_panel` takes
positions whenever an index exists for the prefix, whatever `scan` says. The
first attempt at this check therefore measured positions twice and agreed with
itself. The counts route is reached by naming a prefix that has a table and no
index, which is what the `eval/` prefixes are. The script now refuses a
counts-only reference whose prefix has a position index, rather than reporting a
route it did not take.

Cost of the check: 2.1 s, 132.5 MiB peak.

## Identical sets

Several designs deliver the same primers. Verified by comparing the sorted
primer lists, not by comparing figures: identical sets are ONE observation
reported several times, and reading them as independent would overstate the
evidence.

| Primers | Designs delivering exactly this set |
|---|---|
| 6 | D1 wMel vs *Drosophila*; D3 preliminary |
| 12 | D1 wMel vs *Drosophila*; D3 preliminary |
| 6 | D2 pooled vs *Drosophila*; C1 without wMel; C1 without wRi |
| 12 | C1 without wMel; C1 without wRi |

At n=24 the D1-vs-*Drosophila* and D3 sets differ (host site densities on
*Drosophila* 1.837 against 1.934 sites/Mb), so the D1/D3 identity holds at 6
and 12 primers only. This confirms the earlier report.

The D1/D3 identity has a plain cause: D3 adds *L. plantarum* to the background
of a design that already had *Drosophila*, and at 6 and 12 primers that did not
change what was selected. Every D3 figure at those two sizes is therefore the
corresponding D1 figure and not a second measurement. Likewise, dropping wMel
or wRi from the pooled foreground left the delivered set unchanged at n=6,
which is the strongest form of the Q3 result for supergroup A.

## Against the project's recorded figure for this pool

The recorded figure for the *Wolbachia* pool is selectivity density 60.112 at
coverage 0.6535 for 12 primers, from `plan-pool` with a selectivity-density
floor of 60 and occupancy-weighted coverage
(`docs/validation/frontier_refill_2026-09-17.md`,
`docs/validation/no_search_headroom_on_this_pool_2026-09-18.md`).

**The pool is the same.** The D1 wMel-vs-*Drosophila* candidate pool was
compared with `examples/wolbachia_pool_design/work/step3_df.csv` as CSV: 2,000
candidates on both sides, 2,000 shared, same set and same order.

**The path is not the same, so the figures are not comparable.** This run at
n=12 delivers occupancy-weighted density 26.511 at coverage 0.7381, against the
recorded 60.112 at 0.6535. The recorded figure was produced by `plan-pool`
under a density floor of 60, which is a constraint that held the density up;
this run is `optimize` with `hybrid` and no panel limit set, so nothing stopped
the search from trading density for coverage, and it did: 0.085 more coverage
for 33.6 less density. The two numbers are the same quantity (the
occupancy-weighted density) on different paths with different constraints, and
the difference between them is not a regression. An achievability figure
measured with a constraint in force bounds nothing on a path that has none.

The exact-site density for the same panel is 414.075, which is a different
quantity again and must not be set beside 60.112.

## Decision gate A

The plan states two criteria. They are read here against the figures above and
nothing else.

**Target diversity: the criterion for skipping Phase 6 is not met.** The plan
skips the design-time target work if the worst strain under D2 is within a few
points of the pooled coverage and C1 shows the same. Under D2 against
*Drosophila* the worst strain (wBm) is 0.31, 0.28 and 0.24 below the reference
strain at 6, 12 and 24 primers, and 0.41, 0.35 and 0.32 below the best strain.
Under C1 a held-out strain loses between 0.004 and 0.32 at 12 primers, and the
loss is not predicted by supergroup membership. Neither is within a few
points. So a design made the current way does lose coverage across these
strains, pooling the foregrounds reduces the loss without removing it, and the
levers in Phase 6 (a per-genome candidate gate, and a worst-target figure that
selection can see) have a measured loss to act on. Which lever is built first
is not decided by this record; the shape of the loss says the worst target,
not the pooled figure, is the quantity to act on.

**Host panel: not decided.** The plan skips Phase 7 if per-host density under
D3 is within the spread of C2 across hosts. The D3 that criterion names pools
hosts of very different size, and no design here had hg38 as a background, so
the criterion cannot be read. What is measured bears on it without settling
it: the host a design is given changes its exact-site density on hg38 by a
factor of 4.1, a panel designed against the 3.3 Mb host is inside the control
spread on the 144 Mb host, and adding the 3.3 Mb host to the 144 Mb one changed
no delivered set at 6 or 12 primers. Those are reasons to expect that a pooled
background is decided by its largest member. They are not a measurement of it.
Settling it needs one design with hg38 in the background, which is one `filter`
run of about 8.5 GB.

The size of the claim for both halves: five strains of one genus, three hosts,
one chemistry in a regime the repository describes as unfavourable for
discrimination, three panel sizes, one optimizer seed, 20 control seeds.

## What is not established

- **No design has had hg38, or any host-sized reference, as a background.**
  Every hg38 figure here is an after-the-fact measurement of a panel selected
  against a 3.3 Mb or 144 Mb host. What a design would select with hg38 in its
  background is not measured, and the worst-host table does not stand in for
  it.
- **No coverage or gap figure on hg38.** The counts route cannot produce one.
  So the Q1 to Q3 coverage results rest on the five targets only, and nothing
  here says how evenly a panel tiles a mammalian genome.
- **Size against divergence for wBm.** Confounded in this panel; see Q2.
- **One seed.** Each design was optimized once, at seed 0. No claim here
  distinguishes a property of the optimizer from a property of one run, except
  where the C2 controls (20 seeds) are quoted beside it.
- **One chemistry.** phi29 at 30 C throughout. Known Issue 17 records that this
  pool cannot discriminate at this temperature, so the specificity figures here
  are those of a regime the repository already describes as unfavourable.
- **C2 panels are not dimer-screened.** They are uniform draws from
  `step3_df.csv` and answer "what does any panel from this pool do", not "what
  does any orderable panel do".
- **The full D3** was not built. The D3 rows are the two small hosts only.
- **No per-target coverage floor was set**, so `min_per_target_coverage` was
  neither checked nor repaired in any of these runs.
- **Mismatch discrimination is not claimed.** The occupancy-weighted load uses
  an exact-match model with a uniform assumed per-mismatch correction.

## How to reproduce

The designs, from an empty gitignored work directory (about 9 GB, several hours,
one heavy step at a time under the script's watchdog):

```
python scripts/benchmarking/diversity_baseline.py \
    --manifest tests/validation/genomes/diversity_panel.json \
    --workdir tests/validation/genomes/diversity_baseline
```

The counts-route agreement check, which appends to the results file:

```
python scripts/benchmarking/diversity_baseline.py \
    --counts-route-check D1__wolbachia__vs__drosophila \
    --counts-route-check-host drosophila --counts-route-check-size 12
```

Every table in this record, read from `results.json` and recomputing nothing:

```
python scripts/benchmarking/diversity_baseline.py --tables
```

No number in this record was typed from memory; each was read from
`results.json` by that option, which prints Tables 1 to 8 in the order used
here.

The hg38 figures, which is how they were produced here. The first invocation
builds the table and stops so its peak can be read before anything else starts;
the second measures every recorded panel against it. The ceiling of 7.5 GiB is
needed because the step peaks at 6.21 GiB; at 6.0 GiB the watchdog stops it (see
"The hg38 k-mer table").

```
python scripts/benchmarking/diversity_baseline.py \
    --counts-only-hosts human_full --counts-only-pass --counts-only-tables-only \
    --start-free-percent 65 --rss-limit-gb 7.5 --halt-peak-gb 7.5 \
    --stop-free-percent 25

python scripts/benchmarking/diversity_baseline.py \
    --counts-only-hosts human_full --counts-only-pass \
    --start-free-percent 65 --rss-limit-gb 7.5 --halt-peak-gb 7.5 \
    --stop-free-percent 25
```

That pass reads `results.json`, adds figures to it and writes it back; it
builds no design and rebuilds no pool, and it refuses a results file that
records no design. A host longer than 1 Gb is admitted on the counts route and
nothing else, by its length rather than by a flag: it is kept out of the host
list the designs are planned from, `design_params` refuses parameters naming
one, and `reference_plan` forces the counts route and refuses a prefix that has
a position index.

What that pass cost here: 45 evaluation children (30 set-0 panels and 15 C2
batches of 20 seeds), about 42 minutes in total, each child peaking at about
1.7 GiB. The peak is the hg38 k-mer table held as a dictionary of 8,368,476
entries by `mismatch_counts.load_kmer_counts`, which the occupancy-weighted
load needs; the exact site counts alone stream the table and do not hold it.
System free memory stayed between 73 and 76 percent throughout.

The pass only adds. After it, `results.json` still held all ten designs, all
three step records for each, all 30 original set-0 evaluations and all 15
original C2 blocks, alongside the 30 and 15 new counts-only blocks; the pass
asserts this on every write and refuses a results file recording no design. The
C2 draws it re-measured were checked against the panels the earlier pass
actually evaluated, read from that pass's own output, and matched for all 15
batches.

## The owed real-data check of `improve-set`

`improve-set` was merged with its real-data verification deferred. Run once on
the delivered D1 wMel-vs-*Drosophila* 12-primer set 0, `--max-edits 10`, output
to a scratch directory outside the repository.

Cost: 25.0 s wall, 1.34 GiB peak RSS.

Sections, each stating how many it shows of how many it considered:

| Section | Shown of considered |
|---|---|
| Drops that cost no coverage | 0 of 0 |
| Adds that raise the worst target | 10 of 169 |
| Swaps that raise the worst target | 1 of 1 |
| Trade-off drops | 5 of 5 |

Candidate pool: 1,988 examined, 1,461 conflicting with an oligo that stays, 0
outside the Tm window, 0 whose sites could not be established.

**No proposal raises the worst target without costing specificity.** All 11
proposals with a positive worst-target gain (10 adds and 1 swap) lower the
occupancy-weighted worst host selectivity density below its current 414.075,
and all 11 raise the host exact site density above its current 0.9810 sites/Mb.
The range of the density change is -3.116 (add GCATGGTAGTTA) to -61.072 (add
ATAAAGAAGTAG). So the case the plan called a finding about the optimizer -- a
positive worst-target gain whose predicted worst host density is not lower than
the current one -- does not occur on this set. On this panel the delivered set
is not leaving a free coverage gain behind.

Taking the top add (ATCAAGAAGTAG) and re-running `evaluate-set` on the
13-primer set reproduces every predicted figure exactly: target sites 521,
target coverage 0.7475520239284041, target weighted load 1330.6614132904504,
host sites 153, host sites/Mb 1.0645255407577539, host coverage
0.006252953449578317, exact-site pair density 386.04422156742214, and both
reductions. Eleven figures compared, all identical to full precision.

As the report says of itself, that agreement is expected and is not an
independent check: the prediction is produced by the same evaluation code that
`evaluate-set` then runs. What it does establish is that the seam carrying the
prediction (`PanelSources`, the shared cache `improve-set` passes) returns what
the uncached path returns, which is the thing that could have differed.

Side effects on the design directory: none. No `design_failure.json` was
written anywhere under it, and a directory listing before and after the two
commands is identical. The known `genome_gc.json` side effect of `evaluate-set`
and `improve-set` did NOT occur here, because the design's `params_n12.json`
already states `genome_gc` and the caching branch is only reached when it is
absent.
