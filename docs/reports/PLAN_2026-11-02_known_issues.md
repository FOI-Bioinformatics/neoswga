# Plan: the Known Issues, triaged and the open ones scheduled

Written 2 November 2026, on `main` at `030d1e8`, against Known Issues 1-25 in
[KNOWN_ISSUES.md](../development/KNOWN_ISSUES.md).

**The first finding is that "fix the known issues" is the wrong frame for most
of the list.** Of the twenty-five entries, nine are closed records of defects
already fixed, nine are established facts or constraints where the entry IS the
answer, and four are deliberate positions that evidence has not moved. Treating
those as a backlog would mean reopening decisions this project made on
measurement, which its own conventions forbid: *defaults follow measurement*,
*panel limits have no defaults because no threshold has a reference*, and *do
not record this class as closed*.

Three entries are genuinely open engineering work. They are what this plan
schedules.

## Triage

| Entry | Status | This plan |
|---|---|---|
| 2, 8, 9, 10, 12, 13, 15, 19, 22 | **FIXED**, dated, with the defect recorded | Nothing. The entries are history and the numbering is cited by code and tests |
| 4, 7, 23, 24 | **Invariant or environment** — int64 coordinates, the strand vocabulary, Python 3.13, KMC3 preference | Nothing. Each is enforced by a test or a refusal already |
| 5 | **Upstream, unfixed** (pyahocorasick past 2 Gb) | Nothing beyond the canary that already fails on a version change. Not ours to fix |
| 6, 17, 20, 21 | **Established by measurement**; the entry is the finding | Nothing. 17's remedy is the reaction, and gating candidates was measured to make panels worse |
| 11, 16, 18 | **Deliberate positions**, off or unscored after measurement | Nothing without a new measurement. Phase 4 took the one cheap strengthening available |
| **3** | **Open**, and the entry named the wrong path | **Phase 1.1 DONE** — 0.22% of the table is read, ~1.2 GB at host scale; but 1.2 as specified would be 43,392 KMC calls, so it is NOT written. See the record |
| **14** | **Open and explicitly not closed**; reason now CHANGED | **Phase 2 DONE** — 2.1: field costs +48%/evaluation, 2.05x cheaper than the fallback already in use. 2.2: field added, key runs, axis measured — the sign of its effect depends on panel size, so it stays off |
| **25** | **Open by decision**: the symmetric window overstates coverage 1.3-1.9x | **Phase 3**, and it needs your calls |

## Startability in this container

`tests/validation/genomes/` is gitignored and absent, and
`examples/wolbachia_pool_design/work` is not prepared — so **hg38, the
*Drosophila* host and the prepared Wolbachia pool are all unavailable**, and
NCBI is refused by the egress proxy. Both k-mer counters are installed at
`/root/mm/envs/tools/bin` (export it onto `PATH`).

That gates the plan sharply:

- Phase 1's headline figure is peak memory on a host genome, which **cannot be
  measured here**. A synthetic table of comparable cardinality can stand in for
  the mechanism but not for the number.
- Phase 2 is **fully startable**: it is a per-evaluation cost measurement on a
  synthetic pool, which is how 6.1 and 7.2 were done.
- Phase 3 needs a real reference and, for its last step, a BAM. **Not
  startable.**

## Phase 1 — `filter` stops holding the whole background (Known Issue 3)

The entry's own workaround is "use a Bloom filter", which is a different code
path with its own ceiling (Known Issue 22). The direct fix is already named in
the code that causes it. `mismatch_counts.load_kmer_counts` says:

> This still materialises the whole table as a dict, which is its own ceiling
> on a host genome; asking only for each candidate's neighbours
> (`kmer_tables.counts_for`) is the change that removes it, and it is a
> restructuring of the callers rather than of this function.

**1.1** Measure first. Peak RSS and wall time for `filter` against a background
table at three cardinalities, and the count of distinct canonical forms a run
actually looks up against the table's size. If the lookups are a small fraction
of the table, the dict is pure overhead and the fix is justified by that ratio
alone — which is measurable here on a synthetic table.

**1.2** Move the callers to `kmer_tables.counts_for`, batched per k as the "one
door" rule requires (*a count lookup has a fixed cost per call, so batch by
k*). `load_kmer_counts` keeps its `lru_cache` contract for callers that
genuinely want the whole table; the point is that the mismatch path stops being
one of them.

**1.3** Re-measure. The claim to earn is a memory reduction at equal delivered
output — so a control asserting the same shortlist, byte for byte, matters more
than the memory figure.

**Risk**: `counts_for` per-neighbour lookups could be slower than one dict
build when the table is small and the candidate pool large. 1.1's ratio decides
which regime each run is in; if both exist, the threshold belongs in
`kmer_tables` as one decision in one place, like `lazy_dimer.dimer_screen`.

## Phase 2 — what a host-coverage field would cost (Known Issue 14)

This is the highest-value item, because one cheap measurement can close a
long-open decision either way.

The blocker is precise and already recorded. `stage2_background_weight` is
REFUSED above zero because the weighted axis needs a host-coverage fraction and
`PoolMetrics` deliberately carries five fields and no sixth — its docstring:
*"A caller that reaches for a sixth gets an `AttributeError` rather than a
default that would be indistinguishable from a measurement."* The reason for
five is cost: `compute_metrics` was 27.6 s of a 29.7 s run across 4,040 calls.

So the question is not "should the host term be weighted" but **"what does
carrying host coverage through the search cost per evaluation"**. That is
measurable.

**2.1** Measure the marginal cost of a `bg_coverage` field on `PoolMetrics`:
per-evaluation time with and without, over a search of known evaluation count,
on a pool with a real background. Report it as a fraction of total search time,
not in isolation.

**2.2** Then decide, and both outcomes are results:

- **Affordable** → add the field, lift the refusal, and re-run the measurement
  `2026-10-29-stage2-background-weight-does-not-run.md` could not make: does a
  weighted host axis beat the lexicographic tie-break on a reaction that
  discriminates? Known Issue 14 then closes on evidence, or stays open with the
  weight off and a measured "no benefit".
- **Not affordable** → the refusal stands, and entry 14 gains the cost figure
  it currently lacks. That converts "deliberately carries none" from an
  assertion into a number.

**Do not** reach for the cheap substitute. Re-basing the axis on
`total_bg_sites` instead of a fraction is the cross-scale fudge
`2026-10-26-stage2-background-weight.md` claimed to avoid and was already
rejected once.

## Phase 3 — the coverage geometry (Known Issue 25)

**The path is already written** and this plan does not restate it: the four
steps a switch needs are in
[`2026-10-29-the-coverage-geometry-decision.md`](../validation/2026-10-29-the-coverage-geometry-decision.md),
under "What remains for whoever decides to switch".

What this plan adds is that **three of those four steps are decisions, not
tasks**, and they are yours:

1. **Does the default move at all?** Every recorded `fg_coverage` in
   `docs/validation/` was produced under `symmetric`, and the 3 kb reach was
   fitted against it. Switching means re-baselining a research record.
2. **Scope by polymerase?** The directional band (4,354-6,710 bp) was fitted
   for phi29 against one published outcome. There is none for `bst`,
   `equiphi29` or `klenow`. Either fit three more, or scope the switch to phi29
   and keep refusing the rest — which is what the reach refusal does today.
3. **What does `--target-coverage 0.70` mean across the change?** A 0.70 target
   is a different request under a figure 1.3-1.9x smaller, and so are the panel
   limits.

Only step 1 of that record — teaching Stage 1 the directional geometry and
measuring the panel it then delivers — is engineering, and it is blocked on a
real reference genome.

**The status quo is defensible.** The defect is disclosed in every artifact
(`coverage_geometries` in the summary, both figures in the log), a directional
run must name its own reach, and nothing claims the symmetric figure is the
directional one. Doing nothing here is a position, not a backlog item.

## Phase 4 — the cheap strengthening available on a deliberate position

Not a fix; an honesty improvement.

**Corrected 2026-11-02, on starting it.** This section said entries 11 and 20
"both rest on one measurement". That is wrong about entry 11, which rests on
THREE real pools at ten combinations of tier and panel size, where a 0.9
threshold fired 328,846 times and left the Hamming-1 pair count and the
duplicate 3' hexamer count identical in all six comparable cases. That is a
well-sized measurement and needs no strengthening; reading the index line
("measurement did not support enabling it") rather than the entry is what
produced the error.

Entry 20 does rest on one pair: Prevotella against human chr21, mixed k 10-12
against the single lengths, at panel sizes 8 and 20. The project's own
convention calls one pair "no benefit demonstrated", which the entry correctly
claims. A second pair either holds the claim at a larger size or overturns it.

Cheap, startable, and low value — do it only if Phases 1 and 2 are done and
something is waiting on a reference genome.

**DONE 2026-11-02**, and it changed the entry more than expected. The headline
("choosing k matters far more than choosing whether to mix") replicates on a
second pair with a 290x ratio between the two effects. The subsidiary reading
("lands between the single-length designs") does not generalise, because on the
second pair the mixed arm did not mix: the k 10-12 shortlist held 0 ten-mers and
2 eleven-mers out of 300, so `max_primer` plus the candidate ranking had chosen
the length before any optimizer ran. Entry 20 now says a length RANGE is an
option the pipeline may decline, and that whether it was exercised is visible in
the shortlist's composition — which the funnel does not report.
[The record](../validation/2026-11-02-mixed-length-on-a-second-pair.md). No code
changed. The second pair is synthetic, since no reference genome is reachable
here.

## Sequencing

1. **Phase 2.1** first. It is startable today, it is one measurement, and it
   decides an entry that has been open since September.
2. **Phase 1.1** next, for the same reason: a ratio measurable on a synthetic
   table tells you whether 1.2 is worth writing.
3. **Phase 1.2-1.3** if 1.1 justifies it.
4. **Phase 3** only once a reference genome is reachable, and only after you
   have answered its three decisions.
5. **Phase 4** as filler. **Done** — see that section.

## What this plan does not do

- **It does not treat the closed entries as work.** Nine are fixed and dated;
  their numbers are cited by code and tests and must not be renumbered or
  retired.
- **It does not reopen a measured position.** 11, 16, 17, 18 and 20 are off or
  unscored *because* measurement did not support them. Switching one on without
  a new measurement is the thing the conventions forbid.
- **It does not claim Known Issue 8 can be closed.** The entry says not to, and
  the class produced two new instances within the last week — a guard missing a
  door, and a merge reaching one of two paths. The standing defence is the six
  reachability ratchets plus `test_no_adapter_drops_an_inner_parameter.py`, and
  the work there is to keep extending them, not to declare victory.
- **It promises no wet-lab validation.** Nothing here produces a measured
  amplification outcome, and there is still no BAM or CRAM in the repository,
  so no reach has been fitted against sequencing depth.
