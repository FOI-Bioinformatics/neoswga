# Mixing lengths on a second pair: the filter decides, not the optimizer

Measured 2 November 2026, on `main` at `e8f4da4`. Phase 4 of
[the plan](../reports/PLAN_2026-11-02_known_issues.md), against Known Issue 20.
Reproduce with:

```
PYTHONPATH=. python docs/validation/mixed_length_second_pair_2026-11-02/measure.py \
    --params <work>/params.json --out results.json --lengths 10,11,12 --sizes 8,20
```

Results and the pool's configuration in
[`mixed_length_second_pair_2026-11-02/`](mixed_length_second_pair_2026-11-02/).
**No code changed.**

## Why a second pair

Known Issue 20 concludes from one pair — Prevotella against human chr21 — that
a mixed k 10-12 design "lands BETWEEN the single-length designs on both axes"
at panel sizes 8 and 20, and that **"choosing k matters far more than choosing
whether to mix."** One pair is what this project's convention calls "no benefit
demonstrated". This adds a second: a 1.27 Mb target at 35% GC against a 4 Mb
host at 42% GC, equiphi29 at 43 C, `max_primer` 300, the same panel sizes.

## The measurement

| arm | n | candidates | fg_cov | bg_cov | bg sites | sel_dens |
|---|---|---|---|---|---|---|
| k=10 | 8 | 300 | 0.7427 | 0.1988 | 111 | 2.635 |
| k=11 | 8 | 300 | 0.5026 | 0.0420 | 21 | 3.372 |
| k=12 | 8 | 300 | 0.2790 | 0.0100 | 5 | 4.791 |
| **mixed 10-12** | 8 | 300 | **0.2774** | **0.0080** | 4 | **4.799** |
| k=10 | 20 | 300 | 0.9687 | 0.4457 | 299 | 2.655 |
| k=11 | 20 | 300 | 0.8061 | 0.0984 | 52 | 3.473 |
| k=12 | 20 | 300 | 0.5441 | 0.0200 | 10 | 4.875 |
| **mixed 10-12** | 20 | 300 | **0.5527** | **0.0340** | 17 | **4.857** |

## The main claim replicates, emphatically

**Choosing k matters far more than choosing whether to mix.** Across the single
lengths, coverage spans 0.46 at n=8 and 0.42 at n=20, and selectivity density
spans 2.16 and 2.22 — an enormous coverage-for-specificity trade, with k=10
reaching 0.97 coverage at 45% host coverage and k=12 reaching 0.54 at 2%. The
mixed arm differs from the NEAREST single length by 0.0016 and 0.0086 coverage.

So the effect of choosing k is two to three orders of magnitude larger than the
effect of mixing, on a second and very different pair. That is the entry's
headline and it holds.

## The subsidiary claim does not, and the reason is the interesting part

"Lands between the single-length designs on both axes" is **false at n=8** here
— the mixed arm sits just OUTSIDE the single-length range on both axes, 0.0016
below the lowest coverage and 0.008 above the highest selectivity density — and
true at n=20. But the margin at n=8 is smaller than the difference between any
two single lengths by a factor of 290, so "between or outside" is not the useful
reading. This is:

| arm | oligo lengths DELIVERED |
|---|---|
| mixed, n=8 | **{12: 8}** — not mixed at all |
| mixed, n=20 | **{11: 2, 12: 18}** |

And upstream of that, the k 10-12 **shortlist itself** held:

| | 10-mers | 11-mers | 12-mers |
|---|---|---|---|
| candidates available to the optimizer | **0** | 2 | 298 |

**A mixed-length design on this pair is decided at the filter, not at the
optimizer.** Given a 10-12 range, the shortlist came back with no 10-mers and
two 11-mers out of 300. The n=20 panel then took both 11-mers and 18 of the
298 12-mers — so the "mixing" is the panel using the only two non-12-mers it
was offered, not the search preferring a blend. At n=8 there was no reason to
reach past the 12-mers at all, which is why that arm is an all-12-mer design
wearing a mixed label, and why it lands on top of the k=12 arm rather than
between anything.

That also explains the direction of the tiny n=8 discrepancy: the mixed arm is
not a different kind of design from the k=12 arm, it is the k=12 arm drawn from
a slightly different 300-candidate pool, because the `max_primer` cut ranked
across three lengths instead of one.

**The mechanism is the shortlist cap plus the candidate ranking**, both upstream
of any optimizer: `max_primer` bounds the shortlist at 300 and the ranking
fills it with the length that has the best foreground-to-background ratio,
which on a 4 Mb host is the longest. The single-length k=10 arm proves 10-mers
can fill a shortlist when they are the only competitor — it returned 300
candidates and the highest coverage of any arm. They simply lose the ranking
when 12-mers are in the same pool. **Which of the cap and the ranking dominates
was not isolated here**; raising `max_primer` would separate them and was not
run.

## What this changes in the entry

- The headline — choosing k over choosing to mix — goes from one pair to two,
  with the second showing a 290x ratio between the two effects.
- "Lands between the single-length designs" should not be restated as a general
  property. On this pair the mixed design lands ON the longest single length,
  because the filter had already chosen it. The honest general statement is
  that **a length RANGE is an option the pipeline may decline**, and whether it
  is exercised is visible in the shortlist composition before the optimizer
  runs.
- `core/length_occupancy.py` already reports per-length occupancy from `filter`
  and `optimize` and never enforces it. Reporting the shortlist's length
  composition in the funnel would make "did my mixed design actually mix"
  answerable without a separate measurement — not proposed as a change here,
  but it is the gap this measurement had to work around.

## Scope of the claim

One synthetic pair, one length range, two panel sizes, one reaction, one
method, one run per arm. Each arm reset the artifacts a previous arm wrote and
kept only the k-mer tables; running arms in a shared directory without that
reset is a confound measured earlier in this session, where the candidate
inventory carried across and the same arm disagreed with itself.

**This is a synthetic pair, not a second real one.** The first pair was
Prevotella against human chr21, real genomes whose k-mer sites cluster; here
sites are Poisson by construction, and the number of distinct k-mers at each
length is near the maximum for that length. That affects exactly the quantity
driving the finding — how often a short oligo hits the host — so the shortlist
monopoly by 12-mers may be sharper here than on real sequence. The direction
(short oligos lose the ranking on a large host) follows from the background
frequency gate and generalises; the completeness of the monopoly is for this
site geometry. No reference genome is reachable from this container
(`tests/validation/genomes/` is gitignored and absent, NCBI is refused by the
egress proxy). Nothing here is a wet-lab outcome.
