# No search headroom on this pool, with controls this time

Measured 2026-09-18 on the real Wolbachia pair. Closes the question
[the Stage 1 work](stage_one_constraint_awareness_2026-09-18.md) left open:
whether the search is leaving specificity on the table at panel size 12.

**It is not.** No construction that respects the dimer screen beats the search's
density 60.112 at coverage 0.6535, and every density-oriented construction pays
for its density with coverage far below the 0.5 target.

Reproduce with `scripts/benchmarking/selectivity_compatibility_tension.py`.

## The signal I said was worth chasing, retracted

That note observed that only 6 of the 16 most selective candidates are mutually
dimer-free, and suggested selectivity and compatibility might be in tension in a
way that caps achievable specificity. With a control, most of that goes away.

| k | Set | Compatible pairs | Largest dimer-free subset | Mean GC |
|---|---|---|---|---|
| 16 | top by slack | 57.5% | 6 | 0.479 |
| 16 | 12 random sets | 69.1% | 7.1, range 6 to 8 | 0.392 |
| 32 | top by slack | 59.3% | 10 | 0.479 |
| 32 | 12 random sets | 70.9% | 10.9, range 10 to 12 | 0.393 |
| 64 | top by slack | 62.5% | 18 | 0.467 |
| 64 | 12 random sets | 71.4% | 15.8, range 13 to 21 | 0.393 |
| 128 | top by slack | 61.2% | 26 | 0.454 |
| 128 | 12 random sets | 70.0% | 22.8, range 18 to 27 | 0.395 |
| 256 | top by slack | 62.6% | 41 | 0.442 |
| 256 | 12 random sets | 70.3% | 33.5, range 27 to 47 | 0.397 |

Two things are true and they point opposite ways.

**Pairwise compatibility really is lower among selective candidates**, by 8 to 10
points at every k, and the mechanism is visible in the last column: they are
GC-richer, 0.44 to 0.48 against 0.39, and G/C pairs are what a complementary run
is made of.

**It does not limit panel construction.** The largest mutually compatible subset
is no smaller than random at k = 16 and 32, and LARGER at 64, 128 and 256. The
6 of 16 that prompted this sits inside the random range of 6 to 8, so it was
small-sample noise rather than a barrier. Among the top 64 by slack there are 18
mutually compatible primers, comfortably more than a 12-primer panel needs.

## So compatibility is not the barrier. Coverage is.

Nine valid constructions, every pair screened at `max_dimer_bp` 3, taking the
largest compatible subset of the top k by slack and then its 12 most selective
members:

| Built from | Density | Coverage | Dimer-free | Beats the search |
|---|---|---|---|---|
| top 64 by slack at D=65 | 64.008 | 0.3851 | yes | no |
| top 128 | 64.008 | 0.3851 | yes | no |
| top 256 | 64.008 | 0.3851 | yes | no |
| top 64 by slack at D=70 | 63.824 | 0.3851 | yes | no |
| top 128 | 58.325 | 0.3843 | yes | no |
| top 256 | 63.824 | 0.3851 | yes | no |
| top 64 by slack at D=80 | 62.754 | 0.3870 | yes | no |
| top 128 | 56.993 | 0.3756 | yes | no |
| top 256 | 52.372 | 0.3853 | yes | no |
| **the shipped search** | **60.112** | **0.6535** | **yes** | -- |

Every one lands between 0.375 and 0.387 coverage, far below the 0.5 target, and
none exceeds 64.008 density. The search reaches 60.112 while carrying 0.6535
coverage, which is 0.27 more coverage for at most 4 points of density.

**The binding constraint is the specificity against coverage trade-off, not the
search.** A design asked for a density of 65 with coverage above 0.5 is asking
for a point this pool does not appear to contain, and the search declining it is
the correct answer rather than a defect.

## What this settles, and what it cost to settle

Three claims of mine died here, and the sequence is the lesson:

1. A density ceiling of 79.807, which ignored both the coverage target and the
   dimer screen.
2. An existence proof at density 78.178, which carried 30 dimerising pairs out
   of 66.
3. A compatibility barrier at 6 of 16, which sits inside the random range.

Each survived until it met a control or a constraint it had omitted. Three Stage
1 search rules were built on the first two before the third was tested. The
cheap check that would have killed all three at the outset is the same one:
**evaluate a candidate panel through the same acceptance path a delivered panel
goes through, and compare it against a control.**

## What this does not establish

- One pair, one panel size, one reaction. A larger panel, a different target or
  a looser dimer threshold could all move the frontier.
- The nine constructions are greedy and deterministic, so they are lower bounds.
  The search beating them does not prove it optimal; it means no cheap
  construction tried here beats it, and three of my attempts to find one failed.
- The compatibility control uses 12 random sets per k, which bounds the range
  roughly rather than giving a distribution.
- It does not test a looser `max_dimer_bp`. If the GC-richness of selective
  candidates is what costs compatibility, a design willing to accept 4 rather
  than 3 might reach further, and that is untested.
- Coverage and selectivity remain modelled site geometry, not measured
  amplification. See [evidence_matrix_2026-09-15.md](evidence_matrix_2026-09-15.md).
