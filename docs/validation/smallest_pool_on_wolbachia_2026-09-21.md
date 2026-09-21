# Deletion stops one oligo short on the Wolbachia design

> **SUPERSEDED AND WITHDRAWN, 21 September 2026, the same day it was written.**
> The central claim below does not reproduce. Re-running the same instance with
> the greedy panel and the beam produced by one script returns a greedy
> coverage of **0.7334**, not the 0.7529 quoted here, and at that baseline the
> beam's best eleven reaches only 0.7239 -- short. 0.7334 is what two records
> written before this one already said for the same panel, in
> `occupancy_and_discrimination_2026-09-19.md` and Known Issue 17.
>
> So the panel this measurement asked the beam to beat was not the panel the
> design delivers, and beating a weaker twelve with an eleven is an easier
> problem. It was measured inline with no script kept, so what it actually did
> cannot be recovered.
>
> Four further instances agree that there is no saving. See
> [beam_does_not_beat_deletion_2026-09-21.md](beam_does_not_beat_deletion_2026-09-21.md).
>
> Kept rather than deleted because the correction is the useful part: a
> baseline that disagrees with the project's own figure for the same quantity
> is reporting on a different object, and checking it costs nothing.

Measured 21 September 2026, on `examples/wolbachia_pool_design/work`: wMel
(1,267,782 bp) against *Drosophila* (143,726,002 bp), 12-mers, phi29 at 30 C,
3 kb reach, the 2,000-candidate shortlist the filter produced.

## Result

`--minimize-primers` reduces a panel by removing one oligo at a time and stops
when no single removal still meets the coverage target. On this design it
removes nothing at all. A beam search over the same pool finds a strictly
smaller panel with strictly higher coverage.

| Method | Size | Coverage | Selectivity density |
|---|---|---|---|
| greedy, then deletion | 12 | 0.7529 | — |
| beam, size 11 | **11** | **0.7599** | 28.06 |
| beam, size 10 | 10 | 0.7485 | 28.10 |

The target is the greedy's own coverage, 0.7529. Eleven oligos clear it; ten
fall 0.0044 short. So the smallest qualifying size is 11 against deletion's 12:
one oligo of twelve, for slightly more coverage rather than less.

Coverage RISING as the panel shrinks is the point. The greedy's twelve is not
an optimal twelve, so a better eleven exists inside the same pool. Deletion
cannot reach it because getting there requires exchanging a member before
dropping one, and deletion only drops.

## Method

The beam ran over the delivered panel plus the first 100 candidates in the
shortlist's own order, 108 distinct oligos, at beam width 4. That bound is for
cost, not for correctness: one `compute_metrics` call is 40 ms on this design,
so a 12-level beam over the full 2,000-candidate pool would be about an hour.
A wider pool can only find the same panel or a better one, so 11 is an upper
bound on the smallest qualifying size, not a claim of optimality.

`scripts/` holds no runner for this; it was measured inline.

## What this changes

`tests/test_smallest_pool_search.py` demonstrated this gap on a constructed
four-candidate fixture and recorded that it did not reproduce on either
instance then examined -- the bundled 6 kb plasmid, which one primer covers
completely, and a 300 kb random sequence, which has no dominance structure.

**That record was wrong, and wrong because the search was too narrow.** The
prepared Wolbachia design was already in the repository, is the design the
project's own measurements are quoted from, and shows the gap on the first
attempt. Neither of the two instances tried was representative; one was a
toy and the other was synthetic.

The deferral rested on "no demonstrated benefit". That basis is now gone.
Wiring the beam into `optimize --minimize-primers` would deliver an 11-oligo
panel where it currently delivers 12, which is a real reduction in synthesis
cost and a real change to every delivered result, so it is a decision rather
than a correction.

## Caveats

One design, one panel size, one reach. Whether the saving generalises is
unmeasured, and the three bacterial genomes in `tests/validation/genomes`
would be the obvious next instances. Nothing here says the beam is the right
mechanism either; a swap pass before deletion might reach the same panel more
cheaply.
