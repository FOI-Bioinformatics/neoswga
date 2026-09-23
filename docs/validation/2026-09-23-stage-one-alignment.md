# Aligning Stage 1 with the scoring geometry, measured and not taken

Measured 23 September 2026. **The change was implemented, measured, and
reverted.** This records why, because the first half of the measurement argues
for it and the second half against, and only the second half is decisive.

## The mismatch

`BipartiteGraph.add_primer_coverage` holds both geometries and takes the
directional branch whenever the caller supplies the separated strands.
`optimize_greedy` always did, so **every method's Stage-1 set cover selects
directionally at width `r`** -- while the Stage-2 swap and the reported
`fg_coverage` are symmetric at width `2r`.

Nobody chose that. CLAUDE.md's **Coverage reach (important)** records the
realistic reach being threaded into Stage 1 so that "selection and the reported
`fg_coverage` agree (and ensemble comparisons are fair)". The reach was
reconciled; the geometry was not.

## What aligning them does at a FIXED panel size

Prevotella against human chr21, 12-mers, Stage 1 forced onto its symmetric
branch so it matches the scoring:

| n | Stage 1 directional (today) | Stage 1 aligned | delta | Jaccard |
|---|---|---|---|---|
| 4 | 0.053026 | 0.053026 | 0.00% | 0.600 |
| 8 | 0.090124 | 0.090519 | +0.44% | 0.600 |
| 16 | 0.149496 | 0.149805 | +0.21% | 0.778 |
| 24 | 0.185577 | 0.190676 | +2.75% | 0.778 |
| 32 | 0.209814 | 0.222034 | **+5.82%** | 0.561 |

Never negative, growing with panel size. On its own this reads as a free
improvement, and it is close to expected rather than surprising: a greedy that
optimises the metric it will be judged on should beat one optimising a
different metric.

**That measurement is incomplete, and acting on it alone would have been
wrong.**

## What it does when the panel size is not fixed

`num_primers` is a request. Stage 1 stops once the coverage target is met, so
on a design whose stop is coverage-driven rather than size-driven the panel
size is an OUTPUT. Measured on the 80 kb synthetic pair in
`tests/test_the_repair_scan_is_bounded_in_production.py`, asking for 8:

| | delivered panel |
|---|---|
| Stage 1 directional (today) | **5 primers** |
| Stage 1 aligned to symmetric | **3 primers** |

A 40% smaller panel. The mechanism is immediate: a symmetric window is twice as
wide, so each primer appears to cover twice as much and the greedy reaches its
coverage target sooner.

## Why that decides it

The model that justifies stopping at three primers is the one there is reason
to believe **overstates** coverage. `docs/validation/2026-09-23-directional-coverage.md`
records the mechanism -- a site extends one way -- and
`2026-09-23-does-the-panel-move.md` measures the symmetric figure at about
twice the directional one at equal reach.

So the existing mismatch is **conservative in the direction that matters**.
Stage 1 selecting directionally picks more primers than the symmetric score
would demand, and the symmetric score then reports generously on the larger
panel. Aligning them removes the conservatism and makes a design claim more
coverage from fewer oligos, on the weaker of the two models.

Fewer oligos is cheaper, and if the symmetric figure were right it would be a
real saving. Nothing here establishes that it is right, and the mechanism
argues it is not.

## What was not done, and what would settle it

- **Aligning the other way** -- making the scoring directional so Stage 1 does
  not move -- was measured separately in
  `2026-09-23-does-the-panel-move.md`. It changes delivered panels at two of
  four sizes at equal reach and three of four at the refitted reach, so the
  default stayed `symmetric` there too.
- **Which panel amplifies better** is the question neither measurement can
  answer. `calibrate-reach --bam` fits the reach from measured sequencing depth
  and has never been run; there is no BAM in this repository.
- **One synthetic pair for the size effect.** The 5-to-3 figure is from one
  80 kb fixture. The direction follows from the window widths and should hold
  generally; the magnitude is one number.

## The reusable part

The fixed-size measurement and the free-size measurement disagreed about
whether to ship, and only running both showed it. A panel-size-fixed comparison
cannot see a change that acts on panel SIZE, and `num_primers` being a request
rather than a guarantee is documented in CLAUDE.md -- so the second measurement
was available to be thought of, and was not thought of until a test fixture
whose premise depended on panel size went red.

That test, `test_the_repair_scan_is_bounded_in_production.py`, is what caught
it. Its comment reads "Tight enough that the optimizer's first panel misses it,
so the repair runs", and the first attempt at this change was to retune that
constant. Retuning would have buried the finding.
