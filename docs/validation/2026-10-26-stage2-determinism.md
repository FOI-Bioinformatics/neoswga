# Making hybrid's Stage 2 deterministic

Measured 26 October 2026, on `main` at `a526687` plus the one-line change this
document records. Reproduce with:

```
PYTHONPATH=. python docs/validation/stage2_determinism_2026-10-26/measure.py
PYTHONPATH=. python docs/validation/stage2_determinism_2026-10-26/measure.py --unseeded
```

## The defect

`hybrid` is the default method, so it produces most of the panels this tool
ships, and it was not reproducible. Running `optimize` twice on the same input
could deliver different primer sets.

The two preceding records
([wiring the Tm term](2026-10-26-wiring-the-tm-term.md) and the floor change)
both had to measure around this: a single run told you nothing, so every figure
there is a median over repeats or a comparison isolated from the loop. That
work named the nondeterminism and deferred fixing it, because a fix changes
delivered panels for everyone and is its own decision. This is that fix.

## The one source that mattered

Stage-2 network refinement (`_network_refine`) drops one primer at a time,
scoring each candidate removal on the algebraic connectivity of the network left
behind. `nx.algebraic_connectivity(subgraph, method="tracemin_lu")` starts its
eigenvalue iteration from a RANDOM vector, drawn from numpy's global state with
no seed, advancing per call. So the connectivity figure the removal ranking read
wobbled between runs, and a wobble large enough to flip a near-tie changed which
primer was dropped, and the delivered panel with it.

Two candidates were considered and only one was real:

- **The eigensolver's random start.** Real. Measured below.
- **Hash-ordered set iteration in the removal loop.** Not a factor. The loop
  iterates `current_primers`, a list in Stage-1 order, and breaks ranking ties
  by that order. The one set it derives, `remaining_nodes`, feeds
  `graph.subgraph(...)`, whose node order follows the graph's insertion order,
  not the set's; and `_build_network` inserts nodes in primer-then-position
  order, which is deterministic. The Laplacian eigenvalue is invariant to node
  permutation anyway. So no set needed sorting, and the proposed "order the loop's
  set iterations" turned out to be unnecessary once the seed was fixed. Checked,
  not assumed: with the solver seeded, the delivered panel is identical across
  three `PYTHONHASHSEED` values (below), which it could not be if set order
  still fed a decision.

## Measurement

The mixed-length plasmid pool (k=9-12, the pool that wobbled most in the Tm-term
record), one fixed input, distinct delivered panels over 6 repeats in one
process:

| | distinct panels / 6 |
|---|---|
| unseeded (pre-fix) | **2** |
| seeded (post-fix) | **1** |

Across hash seeds, seeded, the delivered panel is not just stable but identical:
the same signature `7fe3bf28` at `PYTHONHASHSEED` 0, 1 and 2, over 6 runs each.
Unseeded the same harness gave 2-3 distinct panels depending on the seed. The
`--unseeded` flag on the script reproduces the pre-fix behaviour by forcing
`seed=None` back onto the solver, so the defect can be seen without checking out
the old commit.

## The fix

One module constant and a keyword on two calls:

```python
_CONNECTIVITY_SEED = 0
...
nx.algebraic_connectivity(subgraph, method="tracemin_lu", seed=_CONNECTIVITY_SEED)
```

The second-smallest Laplacian eigenvalue is a property of the graph, not of the
vector the iteration starts from, so pinning the vector removes the numerical
wobble without changing what connectivity means. A fixed internal constant, not
a configurable one: reproducibility is not a knob a user should have to set, and
the run manifest already records a seed the tool then did not honour here.

## What it changes, and why that is the point

**It changes delivered panels.** A seeded solver walks a different numerical
path from whatever unseeded path a given run happened to take, so the panel a
user gets today can differ from the one they got yesterday. That is the
one-time cost of becoming reproducible, and it is why this was a decision rather
than a repair. From here on the same input gives the same panel, on a re-run, on
another machine, and under `pytest -n 8`.

The full suite is unaffected in pass/fail terms: 6466 passed, 282 skipped, 55
failed, 4 errors, with the failing IDs identical to the pre-change baseline (all
environmental -- no k-mer counter, no reference genomes). The existing hybrid
tests assert invariants (coverage bounds, requested size, no duplicates, the
Stage-1 floor) rather than exact membership, so a shifted-but-still-valid panel
passes them. `tests/test_hybrid_optimizer_run.py` gains two guards: a source
assertion that every `algebraic_connectivity` call passes a seed, and a
behavioural two-run equality check.

## Size of the claim

Determinism itself is now exact and is pinned by a test, not a measurement: one
input gives one panel, full stop. The quantified *amount* of prior
nondeterminism (2-3 panels) is one pool, one panel size, one reaction, as in the
sibling records. What is not claimed is that the seeded panel is BETTER than any
particular unseeded one -- the eigenvalue is the same quantity either way, so
across inputs this neither helps nor hurts panel quality; it only makes the
answer repeatable.

## Not done

The other methods are out of scope. `dominating-set` has no connectivity
eigensolve. `network` and `clique` were not audited for their own
nondeterminism here; if a run of either proves irreproducible that is a separate
record. This change is confined to `hybrid`'s Stage 2.
