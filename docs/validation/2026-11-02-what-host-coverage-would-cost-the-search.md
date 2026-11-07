# A host-coverage field costs +48% per evaluation, and halves the path already in use

Measured 2 November 2026, on `main` at `112159b`. Phase 2.1 of
[the plan](../reports/PLAN_2026-11-02_known_issues.md), against Known Issue 14.
Reproduce with:

```
PYTHONPATH=. python docs/validation/host_coverage_cost_2026-11-02/measure.py \
    --params <work>/params.json --out results.json
```

Results in [`host_coverage_cost_2026-11-02/results.json`](host_coverage_cost_2026-11-02/results.json),
the pool's configuration beside it. **No code changed.** This is the
measurement the decision needs, not the decision.

## The question, restated

`stage2_background_weight` is REFUSED above zero because the weighted host axis
needs a host-coverage fraction and `PoolMetrics` carries five fields and no
sixth ([`2026-10-29-stage2-background-weight-does-not-run.md`](2026-10-29-stage2-background-weight-does-not-run.md)).
The stated reason for five is cost: `plan_pool` spent 29.7 s, of which 27.6 s
was inside `compute_metrics` across 4,040 calls.

So the open question in Known Issue 14 is not whether a weighted host term is a
good idea. It is **what carrying host coverage through the search would cost**,
and that figure did not exist — the entry asserted the cost without a number.

**The alternative is not "nothing".** `panel_refinement` line 47 already falls
back to the full `compute_metrics` whenever `max_mean_gap`, `max_evenness` or
`max_host_coverage` is set:

```python
extended_limits = (max_mean_gap, max_evenness, max_host_coverage)
evaluate = optimizer.compute_metrics
if all(value is None for value in extended_limits):
    evaluate = getattr(optimizer, "compute_pool_metrics", None) or evaluate
```

So a host-coverage limit today already pays the expensive path. That makes this
a three-way comparison.

## The measurement

40 random 12-primer panels from a 300-candidate pool, equiphi29 at 43 C,
1.27 Mb target at 35% GC, 4 Mb host at 42% GC. Median per call, best of three
repeats each:

| arm | what it computes | median | against A |
|---|---|---|---|
| **A** | `compute_pool_metrics`, as shipped | 3.56 ms | — |
| **B** | A + one host-coverage call | **5.26 ms** | **1.48x** |
| **C** | `compute_metrics`, today's fallback | 10.77 ms | 3.03x |

**C is 2.05x B.** Projected over the 4,040 calls the five-field decision was
made against: A 14.4 s, B 21.3 s, C 43.5 s.

Re-running gave A 3.598/3.556 ms and B/A 1.47/1.48, so the figures are stable
to about 1%.

## Why the marginal cost is small: the positions are already in hand

This is the part that changes the answer. `compute_pool_metrics` **already
gathers and sorts every background position** — it needs them for
`total_bg_sites` — and then discards them:

```python
bg_by_prefix = {prefix: set() for prefix in optimizer.bg_prefixes}
...
bg_by_prefix = {p: sorted(v) for p, v in bg_by_prefix.items()}
...
total_bg = sum(len(bg_by_prefix[p]) for p in optimizer.bg_prefixes)
```

So the marginal work is **one `_coverage_over_prefixes` call over a list
already built**, not the position lookups — which are the part that would have
been expensive, and the part the five-field decision was protecting against.
The marginal is timed here on positions gathered OUTSIDE the timed region,
because that gathering is already inside arm A and counting it twice would
overstate B.

Measured marginal: 1.71 ms, which is 32% of B.

## The caveat that bounds this: it scales with the host's LENGTH

`union_coverage` allocates `np.zeros(total_length, dtype=bool)`, so the cost is
set by how long the host is and not by how many sites the panel has in it. At a
fixed 400 sites:

| host length | marginal call |
|---|---|
| 1 Mb | 0.68 ms |
| 4 Mb | 2.10 ms |
| 16 Mb | 8.09 ms |
| 64 Mb | 33.50 ms |

Linear, within measurement noise. Extrapolated to hg38 at 3.3 Gb that is
roughly **1.7 s per evaluation and a 3.3 GB boolean array**, so about two hours
over 4,040 calls.

**That is not an argument against the field.** Arm C already does exactly this
allocation, so the full-metrics path is equally infeasible at host scale — and
it is the path a `max_host_coverage` limit takes today. The ceiling belongs to
`union_coverage`'s representation, which is Known Issues 1 and 3's territory,
not to the choice between A, B and C.

## What this implies for the decision (Phase 2.2)

Reported, not acted on. The plan's two branches resolve as:

- **Affordable at the scale the search runs at.** +48% per evaluation for a
  quantity the search currently cannot see at all, on a bacterial host.
- **Strictly better than the status quo for host-limited runs.** Any design
  setting `max_host_coverage`, `max_mean_gap` or `max_evenness` is on arm C
  today. Giving `PoolMetrics` the field would let those runs use B instead —
  about **half** the per-evaluation cost — which is a speedup, not a tax.

So the cost objection recorded in `pool_metrics.py` holds for the position
lookups it was written about and does not hold for this field, because the
positions are already gathered. What the figure does not settle is whether a
weighted host axis produces a better panel than the lexicographic tie-break;
that is the measurement
[`2026-10-29-stage2-background-weight-does-not-run.md`](2026-10-29-stage2-background-weight-does-not-run.md)
could not make, and it needs the field to exist before it can be made.

## Scope of the claim

One synthetic pool, one panel size, one host size for the three arms, one
reaction, 40 panels. The host is synthetic, so its site geometry is Poisson
where a real genome clusters — which affects how many spans are marked but not
the array allocation that dominates the marginal cost. No real reference genome
is reachable from this container (`tests/validation/genomes/` is gitignored and
absent, NCBI is refused by the egress proxy), so the hg38 figure above is an
extrapolation from the 64 Mb point and is labelled as one. Nothing here is a
wet-lab outcome, and no delivered panel changed: no code was modified.
