# The fg/bg prefilter deleted candidates on a rule it did not apply

Measured 2026-09-17 on the real Wolbachia pair. Phase 4 increment 6 of the plan
for `docs/validation/pipeline_audit_2026-09-16/`, whose decision this confirms:
"The background prefilter becomes an ordering heuristic. It stops deleting
candidates and instead orders the scan."

The plan asked for the change to be gated on a before-and-after measurement.
This is it, and it turned up a second defect the plan had not named.

## `--min-fg-bg-ratio` was inert

`_prefilter_by_background` kept every candidate at or above the ratio, and then,
if that would remove more than `max_removal_fraction` of them, discarded the
threshold and kept the top 80% by ratio instead.

On the 2,000-candidate shortlist the threshold removes 64.8% at its default of
1.0, so the clause fired:

| `min_ratio` | Candidates below it | Clause fired | Actually removed |
|---|---|---|---|
| 1.0 | 1,295 of 2,000 | yes | exactly 400 |
| 2.0 | 1,700 | yes | exactly 400 |
| 5.0 | 1,926 | yes | exactly 400 |
| 20.0 | 1,991 | yes | exactly 400 |

The flag a user sets changed nothing above about 1.0. The rule in force was
"always drop the worst 20% by ratio", and `--min-fg-bg-ratio` is a documented
command-line option on `optimize`. That is the class Known Issue 8 tracks, in a
shape none of its ratchets look for: not a flag nobody reads, but a flag that is
read and then overruled by a second rule on the same decision.

## Which candidates survived depended on the batch

`max_removal_fraction` bounds the fraction of a batch, which is not a rule about
candidates. The same candidate could be kept in one pool and dropped in another
purely because of how many others happened to be below the threshold alongside
it. A test now pins that the decision is per candidate.

## What ordering does instead

Nothing is deleted. Candidates at or above the ratio are searched first, the
rest are searched last, and every one stays reachable.

| `min_ratio` | Deleting: pool | Deleting: dropped | Ordering: pool | Ordering: deprioritised |
|---|---|---|---|---|
| 1.0 | 1,600 | 400 | 2,000 | 1,295 |
| 20.0 | 1,600 | 400 | 2,000 | 1,991 |

The threshold now has an effect, and 400 candidates that the old path made
unreachable at every setting are reachable again. That matters more after
increment 5 than it would have before: a candidate at the back of the scan can
still be selected, because the frontier refills. Deleting one is final, which is
the shape of Known Issue 9, where an evenness gate removed primers the optimizer
had already selected.

The partition is stable, so the order within each group is the order the
candidates arrived in. That order is the inventory's `search_rank`, which
increment 5 established as the traversal, and the optimizers are
order-sensitive; a full re-sort by ratio would have discarded it.

## What it did to the delivered panel

One design at panel size 12, the deleting pool against the ordering pool:

| Quantity | Value |
|---|---|
| Panel size | 12 both ways |
| Primers shared | 11 of 12 |
| Jaccard | 0.846 |

One primer in twelve changed. `bg_max_removal` is retired with the clause; it
was a keyword argument with no command-line flag, no schema entry and no test.

## What this does not establish

- One pool, one pair, one panel size, single runs. The one-primer difference is
  what this configuration produced, not a bound on the change.
- It does not show the new panel is better. Neither panel was evaluated against
  a measured amplification, and the selectivity a design reports is modelled
  site geometry. See
  [evidence_matrix_2026-09-15.md](evidence_matrix_2026-09-15.md).
- It does not measure the cost of carrying 2,000 candidates rather than 1,600
  through selection. The two designs took 13 s and 8 s here, which is noise at
  this scale rather than a saving.
- The ordering is by fg/bg site-count ratio, the same quantity the deleting
  version used. Whether that is the right quantity to order on is untouched: it
  counts exact matches and carries no genome length, which is the concern Known
  Issue 6 records for `selectivity_ratio`.
