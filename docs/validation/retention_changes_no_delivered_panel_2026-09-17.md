# Does retaining more candidates change a delivered panel?

Measured 2026-09-17 on the real Wolbachia pair. Phase 5 of the plan for
`docs/validation/pipeline_audit_2026-09-16/`.

**On this pair at panel size 12, no.** Three candidate universes spanning two
orders of magnitude reach the same selectivity floor, and the two larger ones
deliver identical panels while costing 19 and 130 times the search time.

The [retention benchmark](wolbachia_retention_benchmark_2026-09-16.md) stopped
after `filter`, so it measured what retention costs and could not say what it
buys. Phase 4 made the question answerable: the inventory is now reachable, the
frontier refills, and the scan is bounded.

Reproduce with `scripts/benchmarking/floor_boundary.py`.

## What was run

wMel `NC_002978.6` against the Drosophila `GCF_000001215.4` background, k = 12,
phi29 at 30 C. Two `filter` runs into fresh directories, one per retention mode,
both reproducing the documented funnel exactly: 874,596 distinct 12-mers,
491,836 after thermodynamic QC, 20,670 after the evenness gate, 2,000 after the
`max_primer` cut.

Then one `plan-pool` size row at panel 12 per universe and floor, with every
budget held equal across universes: 10,000 swap evaluations, 10 s, scan width
64. Only the universe and the floor vary. The shortlist rows pin the frontier at
2,000 and allow no refills; the others allow eight, which is enough to reach
either universe by doubling.

## Where each universe stops

| Floor | Shortlist, 2,000 | Post-Gini, 20,670 | All hard QC, 491,836 |
|---|---|---|---|
| 40 | feasible, density 40.01 | feasible, 0 refills | not run |
| 60 | feasible, density 60.11 | feasible, 0 refills | feasible, 0 refills |
| 80 | fails at 28.78 | fails at 19.16 | fails at 19.16 |
| 100 | fails at 28.78 | fails at 19.16 | not run |
| 140 | fails at 28.78 | fails at 19.16 | not run |

All three stop between 60 and 80. Retaining 10 times or 245 times more
candidates did not raise the floor this pair can reach.

At 40 and 60 the larger universes used no refills at all, because the shortlist
already qualified and a row that qualifies never refills. Those rows are
therefore identical to the shortlist's by construction, which is the refill
logic behaving correctly rather than a measurement.

## The two large universes deliver the same panel

At floor 80, where both fail:

| Quantity | Post-Gini | All hard QC |
|---|---|---|
| Candidates examined | 20,670 | 491,836 |
| Refills used | 4 | 8 |
| Ended | inventory exhausted | inventory exhausted |
| Selectivity density | 19.16197258152639 | 19.16197258152639 |
| Occupancy-weighted coverage | 0.7416316009542975 | 0.7416316009542975 |
| Seconds | 41.1 | 771.9 |

Two independent quantities agree to sixteen significant figures. The 471,166
candidates that only the `all_qc` universe holds changed nothing about the
delivered panel, and cost 772 s against 41 s, plus 123 s to build the position
cache and 785 MB of index against 20 MB.

## The larger universe reports a worse panel, and that is a known defect

At every floor it fails, the larger universe reports a LOWER density than the
shortlist: 19.16 against 28.78, while coverage rises from 0.7200 to 0.7416.

That is not retention's doing and this note should not be read as though it
were. It is the drift recorded in
[violation_magnitude_2026-09-17.md](violation_magnitude_2026-09-17.md): with a
failed repair now returning its input, the panel a rejected row reports is the
optimizer's own selection, and Stage 1 does not see the objective at all. More
candidates therefore let the coverage term climb at specificity's expense. Until
Stage 1 is constraint-aware, a bigger universe will report a worse panel on any
row it cannot satisfy.

So the negative result above has two components that this measurement cannot
separate: retention genuinely not helping, and a selection stage that spends
extra candidates on the wrong axis.

## What this changes about the default

Nothing yet, and deliberately. `candidate_retention` stays as it is, because
this measures one pair at one panel size and the honest conclusion is "no
benefit demonstrated here", not "no benefit exists". What it does establish is
that the cost is real and the benefit is not yet: anyone choosing `all_qc` on a
pair like this is paying 785 MB and an order of magnitude of search time for a
panel they could have had from 2,000 candidates.

## What this does not establish

- One pair, one panel size, one floor sweep, single runs, one machine. A
  different target, a larger panel or a background of a different size could all
  change where each universe stops.
- It does not test whether retention helps at panel sizes where the shortlist is
  genuinely too small. Size 12 out of 2,000 candidates is not a demanding ask,
  and the interesting case may be a 96 or 160-oligo panel.
- The `all_qc` sweep covers floors 60 and 80 only, because a full-universe row
  costs about 13 minutes. Floors 100 and 140 were run for the two smaller
  universes and extrapolated for the third on the strength of the identical
  panel at 80.
- It cannot separate retention's contribution from the Stage 1 defect above,
  which is the single most useful thing to fix before repeating this.
- Coverage and selectivity remain modelled site geometry, not measured
  amplification, so "reaches a floor of 60" is a statement about this model. See
  [evidence_matrix_2026-09-15.md](evidence_matrix_2026-09-15.md).
