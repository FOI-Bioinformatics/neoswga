# Where the panel assessment and the post-optimization validator disagree

Measured 23 September 2026, before either record is allowed to gate more than it
already does.

Since `evaluate_panel` was wired into `run_optimization`, every run writes two
verdicts about the same panel into `step4_improved_df_validation.json`:

- `ok` -- the post-optimization validator, `result_validation.validate_result`
- `assessment.qualified` -- the panel assessment, `panel_evaluation.evaluate_panel`

They check different things and neither subsumes the other, so they are expected
to differ. Nothing had enumerated where. This does.

## Command

```bash
python scripts/benchmarking/assessment_vs_validator.py
```

## Result

| case | n | validator `ok` | `assessment.qualified` | assessment faulted |
|---|---|---|---|---|
| plasmid | 1 | True | True | -- |
| plasmid | 2 | True | False | size |
| plasmid | 4 | True | False | size |
| plasmid | 6 | True | False | size |
| plasmid + unreachable `min_selectivity_density` | 6 | False | False | configured limit, size |

Five runs, five artifacts, **three disagreements, and every one of them is the
size rule.**

## The control is what makes this a finding

The n=1 row is not padding. The shipped plasmid pool supports exactly one primer
at its configured `max_dimer_bp` of 3, so n=1 is the only size at which the
delivered panel MEETS the request -- and it is the only row where the two
records agree. Remove it and the table shows a correlation; with it, the cause
is isolated: **when the requested size is met, the two verdicts agree.**

The last row is the second control, from the other side. Both records say False
there, but for reasons that only partly overlap: the validator because
`panel_limit_not_met` now reaches it, the assessment because it faults the limit
AND the size. Agreement is not the same as agreement for the same reason, and a
count alone would have hidden that.

## Why the size rule is the whole of it

The two compare differently, and this was recorded in
`panel_evaluation._panel_violations` when PR #92 established that the three size
checks in this codebase are not duplicates:

| site | comparison | level |
|---|---|---|
| `result_validation.py` | `n != target_size` | `warning` when the status is already PARTIAL, `error` otherwise |
| `panel_evaluation.py` | `len(primers) < requested` | always a violation |

On every plasmid run above the optimizer returned PARTIAL, so the validator
recorded `set_size_mismatch` at **warning** level and `ok` stayed True. The
assessment has no levels: a violation is a violation, and `qualified` is false
exactly when the list is non-empty.

Neither is wrong. They answer different questions. A short panel is a normal and
documented outcome in this codebase -- `num_primers` is a request, not a
guarantee, and selection stops rather than admitting a pair above
`max_dimer_bp`. The validator treats that as an expected degeneracy; the
assessment treats it as not meeting the request, which it did not.

## What this means for letting the assessment gate

The plan for that increment fixed its decision rule in advance: **zero
disagreements, consolidate; any disagreements, each is its own decision.** There
are disagreements, so the answer is the second branch, and there is exactly one
decision to make.

**Recommendation: do not let `assessment.qualified` gate the export as it
stands.** Every design whose panel is shorter than requested would begin
refusing, and on this fixture that is three runs out of four. That is not a
latent defect being surfaced; it is a different question being asked.

Two ways forward, and this measurement does not choose between them:

1. **Gate on a subset.** Exclude the size violation from what blocks, the way
   PR #99 scoped `panel_limit_not_met` to the configured limits alone. Cheap,
   and it keeps the assessment's other checks -- non-finite metrics, the
   delivered-dimer breach, the configured limits -- available to the gate.
2. **Reconcile the size rule first.** Give the assessment the validator's
   notion that a short panel from a PARTIAL run is expected. That is a change
   to what `qualified` MEANS, and it would make the assessment agree with the
   validator by construction rather than by measurement, which is worth less.

## What this does not establish

**One fixture.** The plasmid example is 5.4 kb and delivers one to a few
primers. A larger design that fills its requested size would agree on the size
rule and might disagree elsewhere -- on a non-finite metric, or on a configured
limit the repair resolved. Nothing here rules that out; it was not reachable.

**Wolbachia was attempted and is not runnable.** `optimize` refuses before any
panel exists:

```
ReferenceDataError: Reference data unusable (position index):
  .../drosophila_12mer_positions.h5: no recorded reference digest to compare;
  .../wmel_12mer_positions.h5: no recorded reference digest to compare
```

**This corrects a claim in CLAUDE.md.** Under **Two scanners, one quantity**, the
consequence for the shipped example is given as: *"wMel is one record, so its
index is accepted. Drosophila has 1,870 and is refused until regenerated."* The
actual refusal is earlier and wider. It comes from the reference-digest check,
not the record-geometry check, and it names **both** indexes including wMel. The
prose was right about the direction and wrong about the mechanism and the scope.
Regenerating either index needs a fresh `count-kmers` plus `filter`, which for
the 144 Mb host is the reason this was not done here.

**The disagreement is measured on the verdicts, not on delivered panels.** No
panel moved in producing this table; both records describe the same search
output.
