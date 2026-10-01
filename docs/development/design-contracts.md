# Design contracts

> Moved verbatim from `CLAUDE.md` on 2026-10-01, when that file was reduced to a
> working summary. Entries are dated records: a figure is what was
> measured on the date given, on the references named, and a later entry may
> correct an earlier one. `CLAUDE.md` keeps the short form and links here.

The rules a design run enforces, each with the failure it was written for:
what fails a run, how a request is resolved, which candidate pool is
searched, how positions and coverage are checked against independent
oracles, what the search budget does and does not bound, and how sequencing
data is read. "Known Issue N" refers to [KNOWN_ISSUES.md](KNOWN_ISSUES.md).

## The design-failure contract (2026-09-21)

A required calculation that fails now fails the run. It does not return a
substitute value. Four errors in `core/exceptions.py` carry this, all under
`DesignError`: `InvalidDesignRequest`, `ReferenceDataError`,
`UnsupportedModelError` and `ModelEvaluationError`. `SearchBudgetExhausted` is
deliberately outside the family, because spending an allowance is a recorded
stopping point rather than a failure.

The distinction the family draws is between a measurement and its absence. A
candidate that misses a Tm window has been measured and rejected; a candidate
whose Tm raised has not been measured at all. `DesignError.qc_reason` is always
None, so code asking "was this a QC rejection" gets a definite no.

What changed, and what each substitution used to cost:

| Site | Was | Now |
|---|---|---|
| `thermodynamics.calculate_tm_batch` | NaN for any failure | raises; a non-ACGT base is a named `InvalidSequenceError` with `qc_reason`, which is a QC rejection and stays one |
| `thermodynamic_filter._check_heterodimer_pair` | 0.0 free energy, the most permissive answer the screen has | raises; a positive or infinite duplex energy still returns 0.0, because that is a measurement |
| `coverage.polymerase_extension_reach`, `product_reach` | a default reach for an unknown polymerase | raises `UnsupportedModelError` naming the supported set |
| `coverage._record_starts_for` | None when the getter raised | raises; a cache with NO getter still returns None, which is absence rather than failure |
| `occupancy.discrimination_profile` | skipped a failed primer | raises; a mean over an unknown subset was reported with the authority of a mean over the pool |
| `unified_optimizer` per-target coverage | empty dict | raises; `base_optimizer` gates the floor on a non-empty dict, so a requested `min_per_target_coverage` passed vacuously |
| `unified_optimizer` application weights | debug line, defaults applied | raises; `--application clinical` silently had no effect |
| `unified_optimizer` ensemble winner evaluator | debug line, `optimizer=None` | raises; with None every configured panel limit went unenforced |
| `unified_optimizer` post-optimization validator | skipped | raises; skipping disarms the duplicate, size-drift, zero-coverage, blacklist and delivered-dimer checks at once, and writes no validation file, so `export` prints "ready for ordering" |
| `_reseed` | `pass` | raises; the caller logged "set for reproducibility" either way |
| ensemble member failure | any exception became an `status: "error"` row | a `DesignError` propagates, because the next member computes the same quantity from the same data; an algorithm that cannot run on this pool is still a visible row |

`core/design_result.py` separates three things that were one. **Run state** is
what happened to the process (`finished`, `failed`, `interrupted`). **Termination
reason** is why the search stopped (`qualified`, `budget_exhausted`,
`candidates_exhausted`, `refill_exhausted`, `error`). **Qualification** is a
property of the panel. `recommendation_allowed(run_state, qualified)` needs both,
and refuses an unknown state rather than defaulting either way.

At the command boundary a `DesignError` prints the stage, field, artifact, model
and input, then writes `design_failure.json` into the run directory and exits
nonzero. The record exists because an output directory holding last week's
`step4_improved_df.csv` reads exactly like one holding this morning's. Each
pipeline step re-raises `DesignError` rather than reducing it to "step N
failed"; `cli/_failure.py` owns the record.

## One resolved design request

`core/design_request.py` resolves a params mapping into a frozen `DesignRequest`
carrying references, chemistry, candidate-source identity, fixed and excluded
oligos, panel limits, size policy, search budgets, seed, model identifiers and
the concentration policy. Nested content is tuples, so a stage cannot append to
a list it was handed.

`default_sources` records, per setting, whether the request supplied it or which
default did. `request_hash` is a SHA-256 over a canonical JSON form, so it is
stable across processes and independent of key order; it is recorded in the run
manifest for `optimize`.

`optimize` resolves the request from the params FILE before the search starts,
not from the `parameter` module: `get_params` runs inside `optimize_step4`, so
at that point every reaction global still holds its default. This is the same
ordering trap `warn_on_condition_drift` documents.

Refusals it makes that used to be silent: unknown and retired keys (a leading
underscore marks a comment and is accepted), non-finite values, `coverage_reach`
of 0, negative budgets, an oligo that is both fixed and excluded, a
background-measured panel limit with no background genome, and an unsupported
polymerase. The explicit zero matters on its own: `design_context_from_params`
used `override or params.get("coverage_reach")`, and 0 is falsy, so it silently
became 3 kb and every coverage figure was reported at a reach the request did
not ask for.

All three design commands resolve it from the same file: `optimize`,
`plan-pool` and `expand-primers`. A command that skipped the gate would accept
what the others refuse, which is how one params file came to mean different
chemistry depending on which command was run.
`tests/test_resolved_design_request.py` walks the call path from each handler
rather than checking that a call appears somewhere in the module.

**Not yet done from the plan's Task 2**: evaluator code still reads `parameter`
globals at run time, and `OptimizationRequest.optimizer` still owns the
scientific settings. The request is a validation gate and a provenance record,
not yet the single channel those settings travel through.

One instance of the mutable-global read was found and removed the hard way. An
index-identity check placed inside `run_optimization` read
`parameter.fg_genomes` and paired it with the prefixes the CALL was given;
under `pytest -n 8` that paired a test's own prefix with another test's FASTA
and the design refused its own index. Reference identity now travels on the
request, where a prefix and its genome are named together.

## Which candidate pool a command searches

`open_source_or_list` is gone, replaced by three functions in
`candidate_source.py` that keep absence and failure apart.

- `open_explicit_source(candidates)`: the pool the user named. It always wins.
- `open_inventory_source(...)`: the inventory, or **None** when the directory has
  none, which is a fact about the directory. It raises `ReferenceDataError` when
  the directory HAS an inventory that holds nothing under this reaction.
- `open_design_source(...)`: the rule every command asks, built from those two.

The old function wrapped the inventory open in `except ValueError` and fell back
to the caller's CSV for both cases. A reaction fingerprint mismatch therefore
became a quiet run over the `max_primer` shortlist, at the shortlist's frontier,
with everything the inventory held unreachable and one `logger.info` line to say
so. That is how the occupancy-gate measurement in Known Issue 17 produced an
apparent density improvement that was not real.

So `filter --preset enhanced_equiphi29` followed by a plain `optimize` now
refuses, naming the remedy, where it used to warn and proceed.
`tests/test_optimize_warns_on_condition_drift.py` was inverted to match.

A self-dimer screen that empties a non-empty frontier now raises
`NoCandidatesError` naming the screen and its threshold, instead of handing an
empty list to an optimizer that answered "candidates list cannot be empty".

## Two scanners, one quantity

Fixed 2026-09-21. `string_search` has two position scanners: the
Aho-Corasick `get_all_positions_multi_k`, and the sliding-window
`get_all_positions_per_k` used when that package is absent or one k is being
scanned. Records are concatenated with no separator, so the last k-1 bases of
one record and the first bases of the next form k-mers that occur in neither.
The Aho-Corasick path rejected those matches. The sliding-window path did not.

So the same reference produced different site sets depending on which path
ran, and the sliding-window path stored up to k-1 fabricated sites per record
join. On a two-record fixture `ACGGTA` is absent from both records and was
stored at offset 4. A fabricated foreground site inflates coverage, a
fabricated background site deflates specificity, and nothing downstream can
tell either from a real one.

Single-record references are unaffected: with no joins the two scanners always
agreed, which is why every complete bacterial genome and every plasmid in this
repository reads the same before and after. Draft assemblies and hg38 are where
it bit.

`spans_a_record_join` now holds the rule once and both scanners call it. A
match beginning exactly ON a boundary starts a record and is kept; off by one
here would delete the first k-1 sites of every contig.

`INDEX_FORMAT_VERSION` is 2. A version 1 index is refused for a MULTI-RECORD
reference only, because a single-record one cannot carry the defect and
refusing it would force a recount for something that never applied to it.

**Missing record GEOMETRY is judged the same way, and by a different layer.**
An index with no `#record_starts` lets a coverage window run past a contig
edge into the next record. Whether that matters depends on how many records
the reference holds, and only the resolved request pairs a prefix with a
genome, so `reference_check.verify_index_geometry` decides it against the
manifest rather than `PositionCache` deciding it alone. Deciding it in the
evaluator means reading `parameter.fg_genomes` and pairing it with the
prefixes the call was GIVEN, which is the defect that made a design refuse
its own index under `pytest -n 8`. Record counting reads header lines; the
genome loader would hold 8.5 GB for hg38 to answer it.

Measured on the shipped Wolbachia design: the 12-oligo panel's `bg_coverage`
against *Drosophila* reads 0.00217080 unconfined against 0.00207162 confined,
an inflation of 14,255 bp or **+4.788% relative**, from 52 host sites across
1,870 records. The error overstates host coverage, so it is not flattering,
but `max_host_coverage` is a configurable limit and a panel could be rejected
for coverage it does not have.

The same measurement on Prevotella, two chromosomes and one join, with 724
target sites, gives **exactly zero**: the region either side of the join is
already covered from both directions. So the magnitude is joins times
sparsity, and neither figure generalises alone
([measurement](../validation/record_geometry_on_drosophila_2026-09-21.md)).

**Consequence for the shipped example.** Both its indexes predate record
geometry. wMel is one record, so its index is accepted. *Drosophila* has 1,870
and is refused until regenerated. `plan-pool` already refused both before this
work; what changed is that `optimize` applies the same standard.

`tests/test_positions_agree_with_an_independent_count.py` checks the scan
against a brute-force sliding window written in that file, which calls neither
scanner nor any helper they call. It covers overlapping occurrences,
palindromes, a reverse complement that also occurs forward, ambiguous bases, a
circular origin, record joins, and an exhaustive pass over every window of a
small reference. It also asserts the two scanners agree, which is the check
that would have caught this.

## What the chemistry model supports

`neoswga/core/registry/model_evidence.json` records, per constant, what it is
and over what domain the model supports it. `model_evidence.py` loads it and
`require_model_support(request)` runs inside `resolve_design_request`, so a
computation outside a recorded domain is refused before any index is opened.

Five statuses, and the line that carries the weight is between the first two.

| Status | Meaning |
|---|---|
| `measured` | the cited work reports this value for a case the model applies it to |
| `estimated` | extrapolated from data at another temperature, on longer DNA, or in another buffer |
| `empirical` | chosen so the model behaves plausibly; no source reports it |
| `assumed` | a modelling decision with a stated reason and no measurement |
| `absent` | nothing computes this effect, and the code must not report zero for it |

Most additive coefficients are 37 C figures for PCR-length duplexes applied to
12-mers at 30 C. That may well be fine; it is not a measurement of it. Only
`tm_urea` was chosen because its source concerns short oligos.

**The registry does not claim the literature was re-read.** Every `source` is
the attribution the repository already carried. What the registry adds is the
second judgement the prose ledger never made: whether the cited work covers
this case. A registry that implied verification would break, in the act of
recording it, the rule it exists to enforce.

**What is refused**: an unknown polymerase; an oligo length outside the
enzyme's modelled range, which is the defect where a Bst design was filtered
through phi29's 6-12 bp window; and an additive whose duplex effect nothing
computes while the literature expects one. Estimates are NOT refused. A model
that declined to run on an extrapolated coefficient would decline to run.

Two findings from compiling it, both in
[docs/validation/chemistry_model_evidence.md](../validation/chemistry_model_evidence.md):

- **Glycerol is accepted, range-validated, printed, and changes no Tm.**
  Measured here: a 12-mer at 10% glycerol returns the same effective Tm as at
  0%, to the last digit. The literature expects a real destabilisation, and the
  shipped `q_solution` preset sets 10%. No coefficient is invented to close
  this, because inventing one is the promotion of an assumption the registry
  exists to prevent; a design that sets glycerol is refused instead. BSA and
  PEG also have no Tm term and are `assumed` rather than `absent`: they act on
  the enzyme and on crowding, not on duplex stability.
- **`neoswga.core.registry` was not installed.** It was missing from
  `pyproject.toml`'s explicit `packages` list, so a built wheel contained none
  of it, while `core/parameter.py` imports `registry.views` at module scope.
  Verified by building a wheel. Package data could not have helped:
  `include-package-data` applies to packages that are being installed. A test
  now compares the declared list against the packages on disk.

`docs/SCIENCE_CITATIONS.md` had drifted: it states Klenow processivity as
10,000 bp citing Bambara (1978) while the shipped registry says 40 bp. The
prose was right when written and the code moved. The registry is checked
against `registry/views.as_characteristics()` by a test, so the two cannot
disagree silently; read the prose document as commentary rather than as the
record.

## One panel assessment, and a coverage oracle that is not the code

`core/panel_evaluation.py` holds `evaluate_panel(request, oligos, metrics) ->
PanelAssessment`: one immutable record carrying the panel, the request hash,
every metric WITH its units and the reach and denominator it was computed at,
per-target results, every hard-constraint violation, and qualification as a
boolean that is true exactly when there are none.

Three rules it enforces, each for a failure this project has seen the shape of:

- **A non-finite required quantity fails the run.** NaN compares False against
  every threshold, so a panel carrying one passes no limit and fails no limit.
- **Unavailable and zero have different representations.** A panel that binds
  the host nowhere and a panel whose host index was never opened both read zero
  otherwise. `Measurement` refuses to hold both a value and a reason for not
  having one.
- **A verified zero background is a zero denominator, not a ratio.**
  `base_optimizer` reports `MAX_SELECTIVITY` (1e6) there, deliberately, because
  it is finite and JSON carries it; its own docstring says that means "no
  background binding was detected", not "measured this well", which concedes a
  reader cannot tell them apart from the number. The assessment says undefined
  and why, and keeps the site count. The sentinel is left in place because
  changing it moves every saved summary.

**The REPORT reads it as of 2026-09-27**, so a rendered figure and a stored one
cannot come from different arithmetic. `report/metrics.py` applies the
assessment LAST, after its own estimate and after the summary JSON, and a
quantity the assessment marks unavailable CLEARS the rendered field rather than
falling back to an estimate -- the report's own coverage guess is primers times
30 kb over genome length, which saturates at 1.0 on a small target.

What that caught: the summary carries `selectivity_ratio: 1000000.0` beside
`total_bg_sites: 0`, the MAX_SELECTIVITY sentinel, and the report rendered the
million as a specificity. It now renders nothing there and says why.

`qualified` still gates nothing, and the reason is now MEASURED rather than
feared. It is the absence of every violation including a panel shorter than
requested, and `num_primers` is a request. On the bundled plasmid example at
requested sizes 6 and 40 the delivered panel was short both times, the
validator said ok with a warning, and `qualified` was False -- so a gate on it
would refuse two ordinary runs. Violations are tagged BLOCKING or ADVISORY at
the site that raises them, and `acceptable` is what a gate may consult.

The optimizer and the acceptance path still assemble their own answers, which
is the rest of Task 5.

`tests/test_coverage_independent_oracle.py` checks coverage against a
base-by-base oracle written in that file. It calls neither
`merged_window_intervals` nor `_mark_window` nor anything they call, so
agreement is evidence rather than a restatement. It is deliberately the slow
implementation production replaced: the fast one accumulates log(1 - theta) at
window edges, and an edge-accounting error is invisible from inside that
formulation. Verified load-bearing by mutating the production grouping from
per-primer to per-site, which fails four of its cases.

The window convention, measured rather than assumed: a site at `pos` with reach
`r` covers `[pos - r, pos + r)`, so a window is `2r` wide.

**The two production coverage paths disagree across a record join, measured and
not fixed.** `compute_per_prefix_coverage` marks through `_mark_window` WITH
record starts, so a window stops at a contig edge. `_union_coverage` and the
occupancy path go through `merged_window_intervals`, which takes no record
starts by design. On the fixture in that file a site 2 bases before a join
covers 12 bases confined and 20 unconfined. So a multi-record reference has two
coverage figures and which one a reader sees depends on the code path. The test
asserts both numbers so the gap cannot grow unnoticed.

## The report and the saved result must agree

`tests/test_report_agrees_with_the_saved_result.py` asserts that every quantity
the report renders equals the one the summary holds, for the exact exported
panel. The report computes its own estimates from the results CSV and overrides
them with the optimizer summary where the summary has an opinion, which is the
right order; `from_optimizer` says which a reader is looking at.

**Writing that check found a favourable default.** `mean_gap`, `max_gap`,
`gap_gini` and `gap_entropy` were read with `.get(key, 0.0)`, and zero is the
BEST value for every one of them. A `max_gap` of 0.0 says the panel leaves no
coverage hole anywhere. A summary that did not carry the key therefore rendered
as the best possible measurement rather than as no measurement, and directories
written before these keys existed are explicitly supported, so this was
reachable rather than theoretical.

The four fields are now `Optional[float] = None`, read without a literal
default, and both render sites skip the gap section unless every figure in it
was measured. Rendering it with one missing is what put a "no coverage hole"
verdict beside three real numbers.

This is the silent-zero family again, in its fourth shape: not a scan that
found nothing, an integer that saturated, a cache asked for what it does not
hold, or a guard with the wrong predicate, but a dictionary lookup whose
fallback happens to be the answer everyone wants.

## Two counters, and only one of them bounds the run

`SearchBudget` is the SHARED ledger. It counts uncached evaluations of the
shared objective across every stage and raises when spent, so a later stage
cannot get a fresh allowance. `swap_max_evaluations` is a PER-STAGE allowance
inside the deletion and swap loops, and it starts at zero every time one of
them is entered.

The per-stage one is not a defect; bounding one loop is reasonable. What would
be a defect is believing it bounds the run. It defaults to 10,000 and looks
like a total. The only setting that is a total is `total_search_evaluations`,
and it is None by default, so **by default there is no total bound at all.**

`describe()` now carries `uncounted_scopes`, naming the two kinds of work the
ledger does not see, rather than leaving a reader to infer them from a count
lower than they expected:

- **final assessment**: one `compute_metrics` per stage once the panel is
  decided, deliberately uncharged so reporting cannot consume a search's
  allowance.

**Proposal generation left that list on 2026-09-27.** `clique` scored the top
`max_scored_sets` dimer-free sets through `compute_metrics` directly, which the
allowance could not see, so a user setting `total_search_evaluations` would find
the run spending past it. `run_panel_search` now attaches the ledger through
`attach_search_config` and that loop consumes it; exhaustion keeps the best set
scored so far rather than failing, because spending an allowance is a recorded
stopping point. `max_scored_sets` remains, bounding the method's own cost.

The ratchet's DETECTOR was refined rather than its allowlist reworded: a loop
that evaluates panels is accepted only when the spend is INSIDE that loop, so
one `consume()` elsewhere in the function cannot excuse an unmetered scan. Four
tests drive it on source written in the test file, because a refinement that
turned the ratchet off would otherwise be invisible.

**Alternatives escaped the ledger entirely until 2026-09-27.**
`budgeted_objective` restores the previous binding when it exits, correctly for
a context manager, so once `run_panel_search` returned the evaluator was
unbound and `collect_alternative_sets` searched each alternative with no budget
at all. A run declaring `total_search_evaluations` could spend past it once per
alternative, and that function's own `except SearchBudgetExhausted` clause was
a handler for something that could not happen -- which is the tell, since it
was written believing the ledger was in force. It now takes the run's budget
and binds it around each attempt.

**The binding does not yet bound anything, and the commit that added it
overstated the case.** It said alternatives then "spend the run's allowance
instead of none". Measured afterwards: an alternative search calls
`optimizer.optimize` and nothing else, and neither `dominating-set` nor
`hybrid` evaluates the shared objective inside `optimize` -- 0 objective
evaluations for both with a 5,000 allowance. So the hole is closed and nothing
in production travels through it, and that `except` clause is still unreachable
on the shipped methods. Routing alternatives through `run_panel_search` earned the
claim on 2026-09-28: the same five alternatives cost 0 counted evaluations
through the bare `optimize` and 2,988 through the contract. Panel quality is a
wash at about 5x the wall clock, paid only when `max_sets` exceeds 1
([measurement](../validation/alternatives_through_the_contract_2026-09-28.md)).

**That measurement found a separate defect, recorded and not fixed.** An
alternative violating a configured limit is offered either way -- set 4 sat at
density 17.05 against a floor of 20 -- and `export_is_blocked` takes a directory
with no set index, so the findings it reads describe set 0 while `export --set N`
delivers a different panel. The primary is held to its limits and the sets after
it are never assessed, which is Known Issue 19's family from the other side.

**Half of that is fixed as of 2026-09-28.** `export_is_blocked` takes the set
index and REFUSES a non-zero set, because nothing evaluated it: the findings and
the assessment in the directory describe set 0. Unknown is not success, and here
the unknown is total. `--allow-unqualified` is the deliberate override, as it is
for every other refusal there. Verified end to end on the plasmid example with
five sets: set 0 exports and set 1 is refused with no files written.

What is NOT fixed is the other half -- an alternative that violates a configured
limit is still OFFERED. Suppressing it, marking it, or assessing every set are
three different answers about what `max_sets` is for, and that is a decision
rather than a repair.
`tests/test_the_allowance_holds_on_the_production_path.py` fails when
`optimize` starts consulting the objective, so the claim is corrected in the
same change that makes it true. This is "verify the deciding stage" again: a
budget reaching the code is not a budget that bounds the answer.

`tests/test_search_budget_contract.py` holds the ratchet.
`UNCOUNTED_SEARCH_LOOPS` lists every function that evaluates panels in a loop
outside the objective, with its reason and the bound that does apply, and the
list can only shrink. **It is empty as of 2026-09-27.** A call made ONCE per stage is final assessment and is
not flagged; one inside a `for` or `while` is search work, and search work the
ledger cannot see is what the check is for. Verified load-bearing by wrapping
an existing single call in a loop, which fails it.

## The smallest pool, and what deletion cannot reach

`--minimize-primers` reduces a panel by removing one primer at a time and
keeping it when no single removal still qualifies. That is a local optimum, not
the smallest pool, and the two differ whenever one candidate covers what two
others cover between them.

`tests/test_smallest_pool_search.py` enumerates every subset of a four-candidate
set-cover fixture and compares the search against the enumerated answer. The
fixture is constructed so the structure is visible: one candidate covers the
union of two others, so `{P1,P2,P3}` qualifies at size 3, no single deletion
from it qualifies, and `{P4,P3}` qualifies at size 2. Deletion stops at 3;
`panel_beam.beam_search` asked for 2 finds `{P4,P3}`.

**Minimisation reports what it did, as of 2026-09-30.** It could leave a panel
unchanged without a word, by two routes. The coverage target is compared
against the objective's coverage, which is occupancy weighted when conditions
are attached, while the run prints the unweighted figure: on the plasmid
example at 100 bp reach a 12-primer panel printed 46.7% against a 30% target
and was not reduced, because the compared figure was 19.7%. The run now logs
the primers removed, the stop reason and the figure the target was compared
against, at WARNING when nothing was removed.

**Which coverage it is compared against is now the user's to choose**, as
`--coverage-metric effective|raw` on `optimize` and as `coverage_metric` in
params.json. `plan-pool` has had that flag since it was written; `optimize`
judged panels on a metric no option could name. The default is unchanged and the
flag is a `None` sentinel, so an absent flag leaves
`panel_refinement.objective_for_optimizer`'s rule in place: effective whenever
reaction conditions are attached, raw for a condition-free library evaluator.
An explicit choice also overrides the metric a configured panel limit arrives
with, which is otherwise built with the default. Measured on the plasmid example
at 100 bp reach and a 0.30 target: the default leaves 12 primers and reports
0.197 effective, `raw` removes 8 of them and reports 0.339.

**It is the other `coverage_metric` in this package that makes the name a
hazard.** `coverage.polymerase_extension_reach` takes one too, whose values are
`realistic` and `processivity` -- same word, different question.
`search_control.COVERAGE_METRICS` names the panel pair once, the schema's enum
refuses the reach words with a message listing the permitted values, and
`resolve_search_settings` refuses them again at load time for a caller that
skips the schema. A test drives `realistic` through it for exactly this reason.

The second route was a candidate list handed to `run_optimization` directly.
The self-dimer screen lives with the candidate source, a plain list has none,
so selection could deliver a self-dimerising primer; the reduction stage counts
that as a violation, and every deletion inherited it.
`optimization_service.screen_supplied_candidates` applies the same screen to a
supplied list, keeps fixed primers, and raises `NoCandidatesError` when it
empties the list. `neoswga optimize` was never affected, since it always has a
source. A panel the reduction stage made smaller is also no longer reported as
a `set_size_mismatch` ERROR; it is a warning that says minimisation did it.

**The capability exists and `optimize` cannot reach it.** `beam_search` is
called only from `pool_planner`, so `plan-pool` can escape this local optimum
and `optimize --minimize-primers` cannot.

**On a real design it costs nothing, measured on five instances.** Greedy at
n=12, deletion at the greedy's own coverage, then a beam at width 4 over the
delivered panel plus 100 candidates:

| target | greedy coverage | deletion | beam best at n=11 |
|---|---|---|---|
| E. coli | 0.2551 | 12 | 0.2497 |
| S. aureus | 0.1985 | 12 | 0.1859 |
| M. tuberculosis | 0.1661 | 12 | 0.1619 |
| wMel against chr21 | 0.6652 | 12 | 0.6553 |
| wMel against *Drosophila* | 0.7334 | 12 | 0.7239 |

Deletion removes nothing on any of them, and neither does the beam find
anything smaller: every best-eleven falls short of the greedy's own coverage
([measurement](../validation/beam_does_not_beat_deletion_2026-09-21.md)).

**An earlier entry here claimed the opposite and was wrong.** It reported the
beam finding 11 at 0.7599 against a greedy baseline of 0.7529 on the last row
above. That baseline matches nothing else in this repository: two records
written before the beam work, Known Issue 17 and
`occupancy_and_discrimination_2026-09-19.md`, both put that panel at **0.7334**,
which is what re-running it returns. The panel the beam was asked to beat was
therefore not the panel this design delivers, and beating a weaker twelve with
an eleven is an easier problem. The earlier run was measured inline with no
script kept, so what it actually did is not recoverable.

The check that would have caught it is free: **compare a baseline against the
project's own recorded figure for the same quantity before drawing a
conclusion from the delta.** A measurement whose baseline disagrees is
reporting on a different object.

A regime explanation was tested and refuted on the way. Widening only the reach
put S. aureus at 0.8017 coverage, comparable to the Wolbachia figure, and the
beam still found no smaller panel (0.7773 at n=11). Saturation is not the
difference, because there is no difference.

So it stays unwired, and back to the reason it had before: a capability with no
demonstrated benefit, the resolution Known Issues 11 and 16 reached. Each row
is "no saving found within a ~110-oligo beam pool" rather than "none exists",
and a wider pool can only help the beam -- which is what makes the earlier
positive suspect rather than these negatives weak.

The file also states the claim the plan forbids, as arithmetic: a search that
examined every candidate examined `n` things, while the subsets number `2**n`.
Exhausting candidates is not exhausting candidate subsets, and only a toy case
like this one can produce a minimum certificate.

## A failed run leaves nothing exportable

An output directory is the only thing a later command sees, and nothing in it
carries a timestamp anyone compares. A directory whose most recent run FAILED
therefore looked exactly like one whose run succeeded: the previous run's
`step4_improved_df.csv` was still sitting there, real and stale, and `export`
turned it into an oligo order.

Task 1 wrote `design_failure.json` for exactly this, **and nothing read it.**
That is the Known Issue 8 class in artifact form: the evidence exists and the
check does not.

`export.export_is_blocked(results_dir)` is that check, and `neoswga export`
now consults it before loading anything, exiting nonzero with the recorded
stage and reason. It refuses a `failed` or `interrupted` run, and also a
`finished` one that did not qualify, since finishing is not finding something.

A record that cannot be parsed blocks rather than passes. Unknown is not
success, and defaulting the other way would make a corrupted artifact the most
permissive state available, which is every silent-zero in this file.

`cli/_failure.clear_failure_artifact` removes the record when step 4 finishes.
A record that is never cleared is as wrong as one that is never read: the user
fixes the problem, the next run succeeds, and the export refuses on evidence
that no longer describes anything.

Verified end to end: a directory carrying a failure record exits 1 and writes
no FASTA; the same directory with the record cleared exits 0 and writes one.
`tests/test_design_report_provenance.py` also pins that the CSV, the summary
JSON, the rendered report and the exported FASTA all name the same panel.

## Sequencing feedback must know which sequence it is reading

`bam_coverage.match_contigs` bound a BAM contig to a foreground reference by
**unique sequence length** when no name matched. That fallback is gone as of
2026-09-21.

Equal length is not identity. This repository ships two plasmids of 5,386 bp
each, and two chromosomes from different assemblies routinely agree. When the
lengths agree every coordinate lines up, so the depth profile reads cleanly
against a sequence the design was not made for. Feedback drives redesign --
low-depth regions become targeted additions -- so the additions would aim at
gaps in the wrong genome, and nothing downstream could notice.

Binding is now by name only: explicit alias, exact match, basename, then
chr-prefix normalisation. Each is a name agreeing with a name, which is a
claim somebody made. An unmatched prefix is skipped with a warning naming
`--contig-alias`, which is the same claim made by someone who can check it.

A name match whose LENGTHS disagree still binds, because the names are an
assertion this code should not overrule, but it now warns: that combination
means the BAM was aligned against a different version of the sequence.

`core/sequencing_feedback.require_disjoint_experiments` refuses a validation
set sharing an experiment with the training set, and names which. It also
refuses an empty set on either side: nothing held out is not the same as
nothing overlapping, and an out-of-sample result computed on no samples is an
unsupported claim rather than a weaker one. A repeated identifier within one
set is untidy rather than leakage and is allowed.

## Reading an alignment file (2026-09-24)

`core/bam_coverage.open_alignment` is the only door. Everything else in the
package is forbidden a raw `pysam.AlignmentFile` by
`tests/test_bam_reading_survives_real_aligners.py`, a source check rather than
a behavioural one because a new raw open is invisible to every behavioural
test until someone hits the failure it mishandles. Two such sites existed, and
both were reached only AFTER a successful open elsewhere, which is why neither
showed up in a failing run.

**A CRAM needs a reference and `--reference` supplies it.** CRAM stores
differences from a reference rather than sequence. htslib reports a missing
one as `OSError: truncated file`, which is a claim about the CRAM and is
wrong; `except RuntimeError` caught neither. htslib also resolves a reference
through the `UR` header field, then `REF_PATH`, `REF_CACHE`, then the EBI, so
the same CRAM reads on one machine and not another with nothing about the file
changed. `--reference` is on all three commands taking `--bam` and is threaded
to `compute_bam_depth` through `bam_depth_profile` and `bam_gaps`.

**`--contig-alias` used to mean two things.** There are two binding rules:
`match_contigs` keys on the foreground PREFIX or its basename, and
`ReferenceLayout.bind` keys on the FASTA RECORD name. The CLI documented only
the first, while the second is what every path carrying `fg_genomes` uses --
`expand-primers`, `analyze-coverage`, `iterate` -- so the documented form
failed on the commoner path, by binding nothing. Measured on a single-record
`mygenome.fasta` holding `contig_A` against a BAM contig `BAMNAME`:
prefix-keyed gave 0 gaps, record-keyed gave 1, and the run SUCCEEDED either
way, having ignored the sequencing data `--bam` exists to use.
`_record_keyed_aliases` translates the prefix form on a single-record
reference, where it can only mean one thing, and warns on a multi-record one
where a prefix names a FILE and a contig names a molecule.

**Two limits of `count_coverage`, measured and named rather than fixed**, in
`core/depth_policy.py` beside the two knobs already declared absent there:

- An `N` in a read contributes no depth -- `ACGT` + ten `N` + `ACGT` covers 8
  of the 18 bases it spans. This is the one place that module's reasoning does
  not carry through: it declines a mapping-quality floor precisely because a
  gap is what expansion then designs primers for, and an ambiguous BASE call
  is the same situation. Closing it needs the pileup API.
- `count_secondary` cannot count a bwa mem secondary record, which carries
  `SEQ` set to `*`: measured identical at 20 covered bases with the knob on
  and off. No production path sets it.

**`DepthPolicy` was unit-tested against a stub and that could not pin the
path.** The policy is handed to `count_coverage` as a `read_callback`, and
whether pysam honours it -- for which record kinds, under which of its own
default filters -- is a fact about pysam. Both ends existed and nothing walked
between them. The file now builds real BAMs covering all seven CIGAR shapes,
secondary/supplementary/QC-fail exclusion, duplicates counted, MAPQ floors on
both the bowtie2 (0-42) and bwa (0-60) scales, and records with no base
qualities as minimap2 writes from FASTA input. Everything in it passes today:
it is a ratchet, since a pysam upgrade changing `count_coverage`'s filtering
would move every depth figure here with no test to notice.

**Depth is never reported for bases the BAM cannot answer for.**
`count_coverage` CLAMPS its `stop` to the contig's length, and
`compute_bam_depth` wrote the shorter result into a `length`-sized array of
zeros. A zero there is not "no reads"; it is no sequence to have reads on. On
a fully covered 2 kb contig asked for 5,000 bases, `bam_gaps` reported a gap
of (1950, 5000) -- 3,050 bp, of which 3,000 is invented. `expand-primers`
designs oligos AT gaps, so they would target a region the BAM says nothing
about, and `calibrate-reach` would fit a reach against the same zeros. Only
`match_contigs` could reach it, since that binds on a NAME and warns rather
than refusing when lengths disagree -- which stays, because binding on the
name and inventing the depth are separate decisions and only the second is
wrong. A contig LONGER than the configured length stays silent: that reads a
prefix of it, and every base reported was observed.

**A shipped config reached it, and the tool's own advice was the way in.**
`tests/validation/genomes/params.json` names `prevotella.fna`, two
chromosomes, `fg_seq_lengths: [3168282]` -- the CONCATENATED total, since
`utility.get_seq_length` sums characters across every record. `match_contigs`
compares the prefix `prevotella` against `NC_014370.1` and `NC_014371.1`,
matches nothing, and `calibrate-reach` then said "Map one explicitly with
--contig-alias FG=BAMCONTIG". Following that produced 1,796,408 real values
and 1,371,874 fabricated zeros, **43.3% of the array**. The alias fix could
not have prevented it: `calibrate-reach` is the one command where the
prefix-keyed alias works exactly as documented.

So the message now checks whether the references are multi-record and says
this command cannot describe one, rather than offering a flag that cannot
help. `calibrate-reach` also builds no `DepthProfile`, so the evaluable mask
that keeps "measured zero" apart from "not observed" is absent on the one
path that could fabricate zeros; the refusal is what stands in for it.

**What the fabrication cost the fit**, measured through the production
`fit_reach` on that geometry with synthetic depth at a true reach of 5,000,
honest array against clamped, four seeds
(`scripts/benchmarking/clamped_tail_fit.py`):

| seed | honest rho | clamped rho | clamped plausible reaches |
|---|---|---|---|
| 7 | 0.993 | 0.272 | 2 |
| 11 | 0.992 | 0.328 | 3 |
| 23 | 0.992 | 0.318 | 2 |
| 41 | 0.991 | 0.283 | 3 |

`best_reach` recovered 5,000 in all eight runs, so the headline number
survived. What did not is the confidence: rho falls about 3.5x and the
plausible range widens from one reach to two or three, so `format_reach_table`
reported "the range this data cannot separate" and told the user to repeat the
design at both ends, for data that did separate them. And 0.27-0.33 sits ABOVE
`MIN_INFORMATIVE_CORRELATION` (0.15), so the guard never fired; a sparser real
profile with the same reduction lands below it and the run exits blaming the
BAM, which was fine.

Four seeds, one geometry, synthetic depth: a direction, not a calibrated
magnitude for a real profile. The figures first arrived from a separate audit
of this code as an inline measurement with no script -- the shape recorded as
unrecoverable under **The smallest pool** above -- and were then written out
and re-run here, reproducing digit for digit. That is why the script is in the
repository rather than the numbers alone.

**Depth costs about 72 bytes per base at peak, and that is the figure that
fails rather than slows.** `count_coverage` allocates four `array('L')` of
contig length, the caller copies them to int64 and sums, and only a 4 B/base
int32 array survives -- an 18x transient. Measured with `ru_maxrss`, one size
per process because it is a high-water mark that never falls
(`scripts/benchmarking/count_coverage_rss.py`): 72.1 B/base at 5 Mb, 72.0 at
10 Mb, 72.0 at 20 Mb. 2.5 Mb reads high (80.1) because fixed process overhead
is a larger share of a small delta, and 40 Mb reads LOW (63.3) for a reason
nobody has established -- macOS memory compression is the guess and was not
verified. So a figure for a 250 Mb chromosome is an UPPER BOUND of about
18 GB extrapolated from a range topping out at 40 Mb, not a prediction, and
the only point above 20 Mb undershoots the trend. A foreground is rarely a
chromosome, which is why this has not bitten.

**Each bound record reopens the file**, since `bam_depth_profile` calls
`compute_bam_depth` inside its loop. Measured at the *Drosophila* record count
(`scripts/benchmarking/reopen_overhead.py`): about 1.2 ms per open, so about
2.3 s for 1,870 records. Quote the seconds and not the ratio -- the fixture's
records are 1,000 bp, so counting is trivial and the ratio is inflated by
construction; on real records the ratio collapses while the toll stays. It is
2.3 s against a run that already reads a 144 Mb reference.

**Every fixture is synthesised.** There is still no BAM or CRAM in this
repository, so `calibrate-reach` has never been run against measured
sequencing depth and every reach figure remains fitted to a breadth proxy.

**A failure record no longer creates a directory where a file was asked for.**
`_failure_artifact_path` read `args.data_dir or args.output` and `makedirs`'d
it, but `--output` names a directory for `analyze-coverage` and a FILE for
`calibrate-reach`, `predict` and `report`. So a failure left a DIRECTORY at
the path the user wanted a file, the next successful run could not write
there, and the record landed where `export.export_is_blocked` never looks --
so the record whose whole purpose is to block a stale export blocked nothing.
`--output` is now consulted last and only when it is already a directory.

## Refusals, end to end

`tests/integration/test_strict_design_pipeline.py` runs the four steps over the
packaged 5.4 kb plasmid and then injects one bad configuration at a time. Each
injection is a params file somebody could write, not a monkeypatched internal:
the point is that a refusal survives argument parsing, parameter resolution,
the step's own `except` clause and the command boundary. Known Issue 8's class
is exactly a check that exists and is not reached.

**Two refusal mechanisms, and only one owes a failure record.**

A setting the SCHEMA rejects -- `coverage_reach: 0`, `polymerase: taq` -- never
enters the design path. It exits nonzero naming the permitted values, writes no
record, and wants none: there is no stage to name and no step-4 output for
anyone to mistake for current. That is an earlier and better refusal than the
design-request contract's.

A setting the DESIGN PATH rejects -- an additive with no Tm model, a
background-measured limit with no background -- owes all three: nonzero exit, a
record saying what failed, and nothing the next command will export.

**The record was not being written**, found here and nowhere else. Both
design-path refusals fire inside `resolve_design_request`, which runs BEFORE
`get_params` populates the `parameter` module, so the failure writer looked for
an output directory the module did not yet know about and wrote nowhere. The
run exited nonzero and the directory still looked like its previous success.
`_failure_artifact_path` now falls back to the `data_dir` named in the params
file, resolved relative to that file as every command resolves paths. Same
ordering trap as `warn_on_condition_drift`, reached from the other side.

Also covered: a deleted position index, a FASTA replaced after indexing so the
index is complete and describes another sequence, a reaction the filter never
recorded, and the retired `score` command name.
