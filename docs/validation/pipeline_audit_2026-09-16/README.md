# Follow-up audit: candidate search, pool size and sequencing feedback

Audited 16 September 2026 at commit `09b7d1d`. This is a review of the current
implementation after the condition-aware design changes. Production code was
not modified. Findings below distinguish executed reproductions from source
inspection and proposed improvements.

## Main conclusion

Candidate storage, chemistry consistency, panel repair and reporting have
improved. However, the default application still optimizes the capped stage-3
CSV. Keeping candidates in SQLite and indexing their background sites does not
yet make them available to the normal optimizer. Several useful new components
are tested separately but are not connected to the main command paths.

The next priority is integration and correctness, followed by a benchmark of
delivered panels over the expanded inventory. Increasing every limit or
removing QC would not establish a better design: search budgets should be
explicit, biological constraints should remain enforced, and pool quality
must be defined against coverage and specificity together.

## Status as of 2026-09-19

The conclusion above records what the audit found on 16 September and is left
as written. Each finding below now carries its own verified status; this is the
summary.

| Finding | Status | What remains |
|---|---|---|
| F1 | Partly fixed | All three commands open the shared inventory source, but only `plan-pool` reaches past the initial frontier. `optimize` and `expand-primers` still search a shortlist-sized universe |
| F2 | Mostly fixed | A coverage-target miss now triggers repair. The Stage 1 greedy still takes no occupancy objective (Known Issue 16) |
| F3 | Fixed | -- |
| F4 | Fixed | -- |
| F5 | Fixed | -- |
| F6 | Fixed | -- |
| F7 | Fixed | -- |
| F8 | Fixed | -- |
| F9 | Fixed | The fitted reach stays a model parameter for one dataset; that is a property of the method, not a defect left open |

F1 and F2 were checked against the code rather than against commit titles, and
both turned out narrower than the commits that closed them suggest. The
headline claim of F1 -- that candidates beyond the `max_primer` shortlist
cannot improve a default design -- still holds for `optimize`.

## Verification

- **98 existing focused tests passed**, including the small end-to-end pipeline,
  inventory, retention modes, real-optimizer repair, focused/full metric
  agreement, BAM coverage and reach calibration. Run time: 19.16 seconds.
- [probes.py](probes.py) reproduces nine additional cases using production
  functions and small synthetic inputs. [probe_results.json](probe_results.json)
  records their outputs. These are audit demonstrations, not claims that the
  existing tests failed.
- No full Wolbachia redesign or sequencing-data analysis was performed in this
  audit. Wolbachia figures below are explicitly from the committed benchmarks.

```bash
python -m pytest tests/test_candidate_inventory.py tests/test_pool_design_end_to_end.py tests/test_pool_plan_repair.py tests/test_pool_plan_repair_runs_on_a_real_optimizer.py tests/test_pool_metrics_agree_with_the_full_evaluation.py tests/test_bam_coverage.py tests/test_reach_calibration.py tests/test_candidate_retention_mode.py tests/test_post_gini_retention_mode.py -q
PYTHONPATH=. python docs/validation/pipeline_audit_2026-09-16/probes.py
```

The probe script also requires the optional `pysam` dependency, available in the
audit environment. Temporary references, indexes and BAMs are removed on exit.

## Findings requiring correction

### F1 — P1: retained candidates are not reached by normal optimization

**PARTLY FIXED, verified 2026-09-19.** `core/candidate_source.py` holds the one
rule all three commands ask: `plan-pool` (`cli/plan_pool.py:66`), `optimize`
(`core/unified_optimizer.py:815`) and `expand-primers`
(`core/primer_expansion.py:926`) open the inventory when the directory has one
and the supplied list otherwise, and the run logs eligible, frontier and
unexamined counts.

**What is NOT fixed:** the frontier opens at exactly the size of the list the
command would have read anyway, and only `plan-pool` reaches past it --
`pool_planner.py:446-475` is the sole caller of `source.advance()`. `optimize`
and `expand-primers` take `source.initial()` and stop. So on the Wolbachia
example the 489,836 candidates beyond the 2,000-primer shortlist are still
unreachable by a default `optimize`, which is this finding's stated effect. What
changed for those two commands is the ORDER and the reporting, not the reach.
**Measured 2026-09-19 before doing that work, and it should not be done as
stated.** Giving `optimize` the whole 20,670-candidate inventory instead of the
2,000-primer shortlist changes the delivered panel substantially -- Jaccard
0.500, 0.263 and 0.171 at n=6, 12 and 24, with 2, 7 and 15 of the delivered
primers coming from outside the shortlist -- and the change trades specificity
for coverage: selectivity density falls 12 to 42 per cent while coverage rises
1 to 6 points. Two causes compound. `max_primer` cuts on `bg_count / fg_count`
ascending, so the candidates a refill reaches bind the host twice as often and
the target half as often (median 20/3 against 10/6); and Stage 1's greedy
carries no specificity term to resist them, which is Known Issue 16.

The earlier `plan-pool` refill measurement does not transfer, because there the
refill fired only on a size row that could not satisfy its constraints. A
`optimize` run with no panel limits has no failing row, so an unconditional
refill has nothing holding it to the specificity the shortlist was selected
for. Close this by making Stage 1 specificity-aware, or by triggering a refill
only on an unmet configured panel limit -- not by refilling unconditionally
([measurement](../frontier_refill_on_optimize_2026-09-19.md)).

**Source evidence:** `neoswga/cli/plan_pool.py:78-85` reads `step3_df.csv` unless
an explicit CSV is supplied. `neoswga/core/unified_optimizer.py:953-956` does the
same for standard `optimize`. Neither path constructs a `CandidateProvider` or
requests expansion. The only production construction of `CandidateProvider`
outside its own module is the library-level `design_sweep` callback interface.
The expansion command also defaults to stage 3 (`neoswga/cli/iterate.py:26`).

**Effect:** `all_qc` changes storage/indexing, not the default search universe.
The published Wolbachia retention benchmark stores/indexes 491,836 hard-QC
survivors but still writes 2,000 stage-2 candidates. Default score carries these
2,000 forward. The 489,836 others cannot improve the default result.

**Correction:** Connect one inventory-aware provider to `optimize`, `plan-pool`
and expansion. Use the CSV ranking only as the initial batch. Expand on
unmet targets or stalled improvement, and expose eligible, indexed, examined
and unexamined counts plus budget stop reasons. Test a real CLI design where
the only primer covering a required interval is outside the initial shortlist.
Do not solve this by eagerly loading every pairwise dimer relationship.

### F2 — P1: coverage-target failures do not trigger condition-aware improvement

**MOSTLY FIXED, verified 2026-09-19.** `pool_planner.py:277-278` now reads
`if repair and (reasons or missed is not None)`, where `missed` is
`max(targets)` when the delivered coverage falls short of it, so a coverage
miss triggers repair with no constraint violated. The most demanding target
drives it, because repairing to the lowest would leave every higher row
reporting `not_found` beside a panel never asked to reach it. The beam is no
longer gated on a whole-pool bound it could not meet: `_beam_candidates`
slices the pool to what the remaining budget affords, and the record says
`not affordable within the remaining budget` when even that does not fit.
`refine_hybrid_stage2` receives the objective as of 2026-09-18 (see CLAUDE.md,
"The objective never reached the stage that refines").

**What is NOT fixed:** this finding's remark that "the dominating-set
initializer also receives no occupancy objective" still holds. `optimize_greedy`
accepts an `objective` and no production caller supplies one; that is tracked
as Known Issue 16 in CLAUDE.md, with the measurement of what it costs.

**Source evidence:** `neoswga/core/pool_planner.py:235` runs repair only when
`reasons` is nonempty. These reasons cover specificity, missing coverage and
dimers, but not coverage below a requested target. `refine_hybrid_stage2`
(`neoswga/core/swap_refinement.py:149`) still calls the raw-bin swap objective
without `PoolObjective`. The dominating-set initializer also receives no
occupancy objective.

**Executed reproduction:** A returned singleton has 50% effective coverage and
passes specificity. Another available singleton has 95% coverage and passes
the same limits. For a 90% request, `plan_pool` reports `not_found` and
`repair.attempted=false`. A controlled optimizer stub isolates this planner
behavior; it does not claim that every real optimizer returns the worse pool.

**Effect:** A pool can be accepted as specificity-feasible while the relevant
coverage search never runs. Increasing pool size may compensate for an
objective mismatch rather than identify the smallest satisfactory pool.

**Correction:** Use the shared chosen coverage metric in initial selection and
all refinement, including when specificity already passes. Treat target misses
as reasons to improve or expand search. Retain feasible incumbents and test
raw/effective disagreement against an exhaustive oracle on small instances.

**Related limits:** Repair uses a fixed beam width of four and only runs the
beam when `4 * candidate_count * panel_size` fits its remaining evaluation
budget. At 2,000 candidates and 12 primers, that bound is 96,000 evaluations,
above the default 10,000 even before swap evaluations. Consequently, beam
fallback is generally unavailable at realistic shortlist sizes under default
budgets. Report actual search termination and make budgets controllable across
the whole search; do not present this fallback as exhaustive exploration.

### F3 — P1: inventory eligibility can outlive the QC decision that created it

**FIXED, verified 2026-09-19.** Two changes, and the finding needed both.
`qc_policy_fingerprint` (`core/candidate_inventory.py:49`) digests the RESOLVED
hard-QC thresholds rather than naming a policy version, so two runs with
different GC, Tm or frequency limits no longer share a verdict; an absent
threshold is encoded distinctly from one set to its default, because "not
configured" and "configured to the default" are different statements about what
was enforced. And assessments carry a `generation`, with `iter_eligible`
returning only `MAX(generation)` for the condition and policy, so a survivor of
an earlier, looser run is no longer eligible after a stricter one.
`tests/test_inventory_eligibility_expires.py` pins it.

**Source evidence:** `record_stage2_inventory`
(`neoswga/core/candidate_inventory.py:236`) only upserts current survivors as
passing. It does not invalidate prior survivors absent from the new run.
The assessment key uses reaction fingerprint and a fixed policy-version
string, not the resolved QC thresholds, reference identity and run generation.

**Executed reproduction:** Write two passing candidates; rerun stage-2
recording under the same reaction with only one passing. Both are still
returned by `iter_eligible`.

**Effect:** Tightening a GC, sequence, frequency or exclusion rule can leave
previously admitted candidates eligible. This becomes particularly important
when F1 is fixed and the inventory is used directly.

**Correction:** Publish an atomic current-run inventory view keyed by full
reference, condition and QC-policy identity. Mark assessed failures explicitly
or exclude candidates not assessed in the current generation. Preserve history
separately. Test stricter rules and changed references in the same output
directory, not only independent writes into empty databases.

### F4 — P1: `--design-grid` is accepted but not executed

**FIXED, verified 2026-09-19.** `cli/plan_pool.py:294` reads `args.design_grid`
and dispatches to `_run_design_grid`, which calls `load_grid_file` and
`design_sweep` and writes `design_sweep.json`. A condition the filter never
recorded is reported as having no eligible candidate rather than designed over
an empty pool, which is the "distinguish unassessed from assessed-and-failed"
half of the correction. Its entries are gone from both the inert-option and
unreachable-capability allowlists, so the two ratchets now hold it.

**Source evidence:** `neoswga/cli/plan_pool.py:237` advertises independent
designs across lengths and conditions. `run_plan_pool` never reads
`args.design_grid` or calls `load_grid_file`/`design_sweep`. The helpers exist,
but are not connected to the command handler.

**Effect:** A user can request a condition grid and receive an ordinary
single-condition design. Also, the library sweep reads existing assessments;
it does not itself perform new condition-specific QC. Missing assessments can
be described as no eligible candidates rather than not yet evaluated.

**Correction:** Wire the command end to end or reject the unsupported flag.
Reassess a condition-independent candidate universe for every grid condition,
and distinguish unassessed from assessed-and-failed. Test invocation of the
actual CLI with two conditions and verify two independently designed results.

### F5 — P1: `all_qc` can still fail at Gini; `post_gini` eligibility is inconsistent

**FIXED, verified 2026-09-19.** Both halves, separately.

`core/stage2_recording.py:95-123` indexes the background and records the
inventory BEFORE `check_gini_stage_kept_something` runs. The guard still
raises with the same message, but under `all_qc` a run where evenness happened
to be unmeasurable no longer aborts before a single inventory row is written.
Those candidates had cleared every declared requirement; what they had not done
was bind often enough for their spacing to be measured, which is what
`min_gini_sites` exists to say is not a judgement about them.
`tests/test_all_qc_seeds_before_the_evenness_guard.py` pins it.

The eligible set and the indexed set are now the same set in both modes, pinned
by `tests/test_post_gini_admits_only_what_it_indexes.py`. `post_gini` is an
ADMISSION policy, not an indexing budget: a candidate missing the evenness gate
is recorded with an explicit failed assessment naming it and is not eligible.
The mode is part of the admission-policy digest, so switching it opens a new
generation rather than inheriting the other mode's verdicts.

**Source evidence:** `neoswga/core/pipeline.py:1336` calls the empty-Gini guard
before inventory recording at line 1360, regardless of retention mode.
`record_stage2_inventory` marks every hard-QC survivor eligible even when
`post_gini` indexed only a subset. `iter_eligible` and `CandidateProvider`
do not filter on `passed_gini` or `indexed`.

**Executed reproductions:** The Gini guard raises when one pre-Gini survivor
has no post-Gini row. Separately, a provider returns a candidate explicitly
recorded as unindexed under the post-Gini indexing pattern.

**Effect:** Candidates useful only in local gaps can remain inaccessible even
with `all_qc`. Expanding a `post_gini` inventory may access candidates whose
background sites were not indexed, contrary to the retention documentation.
`ensure_positions` does not resolve this reliably: it is normally unattached
to a cache, and when attached it asks whether any prefix has positions, not
whether every required target/background index contains the candidate.

**Correction:** In `all_qc`, seed from hard-QC survivors even when none has
measurable Gini. Define explicitly whether `post_gini` is a hard admission
policy or an indexing budget. If it is an indexing budget, index missing
background entries before evaluation. Validate index membership separately
from the measured count: a present zero-match entry is valid data.

### F6 — P2: count-table provenance misses interior reference changes

**FIXED 2026-09-19** (`45d4610`). `genome_fingerprint` is the SHA-256 of the
whole file, cached on `(st_size, st_mtime_ns)` so it costs one pass per input
per run rather than one per k value. The algorithm is named in the sidecar as
`digest_algorithm`, and a record written under the old partial hash is treated
as UNKNOWN rather than as a mismatch -- without that distinction every existing
data directory failed step 2 on upgrade with "counted from a different genome",
which is both alarming and untrue (28 tests caught it).

**Source evidence:** `genome_fingerprint`
(`neoswga/core/kmer_counter.py:327`) hashes file size and first/last 1 MiB.

**Executed reproduction:** A base substitution in the middle of a roughly
3 MiB FASTA leaves the fingerprint unchanged.

**Effect:** A same-length reference update can reuse old k-mer counts. This is
relevant to redesign against sample-specific assemblies or consensus changes.

**Correction:** Use a full content digest computed once per input per run and
reuse the digest across k values. Include it in position-index and inventory
provenance, not only count tables. Checking that record-start metadata exists
is not a full reference-identity check. Standard `optimize` also lacks the
explicit `require_record_metadata` call present in `plan-pool`.

### F7 — P1: BAM-guided redesign does not yet optimize recovery of observed gaps

**FIXED 2026-09-19** (`8a83591`, `d2e834e`, `50e7342`). Both halves.
`core/deficit_objective.py` scores a candidate by the depth deficit it
recovers, weighted per base, in place of requiring the binding site itself to
lie inside a gap; the all-or-nothing fallback to the full candidate list is
gone and candidates are matched to gaps dilated by the coverage reach.
`core/design_context.py` is the single params resolution both `plan-pool` and
`expand-primers` use, so the command that ADDS to a panel no longer runs at a
hard-coded 3 kb reach with `conditions=None`.

**Source evidence:** `PrimerExpander._filter_candidates_to_gaps`
(`neoswga/core/primer_expansion.py:516`) requires the binding site itself to lie
inside the gap. `expand` falls back to the full candidate list when this
leaves fewer candidates than `target_new`. After filtering, selection uses
the ordinary global objective, not a depth-deficit objective.

**Executed reproduction:** For a gap from 1,000 to 2,000 bp, a candidate with a
site at 900 bp and a configured reach of 3,000 bp is excluded, although its
modeled window covers the entire gap.

**Other source findings:** Expansion's CLI constructs `PrimerExpander` without
the resolved chemistry or reach. The expander's hybrid builder leaves chemistry
unset and uses its default polymerase; its reach defaults to 3,000 bp. The
heterodimer limit can still be inherited from global parameters, so it would
be incorrect to say that all configured dimer settings are necessarily lost.
There is nevertheless no shared explicit configuration contract here. Its
geometric coverage path does not pass record starts to the coverage graph.
Its reported `gap_coverage` uses the number of gaps, not recovered gap bases;
splitting one long gap into two smaller gaps can worsen this number despite
covering more bases. Reported post-design coverage is a prediction, not a new
BAM measurement.

**Correction:** Reuse the main evaluator/configuration. Rank candidates by
incremental coverage of observed deficit intervals, including nearby sites
whose modeled reach overlaps them. Keep a fixed validated subset if requested,
allow add/swap/drop comparisons, and enforce limits on the complete delivered
pool. Preserve gap source labels instead of silently treating observed and
geometric gaps as equivalent.

The default minimum gap length is 10,000 bp. Smaller observed deficits are
ignored unless the user changes it. Make this resolution explicit in the
report and relate it to the user's breadth requirement, rather than treating
absence of a reported large gap as complete recovery.

### F8 — P2: sequencing coverage ingestion needs stronger semantics

**FIXED 2026-09-19** (`b68248e`, `b952ee7`). `core/reference_layout.py` binds a
prefix's records to BAM records by name and length, so a two-record FASTA maps
to the two BAM records with their offsets instead of returning an empty
mapping; a name match with a conflicting length is reported rather than
accepted. `core/depth_policy.py` makes the read-inclusion rules explicit and
recorded (MAPQ, base quality, duplicates, supplementary, secondary, QC-fail)
rather than implied by pysam's `'all'` default, which skips UNMAP, SECONDARY,
QCFAIL and DUP but NOT SUPPLEMENTARY. It deliberately omits overlapping-mate
and deletion knobs, which `count_coverage` cannot enforce.

**Source/execute evidence:** `match_contigs`
(`neoswga/core/bam_coverage.py:45`) maps each foreground prefix to one BAM
record. A two-record FASTA represented by one prefix/total length does not map
to BAM records of the two individual lengths. The probe returns an empty
mapping. An alias cannot represent the required two records plus offsets.
Same-name matches are accepted without checking matching sequence lengths;
unique-length fallback alone does not prove reference identity.

`compute_bam_depth` at line 127 uses base-quality threshold zero and no
explicit mapping-quality/supplementary filter. A synthetic BAM with one
high-quality primary alignment, one MAPQ-zero alignment and one supplementary
alignment gives depth three. This may be an intentional definition of aligned
depth, but it is not a configurable high-confidence recovery endpoint. The
default callback already excludes unmapped, secondary, QC-failed and
duplicate-flagged reads; do not claim these are all included. See the
[pysam API](https://pysam.readthedocs.io/en/stable/api.html#pysam.AlignmentFile.count_coverage).

**Correction:** Map BAM records to validated FASTA records and offsets. Fail
or require explicit partial-analysis selection for unmatched target records.
Expose and record MAPQ/base-quality, supplementary, duplicate and overlapping
mate policies. Preserve a mask for non-evaluable bases. Report both included
and excluded denominators so removing difficult regions does not inflate
apparent whole-genome recovery.

### F9 — P2: reach fitting is exploratory and fails for short input

**FIXED 2026-09-19** (`7a4a653`). Handled, with the statistical half the more
consequential.

- A target shorter than one bin returns a non-informative result naming the bin
  count, in place of the raw numpy `ValueError`.
- Bins are formed with `np.add.reduceat` from per-record start offsets, so the
  trailing remainder counts and no bin spans a join between two records of a
  concatenated prefix. A 10,500 bp target now yields 11 bins, not 10.
- `predicted_depth(..., circular=True)` wraps the kernel. The clipped form
  loses a share of each near-boundary triangle that GROWS with the reach being
  tested: for a site 5 kb in, 0.0% at reach 3 kb, 12.5% at 10 kb, 36.7% at
  35 kb and 43.1% at 70 kb -- a bias against exactly the large reaches the fit
  exists to weigh. It is refused, with a note, for a multi-record prefix, where
  the end of the last record is not adjacent to the start of the first.
- `correlation` is relabelled in-sample everywhere it is reported, and
  `cv_correlation` reports the same quantity under spatially blocked
  cross-validation: each fold re-runs the whole grid selection on the remaining
  bins and scores the reach it picks on a contiguous block it never saw. Blocks
  are contiguous rather than interleaved because neighbouring bins share
  binding sites and mappability.
- `plausible_reaches` carries every grid point within 0.02 rho of the winner,
  and `at_grid_edge` flags an optimum at the top of the grid.

**Measured, and it qualifies the fix.** On null depth over twelve seeds the
in-sample maximum averages +0.058 where the truth is zero -- about a third of
the 0.15 informativeness threshold, purely from selecting among eleven
candidates. But that bias is not larger than the fold-to-fold noise in the
held-out estimate (sd 0.052), so the held-out figure landed ABOVE the in-sample
one on 2 of 12 seeds. The pair is therefore a stability check, not a corrected
value, and the documentation says so rather than claiming cross-validation
removes the optimism. `tests/test_reach_fitting_is_honest.py` asserts the mean
across seeds for that reason, not an inequality on one.

The remaining model limits below are unchanged and still stand: the kernel is
symmetric and identical for every site, and the result is a model parameter
fitted to one dataset rather than a measured constant of the polymerase.

**Executed reproduction:** `fit_reach` with 500 depth values and the default
1,000-base bins raises `cannot reshape array of size 500 into shape (1,1000)`
(`neoswga/core/reach_calibration.py:130`). Longer non-divisible inputs discard
their final partial bin.

**Model limits from source:** The fit uses identical symmetric triangular
kernels for all positions, selects reach by in-sample Spearman correlation,
does not implement circular wrapping, and the CLI fits the longest matched
contig. It does not estimate an occupancy-weighted sequencing breadth model,
individual oligo efficacy, amplification kinetics or uncertainty in the reach.
Its correlation threshold is a heuristic, not independent validation.

**Correction:** Handle partial/short bins, record boundaries and topology.
Report flat or boundary optima and fit uncertainty; evaluate predictions on
held-out spatial blocks and preferably independent reactions. Treat the fitted
reach as a model parameter for sensitivity analysis until transfer is tested,
not a direct measurement of polymerase extension or expected recovered bases.

## Stage-by-stage assessment

| Stage | Improved/current behavior | Remaining limitation | Recommended action |
|---|---|---|---|
| `count-kmers` | Counts canonical motifs across configured k values; checks cached-table provenance | Only counted lengths are available; partial digest can miss reference edits | Keep exact enumeration, validate all required lengths and full reference identity; counting must not impose chemistry ranks |
| `filter` | Canonical chemistry in primary QC; broad background indexing and durable survivors | Early loader rejects are not recorded in the inventory; Gini/cap still govern the CSV; stale assessments persist | Separate hard eligibility from ranking, record QC reasons/current generation, and keep all hard-QC survivors reachable |
| `score` | Default preserves stage-2 candidates; legacy efficacy model is opt-in | Carries only the capped CSV; inventory `examined` counts are inferred from stage-3 membership, not actual search | Preserve inventory references and record scoring as annotations; instrument real search examination |
| `optimize` / `plan-pool` | Strict delivered-panel dimer checks; focused metrics; optional specificity repair | Default CSV truncation, raw versus effective objective mismatch, finite unreported/prioritized search paths | Connect provider, shared objective and consistent configuration across entry points; expose stop reasons and stable budgets |

Standard `optimize` also applies an exact-count foreground/background prefilter
by default, removing up to a configured fraction (default 20%). It is not the
same as the final per-base occupancy-density limit, and it does not run in the
same way in `plan-pool`. For maximum-quality design, make such pruning an
explicit approximation or ordering heuristic, not an unreported permanent
exclusion. Preserve explicit hard background/exclusion limits.

Frequency limits for longer oligos still scale by `4 ** (10-k)`. The resolved
thresholds are now logged. This is a declared model choice rather than a newly
discovered bug; compare actual site thresholds across lengths before treating
a longer-oligo run as directly comparable.

Inventory counters also need correction: `counted` currently counts candidates
written after hard QC, `hard_qc_passed` aggregates passing assessments across
conditions/policies, and `examined` counts stage-3 assessments. These are not
the full enumerated universe, current-condition pass count, and actual search
evaluations respectively (`candidate_inventory.py:405`).

## How many oligos should be in the final pool?

There is no defensible universal number, and the current default maximum of
24 is a user-interface search bound, not a result derived from the target.
Pool size should be selected from a curve after F1/F2 are corrected.

1. Declare the target/background references, chemistry, specificity limits,
   coverage metric and target breadth. For lab validation, specify depth
   thresholds and sequencing effort separately from modeled coverage.
2. Search increasing sizes at comparable documented budgets and multiple
   deterministic starts. Carry the best feasible panels forward; permit
   deletion and swaps as well as addition. A smaller panel found while testing
   a larger budget should remain available to the overall recommendation.
3. Expand the inventory when a small pool misses the target before concluding
   that more oligos are needed. Distinguish a candidate-limited search,
   computation-limited search and incompatibility-limited search.
4. Report the smallest qualifying found panel, its nearby alternatives,
   marginal gain per added oligo, maximum persistent gap and reach sensitivity.
   Do not choose solely from average genome coverage if required regions remain
   uncovered. Keep per-target results for multiple target genomes.
5. Retain strict pool-wide dimer/background constraints while removing redundant
   oligos. More oligos can create additional incompatibilities and cost without
   sufficient coverage gain.
6. Report concentration per oligo and nominal total pool concentration. The
   present model assumes fixed per-oligo concentration; fixed-total comparisons
   require recomputation and an explicitly supported concentration model.

The committed [Wolbachia search-budget measurement](../wolbachia_search_budget_2026-09-16.md)
used the same **2,000-candidate** pool. At 16 oligos, increasing the evaluation
budget from 1,000 to 100,000 changed modeled coverage from **74.38% to 76.42%**
and density from **26.93 to 21.86**. These figures demonstrate search sensitivity,
not the required pool size or measured recovery. That study did not test the
full retained inventory, nor a stringent binding constraint. The
[retention benchmark](../wolbachia_retention_benchmark_2026-09-16.md) stopped
after `filter`; it cannot establish the effect on final panel quality.

## Using sequencing results to improve a designed pool

Existing building blocks are `analyze-coverage --bam`, `calibrate-reach --bam`,
`expand-primers --bam`, `contract-set` and `rescore-set`. They are useful starting
points, but the current combination is not a validated closed-loop optimizer.
Correct F7-F9 before relying on it for automatic redesign.

### 1. Capture a complete experimental observation

Associate the exact delivered oligo list and concentrations with the reference
digest, reaction conditions, sample/input composition, run identifier and mapped
BAM. Include unamplified controls where available and independent reactions.
Use the target plus relevant background for assessing mapping specificity;
target-only depth cannot establish enrichment against the background.

Derive per-record depth and windows, breadth at declared thresholds (for example
1x and 10x), target read fraction, coverage variability and lengths of persistent
low-depth intervals. Record the mapping/filtering policies from F8. Compare
panels at equal total sequencing effort for practical recovery efficiency;
also compare at equal target-aligned effort to assess coverage uniformity.
These answer different questions and should be reported separately.

### 2. Distinguish recoverable deficits from other causes

Classify low-depth regions using reference identity, mapping ambiguity,
sequence divergence, expected presence and replicate consistency. A missing
region or unmappable repeat is not automatically a primer-design failure.
Keep excluded/non-evaluable regions visible in the whole-reference summary.

Separate observed low-depth intervals from predicted no-site intervals. Their
overlap is informative, but one is not evidence that the other has disappeared.
An updated design has predicted improvement until the next sequencing run.

### 3. Optimize additional coverage where it is needed

A proposed first implementation can use an observed deficit as a **priority
weight**, without claiming a calibrated predictor of read depth. For region r:

```
weight[r] = max(0, 1 - observed_depth[r] / desired_depth)
gain(candidate | pool) = sum_r length[r] * weight[r]
                         * (modeled_coverage[r, pool + candidate]
                            - modeled_coverage[r, pool])
```

Use depth and masks from a declared normalization policy and aggregate
reproducible deficits across reactions. Enforce global coverage, specificity,
full resolved chemistry and selected-pool dimer constraints alongside this
weighted objective. Count nearby candidate sites whose reach overlaps a
deficit; do not require the site to fall inside it. Diminishing marginal gain
prevents repeatedly targeting a region already addressed by other selected
oligos.

Compare small add-only designs (retaining the current pool), fixed-size swaps,
and smaller pools after removal of redundant oligos. Report the number added,
removed and retained and the predicted coverage change for each persistent
deficit. Candidate discovery must draw from the retained inventory, not just
the original selected panel or capped CSV.

### 4. Learn cautiously across independent panels

One panel's depth profile can identify regional deficits, but generally cannot
identify the causal contribution of each oligo: their binding windows overlap
and amplification is coupled. Do not mark individual primers as failed solely
because a nearby region is under-covered. Explicitly measured failures can
remain an exclusion input.

With multiple independently tested pools, fit a regularized model using
regional sequence/site features, panel composition and recorded conditions.
Compare it with simple geometry and baseline-depth models. Validate across
held-out reactions/panels, not randomly partitioned neighboring bases. Use
spatially blocked evaluation for within-genome diagnostics. Reserve new panels
for prospective assessment before advertising calibrated recovery predictions.

Primary SWGA studies assess actual sequencing outcomes, illustrating why a
design metric and a sequencing endpoint need separate validation:
[Clarke et al. 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5870857/) and
[Dwivedi-Yu et al. 2023](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010137).
The repository's existing evidence matrix correctly distinguishes qPCR yield,
read enrichment and per-base breadth. Its held-out discrimination report
should say **no clear transferable advantage was established**, rather than
literally saying no predictor beat baseline: its own table shows 0.67 versus
0.62 for one predictor, without establishing a reliable improvement.

## Recommended next work

1. **Repair integration and provenance first:** F1, F3, F4 and F5. Add actual
   command-level regressions so an unused helper cannot satisfy acceptance.
2. **Align search and choose pool size:** F2, then a full-path Wolbachia benchmark
   comparing initial shortlist, post-Gini inventory and all hard-QC candidates,
   with equal budgets and also documented larger-budget runs. Include a
   specificity constraint that affects selection. Preserve all current results.
3. **Make coverage feedback reliable:** F6-F9, with a small multi-record BAM
   fixture and known deficits. Forward the exact same chemistry and constraints
   through design, expansion, contraction and reporting.
4. **Add deficit-weighted add/swap/drop design:** start with explicit priority
   weights and measured endpoints, then assess prospective sequencing outcomes.
   Do not delay the integration fixes while waiting for calibration data.

The objective is the best pool found under explicit constraints and search
effort. Neither retaining every candidate nor passing these tests proves a
global optimum or a particular fraction of experimentally recovered genome.
