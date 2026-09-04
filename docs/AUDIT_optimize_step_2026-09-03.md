# Audit: the optimize step and the path that delivers an oligo pool

Date: 2026-09-03. Branch `audit-alternatives-and-scaling`, commit `232d683`.

Scope: step 4 (`neoswga optimize`), the five optimizers it dispatches to, and
the commands between it and the oligos a user orders (`interpret`, `report`,
`export`). Every finding was re-derived against this tree. Findings marked
**[verified here]** were reproduced directly in the course of writing this
document; the remainder come from module-level audits that reported their own
reproductions, and are marked accordingly.

Measurements use `examples/plasmid_example` (pcDNA 6157 bp against pLTR
6258 bp) with its shipped `params.json`, which sets `max_dimer_bp: 3`.

Every finding below survived a green suite when it was written: `pytest -k
"optimiz or dimer or interpret or ensemble or coverage" -n 8` reported 823
passed, 3 skipped against the tree that contained all of them. That is the point
of the document as much as any individual defect.

**Status, 2026-09-04.** Twenty of the twenty-one findings are fixed, each with a
guard test written to fail against the old behaviour first; finding 21 is a
list of smaller items and is partly done. Fixed findings carry a note under
their heading naming the test. The suite now stands at 3929 unit tests and 150
integration and CLI tests, with both size ratchets and the lint gate green.

## Summary

The optimizers are deterministic and the surrounding scaffolding (empty-input
guards, blacklist re-injection, factory dispatch, seeding) is sound. The defects
cluster in three areas.

**Configuration does not reach the optimizer.** The selection weights that
`--application` sets, and the `max_dimer_bp` from params.json, are dropped in an
inner constructor on `hybrid` and `background-aware`; the log states they were
applied. `--auto-size` rebinds the polymerase from raw JSON and discards
`--polymerase`, so the panel is designed for the wrong enzyme at the wrong
reach. `background-aware` drops the configured Tm window and throws away two
thirds of the candidate pool. All three are the same defect class as the
already-fixed `max_dimer_bp` routing bug, and two of them sit on the default
path.

**Two optimizer stages do not do what they are named for.** Hybrid's stage-2
coverage guard is inverted by a dataclass equality subtlety, and preferentially
removes the primer that solely covers a region. Background-aware's pruning stage
is disabled on any run whose coverage is at or below 95 per cent, which is the
normal case.

**The dimer constraint is not enforced or reported anywhere on the delivery
path.** A pool exceeding the configured threshold is graded Excellent and
described as ready to order.

| # | Severity | Finding |
|---|---|---|
| 1 | Critical | Selection weights and `max_dimer_bp` never reach the default optimizer |
| 2 | Critical | `--auto-size` silently discards `--polymerase` and designs for the wrong enzyme |
| 3 | Critical | Hybrid stage 2 prefers to remove the sole coverer of a region |
| 4 | Critical | Background pruning is a no-op whenever coverage is at or below 0.95 |
| 5 | Critical | The network greedy's background penalty can become a background reward |
| 6 | High | `background-aware` ignores the configured Tm window and `--max-extension` |
| 7 | High | `--minimize-primers` with `ensemble` crashes and delivers no pool |
| 8 | High | The quality report divides a 0-1 dimer risk score by ten |
| 9 | High | No stage checks the delivered pool against the configured `max_dimer_bp` |
| 10 | High | Multi-target coverage pools cross-genome coordinates |
| 11 | High | The amplification network joins sites on different molecules |
| 12 | High | The clique optimizer never receives its own configuration |
| 13 | Medium | `interpret` reads neither the summary nor the validation file |
| 14 | Medium | The interpreter's dimer assessment looks for absent column names |
| 15 | Medium | `--validate-simulation` cannot work and fails as a warning |
| 16 | Medium | The ensemble tie-break is list order, and ties are the common case |
| 17 | Medium | Three more inert options: `--min-fg-bg-ratio`, `--no-position-cache`, `iterations` |
| 18 | Medium | The report prints the no-background sentinel as a measured 1000000x |
| 19 | Medium | `coverage_uniformity` rewards clustering |
| 20 | Medium | One lowercase primer defeats every dimer check |
| 21 | Low | Six smaller items, listed at the end |

---

## 1 (CRITICAL) Selection weights and `max_dimer_bp` never reach the default optimizer -- FIXED

> **Fixed.** `hybrid` and `background-aware` now receive all five values; the audit's own table reads correctly for every method. Test: the parameterised matrix in `tests/test_optimizer_config_reaches_optimizers.py`, which was red in nine places.

**[verified here]**

`HybridOptimizer` builds its inner `NetworkOptimizer` from an explicit argument
list that omits `tm_weight`, `dimer_penalty`, `max_dimer_bp` and `template_gc`
(`hybrid_optimizer.py:310-325`). `BackgroundAwareBaseOptimizer` builds its
`HybridOptimizer` without `max_extension` or `uniformity_weight`
(`background_aware_optimizer.py:692-711`).

Constructing all three through `OptimizerFactory` with the clinical profile
(tm 0.30, dimer_penalty 0.45, uniformity 0.05), `max_dimer_bp` 3 and
`--max-extension 25000`, then reading the inner `NetworkOptimizer`:

| | requested | network | hybrid | background-aware |
|---|---|---|---|---|
| `tm_weight` | 0.30 | 0.30 | **0.00** | **0.00** |
| `dimer_penalty` | 0.45 | 0.45 | **0.00** | **0.00** |
| `uniformity_weight` | 0.05 | 0.05 | **0.00** | **0.00** |
| `max_dimer_bp` | 3 | 3 | **4** | **4** |
| `max_extension` | 25000 | 25000 | 25000 | **70000** |

The run log asserts the opposite. `neoswga optimize --application clinical`
prints:

```
INFO: application='clinical' -> selection weights applied: tm=0.30, uniformity=0.05, dimer_penalty=0.45
```

The dimer penalty is the only dimer consideration anywhere in hybrid's stage-2
refinement, so the default optimizer applies no dimer penalty at all whatever
`--application` says, and substitutes 4 for a configured `max_dimer_bp` of 3.
This is the mechanical explanation for finding 9.

`tests/test_optimizer_config_reaches_optimizers.py` already exists for exactly
this defect class and should be extended one constructor deeper.

## 2 (CRITICAL) `--auto-size` silently discards `--polymerase` -- FIXED

> **Fixed.** The rebind is deleted. `--polymerase bst --auto-size` now designs at bst's 1000 bp reach, not phi29's 3000.

**[verified here]**, reported by the factory module audit.

`cli/pipeline.py:570` resolves the polymerase correctly: the `--polymerase` flag
first, then the adapted global, then params.json, with alias normalisation.
Inside the auto-size block, `cli/pipeline.py:635` rebinds the same local from
raw JSON:

```python
polymerase = json_data.get("polymerase", "phi29")
```

That value is what reaches the optimizer at `:771`. So whenever `--auto-size` is
given, the flag is discarded and the params.json value, or a literal `phi29`,
wins.

Measured on the plasmid example with `--polymerase bst`:

| run | extension reach | primers delivered |
|---|---|---|
| without `--auto-size` | 1000 bp (bst) | 6 |
| with `--auto-size` | **3000 bp (phi29)** | **12** |

Both runs print `Polymerase: bst (config applied to optimizer)` at the top. The
second designs a phi29 panel at 30 C while telling the user it is a bst panel
for 60-65 C. The polymerase sets the coverage reach, the Tm window and the
reaction chemistry, so this is not a labelling error: a different pool comes
back, selected against the wrong enzyme's physics.

## 3 (CRITICAL) Hybrid stage 2 prefers to remove the sole coverer of a region -- FIXED

> **Fixed.** `_coverage_bins_by_primer` builds one shared graph, reusing `BipartiteGraph`'s own region de-duplication. Test: `tests/test_stage2_keeps_the_sole_coverer.py`, whose behavioural case reproduced the drop before the fix.

**[verified here]**, reported by the hybrid module audit.

`CoverageRegion` (`dominating_set_optimizer.py:26-37`) is a `@dataclass` that
defines `__hash__` over `(chromosome, start, end)` but inherits the generated
`__eq__`, which compares all four fields including `covered_by`. Two regions for
the same genomic bin therefore hash equal and compare unequal:

```
hash equal : True
eq         : False
set size   : 2
```

`_coverage_bins_by_primer` (`hybrid_optimizer.py:920-947`) constructs a fresh
`BipartiteGraph` per primer, so each primer's regions carry `covered_by = {that
primer}`. The occupancy `Counter` at `:795-798` consequently records 1 for every
bin of every primer, and `unique_bins` at `:831` is the primer's total bin
count, never its uniquely covered bins. Measured on eight plasmid primers:

```
sum of per-primer bin counts: 50
distinct keys in Counter    : 50
keys with count > 1         : 0
```

Every primer's "unique" bin count equalled its total bin count.

The ranking at `:853` reads `norm_keep = 1 - (unique_bins - lo)/span`, so a low
count means cheap to remove. Because the count is now total bins, the primer
with the fewest bins is preferred for removal, which is typically the primer
solely covering a small region. The comment at `:779-793` states this term was
added to prevent a 94.0 to 64.9 per cent coverage regression; it currently
selects for that failure.

The module audit reproduced the end-to-end consequence: given six primers
binding identical positions plus one sole coverer of the right third of a 1 Mb
target, `_network_refine(target_count=5)` drops the sole coverer.

Fix: key the `Counter` on `(chromosome, start, end)`, or declare the dataclass
`eq=False`.

## 4 (CRITICAL) Background pruning is a no-op whenever coverage is at or below 0.95 -- FIXED

**[verified here]**, reported independently by two module audits.

`_prune_background` (`hybrid_optimizer.py:1164`) skips a removal when
`test_coverage < self.min_coverage_threshold`. That threshold is **absolute**
and defaults to 0.95, while the quantity compared is binned coverage at
realistic reach. If the set's coverage is already below 0.95, every candidate
removal is skipped, `best_removal` stays `None`, and the loop breaks on its
first pass.

Driving the real `_prune_background` with a controlled coverage function, 12
primers in, target 9:

| stage-1 coverage | primers kept | background sites |
|---|---|---|
| 0.99 | 9 | 7800 -> 4500 |
| 0.94 | 12 | 7800 -> 7800 |
| 0.80 | 12 | 7800 -> 7800 |
| 0.60 | 12 | 7800 -> 7800 |
| 0.39 | 12 | 7800 -> 7800 |

`--optimization-method=background-aware` sets `background_pruning=True` and
inherits the 0.95 default (`background_aware_optimizer.py:703`), so the method
documented for clinical use with "10-20x background reduction" performs no
pruning at all on a typical run, and degenerates to plain hybrid.

The regression is recent and self-documented. The comment at
`hybrid_optimizer.py:1005-1013` records that `_calculate_coverage` was corrected
from bin occupancy to realistic reach, and that the same set which read 100 per
cent under the old measure scores 39 per cent under the new one. That correction
was right, and it moved typical runs into the dead region of an absolute floor
that was not revisited. The comment even notes the floor "was not measuring what
it protected".

There is no escape hatch: `min_coverage_threshold` is settable only from
`neoswga iterate` (`cli/iterate.py:643`), not from `optimize`, and has no
params.json key.

**Correction to this finding.** An earlier draft said the fix was to make the
floor relative "which is what the docstring describes". The docstring does not
describe that. It says "Stops when coverage drops below
`min_coverage_threshold`", which matches the code exactly. Code and
documentation agreed; what was wrong was the default, which sits above the range
real runs occupy. That makes the fix a design decision rather than a repair of a
mismatch, and it is worth being precise about which of the two it is.

> **Fixed.** The floor is now relative to the coverage the stage was handed,
> with `min_coverage_drop` (default 0.05) as the budget, and
> `absolute_coverage_floor=True` still available for a caller that wants a hard
> bound. Across the five stage-1 coverage levels above, pruning now removes the
> three highest-background primers at every one of them (7800 background sites
> to 4500) instead of stopping dead at four of the five. Tests:
> `tests/test_hybrid_optimizer_run.py::test_the_coverage_floor_is_relative_to_where_pruning_started`
> and the strengthened `test_background_pruning_reduces_background_binding`,
> which now asserts a strict reduction on the pruning stage rather than a `<=`
> on the finished pipeline that a no-op also satisfied.

## 5 (CRITICAL) The network greedy's background penalty can become a reward -- FIXED

> **Fixed.** The background term counts host sites directly instead of diffing an average component size. The score now falls monotonically with host load (+2.0 at zero sites to +0.095 at twenty) where it previously rose to +42. Test: `tests/test_network_background_penalty_is_a_penalty.py`.

Reported by the network module audit, CONFIRMED there with a stubbed cache. Not
independently re-verified here.

`network_optimizer.py:920-923` computes `base_score = fg_improvement /
(bg_penalty + 1.0)`, where `bg_penalty` is the change in the background
network's **average** component size. Adding a primer with isolated host sites
lowers that average, so `bg_penalty` falls below -1 and the denominator turns
negative:

| isolated background sites added | `bg_penalty` | `base_score` |
|---|---|---|
| 0 | +0.00 | +2.000 |
| 1 | -2.50 | -1.333 |
| 13 | -4.58 | -0.558 |

Three consequences follow. Once the denominator is negative the score rises with
background load, so the greedy prefers the primer with more host binding. A
`bg_penalty` of exactly -1.0 is reachable, giving `+inf`, which is selected
immediately and can never be displaced. And because the Tm and dimer modifiers
are applied multiplicatively at `:842` and `:848`, a negative base score inverts
them, so the dimer penalty becomes a dimer bonus.

This is a sign error, distinct from the documented limitation that a
multiplicative soft penalty cannot cut more than its weight.

## 6 (HIGH) `background-aware` ignores the configured Tm window and `--max-extension` -- FIXED

> **Fixed.** The wrapper forwards the six dropped values. On the equiphi29 pool it now passes 500/500 candidates, matching `hybrid`. Test: `tests/test_background_aware_live_path.py`.

**[verified here]**, reported by the factory module audit.

`hybrid_optimizer.py:1284` forwards `min_tm`, `max_tm`, `max_extension`,
`genome_gc_content` and `bin_size` into the inner `HybridOptimizer`. The
background-aware wrapper at `background_aware_optimizer.py:692` wraps the same
class and forwards none of them, so unset bounds fall back to the polymerase
preset at `hybrid_optimizer.py:1064`.

One params file setting `min_tm: 20, max_tm: 70` on equiphi29, same candidate
pool:

| method | Tm window used | candidates passing the pre-stage |
|---|---|---|
| hybrid | 20.0-70.0 C | 500 / 500 |
| background-aware | **37.0-62.0 C** | **167 / 500** |

Two thirds of the user-approved pool is discarded before selection by the method
documented as the clinical one. `--max-extension` is dropped the same way, per
the table in finding 1. This bites only for equiphi29 and bst, the two
polymerases whose preset enables the stage-0 thermodynamic screen.

## 7 (HIGH) `--minimize-primers` with `ensemble` crashes and delivers no pool -- FIXED

> **Fixed.** `optimizer` is bound before the branch AND the ensemble builds a measuring optimizer, since passing None would have made the trim a silent no-op. Test: `tests/test_ensemble_survives_minimize_primers.py`.

**[verified here]**

`unified_optimizer.py:822` passes `optimizer=optimizer` to
`_minimize_primer_count`. The local `optimizer` is assigned only in the
single-method branch at `:791`; the ensemble branch at `:760-787` never assigns
it. Run from the CLI:

```
File "neoswga/core/unified_optimizer.py", line 822, in run_optimization
    optimizer=optimizer,
UnboundLocalError: cannot access local variable 'optimizer' where it is not associated with a value
```

Exit code 1, and no `step4_improved_df.csv` or summary JSON is written. Every
ensemble method has already run to completion at that point and all results are
discarded. The `auto` and `all` aliases reach the same branch.

## 8 (HIGH) The quality report divides a 0-1 dimer risk score by ten -- FIXED

> **Fixed.** Uses `max()` on the 0-1 scale with no division. The plasmid set now grades poor, and the report's overall grade moves from B to C. Tests rewritten to the real scale in `tests/report/test_quality.py`.

**[verified here]**

`report/quality.py:325-331`:

```python
worst_dimer = min((p.dimer_score for p in metrics.primers), default=0)
dimer_risk = min(abs(worst_dimer) / 10.0, 1.0)
```

The comment above it reads "dimer scores are typically negative (kcal/mol), more
negative = worse". The value received is not a free energy:
`report/metrics.py:169` fills `dimer_score` from the `dimer_risk_score` column,
a positive 0-1 risk score. Two errors compound, `min()` selecting the least
risky primer rather than the worst, and a division by ten applied to a value
already on the thresholds' own scale (`{'excellent': 0.1, 'good': 0.2,
'acceptable': 0.35, 'poor': 0.5}`).

| quantity | value |
|---|---|
| `dimer_risk_score` written by step 4 | 0.4667 |
| rating that value earns on its own scale | poor |
| value after `abs(worst)/10` | 0.0467 |
| rating the report prints | Excellent |

The error is systematic and always in the same direction, so every report
understates dimer risk.

## 9 (HIGH) No stage checks the delivered pool against `max_dimer_bp` -- FIXED

> **Fixed.** `optimize` now measures the delivered pool, records the worst pair in the summary, and raises a validation issue when it exceeds the threshold; `interpret` and `export` refuse to call such a pool ready. Test: `tests/test_delivered_pool_is_screened_for_dimers.py`.

**[verified here]**

Recomputing the longest contiguous complementary stretch over every pair of the
delivered set, with `max_dimer_bp: 3` configured:

| method | oligos | worst heterodimer |
|---|---|---|
| network | 12 | **10 bp** |
| network | 6 | **10 bp** |
| clique | 12 | 3 bp |
| clique | 6 | 3 bp |
| hybrid | 2 | 3 bp |
| dominating-set | 2 | 3 bp |
| background-aware | 2 | 3 bp |

On that same run,

- `neoswga interpret` prints "Overall rating: EXCELLENT" and "Ready for
  synthesis", with "Order primers for experimental validation" as step 1,
- `neoswga export --format fasta` prints "Primers ready for ordering!",
- `neoswga report` grades Dimer Risk "Low / Excellent", for the reason in
  finding 8.

The capability to catch this exists and works. `neoswga analyze-dimers` on the
identical six oligos returns `FAIL (Severe interaction detected (severity
1.00))`. It is a separate command taking sequences as arguments, and nothing in
the pipeline invokes it.

That a soft penalty cannot exclude a strongly scoring pair is documented and is
a design decision. The finding here is that nothing downstream reports the
result, and per finding 1 the penalty weight is zero on the default optimizer
anyway.

## 10 (HIGH) Multi-target coverage pools cross-genome coordinates -- FIXED

> **Fixed.** Coverage and gaps are measured per sequence and combined by length. Test: `tests/test_multi_target_coverage_is_not_pooled.py`.

Reported by the dispatch module audit, CONFIRMED there with two reproduction
scripts against the real `compute_metrics`.

`base_optimizer.py:1050-1084` extends one flat `fg_positions` list across every
foreground prefix, then takes `sorted(set(...))` and measures against
`sum(fg_seq_lengths)`. Each prefix is a separate HDF5 with its own 0-based
coordinate space, so position 5000 in target A and position 5000 in target B
collapse, and no position can exceed the longest single genome while the
denominator is the sum of all of them.

Two 100 kb targets, each individually 100 per cent covered:

| | one target | two targets |
|---|---|---|
| `fg_coverage` | 1.0000 | 0.5145 |
| `max_gap` | - | 102,000 on a 100,000 bp target |

`compute_per_prefix_coverage` on the same inputs returns `{'fgA': 1.0, 'fgB':
1.0}`, so the summary carries two contradictory figures and the incorrect one is
authoritative: `report/metrics.py:897` assigns `fg_coverage` into the report
headline. `--minimize-primers` gates on this value, so with a ceiling near 1/N
the 0.70 default is unreachable and the trim silently does nothing.

`hybrid_optimizer._calculate_coverage` gets this right by passing `prefix` into
`add_primer_coverage`; only the base class pools.

## 11 (HIGH) The amplification network joins sites on different molecules -- FIXED

> **Fixed.** Sequences occupy disjoint slots separated by more than the extension reach, in both the network builder and the scoring simulation. Test: `tests/test_amplification_network_respects_molecules.py`.

Reported by the network module audit, CONFIRMED there.

`network_optimizer.py:1216-1224` adds every prefix's positions into one
`AmplificationNetwork`, and `PositionCache` returns genome-local coordinates
with no offsetting. A forward site at 500 on contig C1 and a reverse site at 900
on contig C2 produce 2 nodes and 1 edge, an amplification edge across two
molecules that cannot share an amplicon.

`largest_component_size` is the entire `fg_improvement` term of the greedy
objective, so this changes selection, not only reporting, on any draft assembly,
plasmid-plus-chromosome target, or multi-target run.

## 12 (HIGH) The clique optimizer never receives its own configuration -- FIXED

> **Fixed.** The factory builds the declared `CONFIG_CLASS` and carries base fields across without mutating the caller's config. The four fields that were read by nothing are removed rather than made settable.

Reported by the network module audit, CONFIRMED there.

`_build_optimizer_config` (`unified_optimizer.py:380`) always returns a plain
`OptimizerConfig`, never a `CliqueOptimizerConfig`, so every clique-specific
field sits at its default. `max_pool_size` is 200, so with the default
`max_primer` of 500, 300 candidates are silently discarded, and
`_rank_and_truncate` ranks on raw foreground site count with no background term.
`CliqueOptimizerConfig` is instantiated only in
`tests/test_clique_optimizer.py`, so no test covers the production path.

This weakens the pool the guarantee is applied to, not the guarantee itself.

### The clique guarantee itself is real

Worth recording, since it is the recommended mitigation for finding 9. The
threshold is inclusive-allowed: `is_dimer_fast` sets `threshold = max_dimer_bp +
1` and reports a dimer only when the run reaches it (`dimer.py:95`). So a worst
heterodimer of exactly 3 bp at `max_dimer_bp=3` is the intended output, not a
near miss. `build_compatibility_graph` (`clique_optimizer.py:104-137`) drops
self-dimering primers as nodes and adds an edge only where the heterodimer check
returns 0, so every clique is mutually compatible by construction. When no
clique of the requested size exists it returns a clean failure with the message
"No dimer-free sets of size N found", rather than silently relaxing the
threshold. The two caveats are finding 12 above and finding 20 below.

## 13 (MEDIUM) `interpret` reads neither the summary nor the validation file -- FIXED

> **Fixed.** It reads the validation file, surfaces the findings verbatim, and will not recommend ordering a pool that breaks its own dimer threshold.

**[verified here]**

`results_interpreter.py:182-183` and `:410` open exactly three files:
`step4_improved_df.csv`, `step3_df.csv` and `params.json`. It never opens
`step4_improved_df_summary.json`, which CLAUDE.md names as authoritative, nor
`step4_improved_df_validation.json`.

On the plasmid run the optimizer detected that its headline metric was not
trustworthy and recorded it:

```json
{"level": "warning", "code": "coverage_saturated_on_small_genome",
 "detail": "pcDNA: genome=6157 bp, 6 primers x 2x 3000 bp reach = 36000 bp; coverage=100.0% is saturation-bounded (metric unreliable)"}
```

`neoswga report` renders that warning at the top of the page. `neoswga
interpret`, on the same directory, prints "Genome Coverage: 100.0% / Rating:
EXCELLENT" and recommends ordering.

## 14 (MEDIUM) The interpreter's dimer assessment looks for absent column names -- FIXED

> **Fixed.** The real column name is in the lookup. The fixture that fabricated a column the pipeline never writes is corrected.

**[verified here]**

`results_interpreter.py:384-390` searches for `dimer_score`, `dimer_risk` and
`heterodimer_score`. The column step 4 writes is `dimer_risk_score`, and the
lookup is an exact match. Called on the real CSV the function returns `None`
while the data holds 0.4667, so the assessment block at `:257-272`, including
the warning "High dimer risk may reduce amplification efficiency", is skipped
for every run.

This survived because the two relevant tests use different schemas.
`tests/test_results_interpreter.py:198` builds fixture rows carrying a
`dimer_score` column that the pipeline never produces, while
`tests/integration/test_output_column_parity.py:155` pins `dimer_risk_score` as
canonical. Both pass, and nothing compares them.

## 15 (MEDIUM) `--validate-simulation` cannot work -- FIXED

> **Fixed.** All three defects fixed: the two non-existent methods and the genome length read from a global that does not exist. The flag now completes.

**[verified here]**

`cli/pipeline.py:1018` and `:1033` call `AmplificationNetwork.build_graph()`.
The class has no such method; its builders are `add_primer_sites`,
`add_edges_for_sites` and `build_edges`. The `except Exception` at `:1060`
converts the `AttributeError` into a warning and the command exits 0:

```
WARNING: Validation failed: 'AmplificationNetwork' object has no attribute 'build_graph'
```

Two further defects in the same block are masked by the first. `add_site()` at
`:1030` sits inside a bare `except Exception: pass` and also does not exist, so
the networks would be empty regardless. And `:990-991` reads
`getattr(parameter, "fg_lengths", [1000000])`; the module global is
`fg_seq_lengths`, and `hasattr(parameter, "fg_lengths")` is False, so the
genome length would be a fabricated 1 Mb whatever the target.

`validate_simulation` is not in `cli/_common.py:UNIMPLEMENTED_OPTIONS`, the
registry the project already uses to warn about flags it cannot honour.

## 16 (MEDIUM) The ensemble tie-break is list order -- FIXED

> **Fixed.** Ranked by score, then the smaller set, then the method name. Reordering `--ensemble-methods` no longer moves the winner. Test: `tests/test_ensemble_tie_break_is_deliberate.py`.

**[verified here]**

`unified_optimizer.py:287-291` picks the winner with `max()` over the results
dict, which returns the first key attaining the maximum; insertion order is the
methods list order, beginning with `hybrid`. There is no tie-break on cost or
set size.

Ties are not a corner case. Three of the four default methods returned the
identical two-oligo set on the plasmid example, so their scores are equal to the
last digit in every profile:

| application | hybrid | dominating-set | background-aware | network | winner |
|---|---|---|---|---|---|
| balanced | 0.9003 | 0.9003 | 0.9003 | 0.8018 | hybrid |
| clinical | 0.9301 | 0.9301 | 0.9301 | 0.7981 | hybrid |
| discovery | 0.8795 | 0.8795 | 0.8795 | 0.7807 | hybrid |
| metagenomics | 0.9076 | 0.9076 | 0.9076 | 0.8003 | hybrid |

`--application` moves every score but never breaks the tie, and the log reports
"Ensemble winner: 'hybrid'" with no indication the choice was arbitrary.

That the tie-break is list order is directly demonstrable. Reordering
`--ensemble-methods` over the same three tied methods, same seed, same params,
moves the winner every time:

| `--ensemble-methods` order | winner |
|---|---|
| `hybrid dominating-set background-aware` | hybrid |
| `dominating-set background-aware hybrid` | dominating-set |
| `background-aware hybrid dominating-set` | background-aware |

All three rows score 0.9003 in every run. The winning name propagates into
`optimizer_name` in the summary, which the report reads, so the recorded
provenance of a design is set by flag order.

On this example the tied methods return the identical two primers, so only
provenance changes. The tie is on the composite rather than on the set, though,
so two methods returning *different* primers can tie and then list order picks
the delivered pool. CLAUDE.md's claim that the ensemble is order-independent is
false as written, and `_reseed` does not address it. Given the earlier
measurement that hybrid costs up to 260x dominating-set for an identical set at
128 primers, the default order also credits the most expensive member.

## 17 (MEDIUM) Three further inert options -- FIXED

> **Fixed.** `--min-fg-bg-ratio` now drives the background pre-filter threshold; `--no-position-cache` selects `StreamingPositionCache`, whose own k-range bug (`range(6, 13)`, invisible to every equiphi29 and bst primer) is fixed too; and `iterations` reaches the config, where the lookup used a name -- `max_iterations` -- that no global has. Tests: `tests/test_min_fg_bg_ratio_reaches_the_prefilter.py`, `tests/test_no_position_cache_uses_the_streaming_cache.py`, `tests/test_alternative_primer_sets.py`.

Reported by the dispatch module audit, CONFIRMED there by grep traces from
argparse to every read site.

- **`--min-fg-bg-ratio`** (`_optimize_parser.py:133-138`) is read only under
  `--auto-size` and `--show-frontier`. It is never forwarded into
  `optimize_step4`. The knob it should drive exists and cannot be set:
  `run_optimization` reads `kwargs.get("bg_min_ratio", 1.0)` at `:698` and
  `run_step4` never passes it.
- **`--no-position-cache`** sets `use_cache`, which `optimize_step4` accepts at
  `:1215` and never references. `run_step4` logs "Position cache: False" before
  forwarding it, which reads as confirmation the request was honoured.
- **`iterations` and `max_sets`** from params.json reach no optimizer.
  `max_iterations` is a bare `kwargs.get` at `:382` with no parameter fallback,
  and no optimizer reads `config.max_iterations` at all. `max_sets` is consumed
  only by `search_context.py`, which has no callers.

None of the three are in `UNIMPLEMENTED_OPTIONS`, so nothing warns.

## 18 (MEDIUM) The report prints the no-background sentinel as a measurement -- FIXED

> **Fixed.** Rendered as "no bg binding", matching `condition_sweep`.

**[verified here]**

`base_optimizer.py:33` defines `MAX_SELECTIVITY = 1e6`, and the docstring at
`:76` says it "stands in so the value stays finite" when there is no background
binding. `report/metrics.py:927` assigns `selectivity_density` straight to
`enrichment_ratio`, and the report renders "Enrichment 1000000x Excellent".
`condition_sweep.py:222` handles the same sentinel correctly, printing "no bg
binding".

## 19 (MEDIUM) `coverage_uniformity` rewards clustering -- FIXED

> **Fixed.** Measures the evenness of the gaps between components. Spread now scores 0.3156 against a cluster's 1.7680; it was 0.6546 against ~0. Test: `tests/test_coverage_uniformity_rewards_spread.py`.

Reported by the network module audit, CONFIRMED there on a 4 Mb synthetic
genome.

`network_optimizer.py:369-394` returns the coefficient of variation of component
centroid positions, and lower is treated as better. Eight binding islands spread
evenly across the genome score 0.6546; all eight packed inside one 5 kb island
score 0.0000. The `uniformity_weight` term therefore pushes selection toward
concentrating binding in one place. Per finding 1 it reaches selection only on
the network optimizer.

## 20 (MEDIUM) One lowercase primer defeats every dimer check -- FIXED

> **Fixed.** Candidates are normalised to uppercase at validation. Test: `tests/test_primer_case_is_normalised.py`.

**[verified here]**

`base_optimizer.py:912` validates with `primer.upper()` but appends the primer
unchanged, so lowercase survives into the pool. The dimer scan compares
characters directly, so a lowercase base never matches an uppercase one:

| pair | longest complementary run | `is_dimer_fast(thr=3)` |
|---|---|---|
| both uppercase | 10 | True |
| both lowercase | 10 | True |
| one lowercased | **0** | **False** |

A perfect 10 bp duplex reads as no dimer at all. This is the one input that
breaks the clique guarantee outright rather than weakening it. It is latent for
pipeline-generated candidates, which are uppercase, and reachable through fixed
primers, `expand-primers` and the library API.

Fix: normalise to uppercase in `_validate_candidates`.

## 21 (LOW) Smaller items -- MOSTLY FIXED

> **Mostly fixed.** Trimmed sets now report their own score; `--target-coverage 0` is honoured; `run_optimization_from_config` forwards the five fields it has a home for; the audit-trail hash covers the optimizer and the thresholds that change a set; `analyze-dimers`' contradictory summary and the three unexercised integration scenarios remain open.

All CONFIRMED by the module audits unless noted.

- **Trimmed sets keep the untrimmed score.** `unified_optimizer.py:528` replaces
  primers and metrics but leaves `score`, so a 4-primer-to-1-primer trim still
  reports the 4-primer score in the CSV and summary.
- **`--target-coverage 0` silently becomes 0.70**, via `or 0.70` at `:747`. Same
  shape for `-n 0` at `cli/pipeline.py:584`.
- **`run_optimization_from_config` drops nine of eighteen fields**
  (`:977-996`), including `minimize_primers`, `target_coverage` and
  `uniformity_weight`, contradicting the comment at `:740-745`.
- **`normalized_score` is absent from the summary JSON.**
  `base_optimizer.py:257-293` omits it, so the only cross-optimizer-comparable
  value is missing from the file the report reads. It is present in the CSV.
  Raw `score` across five methods on one dataset, all at `fg_coverage` 1.0:
  15.0, 15.0, 1.0, 1.0, 0.888. **[verified here]**
- **The audit trail cannot distinguish runs.** `_write_audit_trail` hashes a
  fixed 13-key list at `:1087-1101` that omits `optimization_method`, `seed`,
  `application` and `coverage_reach`.
- **`analyze-dimers` contradicts itself.** The same six-oligo set yields "FAIL
  (severity 1.00)" alongside "Total interactions: 0" and zero problematic
  partners for every primer, so the verdict is right and the detail cannot say
  which pair to break. `--output` is also treated as a directory.
  **[verified here]**
- **The three named integration scenarios never run the pipeline.**
  `tests/integration/{phi29_baseline,equiphi29_baseline,phi29_with_bg}/` contain
  only a `params.json`; the genomes they name are absent, and
  `check_genome_available()` is defined and never called. The two tests that run
  assert JSON keys, one of them pinning `optimization_method`, which Known Issue
  8 documents as inert. **[verified here]**

## Status of earlier findings

| Finding | Status |
|---|---|
| F0, `score` does not change the delivered set | Open. `unified_optimizer.py:641` reads the primer column only |
| F1b, `optimization_method` inert in params.json | Open. `_apply_params_only_keys` assigns five keys, not this one |
| F1, `coverage_reach` inert | Fixed. `parameter.py:850` assigns it; the optimizer resolves it |
| F2, `num_primers` above 50 rejected from params.json | Fixed. 12 from params.json and 12 from the flag both deliver 12 |
| `--use-mechanistic-model` silently ignored | Fixed. Wired at `cli/pipeline.py:704-776` and announced in the log |
| `--seed` not reproducible | Not reproduced. Three seeded and three unseeded runs returned identical sets |

## Recommendations, in order of payoff

1. **Close the configuration-routing hole** (1, 2, 6). Forward the four dropped
   parameters through the inner constructors; delete the polymerase rebind at
   `cli/pipeline.py:635`; forward the Tm window and `max_extension` through the
   background-aware wrapper. Extend
   `tests/test_optimizer_config_reaches_optimizers.py` to assert on the
   innermost object rather than the outermost, and add a case that runs the same
   params through every method and compares the effective configuration. Three
   of the audit's five critical and high routing defects would have been caught
   by that one test.
2. **Fix the two broken optimizer stages** (3, 4). Key the bin counter on
   coordinates, or set `eq=False` on `CoverageRegion`; make the pruning floor
   relative to starting coverage. Both are small changes that restore a stage
   which currently does nothing or the opposite of its intent.
3. **Repair the network background term** (5). A denominator that can be
   negative or infinite is a sign error, not a calibration question.
4. **Screen the delivered pool and say so** (8, 9, 14). Fix the `/10.0` and the
   `min()` in the report; fix the key lookup at `results_interpreter.py:386`;
   run the check `analyze-dimers` already performs at the end of `optimize` and
   record the worst pair in the summary. Order-readiness messages should not
   appear when the pool exceeds the configured threshold.
5. **Namespace coordinates per prefix** (10, 11). Both the coverage computation
   and the amplification network treat separate molecules as one.
6. **Fix the ensemble crash** (7), then break ties deliberately (16): prefer the
   smaller set, then the cheaper method, and report that it was a tie.
7. **Give `interpret` the same inputs as `report`** (13), and normalise primer
   case (20).
8. **Decide about the inert options** (15, 17). Either implement them or add
   them to `UNIMPLEMENTED_OPTIONS`, which is the mechanism the project already
   uses for this. Four more instances is enough to argue for a schema-driven
   check rather than a hand-maintained list.

## Method

Commands are relative to a copy of `examples/plasmid_example`.

```bash
neoswga optimize -j params.json --optimization-method=<method> --seed 42
neoswga optimize -j params.json --optimization-method=ensemble --application <profile> --seed 42
neoswga optimize -j params.json --optimization-method=ensemble --minimize-primers --seed 42
neoswga interpret -d ./
neoswga report -d ./
neoswga export -d ./ --format fasta
neoswga analyze-dimers --primers <set> -o <dir> --threshold 3
pytest -k "optimiz or dimer or interpret or ensemble or coverage" -n 8
```

Heterodimer stretches were recomputed independently of the package, by sliding
each primer against the reverse complement of its partner and taking the longest
run of matching bases. Constructor routing was checked by building each
optimizer through `OptimizerFactory` and reading the innermost
`NetworkOptimizer`'s attributes. Background pruning was exercised by calling the
real `_prune_background` with a controlled coverage function.
