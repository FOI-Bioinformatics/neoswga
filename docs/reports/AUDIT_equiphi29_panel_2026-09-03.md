# EquiPhi29 panel design audit, 2026-09-03

Scope: the pipeline path that produced the three GC-tiered EquiPhi29 12-mer
panels in `runs/gc_tiers/` (S. aureus 32.9% GC, E. coli 50.8%, M. tuberculosis
65.6%), scored against whole hg38. Each finding below was re-derived from the
source or the run outputs rather than carried over from notes.

## F1. The dimer penalty saturates, so the optimizer cannot avoid the worst pair

`network_optimizer.calculate_dimer_score` (neoswga/core/network_optimizer.py:40)
returns 1.0 when the longest complementary stretch exceeds `max_dimer_bp` and
0.0 otherwise. It is a step function, not a graded measure. The value feeds
`_dimer_penalty` (neoswga/core/network_optimizer.py:745) as the only dimer term
in selection.

At the M. tuberculosis panel's configured `max_dimer_bp` of 3, the penalty is
saturated across most of the pool:

| max_dimer_bp | pairs scoring 1.0, of 630 |
|---|---|
| 3 (configured) | 328 (52%) |
| 4 | 185 (29%) |
| 5 | 82 (13%) |
| 6 | 48 (8%) |

A pair sharing 4 complementary bases and a pair sharing 11 receive the same
penalty of 1.0. Over half the pairwise space sits at the ceiling, so the term
carries no gradient in the region where it would matter, and the optimizer has
no signal distinguishing an incidental overlap from a near-complete duplex.

This is not the inert-parameter defect of Known Issue #8. `max_dimer_bp` does
reach the optimizer, through `unified_optimizer.py:386`, and
`tests/test_optimizer_config_reaches_optimizers.py` pins that routing. The
parameter arrives and is then applied as a threshold on a binary score.

Fixed. `calculate_dimer_score` now grades with the length of the complementary
stretch: zero at or below `max_dimer_bp`, 1.0 for a full-length duplex, and a
linear ramp between. On the same panel the ceiling is now empty, the penalty
takes nine distinct values, and the worst pair scores 0.889 against 0.111 for a
four-base overlap. Test: `tests/test_dimer_penalty_is_graded.py`.

### F1b. The configured threshold never reached the network optimizer

Found while verifying the fix above. `NetworkBaseOptimizer` builds its inner
`NetworkOptimizer` from an explicit argument list that omitted `max_dimer_bp`,
so the object that scores pairs used its own default of 4 whatever params.json
said. The panel was configured for 3 and scored at 4, which is the looser
threshold, so pairs the user asked to treat as dimers were counted as clean.
The run log said `max bp: 4` next to a params.json reading 3.

This is the same class as Known Issue #8, one constructor deeper than the
existing routing test reached: the value did arrive on the `OptimizerConfig`
and was then dropped on the way in. Fixed, and pinned by a new case in
`tests/test_optimizer_config_reaches_optimizers.py`, which is where this class
of defect is already tracked.

## F2. The M. tuberculosis panel contains an 11 bp heterodimer

The 36-primer panel in `runs/gc_tiers/high_mtb/step4_improved_df_summary.json`
was produced by the `network` optimizer. Recomputing the longest contiguous
complementary stretch over all 630 pairs gives:

| stretch (bp) | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
|---|---|---|---|---|---|---|---|---|---|---|
| pairs | 71 | 231 | 143 | 103 | 34 | 32 | 10 | 4 | 1 | 1 |

The worst pair is `TCGTCAGACCCA` against `GGGTCTGACGAC`, sharing an 11 bp duplex
of a 12 bp oligo. The two are close to reverse complements of each other. The
next is `GTCGGCGACGAC` against `CGTCGCCGACCA` at 10 bp. Worst self-dimer in the
panel is 4 bp, so the self-dimer filter at step 2 is working as intended; the
gap is specific to cross-primer pairs.

This is the expected consequence of F1 rather than a separate fault. The panel
should be re-screened and the 11 bp pair broken before the oligos are ordered.

### The fix does not remove the pair, and cannot

Re-running selection on a copy of the run, with both F1 and F1b fixed, returns
the identical 36 primers. Coverage, selectivity and the duplex distribution are
unchanged, and the 11 bp pair is still there. Raising the dimer weight to the
highest profile value (`--application clinical`, 0.45) also changes nothing.

The reason is the shape of the penalty rather than its calibration. It enters
as `score * (1 - weight * factor)`, so the largest haircut it can ever apply is
the weight itself. At the decisive step the offending primer scored 23.77
before any penalty, and the graded factor of 0.296 at weight 0.30 reduced that
by 8.9 per cent. No weight in the permitted range closes a gap of that size
against the runner-up.

So grading the penalty restores the gradient and makes the term informative,
which is worth having, but a multiplicative soft penalty cannot exclude a
strongly scoring pair on its own. Excluding one requires a hard constraint: the
clique optimizer, which guarantees a dimer-free set, or a screen applied to the
finished panel. That is a design decision rather than a defect, and it is not
made here.

## F3. The multi-genome path disables heterodimer checking outright

`thermodynamic_filter.filter_candidates` accepts `check_heterodimers`, default
True. `neoswga/core/multi_genome_pipeline.py:370` calls it with False. The
pairwise cost is quadratic in pool size, so the choice is understandable, but it
means the multi-genome path has no heterodimer screen at any stage.

## F4. Step 3 scoring assumes one oligo length per batch

`rf_preprocessing` derives k from the first element of the primer list and uses
it to preload k-mer files for every worker. The code says so in a comment,
"assumes uniform length in batch". A mixed-length batch would score every primer
whose length differs from the first against the wrong k-mer dictionary, with no
error raised. Step 3 passes the whole of step2_df as one batch, so any future
mixed-length design would need the batch split by length first.

This is latent rather than active: the shipped runs are single-length, and
mixed-length pools are not viable at EquiPhi29 anyway, because a 10-mer melts at
39-42 C and cannot clear a 42 C reaction temperature.

## F5. The heuristic scoring fallback is biased against high-GC targets

`_heuristic_primer_score` awards its maximum for GC between 0.40 and 0.60 and
for Tm between 30 and 42 C, and penalises outside those bands. Optimal EquiPhi29
12-mers for a 65% GC genome fall outside both. The heuristic runs only when the
random forest model fails to load, but in that case it would rank the better
candidates lowest on exactly the high-GC designs where the margin is tightest.
The windows should follow the genome GC and the reaction conditions the rest of
the pipeline already computes.

## F6. Three GC classifiers still disagree

`gc_adaptive_strategy` splits AT-rich from balanced at 0.35 and balanced from
GC-rich at 0.65. `condition_suggester` and `genome_analysis` use 0.40 and 0.60.
A genome between 0.35 and 0.40, or between 0.60 and 0.65, is classified two ways
within one run. Already recorded in docs/TECH_DEBT_2026-09.md. It does not
affect the three panels here, since all three sit outside the disputed bands.

## Checked and not reproduced

The GC-adaptive strategy was suspected of silently widening the k-mer range past
what params.json asks for. It does not. `pipeline.py:522` tests whether min_k or
max_k came from the JSON, preserves them when they did, and logs the decision
either way. When neither is set it clamps the strategy's range to the
polymerase registry range rather than replacing it.
