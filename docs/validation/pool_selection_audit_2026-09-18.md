# Auditing how a pool gets selected

Audited 2026-09-18. Two questions: whether the selection logic is as informed as
the measurements a design is accepted on, and whether an additive can be made to
admit more candidates that are both more specific and better amplified.

**The additive mechanism works as hypothesised.** Measured on the model, a
modest additive moves GC-rich candidates into the band where they both prime and
discriminate, and those are the same candidates that are compositionally rare in
an AT-rich host. **Nothing in selection can exploit it**, for two specific
reasons that are each a few lines of wiring.

Reproduce with `scripts/benchmarking/candidate_gate_audit.py`.

## Where reaction conditions reach, and where they stop

The reaction reaches more of this pipeline than of any published SWGA tool, and
it does not reach the stage that chooses the panel.

| Stage | Decides | Condition-aware |
|---|---|---|
| Tm gate (`filter`) | admission | yes, additive-corrected effective Tm, but see below |
| shortlist cut to `max_primer` | admission | yes, occupancy-weighted background load, `occupancy_ranking` defaults true |
| evenness gate (`max_gini`) | admission | no, site positions only |
| **Stage 1 set cover** | **panel composition** | **no, unweighted coverage bins** |
| Stage 2 `_network_refine` (default) | panel composition | partly: Tm-weighted edges, gap CV at weight 0.10, background as a count divisor |
| Stage 2 `_swap_refine` (`plan-pool` only) | panel composition | yes, occupancy-weighted coverage under an occupancy-weighted specificity floor |
| `clique` method | panel composition | no, raw exact-match counts until the final re-rank |
| acceptance and report | what the user is told | yes |

## What an additive actually does

An additive lowers effective Tm. Occupancy at the reaction temperature and
mismatch discrimination move in opposite directions along that axis, so the
useful primer is the one near its transition: bound often enough to prime, not
so tightly that a mismatch costs nothing.

Per GC class, 8,000 random 12-mers at 35% GC, equiphi29 at 42 C, one mismatch
costed at 4 C. Classes below G+C 2 and above 8 are omitted as too small to read. `theta` is occupancy, `disc` is the per-candidate ratio
`theta(perfect) / theta(one mismatch)`, and the last column is their product,
a crude joint figure of merit rather than a validated one.

| G+C of 12 | n | Tm plain | theta | disc | Tm +additive | theta | disc | product plain to +additive |
|---|---|---|---|---|---|---|---|---|
| 2 | 886 | 32.7 | 0.025 | 5.83 | 30.2 | 0.009 | 6.08 | 0.14 to 0.05 |
| 3 | 1,599 | 36.1 | 0.096 | 5.44 | 33.3 | 0.031 | 5.95 | 0.52 to 0.18 |
| 4 | 1,882 | 40.0 | 0.321 | 4.32 | 36.8 | 0.117 | 5.50 | 1.39 to 0.65 |
| 5 | 1,641 | 43.7 | 0.658 | 2.67 | 40.1 | 0.330 | 4.42 | 1.75 to 1.46 |
| **6** | 1,022 | 47.4 | 0.888 | 1.54 | 43.4 | 0.627 | 2.91 | **1.37 to 1.82** |
| **7** | 431 | 51.1 | 0.975 | 1.12 | 46.7 | 0.864 | 1.69 | **1.09 to 1.46** |
| 8 | 164 | 54.9 | 0.996 | 1.02 | 50.2 | 0.967 | 1.17 | 1.02 to 1.13 |

The additive shifts Tm down 2.5 to 4.7 C. That degrades every class at or below
GC 5 and improves every class at or above GC 6, moving the best class up one
step and making the GC-rich tail usable for the first time. GC 7 gains 34%.

This is the user-facing claim, confirmed: **an additive makes GC-rich candidates
more specific and better amplified at once.** GC-richness is also how
compositional selectivity presents in a real pool, where the most selective
Wolbachia candidates run GC 0.44 to 0.48 against 0.39 for random ones
([no_search_headroom_on_this_pool_2026-09-18.md](no_search_headroom_on_this_pool_2026-09-18.md)).

### The two routes to specificity conflict, and only one is modelled in selection

They are anti-correlated, strongly.

| Route | Favours | Where it appears |
|---|---|---|
| Compositional: be rare in the host | GC-rich against an AT-rich host | `selectivity_density`, `occupancy_ranking` |
| Thermodynamic: sit near the transition | AT-rich at fixed length and temperature | occupancy-weighted metrics only |

Pearson correlation between GC and thermodynamic discrimination, over the same
same 8,000 12-mers: **about -0.89**. A single scalar gate cannot serve both, which is why
the additive matters: it is the one lever that moves a candidate along the
thermodynamic axis without changing its composition. That is the whole of the
mechanism, and it is why the best design measured in
[additive_specificity.md](additive_specificity.md) is an additive design at
k = 12 rather than a longer-primer one.

### The Tm gate is the wrong parameterisation

The useful band above is GC 4 to 6, roughly theta 0.3 to 0.9. The shipped window
admits GC 3 through 8, which spans theta 0.10 to 0.995. It is wide in both
directions and for the same reason: it measures Tm against a per-polymerase
constant rather than against the temperature the reaction runs at.

| Polymerase | Default Tm window | Typical reaction temp | Floor minus temp |
|---|---|---|---|
| phi29 | [20, 50] | 30 | -10 |
| equiphi29 | [37, 62] | 42 | -5 |
| bst | [50, 75] | 63 | -13 |
| klenow | [20, 55] | 37 | -17 |

A primer on the floor has occupancy 0.002 to 0.014 at phi29 and 0.054 to 0.128
at equiphi29. A primer well below the ceiling is saturated and discriminates
nothing. Neither is excluded, and there is no occupancy gate anywhere;
`occupancy_ranking` is the closest thing and it ranks on the background side,
not on whether the candidate is bound to the target.

Gating on an occupancy band instead would be additive-aware and
temperature-aware by construction, with no per-polymerase constant to maintain.
Measured against the shipped gate over the same samples, a band of
`0.3 <= theta <= 0.8` raises mean per-candidate discrimination at every length
(3.16 against 2.99 at k = 12, and 4.29 against 1.45 at k = 16) and, unlike the Tm
gate, admits MORE candidates under the additive than without it at every length
from k = 14 up. It also admits far fewer candidates in absolute terms at long k,
which is a real cost rather than a detail.

### Why longer primers are not the answer

Under the shipped Tm window, whether an additive enlarges the candidate pool or
shrinks it is decided by which bound is active, and the sign flips at k = 15.

| k | Admitted, plain | Admitted, +additive | Change | Fraction above ceiling |
|---|---|---|---|---|
| 12 | 0.699 | 0.514 | -0.185 | 0.001 |
| 14 | 0.954 | 0.913 | -0.040 | 0.009 |
| **15** | 0.970 | 0.980 | **+0.010** | 0.022 |
| 16 | 0.957 | 0.992 | +0.035 | 0.043 |
| 18 | 0.863 | 0.970 | +0.108 | 0.137 |
| 20 | 0.707 | 0.902 | +0.196 | 0.293 |

At k = 18 the additive costs nothing and admits candidates at mean GC
0.506 against 0.320 for those admitted either way. That looks like the goal and
it is not: at k = 18 the admitted pool is saturated, mean discrimination is 1.09
plain and 1.19 with the additive, and occupancy spread across the pool is 1.0.
**The regime that gains candidates is the regime with no discrimination left to
buy.** An earlier draft of this note recommended k >= 15 on the strength of the
table above, before the discrimination column was measured.

So the tool's short-primer defaults are right for the lever, and the useful
change is to the gate rather than to the length.

## Stage 1 is the least informed stage, and it is the deciding one

`hybrid_optimizer.py:711` calls `optimize_greedy` with no `objective`. So do
`dominating_set_adapter.py:164` and `primer_expansion.py:652`. The parameter is
threaded only inside `dominating_set_optimizer.py` and into
`swap_refinement.py`. **No production path passes it.** The only caller that
does is `tests/test_partial_panel_pruning.py:214`.

Without it, `_select_next_primer` scores a candidate by
`len(primer_regions - covered_regions)`: coverage bins nothing selected already
touches, every site at weight 1.0. Occupancy, background load and evenness are
absent, and so is the background tie-break, which is gated on the same
parameter.

Commit `59a4ee3` added the parameter and opens "The greedy now chooses on the
quantity the design is judged on", with the failure mode stated exactly: "a
primer with many sites and a melting temperature well below the reaction
temperature touches many bins and contributes little amplification". That is the
defect. The remedy was written and never connected.

This is the class of Known Issues 8 and 14 reached a fourth way, and the
ratchets cannot see it: `test_no_capability_is_unreachable.py` walks reach to
FUNCTIONS, and `optimize_greedy` is reachable. A parameter no caller supplies is
invisible to it.

### How much the blindness can cost

Occupancy depends only on the primer, so an unweighted bin count misranks two
candidates by exactly the ratio of their occupancies. Over 12,000 random 12-mers
at 35% GC, restricted to those the Tm gate admits:

| Platform | Admitted | theta p5 | median | p95 | p95/p5 | Fraction below 0.5 |
|---|---|---|---|---|---|---|
| phi29 30 C | 11,049 | 0.544 | 0.988 | 1.000 | 1.8 | 0.042 |
| equiphi29 42 C | 8,354 | 0.128 | 0.613 | 0.991 | 7.8 | 0.415 |
| equiphi29 + DMSO 5% + betaine 1 M | 6,156 | 0.117 | 0.408 | 0.965 | 8.3 | 0.586 |
| bst 63 C, k = 18 | 10,341 | 0.0007 | 0.020 | 0.787 | 1149 | 0.898 |

The cost tracks the platform where additives work. On phi29 occupancy is
saturated and the unweighted count is nearly correct, which is the same reason
phi29 offers no discrimination. On equiphi29 at k = 12, the configuration the
additive work settled on, 42% of admitted candidates are bound less than half
the time and the ranking can be wrong by a factor of 8. Under the additive that
rises to 59%. At k = 18 the spread is 1.0 and the blindness costs nothing,
which is the same saturation that removes the lever.

## Specificity and evenness are never both in force

Two selection regimes, neither carrying both.

`optimize` and `design` construct no `PoolObjective` at all: the only
constructors are `pool_planner.py:404` and `cli/plan_pool.py:143`. Their Stage 2
is `_network_refine`, scoring `fg_improvement / (1 + bg_added)` blended with a
gap-CV improvement at `uniformity_weight`, which `--application` sets to 0.10 by
default. Evenness is present at a small weight; there is no specificity
constraint, only a divisor.

`plan-pool` sets `refinement_method="swap"`, where `_score` is
`(-shortfall, coverage, -total_bg_sites)`. The specificity floor is enforced and
evenness appears nowhere.

Two different quantities are both called evenness: `gap_gini` in the post-hoc
`normalized_score`, and `coverage_uniformity`, a coefficient of variation of
gaps, in selection. `max_gap` is computed, interpreted in prose by
`coverage.interpret_gap_metrics`, and never constrained or optimized.

**This is less damaging than it looks.** Maximising a union of fixed-width
windows implicitly rewards spreading, because overlap is wasted, and it keeps
rewarding it until the genome saturates. At realistic panel sizes saturation is
far off: 32 primers at about 8 sites each with 3 kb reach can reach at most half
of a 3.2 Mb target. The implicit mechanism is visible in the measured numbers,
where max gap falls 108 kb to 46 kb as a panel grows from 8 to 32 primers while
coverage rises 0.382 to 0.539. What union maximisation cannot control is the
tail: one large gap costs little union coverage once the rest is covered, and
max gap decides whether a region is recoverable at all.

## Background structure is measured and never used

`selectivity_density` and `total_bg_sites` are additive in per-primer site
counts. `occupancy.weighted_site_load` sums `count * theta` per mismatch class
and no position enters it. Two backgrounds with identical per-primer counts
therefore score identically whether their sites are clustered or dispersed.

That distinction is most of off-target amplification. SWGA amplifies where two
sites sit in convergent orientation within the polymerase's reach; 100 host
sites spread over 3 Gb support almost none, and 20 inside one 10 kb window
support a great deal. The accounting scores the first as worse.

`bg_coverage` is the one computed quantity that sees position, and nothing reads
it. It is written to the metrics dataclass and serialised into
`step4_improved_df_summary.json`, and it is absent from `normalized_score`, from
`PoolObjective` and from every optimizer's scoring. No background amplification
network is built anywhere, in contrast to the foreground network the hybrid and
network methods build at about 70 kb processivity.

## Multi-target selection does not balance across targets

`per_target_coverage` is populated after optimization
(`unified_optimizer.py:1121`) and consumed only as a post-hoc validation warning
when `min_per_target_coverage` is set. No selection stage sees it and there is
no minimax term, so on aggregate coverage a panel covering two targets 0.9 and
0.1 beats one covering them 0.5 and 0.5.

## The clique method is the least condition-aware

It is the only method that guarantees a dimer-free set, which is a real property
the others lack. Its scoring is also the weakest here:

- `_rank_and_truncate` cuts the candidate list to 200 by **raw foreground site
  count alone**, with no specificity, Tm or occupancy term.
- `_score_set` ranks enumerated sets by `fg_sites / (bg_sites + 1)`, exact
  matches only, no genome length, no occupancy, no position.
- `normalized_score` re-ranks the top 100, the first point at which any
  condition enters.

An additive cannot move either of the first two by any amount.

## Compared with the published tools, on selection logic

[tool_comparison.md](tool_comparison.md) compares capabilities. This compares
the selection rule, which is a narrower and less flattering question. Sources
are listed at the end of this section.

| Criterion | swga 1.0 (2017) | swga 2.0 / SoapSWGA (2023) | COATswga (2025) | NeoSWGA |
|---|---|---|---|---|
| Set search | max clique, `cliquer` branch and bound | breadth-first greedy with a drop-out step | greedy interval tiling via bedtools | greedy set cover, then swap or network refinement; clique available |
| Dimer-free set | hard, by construction | pairwise check while growing | hard, by construction | hard in `clique` only; penalised in the others |
| Foreground site spacing | **hard**, `max_fg_bind_dist`, and "even binding site spacing" is a stated set criterion | **fitted**, `on_gap_gini` | implicit in tiling, with `target_coverage` default 0.95 | reported only; gap CV at weight 0.10 in one refinement mode |
| Background site spacing | **hard**, `min_bg_bind_dist`, and the clique vertex weight is the mean background binding distance | **fitted**, `off_gap_gini` | not modelled | **not modelled** |
| Convergent-orientation amplicon geometry | no | yes, `coverage_ratio` at 70 kb, computed for both genomes and taken as a ratio | `fragment_length`, default 10,000, target only | foreground network at about 70 kb; **no background network** |
| Set score fitted to wet-lab data | no | yes, ridge regression on 46 published SWGA and sequencing sets | no | no, hand-chosen weights |
| Reaction conditions | mono- and divalent cation Tm correction | Tm window, default [15, 45], no correction stated | Tm window, default [15, 45], none stated | additives, salt, and site occupancy at the reaction temperature |

Three readings follow.

**NeoSWGA is alone in modelling the reaction, and alone in constraining neither
spacing criterion.** swga 1.0 makes both hard: a floor on background site
spacing and a ceiling on foreground site spacing, with cliques ranked by total
background binding distance. swga 2.0 replaces both with fitted terms. NeoSWGA
computes the foreground quantities and constrains none of them, and does not
compute the background one at all. The oldest tool here is stricter about
exactly the two things this audit found unconstrained.

**The one externally fitted opinion about what makes an SWGA pool work assigns
background evenness a large weight.** swga 2.0's five set-level terms are
`freq_ratio`, `mean_gap_ratio`, `coverage_ratio`, `on_gap_gini` and
`off_gap_gini`, fitted by ridge regression with 10-fold cross-validation against
the proportion of the target genome reaching 1x sequencing coverage, over 46
published sets. `off_gap_gini` is the Gini index of distances between background
binding sites, per strand. The in-repo record of the coefficients puts it second
largest at +0.281, behind `freq_ratio` at +0.321; a positive weight means an
UNEVEN background is better, which is what the amplicon geometry predicts, since
clustered host sites leave most of the host unreachable. I confirmed the term
list, the response variable, the sample size and the method from the paper, and
could not re-extract the coefficient table, so treat the individual values as
recorded rather than verified. One secondary reading of the table disagreed with
the recorded sign of `mean_gap_ratio`.

**On coverage geometry NeoSWGA and COATswga converge, and COATswga is stricter.**
Both union fixed-width windows around binding sites and select to tile the
target: `fragment_length` against `coverage_reach`. COATswga takes dimer-freedom
as a hard constraint on every set it forms and carries an explicit
`target_coverage` of 0.95. NeoSWGA's equivalent guarantee exists in one
non-default method.

### Where NeoSWGA is genuinely ahead

These are not restatements of the capability table; each is a selection-logic
difference.

- **Occupancy-weighted coverage and selectivity.** No other tool weights a
  binding site by how much of the time it is bound at the reaction temperature.
  Every one of the three uses a Tm window and then counts sites at weight 1.0,
  which is the defect this audit finds in NeoSWGA's Stage 1 and which the
  others have everywhere. Two of the three use the same [15, 45] window
  irrespective of the enzyme.
- **Additives at all.** swga 1.0 corrects Tm for mono- and divalent cations.
  Neither of the later two states any correction. None models DMSO, betaine or
  formamide, so none can express the lever measured above.
- **Separating coverage reach from processivity**, and fitting the reach to
  sequencing depth.
- **A specificity floor as a constraint** rather than a scoring term, in
  `plan-pool`.
- **Multi-genome targets**, though without a worst-target term.

### Sources

- swga 1.0: [Bioinformatics 33(14):2071](https://academic.oup.com/bioinformatics/article/33/14/2071/3056525),
  [github.com/eclarke/swga](https://github.com/eclarke/swga)
- swga 2.0 / SoapSWGA: [PLOS Comput Biol 19(4):e1010137](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1010137),
  [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2022.04.27.489632v1.full)
- COATswga: [bioRxiv 2025.11.26.688640](https://www.biorxiv.org/content/10.1101/2025.11.26.688640.full.pdf),
  [github.com/bailey-lab/coatswga](https://github.com/bailey-lab/coatswga)

## What would change a delivered panel

Ordered by expected effect per unit of work, and none of it is measured against
a re-derived panel yet.

1. **Pass the objective to Stage 1.** The parameter, its documentation and its
   tie-breaks already exist. On equiphi29 at k = 12 it is choosing between
   candidates whose contributions differ by up to 8x.
2. **Add an occupancy gate beside the Tm window.** A floor alone removes the
   unbound tail; a band also removes the saturated head, which is where the
   additive lever is spent. This is the change that makes "more candidates that
   are more specific" expressible, and it is the only one of these five that no
   published tool has, since all three of them gate on a Tm window and two use
   the same [15, 45] window whatever the enzyme.
3. **Let the two regimes share their terms.** A specificity floor on `optimize`,
   or an evenness term in the swap score, whichever the design calls for.
4. **Compute background gap evenness, or read `bg_coverage`.** This is the
   only item here with external support: `off_gap_gini` carries the second
   largest weight in the one set-level model fitted against measured sequencing
   breadth, and swga 1.0 made background site spacing a hard filter in 2017.
   `bg_coverage` already exists and is unread, so the cheap version is to use
   it; the faithful version is a Gini of background inter-site distances.
5. **A worst-target term for multi-genome designs.**

## What this does not establish

- Every candidate-gate measurement is over random sequence at fixed GC, not over
  a genome's distinct k-mers. It answers what fraction of sequence space a gate
  admits, which is a property of the gate and not a pool composition.
- No delivered panel was re-derived. The occupancy spread bounds how far Stage
  1's ranking CAN be wrong; it does not show a shipped panel is worse than an
  objective-aware Stage 1 would return. An earlier attempt to improve Stage 1 on
  the Wolbachia pool failed four ways
  ([stage_one_constraint_awareness_2026-09-18.md](stage_one_constraint_awareness_2026-09-18.md)),
  and at panel size 12 no dimer-respecting construction beat the shipped search.
- The mismatch cost is a flat 4 C per mismatch, and discrimination is linear in
  that constant. A different penalty moves every `disc` figure here.
- The occupancy band 0.3 to 0.8 is a plausible range, not a fitted one. Nothing
  measured where the band should sit.
- `theta x disc` is a product of two modelled quantities offered as a reading
  aid. It is not a validated figure of merit and should not be optimized against
  without one.
- The k = 15 crossover is specific to equiphi29's [37, 62] window, 35% GC and
  this additive pair. The mechanism generalises; the number does not.
- Background clustering blindness is established from the formula, not from a
  case where it changed a panel.
- The comparison rests on published methods sections and repository
  documentation, not on running the other three tools. Defaults quoted for them
  are what their documentation states. The swga 2.0 coefficient values are
  carried over from [tool_comparison.md](tool_comparison.md) and were not
  re-verified here; the term list, response variable, sample size and regression
  method were.
- Nothing here compares a NeoSWGA panel with a panel from another tool on the
  same input. The comparison is of selection rules, not of outputs.
- Coverage, occupancy and selectivity all remain modelled site geometry rather
  than measured amplification. See
  [evidence_matrix_2026-09-15.md](evidence_matrix_2026-09-15.md).
