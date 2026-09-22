# What each capability's evidence actually is

Generated from `capability_evidence.json`. Do not edit this page by hand; run
`scripts/generate_evidence_page.py`.

Every capability this tool advertises is placed on one of five tiers. The tiers
are ordered, and the line a reader should look for first is the one at the top:
**nothing here is prospectively validated**, because no pool this tool designed
has been synthesised, run and sequenced. [What NeoSWGA does not
establish](LIMITATIONS.md) says that first and in more detail.

The point of a tier is to stop "we built it" and "we measured it" from reading
alike. This repository has shipped a capability that was built, unit-tested,
documented and reachable by no command at all -- six of them at once -- which is
why `connected` is a tier of its own rather than an assumption.


## The tiers

**implemented** -- The code exists and its own unit tests pass. Nothing establishes that any command can reach it.

**connected** -- A command reaches it on a real path. `tests/test_no_capability_is_unreachable.py` and `tests/test_every_cli_option_has_an_effect.py` hold this: a capability built and tested but reachable by nothing has happened here six times at once.

**oracle-tested** -- Its output is checked against an independent implementation written for the purpose, not against itself. Coverage and the position scan both have one.

**retrospectively evaluated** -- Measured on a real reference, or against published primer sets with known outcomes, and the measurement is recorded with its caveats.

**prospectively validated** -- A pool this tool designed was synthesised, run and sequenced. NOTHING here is at this tier, and docs/LIMITATIONS.md says so first.


## Where each capability sits

| Capability | Tier | Reached by |
|---|---|---|
| Background-aware selection | retrospectively evaluated | `neoswga optimize --optimization-method background-aware` |
| Candidate filtering and the k-mer funnel | retrospectively evaluated | `neoswga count-kmers, neoswga filter` |
| Reaction chemistry and additives | retrospectively evaluated | `every command that resolves conditions` |
| Several optimization methods | retrospectively evaluated | `neoswga optimize --optimization-method` |
| Specificity against a host | retrospectively evaluated | `neoswga optimize with bg_genomes configured` |
| Binding position scan | oracle-tested | `neoswga filter` |
| Foreground coverage | oracle-tested | `neoswga optimize` |
| Adaptive GC filtering | connected | `neoswga filter, neoswga design` |
| Bloom filter for large backgrounds | connected | `neoswga build-filter, neoswga filter` |
| Export for ordering | connected | `neoswga export` |
| Position cache | connected | `every design command` |
| Quality reports | connected | `neoswga report, neoswga interpret` |
| Iterative design from sequencing depth | implemented | `neoswga analyze-coverage, neoswga expand-primers --bam` |
| Panel spacing statistics | implemented | `neoswga optimize, neoswga report` |

## The evidence, per capability

### Background-aware selection

*Adds a host-binding term to the pruning and to the refinement that chooses the panel.*

**Tier:** retrospectively evaluated. **Reached by:** `neoswga optimize --optimization-method background-aware`

Tests: `tests/test_background_aware_live_path.py`, `tests/test_expansion_uses_the_background.py`

Measurements: [pool_selection_audit_2026-09-18.md](validation/pool_selection_audit_2026-09-18.md)

Measured against hg38 on three GC-tier designs: host sites in the delivered panel fall 7-35% against hybrid and coverage falls 0.1-3.1 points. At n=12 it returns the same panel as hybrid. Reading the background is not acting on it, and three seams had to be closed before it did.

### Candidate filtering and the k-mer funnel

*Counts k-mers, applies frequency, thermodynamic, exclusion and evenness gates, and records the per-stage counts.*

**Tier:** retrospectively evaluated. **Reached by:** `neoswga count-kmers, neoswga filter`

Tests: `tests/test_gini_needs_enough_sites.py`, `tests/test_filter_funnel_stages.py`

Measurements: [occupancy_and_discrimination_2026-09-19.md](validation/occupancy_and_discrimination_2026-09-19.md)

The evenness gate's threshold was re-derived against delivered coverage. The discrimination profile the filter reports is a measurement, not a gate: an occupancy gate was tried and made panels worse.

### Reaction chemistry and additives

*SantaLucia nearest-neighbour melting temperatures with salt and additive corrections.*

**Tier:** retrospectively evaluated. **Reached by:** `every command that resolves conditions`

Tests: `tests/test_an_additive_stays_inside_its_recorded_temperature.py`, `tests/test_model_evidence_contract.py`

Measurements: [chemistry_model_evidence.md](validation/chemistry_model_evidence.md), [additive_specificity.md](validation/additive_specificity.md)

Per-constant evidence is recorded in neoswga/core/registry/model_evidence.json. Most additive coefficients are 37 C figures for PCR-length duplexes applied to 12-mers at 30 C, and only one of twenty-three records states a temperature range that can be checked mechanically. Glycerol has no melting-temperature term at all and a design that sets it is refused rather than given an invented coefficient.

### Several optimization methods

*hybrid, dominating-set, network, background-aware, clique, and ensemble.*

**Tier:** retrospectively evaluated. **Reached by:** `neoswga optimize --optimization-method`

Tests: `tests/integration/test_optimizer_method_coverage.py`, `tests/test_ensemble_tie_break_is_deliberate.py`

Measurements: [optimizer_cost_2026-09.md](validation/optimizer_cost_2026-09.md), [pool_selection_audit_2026-09-18.md](validation/pool_selection_audit_2026-09-18.md)

Measured against each other on real pools. hybrid returned a set identical to dominating-set at 7.8x the cost at 32 primers and 260x at 128. No method has been shown to produce a better experimentally performing pool, because no pool from any of them has been tested.

### Specificity against a host

*Reports selectivity ratio, selectivity density, host site counts and host coverage.*

**Tier:** retrospectively evaluated. **Reached by:** `neoswga optimize with bg_genomes configured`

Tests: `tests/test_expansion_counts_background.py`, `tests/validation/test_prevotella_specificity.py`

Measurements: [additive_specificity.md](validation/additive_specificity.md), [published_primer_sets.md](validation/published_primer_sets.md)

Read selectivity_density rather than selectivity_ratio when comparing designs scored against different backgrounds: the ratio has no genome length in it and moves about 66x between chr21 and whole hg38 with nothing about the primers changed. Only bg_coverage sees where host sites fall, and nothing selects on it.

### Binding position scan

*Finds every site of every candidate in the reference, on both strands.*

**Tier:** oracle-tested. **Reached by:** `neoswga filter`

Tests: `tests/test_positions_agree_with_an_independent_count.py`, `tests/test_position_cache.py`

Measurements: none.

Checked against a brute-force sliding window written in the test file, over overlapping occurrences, palindromes, ambiguous bases, a circular origin and record joins. The two scanners disagreed across a record join until 2026-09-21; the test that asserts they agree is the one that would have caught it.

### Foreground coverage

*Reports the fraction of the target within reach of a binding site.*

**Tier:** oracle-tested. **Reached by:** `neoswga optimize`

Tests: `tests/test_coverage_independent_oracle.py`, `tests/test_positions_agree_with_an_independent_count.py`

Measurements: [occupancy_coverage_rewrite_2026-09-17.md](validation/occupancy_coverage_rewrite_2026-09-17.md), [record_geometry_on_drosophila_2026-09-21.md](validation/record_geometry_on_drosophila_2026-09-21.md)

Checked base by base against an oracle written in the test file that calls nothing production calls. It is a geometric proxy at a modelled reach, not sequencing breadth, and the two production paths still disagree across a contig join by up to 4.8% on a fragmented assembly.

### Adaptive GC filtering

*Supports extreme GC genomes, 32-68% GC.*

**Tier:** connected. **Reached by:** `neoswga filter, neoswga design`

Tests: `tests/test_adaptive_gc_auto.py`, `tests/test_gc_adaptive_pipeline.py`, `tests/test_long_primer_gc_extremes.py`

Measurements: none.

Connected and tested, not measured against a design outcome. The window matters: a fixed 0.20 lower bound excluded exactly the zero-GC primers published AT-rich designs are built from, on a 19% GC target.

### Bloom filter for large backgrounds

*Handles a human-sized background genome.*

**Tier:** connected. **Reached by:** `neoswga build-filter, neoswga filter`

Tests: `tests/test_bloom_background_gate.py`

Measurements: none.

Often unnecessary: exact jellyfish counting of hg38 at k=12 costs about 7 minutes and a 138 MB table. The Bloom path matters at longer k. The sampled-index path has its own resolution trap.

### Export for ordering

*FASTA, vendor CSV, BED and BedGraph.*

**Tier:** connected. **Reached by:** `neoswga export`

Tests: `tests/test_export.py`, `tests/test_export_refuses_before_it_writes.py`

Measurements: none.

Refuses before it writes, as of 2026-09-22. It used to write every order file and then print that the pool was not ready, leaving a complete FASTA on disk beside the warning.

### Position cache

*About 1000x faster than re-reading the HDF5 files.*

**Tier:** connected. **Reached by:** `every design command`

Tests: `tests/test_position_cache.py`, `tests/test_positions_arrive_on_demand.py`

Measurements: [positions_on_demand_2026-09-17.md](validation/positions_on_demand_2026-09-17.md)

The speed claim is a timing, not a design outcome. Two defects here were silent zeros: int32 positions saturated past 2.1 Gb, and an unindexed prefix answered with an empty array.

### Quality reports

*A technical report badging each value MEASURED or ESTIMATED.*

**Tier:** connected. **Reached by:** `neoswga report, neoswga interpret`

Tests: `tests/test_report_agrees_with_the_saved_result.py`, `tests/test_design_report_provenance.py`

Measurements: none.

A test asserts every quantity the report renders equals the one the saved summary holds, for the exact exported panel. Writing that check found four gap metrics read with a default of 0.0, which is the BEST value for all four: a summary missing the key rendered as 'no coverage hole anywhere'.

### Iterative design from sequencing depth

*Adds oligos to a validated set using real BAM coverage gaps.*

**Tier:** implemented. **Reached by:** `neoswga analyze-coverage, neoswga expand-primers --bam`

Tests: `tests/test_bam_coverage.py`, `tests/test_bam_depth_lands_in_the_right_place.py`

Measurements: none.

The lowest tier here and the claim a reader is most likely to over-read. There is no BAM or CRAM in this repository, so nothing has run against measured depth: the tests use constructed alignments. Contig binding is by name only as of 2026-09-21, because equal length is not identity and this repository ships two plasmids of 5,386 bp.

### Panel spacing statistics

*Reports worst gap, mean gap, gap evenness and strand balance.*

**Tier:** implemented. **Reached by:** `neoswga optimize, neoswga report`

Tests: `tests/test_strand_scores_say_when_they_are_unmeasurable.py`

Measurements: [getting_ahead_on_spacing_2026-09-18.md](validation/getting_ahead_on_spacing_2026-09-18.md), [published_primer_sets.md](validation/published_primer_sets.md)

Computed and reported; no optimizer ranks on them, deliberately. No threshold derived from the polymerase reach separates the 18 published sets with wet-lab outcomes, the winners included, so the tool will not pick one. You can impose your own through the panel limits.
