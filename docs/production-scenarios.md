# Production Scenarios

Reference configurations for the scenario matrix that NeoSWGA supports
end-to-end. Each scenario ships as a ready-to-edit template in
`examples/`. Copy a template, replace the FASTA paths with your own
genomes, and run the four-step pipeline.

## Quick lookup

| Scenario                                    | Example template                      | Polymerase  | Typical primer length | Additives                 |
|--------------------------------------------|---------------------------------------|-------------|----------------------|---------------------------|
| Small plasmid (5-10 kb), phi29              | `examples/plasmid_example/`           | phi29       | 6-12 bp              | none                      |
| Bacterial chromosome, EquiPhi29, long primers | `examples/equiphi29_scenario/`       | equiphi29   | 10-18 bp             | DMSO 5%, betaine 1 M      |
| Pan-genome + blacklist                      | `examples/multi_genome_blacklist/`   | phi29       | 8-12 bp              | none                      |
| Extreme GC (AT-rich or GC-rich)             | `examples/gc_extreme/`                | equiphi29   | 12-18 bp             | betaine 1.5 M, DMSO 5%, TMAC 0.05 M |
| Human host background (Gb scale)            | add `use_bloom_filter: true` to any scenario | -   | -                    | -                         |

## Scenario 1: Small plasmid, phi29

Starting point for most users. Short primers, standard phi29 at 30 C,
no additives required. Works for plasmids, viral genomes, and small
bacterial targets up to ~1 Mb.

```bash
cd examples/plasmid_example
neoswga count-kmers -j params.json
neoswga filter -j params.json
neoswga score -j params.json
neoswga optimize -j params.json
neoswga report -d .
```

Expected output: `step4_improved_df.csv` with 5-10 primer sets and an
HTML report in the current directory.

## Scenario 2: Bacterial chromosome with EquiPhi29

For targets where phi29's 30 C is too permissive (e.g. GC-rich bacteria,
strong secondary structure). EquiPhi29 at 42-45 C with longer primers
improves specificity.

Key differences from scenario 1:
- `polymerase: "equiphi29"`
- `reaction_temp: 43.0` (anywhere in 42-45 C works)
- `min_k: 10, max_k: 18`
- `dmso_percent: 5.0, betaine_m: 1.0`

Template: `examples/equiphi29_scenario/params.json`

## Scenario 3: Pan-genome design with blacklist

Primers that amplify multiple target strains (`fg_genomes: [a, b, c]`),
avoid a single host (`bg_genomes: [host]`), and have zero tolerance for
a known contaminant (`bl_genomes: [contam]`).

Key parameters:
- `fg_genomes`: list of targets (pan-primer design)
- `bg_genomes`: host or off-target sequences (frequency tolerated)
- `bl_genomes`: zero-tolerance avoid list
- `max_bl_freq: 0.0` for strict rejection, or `max_bl_freq: 1e-6` for
  low-penalty weighted filtering

Template: `examples/multi_genome_blacklist/params.json`

## Scenario 4: Extreme-GC genomes

Targets outside the 35-65% GC range (*Plasmodium* ~20%, *Mycobacterium*
~65%, *Burkholderia* ~67%). The adaptive GC filter engages automatically
when `genome_gc` is outside that band and narrows the primer GC window to
`genome_gc +/- gc_tolerance`. Additives help further:
- Betaine equalises AT/GC kinetics
- TMAC tightens Tm across GC-extreme sequences
- DMSO reduces secondary structure

Template: `examples/gc_extreme/params.json`

## Scenario 5: Human / Gb-scale host background

For targets from a host background too large to load into memory (e.g.
human DNA in a clinical sample). Pre-build a Bloom filter once and reuse
it across runs:

```bash
neoswga build-filter GRCh38.fna ./bloom
```

Add to params.json:
```json
{
  "use_bloom_filter": true,
  "bloom_filter_path": "./bloom/background.pkl",
  "sampled_index_path": "./bloom/sampled.pkl"
}
```

The Bloom filter estimates background k-mer counts without storing the
full genome, trading a small false-positive rate for a large memory
reduction.

## Validating your params.json before running

Before starting a long pipeline, validate the configuration:

```bash
neoswga validate-params -j params.json
```

This runs type, range, and interdependency checks, including the new
canonical JSON schema, and tells you exactly which parameters are wrong.

To inspect the schema directly (useful for IDE autocomplete):

```bash
neoswga schema --dump > params.schema.json
```

## Picking a polymerase

| Polymerase | Temp      | Primer length | When to use                                        |
|-----------|-----------|---------------|---------------------------------------------------|
| phi29     | 30 C      | 6-12 bp       | Standard SWGA; small genomes, low complexity       |
| equiphi29 | 42-45 C   | 10-18 bp      | GC-rich or structured targets; need specificity    |
| bst       | 60-65 C   | 15-25 bp      | LAMP-like applications; thermostable               |
| klenow    | 25-40 C   | 8-15 bp       | Room temperature; lower processivity               |

If `polymerase` is set in params.json but `min_k`/`max_k`/`mg_conc` are
not, NeoSWGA uses polymerase-appropriate defaults (6-12 bp / 10 mM for
phi29, 10-18 bp / 10 mM for equiphi29, etc.).

## Which optimizers see my additives?

NeoSWGA pipelines apply ReactionConditions (polymerase, temperature,
DMSO / betaine / etc.) at two points:

1. **At filter time (always).** `filter` applies Tm windows, adaptive GC
   ranges, and the additive-adjusted Tm when deciding which k-mers make
   it into `step2_df.csv`. Every downstream optimizer picks from this
   already-screened pool, so the reaction buffer shapes the universe of
   possible primers for every method.

2. **At selection time (some optimizers).** Five optimizers re-apply
   ReactionConditions inside their selection loop to weight Tm score and
   edge stability under the actual buffer:

   - `network`
   - `hybrid` (default; the inner network stage)
   - `background-aware` (three-stage; its final network stage)
   - `genetic`
   - `equiphi29` (delegates to hybrid)

   The remaining registered optimizers pick from the filter-screened pool
   by their own set-cover / coverage objective without re-applying
   conditions:

   - `greedy`
   - `dominating-set`
   - `weighted-set-cover`
   - `tiling`
   - `clique`
   - `milp`

   This is by design. Pure coverage optimisers trust that the filter step
   already respected the user's reaction buffer, and keep their
   coverage-optimality guarantees intact. Run `neoswga doctor` to see
   the live capability matrix — the `additive-aware: yes/no` column is
   declared on each optimizer class (it is not a heuristic) and
   accurately reflects which path shapes primer selection.

## Multi-target coverage: union vs. intersection

When `fg_prefixes` lists more than one target genome, the coverage
metric `fg_coverage` aggregates as the **union** across prefixes —
i.e. a primer that binds only target A still contributes to the
global coverage number, because SWGA amplifies any target it touches.
`per_target_coverage` (populated by every optimizer since Phase 15A)
is the per-target dict `{prefix: covered_fraction}` so users can see
the imbalance directly.

| Scenario | What you see | What it means |
|---|---|---|
| `fg_coverage = 0.95`, `per_target = {A: 0.95, B: 0.95}` | Balanced | Primer set amplifies both targets well. |
| `fg_coverage = 0.90`, `per_target = {A: 0.95, B: 0.40}` | Target B underserved | Aggregate is misleading; target B is under-covered. A `per_target_coverage_below_threshold` warning fires when `min_per_target_coverage` > B's coverage. |
| `fg_coverage = 1.00`, `per_target = {plasmid: 1.00}` on a 5 kb plasmid | Saturation | Phase 17B emits `coverage_saturated_on_small_genome` because N primers x 2 x per-primer-reach (~6-18 kb) exceeds the genome length; the 100% number is not a quality signal. |

Neoswga does **not** compute an intersection coverage metric (how much
of the genome is covered in *every* target simultaneously). For
pan-primer designs that require "this primer must work on all targets",
inspect `per_target_coverage` directly and reject sets where any entry
is below your threshold. The HTML report surfaces the
`per_target_coverage_below_threshold` warning visibly (Phase 17A) so
this does not require reading the JSON.

## Diagnosing "no primers selected"

The most common failure is empty filter output. In order of likelihood:

1. `max_bg_freq` too strict — relax to `1e-5`.
2. Adaptive GC window too narrow for an extreme-GC genome — widen
   `gc_tolerance` or set `adaptive_gc: false` and pin your own window.
3. `min_k`/`max_k` not in polymerase's sweet spot — remove these keys
   to use polymerase-aware defaults.
4. Blacklist too strict — relax `max_bl_freq` if legitimate primers are
   being rejected.

Step 2 now fails fast with actionable suggestions when filtering
produces zero candidates.
