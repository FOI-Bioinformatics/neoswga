# Migration 3.6 -> 3.7

NeoSWGA 3.7 is a production-readiness release. Existing 3.6 configurations
continue to run unchanged; this note documents the new behaviour you
should be aware of.

## What changed

### 1. Polymerase-aware defaults

If you omit `min_k` / `max_k` / `mg_conc` from `params.json`, NeoSWGA now
fills them in based on `polymerase`:

| Polymerase | `min_k, max_k` | `mg_conc` |
|-----------|----------------|-----------|
| phi29     | 6, 12          | 10 mM     |
| equiphi29 | 10, 18         | 10 mM     |
| bst       | 15, 25         | 8 mM      |
| klenow    | 8, 15          | 10 mM     |

**Action**: if your 3.6 run relied on the old `mg_conc` default of 2 mM,
pin `"mg_conc": 2.0` explicitly. The 3.6 value was suboptimal for every
supported polymerase; the new 10 mM default follows vendor guidance.

### 2. Adaptive GC filter auto-engages

When `genome_gc` is set (or auto-computed) and sits outside `[0.35, 0.65]`
(i.e. extreme-GC genomes like *Plasmodium* or *Mycobacterium*), the
pipeline now switches to `gc_min = genome_gc - gc_tolerance` and
`gc_max = genome_gc + gc_tolerance` automatically. A single INFO log
line announces this.

**Action to opt out**: set `"adaptive_gc": false` in `params.json`, or
pin your own `gc_min` / `gc_max`.

### 3. Blacklist length auto-computation

`bl_seq_lengths` is now auto-computed when `bl_genomes` is present but
`bl_seq_lengths` is missing or wrong. Prior runs that set `bl_genomes`
via `--blacklist` on the CLI had no `bl_seq_lengths` and produced
incorrect blacklist frequencies in edge cases.

**Action**: none. If you had an explicit `bl_seq_lengths` that matched
your genomes it continues to be used.

### 4. Additives now affect scoring

`IntegratedQualityScorer` previously used Wallace's rule for Tm, ignoring
DMSO / betaine / trehalose / formamide / ethanol / urea / TMAC. Scores
(and therefore top-K rankings) now reflect the effective Tm your reaction
will actually see. Same primer pool + same targets + different additives
now produces different top-K rankings.

**Action**: none for new runs. If you want to reproduce a 3.6 ranking
exactly, set all additive concentrations to 0 in `params.json`.

### 5. ParamValidator is stricter

`formamide_percent`, `ethanol_percent`, `urea_m`, `tmac_m`,
`glycerol_percent`, `peg_percent`, `bsa_ug_ml`, `mg_conc`, `gc_tolerance`,
`genome_gc`, `bl_penalty`, and `max_bl_freq` now have explicit ranges.
Previously-silent out-of-range values now emit validation errors.

A canonical JSON schema ships at `neoswga/core/schema/params.schema.json`
and is now loaded automatically when `jsonschema` is installed (it is a
core dependency as of 3.7).

**Action**:
- Dump the schema for IDE autocomplete:
  ```bash
  neoswga schema --dump > params.schema.json
  ```
- Run `neoswga validate-params -j params.json` if you see unexpected
  errors; the message will point you at the offending field.

### 6. Multi-genome `[0]` hardcoding fixed

Two places in the pipeline only used the first foreground / background
genome in multi-genome runs:
- `pipeline.py` position-file cache check (now all `fg_prefixes` are
  required to have caches before the optimization skips position files).
- `background_aware_optimizer.compare_optimizers()` (now aggregates hits
  across all `bg_prefixes`).

**Action**: none. If you run with one target and one background, results
are unchanged.

### 7. New CLI commands

- `neoswga schema [--dump] [-o FILE]` — dump the canonical JSON schema.

## New example templates

- `examples/equiphi29_scenario/` — long primers + DMSO + betaine
- `examples/multi_genome_blacklist/` — pan-primer with zero-tolerance
  contaminant list
- `examples/gc_extreme/` — AT-rich or GC-rich targets

See `docs/production-scenarios.md` for a scenario-by-scenario walkthrough.

## Release infrastructure

3.7 also adds:
- `.github/workflows/publish.yml` — PyPI trusted publishing on release
- `.github/workflows/nightly.yml` — nightly end-to-end pipeline tests
- `.pre-commit-config.yaml` — ruff / black / isort / yaml / toml hooks
- `ruff` and `mypy` baselines in CI (non-blocking) and `pyproject.toml`

## Upgrading

```bash
pip install --upgrade neoswga
```

No action is required for pipelines that:
- Explicitly set `min_k`, `max_k`, `mg_conc`, `gc_min`, `gc_max`
- Use only phi29 or equiphi29 without additives
- Use a single target and single background
