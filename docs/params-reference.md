# params.json Reference

This document is generated from `neoswga/core/schema/params.schema.json`.
Do not edit by hand — run `python scripts/render_schema.py` after
changing the schema.

To dump the schema JSON directly:

```bash
neoswga schema --dump > params.schema.json
```

## Required parameters

| Parameter | Type | Range / allowed | Default | Description |
|---|---|---|---|---|
| `data_dir` | string | - | - | Directory for pipeline outputs (CSVs, HDF5 caches, logs). |
| `fg_genomes` | array of string | minItems: 1 | - | Foreground (target) genome FASTA paths, one per target. |
| `fg_prefixes` | array of string | - | - | K-mer file prefixes for each foreground genome. Populated by count-kmers. |

## Optional parameters

| Parameter | Type | Range / allowed | Default | Description |
|---|---|---|---|---|
| `adaptive_gc` | boolean | - | `True` | Enable genome-GC-aware filtering. Set false to force gc_min/gc_max to user values. |
| `betaine_m` | number | min: 0.0; max: 2.5 | `0.0` | - |
| `bg_circular` | boolean | - | `False` | - |
| `bg_genomes` | array of string | - | - | Background (host / off-target) genome FASTA paths. |
| `bg_prefixes` | array of string | - | - | - |
| `bg_seq_lengths` | array of integer | - | - | - |
| `bl_genomes` | array of string | - | - | Blacklist genomes; primers matching these are penalized or rejected. |
| `bl_penalty` | number | min: 0.0; max: 100.0 | `5.0` | Penalty weight applied to blacklist matches. |
| `bl_prefixes` | array of string | - | - | - |
| `bl_seq_lengths` | array of integer | - | - | - |
| `bloom_filter_path` | string or null | - | - | - |
| `bsa_ug_ml` | number | min: 0.0; max: 400.0 | `0.0` | - |
| `cpus` | integer | min: 1; max: 128 | - | - |
| `dmso_percent` | number | min: 0.0; max: 10.0 | `0.0` | - |
| `drop_iterations` | integer | min: 0; max: 100 | - | - |
| `ethanol_percent` | number | min: 0.0; max: 5.0 | `0.0` | - |
| `excl_genomes` | array of string | - | - | Zero-tolerance exclusion genomes (e.g., mtDNA). |
| `excl_prefixes` | array of string | - | - | - |
| `excl_threshold` | integer | min: 0 | `0` | - |
| `fg_circular` | boolean | - | `True` | - |
| `fg_seq_lengths` | array of integer | - | - | Per-genome foreground lengths in bp. Auto-computed if missing. |
| `formamide_percent` | number | min: 0.0; max: 10.0 | `0.0` | - |
| `gc_max` | number | min: 0.0; max: 1.0 | - | - |
| `gc_min` | number | min: 0.0; max: 1.0 | - | - |
| `gc_tolerance` | number | min: 0.0; max: 0.5 | `0.15` | - |
| `genome_gc` | number | min: 0.0; max: 1.0 | - | - |
| `glycerol_percent` | number | min: 0.0; max: 15.0 | `0.0` | - |
| `iterations` | integer | min: 1; max: 100 | `8` | - |
| `long_primer_mode` | boolean | - | `False` | - |
| `max_bg_freq` | number | min: 0.0; max: 1.0 | - | - |
| `max_bl_freq` | number | min: 0.0; max: 1.0 | `0.0` | Maximum permissible blacklist frequency; 0 = zero tolerance. |
| `max_dimer_bp` | integer | min: 1; max: 15 | - | - |
| `max_gini` | number | min: 0.0; max: 1.0 | - | - |
| `max_k` | integer | min: 4; max: 30 | - | Maximum primer length (bp). Polymerase-aware default is used if absent. |
| `max_primer` | integer | min: 1; max: 10000 | - | - |
| `max_self_dimer_bp` | integer | min: 1; max: 15 | - | - |
| `max_sets` | integer | min: 1; max: 100 | `5` | - |
| `max_tm` | number | min: 0.0; max: 100.0 | `45.0` | - |
| `mg_conc` | number | min: 0.0; max: 20.0 | - | Mg2+ concentration (mM). Polymerase-aware default is used if absent. |
| `min_amp_pred` | number | - | - | - |
| `min_fg_freq` | number | min: 0.0; max: 1.0 | - | - |
| `min_k` | integer | min: 4; max: 30 | - | Minimum primer length (bp). Polymerase-aware default is used if absent. |
| `min_sample_count` | integer | min: 1 | - | - |
| `min_tm` | number | min: 0.0; max: 100.0 | `15.0` | - |
| `na_conc` | number | min: 0.0; max: 1000.0 | `50.0` | - |
| `num_primers` | integer | min: 1; max: 50 | `6` | - |
| `optimization_method` | string | one of: hybrid, greedy, network, genetic, milp, dominating-set, background-aware, moea | - | - |
| `peg_percent` | number | min: 0.0; max: 15.0 | `0.0` | - |
| `polymerase` | string | one of: phi29, equiphi29, bst, klenow | `phi29` | - |
| `primer_conc` | number | min: 1e-09; max: 0.0001 | `5e-07` | - |
| `reaction_temp` | number | min: 20.0; max: 70.0 | - | - |
| `retries` | integer | min: 0; max: 100 | - | - |
| `sample_rate` | number or null | min: 0.001; max: 1.0 | - | - |
| `schema_version` | integer | min: 1; max: 1 | `1` | Version of this schema the file was written for. |
| `selection_metric` | string | one of: deterministic, random, stochastic | - | - |
| `src_dir` | string | - | - | Source directory; usually equal to data_dir. |
| `target_set_size` | integer | min: 1; max: 50 | - | - |
| `tmac_m` | number | min: 0.0; max: 0.1 | `0.0` | - |
| `top_set_count` | integer | min: 1; max: 100 | - | - |
| `trehalose_m` | number | min: 0.0; max: 1.0 | `0.0` | - |
| `urea_m` | number | min: 0.0; max: 2.0 | `0.0` | - |
| `use_bloom_filter` | boolean | - | `False` | - |

---

## Additional guidance

- See `docs/production-scenarios.md` for scenario-specific recipes
  (phi29 plasmid, equiphi29 bacterium, extreme GC, multi-genome, etc.).
- See `docs/SWGA_SCIENCE.md` for theoretical background.
- Use `neoswga validate-params -j params.json` to check a configuration
  against all validation layers (schema + ranges + interdependencies).

