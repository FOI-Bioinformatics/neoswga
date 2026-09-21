# Training the amp_pred random forest on real lab data

The shipped `random_forest_filter.skops` was trained on synthetic data via
`scripts/retrain_rf_model.py`: features sampled at random, labels supplied by
`compute_target_score`, a hand-written rule over Tm, GC, GC clamp and
homopolymer runs. To replace it with a model fit to measured SWGA outcomes, use
the workflow below.

## Overview

```
lab CSV (per-primer or per-set measurements)
       |
       v   scripts/lab_csv_to_training_data.py
       |   (computes 120+ features per primer using AdvancedFeatureEngineer)
       v
training_data.csv
       |
       v   scripts/train_enhanced_rf.py
       v
enhanced_rf_model.pkl
       |
       v   predict_new_primers_enhanced(model_path=...) in Python -- see Step 4;
           there is no CLI flag that wires this into `neoswga prepare-candidates` yet
```

## Step 1: collect lab measurements

The bridge accepts two CSV layouts.

**Per-primer** (one primer measured at a time):

```csv
primer,enrichment_fold
ATCGATCGATCGAT,12.3
GCTAGCTAGCTAGC,4.1
ATCGATCGATCGAT,15.0
```

Repeated rows for the same primer are averaged.

**Per-set** (one primer set measured as a unit; primers semicolon-separated):

```csv
primers,enrichment_fold
ATCGATCGATCGAT;GCTAGCTAGCTAGC,8.5
GCTAGCTAGCTAGC;AAAGGGCCCTTTAA,11.2
```

The set's enrichment is attributed to each primer in the set, then averaged
per primer across sets. The signal is noisier than per-primer measurements;
plan for more samples to compensate.

`enrichment_fold` is the measured target/background DNA ratio after SWGA
(e.g., from sequencing read counts or qPCR). Replicates can be encoded as
duplicate rows.

## Step 2: compute features

```bash
python scripts/lab_csv_to_training_data.py \
    --lab-csv lab_results.csv \
    --fg-genome ecoli_reference.fna \
    --params params.json \
    --output training_data.csv \
    --mode per-primer
```

`--params` reuses your pipeline params.json so reaction conditions
(`reaction_temp`, `na_conc`, `mg_conc`, additives, polymerase) match the
conditions under which the lab data was generated. Individual flags
(`--polymerase`, `--reaction-temp`, ...) override params.json fields.

If lab data was collected under multiple reaction conditions, generate a
training_data.csv for each condition and concatenate, or extend the script
to read condition columns from the lab CSV.

## Step 3: train

```bash
python scripts/train_enhanced_rf.py \
    --training-data training_data.csv \
    --output-model enhanced_rf_model.pkl \
    --target-variable log_enrichment \
    --cross-validate
```

Useful flags:

- `--n-features 50`: select top features by mutual information. Recommended
  for small training sets (<200 primers) to reduce overfitting.
- `--tune-hyperparams`: grid search over RF hyperparameters. Adds runtime;
  worth using once when the dataset is reasonably large.

The trainer reports R-squared, RMSE, and feature importance. Cross-validated
R-squared below ~0.2 indicates the model is not learning a useful signal;
consider collecting more data or revisiting feature selection before using
the model in production.

## Step 4: use the trained model

There is currently no CLI path to this: `neoswga prepare-candidates`'s `--use-enhanced-features`
and `--enhanced-model-path` flags are accepted but explicitly unimplemented
(`neoswga/cli/_common.py:UNIMPLEMENTED_OPTIONS`; the `score` step's help text
says so too), and the real `score` step (`core/pipeline.py`) calls only
`rf_preprocessing.predict_new_primers`, never the enhanced path. The
`is_enhanced_model_available()` / `get_enhanced_model_info()` helpers exist
but are wired only into `rf_preprocessing.py`'s own `__main__` block -- a
manual status check, not the prediction call the `score` step makes.

The underlying prediction function works and takes a model path directly; use
it from Python instead of the CLI:

```python
from neoswga.core.rf_preprocessing import predict_new_primers_enhanced

result = predict_new_primers_enhanced(
    primer_list=primers,              # list[str]
    fg_genome_sequence=genome_seq,    # str, the foreground genome
    conditions=reaction_conditions,   # ReactionConditions
    primer_positions=positions,       # dict, as produced by PositionCache
    model_path="enhanced_rf_model.pkl",
)
```

This falls back to the standard (synthetic-trained) model if the enhanced
model can't be loaded, and to `predict_new_primers(df)` if `df` is supplied
and `use_enhanced=False`. Wiring `--use-enhanced-features` /
`--enhanced-model-path` through to this function in the `score` step is
outstanding work, not yet done.

## Sample-size guidance

With 120+ features, honest evaluation needs on the order of 500+ labeled
primers (or sets). With fewer than ~100 samples, prefer:

- `--n-features 30` to `--n-features 50` to limit overfitting.
- Per-primer mode over per-set mode (lower label noise).
- Reporting cross-validated R-squared, not training R-squared.

## When to retrain

Retrain when any of the following change materially:

- Reaction conditions (polymerase, temperature, additive concentrations).
- Target genome class (e.g., switching from gram-negative to mycobacterial
  GC-rich targets).
- Primer length range (the synthetic-trained model is unreliable above
  ~12 bp; a real-data-trained model inherits the length range of its
  training set).

## Limitations

- Per-set labels assign the same enrichment value to every primer in the
  set. This conflates the contribution of individual primers; the model
  will learn a "primer that tends to appear in good sets" signal, not a
  per-primer enrichment signal.
- The bridge computes features on the foreground genome only. If you want
  background-aware features, extend the bridge to also pass background
  positions to `AdvancedFeatureEngineer`.
- Reaction conditions are taken as constant across the training set. Mixing
  conditions in one training run requires extending the feature set with
  condition columns and is not handled by the current bridge.
