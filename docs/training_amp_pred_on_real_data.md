# Training the amp_pred random forest on real lab data

The shipped `random_forest_filter.p` was trained on synthetic data via
`scripts/retrain_rf_model.py`. To replace it with a model fit to measured
SWGA outcomes, use the workflow below.

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
enhanced_rf_model.pkl  ->  drop into neoswga/core/models/
       |
       v   neoswga score -j params.json --use-enhanced-model
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

## Step 4: install and use

```bash
cp enhanced_rf_model.pkl neoswga/core/models/
neoswga score -j params.json --use-enhanced-model
```

`--use-enhanced-model` switches the score step to the trained model. The
default model path is resolved by `is_enhanced_model_available()` in
`neoswga/core/rf_preprocessing.py`. To use a model from a non-default
location, pass `--enhanced-model-path /path/to/model.pkl`.

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
