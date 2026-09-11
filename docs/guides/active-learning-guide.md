# Active learning with lab results

This guide explains how to feed wet-lab measurements back into NeoSWGA so
subsequent primer-set designs learn from real amplification performance.

The active-learn subsystem uses Gaussian Process regression over primer-set
features (network connectivity, coverage, Tm statistics) trained on
measured outcomes (enrichment fold, coverage uniformity) to propose the
next set most likely to succeed.

## The data format

Each experimental run is an `ExperimentalResult` record. The JSON payload
used by `neoswga active-learn --experimental-results` is an array of
these objects:

```json
[
  {
    "primer_set": ["ATCGATCGATCG", "GCGTAGCATAGC", "TATACGCATGGA"],
    "timestamp": "2026-04-15T14:30:00Z",
    "enrichment_fold": 12.4,
    "uniformity_score": 0.18,

    "total_amplification": 1.8,
    "off_target_fraction": 0.07,

    "temperature": 30.0,
    "time_hours": 4.0,
    "primer_concentration": 1.0,

    "replicate_cv": 0.12,
    "passed_qc": true,
    "notes": "Standard phi29, DMSO 5%, betaine 1 M"
  }
]
```

### Field reference

**Required**

| Field | Type | Description |
|---|---|---|
| `primer_set` | list[str] | Primer sequences tested as a set. |
| `timestamp` | str | ISO-8601 timestamp of the experiment. |
| `enrichment_fold` | float | Measured target / background read ratio (higher is better). |
| `uniformity_score` | float | Coefficient of variation across target regions (lower is better, usually 0-1). |

**Optional**

| Field | Type | Default | Description |
|---|---|---|---|
| `total_amplification` | float | null | DNA yield (ng/uL or fold). |
| `off_target_fraction` | float | null | Fraction of reads mapping off-target. |
| `temperature` | float | 30.0 | Amplification temperature (C). |
| `time_hours` | float | 4.0 | Amplification time. |
| `primer_concentration` | float | 1.0 | Primer concentration (uM). |
| `replicate_cv` | float | null | Coefficient of variation across replicates. |
| `passed_qc` | bool | true | Whether the run passed sequencing QC. |
| `notes` | str | "" | Free-text. Include buffer/additive summary if relevant. |

## From lab CSV to active-learn JSON

Most wet-labs keep results in a spreadsheet or CSV, not JSON. The
`scripts/csv_to_experimental_results_json.py` helper converts a CSV like

```csv
primer_set_id,primers,enrichment_fold,uniformity,coverage,temperature,notes
set_A,ATCGATCGATCG;GCGTAGCATAGC;TATACGCATGGA,12.4,0.18,0.72,30.0,phi29 + 1 M betaine
set_B,ATCGATCGATCG;GCGTAGCATAGC;CAGCTAGTCATT,6.1,0.31,0.58,30.0,phi29 no additive
```

into the ExperimentalResult JSON accepted by `active-learn`:

```bash
python scripts/csv_to_experimental_results_json.py \
    --input results.csv \
    --output experimental_results.json
```

The `primers` column is a semicolon-separated list. Empty optional
columns become `null`.

## End-to-end iterative workflow

```bash
# 1. Design initial set
neoswga design -j params.json

# 2. Synthesize the top primer set; run SWGA + sequencing in the lab.

# 3. Tally outcomes into a CSV (one row per primer set tested).

# 4. Convert CSV to the JSON format active-learn expects.
python scripts/csv_to_experimental_results_json.py \
    --input results.csv \
    --output experimental_results.json

# 5. Ask active-learn to propose the next experiment, informed by
#    the lab measurements.
neoswga active-learn -j params.json \
    --experimental-results experimental_results.json \
    --num-candidates 10 \
    --output round_2/

# 6. Repeat. Each round tightens the GP model's predictions and lets
#    later rounds converge toward sets the model is confident about.
```

## Calibration diagnostics

`ExperimentalTracker` (wired into `neoswga predict-efficiency --track`)
records the model's prediction for each primer set before synthesis and
the measured outcome after sequencing. After a few rounds, call the
tracker's calibration report to see whether the model's enrichment and
coverage estimates correlate with reality:

```python
from neoswga.core.experimental_tracker import ExperimentalTracker
tracker = ExperimentalTracker.load("./tracker.json")
report = tracker.get_calibration_report()
print(report)  # MAE for each metric + recommended scaling factor
```

Systematic over- or under-prediction is a signal to adjust either the
scoring weights or the mechanistic-model activation energies (see
`mechanistic_params.py`).

## Common pitfalls

- **Primer ordering** inside `primer_set` matters only insofar as
  downstream analyzers use set-level features (sorted internally), but
  keeping order stable across rounds aids reproducibility.
- **`passed_qc: false`** records are still used by the GP, just
  weighted lower. Remove them manually only if the experimental
  failure is unrelated to primer choice (e.g., sample handling).
- **Too few results** (< 5) give an under-determined GP; start with a
  2^k factorial over key additives to build a diverse training set.
- **Single-primer tests** cannot be used; active-learn works only on
  sets of the same target size you plan to keep using.
