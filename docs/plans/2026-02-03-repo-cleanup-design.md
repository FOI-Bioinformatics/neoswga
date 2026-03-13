**Note: This planning document is historical. The cleanup described here has been completed.**

# Repo Cleanup and Test Streamlining

## Goal

Reduce repo from 335 files to ~180 files by removing redundant tests, historical docs, and committed fixtures.

## Deletions

### Directories to Remove

| Path | Reason |
|------|--------|
| `validation_tests/` | Local scripts with hardcoded paths |
| `docs/archive/` | Historical QA reports |
| `docs/validation/` | Validation reports |
| `tests/integration/francisella_clinical/` | Factorial experiment data |
| `tests/integration/francisella_e2e/` | Redundant e2e tests |
| `tests/integration/equiphi29_full_factorial/` | Factorial experiment data |
| `tests/integration/equiphi29_long/` | Covered by parameterized test |
| `tests/integration/equiphi29_validation/` | Duplicate of baseline |
| `tests/integration/phi29_ga/` | Covered by parameterized test |
| `tests/integration/phi29_minimal/` | Duplicate of baseline |
| `tests/integration/phi29_uniformity/` | Duplicate of baseline |

### Files to Remove from Kept Directories

- `tests/integration/*/results/` - CSV fixtures (should be generated)
- `tests/integration/run_integration_tests.py` - Replaced by pytest
- `tests/integration/test_plan.md` - Outdated

## Consolidation

### New Integration Test Structure

```
tests/integration/
├── phi29_baseline/params.json
├── equiphi29_baseline/params.json
├── phi29_with_bg/params.json
└── test_integration.py        # Parameterized pytest
```

### test_integration.py Design

```python
import pytest
from pathlib import Path

SCENARIOS = [
    ("phi29_baseline", "phi29", 10),
    ("equiphi29_baseline", "equiphi29", 12),
    ("phi29_with_bg", "phi29", 10),
]

@pytest.mark.integration
@pytest.mark.parametrize("scenario,polymerase,primer_len", SCENARIOS)
def test_pipeline_scenario(scenario, polymerase, primer_len):
    """Run pipeline for each scenario, validate outputs."""
    # Skip if test genome not available
    # Load params.json
    # Run filter/score/optimize
    # Assert: primers returned, coverage > 0, no errors
```

## Verification

After cleanup:
1. `pytest tests/` passes
2. `git status` shows no untracked important files
3. File count reduced to ~180

## Result

| Metric | Before | After |
|--------|--------|-------|
| Total files | 335 | ~180 |
| Markdown docs | 50 | 24 |
| Integration test dirs | 12 | 3 |
| CSV fixtures | 39 | 0 |
