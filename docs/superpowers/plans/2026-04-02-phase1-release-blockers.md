# Phase 1: Release Blockers Remediation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the 10 release-blocking issues identified in the production readiness audit so NeoSWGA can be published to PyPI.

**Architecture:** Targeted fixes to existing modules — no new modules created. Each task is independent and can be committed separately. Tasks are ordered by risk: legal/security first, then correctness, then UX.

**Tech Stack:** Python 3.11+, pytest, pyproject.toml, argparse

---

### Task 1: Fix License Mismatch

The LICENSE file contains AGPL-3.0 but all metadata says MIT. The LICENSE file is the legally binding document, so update all metadata to match it.

**Files:**
- Modify: `pyproject.toml:11,34`
- Modify: `neoswga/__init__.py:16`

- [ ] **Step 1: Update pyproject.toml license field**

In `pyproject.toml`, change line 11:
```python
# Old:
license = {text = "MIT"}
# New:
license = {text = "AGPL-3.0-or-later"}
```

And change the classifier at line 34:
```python
# Old:
"License :: OSI Approved :: MIT License",
# New:
"License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)",
```

- [ ] **Step 2: Update __init__.py license**

In `neoswga/__init__.py`, change line 16:
```python
# Old:
__license__ = "MIT"
# New:
__license__ = "AGPL-3.0-or-later"
```

- [ ] **Step 3: Verify consistency**

Run:
```bash
grep -rn "MIT" pyproject.toml neoswga/__init__.py
```
Expected: No matches.

- [ ] **Step 4: Commit**

```bash
git add pyproject.toml neoswga/__init__.py
git commit -m "fix: reconcile license metadata to match AGPL-3.0 LICENSE file"
```

---

### Task 2: Remove Pickle Deserialization Fallbacks

Four sites silently fall back to unrestricted `pickle.load()` when `safe_load()` fails. This defeats the entire purpose of the safe_pickle module and is exploitable. Remove all fallbacks — if safe_load fails, let the error propagate.

**Files:**
- Modify: `neoswga/core/background_filter.py:289-295`
- Modify: `neoswga/core/background_filter.py:409-412`
- Modify: `neoswga/core/efficiency_predictor.py:389-395`
- Modify: `neoswga/core/rf_preprocessing.py:132-137`
- Modify: `neoswga/core/safe_pickle.py:112-124` (remove `unsafe_load_with_warning`)

- [ ] **Step 1: Write test for safe_load enforcement**

Create test at `tests/test_safe_pickle.py`:
```python
"""Tests for safe_pickle module — verifies no unsafe fallbacks exist."""
import pickle
import tempfile
import os
import pytest
from unittest.mock import patch


def test_safe_load_raises_on_disallowed_class():
    """safe_load must raise, not fall back, when encountering disallowed classes."""
    from neoswga.core.safe_pickle import safe_load

    # Create a pickle containing a disallowed class (os.system)
    payload = pickle.dumps({"action": os.getcwd})

    with tempfile.NamedTemporaryFile(suffix='.p', delete=False) as f:
        f.write(payload)
        tmp_path = f.name

    try:
        with pytest.raises(pickle.UnpicklingError):
            safe_load(tmp_path, context='sklearn_model')
    finally:
        os.unlink(tmp_path)


def test_no_unsafe_load_with_warning_function():
    """The unsafe_load_with_warning function should not exist."""
    from neoswga.core import safe_pickle
    assert not hasattr(safe_pickle, 'unsafe_load_with_warning'), \
        "unsafe_load_with_warning must be removed — it defeats safe_load"


def test_background_bloom_filter_load_no_fallback():
    """BackgroundBloomFilter.load must not fall back to raw pickle.load."""
    import inspect
    from neoswga.core.background_filter import BackgroundBloomFilter
    source = inspect.getsource(BackgroundBloomFilter.load)
    assert 'pickle.load' not in source, \
        "BackgroundBloomFilter.load must not contain pickle.load fallback"


def test_sampled_genome_index_load_no_fallback():
    """SampledGenomeIndex.load must not fall back to raw pickle.load."""
    import inspect
    from neoswga.core.background_filter import SampledGenomeIndex
    source = inspect.getsource(SampledGenomeIndex.load)
    assert 'pickle.load' not in source, \
        "SampledGenomeIndex.load must not contain pickle.load fallback"


def test_rf_preprocessing_load_no_fallback():
    """load_model must not fall back to raw pickle.load."""
    import inspect
    from neoswga.core.rf_preprocessing import load_model
    source = inspect.getsource(load_model)
    assert 'pickle.load' not in source and 'unsafe_load' not in source, \
        "load_model must not contain pickle.load or unsafe_load fallback"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_safe_pickle.py -v`
Expected: 3 of 5 tests FAIL (the source-inspection tests for background_filter, rf_preprocessing).

- [ ] **Step 3: Fix BackgroundBloomFilter.load in background_filter.py**

In `neoswga/core/background_filter.py`, replace lines 286-295:
```python
        try:
            from neoswga.core.safe_pickle import safe_load
            data = safe_load(path, context='bloom_filter')
        except Exception as e:
            logger.warning(
                f"Restricted unpickling failed ({e}), falling back to standard load. "
                "Only load files from trusted sources."
            )
            with open(path, 'rb') as f:
                data = pickle.load(f)
```

With:
```python
        from neoswga.core.safe_pickle import safe_load
        data = safe_load(path, context='bloom_filter')
```

- [ ] **Step 4: Fix SampledGenomeIndex.load in background_filter.py**

In `neoswga/core/background_filter.py`, replace lines 406-412:
```python
        try:
            from neoswga.core.safe_pickle import safe_load
            data = safe_load(path, context='bloom_filter')
        except Exception:
            logger.warning("Restricted unpickling failed, using standard load")
            with open(path, 'rb') as f:
                data = pickle.load(f)
```

With:
```python
        from neoswga.core.safe_pickle import safe_load
        data = safe_load(path, context='bloom_filter')
```

- [ ] **Step 5: Fix efficiency_predictor.py**

In `neoswga/core/efficiency_predictor.py`, replace lines 386-395:
```python
            try:
                from neoswga.core.safe_pickle import safe_load
                model = safe_load(model_path, context='sklearn_model')
            except Exception as e:
                logger.warning(
                    f"Restricted unpickling failed ({e}), falling back to standard load. "
                    "Only load files from trusted sources."
                )
                with open(model_path, 'rb') as f:
                    model = pickle.load(f)
```

With:
```python
            from neoswga.core.safe_pickle import safe_load
            model = safe_load(model_path, context='sklearn_model')
```

- [ ] **Step 6: Fix rf_preprocessing.py**

In `neoswga/core/rf_preprocessing.py`, replace lines 129-137:
```python
    try:
        from neoswga.core.safe_pickle import safe_load
        return safe_load(model_path, context='sklearn_model')
    except Exception as e:
        logging.getLogger(__name__).warning(
            f"Restricted unpickling failed ({e}), falling back to standard load"
        )
        from neoswga.core.safe_pickle import unsafe_load_with_warning
        return unsafe_load_with_warning(model_path, purpose='sklearn RF model')
```

With:
```python
    from neoswga.core.safe_pickle import safe_load
    return safe_load(model_path, context='sklearn_model')
```

- [ ] **Step 7: Remove unsafe_load_with_warning from safe_pickle.py**

In `neoswga/core/safe_pickle.py`, delete lines 112-124 (the entire `unsafe_load_with_warning` function).

- [ ] **Step 8: Run tests to verify they pass**

Run: `pytest tests/test_safe_pickle.py -v`
Expected: All 5 tests PASS.

- [ ] **Step 9: Run full test suite to check for regressions**

Run: `pytest tests/ -x --timeout=60`
Expected: No new failures.

- [ ] **Step 10: Commit**

```bash
git add neoswga/core/background_filter.py neoswga/core/efficiency_predictor.py neoswga/core/rf_preprocessing.py neoswga/core/safe_pickle.py tests/test_safe_pickle.py
git commit -m "security: remove all unsafe pickle.load() fallbacks

safe_load() must be the only deserialization path. Fallbacks
defeated the RestrictedUnpickler by catching all exceptions
and falling back to unrestricted pickle.load()."
```

---

### Task 3: Fix Incorrect Tm Calculation in adaptive_filters.py

The `ThermodynamicFilter._nn_tm()` method reimplements nearest-neighbor Tm calculation with an incorrect simplification: it drops the `R * ln(C/4)` primer concentration term, producing systematically wrong Tm values. Replace with the validated `thermodynamics.calculate_tm_with_salt()`.

**Files:**
- Modify: `neoswga/core/adaptive_filters.py:121-185`

- [ ] **Step 1: Write test for Tm consistency**

Create test at `tests/test_adaptive_filter_tm.py`:
```python
"""Test that adaptive_filters uses the canonical thermodynamics module for Tm."""
import pytest
import numpy as np


def test_thermodynamic_filter_tm_matches_canonical():
    """ThermodynamicFilter Tm must agree with thermodynamics.calculate_tm_with_salt."""
    from neoswga.core.adaptive_filters import ThermodynamicFilter
    from neoswga.core.thermodynamics import calculate_tm_with_salt

    tf = ThermodynamicFilter(min_tm=15.0, max_tm=60.0, na_conc=50.0)

    test_seqs = [
        'ATCGATCGATCG',
        'GCGCGCGCGC',
        'AAAATTTTCCCC',
        'ATATATATATAT',
        'GCTAGCTAGCTA',
    ]

    for seq in test_seqs:
        adaptive_tm = tf._calculate_tm(seq)
        canonical_tm = calculate_tm_with_salt(seq, na_conc=50.0)
        assert abs(adaptive_tm - canonical_tm) < 2.0, \
            f"Tm mismatch for {seq}: adaptive={adaptive_tm:.1f}, canonical={canonical_tm:.1f}"


def test_thermodynamic_filter_accepts_valid_primers():
    """Primers with Tm in range should pass the filter."""
    from neoswga.core.adaptive_filters import ThermodynamicFilter

    tf = ThermodynamicFilter(min_tm=15.0, max_tm=60.0, na_conc=50.0)
    # A typical 12-mer should have Tm in a reasonable range
    assert tf.passes_filter('ATCGATCGATCG')


def test_thermodynamic_filter_rejects_extreme_primers():
    """Very AT-rich short primers should have low Tm and potentially be rejected."""
    from neoswga.core.adaptive_filters import ThermodynamicFilter

    # Narrow Tm window that excludes low-Tm sequences
    tf = ThermodynamicFilter(min_tm=50.0, max_tm=70.0, na_conc=50.0)
    # Short AT-rich sequence has very low Tm
    result = tf.passes_filter('AAAAAA')
    assert result is False
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_adaptive_filter_tm.py -v`
Expected: `test_thermodynamic_filter_tm_matches_canonical` FAILS due to Tm divergence.

- [ ] **Step 3: Replace _nn_tm and _calculate_tm with canonical implementation**

In `neoswga/core/adaptive_filters.py`, replace the `_calculate_tm` and `_nn_tm` methods (lines ~121-185) with:

```python
    def _calculate_tm(self, seq: str) -> float:
        """
        Calculate melting temperature using the canonical thermodynamics module.

        Uses SantaLucia nearest-neighbor model with Owczarzy salt correction.
        """
        from neoswga.core.thermodynamics import calculate_tm_with_salt
        try:
            return calculate_tm_with_salt(seq, na_conc=self.na_conc)
        except (ValueError, ZeroDivisionError):
            # Fallback for very short or degenerate sequences
            return 0.0
```

Delete the `_nn_tm` method entirely.

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_adaptive_filter_tm.py -v`
Expected: All 3 tests PASS.

- [ ] **Step 5: Run full test suite for regressions**

Run: `pytest tests/ -x --timeout=60`
Expected: No new failures.

- [ ] **Step 6: Commit**

```bash
git add neoswga/core/adaptive_filters.py tests/test_adaptive_filter_tm.py
git commit -m "fix: replace incorrect Tm calculation in adaptive_filters

The _nn_tm() method omitted the R*ln(C/4) primer concentration
term and used different entropy values than the canonical
thermodynamics module, producing systematically wrong Tm values.
Now delegates to thermodynamics.calculate_tm_with_salt()."
```

---

### Task 4: Fix GA Optimizer ProcessPoolExecutor Crash

`genetic_algorithm.py` submits `self._evaluate_individual` (a bound method) to `ProcessPoolExecutor`. Bound methods of objects containing numpy arrays, caches, and complex objects fail to pickle on macOS (spawn start method). Switch to `ThreadPoolExecutor` since the GIL is released during numpy operations and file I/O.

**Files:**
- Modify: `neoswga/core/genetic_algorithm.py:341`

- [ ] **Step 1: Write test for GA population evaluation**

Create test at `tests/test_genetic_algorithm_eval.py`:
```python
"""Test that GA population evaluation does not crash due to pickling."""
import pytest
import numpy as np
from unittest.mock import MagicMock


def test_evaluate_population_does_not_use_process_pool():
    """GA must not use ProcessPoolExecutor with bound methods."""
    import inspect
    from neoswga.core.genetic_algorithm import GeneticAlgorithmOptimizer
    source = inspect.getsource(GeneticAlgorithmOptimizer._evaluate_population)
    assert 'ProcessPoolExecutor' not in source, \
        "ProcessPoolExecutor with bound methods fails on macOS spawn start method"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_genetic_algorithm_eval.py -v`
Expected: FAIL — source currently contains `ProcessPoolExecutor`.

- [ ] **Step 3: Replace ProcessPoolExecutor with ThreadPoolExecutor**

In `neoswga/core/genetic_algorithm.py`, find the import at the top of the file:
```python
from concurrent.futures import ProcessPoolExecutor, as_completed
```

Change to:
```python
from concurrent.futures import ThreadPoolExecutor, as_completed
```

Then in `_evaluate_population()` (around line 341), change:
```python
        with ProcessPoolExecutor(max_workers=self.config.n_processes) as executor:
```

To:
```python
        with ThreadPoolExecutor(max_workers=self.config.n_processes) as executor:
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_genetic_algorithm_eval.py tests/test_genetic_algorithm_integration.py -v`
Expected: All PASS.

- [ ] **Step 5: Commit**

```bash
git add neoswga/core/genetic_algorithm.py tests/test_genetic_algorithm_eval.py
git commit -m "fix: replace ProcessPoolExecutor with ThreadPoolExecutor in GA

Bound methods of complex objects (with numpy arrays, caches)
fail to pickle on macOS spawn start method, crashing the GA
optimizer. ThreadPoolExecutor avoids pickling entirely."
```

---

### Task 5: Fix Placeholder Emails

Replace `example.com` placeholder emails with real contact info.

**Files:**
- Modify: `pyproject.toml:13,16`
- Modify: `neoswga/__init__.py:15`

- [ ] **Step 1: Update pyproject.toml**

In `pyproject.toml`, replace `andreas.sjodin@example.com` on lines 13 and 16 with the actual email address. Ask the user for their preferred contact email if unknown.

**NOTE TO IMPLEMENTER**: Ask the user what email to use before making this change. Do not guess.

- [ ] **Step 2: Update __init__.py**

In `neoswga/__init__.py`, update line 15 to match the email set in pyproject.toml.

- [ ] **Step 3: Update GitHub URL in __init__.py**

In `neoswga/__init__.py`, line 31, change:
```python
    print("Documentation: https://github.com/andreassjodin/neoswga")
```
To:
```python
    print("Documentation: https://github.com/FOI-Bioinformatics/neoswga")
```

- [ ] **Step 4: Commit**

```bash
git add pyproject.toml neoswga/__init__.py
git commit -m "fix: replace placeholder email and fix GitHub URL"
```

---

### Task 6: Fix MANIFEST.in

Remove references to non-existent files. Remove CLAUDE.md from distribution. Add LICENSE (not LICENSE.txt).

**Files:**
- Modify: `MANIFEST.in`

- [ ] **Step 1: Rewrite MANIFEST.in**

Replace the entire contents of `MANIFEST.in` with:
```
# Include documentation
include README.md
include LICENSE

# Include configuration files
include pyproject.toml

# Include Python source and model files
recursive-include neoswga *.py
recursive-include neoswga/core/models *.p

# Exclude compiled files
global-exclude *.pyc
global-exclude *.pyo
global-exclude __pycache__
global-exclude .DS_Store

# Exclude test files and caches
prune tests
prune validation_tests
prune examples
prune docs
prune scripts
prune .pytest_cache
prune .mypy_cache
prune htmlcov
prune dist
prune build
prune *.egg-info

# Exclude version control and IDE files
prune .git
prune .vscode
prune .idea
exclude .gitignore
exclude .gitattributes
exclude CLAUDE.md
```

- [ ] **Step 2: Verify with a test build**

Run:
```bash
python -m build --sdist 2>&1 | tail -5
```
Expected: No warnings about missing files.

- [ ] **Step 3: Commit**

```bash
git add MANIFEST.in
git commit -m "fix: clean up MANIFEST.in — remove missing file references, prune tests/docs"
```

---

### Task 7: Verify Model File Inclusion in Wheel

The pickle model at `neoswga/core/models/random_forest_filter.p` must be in the wheel. The `models/` directory has no `__init__.py` but `package-data` should still include it.

**Files:**
- No files modified if test passes; otherwise create `neoswga/core/models/__init__.py`

- [ ] **Step 1: Build and check wheel contents**

Run:
```bash
python -m build --wheel 2>&1 | tail -3 && unzip -l dist/neoswga-*.whl | grep -E '(random_forest|models)'
```
Expected: `neoswga/core/models/random_forest_filter.p` appears in the listing.

- [ ] **Step 2: If model is missing, add __init__.py**

Only if the model file is NOT in the wheel:
```bash
touch neoswga/core/models/__init__.py
```
Then rebuild and verify.

- [ ] **Step 3: Commit (only if __init__.py was needed)**

```bash
git add neoswga/core/models/__init__.py
git commit -m "fix: add __init__.py to models/ to ensure pickle model is included in wheel"
```

---

### Task 8: Wire ParamValidator Into Pipeline Commands

Currently `ParamValidator` only runs when the user explicitly runs `neoswga validate-params`. It should run automatically (at ERROR level only) before pipeline commands.

**Files:**
- Modify: `neoswga/cli_unified.py` (in `validate_params_json_file` function, around line 301)

- [ ] **Step 1: Write test for automatic validation**

Add to `tests/test_cli_smoke.py` or create `tests/test_param_validation_integration.py`:
```python
"""Test that pipeline commands auto-validate params.json."""
import json
import tempfile
import os
import pytest


def test_invalid_param_type_caught_before_pipeline():
    """A string value for min_k should produce a clear error, not a traceback."""
    from neoswga.cli_unified import validate_params_json_file

    params = {
        'data_dir': '/tmp/test_neoswga',
        'fg_genome': '/tmp/test.fna',
        'fg_prefixes': ['/tmp/test'],
        'min_k': 'six',  # Invalid type
    }

    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        json.dump(params, f)
        tmp_path = f.name

    try:
        # Should not crash — validate_params_json_file handles this gracefully
        # (it calls sys.exit on errors, so we need to catch SystemExit)
        with pytest.raises(SystemExit):
            validate_params_json_file(tmp_path)
    finally:
        os.unlink(tmp_path)
```

- [ ] **Step 2: Add ParamValidator call to validate_params_json_file**

In `neoswga/cli_unified.py`, at the end of `validate_params_json_file()` (after line 349), add:
```python
    # Run automatic parameter validation (ERROR-level issues only)
    try:
        from neoswga.core.param_validator import ParamValidator
        validator = ParamValidator(data)
        results = validator.validate()
        errors = [r for r in results if r.level == 'error']
        if errors:
            print(f"\nParameter validation found {len(errors)} error(s) in '{path}':",
                  file=sys.stderr)
            for err in errors:
                print(f"  - {err.message}", file=sys.stderr)
            print(f"\nRun 'neoswga validate-params -j {path}' for full details.",
                  file=sys.stderr)
            sys.exit(1)
    except ImportError:
        pass  # ParamValidator not available — skip validation
    except Exception as e:
        # Don't block pipeline on validator bugs
        logger.debug(f"Parameter validation skipped: {e}")
```

- [ ] **Step 3: Run the test**

Run: `pytest tests/test_param_validation_integration.py -v`
Expected: PASS.

- [ ] **Step 4: Run full test suite**

Run: `pytest tests/ -x --timeout=60`
Expected: No new failures.

- [ ] **Step 5: Commit**

```bash
git add neoswga/cli_unified.py tests/test_param_validation_integration.py
git commit -m "feat: auto-validate params.json before pipeline commands

Runs ParamValidator at ERROR level during validate_params_json_file().
Users get clear messages for type errors and invalid values instead
of raw Python tracebacks downstream."
```

---

### Task 9: Add Empty-Result Detection at End of Step 2

When all primers are filtered out, step2 writes an empty CSV. The user doesn't find out until step3 starts. Detect and report immediately.

**Files:**
- Modify: `neoswga/core/pipeline.py` (in `step2()`, after the DataFrame is written)

- [ ] **Step 1: Locate the step2 CSV write**

Read `neoswga/core/pipeline.py` and find where `step2_df.csv` is written (the `to_csv` call after filtering). The fix goes immediately after that write.

- [ ] **Step 2: Add empty-result check**

After the `filtered_rate_df.to_csv(...)` call in `step2()`, add:
```python
    if len(filtered_rate_df) == 0:
        logger.error(
            "No primers passed filtering. The output file is empty.\n"
            "Suggestions to get more candidates:\n"
            "  - Increase max_bg_freq (e.g., 1e-5)\n"
            "  - Decrease min_fg_freq (e.g., 1e-6)\n"
            "  - Increase max_gini (e.g., 0.8)\n"
            "  - Widen the k-mer range (min_k/max_k)\n"
            "  - Check that foreground and background genomes are not swapped"
        )
        raise ValueError(
            f"Filtering produced 0 candidate primers. "
            f"Started with {len(rate_df)} k-mers, all were removed by filters. "
            f"Adjust filtering thresholds in params.json."
        )
```

- [ ] **Step 3: Run existing tests**

Run: `pytest tests/test_pipeline.py -v`
Expected: No regressions.

- [ ] **Step 4: Commit**

```bash
git add neoswga/core/pipeline.py
git commit -m "fix: detect and report empty filter results at end of step2

Instead of silently writing an empty CSV that causes a confusing
error at step3, raise immediately with actionable suggestions
for adjusting filter thresholds."
```

---

### Task 10: Fix Division by Zero in filter.py

`get_bg_rates_via_bloom()` divides by `total` (length of `primer_list`) at line 170 without checking for zero.

**Files:**
- Modify: `neoswga/core/filter.py:169-170`

- [ ] **Step 1: Write test for empty primer list**

Add to `tests/test_filter.py`:
```python
def test_get_bg_rates_via_bloom_empty_list():
    """get_bg_rates_via_bloom should handle an empty primer list without crashing."""
    from neoswga.core.filter import get_bg_rates_via_bloom
    from unittest.mock import MagicMock

    mock_bloom = MagicMock()
    result = get_bg_rates_via_bloom([], mock_bloom)
    assert result == {}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_filter.py::test_get_bg_rates_via_bloom_empty_list -v`
Expected: FAIL with `ZeroDivisionError`.

- [ ] **Step 3: Add guard for empty list**

In `neoswga/core/filter.py`, before line 169 (the `logger.info` with division), add:
```python
    if total == 0:
        logger.info("Bloom filter: no primers to check")
        return primer_to_count
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_filter.py::test_get_bg_rates_via_bloom_empty_list -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add neoswga/core/filter.py tests/test_filter.py
git commit -m "fix: guard against division by zero in get_bg_rates_via_bloom

Empty primer list caused ZeroDivisionError at the logging
statement. Now returns early with empty dict."
```

---

## Summary

| Task | Issue | Severity | Est. Time |
|------|-------|----------|-----------|
| 1 | License mismatch | Critical | 10 min |
| 2 | Pickle fallbacks | Critical | 30 min |
| 3 | Incorrect Tm calculation | Critical | 20 min |
| 4 | GA ProcessPoolExecutor crash | Critical | 15 min |
| 5 | Placeholder emails | Critical | 5 min |
| 6 | MANIFEST.in broken | High | 10 min |
| 7 | Model file in wheel | High | 10 min |
| 8 | ParamValidator wiring | High | 20 min |
| 9 | Empty filter results | High | 15 min |
| 10 | Division by zero | High | 10 min |

**Total estimated: ~2.5 hours**

All tasks are independent and can be executed in any order. Each produces a clean, atomic commit.
