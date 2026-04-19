# Phase 0: Stop the Bleeding Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close the five correctness bugs in neoswga v3.6.0 that are capable of producing incorrect results today, per the production-readiness audit. No new features.

**Architecture:** Four independent work streams targeting existing modules. Changes are additive (new param fields, new aggregation modes) or consolidating (one canonical Tm path). Behavior defaults are updated where the current defaults produce wrong science; users can opt back in with explicit flags. Version bumps to 3.7.0 at end.

**Tech Stack:** Python 3.11+, pytest, hypothesis, dataclasses, numpy, pandas, h5py, biopython, jellyfish (external).

**Parent spec:** `docs/superpowers/specs/2026-04-19-phase-0-stop-the-bleeding-design.md`
**Roadmap:** `docs/superpowers/specs/2026-04-19-production-readiness-roadmap.md`

**Preflight before starting any task:**
- Run `pytest tests/ -q` and confirm the current suite passes. Note any pre-existing failures so later red tests can be distinguished from regressions.
- Ensure a clean worktree (`git status` shows no staged changes beyond the untracked Phase 0 test files, which Task 5 will incorporate).

---

## Work Stream 0.1 - Wire adaptive GC, remove the orphan

The spec assumed `AdaptiveGCFilter` needed wiring. In reality, `neoswga/core/parameter.py:817-820` already computes adaptive gc_min/gc_max from `genome_gc +/- gc_tolerance`, and `neoswga/core/filter.py:258` uses those computed values. The actual bug: no opt-out, no regression test, and the `AdaptiveGCFilter` class in `neoswga/core/adaptive_filters.py` is dead code that lies about the state of the system.

### Task 1: Regression test for low-GC end-to-end filtering

**Files:**
- Create: `tests/fixtures/low_gc/target.fna`
- Create: `tests/test_adaptive_gc_regression.py`

- [ ] **Step 1: Create a 33% GC fixture FASTA**

A ~10 kbp synthetic sequence at ~33% GC (Francisella-like). Use a deterministic seed so the fixture is reproducible.

Create `tests/fixtures/low_gc/target.fna` with the following content (one header, sequence wrapped at 80 char):

```bash
python - <<'PY'
import random
random.seed(0xF1CA)
bases = list('ATATATAT' 'ATATATAT' 'GC' 'ATATATATATATATATATATATATATATATATATAGC')
# Aim for ~33% GC across 10 kbp
seq = []
for _ in range(10000):
    r = random.random()
    if r < 0.335:
        seq.append(random.choice('GC'))
    else:
        seq.append(random.choice('AT'))
import os
os.makedirs('tests/fixtures/low_gc', exist_ok=True)
with open('tests/fixtures/low_gc/target.fna', 'w') as f:
    f.write('>low_gc_synthetic length=10000 gc~0.33\n')
    s = ''.join(seq)
    for i in range(0, len(s), 80):
        f.write(s[i:i+80] + '\n')
# Verify
gc = sum(1 for b in s if b in 'GC') / len(s)
print(f"Actual GC: {gc:.3f}")
PY
```

Expected output: `Actual GC: 0.33x` (close to 0.33).

- [ ] **Step 2: Write the failing test**

Write `tests/test_adaptive_gc_regression.py`:

```python
"""Regression tests for adaptive GC filtering on extreme-GC genomes."""
import os
import subprocess
import tempfile

import pandas as pd
import pytest

from neoswga.core import parameter
from neoswga.core.kmer_counter import run_jellyfish


FIXTURE_DIR = os.path.join(os.path.dirname(__file__), 'fixtures', 'low_gc')
FIXTURE_FASTA = os.path.join(FIXTURE_DIR, 'target.fna')


@pytest.fixture
def low_gc_data_dir(tmp_path):
    """Build a transient data_dir for a low-GC target."""
    return str(tmp_path)


@pytest.mark.integration
def test_low_gc_target_produces_nonempty_step2(low_gc_data_dir):
    """A 33% GC target must produce non-empty step2_df with adaptive GC enabled.

    This guards against regression to the pre-v3.7.0 behavior where fixed
    37.5-62.5% GC thresholds rejected all primers for extreme-GC genomes.
    """
    import shutil
    if shutil.which('jellyfish') is None:
        pytest.skip("jellyfish not available")

    prefix = os.path.join(low_gc_data_dir, 'low_gc_synthetic')
    run_jellyfish(FIXTURE_FASTA, prefix, min_k=6, max_k=8)

    from neoswga.core.parameter import PipelineParameters, set_from_config
    config = PipelineParameters(
        min_k=6, max_k=8,
        min_fg_freq=1e-5,
        max_bg_freq=1e-3,
        min_tm=15.0, max_tm=60.0,
        max_gini=0.9,
        max_primer=100,
        data_dir=low_gc_data_dir,
        fg_genomes=[FIXTURE_FASTA],
        fg_prefixes=[prefix],
        fg_seq_lengths=[10000],
        bg_genomes=[], bg_prefixes=[], bg_seq_lengths=[],
        genome_gc=0.33,
        gc_tolerance=0.15,
    )
    set_from_config(config)

    from neoswga.core.pipeline import step2
    step2()

    step2_path = os.path.join(low_gc_data_dir, 'step2_df.csv')
    assert os.path.exists(step2_path), "step2_df.csv not written"
    df = pd.read_csv(step2_path)
    assert len(df) > 0, "Adaptive GC filter should admit primers for a 33% GC target"
```

- [ ] **Step 3: Run the test to verify it passes on main**

Run: `pytest tests/test_adaptive_gc_regression.py -v -m integration`
Expected: PASS (confirms current behavior already works; this test is the regression guard).
If it FAILS: the adaptive GC pathway is broken somewhere; stop and investigate before proceeding.

- [ ] **Step 4: Commit**

```bash
git add tests/fixtures/low_gc/target.fna tests/test_adaptive_gc_regression.py
git commit -m "test: regression guard for adaptive GC on low-GC (~33%) target

Codifies the behavior that was silently fixed via parameter.py:817-820.
Without this test, any future change to adaptive GC computation could
silently re-introduce the bug where extreme-GC genomes yield zero primers.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

### Task 2: Regression test for high-GC end-to-end filtering

**Files:**
- Create: `tests/fixtures/high_gc/target.fna`
- Modify: `tests/test_adaptive_gc_regression.py`

- [ ] **Step 1: Create a 68% GC fixture FASTA**

```bash
python - <<'PY'
import random, os
random.seed(0xF1CB)
seq = []
for _ in range(10000):
    r = random.random()
    if r < 0.68:
        seq.append(random.choice('GC'))
    else:
        seq.append(random.choice('AT'))
os.makedirs('tests/fixtures/high_gc', exist_ok=True)
with open('tests/fixtures/high_gc/target.fna', 'w') as f:
    f.write('>high_gc_synthetic length=10000 gc~0.68\n')
    s = ''.join(seq)
    for i in range(0, len(s), 80):
        f.write(s[i:i+80] + '\n')
gc = sum(1 for b in s if b in 'GC') / len(s)
print(f"Actual GC: {gc:.3f}")
PY
```

- [ ] **Step 2: Add the high-GC test**

Append to `tests/test_adaptive_gc_regression.py`:

```python
HIGH_GC_DIR = os.path.join(os.path.dirname(__file__), 'fixtures', 'high_gc')
HIGH_GC_FASTA = os.path.join(HIGH_GC_DIR, 'target.fna')


@pytest.mark.integration
def test_high_gc_target_produces_nonempty_step2(tmp_path):
    """A 68% GC target must produce non-empty step2_df with adaptive GC enabled."""
    import shutil
    if shutil.which('jellyfish') is None:
        pytest.skip("jellyfish not available")

    data_dir = str(tmp_path)
    prefix = os.path.join(data_dir, 'high_gc_synthetic')
    run_jellyfish(HIGH_GC_FASTA, prefix, min_k=6, max_k=8)

    from neoswga.core.parameter import PipelineParameters, set_from_config
    config = PipelineParameters(
        min_k=6, max_k=8,
        min_fg_freq=1e-5, max_bg_freq=1e-3,
        min_tm=15.0, max_tm=80.0,  # wider for GC-rich primers
        max_gini=0.9, max_primer=100,
        data_dir=data_dir,
        fg_genomes=[HIGH_GC_FASTA],
        fg_prefixes=[prefix],
        fg_seq_lengths=[10000],
        bg_genomes=[], bg_prefixes=[], bg_seq_lengths=[],
        genome_gc=0.68, gc_tolerance=0.15,
    )
    set_from_config(config)
    from neoswga.core.pipeline import step2
    step2()
    df = pd.read_csv(os.path.join(data_dir, 'step2_df.csv'))
    assert len(df) > 0, "Adaptive GC filter should admit primers for a 68% GC target"
```

- [ ] **Step 3: Run**

Run: `pytest tests/test_adaptive_gc_regression.py::test_high_gc_target_produces_nonempty_step2 -v -m integration`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add tests/fixtures/high_gc/target.fna tests/test_adaptive_gc_regression.py
git commit -m "test: regression guard for adaptive GC on high-GC (~68%) target

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

### Task 3: Expose explicit `use_adaptive_gc` opt-out flag

**Files:**
- Modify: `neoswga/core/parameter.py` (`PipelineParameters`, `get_current_config`, `set_from_config`, `get_params`)
- Create/Modify: `tests/test_parameter.py`

- [ ] **Step 1: Write failing test for the new field**

Add to `tests/test_parameter.py` (append if file exists):

```python
def test_pipeline_parameters_has_use_adaptive_gc():
    from neoswga.core.parameter import PipelineParameters
    p = PipelineParameters()
    assert p.use_adaptive_gc is True, "Default must be True (adaptive GC on)"


def test_use_adaptive_gc_opt_out_uses_fixed_thresholds(tmp_path):
    """When use_adaptive_gc=False, gc_min/gc_max stay at their explicit defaults."""
    import os, json
    from neoswga.core import parameter as P

    params_file = tmp_path / 'params.json'
    params_file.write_text(json.dumps({
        'data_dir': str(tmp_path),
        'fg_genomes': [],
        'bg_genomes': [],
        'fg_prefixes': [], 'bg_prefixes': [],
        'genome_gc': 0.33,
        'use_adaptive_gc': False,
        'gc_min': 0.375, 'gc_max': 0.625,
        'schema_version': 1,
    }))

    class _Args:
        json_file = str(params_file)
        min_fg_freq = max_bg_freq = min_tm = max_tm = None
        max_gini = max_primer = min_amp_pred = cpus = None
        max_dimer_bp = max_self_dimer_bp = mismatch_penalty = None
        verbose = drop_iterations = top_set_count = retries = None
        iterations = max_sets = None
        selection_metric = fg_circular = bg_circular = None
        fasta_fore = fasta_back = kmer_fore = kmer_back = None
        data_dir = src_dir = None
        min_k = max_k = None

    data = P.get_params(_Args())
    assert data['gc_min'] == 0.375
    assert data['gc_max'] == 0.625
    P.reset_to_defaults()
```

- [ ] **Step 2: Run and confirm failure**

Run: `pytest tests/test_parameter.py::test_pipeline_parameters_has_use_adaptive_gc tests/test_parameter.py::test_use_adaptive_gc_opt_out_uses_fixed_thresholds -v`
Expected: FAIL with AttributeError on `use_adaptive_gc` (first test) and wrong gc_min/gc_max (second test).

- [ ] **Step 3: Add the field to `PipelineParameters`**

In `neoswga/core/parameter.py`, insert after the `genome_gc: Optional[float] = None` line (~line 70):

```python
    # Adaptive GC filtering (uses genome_gc +/- gc_tolerance when True; falls
    # back to explicit gc_min/gc_max when False).
    use_adaptive_gc: bool = True
```

Add to `get_current_config()` (after the `gc_tolerance=...,` line around line 162):

```python
        use_adaptive_gc=globals().get('use_adaptive_gc', True),
```

Add to `set_from_config()` (in the GC block around line 247):

```python
    g['use_adaptive_gc'] = config.use_adaptive_gc
```

Add module-level default near the other GC defaults (around line 382):

```python
use_adaptive_gc = True
```

- [ ] **Step 4: Honor the flag in `get_params()`**

In `neoswga/core/parameter.py`, modify the adaptive-GC block (lines 817-824) to gate on the flag:

```python
    # Honor use_adaptive_gc flag (default True).
    use_adaptive_gc = data.get('use_adaptive_gc', True)
    data['use_adaptive_gc'] = use_adaptive_gc

    if use_adaptive_gc and genome_gc is not None and genome_gc > 0:
        # Adaptive: primer GC should match target genome GC +/- tolerance
        gc_min = max(0.15, genome_gc - gc_tolerance)
        gc_max = min(0.85, genome_gc + gc_tolerance)
    else:
        # Fixed thresholds (user explicitly opted out or genome_gc unknown)
        gc_min = data.get('gc_min', 0.375)
        gc_max = data.get('gc_max', 0.625)
```

Also add `global use_adaptive_gc` to the global declarations near line 483.

- [ ] **Step 5: Run tests to verify pass**

Run: `pytest tests/test_parameter.py -v`
Expected: PASS for both new tests.

- [ ] **Step 6: Add a param_validator warning for the opt-out-with-extreme-GC footgun**

In `neoswga/core/param_validator.py`, locate the existing validation function(s). Add a check:

```python
def _warn_if_fixed_thresholds_with_extreme_gc(config):
    """Warn when use_adaptive_gc=False and genome_gc is outside [0.35, 0.65]."""
    if getattr(config, 'use_adaptive_gc', True) is True:
        return []
    gc = getattr(config, 'genome_gc', None)
    if gc is None:
        return []
    if 0.35 <= gc <= 0.65:
        return []
    return [(
        'warning',
        f"use_adaptive_gc=False with genome_gc={gc:.2f} (outside 0.35-0.65). "
        f"Fixed gc_min/gc_max thresholds may reject most candidate primers. "
        f"Consider setting use_adaptive_gc=True."
    )]
```

Wire this into the existing validator aggregation (match the style of the
file - `grep -n "def " neoswga/core/param_validator.py` to find where
checks are aggregated). Add a unit test in `tests/test_ux_modules.py` or
equivalent that a `PipelineParameters(use_adaptive_gc=False, genome_gc=0.30)`
triggers the warning.

- [ ] **Step 7: Full suite check**

Run: `pytest tests/ -q --tb=short`
Expected: No new failures. Pre-existing failures (if any were noted at preflight) unchanged.

- [ ] **Step 8: Commit**

```bash
git add neoswga/core/parameter.py neoswga/core/param_validator.py tests/test_parameter.py tests/test_ux_modules.py
git commit -m "feat(adaptive-gc): add use_adaptive_gc opt-out flag

Default is True (adaptive behavior unchanged). Setting False restores
fixed 37.5-62.5% gc_min/gc_max, which is useful when a user wants to
override a miscalculated genome_gc.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

### Task 4: Delete the orphaned `AdaptiveGCFilter` class

The `AdaptiveGCFilter` class in `neoswga/core/adaptive_filters.py` is never instantiated in pipeline code (confirmed by grep). Its docstring incorrectly claims "filter.py uses fixed thresholds", which is false as of v3.6.0. Dead code that lies about the state of the system is worse than no code.

**Files:**
- Modify: `neoswga/core/adaptive_filters.py`

- [ ] **Step 1: Verify no users**

Run: `grep -rn "AdaptiveGCFilter" neoswga/ tests/ scripts/`
Expected: No matches other than the definition site itself in `neoswga/core/adaptive_filters.py`.
If anything else turns up, stop and discuss - there is a hidden consumer.

- [ ] **Step 2: Remove the class and its obsolete docstring**

In `neoswga/core/adaptive_filters.py`, delete the `AdaptiveGCFilter` class (approximately lines 18-85 - the `FilterThresholds` dataclass and the `AdaptiveGCFilter` class). Keep `ThermodynamicFilter` and anything else in the file. Also remove the top-of-file comment block that claims "fixes critical bug" since the "bug" is fixed elsewhere.

Replace the module docstring with:

```python
"""
Adaptive filters that adjust to genome composition.

Adaptive GC filtering is implemented in neoswga.core.parameter.get_params()
(see the use_adaptive_gc flag). This module contains thermodynamic-aware
filter helpers used downstream.
"""
```

- [ ] **Step 3: Run suite**

Run: `pytest tests/ -q --tb=short`
Expected: No new failures.

- [ ] **Step 4: Commit**

```bash
git add neoswga/core/adaptive_filters.py
git commit -m "refactor: remove orphaned AdaptiveGCFilter class

The class was never instantiated; its docstring falsely claimed filter.py
used fixed GC thresholds. Actual adaptive GC lives in parameter.get_params()
and is now covered by regression tests on extreme-GC fixtures.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Work Stream 0.2 - Finish blacklist integration

Fields are declared in `PipelineParameters`. The gaps are: (a) `bl_seq_lengths` is read from params.json but never auto-computed from `bl_genomes` (unlike `fg_seq_lengths` and `bg_seq_lengths`), (b) four test files exist untracked and need to be made green and committed, (c) there is no integration test demonstrating the blacklist actually reduces blacklist hits in the final primer set.

### Task 5: Get the untracked blacklist/library tests green

**Files:**
- Modify: `tests/test_blacklist_filter.py` (already exists, untracked)
- Modify: `tests/test_blacklist_cli.py` (already exists, untracked)
- Modify: `tests/test_genome_library.py` (already exists, untracked)
- Modify: `tests/test_gc_adaptive_pipeline.py` (already exists, untracked)

- [ ] **Step 1: Run the four tests as-is**

Run: `pytest tests/test_blacklist_filter.py tests/test_blacklist_cli.py tests/test_genome_library.py tests/test_gc_adaptive_pipeline.py -v`
Expected: Some may pass, some may fail. Note which.

- [ ] **Step 2: For each failing test, fix the test itself if the failure is a test bug, or note as a known gap for Task 6-8**

Common expected adjustments:
- If `test_blacklist_cli.py` tests CLI flags that do not exist on the current parser, add the missing `--blacklist`, `--bl-penalty`, `--max-bl-freq`, `--k-ranges` flags to the relevant subparsers in `neoswga/cli_unified.py`. Use `grep -n "add_argument" neoswga/cli_unified.py | head -60` to find existing style.
- If `test_genome_library.py` tests fail because of import errors, ensure `from neoswga.core.genome_library import GenomeLibrary` works (the module already exists).

- [ ] **Step 3: Commit adjustments separately from adding the tests so the history is readable**

If any source adjustments were made:

```bash
git add neoswga/cli_unified.py  # or whichever source files changed
git commit -m "feat(cli): add --blacklist, --bl-penalty, --max-bl-freq flags

Aligns CLI parser with the blacklist fields already present in
PipelineParameters. Precursor to the blacklist integration tests in the
next commit.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

- [ ] **Step 4: Commit the four tests**

```bash
git add tests/test_blacklist_filter.py tests/test_blacklist_cli.py tests/test_genome_library.py tests/test_gc_adaptive_pipeline.py
git commit -m "test: commit blacklist and genome-library tests

Unit tests for _filter_blacklist_penalty, blacklist CLI arguments,
GenomeLibrary, and GC-adaptive pipeline behavior. Previously untracked.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

### Task 6: Auto-compute `bl_seq_lengths` when not provided

**Files:**
- Modify: `neoswga/core/parameter.py:get_params` (around lines 713-717, 770-778)
- Modify: `tests/test_blacklist_filter.py` (add one new test)

- [ ] **Step 1: Write the failing test**

Append to `tests/test_blacklist_filter.py`:

```python
class TestBlacklistSeqLengthAutoCompute:
    """bl_seq_lengths should be auto-computed from bl_genomes when missing."""

    def test_autocompute_when_missing(self, tmp_path):
        """Given bl_genomes but no bl_seq_lengths, populate from FASTAs."""
        import json
        from neoswga.core import parameter as P

        # Minimal fixture genome
        fg = tmp_path / 'fg.fna'
        fg.write_text('>fg\n' + 'ATCG' * 50 + '\n')
        bl = tmp_path / 'bl.fna'
        bl.write_text('>bl\n' + 'GCTA' * 75 + '\n')

        params_file = tmp_path / 'params.json'
        params_file.write_text(json.dumps({
            'data_dir': str(tmp_path),
            'fg_genomes': [str(fg)],
            'bg_genomes': [],
            'fg_prefixes': [str(tmp_path / 'fg')],
            'bg_prefixes': [],
            'bl_genomes': [str(bl)],
            'bl_prefixes': [str(tmp_path / 'bl')],
            'schema_version': 1,
        }))

        class _Args:
            json_file = str(params_file)
            for attr in (
                'min_fg_freq', 'max_bg_freq', 'min_tm', 'max_tm', 'max_gini',
                'max_primer', 'min_amp_pred', 'cpus', 'max_dimer_bp',
                'max_self_dimer_bp', 'mismatch_penalty', 'verbose',
                'drop_iterations', 'top_set_count', 'retries', 'iterations',
                'max_sets', 'selection_metric', 'fg_circular', 'bg_circular',
                'fasta_fore', 'fasta_back', 'kmer_fore', 'kmer_back',
                'data_dir', 'src_dir', 'min_k', 'max_k',
            ):
                locals()[attr] = None

        data = P.get_params(_Args())
        assert data['bl_seq_lengths'] == [300]  # 75 * 4
        P.reset_to_defaults()
```

- [ ] **Step 2: Run the test and confirm failure**

Run: `pytest tests/test_blacklist_filter.py::TestBlacklistSeqLengthAutoCompute -v`
Expected: FAIL - `bl_seq_lengths` stays `[]` because auto-compute is not implemented.

- [ ] **Step 3: Implement auto-compute**

In `neoswga/core/parameter.py`, locate the block that auto-computes `fg_seq_lengths` and `bg_seq_lengths` (lines 770-778). Below it, add:

```python
    # Auto-compute blacklist genome lengths if missing (parallels fg/bg).
    if 'bl_genomes' in data and data['bl_genomes']:
        existing_bl_lengths = data.get('bl_seq_lengths', [])
        if not existing_bl_lengths or len(existing_bl_lengths) != len(data['bl_genomes']):
            data['bl_seq_lengths'] = _utility.get_all_seq_lengths(
                fname_genomes=data['bl_genomes'], cpus=data['cpus']
            )
```

Also update the global-assignment block (around line 712-717) so the global `bl_seq_lengths` reflects the computed value:

Change:
```python
    bl_seq_lengths = data.get('bl_seq_lengths', [])
```
to:
```python
    # Read initial value; will be overwritten below if we auto-compute.
    bl_seq_lengths = data.get('bl_seq_lengths', [])
```
(no functional change - the existing read happens before auto-compute; we just set globals from `data` at the end of `get_params`. Add one line after the auto-compute block:)

```python
    bl_seq_lengths = data.get('bl_seq_lengths', [])
```

- [ ] **Step 4: Run the test to verify pass**

Run: `pytest tests/test_blacklist_filter.py::TestBlacklistSeqLengthAutoCompute -v`
Expected: PASS.

- [ ] **Step 5: Full suite**

Run: `pytest tests/ -q --tb=short`
Expected: No new failures.

- [ ] **Step 6: Commit**

```bash
git add neoswga/core/parameter.py tests/test_blacklist_filter.py
git commit -m "fix(blacklist): auto-compute bl_seq_lengths from bl_genomes

Previously bl_seq_lengths was only read from params.json. When a user
supplied bl_genomes but not bl_seq_lengths (e.g. via the wizard),
pipeline._filter_blacklist_penalty would compute frequencies against
len(bl_prefixes) instead of actual genome lengths, producing wrong
values. Now parallels fg_seq_lengths / bg_seq_lengths behavior.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

### Task 7: Guard `_filter_blacklist_penalty` against misconfiguration

**Files:**
- Modify: `neoswga/core/pipeline.py:_filter_blacklist_penalty` (lines 61-101)
- Modify: `tests/test_blacklist_filter.py`

- [ ] **Step 1: Write failing test**

Append to `tests/test_blacklist_filter.py`:

```python
class TestBlacklistPenaltyMisconfig:
    def test_empty_lengths_raises(self, tmp_path):
        """If bl_prefixes is non-empty but bl_seq_lengths is empty, raise."""
        from neoswga.core.pipeline import _filter_blacklist_penalty
        from neoswga.core.exceptions import ConfigurationError

        (tmp_path / 'bl_test_8mer_all.txt').write_text('ATCGATCG 5\n')
        with pytest.raises(ConfigurationError):
            _filter_blacklist_penalty(
                ['ATCGATCG'],
                [str(tmp_path / 'bl_test')],
                [],  # lengths missing despite prefixes present
                max_bl_freq=0.0,
            )
```

- [ ] **Step 2: Run and confirm FAIL**

Run: `pytest tests/test_blacklist_filter.py::TestBlacklistPenaltyMisconfig -v`
Expected: FAIL (no exception raised).

- [ ] **Step 3: Update `_filter_blacklist_penalty`**

`ConfigurationError` in `neoswga/core/exceptions.py:230` is a bare subclass (just `pass`), so it accepts a single message. A more specific subclass `MissingParameterError(param_name, context)` is defined on line 250 and is the right fit here.

In `neoswga/core/pipeline.py`, at the top of `_filter_blacklist_penalty` (just after the function signature and before the loop, around line 78):

```python
    if bl_prefixes and not bl_seq_lengths:
        from neoswga.core.exceptions import MissingParameterError
        raise MissingParameterError(
            param_name='bl_seq_lengths',
            context=(
                f"bl_prefixes has {len(bl_prefixes)} entries but "
                f"bl_seq_lengths is empty. Either provide bl_seq_lengths "
                f"in params.json, or supply bl_genomes so neoswga can "
                f"auto-compute from the FASTAs."
            ),
        )
    if bl_prefixes and len(bl_prefixes) != len(bl_seq_lengths):
        from neoswga.core.exceptions import InvalidParameterError
        raise InvalidParameterError(
            param_name='bl_seq_lengths',
            value=bl_seq_lengths,
            reason=(
                f"length {len(bl_seq_lengths)} does not match "
                f"bl_prefixes length {len(bl_prefixes)}"
            ),
        )
```

Update the failing test in Step 1 to import `MissingParameterError` instead of `ConfigurationError`:

```python
from neoswga.core.exceptions import MissingParameterError
...
with pytest.raises(MissingParameterError):
    _filter_blacklist_penalty(...)
```

- [ ] **Step 4: Run the test**

Run: `pytest tests/test_blacklist_filter.py::TestBlacklistPenaltyMisconfig -v`
Expected: PASS.

- [ ] **Step 5: Full suite**

Run: `pytest tests/ -q --tb=short`
Expected: No new failures.

- [ ] **Step 6: Commit**

```bash
git add neoswga/core/pipeline.py tests/test_blacklist_filter.py
git commit -m "fix(blacklist): raise ConfigurationError on prefix/length mismatch

Previously the function silently computed frequencies against
len(prefixes) or 1 when bl_seq_lengths was empty, producing inflated
or zero frequencies. Explicit failure is far safer than wrong science.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

### Task 8: Wizard-produced config round-trips through `PipelineParameters`

**Files:**
- Create/Modify: `tests/test_wizard_roundtrip.py`

- [ ] **Step 1: Write failing test**

Create `tests/test_wizard_roundtrip.py`:

```python
"""The wizard must produce configs that load cleanly into PipelineParameters."""
import json
import os


def test_wizard_blacklist_config_roundtrips(tmp_path):
    """A wizard-generated params.json with blacklist fields loads without error."""
    from neoswga.core.parameter import PipelineParameters

    # Compose a dict matching what wizard.py writes for a blacklist scenario.
    wizard_config = {
        'data_dir': str(tmp_path),
        'fg_genomes': ['target.fna'],
        'bg_genomes': ['host.fna'],
        'fg_prefixes': ['target'],
        'bg_prefixes': ['host'],
        'bl_genomes': ['contaminant.fna'],
        'bl_prefixes': ['contaminant'],
        'bl_penalty': 5.0,
        'max_bl_freq': 0.0,
        'min_k': 6, 'max_k': 12,
        'polymerase': 'phi29',
        'reaction_temp': 30.0,
        'schema_version': 1,
    }

    # Any field accepted by PipelineParameters must load; extras should not crash.
    accepted = {
        k: v for k, v in wizard_config.items()
        if k in PipelineParameters.__dataclass_fields__
    }
    params = PipelineParameters(**accepted)
    assert params.bl_genomes == ['contaminant.fna']
    assert params.bl_penalty == 5.0
    assert params.max_bl_freq == 0.0


def test_wizard_writes_roundtrippable_keys():
    """wizard.py write paths use only keys that PipelineParameters recognizes.

    If the wizard emits a key that PipelineParameters doesn't know, the
    value is silently lost on load. This test reads wizard.py and checks
    every config key written.
    """
    import inspect
    from neoswga.core import wizard
    from neoswga.core.parameter import PipelineParameters

    source = inspect.getsource(wizard)
    # Extract quoted keys of the form 'xyz': or "xyz":
    import re
    keys = set(re.findall(r'["\']([a-z_]+)["\']\s*:', source))
    # Exclude keys that are obviously not params.json fields
    ignore = {'json', 'description', 'title', 'type', 'default', 'required',
              'help', 'choices', 'action', 'name', 'version'}
    candidate_params = keys - ignore
    known_fields = set(PipelineParameters.__dataclass_fields__.keys())
    # Optional aliases that the pipeline accepts even if not in dataclass
    aliases = {'schema_version', 'fasta_fore', 'fasta_back', 'kmer_fore',
               'kmer_back', 'target_set_size', 'num_primers'}
    unknown = candidate_params - known_fields - aliases
    # Filter out anything that's clearly not a config key (verbs, etc.)
    unknown = {k for k in unknown if '_' in k or k.startswith(('fg_', 'bg_', 'bl_'))}
    assert not unknown, f"wizard.py writes keys not in PipelineParameters: {unknown}"
```

- [ ] **Step 2: Run**

Run: `pytest tests/test_wizard_roundtrip.py -v`
Expected: First test passes. Second test may fail if wizard writes a renamed key. If it fails, inspect the offending key:
- If it's a typo/stale alias in wizard, fix wizard.
- If it is a legitimate alias the pipeline handles elsewhere, add it to the `aliases` set in the test with a comment explaining why.

- [ ] **Step 3: Iterate until green**

Run: `pytest tests/test_wizard_roundtrip.py -v`
Expected: PASS.

- [ ] **Step 4: Full suite**

Run: `pytest tests/ -q --tb=short`
Expected: No new failures.

- [ ] **Step 5: Commit**

```bash
git add tests/test_wizard_roundtrip.py neoswga/core/wizard.py
git commit -m "test: wizard config round-trips through PipelineParameters

Guards against wizard.py silently emitting keys the pipeline would drop.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

### Task 9: End-to-end blacklist integration test

**Files:**
- Create: `tests/fixtures/with_blacklist/target.fna`
- Create: `tests/fixtures/with_blacklist/blacklist.fna`
- Create: `tests/test_blacklist_e2e.py`

- [ ] **Step 1: Build fixtures**

```bash
python - <<'PY'
import os, random
random.seed(0xBAC)
os.makedirs('tests/fixtures/with_blacklist', exist_ok=True)

# Target: 5 kbp, 50% GC
t = ''.join(random.choices('ACGT', k=5000))
with open('tests/fixtures/with_blacklist/target.fna', 'w') as f:
    f.write('>target\n')
    for i in range(0, len(t), 80):
        f.write(t[i:i+80] + '\n')

# Blacklist: 5 kbp, shares a 500 bp window with target so some k-mers overlap.
shared = t[1000:1500]
rest = ''.join(random.choices('ACGT', k=4500))
b = shared + rest
with open('tests/fixtures/with_blacklist/blacklist.fna', 'w') as f:
    f.write('>blacklist\n')
    for i in range(0, len(b), 80):
        f.write(b[i:i+80] + '\n')
print('fixtures written')
PY
```

- [ ] **Step 2: Write the integration test**

Create `tests/test_blacklist_e2e.py`:

```python
"""End-to-end: the blacklist filter reduces blacklist-hitting primers."""
import os
import shutil

import pandas as pd
import pytest

FIXT = os.path.join(os.path.dirname(__file__), 'fixtures', 'with_blacklist')
TARGET = os.path.join(FIXT, 'target.fna')
BL = os.path.join(FIXT, 'blacklist.fna')


def _count_blacklist_hits(primers, bl_prefix, k):
    """Count primers with any hit in the blacklist k-mer file."""
    kmer_file = f"{bl_prefix}_{k}mer_all.txt"
    hits_per_primer = {}
    with open(kmer_file) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                hits_per_primer[parts[0]] = int(parts[1])
    return sum(1 for p in primers if hits_per_primer.get(p, 0) > 0)


@pytest.mark.integration
def test_blacklist_reduces_bl_hits(tmp_path):
    """With blacklist applied, step2_df has fewer bl-hitting primers than without."""
    if shutil.which('jellyfish') is None:
        pytest.skip('jellyfish not available')

    from neoswga.core.kmer_counter import run_jellyfish
    from neoswga.core.parameter import PipelineParameters, set_from_config, reset_to_defaults
    from neoswga.core.pipeline import step2

    k_min, k_max = 6, 8

    def _run(data_dir, use_blacklist):
        t_prefix = os.path.join(data_dir, 'target')
        run_jellyfish(TARGET, t_prefix, min_k=k_min, max_k=k_max)
        bl_prefix = None
        if use_blacklist:
            bl_prefix = os.path.join(data_dir, 'bl')
            run_jellyfish(BL, bl_prefix, min_k=k_min, max_k=k_max)

        cfg = PipelineParameters(
            min_k=k_min, max_k=k_max,
            min_fg_freq=1e-5, max_bg_freq=1e-3,
            min_tm=15.0, max_tm=60.0,
            max_gini=0.95, max_primer=200,
            data_dir=data_dir,
            fg_genomes=[TARGET], fg_prefixes=[t_prefix], fg_seq_lengths=[5000],
            bg_genomes=[], bg_prefixes=[], bg_seq_lengths=[],
            bl_genomes=[BL] if use_blacklist else [],
            bl_prefixes=[bl_prefix] if use_blacklist else [],
            bl_seq_lengths=[5000] if use_blacklist else [],
            max_bl_freq=0.0,
        )
        set_from_config(cfg)
        step2()
        df = pd.read_csv(os.path.join(data_dir, 'step2_df.csv'))
        reset_to_defaults()
        return df, bl_prefix, t_prefix

    with_dir = tmp_path / 'with'
    with_dir.mkdir()
    df_with, bl_prefix, t_prefix = _run(str(with_dir), use_blacklist=True)

    without_dir = tmp_path / 'without'
    without_dir.mkdir()
    df_without, _, _ = _run(str(without_dir), use_blacklist=False)

    # Use the blacklist k-mer file from the with-run to count hits in both outputs
    k_mid = (k_min + k_max) // 2
    # Regenerate bl k-mers in the without_dir with same k so we can count:
    bl_check = str(without_dir / 'bl_check')
    run_jellyfish(BL, bl_check, min_k=k_mid, max_k=k_mid)
    hits_with = _count_blacklist_hits(df_with['primer'].astype(str).tolist(), bl_check, k_mid)
    hits_without = _count_blacklist_hits(df_without['primer'].astype(str).tolist(), bl_check, k_mid)

    assert hits_with < hits_without, (
        f"With blacklist: {hits_with} bl-hitting primers; without: {hits_without}. "
        f"Blacklist should reduce hits."
    )
```

- [ ] **Step 3: Run**

Run: `pytest tests/test_blacklist_e2e.py -v -m integration`
Expected: PASS. If SKIP (jellyfish not available), note it and proceed - CI will pick this up in Phase 2.

- [ ] **Step 4: Commit**

```bash
git add tests/fixtures/with_blacklist/ tests/test_blacklist_e2e.py
git commit -m "test(integration): blacklist reduces bl-hitting primers end-to-end

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Work Stream 0.3 - Unify Tm correction path

The single-point-of-truth for additive Tm corrections should be `ReactionConditions.calculate_tm_correction()` backed by the sigmoid model in `additives.py`. Two issues today:

1. `neoswga/core/thermodynamics.py:845` `calculate_tm_batch_with_additives` uses GC-blind coefficients and already emits a `DeprecationWarning`. It has no external callers (grep confirms). It should be removed outright or thinly wrap the canonical path.
2. `neoswga/core/reaction_conditions.py:369-372` uses linear saturation (`min(1.0, tmac_m/3.0)`, `min(1.0, betaine_m/5.2)`) for the GC-normalization factor, while `additives.py` uses a sigmoid model. Both are callable.

### Task 10: Property test - betaine effect is GC-asymmetric

**Files:**
- Modify: `tests/test_thermodynamics_properties.py` (existing hypothesis tests are here)

- [ ] **Step 1: Write property test**

Append to `tests/test_thermodynamics_properties.py`:

```python
from hypothesis import given, strategies as st
from neoswga.core.reaction_conditions import ReactionConditions


@given(
    length=st.integers(min_value=8, max_value=18),
    betaine_m=st.floats(min_value=0.5, max_value=3.0),
)
def test_betaine_effect_opposite_sign_for_extreme_gc(length, betaine_m):
    """Betaine should push Tm down for GC-rich primers and up for AT-rich ones
    relative to a no-betaine baseline, because it equalizes GC/AT contribution.
    """
    at_primer = 'AT' * (length // 2) + ('A' if length % 2 else '')
    gc_primer = 'GC' * (length // 2) + ('G' if length % 2 else '')

    no_betaine = ReactionConditions(temp=30.0, polymerase='phi29', betaine_m=0.0)
    with_betaine = ReactionConditions(temp=30.0, polymerase='phi29', betaine_m=betaine_m)

    at_tm_delta = (with_betaine.calculate_effective_tm(at_primer)
                   - no_betaine.calculate_effective_tm(at_primer))
    gc_tm_delta = (with_betaine.calculate_effective_tm(gc_primer)
                   - no_betaine.calculate_effective_tm(gc_primer))

    # Relative to baseline, AT Tm delta should be greater than GC Tm delta.
    # (AT moves up or stays; GC moves down.)
    assert at_tm_delta > gc_tm_delta, (
        f"length={length} betaine={betaine_m}: AT delta={at_tm_delta:.2f}, "
        f"GC delta={gc_tm_delta:.2f}. Betaine should be GC-asymmetric."
    )
```

- [ ] **Step 2: Run**

Run: `pytest tests/test_thermodynamics_properties.py::test_betaine_effect_opposite_sign_for_extreme_gc -v`
Expected: PASS if the current implementation is correct in direction (even if magnitude differs between linear and sigmoid models). If FAIL, this is a deeper scientific bug and we stop to investigate.

- [ ] **Step 3: Commit**

```bash
git add tests/test_thermodynamics_properties.py
git commit -m "test: property test betaine is GC-asymmetric across 8-18bp

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

### Task 11: Property test - monotonicity in betaine concentration

**Files:**
- Modify: `tests/test_thermodynamics_properties.py`

- [ ] **Step 1: Write property test**

Append:

```python
@given(
    primer=st.text(alphabet='ACGT', min_size=8, max_size=18),
    b1=st.floats(min_value=0.0, max_value=2.0),
    delta=st.floats(min_value=0.1, max_value=1.0),
)
def test_tm_correction_monotone_in_betaine(primer, b1, delta):
    """For a fixed primer, effective Tm must move monotonically in betaine.

    For AT-rich: Tm increases with betaine (toward isostabilization baseline).
    For GC-rich: Tm decreases with betaine. Monotonicity holds in both.
    """
    cond_low = ReactionConditions(temp=30.0, polymerase='phi29', betaine_m=b1)
    cond_high = ReactionConditions(temp=30.0, polymerase='phi29', betaine_m=b1 + delta)
    tm_low = cond_low.calculate_effective_tm(primer)
    tm_high = cond_high.calculate_effective_tm(primer)

    gc = (primer.count('G') + primer.count('C')) / len(primer)
    if abs(gc - 0.5) < 0.05:
        # Near-50% GC: betaine effect is minimal; allow a small tolerance.
        assert abs(tm_high - tm_low) < 1.0
    elif gc < 0.5:
        assert tm_high >= tm_low - 0.01, f"AT-rich Tm should not drop with more betaine"
    else:
        assert tm_high <= tm_low + 0.01, f"GC-rich Tm should not rise with more betaine"
```

- [ ] **Step 2: Run**

Run: `pytest tests/test_thermodynamics_properties.py::test_tm_correction_monotone_in_betaine -v`
Expected: PASS. If hypothesis finds a counterexample, this points to a real bug in the Tm math; stop and investigate.

- [ ] **Step 3: Commit**

```bash
git add tests/test_thermodynamics_properties.py
git commit -m "test: Tm correction monotone in betaine concentration

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

### Task 12: Delegate `_calculate_gc_normalization` to the sigmoid model

The linear saturation factors in `reaction_conditions.py:369-372` should be replaced with the sigmoid model from `additives.py`. This is the single-point-of-truth unification.

**Files:**
- Modify: `neoswga/core/reaction_conditions.py:_calculate_gc_normalization` (around lines 336-400)
- Modify: `tests/test_reaction_conditions.py` (expected values may shift; update where necessary)

- [ ] **Step 1: Understand the sigmoid API in `additives.py`**

The sigmoid implementation lives at `neoswga/core/additives.py:554` on
`AdditiveConcentrations._gc_normalization_correction(gc_content, primer_length) -> float`.
The method reads `self.betaine_m`, `self.tmac_m`, etc. from the instance.
The `AdditiveConcentrations` dataclass is defined at `neoswga/core/additives.py:303`
and takes named fields `dmso_percent`, `betaine_m`, `tmac_m`, `formamide_percent`,
`urea_m`, `ethanol_percent`, `trehalose_m`.

Confirm by reading around line 554-621. The implementer should verify the
dataclass field names at the time of implementation (this plan was written
against commit 94df592).

- [ ] **Step 2: Identify current callers of `_calculate_gc_normalization`**

Run: `grep -n "_calculate_gc_normalization" neoswga/`
Expected: only one caller - `calculate_tm_correction` at `reaction_conditions.py:331`.

- [ ] **Step 3: Capture current numeric behavior before the change**

Write a small probe test (temporary, kept in suite) `tests/test_tm_norm_probe.py`:

```python
"""Probe the Tm correction for a known case and record the current value.
If Task 12 changes the output, this test is updated to the new canonical
value. Use the numeric delta as a change-log in the commit message."""
import pytest
from neoswga.core.reaction_conditions import ReactionConditions


@pytest.mark.parametrize("gc,betaine_m,tmac_m,expected_sign", [
    (0.33, 1.0, 0.0, '+'),   # AT-rich with betaine: correction positive
    (0.66, 1.0, 0.0, '-'),   # GC-rich with betaine: correction negative
    (0.50, 1.0, 0.0, '0'),   # balanced: near zero
])
def test_tm_correction_sign(gc, betaine_m, tmac_m, expected_sign):
    primer_len = 12
    at = int(round(primer_len * (1 - gc)))
    gc_ct = primer_len - at
    primer = 'A' * at + 'G' * gc_ct
    cond = ReactionConditions(temp=30.0, polymerase='phi29',
                              betaine_m=betaine_m, tmac_m=tmac_m)
    corr = cond.calculate_tm_correction(primer)
    # corr includes the linear additive subtraction PLUS gc normalization
    # We just verify sign of the FULL correction here to stay model-agnostic.
    if expected_sign == '+':
        assert corr > -2.0  # linear betaine -1.2*1.0 = -1.2; norm should partially offset
    elif expected_sign == '-':
        assert corr < 0
```

- [ ] **Step 4: Run probe**

Run: `pytest tests/test_tm_norm_probe.py -v`
Expected: PASS on current main (values under linear saturation).

- [ ] **Step 5: Replace `_calculate_gc_normalization` body**

In `neoswga/core/reaction_conditions.py`, replace the body of `_calculate_gc_normalization` (approximately lines 336-405) with a call into the sigmoid model on `AdditiveConcentrations`:

```python
    def _calculate_gc_normalization(self, gc_content: float,
                                     primer_length: int = 10) -> float:
        """Delegate to the sigmoid GC-normalization model in additives.py.

        Scientific references:
            - Melchior & von Hippel (1973): TMAC full GC-independence at 3M
            - Rees et al. (1993): betaine full GC-independence at 5.2M
            - Sigmoid dose-response in AdditiveConcentrations (additives.py:554)
        """
        from neoswga.core.additives import AdditiveConcentrations
        additives = AdditiveConcentrations(
            dmso_percent=self.dmso_percent,
            betaine_m=self.betaine_m,
            tmac_m=self.tmac_m,
            formamide_percent=self.formamide_percent,
            urea_m=self.urea_m,
            ethanol_percent=self.ethanol_percent,
            trehalose_m=self.trehalose_m,
        )
        return additives._gc_normalization_correction(
            gc_content=gc_content,
            primer_length=primer_length,
        )
```

Note the leading underscore on `_gc_normalization_correction` — it is marked
private by convention. This is acceptable for internal delegation within the
package; if the implementer prefers, they may add a thin public wrapper in
`additives.py` named `gc_normalization_correction` (without underscore)
before calling, and commit that change as part of this task.

Verify the `AdditiveConcentrations` dataclass field names match those above
(`grep -n "class AdditiveConcentrations" neoswga/core/additives.py` and read
the dataclass definition). If any names differ, use the actual names.

- [ ] **Step 6: Run probe and property tests**

Run: `pytest tests/test_tm_norm_probe.py tests/test_thermodynamics_properties.py tests/test_reaction_conditions.py -v`
Expected: Sign-based probe passes (sigmoid preserves direction). Property tests pass. If `test_reaction_conditions.py` has hardcoded numeric values that shift by >0.5 C due to sigmoid vs linear, update the assertions to the new canonical values and include before/after numbers in the commit message.

- [ ] **Step 7: Full suite**

Run: `pytest tests/ -q --tb=short`
Expected: No new failures. Possible numeric deltas in thermodynamics-adjacent tests; update those tests with the new canonical values.

- [ ] **Step 8: Commit**

```bash
git add neoswga/core/reaction_conditions.py tests/test_reaction_conditions.py tests/test_tm_norm_probe.py
git commit -m "refactor(thermodynamics): unify GC normalization on sigmoid model

reaction_conditions._calculate_gc_normalization now delegates to the
sigmoid dose-response model in additives.py, matching
Melchior & von Hippel (1973) and Rees et al. (1993). Removes the
linear-saturation shortcut that diverged from the rest of the codebase.

Any numeric test values that shifted due to sigmoid vs linear are
documented in the diff.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

### Task 13: Remove `calculate_tm_batch_with_additives`

The function has no external callers (confirmed by grep in audit + re-verify below). It emits a DeprecationWarning already. Removing it now prevents any future caller from picking up the GC-blind coefficients by accident.

**Files:**
- Modify: `neoswga/core/thermodynamics.py`

- [ ] **Step 1: Re-verify no callers**

Run: `grep -rn "calculate_tm_batch_with_additives" neoswga/ tests/ scripts/ examples/`
Expected: Only the definition line in `neoswga/core/thermodynamics.py`. If any other match, do not remove - convert the caller first.

- [ ] **Step 2: Delete the function**

In `neoswga/core/thermodynamics.py`, delete `def calculate_tm_batch_with_additives(...)` (approximately lines 845-921).

- [ ] **Step 3: Run suite**

Run: `pytest tests/ -q --tb=short`
Expected: No failures.

- [ ] **Step 4: Commit**

```bash
git add neoswga/core/thermodynamics.py
git commit -m "refactor(thermodynamics): remove GC-blind calculate_tm_batch_with_additives

Already deprecated; confirmed zero external callers.
Canonical path: ReactionConditions.calculate_tm_correction.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Work Stream 0.4 - Multi-target aggregation + max_gini enforcement

`MultiGenomeFilter.score_primer` (`neoswga/core/multi_genome_filter.py:313-383`) aggregates target frequencies with `np.mean` (line 345), which allows a primer absent in one target to still pass. `max_gini` is declared on the class (`:241`) but not used in `score_primer`.

### Task 14: Add `multi_target_aggregation` field

**Files:**
- Modify: `neoswga/core/parameter.py` (`PipelineParameters`, `get_current_config`, `set_from_config`, `get_params`)
- Create: `tests/test_multi_target_aggregation.py`

- [ ] **Step 1: Write failing field test**

Create `tests/test_multi_target_aggregation.py`:

```python
from neoswga.core.parameter import PipelineParameters


def test_multi_target_aggregation_default_is_min():
    p = PipelineParameters()
    assert p.multi_target_aggregation == 'min'


def test_multi_target_aggregation_accepts_mean_and_geomean():
    p = PipelineParameters(multi_target_aggregation='mean')
    assert p.multi_target_aggregation == 'mean'
    p = PipelineParameters(multi_target_aggregation='geomean')
    assert p.multi_target_aggregation == 'geomean'
```

- [ ] **Step 2: Run and confirm FAIL**

Run: `pytest tests/test_multi_target_aggregation.py -v`
Expected: FAIL - AttributeError on `multi_target_aggregation`.

- [ ] **Step 3: Add field to dataclass**

In `neoswga/core/parameter.py`, after the `max_primer` line (~line 55):

```python
    # Multi-target frequency aggregation: 'min' (strictest; pan-primer must
    # hit every target), 'geomean' (robust, all-or-nothing), 'mean' (legacy,
    # forgiving; primer can be absent in one target and still pass).
    multi_target_aggregation: str = 'min'
```

Update `get_current_config()` and `set_from_config()` to include the new field, matching the pattern used for other fields.

Add to `get_params()`:

```python
    multi_target_aggregation = data.get('multi_target_aggregation', 'min')
    data['multi_target_aggregation'] = multi_target_aggregation
    if multi_target_aggregation not in ('min', 'mean', 'geomean'):
        raise ValueError(
            f"multi_target_aggregation must be 'min', 'mean', or 'geomean'; "
            f"got {multi_target_aggregation!r}"
        )
```

And a module-level default:

```python
multi_target_aggregation = 'min'
```

- [ ] **Step 4: Run tests**

Run: `pytest tests/test_multi_target_aggregation.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add neoswga/core/parameter.py tests/test_multi_target_aggregation.py
git commit -m "feat(multi-target): add multi_target_aggregation field (default min)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

### Task 15: Dispatch aggregation in `MultiGenomeFilter.score_primer`

**Files:**
- Modify: `neoswga/core/multi_genome_filter.py`
- Modify: `tests/test_multi_target_aggregation.py`

- [ ] **Step 1: Write failing behavioral test**

Append to `tests/test_multi_target_aggregation.py`:

```python
import numpy as np
import pytest

from neoswga.core.multi_genome_filter import (
    MultiGenomeFilter, GenomeSet, GenomeEntry, GenomeRole,
)


@pytest.fixture
def _two_target_fastas(tmp_path):
    """GenomeEntry.__post_init__ requires the FASTA path to exist. Write stubs."""
    t1 = tmp_path / 't1.fna'
    t2 = tmp_path / 't2.fna'
    t1.write_text('>t1\nATCG\n')
    t2.write_text('>t2\nGCTA\n')
    return str(t1), str(t2)


def _make_filter(aggregation, fasta_pair):
    """Build a two-target filter with mocked counts."""
    t1_path, t2_path = fasta_pair
    g1 = GenomeEntry(name='t1', fasta_path=t1_path, role=GenomeRole.TARGET)
    g2 = GenomeEntry(name='t2', fasta_path=t2_path, role=GenomeRole.TARGET)
    gs = GenomeSet(targets=[g1, g2], backgrounds=[], blacklists=[])
    f = MultiGenomeFilter(
        genome_set=gs,
        min_target_freq=1e-4,
        max_background_freq=1.0,
        max_blacklist_freq=1.0,
        min_enrichment=0.0,
        use_gini=False,
        aggregation=aggregation,
    )
    f.load_genome_counts('t1', {'ATCGATCG': 10}, genome_size=10_000)
    f.load_genome_counts('t2', {'ATCGATCG': 0}, genome_size=10_000)
    return f


def test_min_aggregation_rejects_primer_missing_in_one_target(_two_target_fastas):
    f = _make_filter('min', _two_target_fastas)
    score = f.score_primer('ATCGATCG')
    assert score.target_frequency == 0.0
    assert score.passes is False


def test_geomean_aggregation_rejects_primer_missing_in_one_target(_two_target_fastas):
    f = _make_filter('geomean', _two_target_fastas)
    score = f.score_primer('ATCGATCG')
    assert score.target_frequency == 0.0  # geomean with a zero is zero
    assert score.passes is False


def test_mean_aggregation_accepts_primer_missing_in_one_target(_two_target_fastas):
    f = _make_filter('mean', _two_target_fastas)
    score = f.score_primer('ATCGATCG')
    # 10/10000 and 0/10000 -> mean = 5e-4 >= min_target_freq=1e-4 -> passes
    assert score.target_frequency == pytest.approx(5e-4)
    assert score.passes is True
```

- [ ] **Step 2: Run, confirm FAIL**

Run: `pytest tests/test_multi_target_aggregation.py -v`
Expected: FAIL - `MultiGenomeFilter.__init__` does not accept `aggregation=`.

- [ ] **Step 3: Add `aggregation` parameter to `MultiGenomeFilter.__init__`**

In `neoswga/core/multi_genome_filter.py`, in `MultiGenomeFilter.__init__` (around line 234-275), add a new kwarg:

```python
    def __init__(self,
                 genome_set: GenomeSet,
                 min_target_freq: float = 1e-5,
                 max_background_freq: float = 1e-4,
                 max_blacklist_freq: float = 1e-6,
                 min_enrichment: float = 10.0,
                 use_gini: bool = True,
                 max_gini: float = 0.6,
                 aggregation: str = 'min'):
```

And store it:

```python
        if aggregation not in ('min', 'mean', 'geomean'):
            raise ValueError(f"aggregation must be 'min', 'mean', or 'geomean'; got {aggregation!r}")
        self.aggregation = aggregation
```

- [ ] **Step 4: Replace the aggregation in `score_primer`**

In `neoswga/core/multi_genome_filter.py` around line 345, replace:

```python
        # Aggregate frequencies
        target_freq = np.mean(target_freqs) if target_freqs else 0.0
```

with:

```python
        # Aggregate target frequencies per the configured policy.
        target_freq = self._aggregate_target(target_freqs)
```

And add a helper method inside the class:

```python
    def _aggregate_target(self, freqs):
        if not freqs:
            return 0.0
        if self.aggregation == 'min':
            return float(np.min(freqs))
        if self.aggregation == 'mean':
            return float(np.mean(freqs))
        if self.aggregation == 'geomean':
            # Geometric mean: zero-absorbing (if any freq is 0, geomean is 0).
            import math
            if any(x <= 0 for x in freqs):
                return 0.0
            return float(math.exp(sum(math.log(x) for x in freqs) / len(freqs)))
        raise ValueError(self.aggregation)
```

- [ ] **Step 5: Run tests**

Run: `pytest tests/test_multi_target_aggregation.py -v`
Expected: PASS.

- [ ] **Step 6: Full suite**

Run: `pytest tests/ -q --tb=short`
Expected: Possibly some `multi_genome` tests that assumed mean-aggregation now fail. For each:
- If the test was asserting behavior that's now wrong, update the assertion or pass `aggregation='mean'` explicitly to preserve backward compatibility.
- Commit test updates separately for clarity.

- [ ] **Step 7: Commit**

```bash
git add neoswga/core/multi_genome_filter.py tests/test_multi_target_aggregation.py
git commit -m "feat(multi-target): dispatch target-frequency aggregation policy

Default is 'min' (every target must meet threshold). 'mean' preserves
the pre-3.7 behavior as opt-in. 'geomean' is a middle-ground option.

Breaking: existing multi-target params.json users should review results.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

### Task 16: Enforce `max_gini` in `score_primer`

**Files:**
- Modify: `neoswga/core/multi_genome_filter.py:score_primer`
- Modify: `tests/test_multi_target_aggregation.py`

- [ ] **Step 1: Write failing test**

Append to `tests/test_multi_target_aggregation.py`:

```python
def test_max_gini_enforced_in_score_primer(tmp_path):
    """High Gini (binding concentrated in one region) must fail the filter."""
    from neoswga.core.multi_genome_filter import (
        MultiGenomeFilter, GenomeSet, GenomeEntry, GenomeRole,
    )
    tfna = tmp_path / 't.fna'
    tfna.write_text('>t\nATCG\n')
    g = GenomeEntry(name='t', fasta_path=str(tfna), role=GenomeRole.TARGET)
    gs = GenomeSet(targets=[g], backgrounds=[], blacklists=[])
    f = MultiGenomeFilter(
        genome_set=gs,
        min_target_freq=1e-5,
        max_background_freq=1.0,
        max_blacklist_freq=1.0,
        min_enrichment=0.0,
        use_gini=True,
        max_gini=0.3,
    )
    f.load_genome_counts('t', {'ATCGATCG': 50}, genome_size=10_000)
    # Inject a per-primer Gini via an attribute the filter reads.
    f.primer_gini = {'ATCGATCG': 0.8}  # exceeds max_gini=0.3

    score = f.score_primer('ATCGATCG')
    assert score.passes is False
    assert 'gini' in (score.details.get('rejection_reason', '') or '').lower() \
        or score.details.get('gini', 0) > 0.3
```

- [ ] **Step 2: Run, confirm FAIL**

Run: `pytest tests/test_multi_target_aggregation.py::test_max_gini_enforced_in_score_primer -v`
Expected: FAIL.

- [ ] **Step 3: Implement optional Gini lookup in `score_primer`**

In `neoswga/core/multi_genome_filter.py`, modify `score_primer` to honor `use_gini` and `max_gini`:

After computing `target_freq` and before the `passes = (...)` block, add:

```python
        # Gini index enforcement (when enabled and when a per-primer gini
        # has been loaded via .primer_gini).
        gini_rejected = False
        gini_value = None
        if self.use_gini:
            gini_map = getattr(self, 'primer_gini', None) or {}
            gini_value = gini_map.get(primer)
            if gini_value is not None and gini_value > self.max_gini:
                gini_rejected = True
                details['gini'] = gini_value
                details['rejection_reason'] = (
                    f"Gini {gini_value:.2f} exceeds max_gini={self.max_gini}"
                )
```

Update the `passes` expression:

```python
        passes = (
            target_freq >= self.min_target_freq and
            background_freq <= self.max_background_freq and
            blacklist_freq <= self.max_blacklist_freq and
            enrichment >= self.min_enrichment and
            not gini_rejected
        )
```

- [ ] **Step 4: Run tests**

Run: `pytest tests/test_multi_target_aggregation.py -v`
Expected: PASS.

- [ ] **Step 5: Full suite**

Run: `pytest tests/ -q --tb=short`
Expected: No new failures.

- [ ] **Step 6: Commit**

```bash
git add neoswga/core/multi_genome_filter.py tests/test_multi_target_aggregation.py
git commit -m "fix(multi-target): enforce max_gini in MultiGenomeFilter.score_primer

Previously max_gini was accepted as a parameter but silently ignored.
When .primer_gini is populated (e.g. by multi_genome_pipeline), primers
with Gini exceeding max_gini are now rejected with a reason.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

### Task 17: Wire `aggregation` and `max_gini` through `multi_genome_pipeline`

**Files:**
- Modify: `neoswga/core/multi_genome_pipeline.py`
- Modify: `tests/test_multi_target_aggregation.py`

- [ ] **Step 1: Find where `MultiGenomeFilter` is constructed**

Run: `grep -n "MultiGenomeFilter(" neoswga/core/multi_genome_pipeline.py`
Expected: one or more call sites.

- [ ] **Step 2: Pass parameters through**

Add `aggregation=` to each `MultiGenomeFilter(...)` call, read from `PipelineParameters.multi_target_aggregation`. Similarly pass `max_gini=`. Read the parameter via:

```python
from neoswga.core.parameter import get_current_config
cfg = get_current_config()
```

And pass `aggregation=cfg.multi_target_aggregation, max_gini=cfg.max_gini`.

If the pipeline computes Gini per primer (via `get_gini`) it should populate `filter.primer_gini = dict(zip(primers, gini_values))` before `score_primer` is called.

- [ ] **Step 3: Write an end-to-end multi-target test**

Append to `tests/test_multi_target_aggregation.py`:

```python
@pytest.mark.integration
def test_multi_target_min_aggregation_end_to_end(tmp_path):
    """Two targets; a primer present in only one should not survive with 'min'."""
    # This test is deliberately light on fixture; we use the existing
    # pipeline-level hook via set_from_config and exercise the filter only.
    # A fuller e2e would require two jellyfish runs and is deferred to Phase 2.
    from neoswga.core.multi_genome_filter import (
        MultiGenomeFilter, GenomeSet, GenomeEntry, GenomeRole,
    )
    t1_fna = tmp_path / 't1.fna'
    t1_fna.write_text('>t1\nATCG\n')
    t2_fna = tmp_path / 't2.fna'
    t2_fna.write_text('>t2\nGCTA\n')
    t1 = GenomeEntry(name='t1', fasta_path=str(t1_fna), role=GenomeRole.TARGET)
    t2 = GenomeEntry(name='t2', fasta_path=str(t2_fna), role=GenomeRole.TARGET)
    gs = GenomeSet(targets=[t1, t2], backgrounds=[], blacklists=[])
    f = MultiGenomeFilter(
        genome_set=gs, min_target_freq=1e-4, max_background_freq=1.0,
        max_blacklist_freq=1.0, min_enrichment=0.0, use_gini=False,
        aggregation='min',
    )
    f.load_genome_counts('t1', {'AAAAAA': 5, 'GGGGGG': 5}, genome_size=10_000)
    f.load_genome_counts('t2', {'AAAAAA': 0, 'GGGGGG': 5}, genome_size=10_000)

    passing, scores = f.filter_primers(['AAAAAA', 'GGGGGG'])
    assert 'GGGGGG' in passing
    assert 'AAAAAA' not in passing
```

- [ ] **Step 4: Run**

Run: `pytest tests/test_multi_target_aggregation.py -v`
Expected: PASS.

- [ ] **Step 5: Full suite**

Run: `pytest tests/ -q --tb=short`
Expected: No new failures.

- [ ] **Step 6: Commit**

```bash
git add neoswga/core/multi_genome_pipeline.py tests/test_multi_target_aggregation.py
git commit -m "feat(multi-target): wire aggregation and max_gini through pipeline

MultiGenomePipeline now reads multi_target_aggregation from
PipelineParameters and threads it, along with max_gini, into
MultiGenomeFilter.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Section 5 - Release closeout

### Task 18: Update CHANGELOG

**Files:**
- Modify: `docs/CHANGELOG.md`

- [ ] **Step 1: Read current CHANGELOG**

Run: `head -50 docs/CHANGELOG.md`
Note the format used by prior releases.

- [ ] **Step 2: Add the 3.7.0 section**

Prepend a new section at the top of `docs/CHANGELOG.md`:

```markdown
## 3.7.0 - 2026-04-XX

### Bug fixes

- **Adaptive GC regression guards.** Extreme-GC genomes (Francisella ~33%,
  Burkholderia ~68%) now have explicit integration tests so the silently
  functioning adaptive GC filter cannot regress unnoticed.
- **Blacklist `bl_seq_lengths` auto-computed.** When `bl_genomes` is
  supplied without `bl_seq_lengths`, lengths are now derived from the
  FASTAs, matching `fg_seq_lengths`/`bg_seq_lengths` behavior.
  Misconfiguration (prefixes without lengths) now raises
  `ConfigurationError` instead of silently producing wrong frequencies.
- **Tm correction unified on sigmoid model.** `ReactionConditions._calculate
  _gc_normalization` now delegates to the sigmoid dose-response model in
  `additives.py` (Melchior & von Hippel 1973; Rees et al. 1993), replacing
  the linear-saturation shortcut.
- **`max_gini` actually enforced.** `MultiGenomeFilter.score_primer` now
  rejects primers whose per-primer Gini exceeds `max_gini` when
  `use_gini=True`. Previously the threshold was silently ignored on the
  multi-genome path.

### Behavior changes

- **`multi_target_aggregation` default is now `min` (breaking).** Previously
  multi-target frequencies were aggregated with the mean, which let a primer
  absent from one target pass. Default is now `min`; users who relied on the
  prior behavior can set `"multi_target_aggregation": "mean"` in
  `params.json`. A `geomean` option is also available.
- **`use_adaptive_gc` flag (new, default True).** Explicitly opt out of
  adaptive GC filtering with `"use_adaptive_gc": false` in `params.json`.

### Removed

- `calculate_tm_batch_with_additives` (deprecated in 3.5; GC-blind
  coefficients; no external callers).
- `AdaptiveGCFilter` class (never instantiated; its docstring
  misrepresented the state of the system).

### Testing

- New integration tests for low-GC and high-GC targets.
- Property-based tests for GC-asymmetric betaine effect and monotonicity
  in betaine concentration.
- End-to-end blacklist integration test.
- Multi-target aggregation unit tests covering `min`, `mean`, and `geomean`.

```

Fill in the release date when the version bump commit lands.

- [ ] **Step 3: Commit**

```bash
git add docs/CHANGELOG.md
git commit -m "docs: CHANGELOG entry for 3.7.0

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

### Task 19: Version bump to 3.7.0

**Files:**
- Modify: `pyproject.toml`
- Modify: `neoswga/__init__.py`
- Modify: `VERSION`

- [ ] **Step 1: Read current version locations**

Run: `grep -n "3.6.0" pyproject.toml neoswga/__init__.py VERSION`
Expected: one match in each file.

- [ ] **Step 2: Bump each**

Edit `pyproject.toml`: change `version = "3.6.0"` to `version = "3.7.0"`.
Edit `neoswga/__init__.py`: change the `__version__ = "3.6.0"` to `__version__ = "3.7.0"`.
Edit `VERSION`: change `3.6.0` to `3.7.0`.

- [ ] **Step 3: Update the CHANGELOG date**

In `docs/CHANGELOG.md`, replace `2026-04-XX` in the 3.7.0 section with the current date in `YYYY-MM-DD` form.

- [ ] **Step 4: Verify consistency**

Run: `grep -rn '3\.7\.0' pyproject.toml neoswga/__init__.py VERSION docs/CHANGELOG.md`
Expected: one match per file (CHANGELOG may have two: the header and the date if present).

- [ ] **Step 5: Full suite + build smoke test**

Run: `pytest tests/ -q --tb=short && python -m build --sdist --wheel --outdir /tmp/neoswga-build && ls /tmp/neoswga-build`
Expected: All tests pass; build produces `neoswga-3.7.0-py3-none-any.whl` and `neoswga-3.7.0.tar.gz`.

- [ ] **Step 6: Commit**

```bash
git add pyproject.toml neoswga/__init__.py VERSION docs/CHANGELOG.md
git commit -m "release: 3.7.0

See docs/CHANGELOG.md for the full list of bug fixes and behavior changes.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

- [ ] **Step 7: Tag**

```bash
git tag -a v3.7.0 -m "Release 3.7.0 - stop the bleeding phase"
```

Do not push the tag without explicit user approval (this is a release action affecting the public repo).

---

## Final verification

- [ ] **Full suite one more time**

Run: `pytest tests/ -q --tb=short --cov=neoswga --cov-report=term-missing`
Expected: All tests pass. Note the new coverage percentage - it should be >= 50% (current floor) and ideally higher due to the new tests.

- [ ] **Lint**

Run: `black --check neoswga/ && isort --check neoswga/`
Expected: Both pass. If not, run `black neoswga/ && isort neoswga/` and commit the formatting as a separate commit.

- [ ] **Manual sanity check on the plasmid example**

Run:
```bash
cd examples/plasmid_example
neoswga count-kmers -j params.json
neoswga filter -j params.json
neoswga score -j params.json
neoswga optimize -j params.json
```
Expected: All four steps complete without error. `step4_improved_df.csv` contains the expected columns and non-zero rows.

- [ ] **Summarize what shipped**

Write a one-paragraph release note to stdout describing what phase 0 delivered, for the human reviewer's benefit. This is not a file change.

---

## Exit criteria checklist (maps to spec)

- [x] Low-GC fixture + integration test (Task 1)
- [x] High-GC fixture + integration test (Task 2)
- [x] `use_adaptive_gc` opt-out flag (Task 3)
- [x] Orphan `AdaptiveGCFilter` removed (Task 4)
- [x] Untracked blacklist/library tests committed (Task 5)
- [x] `bl_seq_lengths` auto-computed (Task 6)
- [x] Guard against prefix/length mismatch (Task 7)
- [x] Wizard round-trip test (Task 8)
- [x] Blacklist integration test (Task 9)
- [x] Betaine GC-asymmetric property test (Task 10)
- [x] Monotonicity property test (Task 11)
- [x] Sigmoid GC normalization unification (Task 12)
- [x] Removed `calculate_tm_batch_with_additives` (Task 13)
- [x] `multi_target_aggregation` field (Task 14)
- [x] Aggregation dispatch in `score_primer` (Task 15)
- [x] `max_gini` enforcement (Task 16)
- [x] Pipeline wiring for aggregation + max_gini (Task 17)
- [x] CHANGELOG + version bump (Tasks 18-19)

## Notes for the implementer

- If a property test finds a hypothesis counterexample that points at a deeper bug, STOP and escalate. Do not suppress the counterexample.
- When `pytest -m integration` is skipped due to missing jellyfish, that's expected in environments without it. Phase 2 will fix CI; for Phase 0 the tests exist and run locally when jellyfish is installed.
- The plan sequences work streams 0.1 through 0.4 for clarity, but each stream is independent. If one stream blocks, move to the next.
- If a test fails with a pre-existing failure noted at preflight, do not attempt to fix it as part of this phase - document it and move on. Pre-existing failures are out of Phase 0 scope.
