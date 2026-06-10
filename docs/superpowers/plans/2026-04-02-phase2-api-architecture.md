# Phase 2: API & Architecture Remediation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix API stability, version management, logging hygiene, and two targeted code bugs before PyPI publication.

**Architecture:** Targeted fixes to existing modules. Tasks are independent. The print-to-logger conversion is the largest task (1028 print() calls across 64 files) but most are in __main__ guards or interactive UX modules where print() is intentional.

**Tech Stack:** Python 3.11+, pytest, logging

---

### Task 1: Single-Source Version Number

Version "3.0.0" is hardcoded in 5+ places. Use importlib.metadata for report modules. Fix changelog/version mismatch.

**Files:**
- Modify: `neoswga/core/report/metrics.py` (~line 251)
- Modify: `neoswga/core/report/executive_summary.py` (~line 51)
- Modify: `neoswga/core/report/technical_report.py` (~line 103)

- [ ] **Step 1: Find all hardcoded version strings in report modules**

Run: `grep -rn '"3.0.0"' neoswga/core/report/`

- [ ] **Step 2: Replace hardcoded versions with dynamic import**

In each report module that has `version: str = "3.0.0"` as a dataclass default, change it to use a helper. Add this function to `neoswga/core/report/utils.py`:

```python
def get_version() -> str:
    """Get NeoSWGA version from package metadata."""
    try:
        from neoswga import __version__
        return __version__
    except ImportError:
        return "unknown"
```

Then in each dataclass, use `field(default_factory=...)` or set it in `__post_init__`. The simplest approach: change the default to call the function where the version is assigned (not in the dataclass field default, but where instances are created).

Read each file first to determine the best insertion point.

- [ ] **Step 3: Fix changelog version — set __version__ to match CHANGELOG**

The CHANGELOG documents through 3.6.0 but __version__ is 3.0.0. Update `neoswga/__init__.py` and `pyproject.toml` to `3.6.0` to match reality.

- [ ] **Step 4: Verify**

Run: `grep -rn '"3.0.0"' neoswga/` — should only appear in test fixtures, not production code.
Run: `python -c "import neoswga; print(neoswga.__version__)"` — should print 3.6.0.

- [ ] **Step 5: Run tests and commit**

Run: `pytest tests/ -x --timeout=120 -q`
Commit: "fix: single-source version number, bump to 3.6.0 to match changelog"

---

### Task 2: Add params.json Schema Versioning

Add a `schema_version` field to params.json. Emit a warning when loading unversioned configs.

**Files:**
- Modify: `neoswga/core/parameter.py` (in `get_params()`)
- Modify: `neoswga/core/wizard.py` (in init command output)

- [ ] **Step 1: Read parameter.py get_params() to find where JSON is loaded**

Read `neoswga/core/parameter.py` around the `get_params()` function.

- [ ] **Step 2: Add schema version check after JSON load**

After `data = json.load(...)` in `get_params()`, add:

```python
    # Schema versioning
    CURRENT_SCHEMA_VERSION = 1
    schema_version = data.get('schema_version', None)
    if schema_version is None:
        logger.warning(
            "params.json has no 'schema_version' field. "
            "Defaults may differ between NeoSWGA versions. "
            "Add '\"schema_version\": 1' to your params.json for reproducibility."
        )
    elif schema_version > CURRENT_SCHEMA_VERSION:
        logger.warning(
            f"params.json schema_version {schema_version} is newer than "
            f"this NeoSWGA version supports (max: {CURRENT_SCHEMA_VERSION}). "
            f"Some parameters may not be recognized."
        )
```

- [ ] **Step 3: Update wizard.py to include schema_version in generated configs**

Find where the wizard writes params.json and add `"schema_version": 1` to the output dict.

- [ ] **Step 4: Run tests and commit**

Run: `pytest tests/ -x --timeout=120 -q`
Commit: "feat: add params.json schema versioning for reproducibility"

---

### Task 3: Fix Stale Reaction Conditions Cache

`filter.py` caches `ReactionConditions` at module level. Adaptive GC filtering changes parameters but doesn't reset the cache.

**Files:**
- Modify: `neoswga/core/adaptive_filters.py` (in `run_step2_with_adaptive_gc()`)

- [ ] **Step 1: Read adaptive_filters.py to find run_step2_with_adaptive_gc**

Find where parameters are modified and where step2() is called.

- [ ] **Step 2: Add reset call before step2()**

Before the call to `step2()` (or equivalent filter function), add:

```python
from neoswga.core.filter import reset_reaction_conditions
reset_reaction_conditions()
```

Also add it in the `finally` block after parameters are restored.

- [ ] **Step 3: Write test**

Create `tests/test_reaction_conditions_cache.py`:
```python
"""Test that reaction conditions cache is properly invalidated."""

def test_reset_reaction_conditions_clears_cache():
    from neoswga.core.filter import _get_reaction_conditions, reset_reaction_conditions
    
    # First call creates cache
    rc1 = _get_reaction_conditions()
    # Same object returned
    rc2 = _get_reaction_conditions()
    assert rc1 is rc2
    
    # Reset clears cache
    reset_reaction_conditions()
    rc3 = _get_reaction_conditions()
    assert rc3 is not rc1
```

- [ ] **Step 4: Run tests and commit**

Run: `pytest tests/test_reaction_conditions_cache.py tests/ -x --timeout=120 -q`
Commit: "fix: reset reaction conditions cache during adaptive GC filtering"

---

### Task 4: Fix BackgroundAwareOptimizer Multi-Genome Coverage

`_calculate_coverage()` only uses `fg_prefixes[0]`, ignoring additional foreground genomes.

**Files:**
- Modify: `neoswga/core/background_aware_optimizer.py` (~line 358-361)

- [ ] **Step 1: Read the _calculate_coverage method**

Read `neoswga/core/background_aware_optimizer.py` around line 350-380.

- [ ] **Step 2: Fix to iterate over all foreground prefixes**

Replace:
```python
        all_positions = []
        for primer in primers:
            positions = self.cache.get_positions(self.fg_prefixes[0], primer, 'both')
            all_positions.extend(positions)
```

With:
```python
        all_positions = []
        for primer in primers:
            for prefix in self.fg_prefixes:
                positions = self.cache.get_positions(prefix, primer, 'both')
                all_positions.extend(positions)
```

- [ ] **Step 3: Run tests and commit**

Run: `pytest tests/ -x --timeout=120 -q`
Commit: "fix: BackgroundAwareOptimizer now uses all foreground genomes for coverage"

---

### Task 5: Replace print() with logger in Pipeline-Critical Modules

1028 print() calls across 64 files. However, many are in `__main__` guards, interactive UX (wizard, workflow_selector, condition_suggester, results_interpreter), or simulation plot output where print() is appropriate. Focus on the pipeline-critical modules where print() pollutes stdout for library users.

**Target files** (modules that run during normal pipeline execution):
- `neoswga/core/genetic_algorithm.py` (19 prints — GA optimizer output)
- `neoswga/core/dimer.py` (3 prints — debug output in production path at lines 29-30, 327-329)
- `neoswga/core/improved_pipeline.py` (28 prints)
- `neoswga/core/milp_optimizer.py` (19 prints)
- `neoswga/core/dominating_set_optimizer.py` (check count)
- `neoswga/core/network_optimizer.py` (check count)
- `neoswga/core/hybrid_optimizer.py` (check count)
- `neoswga/core/background_aware_optimizer.py` (check count)
- `neoswga/core/equiphi29_optimizer.py` (check count)
- `neoswga/core/moea_optimizer.py` (check count)
- `neoswga/core/rf_preprocessing.py` (check count)
- `neoswga/core/genome_io.py` (30 prints — but check if in __main__)

**DO NOT convert** print() in these intentional-UI modules:
- `wizard.py` (interactive prompts)
- `workflow_selector.py` (interactive menu)
- `condition_suggester.py` (user-facing suggestions)
- `results_interpreter.py` (user-facing output)
- `param_validator.py` (user-facing validation output)
- `progress.py` (progress bars)
- Any `if __name__ == "__main__":` blocks

**Rules:**
- `print(f"...")` with status/progress info → `logger.info(f"...")`
- `print(f"...")` with debug/diagnostic info → `logger.debug(f"...")`
- `print(...)` that is clearly debug leftovers (like dimer.py lines 29-30) → delete entirely or `logger.debug()`
- If the module doesn't have `logger = logging.getLogger(__name__)`, add it near the top after imports

- [ ] **Step 1: Convert genetic_algorithm.py**

Replace all 19 print() calls with logger.info()/logger.debug(). The module imports logging but doesn't create a logger — add `logger = logging.getLogger(__name__)` after imports.

- [ ] **Step 2: Fix dimer.py debug prints**

Lines 29-30: delete or convert to logger.debug(). Lines 327-329: these are in a `__main__` block — leave them.

- [ ] **Step 3: Convert remaining optimizer modules**

For each optimizer file, replace print() in production code paths with logger calls. Leave __main__ blocks alone.

- [ ] **Step 4: Convert improved_pipeline.py**

Replace print() in production code paths.

- [ ] **Step 5: Convert genome_io.py**

Only convert print() outside __main__ blocks.

- [ ] **Step 6: Convert rf_preprocessing.py**

Replace any remaining print() in production paths.

- [ ] **Step 7: Verify no print() in pipeline path**

Run: `grep -rn 'print(' neoswga/core/genetic_algorithm.py neoswga/core/dimer.py neoswga/core/improved_pipeline.py neoswga/core/milp_optimizer.py neoswga/core/dominating_set_optimizer.py neoswga/core/network_optimizer.py neoswga/core/hybrid_optimizer.py neoswga/core/background_aware_optimizer.py neoswga/core/equiphi29_optimizer.py neoswga/core/moea_optimizer.py neoswga/core/rf_preprocessing.py neoswga/core/genome_io.py | grep -v '__main__' | grep -v '# print'`

Should return minimal results (only __main__ blocks and intentional user output).

- [ ] **Step 8: Run tests and commit**

Run: `pytest tests/ -x --timeout=120 -q`
Commit: "refactor: replace print() with logger in pipeline-critical modules"

---

### Task 6: Define Public API Surface

Add `__all__` to key package init files to declare the stable public API.

**Files:**
- Modify: `neoswga/__init__.py`
- Modify: `neoswga/core/__init__.py`

- [ ] **Step 1: Define top-level __all__**

In `neoswga/__init__.py`, add:
```python
__all__ = [
    "__version__",
    "get_version",
]
```

- [ ] **Step 2: Define core __all__ with key public classes**

In `neoswga/core/__init__.py`, add lazy imports and __all__ for the most important public symbols:

```python
"""NeoSWGA core modules for primer design and optimization."""

__all__ = [
    "ReactionConditions",
    "PositionCache",
    "MechanisticModel",
    "MechanisticEffects",
    "get_standard_conditions",
    "get_enhanced_conditions",
    "calculate_tm_with_salt",
    "calculate_tm_basic",
]


def __getattr__(name):
    """Lazy imports for public API symbols."""
    if name in ("ReactionConditions", "get_standard_conditions", "get_enhanced_conditions"):
        from neoswga.core.reaction_conditions import ReactionConditions, get_standard_conditions, get_enhanced_conditions
        return locals()[name]
    if name in ("PositionCache",):
        from neoswga.core.position_cache import PositionCache
        return PositionCache
    if name in ("MechanisticModel", "MechanisticEffects"):
        from neoswga.core.mechanistic_model import MechanisticModel, MechanisticEffects
        return locals()[name]
    if name in ("calculate_tm_with_salt", "calculate_tm_basic"):
        from neoswga.core.thermodynamics import calculate_tm_with_salt, calculate_tm_basic
        return locals()[name]
    raise AttributeError(f"module 'neoswga.core' has no attribute {name!r}")
```

- [ ] **Step 3: Write test**

Create `tests/test_public_api.py`:
```python
"""Test that the public API is importable and stable."""

def test_top_level_exports():
    import neoswga
    assert hasattr(neoswga, '__version__')
    assert hasattr(neoswga, 'get_version')

def test_core_lazy_imports():
    from neoswga.core import ReactionConditions
    from neoswga.core import PositionCache
    from neoswga.core import MechanisticModel
    from neoswga.core import calculate_tm_with_salt
    assert callable(calculate_tm_with_salt)

def test_core_all_is_defined():
    import neoswga.core
    assert hasattr(neoswga.core, '__all__')
    assert 'ReactionConditions' in neoswga.core.__all__
```

- [ ] **Step 4: Run tests and commit**

Run: `pytest tests/test_public_api.py tests/ -x --timeout=120 -q`
Commit: "feat: define public API surface with __all__ and lazy imports"

---

### Task 7: Document Deprecation Policy

Add a deprecation policy section to README.md.

**Files:**
- Modify: `README.md` (add section before or after Contributing)

- [ ] **Step 1: Add deprecation policy**

Add this section to README.md:

```markdown
## Deprecation Policy

NeoSWGA follows semantic versioning. When features are deprecated:

- **Deprecated features** emit a `DeprecationWarning` for at least one minor release before removal.
- **Removed features** are documented in the CHANGELOG with migration guidance.
- **params.json changes** are backwards compatible within the same major version. New parameters use sensible defaults.
- **CLI flag changes** follow the same deprecation cycle: warning first, removal in next major.
```

- [ ] **Step 2: Commit**

Commit: "docs: add deprecation policy to README"

---

## Summary

| Task | Issue | Severity | Est. Time |
|------|-------|----------|-----------|
| 1 | Single-source version + changelog fix | High | 20 min |
| 2 | params.json schema versioning | Critical | 15 min |
| 3 | Stale reaction conditions cache | High | 10 min |
| 4 | BackgroundAwareOptimizer multi-genome | High | 10 min |
| 5 | print() → logger in pipeline modules | Medium | 45 min |
| 6 | Define public API surface | Critical | 20 min |
| 7 | Deprecation policy docs | High | 5 min |

**Total estimated: ~2 hours**
