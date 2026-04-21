"""Phase 16 critical gap #3: library users can import the public API
directly from ``neoswga``, not just from deep submodules.

Before the fix, a Python user had to write
``from neoswga.core.unified_optimizer import run_optimization`` which
leaked implementation layout. Now ``neoswga.__all__`` lists the public
API and ``from neoswga import run_optimization`` works via PEP 562
lazy imports.
"""

import pytest


PUBLIC_API = [
    "run_optimization",
    "optimize_step4",
    "OptimizationResult",
    "OptimizationStatus",
    "PrimerSetMetrics",
    "OptimizerFactory",
    "ReactionConditions",
    "get_polymerase_processivity",
    "get_typical_amplicon_length",
    "compute_per_prefix_coverage",
    "polymerase_extension_reach",
    "PositionCache",
]


@pytest.mark.parametrize("name", PUBLIC_API)
def test_public_api_importable_from_top_level(name):
    """`from neoswga import <name>` must resolve for every public symbol."""
    import importlib
    neoswga = importlib.import_module("neoswga")
    symbol = getattr(neoswga, name)
    assert symbol is not None, f"neoswga.{name} resolved to None"


def test_all_dunder_contains_public_api():
    import neoswga
    missing = [n for n in PUBLIC_API if n not in neoswga.__all__]
    assert not missing, (
        f"neoswga.__all__ missing public symbols: {missing}"
    )


def test_dir_includes_public_api():
    """`dir(neoswga)` should include the lazy exports for IDE autocomplete."""
    import neoswga
    visible = set(dir(neoswga))
    missing = [n for n in PUBLIC_API if n not in visible]
    assert not missing, (
        f"dir(neoswga) missing public symbols: {missing}"
    )


def test_lazy_import_actually_imports():
    """The lazy __getattr__ must hit the real module, not a stub. Use the
    run_optimization function as the canary: it's the highest-level API."""
    import neoswga
    run_optimization = neoswga.run_optimization
    # Confirm the resolved object is the real function from
    # unified_optimizer, not something else.
    import neoswga.core.unified_optimizer as _uo
    assert run_optimization is _uo.run_optimization


def test_get_typical_amplicon_length_reachable_from_top_level():
    """The new Phase 16 helper must be reachable from the top-level
    package so library users can use it without knowing submodule layout."""
    from neoswga import get_typical_amplicon_length
    assert get_typical_amplicon_length("phi29") == 3000
    assert get_typical_amplicon_length("equiphi29") == 4000


def test_attribute_error_for_nonexistent_symbol():
    """Typos should raise AttributeError with a clear message, not
    silently succeed."""
    import neoswga
    with pytest.raises(AttributeError, match="has no attribute 'does_not_exist'"):
        _ = neoswga.does_not_exist


def test_import_neoswga_is_cheap():
    """Importing the top-level package must not trigger the heavy
    optimizer import graph. Check by importing it in a clean subprocess
    and measuring that a handful of known-heavy modules are NOT loaded
    until a public symbol is requested."""
    import subprocess
    import sys

    script = (
        "import sys\n"
        "import neoswga\n"
        "heavy = [m for m in sys.modules if m.startswith('neoswga.core.optimizer_factory')]\n"
        "assert not heavy, f'optimizer_factory loaded at import time: {heavy}'\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True, text=True, timeout=30,
    )
    assert result.returncode == 0, (
        f"Import was not lazy: {result.stderr}"
    )


def test_lazy_imports_are_cached():
    """Second access to a lazy export must not re-import the module."""
    import neoswga
    first = neoswga.run_optimization
    second = neoswga.run_optimization
    assert first is second, "Lazy export should be cached after first access"
    # After first access the symbol should be a plain module-level attr.
    assert "run_optimization" in neoswga.__dict__
