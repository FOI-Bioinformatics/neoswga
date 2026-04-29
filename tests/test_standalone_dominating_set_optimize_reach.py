"""Regression test: dominating_set_optimizer.optimize() (standalone) must
construct DominatingSetOptimizer with a non-zero extension_reach.

This standalone function is the legacy CLI-integration entry. It is no
longer called from the unified CLI path, but external scripts and tests
may still import it. Without an extension_reach argument the constructor
defaults to 0, meaning each binding site covers exactly one 10 kb bin
and reported coverage understates amplicon-extended coverage.
"""

from pathlib import Path


def _read_function_source(name: str) -> str:
    src = (
        Path(__file__).resolve().parents[1]
        / "neoswga" / "core" / "dominating_set_optimizer.py"
    ).read_text()
    # Pull the chunk between "def {name}(" and the next top-level "def" or end.
    needle = f"def {name}("
    start = src.index(needle)
    rest = src[start + len(needle):]
    next_def = rest.find("\ndef ")
    end = start + len(needle) + (next_def if next_def != -1 else len(rest))
    return src[start:end]


def test_standalone_optimize_does_not_default_extension_reach_to_zero():
    body = _read_function_source("optimize")

    # The function must either pass extension_reach explicitly or import
    # polymerase_extension_reach to derive it from the polymerase.
    has_explicit_reach = "extension_reach=" in body
    has_helper = "polymerase_extension_reach" in body

    assert has_explicit_reach or has_helper, (
        "Standalone optimize() must set extension_reach when constructing "
        "DominatingSetOptimizer; otherwise it defaults to 0 and reports "
        "bin-only coverage with no amplicon credit."
    )


def test_standalone_optimize_uses_realistic_reach_helper():
    """Match the CLI dispatch convention: derive reach from
    coverage.polymerase_extension_reach with coverage_metric='realistic',
    not a hardcoded number that could drift from the rest of the codebase.
    """
    body = _read_function_source("optimize")
    assert "polymerase_extension_reach" in body, (
        "Use coverage.polymerase_extension_reach(polymerase, coverage_metric='realistic') "
        "to derive per-primer reach so this entry stays consistent with "
        "the unified_optimizer / dominating_set_adapter path."
    )
