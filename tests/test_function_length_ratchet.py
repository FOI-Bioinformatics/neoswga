"""Function-length ratchet — stops new giant functions landing silently.

Mirrors test_module_size_ratchet at function scope. A new (or unlisted) function
must stay under DEFAULT_BUDGET source lines; the handful of existing oversized
functions are pinned at their current size in _BUDGETS so they may shrink but
never grow. Lowering an entry after a refactor is encouraged; raising one (or
adding a new entry) should be a deliberate, reviewed decision.

Length is the AST span (def line .. last body line), which ignores trailing
module-level constants — so a big inline template after a function is not blamed
on it.
"""

import ast
from pathlib import Path

import pytest

_PKG_ROOT = Path(__file__).resolve().parent.parent / "neoswga"

# Any function NOT listed here must stay under this many source lines.
DEFAULT_BUDGET = 200

# Pinned ceilings for the existing oversized functions ("relpath::qualname").
# These are the only functions currently over DEFAULT_BUDGET; each is a known
# refactor candidate (e.g. create_parser -> per-group builders).
_BUDGETS = {
    "cli_unified.py::create_parser": 1460,
    "cli/pipeline.py::run_step4": 550,
    "core/parameter.py::get_params": 540,
    "core/unified_optimizer.py::run_optimization": 470,
    "core/hybrid_optimizer.py::optimize": 350,
    "core/report/technical_report.py::render_technical_report": 325,
    "core/workflow_selector.py::run_workflow_selector": 260,
    "cli/iterate.py::run_expand_primers": 245,
    "cli/pipeline.py::run_step2": 215,
    "core/additive_interactions.py::load_defaults": 215,
    "core/pipeline.py::step2": 210,
}


def _function_spans():
    """Yield (relpath, qualname, line_span) for every def in the package."""
    for path in sorted(_PKG_ROOT.rglob("*.py")):
        rel = path.relative_to(_PKG_ROOT).as_posix()
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"))
        except SyntaxError:  # pragma: no cover - defensive
            continue
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                span = (node.end_lineno or node.lineno) - node.lineno + 1
                yield rel, node.name, span


def test_no_function_exceeds_its_budget():
    offenders = []
    for rel, name, span in _function_spans():
        budget = _BUDGETS.get(f"{rel}::{name}", DEFAULT_BUDGET)
        if span > budget:
            offenders.append(f"{rel}::{name} = {span} lines (budget {budget})")
    assert not offenders, (
        "Function(s) over budget — split them (do not raise the budget without "
        f"review). See tests/test_function_length_ratchet.py:\n  " + "\n  ".join(offenders)
    )


def test_budget_entries_are_not_stale():
    """Every pinned entry must still exist and still be over DEFAULT_BUDGET.

    Keeps the allowlist honest: once a giant is refactored below the default it
    should be removed from _BUDGETS so it can't silently regrow.
    """
    # A file may hold two functions with the same name (e.g. an optimize() on
    # two classes); keep the largest span per key so a small namesake does not
    # mask the pinned giant.
    spans = {}
    for rel, name, span in _function_spans():
        key = f"{rel}::{name}"
        spans[key] = max(spans.get(key, 0), span)
    stale = [key for key in _BUDGETS if key not in spans or spans[key] <= DEFAULT_BUDGET]
    assert not stale, f"Remove stale/under-budget _BUDGETS entries: {stale}"
