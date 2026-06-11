"""Module-size ratchet to stop the god-file from growing.

`cli_unified.py` is the project's largest module by a wide margin (it holds the
argparse surface and every CLI handler). This guard does not force an immediate
split, but it pins the current sizes as a CEILING: a module may shrink, never
grow past its budget. Lowering a budget here after a refactor is encouraged;
raising one should be a deliberate, reviewed decision.

To tighten the ratchet after extracting code, lower the matching budget.
"""

from pathlib import Path

import pytest

_PKG_ROOT = Path(__file__).resolve().parent.parent / "neoswga"

# Per-module ceilings (current size, rounded up slightly for trailing edits).
# A module NOT listed here must stay under DEFAULT_BUDGET.
_BUDGETS = {
    # Being split into neoswga/cli/ submodules; budget ratchets down per slice.
    "cli_unified.py": 2400,
    "core/report/technical_report.py": 2200,
}

# Any other single module should stay below this. The current second-largest
# non-listed module is reaction_conditions.py at ~1416 LOC.
DEFAULT_BUDGET = 1600


def _rel(path: Path) -> str:
    return path.relative_to(_PKG_ROOT).as_posix()


def _all_modules():
    return sorted(_PKG_ROOT.rglob("*.py"))


@pytest.mark.parametrize("path", _all_modules(), ids=_rel)
def test_module_within_size_budget(path):
    loc = sum(1 for _ in path.open(encoding="utf-8"))
    rel = _rel(path)
    budget = _BUDGETS.get(rel, DEFAULT_BUDGET)
    assert loc <= budget, (
        f"{rel} is {loc} LOC, over its {budget} budget. Split it into smaller "
        f"modules (do not raise the budget without review). See "
        f"tests/test_module_size_ratchet.py."
    )


def test_no_unlisted_module_exceeds_default():
    """Catch a brand-new oversized module that isn't in _BUDGETS yet."""
    offenders = []
    for path in _all_modules():
        rel = _rel(path)
        if rel in _BUDGETS:
            continue
        loc = sum(1 for _ in path.open(encoding="utf-8"))
        if loc > DEFAULT_BUDGET:
            offenders.append(f"{rel} ({loc} LOC)")
    assert not offenders, (
        "New oversized module(s) exceed the default budget "
        f"({DEFAULT_BUDGET} LOC): {offenders}. Split them or add a reviewed "
        "budget entry in tests/test_module_size_ratchet.py."
    )
