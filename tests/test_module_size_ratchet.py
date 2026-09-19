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
    # cli_unified is now just create_parser orchestration + main dispatch; the
    # per-group subparser builders live in neoswga/cli/<group>.py next to their
    # handlers. pipeline.py holds the count-kmers/filter/score/optimize/build
    # handlers + their argparse (the optimize command alone is ~250 arg lines).
    "cli_unified.py": 500,
    "cli/pipeline.py": 1800,
    # HTML template extracted to report/templates/technical_report.html.
    "core/report/technical_report.py": 1300,
    # Reviewed 2026-09-14. Both grew carrying the dimer-policy work that settled
    # `num_primers` as a request: the relaxation flag threaded through every
    # stage, the swap refinement behind `refinement_method`, and the reporting
    # that makes a short panel say why it is short. hybrid_optimizer gave back
    # 26 lines when the Stage-2 swap glue moved to core/swap_refinement.py;
    # these ceilings cover what is left, and splitting `HybridOptimizer.optimize`
    # into its three stages is still the right answer for the residue.
    "core/hybrid_optimizer.py": 1650,
    "core/unified_optimizer.py": 1650,
    # Reviewed 2026-09-19, +1 for `stage1_objective_width`. parameter.py sat at
    # exactly 1600, so ANY new params.json key costs a line: the key needs a
    # module global because `unified_optimizer.pick` resolves through
    # `getattr(parameter, name)`. One `global` declaration was folded into an
    # adjacent one to pay for half of it. This ceiling is a holding action --
    # the module is a flat list of per-key assignments and wants splitting by
    # concern (thermodynamics, dimers, search budgets), not another +1.
    "core/parameter.py": 1601,
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
