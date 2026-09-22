"""A size ceiling with no headroom fails on the merge, not on the change.

`main` went red on 2026-09-22 with `base_optimizer.py` at 1603 lines against a
1600 budget, and no single pull request put it there. Three did: each added a
few lines, each measured 1600 exactly on its own branch, each passed CI, and
the sum arrived only once they were merged.

That is a property of the ceiling rather than of any change. A module sitting
ON its budget has no room for a correct small edit, so the failure surfaces
where it is hardest to attribute -- on a merge commit, in someone else's pull
request, naming a module they did not touch.

The remedy applied was the one the size ratchet's own docstring asks for:
extract, do not raise. `OptimizationResult.validate` moved to
`core/result_validation.py`, which is its own subject anyway.

This test is the early warning. It does not forbid a module from approaching
its budget -- that is what the budget is for -- but a module with almost no
room left is a merge conflict waiting to happen, and saying so before it
happens is cheaper than attributing it afterwards.
"""

import pathlib

import pytest

from tests.test_module_size_ratchet import _BUDGETS, DEFAULT_BUDGET, _PKG_ROOT

#: How many lines a module must be able to gain without breaking its budget.
#: Three, because three is what it took: the smallest useful edits here are a
#: deferred import plus a two-line comment saying why.
REQUIRED_HEADROOM = 4

#: Modules allowed to sit close to their ceiling, with the reason. This list
#: can only shrink; an entry is a promise to extract, not a permanent excuse.
ACCEPTED: dict = {}


def _budget(relative):
    return _BUDGETS.get(relative, DEFAULT_BUDGET)


def _measured():
    sizes = {}
    for path in _PKG_ROOT.rglob("*.py"):
        relative = path.relative_to(_PKG_ROOT).as_posix()
        sizes[relative] = len(path.read_text(encoding="utf-8").splitlines())
    return sizes


def test_no_module_sits_against_its_ceiling():
    tight = {
        relative: (size, _budget(relative))
        for relative, size in _measured().items()
        if relative not in ACCEPTED and _budget(relative) - size < REQUIRED_HEADROOM
    }

    assert not tight, (
        "these modules have fewer than "
        f"{REQUIRED_HEADROOM} lines of headroom, so the next correct small edit "
        "fails -- or, if two land separately, the merge fails and names a module "
        "nobody touched:\n  "
        + "\n  ".join(f"{p}: {size} of {budget}" for p, (size, budget) in tight.items())
        + "\n\nExtract; do not raise the budget."
    )


@pytest.mark.parametrize("relative", sorted(ACCEPTED))
def test_the_accepted_list_has_no_stale_entries(relative):
    """An entry is a promise to extract. Once kept, it goes."""
    size = _measured().get(relative)

    assert size is not None, f"{relative} no longer exists; drop its entry"
    assert (
        _budget(relative) - size < REQUIRED_HEADROOM
    ), f"{relative} now has room, so its entry in ACCEPTED can go"


def test_the_validator_left_the_optimizer():
    """The extraction that fixed it, pinned so it is not folded back in.

    Driven through a real result rather than by reading source text: a check
    that `base_optimizer` no longer contains the code would pass on a copy that
    computes the wrong thing.
    """
    from neoswga.core.base_optimizer import (
        OptimizationResult,
        OptimizationStatus,
        PrimerSetMetrics,
    )
    from neoswga.core.result_validation import validate_result

    result = OptimizationResult(
        primers=("ACGTACGTACGT", "ACGTACGTACGT"),
        score=1.0,
        status=OptimizationStatus.SUCCESS,
        metrics=PrimerSetMetrics.empty(),
        iterations=1,
        optimizer_name="test",
    )

    assert result.validate(target_size=2) == validate_result(result, 2)
    assert result.validate(target_size=2)["ok"] is False, "a duplicate is an error"


def test_the_extracted_module_is_not_where_the_lines_went():
    """Extraction must move a subject, not hide a module behind a re-export."""
    moved = pathlib.Path(_PKG_ROOT / "core" / "result_validation.py")

    assert moved.exists()
    assert len(moved.read_text().splitlines()) < 400, "this is becoming a second god-file"
