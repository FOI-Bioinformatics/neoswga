"""Every repair path reaches one implementation, so two cannot disagree.

CLAUDE.md states that `optimize` and `plan-pool` use "the same repair ... so
there is one in the codebase rather than two that can disagree". That was true
and unenforced: `optimization_service.repair_result` delegates to
`pool_planner.repair_panel`, and nothing stopped a second implementation
appearing beside it.

This repository has paid for exactly that twice. `string_search` carried two
position scanners that disagreed across a record join, and the delivered site set
depended on which ran. `panel_evaluation` measured complementary runs itself
while `dimer.dimer_validation_issue` did too, and the two agreed on every panel a
resolved request could produce, which is what made the duplication invisible
until they drifted.

So the check is structural rather than behavioural: a second repair would agree
with the first on most panels, and a behavioural test would pass right up to the
day it mattered.

Also recorded here: the late repair is a FINAL pass, not a duplicate of the
in-search repair stage. `run_panel_search` runs `repair` before `refinement` and
`reduction`, either of which can move the panel afterwards, so a violation can
appear after the in-search repair has already run. Measured on the plasmid
example with an unsatisfiable density floor: the late repair ran, changed nothing
and reported the violation rather than claiming compliance.
"""

import ast
import inspect
import pathlib

from neoswga.core import optimization_service, panel_acceptance, pool_planner

#: The one bounded repair. Everything else must reach it rather than reimplement it.
REPAIR = "repair_panel"


def _called_names(module):
    tree = ast.parse(inspect.getsource(module))
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            if isinstance(node.func, ast.Name):
                names.add(node.func.id)
            elif isinstance(node.func, ast.Attribute):
                names.add(node.func.attr)
    return names


def test_the_service_repair_delegates_rather_than_reimplementing():
    assert REPAIR in _called_names(optimization_service), (
        "repair_result must call pool_planner.repair_panel; a second bounded "
        "repair would agree with the first on most panels and diverge silently"
    )


def test_both_acceptance_paths_reach_the_one_repair():
    """`panel_acceptance` has two entry points and neither writes its own loop.

    `enforce_constraints` calls `repair_panel` directly and
    `apply_configured_limits` reaches it through `repair_result`. Calling the one
    implementation from two places is use, not duplication -- an earlier version
    of this test forbade the direct call and was simply wrong about that. What
    matters is that neither path contains a swap loop of its own.
    """
    names = _called_names(panel_acceptance)
    assert {
        "repair_panel",
        "repair_result",
    } & names, "neither acceptance path reaches the bounded repair"


def test_only_one_module_defines_a_bounded_repair():
    """The ratchet. A new `def repair_*` elsewhere is a second implementation.

    `repair_result` is the service's thin adapter and is named explicitly;
    anything else calling itself a repair has to justify itself here.
    """
    package = pathlib.Path("neoswga")
    allowed = {
        "core/pool_planner.py": {"repair_panel"},
        "core/optimization_service.py": {"repair_result"},
    }
    found = {}
    for path in sorted(package.rglob("*.py")):
        tree = ast.parse(path.read_text())
        defined = {
            node.name
            for node in ast.walk(tree)
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
            and node.name.startswith("repair")
        }
        if defined:
            found[path.relative_to(package).as_posix()] = defined

    unexpected = {rel: names - allowed.get(rel, set()) for rel, names in found.items()}
    unexpected = {rel: names for rel, names in unexpected.items() if names}
    assert not unexpected, (
        f"these define a repair outside the one implementation: {unexpected}. "
        f"Route them through pool_planner.repair_panel, or add them here with a "
        f"reason they are a different operation."
    )
    assert set(found) == set(allowed), (
        f"the allowlist names a module that no longer defines a repair: "
        f"{set(allowed) - set(found)}"
    )


def test_the_planner_repair_is_the_one_that_takes_the_objective():
    """It is the implementation, so it carries the parameters a repair needs."""
    signature = inspect.signature(pool_planner.repair_panel)
    for parameter in ("primers", "pool", "objective", "reasons", "config"):
        assert parameter in signature.parameters, (
            f"repair_panel lost {parameter!r}; callers pass the run's evaluator "
            f"and chemistry through it"
        )


def test_the_late_repair_runs_after_the_stages_that_can_break_a_panel():
    """Why the late pass is not a duplicate of the in-search repair stage.

    `run_panel_search` orders repair before refinement and reduction, and both
    of those move the panel, so a violation can appear after the in-search
    repair has already run. Pinned by reading the order the service records.
    """
    source = inspect.getsource(optimization_service)
    repair_at = source.index('record("repair"')
    refine_at = source.index('record("refinement"')
    assert repair_at < refine_at, (
        "repair now runs after refinement; if the order changed deliberately, "
        "the late repair's justification changes with it"
    )
