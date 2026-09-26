"""A module nothing calls is not a feature.

Phase 1 of the plan in `docs/validation/pipeline_audit_2026-09-16/`, the second
half. `tests/test_every_cli_option_has_an_effect.py` catches an option that is
declared and unread; this catches the same defect one layer down, where the
capability is built and tested and no command can reach it.

The audit found six at once, all from the condition-aware pool design work, all
with passing tests:

    CandidateProvider          built only by design_sweep, which no command calls
    design_sweep               no production caller
    load_design_grid           no production caller
    load_grid_file             no production caller
    ensure_positions           called nowhere, and a no-op if it were
                               (wired in Phase 4 increment 3; see
                               tests/test_the_frontier_is_vouched_for_before_it_is_scored.py)
    beam_search                reachable in principle, unreachable in practice

Unit tests cannot see this. Each of those has tests that construct the thing
directly and assert it behaves, which is what a caller would do, so they pass
whether or not anything calls it.

The check is transitive reach from the command dispatch table, not "referenced
somewhere". The difference is the whole point: `design_sweep` calls
`provider.expand`, so a check that accepted any reference would call `expand`
reachable on the strength of a caller that is itself unreachable. Reachability
has to start where a user starts.

Scope is deliberately the modules that defect came from rather than the whole
package: a package-wide version needs an allowlist longer than the check, and a
check nobody reads is the same failure in a different place. Add a module here
when it grows a capability a command is supposed to reach.
"""

import ast
import pathlib
import sys

import pytest

target_module = sys.modules[__name__]

ROOT = pathlib.Path(__file__).resolve().parent.parent
PACKAGE = ROOT / "neoswga"

# Modules whose public surface is supposed to be reachable from a command.
WATCHED = (
    "core/candidate_source.py",
    "core/candidate_provider.py",
    "core/pool_design_sweep.py",
    "core/panel_beam.py",
    "core/pool_objective.py",
    "core/pool_metrics.py",
    "core/partial_panel.py",
    "core/lazy_dimer.py",
    "core/panel_evaluation.py",
    "core/query_scan.py",
)

# Public names no command can reach today, and why. Shrink this list.
# `test_the_unreachable_list_has_no_stale_entries` fails on an entry that has
# since become reachable, so wiring something is not finished until its excuse
# is deleted.
# Empty, and meant to stay that way. Every capability the 2026-09-16 audit found
# unreachable now has a command behind it: `design_sweep` and `load_design_grid`
# through `plan-pool --design-grid` in Phase 4 increment 6, `CandidateProvider`
# and `ensure_positions` in increments 1 and 3, `beam_search` in Phase 2.
KNOWN_UNREACHABLE: dict[str, str] = {}


def _parse_all():
    """Every function in the package, by name, plus the module each came from."""
    functions, modules = {}, {}
    for path in sorted(PACKAGE.rglob("*.py")):
        try:
            tree = ast.parse(path.read_text())
        except SyntaxError:  # pragma: no cover
            continue
        modules[path] = tree
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                functions.setdefault(node.name, []).append(node)
    return functions, modules


FUNCTIONS, MODULES = _parse_all()


def _names_used_in(node):
    """Names a node references, as a call, an attribute, or a bare name.

    All three count. A helper is reached as `thing(...)`, a method as
    `obj.thing(...)`, and a class as a bare `CandidateProvider` handed to
    something else. Counting all three is the conservative direction for a test
    that fails on unreachable names.
    """
    used = set()
    for child in ast.walk(node):
        if isinstance(child, ast.Name):
            used.add(child.id)
        elif isinstance(child, ast.Attribute):
            used.add(child.attr)
        elif isinstance(child, ast.Call) and isinstance(child.func, ast.Name):
            # A dynamic reference is still a reference. `plan_pool` reaches
            # `compute_pool_metrics` as getattr(optimizer, "compute_pool_metrics",
            # None), so a walk that read only attribute access would call the
            # focused evaluator unreachable while it runs on every design.
            if child.func.id in ("getattr", "hasattr", "setattr"):
                if len(child.args) >= 2 and isinstance(child.args[1], ast.Constant):
                    if isinstance(child.args[1].value, str):
                        used.add(child.args[1].value)
    return used


def _dispatch_handlers():
    """The handler names in `main()`'s dispatch dict: where a user starts."""
    tree = MODULES[PACKAGE / "cli_unified.py"]
    main = next(
        node for node in ast.walk(tree) if isinstance(node, ast.FunctionDef) and node.name == "main"
    )
    handlers = set()
    for node in ast.walk(main):
        if isinstance(node, ast.Assign) and isinstance(node.value, ast.Dict):
            for key, value in zip(node.value.keys, node.value.values):
                if isinstance(key, ast.Constant) and isinstance(value, ast.Name):
                    handlers.add(value.id)
    return handlers


def reachable_names():
    """Every name transitively reachable from a command handler.

    Module-level code counts as reachable too: importing a module runs it, and a
    registration performed at import is a real edge. `unified_optimizer`'s
    optimizer registry works that way.
    """
    frontier = set(_dispatch_handlers())
    assert frontier, "no dispatch handlers found; the dispatch table moved"

    for path, tree in MODULES.items():
        for node in tree.body:
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                frontier |= _names_used_in(node)

    reached, pending = set(), list(frontier)
    while pending:
        name = pending.pop()
        if name in reached:
            continue
        reached.add(name)
        for function in FUNCTIONS.get(name, []):
            for used in _names_used_in(function):
                if used not in reached:
                    pending.append(used)
    return reached


def _public_names(path):
    """Public module-level functions and classes in one module.

    Methods are deliberately not tracked. Reachability here is by name, and a
    method name is rarely distinctive enough to carry one: `expand` is both
    `CandidateProvider.expand`, which no command reaches, and
    `PrimerExpander.expand`, which `expand-primers` reaches every run, and
    `initial` and `frontier` occur as ordinary local variables. Nothing short of
    type inference separates those, and a check that guesses is worse than one
    with a stated limit.

    Little is lost, because an unreachable class takes its methods with it.
    `CandidateProvider` is tracked, so `ensure_positions` and the batch methods
    are covered by the finding that matters, one level up.
    """
    names = {}
    for node in MODULES[path].body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            if not node.name.startswith("_"):
                names[node.name] = path
    return names


def unreachable_names():
    """Public names in the watched modules that no command can reach."""
    reached = reachable_names()
    orphans = {}
    for relative in WATCHED:
        path = PACKAGE / relative
        assert path.is_file(), f"{relative} is watched but does not exist"
        for name in _public_names(path):
            if name not in reached:
                orphans[name] = relative
    return orphans


def test_no_watched_capability_is_unreachable():
    """The ratchet. A new capability no command can reach fails here."""
    orphans = unreachable_names()
    surprises = {n: m for n, m in orphans.items() if n not in KNOWN_UNREACHABLE}

    assert not surprises, (
        "No command can reach these public names:\n"
        + "\n".join(
            f"  {module}: {name}"
            for name, module in sorted(surprises.items(), key=lambda kv: (kv[1], kv[0]))
        )
        + "\n\nConnect it to a command path, or add it to KNOWN_UNREACHABLE with the reason."
    )


def test_the_unreachable_list_has_no_stale_entries():
    """What makes this a ratchet rather than a suppression list."""
    orphans = unreachable_names()
    stale = sorted(set(KNOWN_UNREACHABLE) - set(orphans))

    assert (
        not stale
    ), "These KNOWN_UNREACHABLE entries are now reachable. Remove them:\n" + "\n".join(
        f"  {name}" for name in stale
    )


def test_the_capabilities_this_ratchet_was_written_for_stay_reachable():
    """The six findings this file was created for, asserted as fixed.

    Each was built, tested, merged and callable by nothing. Pinned by name so
    that losing a caller fails loudly rather than quietly restoring the defect.
    """
    orphans = unreachable_names()

    assert "design_sweep" not in orphans, (
        "design_sweep has lost its caller. It was wired to "
        "`plan-pool --design-grid` in Phase 4 increment 6; if that was "
        "deliberately reverted, restore its KNOWN_UNREACHABLE entry."
    )
    assert "CandidateProvider" not in orphans, (
        "CandidateProvider was wired in Phase 4 increment 1, through "
        "candidate_source.open_candidate_source. If it has gone out of reach "
        "again, the inventory is write-only once more."
    )


def test_reachability_does_not_flow_through_an_unreachable_caller(monkeypatch):
    """The property that made the first version of this test wrong.

    A check that counted any reference anywhere would call a name reachable on
    the strength of a caller no command can reach, and every finding would read
    as fine. `load_design_grid` used to be the live example: `cli/plan_pool.py`
    imported it inside `load_grid_file`, which nothing called. Phase 4 increment
    6 wired that path, so the example is gone and the property is asserted on a
    synthetic graph instead. A test that depends on a particular bug still
    existing stops testing anything the moment the bug is fixed.
    """
    import ast as _ast
    import textwrap as _textwrap

    module = _ast.parse(_textwrap.dedent("""
            def run_thing(args):
                reached_helper()

            def reached_helper():
                pass

            def nobody_calls_this():
                only_referenced_here()

            def only_referenced_here():
                pass
            """))
    functions = {}
    for node in _ast.walk(module):
        if isinstance(node, (_ast.FunctionDef, _ast.AsyncFunctionDef)):
            functions.setdefault(node.name, []).append(node)

    monkeypatch.setattr(target_module, "FUNCTIONS", functions)
    monkeypatch.setattr(target_module, "MODULES", {})
    monkeypatch.setattr(target_module, "_dispatch_handlers", lambda: {"run_thing"})

    reached = target_module.reachable_names()

    assert "reached_helper" in reached, "the walk did not follow a real call edge"
    assert "only_referenced_here" not in reached, (
        "a name referenced only inside an uncalled function was reported "
        "reachable, so the walk has regressed to counting bare references"
    )


def test_a_reached_capability_is_not_reported():
    """Guard the guard.

    `PoolObjective` and `compute_pool_metrics` are reached from `plan-pool`. If
    the walk ever breaks, everything would look unreachable, the list would be
    too long to act on, and the ratchet would be switched off rather than fixed.
    """
    orphans = unreachable_names()

    for reached in ("PoolObjective", "PoolConstraints", "compute_pool_metrics", "can_prune"):
        assert reached not in orphans, f"{reached} is reachable but was reported otherwise"


@pytest.mark.parametrize("relative", WATCHED)
def test_every_watched_module_is_scanned(relative):
    """A watched module that has been deleted or renamed must not pass silently."""
    assert (PACKAGE / relative).is_file()
