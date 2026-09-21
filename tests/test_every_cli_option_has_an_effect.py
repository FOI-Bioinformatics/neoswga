"""A CLI option that nothing reads is a lie in the help text.

Phase 1 of the plan in `docs/validation/pipeline_audit_2026-09-16/`.

Known Issue 8 in CLAUDE.md records this class -- an option that is declared,
documented, accepted, and read by nothing -- and twice recorded it as closed.
It was not. The audit of 2026-09-16 found `--design-grid` on `plan-pool`,
added after the second closure, parsed and never read. Checking for more found
`--data-dir` and `--min-gini-sites` inert on the pipeline commands as well, the
second of which CLAUDE.md explicitly documents as working.

The three ratchets that were supposed to hold this line cannot see a CLI flag.
`tests/test_no_schema_key_is_inert.py` and
`tests/test_params_json_routes_optional_keys.py` iterate params.json schema
keys; `tests/test_cli_defaults_do_not_beat_params_json.py` checks four argparse
defaults and asserts nothing about whether a flag is read;
`tests/test_design_options_have_effect.py` calls `run_optimization` directly and
so covers the optimize path only.

This walks from each command's dispatch entry through the functions it calls,
collecting every attribute read off the argparse namespace, and compares that
against the options the parser declares. The dispatch table is read from
`main()` rather than restated here, so a command cannot be added without also
being checked.

An option that genuinely cannot be read yet needs an entry in `KNOWN_INERT`
naming the reason. The list can only shrink: `test_the_known_inert_list_has_no_
stale_entries` fails on an entry that is no longer needed, so a fix is not
finished until it also removes its excuse.
"""

import ast
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
PACKAGE = ROOT / "neoswga"

# How far to follow the argparse namespace through call chains. Six is enough
# for every current handler; raising it to twelve changed nothing.
MAX_DEPTH = 8

# Options that no code path reads today, with the reason. Shrink this list;
# do not grow it. A new option belongs in a handler, not here.
#
# "EmptyOptions" below is `core/pipeline.py:_initialize`, which builds an object
# returning None for every attribute except `json_file` and hands THAT to
# `parameter.get_params`. So a resolver reading `args.<name>` there receives
# None on every pipeline run, and the params.json value always wins. Verified
# for `data_dir` and `min_gini_sites` by resolving both with a flag value that
# differs from the configured one; the flag lost in each case.
KNOWN_INERT = {
    "count-kmers": {
        "data_dir": "EmptyOptions; params.json wins",
        "fasta_back": "EmptyOptions; params.json wins",
        "fasta_fore": "EmptyOptions; params.json wins",
        "kmer_back": "EmptyOptions; params.json wins",
        "kmer_fore": "EmptyOptions; params.json wins",
        "polymerase": "EmptyOptions; params.json wins",
    },
    "filter": {
        "data_dir": "EmptyOptions; params.json wins (verified)",
        "min_gini_sites": "EmptyOptions; params.json wins (verified). "
        "CLAUDE.md documents this flag as working; it does not",
    },
    "prepare-candidates": {
        "data_dir": "EmptyOptions; params.json wins",
        "enhanced_model_path": "no reader on the candidate-preparation path",
        "polymerase": "EmptyOptions; params.json wins",
        "use_enhanced_features": "no reader on the candidate-preparation path",
    },
    "optimize": {
        "background_bloom_path": "no reader on the optimize path",
        "background_sampled_path": "no reader on the optimize path",
        "data_dir": "EmptyOptions; params.json wins",
        "max_optimization_time": "no reader on the optimize path",
        "polymerase": "EmptyOptions; params.json wins",
    },
    "design": {"data_dir": "EmptyOptions; params.json wins"},
    "evaluate-set": {
        "data_dir": "EmptyOptions; params.json wins",
        "enable_qa": "no reader on the evaluate-set path",
        "gpu_device": "GPU helpers are not reached by any stage; see CLAUDE.md",
        "no_gpu": "GPU helpers are not reached by any stage; see CLAUDE.md",
        "reaction_temp": "no reader on the evaluate-set path",
        "use_gpu": "GPU helpers are not reached by any stage; see CLAUDE.md",
    },
    "rescore-set": {
        "betaine_m": "additive flags are not merged on this path",
        "dmso_percent": "additive flags are not merged on this path",
        "ethanol_percent": "additive flags are not merged on this path",
        "formamide_percent": "additive flags are not merged on this path",
        "json_file": "no reader; the handler uses the already-initialised globals",
        "mg_conc": "additive flags are not merged on this path",
        "na_conc": "additive flags are not merged on this path",
        "tmac_m": "additive flags are not merged on this path",
        "trehalose_m": "additive flags are not merged on this path",
        "urea_m": "additive flags are not merged on this path",
    },
    "contract-set": {"json_file": "no reader; uses the already-initialised globals"},
    "swap-primer": {"json_file": "no reader; uses the already-initialised globals"},
    "analyze-dimers": {"visualize": "no reader on the analyze-dimers path"},
    "analyze-set": {"simulate": "no reader on the analyze-set path"},
    "analyze-stability": {"temp": "no reader on the analyze-stability path"},
    "schema": {"dump": "no reader on the schema path"},
}


def _index_functions():
    """Every function definition in the package, by name.

    By name rather than by qualified path because the walk follows call sites,
    and a call site names the function rather than its module. Two functions
    sharing a name are both searched, which can only widen what counts as a
    read -- the conservative direction for a test that fails on unread options.
    """
    functions = {}
    for path in sorted(PACKAGE.rglob("*.py")):
        try:
            tree = ast.parse(path.read_text())
        except SyntaxError:  # pragma: no cover - would fail the suite elsewhere
            continue
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                functions.setdefault(node.name, []).append(node)
    return functions


FUNCTIONS = _index_functions()


def _dispatch_table():
    """Command name -> handler function name, read from `main()`.

    Read rather than restated, so a command added to the dispatch dict is
    checked without anyone remembering to add it here. A lambda maps to None
    and is treated as reading nothing.
    """
    tree = ast.parse((PACKAGE / "cli_unified.py").read_text())
    main = next(
        node for node in ast.walk(tree) if isinstance(node, ast.FunctionDef) and node.name == "main"
    )
    table, read_in_main = {}, set()
    for node in ast.walk(main):
        if isinstance(node, ast.Assign) and isinstance(node.value, ast.Dict):
            for key, value in zip(node.value.keys, node.value.values):
                if isinstance(key, ast.Constant):
                    table[key.value] = value.id if isinstance(value, ast.Name) else None
        read_in_main |= _namespace_reads(node, {"args"})
    return table, read_in_main


def _namespace_reads(node, names):
    """Attribute names read off one of `names` by this single node."""
    found = set()
    if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name):
        if node.value.id in names:
            found.add(node.attr)
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
        if (
            node.func.id in ("getattr", "hasattr")
            and len(node.args) >= 2
            and isinstance(node.args[0], ast.Name)
            and node.args[0].id in names
            and isinstance(node.args[1], ast.Constant)
        ):
            found.add(node.args[1].value)
    return found


def _declared_merges(function):
    """Options a handler merges into the `parameter` module globals.

    `@params_command(merge=[...])` and `merge_args_to_parameter(args, parameter,
    [...])` both read the namespace without ever writing `args.<name>`, so a
    walk that looked only for attribute access would call every chemistry flag
    on `filter` inert. It is a read, by a different route.
    """
    merged = set()
    for decorator in function.decorator_list:
        if isinstance(decorator, ast.Name) and decorator.id == "params_command":
            merged.add("json_file")
        if isinstance(decorator, ast.Call):
            for keyword in decorator.keywords:
                if keyword.arg == "merge" and isinstance(keyword.value, (ast.List, ast.Tuple)):
                    merged |= {
                        element.value
                        for element in keyword.value.elts
                        if isinstance(element, ast.Constant)
                    }
                if keyword.arg == "seed" and isinstance(keyword.value, ast.Constant):
                    if keyword.value.value:
                        merged.add("seed")
    for node in ast.walk(function):
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "merge_args_to_parameter"
            and len(node.args) >= 3
            and isinstance(node.args[2], (ast.List, ast.Tuple))
        ):
            merged |= {
                element.value for element in node.args[2].elts if isinstance(element, ast.Constant)
            }
    return merged


def _forwarded_calls(node, local):
    """(callee name, parameter index) for each call handed one of `local`."""
    if not isinstance(node, ast.Call):
        return
    name = node.func.id if isinstance(node.func, ast.Name) else None
    if name is None and isinstance(node.func, ast.Attribute):
        name = node.func.attr
    if name is None:
        return
    for index, argument in enumerate(node.args):
        if isinstance(argument, ast.Name) and argument.id in local:
            yield name, index, None
    for keyword in node.keywords:
        if isinstance(keyword.value, ast.Name) and keyword.value.id in local:
            yield name, None, keyword.arg


def options_read_by(handler_name, bound=("args",), depth=0, seen=None):
    """Every namespace attribute `handler_name` reads, following the namespace.

    Returns `(names, reads_everything)`. `reads_everything` is set by a
    `vars(args)` or `args.__dict__`, where naming individual options is
    impossible and the honest answer is that all of them are reachable.
    """
    if seen is None:
        seen = set()
    key = (handler_name, tuple(sorted(bound)))
    if depth > MAX_DEPTH or key in seen or handler_name not in FUNCTIONS:
        return set(), False
    seen.add(key)

    found, everything = set(), False
    for function in FUNCTIONS[handler_name]:
        parameters = [a.arg for a in function.args.args] + [a.arg for a in function.args.kwonlyargs]
        local = {p for p in parameters if p in bound}
        if not local:
            continue
        found |= _declared_merges(function)
        for node in ast.walk(function):
            found |= _namespace_reads(node, local)
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Name)
                and node.func.id == "vars"
                and node.args
                and isinstance(node.args[0], ast.Name)
                and node.args[0].id in local
            ):
                everything = True
            for callee, index, keyword in _forwarded_calls(node, local):
                for target in FUNCTIONS.get(callee, []):
                    names = [a.arg for a in target.args.args]
                    if keyword is not None:
                        inner = {keyword}
                    else:
                        position = index + 1 if names and names[0] == "self" else index
                        if position >= len(names):
                            continue
                        inner = {names[position]}
                    deeper, deeper_all = options_read_by(callee, inner, depth + 1, seen)
                    found |= deeper
                    everything = everything or deeper_all
    return found, everything


def unread_options():
    """Command -> the options it declares and never reads."""
    from neoswga.cli_unified import create_parser

    table, read_in_main = _dispatch_table()
    parser = create_parser()
    offenders = {}
    for group in parser._subparsers._group_actions:
        for command, subparser in (group.choices or {}).items():
            handler = table.get(command)
            if handler is None:
                continue
            declared = {
                action.dest
                for action in subparser._actions
                if action.dest not in ("help", "command")
            }
            found, everything = options_read_by(handler)
            if everything:
                continue
            missing = declared - found - read_in_main
            if missing:
                offenders[command] = missing
    return offenders


def test_no_command_declares_an_option_it_never_reads():
    """The ratchet. A new unread option fails here."""
    offenders = unread_options()
    surprises = {}
    for command, missing in offenders.items():
        unexcused = missing - set(KNOWN_INERT.get(command, {}))
        if unexcused:
            surprises[command] = sorted(unexcused)

    assert not surprises, (
        "These CLI options are declared and never read on their command's path:\n"
        + "\n".join(
            f"  {command}: {', '.join(names)}" for command, names in sorted(surprises.items())
        )
        + "\n\nWire the option to something, or add it to KNOWN_INERT with the reason."
    )


def test_the_known_inert_list_has_no_stale_entries():
    """What makes this a ratchet rather than a suppression list.

    An entry that is no longer needed fails, so wiring an option is not finished
    until its excuse is removed. Without this the list would outlive the defects
    and quietly re-admit them.
    """
    offenders = unread_options()
    stale = {}
    for command, excused in KNOWN_INERT.items():
        still_unread = offenders.get(command, set())
        gone = sorted(set(excused) - still_unread)
        if gone:
            stale[command] = gone

    assert (
        not stale
    ), "These KNOWN_INERT entries name options that are now read. Remove them:\n" + "\n".join(
        f"  {command}: {', '.join(names)}" for command, names in sorted(stale.items())
    )


def test_the_walker_finds_a_read_that_goes_through_the_merge_decorator():
    """Guard the guard, part one.

    `filter` declares more than twenty chemistry options and reads none of them
    as `args.<name>`; they reach the code through `merge_args_to_parameter`. A
    walker blind to that route would report them all as inert, the list would be
    too long to act on, and the ratchet would be switched off.
    """
    found, _ = options_read_by("run_step2")

    for option in ("dmso_percent", "betaine_m", "na_conc", "max_gini", "max_primer"):
        assert option in found, f"{option} is merged by run_step2 but the walker missed it"


def test_the_walker_follows_the_namespace_into_a_helper():
    """Guard the guard, part two.

    `run_plan_pool` reads most of its options directly, but a walker that
    stopped at the handler would miss anything resolved in a helper it hands
    `args` to. If this ever passes trivially the interprocedural step has been
    lost and every offender list becomes too long to trust.
    """
    found, _ = options_read_by("run_plan_pool")

    assert "coverage_metric" in found
    assert "min_selectivity_density" in found
    assert "no_repair" in found


def test_the_design_grid_stays_wired():
    """The finding this ratchet was written for, now the other way round.

    `--design-grid` was audit finding F4: parsed, in the help text, and read by
    nothing. Phase 4 increment 6 wired it to `design_sweep`, so this asserts it
    stays read rather than that it stays unread. Pinned by name because a flag
    that silently stops being read is the whole class this file exists for.
    """
    offenders = unread_options()

    assert "design_grid" not in offenders.get("plan-pool", set()), (
        "plan-pool --design-grid is no longer read. It was wired to design_sweep "
        "in Phase 4 increment 6; if that was deliberately reverted, restore its "
        "KNOWN_INERT entry with the reason."
    )


@pytest.mark.parametrize("command", sorted(KNOWN_INERT))
def test_every_excused_command_still_exists(command):
    """An excuse for a command that has gone is dead weight."""
    from neoswga.cli_unified import create_parser

    names = {
        name
        for group in create_parser()._subparsers._group_actions
        for name in (group.choices or {})
    }

    assert command in names, f"KNOWN_INERT names {command!r}, which is not a command"
