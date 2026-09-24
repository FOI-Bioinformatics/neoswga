"""A message that names a remedy must name one the user can actually use.

Found 2026-09-24 while auditing where error messages send people. The costly
instance was not mechanical -- `calibrate-reach` named a real flag that made
things worse, which no test can catch -- but two instances WERE mechanical and
this closes those:

- `wizard.py` told a first-time user "Consider using --cpus to parallelize".
  There is no `--cpus` on any command; `cpus` is a params.json key. Every
  command answers `unrecognized arguments: --cpus`. The setup wizard is the
  worst possible audience for that.
- `dimer.py` named "options -t or --max_dimer_bp" in three docstrings. The
  flag is `--max-dimer-bp` and there is no `-t`.

So this checks the three things that can be checked: a `--flag` named in a
string is declared on some parser, a `neoswga <subcommand>` named is a real
subcommand, and a `scripts/...py` named exists on disk.

It cannot check whether the advice is CORRECT, which is the failure that
actually cost something here. `test_the_advice_does_not_lead_into_the_defect.py`
covers that for the one case where it was measured. Read this as the floor, not
the ceiling.
"""

import ast
import pathlib
import re

PACKAGE = pathlib.Path(__file__).resolve().parent.parent / "neoswga"

FLAG = re.compile(r"--[a-zA-Z][a-zA-Z0-9-]{2,}")
COMMAND = re.compile(r"neoswga\s+([a-z][a-z0-9-]{2,})")
SCRIPT = re.compile(r"(scripts/[\w/]+\.py)")

# Flags belonging to OTHER tools, named deliberately in advice about them.
FOREIGN_FLAGS = {
    "--force-reinstall",  # pip
    "--help",  # argparse supplies it; never declared by hand
}

# Words that follow "neoswga " in ordinary prose. Each is a phrase, not a
# command: "a neoswga process", "two neoswga runs", "the neoswga pipeline".
PROSE_AFTER_NEOSWGA = {
    "command",
    "configuration",
    "features",
    "import",
    "pipeline",
    "primer",
    "process",
    "run",
    "versions",
}


def _string_constants():
    for path in sorted(PACKAGE.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and isinstance(node.value, str):
                yield path, node.lineno, node.value


def _declared_options():
    declared = set()
    for path in PACKAGE.rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "add_argument"
            ):
                for arg in node.args:
                    if isinstance(arg, ast.Constant) and str(arg.value).startswith("-"):
                        declared.add(arg.value)
    return declared


def test_every_flag_named_in_a_message_is_declared_somewhere():
    """`--cpus` was the live one: advice a first-time user cannot follow."""
    declared = _declared_options() | FOREIGN_FLAGS
    offenders = []
    for path, lineno, text in _string_constants():
        for flag in FLAG.findall(text):
            if flag in declared:
                continue
            # A double hyphen mid-prose ("min--max", "well--basically") is not
            # a flag. A real one starts the string or follows a space, quote or
            # bracket.
            if not re.search(rf"(^|[\s\'\"`(]){re.escape(flag)}\b", text):
                continue
            offenders.append(f"{path.relative_to(PACKAGE.parent)}:{lineno} names {flag}")

    assert not offenders, (
        "these name a flag no parser declares; name the params.json key or the "
        "real flag instead:\n  " + "\n  ".join(sorted(set(offenders)))
    )


def test_every_command_named_in_a_message_exists():
    from neoswga.cli_unified import create_parser

    real = set()
    for group in create_parser()._subparsers._group_actions:
        real |= set((group.choices or {}).keys())
    assert real, "no subcommands found; this test has gone stale"

    offenders = []
    for path, lineno, text in _string_constants():
        for command in COMMAND.findall(text):
            if command in real or command in PROSE_AFTER_NEOSWGA:
                continue
            offenders.append(
                f"{path.relative_to(PACKAGE.parent)}:{lineno} names 'neoswga {command}'"
            )

    assert not offenders, (
        "these name a subcommand that does not exist. `score` was renamed to "
        "`prepare-candidates` with no alias, so a stale name is reachable "
        "advice that fails:\n  " + "\n  ".join(sorted(set(offenders)))
    )


def test_every_script_named_in_a_message_exists():
    repo = PACKAGE.parent
    offenders = [
        f"{path.relative_to(repo)}:{lineno} names {script}"
        for path, lineno, text in _string_constants()
        for script in SCRIPT.findall(text)
        if not (repo / script).exists()
    ]

    assert not offenders, "these name a script that is not in the repository:\n  " + "\n  ".join(
        offenders
    )
