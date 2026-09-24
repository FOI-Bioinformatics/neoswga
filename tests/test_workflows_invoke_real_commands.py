"""Every `neoswga` command a workflow runs must exist.

`neoswga score` was renamed to `prepare-candidates` on 2026-09-21,
deliberately with no alias, so that stale usages would surface rather than
keep working. One surfaced in this repository's own Nightly E2E workflow and
ran red for three nights before anyone read it.

A ratchet for that rename already existed. It opens with "No runnable example
may still use the old name" and then checks `README.md` and `CLAUDE.md` only.
A workflow is the one runnable example that runs itself, and it was out of
scope, so the check covered the copies and missed the original.

This is the general form rather than another rename-specific check: whatever
the argument list is, a command name a workflow invokes must be one the CLI
registers. It would have failed on the day `score` was retired.

Argparse is the authority for what exists, because a list written here would
be a second place to forget.
"""

import pathlib
import re

import pytest

yaml = pytest.importorskip("yaml")

_WORKFLOWS = pathlib.Path(__file__).resolve().parent.parent / ".github" / "workflows"

# `neoswga` where it STARTS a command, followed by the first bare word.
#
# Anchoring matters: `ruff check neoswga tests` passes the package as a path
# argument, and a bare \bneoswga\b would read `tests` as a subcommand. So the
# name must begin a line or follow a shell separator, optionally with a path
# prefix such as /tmp/wheel-venv/bin/. A leading `-` is an option, so
# `neoswga --version` carries no subcommand and is skipped.
_INVOCATION = re.compile(r"(?m)(?:^|[;&|]\s*)\s*(?:[\w./-]*/)?neoswga\s+(?!-)([a-z][a-z0-9-]*)")


def _registered_subcommands():
    from neoswga.cli_unified import create_parser

    for action in create_parser()._actions:
        if getattr(action, "choices", None) and isinstance(action.choices, dict):
            return set(action.choices)
    raise AssertionError("no subparser action found on the neoswga parser")


def _run_scripts(path):
    """Every shell script a workflow's steps run, with the step name."""
    document = yaml.safe_load(path.read_text())
    for job_name, job in (document.get("jobs") or {}).items():
        for step in job.get("steps") or []:
            script = step.get("run")
            if script:
                yield f"{path.name}:{job_name}:{step.get('name', '<unnamed>')}", script


def _workflow_files():
    return sorted(_WORKFLOWS.glob("*.yml")) + sorted(_WORKFLOWS.glob("*.yaml"))


def test_there_are_workflows_to_check():
    """A glob that matches nothing would make every assertion below vacuous."""
    assert _workflow_files(), f"no workflow files found under {_WORKFLOWS}"


def test_every_command_a_workflow_runs_is_registered():
    registered = _registered_subcommands()
    unknown = []

    for path in _workflow_files():
        for where, script in _run_scripts(path):
            for name in _INVOCATION.findall(script):
                if name not in registered:
                    unknown.append((where, name))

    assert not unknown, "a workflow invokes a neoswga command that does not exist: " + "; ".join(
        f"{where} runs 'neoswga {name}'" for where, name in unknown
    )


def test_the_check_can_see_a_retired_command():
    """Verified load-bearing: the pattern must match the shape that broke.

    `score` was retired, so it must not be registered, and the regex must find
    it in a line of the form the nightly workflow actually carried.
    """
    assert "score" not in _registered_subcommands()
    assert _INVOCATION.findall("          neoswga score -j params.json") == ["score"]


def test_options_are_not_mistaken_for_commands():
    assert _INVOCATION.findall("neoswga --version") == []
    assert _INVOCATION.findall("/tmp/venv/bin/neoswga --help > /dev/null") == []
    assert _INVOCATION.findall("neoswga show-presets > /dev/null") == ["show-presets"]


def test_the_package_as_a_path_argument_is_not_an_invocation():
    """`ruff check neoswga tests` is not a run of `neoswga tests`.

    The first version of this check reported exactly that, which is how the
    anchoring came to be here rather than a bare word boundary.
    """
    assert _INVOCATION.findall("          ruff check neoswga tests") == []
    assert _INVOCATION.findall("mypy neoswga --ignore-missing-imports") == []
    assert _INVOCATION.findall("          black --check neoswga/") == []


def test_a_command_after_a_shell_separator_is_still_seen():
    assert _INVOCATION.findall("cd examples && neoswga optimize -j p.json") == ["optimize"]
