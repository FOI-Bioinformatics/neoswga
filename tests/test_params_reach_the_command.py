"""`-j params.json` must actually load the params.

`pipeline._initialize()` does not take a path. It reads `parameter.json_file`,
a module global, and every CLI command is expected to set that from
`args.json_file` before calling it. Miss the assignment and `_initialize()`
still runs, still succeeds, and loads nothing -- so the command proceeds with
an empty configuration.

`evaluate-set` missed it. The symptom was a report of a missing `fg_prefixes`
on a params.json that plainly contained one, preceded by a wall of "Missing
parameter" warnings for keys that were also present. Nothing pointed at the
real cause, because the file was never opened: the command was reporting what
it found in an unpopulated global, accurately.

Two tests here. One drives `_resolve_sources` against a real file, which is
the specific bug. The other walks every CLI function that calls
`_initialize()` and checks it sets `json_file` first, which is the class --
the assignment is easy to leave out precisely because omitting it produces no
error at the call site.
"""

import ast
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

CLI_DIR = Path(__file__).resolve().parent.parent / "neoswga" / "cli"


@pytest.fixture
def params_file(tmp_path):
    """A params.json carrying foreground genome configuration."""
    genome = tmp_path / "target.fasta"
    genome.write_text(">t\n" + "ACGT" * 500 + "\n")

    path = tmp_path / "params.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 2,
                "data_dir": str(tmp_path),
                "fg_genomes": [str(genome)],
                "fg_prefixes": [str(tmp_path / "target")],
                "fg_seq_lengths": [2000],
                "polymerase": "phi29",
                "reaction_temp": 30.0,
            }
        )
    )
    return path


@pytest.fixture(autouse=True)
def restore_pipeline_state():
    """`_initialize` caches on module globals; leave them as we found them."""
    from neoswga.core import parameter
    from neoswga.core import pipeline as pipeline_mod

    previous = getattr(parameter, "json_file", None)
    yield
    parameter.json_file = previous
    pipeline_mod._initialized = False


# ----------------------------------------------------------------------
# The specific bug
# ----------------------------------------------------------------------


def test_evaluate_set_reads_the_params_file_it_was_given(params_file, tmp_path):
    """`-j params.json` has to reach `parameter.json_file`.

    Before the fix this raised "No foreground genome information available"
    against a file that lists fg_prefixes, fg_seq_lengths and fg_genomes.
    """
    from neoswga.cli.evaluate import _resolve_sources

    args = SimpleNamespace(
        json_file=str(params_file),
        genome=None,
        output=str(tmp_path / "out"),
        linear=False,
    )

    prefixes, genomes, lengths = _resolve_sources(args)

    assert prefixes == [str(tmp_path / "target")]
    assert lengths == [2000]
    assert genomes, "fg_genomes did not come through either"


def test_a_params_file_without_foreground_still_reports_clearly(tmp_path):
    """The error message is only useful once the file is genuinely read."""
    from neoswga.cli.evaluate import _resolve_sources

    empty = tmp_path / "empty.json"
    empty.write_text(json.dumps({"schema_version": 2, "data_dir": str(tmp_path)}))

    args = SimpleNamespace(
        json_file=str(empty), genome=None, output=str(tmp_path / "o"), linear=False
    )

    with pytest.raises(ValueError, match="--genome"):
        _resolve_sources(args)


# ----------------------------------------------------------------------
# The class of bug
# ----------------------------------------------------------------------


def _sets_json_file(node: ast.AST) -> bool:
    """Does this function set `parameter.json_file` before initializing?

    Three accepted spellings: the `@params_command` decorator, whose whole job
    is to validate the path and merge it onto `parameter` (this is the one
    `evaluate-set` was missing); the direct assignment used in
    `cli/pipeline.py`; and an explicit
    `merge_args_to_parameter(args, parameter, [..., "json_file", ...])`.
    """
    for decorator in getattr(node, "decorator_list", []):
        name = (
            getattr(decorator, "id", None)
            or getattr(decorator, "attr", None)
            or getattr(getattr(decorator, "func", None), "id", None)
            or getattr(getattr(decorator, "func", None), "attr", None)
        )
        if name == "params_command":
            return True

    for child in ast.walk(node):
        if isinstance(child, ast.Assign):
            for target in child.targets:
                if isinstance(target, ast.Attribute) and target.attr == "json_file":
                    return True
        if isinstance(child, ast.Call):
            name = getattr(child.func, "id", None) or getattr(child.func, "attr", None)
            if name == "merge_args_to_parameter":
                if "json_file" in ast.dump(child):
                    return True
    return False


def _calls_initialize(node: ast.AST) -> bool:
    for child in ast.walk(node):
        if isinstance(child, ast.Call):
            if getattr(child.func, "attr", None) == "_initialize":
                return True
    return False


def _cli_functions():
    for path in sorted(CLI_DIR.glob("*.py")):
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                yield path.name, node


CANDIDATES = [(f, n) for f, n in _cli_functions() if _calls_initialize(n)]


def test_there_are_functions_to_check():
    """Guards the guard: an AST walk that matches nothing passes vacuously."""
    assert CANDIDATES, "no CLI function calling _initialize() was found"


@pytest.mark.parametrize("filename,func", CANDIDATES, ids=[f"{f}::{n.name}" for f, n in CANDIDATES])
def test_initialize_callers_set_the_params_path(filename, func):
    """`_initialize()` reads a global, so forgetting to set it is silent.

    It does not raise, does not warn about the real problem, and loads an
    empty configuration that the command then reports on as though it were
    the user's.
    """
    assert _sets_json_file(func), (
        f"{filename}::{func.name} calls pipeline._initialize() without setting "
        f"parameter.json_file first, so -j would be ignored"
    )
