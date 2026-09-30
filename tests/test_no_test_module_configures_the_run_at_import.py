"""A test module must not assign to `neoswga.core.parameter` when imported.

pytest imports every test module during collection, before any test runs. A
module-level `parameter.max_self_dimer_bp = 3` therefore configures the whole
session, not the file it sits in, and what a later test measures depends on
which files were collected with it.

That is how `test_option_changes_the_design` came to fail for
`minimize_primers` and `target_coverage` in a full run and pass alone:
`tests/test_dimer.py` lowered the self-dimer limit from 4 to 3 at import, the
delivered panel held a self-dimer at 3, and every deletion from it counted as
a violation. Set such a value in a fixture with `monkeypatch.setattr`, which
restores it on teardown.
"""

import ast
from pathlib import Path

TESTS = Path(__file__).resolve().parent


def _import_time_parameter_assignments(source: str) -> list[int]:
    """Lines assigning `parameter.<name>` outside any function or class body."""
    found = []

    def visit(node):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda)):
            return
        if isinstance(node, (ast.Assign, ast.AugAssign, ast.AnnAssign)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            for target in targets:
                owner = getattr(target, "value", None)
                if isinstance(target, ast.Attribute) and (
                    (isinstance(owner, ast.Name) and owner.id == "parameter")
                    or (isinstance(owner, ast.Attribute) and owner.attr == "parameter")
                ):
                    found.append(node.lineno)
        for child in ast.iter_child_nodes(node):
            visit(child)

    visit(ast.parse(source))
    return found


def test_no_test_module_assigns_a_parameter_at_import():
    offenders = {}
    for path in sorted(TESTS.rglob("*.py")):
        lines = _import_time_parameter_assignments(path.read_text(encoding="utf-8"))
        if lines:
            offenders[path.relative_to(TESTS).as_posix()] = lines

    assert not offenders, (
        "these test modules assign to neoswga.core.parameter when imported, which "
        f"configures every test collected with them: {offenders}. Use a fixture "
        "with monkeypatch.setattr."
    )


def test_the_detector_sees_a_module_level_assignment_and_ignores_a_fixture():
    """Driven on source written here, so the detector cannot go blind quietly."""
    source = (
        "import neoswga.core.parameter as parameter\n"
        "parameter.max_dimer_bp = 3\n"
        "if True:\n"
        "    neoswga.core.parameter.max_self_dimer_bp = 3\n"
        "def fixture():\n"
        "    parameter.max_dimer_bp = 3\n"
    )

    assert _import_time_parameter_assignments(source) == [2, 4]
