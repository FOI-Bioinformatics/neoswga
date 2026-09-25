"""Every read of a position index goes through `core/position_index.py`.

The sorted-blocks layout stores no dataset per k-mer, so a reader that asks
`primer in h5py_file` finds nothing and reports every primer as never scanned
-- or, where a missing key is read as zero sites, as binding nowhere. That is
the silent zero of Known Issues 5, 6, 13 and 15, and no behavioural test can
see a new raw reader until someone runs it against a converted index. So this
is a source check, in the manner of the one-door rule for alignment files.

What remains outside the layer is listed with its reason. The list can only
shrink.
"""

import ast
from pathlib import Path

PACKAGE = Path(__file__).resolve().parent.parent / "neoswga"

# file -> (number of h5py.File calls, the only modes they may use, reason)
ALLOWED = {
    "core/position_index.py": (None, None, "the layer itself"),
    "core/string_search.py": (
        4,
        {"a"},
        "creates an empty index before a scan, and writes or clears the "
        "provenance ROOT ATTRIBUTES, which both layouts store identically",
    ),
    "core/swga_simulator.py": (
        1,
        {"r"},
        "detects a legacy nested genome/primer layout, and hands the flat "
        "case to PositionIndex over the same handle",
    ),
}


def _h5py_file_calls(tree):
    for node in ast.walk(tree):
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "File"
            and isinstance(node.func.value, ast.Name)
            and node.func.value.id == "h5py"
        ):
            mode = None
            if len(node.args) > 1 and isinstance(node.args[1], ast.Constant):
                mode = node.args[1].value
            for keyword in node.keywords:
                if keyword.arg == "mode" and isinstance(keyword.value, ast.Constant):
                    mode = keyword.value.value
            yield node.lineno, mode


def test_no_module_opens_an_index_behind_the_layers_back():
    offenders = []
    counts = {}
    for path in sorted(PACKAGE.rglob("*.py")):
        rel = path.relative_to(PACKAGE).as_posix()
        calls = list(_h5py_file_calls(ast.parse(path.read_text())))
        if not calls:
            continue
        if rel not in ALLOWED:
            offenders.extend(f"{rel}:{line}" for line, _ in calls)
            continue
        expected, modes, _reason = ALLOWED[rel]
        counts[rel] = len(calls)
        if modes is not None:
            offenders.extend(
                f"{rel}:{line} mode={mode!r}" for line, mode in calls if mode not in modes
            )
        if expected is not None and len(calls) > expected:
            offenders.append(f"{rel}: {len(calls)} h5py.File calls, allowed {expected}")
    assert not offenders, (
        "Open position indexes through neoswga.core.position_index.open_index; "
        "a raw h5py reader cannot see the sorted-blocks layout: " + ", ".join(offenders)
    )


def test_the_allowlist_holds_nothing_stale():
    for rel, (expected, _modes, _reason) in ALLOWED.items():
        path = PACKAGE / rel
        assert path.exists(), rel
        if expected is None:
            continue
        found = len(list(_h5py_file_calls(ast.parse(path.read_text()))))
        assert found == expected, f"{rel} now has {found} h5py.File calls; lower its allowance"
