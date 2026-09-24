"""Guard: every test path referenced in a CI workflow must exist.

The nightly workflow once invoked `pytest tests/test_e2e_integration.py`, a file
that never existed, so the job exited with pytest's usage error (code 4) every
night. This test scans `.github/workflows/*.yml` for pytest path arguments and
asserts each referenced file/directory is present, catching that typo class at
PR time instead of at 04:00 UTC.
"""

import re
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parent.parent
_WORKFLOWS = sorted((_ROOT / ".github" / "workflows").glob("*.yml"))


def _pytest_path_args(text):
    """Yield path-like tokens passed to `pytest` invocations in a workflow.

    Comment lines are skipped. They are not executed in either context a
    workflow has -- `#` starts a YAML comment at the top level and a shell
    comment inside a `run:` block -- and prose about the suite legitimately
    mentions commands. A comment reading "a subset of the `pytest tests/` that
    every pull request runs" otherwise yielded the token "tests/`", backtick
    included, and reported it as a missing path.
    """
    for line in text.splitlines():
        if line.lstrip().startswith("#"):
            continue
        if "pytest " not in line:
            continue
        # Tokens after 'pytest' up to a flag/pipe; keep ones that look like paths
        # (contain '/' or end in .py) and are not flags or flag values.
        after = line.split("pytest ", 1)[1]
        tokens = after.replace("|", " ").split()
        skip_next = False
        for tok in tokens:
            if skip_next:
                skip_next = False
                continue
            if tok.startswith("-"):
                # flags like -m take a following value (e.g. -m "not scale")
                if tok in ("-m", "-k", "-p", "-o", "--timeout", "-n"):
                    skip_next = True
                continue
            if tok.startswith('"') or tok.startswith("'"):
                continue
            if tok.endswith(".py") or "/" in tok:
                yield tok


@pytest.mark.skipif(not _WORKFLOWS, reason="no workflow files")
@pytest.mark.parametrize("wf", _WORKFLOWS, ids=lambda p: p.name)
def test_workflow_pytest_paths_exist(wf):
    text = wf.read_text()
    missing = [tok for tok in _pytest_path_args(text) if not (_ROOT / tok).exists()]
    assert not missing, (
        f"{wf.name} references test paths that do not exist: {missing}. "
        f"Fix the path (pytest exits with code 4 on a missing file)."
    )


def test_guard_actually_inspects_a_pytest_line():
    """Sanity: the parser extracts a real path from a real workflow.

    This used to assert that the nightly workflow referenced `tests/integration`
    specifically, which pinned that workflow's CONTENT rather than the parser's
    behaviour. When nightly stopped re-running the integration suite that every
    pull request already covers, the guard failed for a change that was correct.
    A sanity check must break when the parser breaks, not when the thing it
    parses changes legitimately.
    """
    nightly = _ROOT / ".github" / "workflows" / "nightly.yml"
    if not nightly.exists():
        pytest.skip("no nightly workflow")
    toks = list(_pytest_path_args(nightly.read_text()))
    assert toks, "the parser found no pytest path at all in the nightly workflow"
    assert all((_ROOT / tok).exists() for tok in toks), toks


def test_the_parser_ignores_prose_about_pytest():
    """A comment is not an invocation.

    Verified against the exact shape that broke it: a backticked command inside
    a YAML comment.
    """
    text = "# files -- is a subset of the `pytest tests/` that every PR runs\n"
    assert list(_pytest_path_args(text)) == []


def test_the_parser_still_reads_a_real_invocation():
    """The complement of the test above: skipping comments must not skip code."""
    text = '      - name: x\n        run: pytest tests/integration/ -m "not scale" -q\n'
    assert list(_pytest_path_args(text)) == ["tests/integration/"]


def test_a_missing_path_is_still_caught():
    """The guard's whole purpose, pinned against a synthetic workflow."""
    text = "        run: pytest tests/test_does_not_exist.py -q\n"
    toks = list(_pytest_path_args(text))
    assert toks == ["tests/test_does_not_exist.py"]
    assert not (_ROOT / toks[0]).exists()
