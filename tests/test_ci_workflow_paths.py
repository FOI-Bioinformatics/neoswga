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
    """Yield path-like tokens passed to `pytest` invocations in a workflow."""
    for line in text.splitlines():
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
    """Sanity: the parser finds the integration path in the nightly workflow."""
    nightly = _ROOT / ".github" / "workflows" / "nightly.yml"
    if not nightly.exists():
        pytest.skip("no nightly workflow")
    toks = list(_pytest_path_args(nightly.read_text()))
    assert any("tests/integration" in t for t in toks), toks
