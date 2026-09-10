"""Tests for scripts/render_schema.py and the generated params-reference.md.

These lock in that the rendered reference stays in sync with the schema,
so schema changes require regenerating the doc (and the CI can fail
loudly if someone forgets).
"""

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def test_render_schema_script_exists():
    assert (ROOT / "scripts" / "render_schema.py").is_file()


def test_params_reference_is_generated():
    """docs/params-reference.md must be exactly what the renderer emits.

    This used to RUN scripts/render_schema.py, which WRITES the tracked file,
    and then assert that a few strings appeared in the result. It therefore
    asserted against content it had just written and could not fail: a stale
    doc passed, and the run left the repository dirty instead of reporting
    anything. Verified on 2026-09-10 by deleting a parameter row -- the old
    test passed and put the row back.

    It cost more than tidiness. A hand-written addition to the "Additional
    guidance" section was silently reverted by a suite run, so the commit
    meant to carry it carried nothing. Prose for that file belongs in the
    renderer's FOOTER, which is where it now lives.

    Comparing in memory instead: the renderer exposes `render()` as a pure
    function, so nothing needs to be written to check the committed file.
    """
    import importlib.util

    ref = ROOT / "docs" / "params-reference.md"
    assert ref.is_file(), "docs/params-reference.md missing; run scripts/render_schema.py"

    spec = importlib.util.spec_from_file_location(
        "_render_schema", ROOT / "scripts" / "render_schema.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    with (ROOT / "neoswga" / "core" / "schema" / "params.schema.json").open() as fh:
        schema = json.load(fh)
    expected = module.render(schema)

    content = ref.read_text()
    assert content == expected, (
        "docs/params-reference.md is out of date with params.schema.json. "
        "Regenerate it: python scripts/render_schema.py"
    )

    # Structural assertions that both schema and renderer emit
    assert "## Required parameters" in content
    assert "## Optional parameters" in content
    assert "`data_dir`" in content
    assert "`fg_genomes`" in content
    assert "`polymerase`" in content


def test_params_reference_in_sync_with_schema():
    """Regenerate and confirm no diff — stale docs are a signal error."""
    ref = ROOT / "docs" / "params-reference.md"
    before = ref.read_text() if ref.is_file() else ""

    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "render_schema.py")],
        capture_output=True,
        text=True,
        timeout=30,
        cwd=ROOT,
    )
    assert result.returncode == 0

    after = ref.read_text()
    if before != after:
        pytest.fail(
            "docs/params-reference.md is out of sync with the schema. "
            "Run `python scripts/render_schema.py` and commit the result."
        )


def test_lock_file_exists():
    """requirements-dev.lock should be committed for reproducible CI envs."""
    lock = ROOT / "requirements-dev.lock"
    assert lock.is_file(), (
        "requirements-dev.lock missing; "
        "run `pip-compile --extra dev --output-file requirements-dev.lock pyproject.toml`"
    )
    content = lock.read_text()
    # Must pin some expected packages
    assert "numpy==" in content
    assert "pytest==" in content
    assert "jsonschema==" in content
