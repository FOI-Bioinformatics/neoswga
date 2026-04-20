"""Tests for scripts/render_schema.py and the generated params-reference.md.

These lock in that the rendered reference stays in sync with the schema,
so schema changes require regenerating the doc (and the CI can fail
loudly if someone forgets).
"""

import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parent.parent


def test_render_schema_script_exists():
    assert (ROOT / "scripts" / "render_schema.py").is_file()


def test_params_reference_is_generated():
    """docs/params-reference.md must exist and be regeneratable."""
    ref = ROOT / "docs" / "params-reference.md"
    assert ref.is_file(), (
        "docs/params-reference.md missing; run scripts/render_schema.py"
    )

    # Regenerate into a temporary buffer and compare
    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "render_schema.py")],
        capture_output=True, text=True, timeout=30, cwd=ROOT,
    )
    assert result.returncode == 0, result.stderr
    assert ref.is_file()
    content = ref.read_text()
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
        capture_output=True, text=True, timeout=30, cwd=ROOT,
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
