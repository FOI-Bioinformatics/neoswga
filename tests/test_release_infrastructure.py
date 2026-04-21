"""Tests that the release infrastructure files are wired correctly.

These are smoke tests that catch dumb mistakes (broken YAML, missing workflow
entries) before they land in CI or on PyPI.
"""

import json
import os
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def _read_yaml(path: Path):
    try:
        import yaml  # type: ignore
    except ImportError:
        pytest.skip("PyYAML not installed")
    with path.open("r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def test_ci_workflow_exists_and_parses():
    path = ROOT / ".github" / "workflows" / "ci.yml"
    assert path.is_file(), f"{path} missing"
    data = _read_yaml(path)
    # Must have lint, test, build jobs
    assert set(data["jobs"].keys()) >= {"lint", "test", "build"}


def test_no_pypi_publish_workflow():
    """Phase 16.5: releases are tagged on GitHub only; PyPI publishing is
    intentionally NOT wired up. This test locks that decision in so a
    publish.yml cannot be reintroduced without a corresponding policy
    change. If you need to publish to PyPI in the future, delete this
    test together with the new workflow file."""
    path = ROOT / ".github" / "workflows" / "publish.yml"
    assert not path.exists(), (
        f"{path} should not exist — releases are GitHub-only. "
        f"If intentionally reintroducing PyPI publishing, remove this test."
    )


def test_nightly_workflow_exists_and_parses():
    path = ROOT / ".github" / "workflows" / "nightly.yml"
    assert path.is_file(), f"{path} missing"
    data = _read_yaml(path)
    # The string 'on' is parsed as boolean True by some YAML loaders (YAML 1.1).
    on_block = data.get("on") or data.get(True)
    assert on_block is not None, f"nightly workflow has no triggers; got keys: {list(data.keys())}"
    assert "schedule" in on_block
    assert "workflow_dispatch" in on_block


def test_precommit_config_exists_and_parses():
    path = ROOT / ".pre-commit-config.yaml"
    assert path.is_file(), f"{path} missing"
    data = _read_yaml(path)
    repos = [r["repo"] for r in data["repos"]]
    assert any("ruff" in r for r in repos)
    assert any("black" in r.lower() for r in repos)
    assert any("isort" in r.lower() for r in repos)


def test_pyproject_dev_extras_include_release_tools():
    import tomllib  # Python 3.11+
    path = ROOT / "pyproject.toml"
    with path.open("rb") as fh:
        data = tomllib.load(fh)
    dev = set(data["project"]["optional-dependencies"]["dev"])
    tools = [t.split(">")[0].split("=")[0].strip() for t in dev]
    assert "ruff" in tools
    assert "mypy" in tools
    assert "pre-commit" in tools
    assert "build" in tools
    assert "twine" in tools


def test_ruff_config_present():
    import tomllib
    path = ROOT / "pyproject.toml"
    with path.open("rb") as fh:
        data = tomllib.load(fh)
    assert "ruff" in data.get("tool", {}), "Missing [tool.ruff] config"
    assert "mypy" in data.get("tool", {}), "Missing [tool.mypy] config"
