"""Tests for the canonical params.json JSON Schema and its validator wiring."""

import json
import subprocess
import sys

import pytest


def test_schema_ships_as_valid_json():
    """The shipped schema file is valid JSON and names the required keys."""
    from neoswga.core.schema import load_schema
    schema = load_schema()
    assert schema["$schema"].startswith("https://json-schema.org/")
    assert "data_dir" in schema["required"]
    assert "fg_genomes" in schema["required"]
    assert "fg_prefixes" in schema["required"]
    # Polymerase enum includes both phi29 and equiphi29
    polymerase_enum = schema["properties"]["polymerase"]["enum"]
    assert "phi29" in polymerase_enum
    assert "equiphi29" in polymerase_enum


def test_schema_passes_plasmid_example_params():
    """The real plasmid_example/params.json must pass the schema."""
    import os
    from neoswga.core.schema import load_schema
    try:
        import jsonschema
    except ImportError:
        pytest.skip("jsonschema not installed")

    params_path = os.path.join(
        os.path.dirname(__file__),
        '..', 'examples', 'plasmid_example', 'params.json',
    )
    if not os.path.isfile(params_path):
        pytest.skip(f"Example params not found at {params_path}")

    with open(params_path, encoding='utf-8') as fh:
        params = json.load(fh)

    schema = load_schema()
    # Use iter_errors so we surface all violations if any
    errors = list(jsonschema.Draft202012Validator(schema).iter_errors(params))
    assert not errors, "plasmid_example params.json should match the schema; got: " + \
        "; ".join(e.message for e in errors)


def test_param_validator_flags_schema_violation():
    """An invalid polymerase value should be caught by the schema layer."""
    try:
        import jsonschema  # noqa
    except ImportError:
        pytest.skip("jsonschema not installed")

    from neoswga.core.param_validator import ParamValidator, ValidationLevel
    v = ParamValidator()
    messages = v.validate_params({
        'data_dir': 'results',
        'fg_genomes': ['foo.fna'],
        'fg_prefixes': ['foo'],
        'polymerase': 'not_a_real_polymerase',
    })
    errors = [m for m in messages if m.level == ValidationLevel.ERROR]
    assert any('polymerase' in m.parameter for m in errors), (
        f"Expected polymerase error; got: {[str(m) for m in errors]}"
    )


def test_schema_cli_dump_produces_valid_json():
    """`neoswga schema --dump` emits parseable JSON on stdout."""
    result = subprocess.run(
        [sys.executable, '-m', 'neoswga.cli_unified', 'schema', '--dump'],
        capture_output=True, text=True, timeout=30,
    )
    assert result.returncode == 0, f"schema --dump failed: {result.stderr}"
    schema = json.loads(result.stdout)
    assert schema.get("title") == "neoswga params.json"


def test_schema_cli_output_writes_file(tmp_path):
    """`neoswga schema --dump -o FILE` writes the schema."""
    out = tmp_path / "params.schema.json"
    result = subprocess.run(
        [sys.executable, '-m', 'neoswga.cli_unified', 'schema',
         '--dump', '--output', str(out)],
        capture_output=True, text=True, timeout=30,
    )
    assert result.returncode == 0, f"schema -o failed: {result.stderr}"
    assert out.is_file()
    schema = json.loads(out.read_text())
    assert "properties" in schema
