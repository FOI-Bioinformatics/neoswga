"""Test that pipeline commands auto-validate params.json."""
import json
import tempfile
import os
import pytest


def test_valid_params_pass_validation():
    """A valid params.json should not trigger validation errors."""
    from neoswga.cli_unified import validate_params_json_file

    # Note the plural 'fg_genomes' — the pipeline and schema require a list.
    params = {
        'data_dir': '/tmp/test_neoswga_valid',
        'fg_genomes': ['/tmp/test.fna'],
        'fg_prefixes': ['/tmp/test'],
    }

    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        json.dump(params, f)
        tmp_path = f.name

    try:
        # Should not raise or exit
        validate_params_json_file(tmp_path)
    finally:
        os.unlink(tmp_path)
        # Clean up created directory
        if os.path.isdir('/tmp/test_neoswga_valid'):
            os.rmdir('/tmp/test_neoswga_valid')
