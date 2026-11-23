"""`seed` in params.json seeds the run; request-only keys are not called unknown.

`resolve_design_request` hashes params.json's `seed` into the request, but
`optimize` seeded only from `--seed`, so a file declaring `"seed": 42` recorded
42 in the request identity and ran unseeded. And `seed`, `fixed_oligos` and
`excluded_oligos`, which the design request reads, are not schema keys, so
loading a file that set them warned "Unknown parameter ... it will be ignored",
which was false (2026-11-23 sweep, tier 2).
"""

import argparse
import json
import pytest

from neoswga.cli._common import resolve_seed
from neoswga.core.param_validator import unknown_param_keys


def _args(tmp_path, flag, params):
    path = tmp_path / "params.json"
    path.write_text(json.dumps(params))
    return argparse.Namespace(seed=flag, json_file=str(path))


def test_the_flag_wins(tmp_path):
    assert resolve_seed(_args(tmp_path, 7, {"seed": 42})) == 7


def test_the_file_seeds_without_the_flag(tmp_path):
    assert resolve_seed(_args(tmp_path, None, {"seed": 42})) == 42


def test_neither_means_unseeded(tmp_path):
    assert resolve_seed(_args(tmp_path, None, {})) is None


@pytest.mark.parametrize(
    "key", ["seed", "fixed_oligos", "excluded_oligos", "concentration_mode", "total_primer_molar"]
)
def test_a_request_only_key_is_not_unknown(key):
    assert key not in [found for found, _ in unknown_param_keys({key: 1})]


def test_a_typo_is_still_unknown():
    assert "max_bg_freqency" in [found for found, _ in unknown_param_keys({"max_bg_freqency": 1})]
