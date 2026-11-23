"""`evaluate-set` treats targets as circular or linear as params.json says.

Every coverage call in `evaluate-set` used `not args.linear`, so with `-j` a
params file declaring `fg_circular: false` was evaluated as circular, and its
coverage differed from the `optimize` figure for the same set (0.6497 against
0.6177 on the plasmid example; 2026-11-23 sweep, tier 2).
"""

import argparse
import json

import pytest

from neoswga.cli.evaluate import _targets_circular


def _args(tmp_path, linear=False, **params):
    path = None
    if params is not None:
        path = tmp_path / "params.json"
        path.write_text(json.dumps(params))
    return argparse.Namespace(linear=linear, json_file=str(path) if params else None)


def test_the_flag_makes_targets_linear(tmp_path):
    assert _targets_circular(_args(tmp_path, linear=True, fg_circular=True)) is False


@pytest.mark.parametrize("value", [True, False])
def test_params_json_decides_without_the_flag(tmp_path, value):
    assert _targets_circular(_args(tmp_path, fg_circular=value)) is value


def test_the_default_without_either_is_circular(tmp_path):
    assert _targets_circular(argparse.Namespace(linear=False, json_file=None)) is True


def test_a_params_file_without_the_key_keeps_the_default(tmp_path):
    assert _targets_circular(_args(tmp_path, data_dir="x")) is True
