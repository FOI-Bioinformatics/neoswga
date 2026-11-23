"""An absent filtering or set-count key takes the default the documentation states.

`get_value_or_default` returned None for every absent key outside a short
optional list. For `min_fg_freq`, `max_bg_freq`, `max_gini` and `max_primer`
that made a minimal params.json crash `filter` with a TypeError; for `max_sets`
and `iterations` the readers' own fallbacks delivered one set (the request hash
recorded 5) and a 100-iteration bound (the schema says 8)
(2026-11-23 sweep, tier 1 item 7).
"""

import json
import pathlib

import pytest

from neoswga.core import parameter
from neoswga.core.parameter_config import PipelineParameters

DOCUMENTED = {
    "min_fg_freq": 1e-5,
    "max_bg_freq": 5e-6,
    "max_gini": 0.7,
    "max_primer": 500,
    "max_sets": 5,
    "iterations": 8,
    "cpus": 1,
}


class _Args:
    def __getattr__(self, name):
        return None


@pytest.fixture
def minimal(tmp_path, monkeypatch):
    (tmp_path / "fg.fna").write_text(">fg\n" + "ACGT" * 50 + "\n")
    params_file = tmp_path / "params.json"
    params_file.write_text(
        json.dumps({"fg_genomes": ["fg.fna"], "fg_prefixes": ["fg"], "data_dir": "./"})
    )
    monkeypatch.chdir(tmp_path)
    args = _Args()
    args.json_file = str(params_file)
    return parameter.get_params(args)


@pytest.mark.parametrize("key,value", sorted(DOCUMENTED.items()))
def test_absent_key_takes_its_documented_default(minimal, key, value):
    assert minimal[key] == value
    assert getattr(parameter, key) == value


@pytest.mark.parametrize("key,value", sorted(DOCUMENTED.items()))
def test_the_three_statements_of_the_default_agree(key, value):
    schema_path = pathlib.Path(parameter.__file__).parent / "schema" / "params.schema.json"
    schema = json.loads(schema_path.read_text())["properties"]
    assert PipelineParameters.__dataclass_fields__[key].default == value
    if key != "cpus":
        assert schema[key].get("default") == value
