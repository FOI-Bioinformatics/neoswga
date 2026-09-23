"""The minimum-site threshold for the Gini gate is settable, defaulting to 3.

Audit finding B4 established the threshold. This pins its routing, because this
repository's recurring defect is a config key that is declared in the schema,
accepted by the validator, and read by nothing: `additionalProperties: true`
means such a key produces no warning. See Known Issue 8 in CLAUDE.md.

The value must survive the multiprocessing boundary. `get_gini_from_txt` runs
`get_gini_from_txt_for_one_k` through a Pool, and macOS spawns rather than
forks, so a worker reading `parameter.min_gini_sites` would see the module's
import-time default and the parallel path would disagree with the serial one on
the same primers. The parent resolves it and passes it as data.
"""

import json

import pytest

from neoswga.core import parameter
from neoswga.core import primer_attributes as pa


def test_the_default_is_three():
    assert pa.DEFAULT_MIN_GINI_SITES == 3


def test_the_schema_declares_it():
    import pathlib

    schema = json.loads(
        (pathlib.Path(pa.__file__).parent / "schema" / "params.schema.json").read_text()
    )
    prop = schema["properties"]["min_gini_sites"]
    assert prop["type"] == "integer"
    assert prop["minimum"] == 1


def test_params_json_sets_the_module_global(tmp_path, monkeypatch):
    params = {
        "schema_version": 2,
        "data_dir": str(tmp_path),
        "fg_genomes": [],
        "bg_genomes": [],
        "fg_prefixes": [],
        "bg_prefixes": [],
        "fg_seq_lengths": [],
        "bg_seq_lengths": [],
        "min_k": 12,
        "max_k": 12,
        "min_gini_sites": 5,
    }
    path = tmp_path / "params.json"
    path.write_text(json.dumps(params))

    class _Args:
        json_file = str(path)

        def __getattr__(self, name):
            return None

    parameter.get_params(_Args())
    assert parameter.min_gini_sites == 5


def test_an_absent_key_leaves_the_default(tmp_path):
    params = {
        "schema_version": 2,
        "data_dir": str(tmp_path),
        "fg_genomes": [],
        "bg_genomes": [],
        "fg_prefixes": [],
        "bg_prefixes": [],
        "fg_seq_lengths": [],
        "bg_seq_lengths": [],
        "min_k": 12,
        "max_k": 12,
    }
    path = tmp_path / "params.json"
    path.write_text(json.dumps(params))

    class _Args:
        json_file = str(path)

        def __getattr__(self, name):
            return None

    parameter.get_params(_Args())
    assert parameter.min_gini_sites == pa.DEFAULT_MIN_GINI_SITES


def test_the_helper_carries_the_value_across_the_process_boundary():
    """A 7-element task tuple, so the worker is told rather than asked."""
    PRIMER = "TTGACCATGA"
    cache = {("x", PRIMER): [500, 6_000]}
    task = ([PRIMER], "x", None, 10_000, True, cache, 2)
    result = pa.get_gini_from_txt_for_one_k_helper(task)
    forward, _ = result[PRIMER]
    import math

    assert not math.isnan(forward), "min_sites=2 must make two sites measurable"

    task = ([PRIMER], "x", None, 10_000, True, cache, 3)
    forward, _ = pa.get_gini_from_txt_for_one_k_helper(task)[PRIMER]
    assert math.isnan(forward), "min_sites=3 must make two sites unmeasurable"


def test_a_six_element_task_still_works():
    """The helper's existing 5- and 6-element forms must keep working."""
    PRIMER = "TTGACCATGA"
    cache = {("x", PRIMER): [500, 4_000, 8_000]}
    result = pa.get_gini_from_txt_for_one_k_helper(([PRIMER], "x", None, 10_000, True, cache))
    assert PRIMER in result


def test_the_cli_flag_parses_on_filter():
    """Parse a real command line rather than inspecting a parser's internals.

    The plan's version called `cli_pipeline.add_pipeline_parsers`, which does
    not exist, and reached into `_actions`. Going through the real parser is
    what a user does and is what would notice the flag being registered on the
    wrong subcommand.
    """
    from neoswga import cli_unified

    args = cli_unified.create_parser().parse_args(
        ["filter", "-j", "params.json", "--min-gini-sites", "2"]
    )
    assert args.min_gini_sites == 2


def test_the_flag_beats_params_json(tmp_path):
    """A flag the user typed must win over the file, or the escape hatch this
    exists to provide does not work on the command line."""
    params = {
        "schema_version": 2,
        "data_dir": str(tmp_path),
        "fg_genomes": [],
        "bg_genomes": [],
        "fg_prefixes": [],
        "bg_prefixes": [],
        "fg_seq_lengths": [],
        "bg_seq_lengths": [],
        "min_k": 12,
        "max_k": 12,
        "min_gini_sites": 5,
    }
    path = tmp_path / "params.json"
    path.write_text(json.dumps(params))

    class _Args:
        json_file = str(path)
        min_gini_sites = 2

        def __getattr__(self, name):
            return None

    parameter.get_params(_Args())
    assert parameter.min_gini_sites == 2
