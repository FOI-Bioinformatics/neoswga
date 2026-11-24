"""`design` takes `-n`, `-m` and `--seed` as `optimize` does.

The one-shot command had none of the three, so a size, method or seed could be
set only in params.json, and a script moving from `optimize -n 8` to `design`
lost the flag. Step 3 also seeded from `--seed` alone and ignored a params.json
`seed`, which step 4 has honoured since 2026-11-24 (2026-11-23 sweep, tier 3).
"""

import argparse
import json

import pytest

from neoswga.cli_unified import create_parser
from neoswga.core.design_request import design_request_for_run


def test_the_three_flags_parse():
    args = create_parser().parse_args(
        ["design", "-j", "p.json", "-n", "4", "-m", "dominating-set", "--seed", "7"]
    )
    assert (args.num_primers, args.optimization_method, args.seed) == (4, "dominating-set", 7)


def test_an_absent_flag_is_none():
    args = create_parser().parse_args(["design", "-j", "p.json"])
    assert (args.num_primers, args.optimization_method, args.seed) == (None, None, None)


def test_a_non_positive_size_is_refused():
    with pytest.raises(SystemExit) as info:
        create_parser().parse_args(["design", "-j", "p.json", "-n", "0"])
    assert info.value.code == 2


@pytest.fixture
def params_file(tmp_path):
    path = tmp_path / "params.json"
    path.write_text(
        json.dumps(
            {
                "fg_genomes": ["t.fna"],
                "fg_prefixes": ["t"],
                "fg_seq_lengths": [1000],
                "data_dir": str(tmp_path),
            }
        )
    )
    return str(path)


def test_design_flags_move_the_request_as_optimize_flags_do(params_file):
    def request(*argv):
        return design_request_for_run(create_parser().parse_args(list(argv)), None).request_hash

    for flags in (["-m", "dominating-set"], ["--seed", "7"]):
        assert request("design", "-j", params_file, *flags) == request(
            "optimize", "-j", params_file, *flags
        )
        assert request("design", "-j", params_file, *flags) != request("design", "-j", params_file)


class _Stop(BaseException):
    pass


def test_step_three_seeds_from_the_params_file(tmp_path, monkeypatch):
    from neoswga.cli import pipeline as cli_pipeline
    from neoswga.core import pipeline, rf_preprocessing

    path = tmp_path / "params.json"
    path.write_text(
        json.dumps(
            {"seed": 11, "data_dir": str(tmp_path), "fg_genomes": ["t.fna"], "fg_prefixes": ["t"]}
        )
    )
    seen = []

    def _record(seed):
        seen.append(seed)
        raise _Stop

    monkeypatch.setattr(pipeline, "_initialize", lambda: None)
    monkeypatch.setattr(rf_preprocessing, "set_kmer_sampling_seed", _record)
    with pytest.raises(_Stop):
        cli_pipeline.run_step3(argparse.Namespace(json_file=str(path), seed=None, quiet=True))
    assert seen == [11]
