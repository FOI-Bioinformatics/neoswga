"""BED/BedGraph export reads the run's targets, from the params file given.

`_load_primer_positions` globbed every `*_positions.h5` in the results
directory, so a host index's sites were written under the target's name; it
read params only from `DIR/params.json`, ignoring `-j`; a primer the index did
not hold produced no line and no warning; and `_get_genome_length` fell back to
0, so BedGraph had no windows (2026-11-23 sweep, tier 3).
"""

import json

import pytest

from neoswga.cli.report import _get_genome_length, _load_primer_positions

TARGET = "ACGTACGTACTTTTTTTTTTGTACGTACGTAAAAAAAAAACCCCCC"


@pytest.fixture
def run(tmp_path):
    """A params file outside the results directory, naming one target."""
    config = tmp_path / "config"
    config.mkdir()
    (config / "target.fna").write_text(f">t\n{TARGET}\n")
    results = tmp_path / "results"
    results.mkdir()
    params = {"fg_prefixes": ["target"], "fg_genomes": ["target.fna"], "data_dir": "."}
    (config / "params.json").write_text(json.dumps(params))
    return results, config / "params.json"


def test_the_params_file_given_is_the_one_read(run):
    results, params_file = run
    positions = _load_primer_positions(str(results), ["ACGTACGTAC"], str(params_file))
    assert sorted(positions["ACGTACGTAC"]) == [(0, "forward"), (20, "reverse")]


def test_a_primer_with_no_answer_refuses_rather_than_writing_nothing(run):
    results, params_file = run
    params = json.loads(params_file.read_text())
    params["fg_genomes"] = []
    params_file.write_text(json.dumps(params))
    with pytest.raises(Exception, match="ACGTACGTAC|target"):
        _load_primer_positions(str(results), ["ACGTACGTAC"], str(params_file))


def test_two_targets_are_refused(run):
    results, params_file = run
    params = json.loads(params_file.read_text())
    params["fg_prefixes"] = ["target", "other"]
    params["fg_genomes"] = ["target.fna", "other.fna"]
    params_file.write_text(json.dumps(params))
    with pytest.raises(ValueError, match="one target"):
        _load_primer_positions(str(results), ["ACGTACGTAC"], str(params_file))


def test_no_params_file_is_refused(tmp_path):
    with pytest.raises(ValueError, match="params"):
        _load_primer_positions(str(tmp_path), ["ACGTACGTAC"], None)


def test_the_genome_length_is_read_or_measured_never_zero(run):
    results, params_file = run
    assert _get_genome_length(str(results), str(params_file)) == len(TARGET)
    params = json.loads(params_file.read_text())
    params["fg_seq_lengths"] = [5000]
    params_file.write_text(json.dumps(params))
    assert _get_genome_length(str(results), str(params_file)) == 5000
    params["fg_seq_lengths"], params["fg_genomes"] = [], []
    params_file.write_text(json.dumps(params))
    with pytest.raises(ValueError, match="length"):
        _get_genome_length(str(results), str(params_file))
