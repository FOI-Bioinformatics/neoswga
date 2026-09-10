"""Verbose about absent parameters, silent about misspelled ones.

get_value_or_default logged a debug line for each of twelve OPTIONAL_PARAMS a
params.json omits -- seven of them on the shipped config, one announcing a
default for min_amp_pred, a gate that no longer runs without --amp-model.

The other direction had no code at all: the schema sets additionalProperties
true and param_validator had no unknown-key check, so `max_bg_freqency` was
accepted in silence and the default applied, changing the design.
"""

import json
import logging

from neoswga.core.param_validator import ParamValidator, unknown_param_keys


def test_a_misspelled_key_is_reported_with_the_closest_match():
    unknown = dict(unknown_param_keys({"data_dir": ".", "max_bg_freqency": 5e-6}))

    assert "max_bg_freqency" in unknown
    assert unknown["max_bg_freqency"] == "max_bg_freq"


def test_a_declared_key_is_not_reported():
    assert unknown_param_keys({"data_dir": ".", "max_bg_freq": 5e-6}) == []


def test_an_unrecognisable_key_is_reported_without_a_suggestion():
    unknown = dict(unknown_param_keys({"data_dir": ".", "zzzzzzzzz": 1}))

    assert "zzzzzzzzz" in unknown
    assert unknown["zzzzzzzzz"] is None


def test_the_validator_reports_it_as_a_warning_not_an_error(tmp_path):
    """additionalProperties stays true; this must not reject the file."""
    from neoswga.core.param_validator import ValidationLevel

    params = {
        "data_dir": str(tmp_path),
        "fg_genomes": [],
        "fg_prefixes": [],
        "max_bg_freqency": 5e-6,
    }
    messages = ParamValidator().validate_params(params)

    unknown = [m for m in messages if m.parameter == "max_bg_freqency"]
    assert unknown, [m.parameter for m in messages]
    assert unknown[0].level == ValidationLevel.WARNING
    assert "max_bg_freq" in (unknown[0].suggestion or "")


def test_the_schema_still_allows_additional_properties():
    from neoswga.core.schema import load_schema

    assert (
        load_schema()["additionalProperties"] is True
    ), "the unknown-key check warns; it must not become a rejection"


def test_absent_optional_parameters_are_not_logged(caplog):
    """Seven lines at the head of every step log, for documented static
    defaults."""
    from neoswga.core import parameter

    with caplog.at_level(logging.DEBUG, logger="neoswga.core.parameter"):
        parameter.get_value_or_default(None, {}, "min_amp_pred")
        parameter.get_value_or_default(None, {}, "max_dimer_bp")
        parameter.get_value_or_default(None, {}, "retries")

    assert "Missing optional parameter" not in caplog.text


def test_a_genuinely_missing_required_parameter_is_still_warned(caplog):
    from neoswga.core import parameter

    with caplog.at_level(logging.WARNING, logger="neoswga.core.parameter"):
        parameter.get_value_or_default(None, {}, "max_bg_freq")

    assert "max_bg_freq" in caplog.text


def test_get_params_warns_about_an_unknown_key(tmp_path, caplog):
    """The warning has to reach someone who never runs `validate params`."""
    from neoswga.core import parameter

    fasta = tmp_path / "t.fasta"
    fasta.write_text(">t\n" + "ACGT" * 500 + "\n")
    params_file = tmp_path / "params.json"
    params_file.write_text(
        json.dumps(
            {
                "data_dir": str(tmp_path),
                "fg_genomes": [str(fasta)],
                "fg_prefixes": [str(tmp_path / "t")],
                "bg_genomes": [],
                "bg_prefixes": [],
                "fg_seq_lengths": [2000],
                "bg_seq_lengths": [],
                "max_bg_freqency": 5e-6,
                "cpus": 1,
            }
        )
    )

    class _Args:
        json_file = str(params_file)

        def __getattr__(self, name):
            return None

    with caplog.at_level(logging.WARNING, logger="neoswga.core.parameter"):
        parameter.get_params(_Args())

    assert "max_bg_freqency" in caplog.text
    assert "max_bg_freq" in caplog.text


def test_an_underscore_annotation_key_is_not_reported():
    """JSON has no comment syntax, so `_comment` is the convention people reach
    for, and two shipped example configs use it. Warning about it would train
    readers to ignore the warning."""
    assert unknown_param_keys({"data_dir": ".", "_comment": "why this config"}) == []


def test_every_shipped_example_config_is_clean():
    """The check is only useful if a correct config is silent.

    This found two real problems on first run: `mismatch_penalty` was read by
    the occupancy model and listed in OPTIONAL_PARAMS but absent from the
    schema, and `_comment` was being reported as a typo.
    """
    import glob

    offenders = {}
    for path in sorted(glob.glob("examples/**/params.json", recursive=True)):
        with open(path) as fh:
            unknown = unknown_param_keys(json.load(fh))
        if unknown:
            offenders[path] = unknown

    assert not offenders, offenders
