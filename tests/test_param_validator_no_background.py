"""Tests for the no-background-config validation guard.

When a user runs the pipeline without any background source (no bg_genomes,
no bg_prefixes, no bloom_filter_path), background-related filters such as
max_bg_freq and min_fg_to_bg_ratio are silently inert. The validator should
report this clearly so users do not assume they have specificity filtering
when they do not.
"""

from neoswga.core.param_validator import ParamValidator, ValidationLevel


REQUIRED = {
    "fg_genomes": ["foo.fna"],
    "fg_prefixes": ["foo"],
    "data_dir": "results",
}


def _validate(params):
    return ParamValidator().validate_params({**REQUIRED, **params})


def _info_messages_for(messages, parameter_name):
    return [
        m for m in messages
        if m.level == ValidationLevel.INFO and m.parameter == parameter_name
    ]


def test_absent_bg_keys_emit_no_background_info():
    messages = _validate({})
    msgs = _info_messages_for(messages, "bg_genomes")
    assert len(msgs) == 1
    assert "no-background" in msgs[0].message.lower() or "background-free" in msgs[0].message.lower()


def test_empty_bg_lists_emit_no_background_info():
    messages = _validate({"bg_genomes": [], "bg_prefixes": []})
    msgs = _info_messages_for(messages, "bg_genomes")
    assert len(msgs) == 1


def test_bg_genomes_present_does_not_emit_info():
    messages = _validate({
        "bg_genomes": ["host.fna"],
        "bg_prefixes": ["host"],
    })
    assert _info_messages_for(messages, "bg_genomes") == []


def test_bloom_filter_path_satisfies_check():
    messages = _validate({"bloom_filter_path": "host.bloom"})
    assert _info_messages_for(messages, "bg_genomes") == []
