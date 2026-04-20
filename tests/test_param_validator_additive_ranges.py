"""Tests for the newly added additive and reaction-condition ranges in
PARAM_RANGES. These bounds mirror ReactionConditions._validate() so callers
hit ParamValidator errors rather than a deep ValueError from the pipeline.
"""

import pytest

from neoswga.core.param_validator import ParamValidator, ValidationLevel, PARAM_RANGES


REQUIRED = {
    'fg_genomes': ['foo.fna'],
    'fg_prefixes': ['foo'],
    'data_dir': 'results',
}


def _validate(params):
    validator = ParamValidator()
    return validator.validate_params({**REQUIRED, **params})


def _has_range_error(messages, param):
    return any(
        m.level == ValidationLevel.ERROR and param in m.parameter
        for m in messages
    )


@pytest.mark.parametrize("param", [
    'formamide_percent', 'ethanol_percent', 'urea_m', 'tmac_m',
    'glycerol_percent', 'peg_percent', 'bsa_ug_ml', 'mg_conc',
    'gc_tolerance', 'genome_gc', 'bl_penalty', 'max_bl_freq',
])
def test_new_param_has_range(param):
    """New additive / reaction-condition parameters must have an explicit range."""
    assert param in PARAM_RANGES, (
        f"{param} should be in PARAM_RANGES so invalid values are caught by the validator"
    )


def test_formamide_upper_bound_rejected():
    messages = _validate({'formamide_percent': 11.0})
    assert _has_range_error(messages, 'formamide_percent')


def test_formamide_in_range_accepted():
    messages = _validate({'formamide_percent': 5.0})
    assert not _has_range_error(messages, 'formamide_percent')


def test_urea_upper_bound_rejected():
    messages = _validate({'urea_m': 2.5})
    assert _has_range_error(messages, 'urea_m')


def test_tmac_upper_bound_rejected():
    messages = _validate({'tmac_m': 0.15})
    assert _has_range_error(messages, 'tmac_m')


def test_mg_conc_upper_bound_rejected():
    messages = _validate({'mg_conc': 25.0})
    assert _has_range_error(messages, 'mg_conc')


def test_mg_conc_in_range_accepted():
    messages = _validate({'mg_conc': 10.0})
    assert not _has_range_error(messages, 'mg_conc')


def test_gc_tolerance_out_of_range_rejected():
    messages = _validate({'gc_tolerance': 0.6})
    assert _has_range_error(messages, 'gc_tolerance')


def test_genome_gc_out_of_range_rejected():
    messages = _validate({'genome_gc': 1.2})
    assert _has_range_error(messages, 'genome_gc')
