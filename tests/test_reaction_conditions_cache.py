"""Test that reaction conditions cache is properly invalidated."""
import pytest
from unittest.mock import patch
from neoswga.core import parameter


def _set_default_params():
    """Set minimal valid parameter values for filter module."""
    parameter.reaction_temp = 30.0
    parameter.polymerase = 'phi29'
    parameter.dmso_percent = 0.0
    parameter.betaine_m = 0.0
    parameter.trehalose_m = 0.0
    parameter.formamide_percent = 0.0
    parameter.ethanol_percent = 0.0
    parameter.urea_m = 0.0
    parameter.tmac_m = 0.0
    parameter.na_conc = 50.0
    parameter.mg_conc = 0.0


def test_reset_reaction_conditions_clears_cache():
    from neoswga.core.filter import _get_reaction_conditions, reset_reaction_conditions

    _set_default_params()

    # Ensure clean state
    reset_reaction_conditions()

    # First call creates cache
    rc1 = _get_reaction_conditions()
    # Same object returned (cached)
    rc2 = _get_reaction_conditions()
    assert rc1 is rc2

    # Reset clears cache
    reset_reaction_conditions()
    rc3 = _get_reaction_conditions()
    assert rc3 is not rc1
