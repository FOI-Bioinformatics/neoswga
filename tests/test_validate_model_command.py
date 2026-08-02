"""`neoswga validate-model` must pass its own checks.

The command exists to tell a user whether the mechanistic model still behaves
the way the literature says it should. A standing FAIL is worse than no check:
it either trains the reader to ignore the output, or invites them to "fix" the
model back to whatever the stale expectation asserts.

That is what had happened. The Mg2+ expectations were left at the PCR figures
(optimum 2-4 mM, decline above 6 mM) when the model was recalibrated for
isothermal amplification in schema v2 -- where the phi29 buffer value of 10 mM
is the optimum and 2.5 mM sits below the low-activity threshold. The test
therefore asserted that 8 mM must be worse than 2.5 mM, the opposite of the
corrected behaviour, and the command reported 6/7 against a model that was
right.
"""

import pytest

from neoswga.core.model_validation import validate_mechanistic_model


def test_all_model_validations_pass():
    failed = [t["test"] for t in validate_mechanistic_model() if not t["passed"]]

    assert not failed, f"model validation failing: {failed}"


def test_magnesium_optimum_is_the_isothermal_value():
    """Pins the direction of the correction.

    Reverting the model to the PCR optimum would make this fail, which is the
    point: the two numbers differ because the chemistry differs, not because
    one of them is a typo.
    """
    from neoswga.core.mechanistic_params import MECHANISTIC_MODEL_PARAMS

    enzyme = MECHANISTIC_MODEL_PARAMS["enzyme"]
    assert enzyme["mg_optimal"] == pytest.approx(10.0)
    assert enzyme["mg_low_threshold"] < 10.0 < enzyme["mg_high_threshold"]
