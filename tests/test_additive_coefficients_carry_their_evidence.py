"""Every additive coefficient states how well it is supported.

Audit finding: several exact numerical attributions in `mechanistic_params` are
wrong or could not be verified, and nothing in the code distinguished a
coefficient measured in the relevant domain from one extrapolated into it or one
with no located source at all.

The audit is explicit that an unverified coefficient is not thereby false. It is
insufficiently supported for the certainty currently attached to it, and
substituting a new point value merely because it looks plausible repeats the
error. So this does not change any number: it requires each to declare what
stands behind it.

Categories:

``measured``
    A source reports this quantity, at a concentration and duplex length
    relevant to SWGA.
``extrapolated``
    A source reports the quantity, but outside the domain it is used in --
    long DNA, PCR amplicons, or a concentration range the tool does not reach.
``unsupported``
    No located source reports this value. Retained for continuity, not to be
    presented as a literature constant.
"""

import pytest

from neoswga.core.mechanistic_params import ADDITIVE_TM_PARAMS

CATEGORIES = {"measured", "extrapolated", "unsupported"}


@pytest.mark.parametrize("additive", sorted(ADDITIVE_TM_PARAMS))
def test_each_coefficient_declares_an_evidence_category(additive):
    params = ADDITIVE_TM_PARAMS[additive]

    assert "evidence" in params, f"{additive} has no evidence category; see this module's docstring"
    assert (
        params["evidence"] in CATEGORIES
    ), f"{additive} declares {params['evidence']!r}, not one of {sorted(CATEGORIES)}"


@pytest.mark.parametrize("additive", sorted(ADDITIVE_TM_PARAMS))
def test_each_coefficient_names_a_source(additive):
    """`unsupported` is the one category allowed to name none."""
    params = ADDITIVE_TM_PARAMS[additive]
    source = params.get("source", "")

    if params.get("evidence") == "unsupported":
        return
    assert source, f"{additive} is not marked unsupported but names no source"


def test_ethanol_is_not_presented_as_a_literature_constant():
    """Its cited source is a long-PCR study describing glycerol and DMSO.

    Pinned by name because it is the clearest case the audit found: the value
    is in the code, the citation does not report it, and no replacement source
    was located.
    """
    assert ADDITIVE_TM_PARAMS["ethanol"]["evidence"] == "unsupported"
