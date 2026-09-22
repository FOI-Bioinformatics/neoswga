"""An additive coefficient must not be used outside the temperature it covers.

`require_model_support` refuses an unknown polymerase, an oligo length outside
the enzyme's parameter range, and an additive whose duplex effect nothing
computes. It did not check the TEMPERATURE a coefficient was recorded for.

The DMSO record says, in its own prose, "37 C reference, extrapolated to
30-45 C by an Arrhenius term". A Bst design at 63 C with DMSO was accepted and
returned an effective Tm of 58.72 C, computed from that coefficient forty-five
degrees above its reference and eighteen above the top of its stated range.
Nothing said so.

## Why only one additive is checked

Measured across all 23 records in the registry:

| temperature_domain | records |
|---|---|
| states an explicit range | **1** (`tm_dmso`) |
| states a single reference point | 10 |
| states no number at all | 12 |

So "replace free-text domains with checked ranges" cannot be done from the
evidence the registry holds. Doing it properly means re-reading the primary
literature, which the registry says plainly was not done when it was compiled.

Inventing the other twenty-two would be the promotion of an assumption that the
registry exists to prevent, and is the same decision taken over glycerol: the
literature expects an effect, no coefficient is invented to supply one, and the
design is refused instead.

So one range is enforced because one range is recorded. The structured field
exists for the rest to be filled in as they are verified, and a record without
one is not refused -- absence of a domain is not a domain of zero.

## What this refuses, measured before it was written

Nothing that ships. All four bundled configurations using DMSO run at 42-43 C,
inside the recorded 30-45 C. The combination it refuses is one nothing in the
repository asks for.
"""

import pytest

from neoswga.core.design_request import resolve_design_request
from neoswga.core.exceptions import UnsupportedModelError

BASE = {
    "fg_genomes": ["a.fna"],
    "fg_prefixes": ["a"],
    "num_primers": 12,
}


def design(polymerase, temp, lengths, **chemistry):
    return resolve_design_request(
        {
            **BASE,
            "polymerase": polymerase,
            "reaction_temp": temp,
            "min_k": lengths[0],
            "max_k": lengths[1],
            **chemistry,
        }
    )


# ---------------------------------------------------------------------------
# The refusal
# ---------------------------------------------------------------------------


def test_dmso_above_its_recorded_range_is_refused():
    """Bst runs at 60-65 C; the coefficient is recorded for 30-45 C."""
    with pytest.raises(UnsupportedModelError) as excinfo:
        design("bst", 63.0, (18, 20), dmso_percent=5.0)

    message = str(excinfo.value)
    assert "dmso" in message.lower()
    assert "63" in message
    assert "45" in message, "the message must name the recorded limit"


def test_the_refusal_names_the_record_it_comes_from():
    """A refusal on evidence has to say which evidence."""
    with pytest.raises(UnsupportedModelError) as excinfo:
        design("bst", 63.0, (18, 20), dmso_percent=5.0)

    assert "30" in str(excinfo.value) and "45" in str(excinfo.value)


# ---------------------------------------------------------------------------
# What must still be accepted
# ---------------------------------------------------------------------------


def test_dmso_inside_its_range_is_accepted():
    """Guard the guard: a check that refuses everything passes the tests above.

    42 C is where all four bundled DMSO configurations run.
    """
    assert design("equiphi29", 42.0, (12, 16), dmso_percent=5.0) is not None


def test_dmso_at_the_boundaries_is_accepted():
    """The recorded range is inclusive; a coefficient valid to 45 C is valid at 45."""
    assert design("phi29", 30.0, (10, 12), dmso_percent=5.0) is not None


def test_a_hot_design_without_dmso_is_accepted():
    """The refusal is about the additive, not about the temperature."""
    assert design("bst", 63.0, (18, 20)) is not None


def test_an_additive_with_no_recorded_range_is_not_refused():
    """Absence of a domain is not a domain of zero.

    Twenty-two of twenty-three records state no range. Refusing on that would
    turn a gap in the evidence into a gap in the tool, and would stop every
    additive design this project ships.
    """
    assert design("bst", 63.0, (18, 20), betaine_m=1.0) is not None


# ---------------------------------------------------------------------------
# The registry side
# ---------------------------------------------------------------------------


def test_the_recorded_range_agrees_with_its_own_prose():
    """The structured field must not drift from the sentence it came from.

    It was derived from that sentence rather than from a source, so if someone
    edits one and not the other the registry is asserting two things.
    """
    import json
    import pathlib

    records = json.loads(
        (
            pathlib.Path(__file__).resolve().parent.parent
            / "neoswga"
            / "core"
            / "registry"
            / "model_evidence.json"
        ).read_text()
    )["records"]

    dmso = next(r for r in records if r["quantity"] == "tm_dmso")

    assert dmso["temperature_range_c"] == [30.0, 45.0]
    assert "30-45 C" in dmso["temperature_domain"]


def test_no_range_was_invented_for_the_other_records():
    """The registry's own rule: recording a domain nobody stated would break,
    in the act of recording it, the thing the registry exists to enforce."""
    import json
    import pathlib

    records = json.loads(
        (
            pathlib.Path(__file__).resolve().parent.parent
            / "neoswga"
            / "core"
            / "registry"
            / "model_evidence.json"
        ).read_text()
    )["records"]

    with_range = [r["quantity"] for r in records if r.get("temperature_range_c")]

    assert with_range == ["tm_dmso"], (
        "a temperature range appeared for a record whose prose does not state "
        f"one: {with_range}. Verify it against the source first."
    )
