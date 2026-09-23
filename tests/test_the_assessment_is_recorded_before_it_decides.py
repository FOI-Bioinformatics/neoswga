"""The assessment is written, and is deliberately not yet the gate.

`evaluate_panel` is described as the authoritative record of whether a panel is
acceptable and had no production caller: only its own two test files imported
it. It escaped `test_no_capability_is_unreachable` because that check walks an
eight-module `WATCHED` tuple which did not list it, not because anything
allowed it. Task 1 gave `run_optimization` the request it needs; this records
what it produces.

Recording it without letting it decide is the whole of this increment, and the
reason is that the two records check different things. `validate_result` owns
duplicates, blacklist re-injection and size drift in either direction;
`evaluate_panel` owns the configured panel limits, the units and basis of every
metric, a non-finite refusal, and the distinction between an unavailable
quantity and a measured zero. Neither subsumes the other, so letting the
assessment decide today would change what `export` refuses in cases nobody has
enumerated. Measuring which cases is the next increment.

The gate is `BLOCKING_VALIDATOR_CODES`, not `ok`. `export_is_blocked` selects
issues whose `code` is in that set; `ok` has one reader, which stores it on
`metrics.validation_ok` and never acts on it. So this increment is
verdict-neutral exactly while it adds no issue at all, which is what the second
test asserts.
"""

import json

import pytest

from neoswga.core.design_result import VALIDATION_FILENAME
from tests.conftest import plasmid_example_ready


@pytest.mark.skipif(not plasmid_example_ready(), reason="needs jellyfish and the plasmid example")
def test_the_artifact_carries_the_assessment(plasmid_run):
    payload = json.loads((plasmid_run / VALIDATION_FILENAME).read_text())

    assert "assessment" in payload, "the assessment was not recorded"
    assessment = payload["assessment"]
    assert assessment["request_hash"], "the assessment names no request"
    assert isinstance(assessment["qualified"], bool)
    assert assessment["primers"], "the assessment names no panel"


@pytest.mark.skipif(not plasmid_example_ready(), reason="needs jellyfish and the plasmid example")
def test_the_assessment_describes_the_panel_that_was_delivered(plasmid_run):
    """A record of a different panel is worse than no record.

    The assessment is built from `result.primers` and the CSV is written from
    the same result, so this is cheap to assert and expensive to get wrong:
    four commands once read every alternative set as one panel, and the export
    named eighteen oligos when the orderable set was eight.
    """
    import pandas as pd

    payload = json.loads((plasmid_run / VALIDATION_FILENAME).read_text())
    frame = pd.read_csv(plasmid_run / "step4_improved_df.csv")
    if "set_index" in frame.columns:
        frame = frame[frame["set_index"] == 0]

    assert sorted(payload["assessment"]["primers"]) == sorted(frame["primer"].tolist())


@pytest.mark.skipif(not plasmid_example_ready(), reason="needs jellyfish and the plasmid example")
def test_every_metric_says_what_it_is(plasmid_run):
    """The reason for the record: a coverage figure without its reach carries
    almost no information. One saved 26-oligo panel reads 41.3% at 1 kb and
    93.5% at 5 kb."""
    metrics = json.loads((plasmid_run / VALIDATION_FILENAME).read_text())["assessment"]["metrics"]

    assert "fg_coverage" in metrics
    for name, measurement in metrics.items():
        assert measurement["units"], f"{name} carries no units"
        # Either a value or a stated reason for not having one, never both and
        # never neither -- `Measurement.__post_init__` enforces it and this
        # asserts it survives the round trip through JSON.
        has_value = measurement["value"] is not None
        has_reason = bool(measurement["unavailable"])
        assert has_value != has_reason, f"{name} is {measurement}"


@pytest.mark.skipif(not plasmid_example_ready(), reason="needs jellyfish and the plasmid example")
def test_recording_it_adds_no_blocking_code(plasmid_run):
    """Verdict-neutral, asserted against the gate rather than against `ok`."""
    from neoswga.core.design_result import blocking_validator_findings

    payload = json.loads((plasmid_run / VALIDATION_FILENAME).read_text())
    codes = [issue.get("code") for issue in payload.get("issues") or []]

    assert "panel_assessment_violation" not in codes
    assert blocking_validator_findings(str(plasmid_run)) == []
