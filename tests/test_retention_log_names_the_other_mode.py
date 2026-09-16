"""The retention log must not describe what happened as a hypothetical.

The message names the mode that was not used, so a reader can see what the
alternative would have cost. It used to name `post_gini` unconditionally, which
under `post_gini` produced "'post_gini' would index the 20670" on a line that
had just reported indexing exactly 20,670. Read literally it says the run might
do what it did.

Cheap to get wrong again, because the sentence is only wrong in one of the two
modes and the other reads fine.
"""

import logging

import pandas as pd
import pytest

from neoswga.core import pipeline, stage2_recording

CLEARED = ["AAAAAAAAAAAA", "CCCCCCCCCCCC", "GGGGGGGGGGGG", "TTTTTTTTTTTT"]
AFTER_GINI = CLEARED[:3]
SHORTLIST = CLEARED[:1]


@pytest.fixture
def indexed(monkeypatch, caplog):
    """Run the retention step with the background scan stubbed out."""

    def run(retention):
        monkeypatch.setattr(
            stage2_recording.parameter, "candidate_retention", retention, raising=False
        )
        monkeypatch.setattr(pipeline, "_scan_background_positions", lambda *a, **k: None)
        caplog.clear()
        with caplog.at_level(logging.INFO, logger=stage2_recording.logger.name):
            selected = stage2_recording._index_background_by_retention(
                pd.DataFrame({"primer": CLEARED}),
                pd.DataFrame({"primer": AFTER_GINI}),
                pd.DataFrame({"primer": SHORTLIST}),
                ["bg"],
                ["bg.fna"],
            )
        message = next(
            r.getMessage() for r in caplog.records if "Background index" in r.getMessage()
        )
        return selected, message

    return run


def test_all_qc_names_post_gini_as_the_alternative(indexed):
    selected, message = indexed("all_qc")

    assert sorted(selected) == sorted(CLEARED)
    assert f"indexed under candidate_retention='all_qc'" in message
    assert f"'post_gini' would index the {len(AFTER_GINI)}" in message


def test_post_gini_names_all_qc_as_the_alternative(indexed):
    selected, message = indexed("post_gini")

    assert sorted(selected) == sorted(AFTER_GINI)
    assert f"indexed under candidate_retention='post_gini'" in message
    assert f"'all_qc' would index all {len(CLEARED)}" in message


@pytest.mark.parametrize("retention", ["all_qc", "post_gini"])
def test_the_message_never_offers_the_mode_it_just_ran(indexed, retention):
    """The defect, stated directly rather than through either example."""
    _, message = indexed(retention)

    assert (
        f"'{retention}' would index" not in message
    ), f"the log offers {retention!r} as an alternative to itself"


@pytest.mark.parametrize("retention", ["all_qc", "post_gini"])
def test_the_count_reported_is_the_count_indexed(indexed, retention):
    """Guard the guard: a message naming the right mode with the wrong number
    would still be wrong, and both modes index a different amount here."""
    selected, message = indexed(retention)

    assert message.startswith(f"Background index: {len(selected)} of {len(CLEARED)}")
