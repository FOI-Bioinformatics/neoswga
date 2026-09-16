"""A third retention point: keep the evenness gate, drop the arbitrary cut.

Plan step 249 of the condition-aware pool design plan asks for three modes to be
compared, and only two existed. This is the middle one.

`all_qc` indexes all 491,836 candidates that cleared hard QC on the Wolbachia
design, costing 443 MB of background index
(`docs/validation/wolbachia_retention_benchmark_2026-09-16.md`).

`post_gini` sits below it at the 20,670 that cleared the Gini gate, about 19 MB.
The distinction it draws is between a GATE and a RANKING. The Gini gate is a
declared requirement on binding evenness that a candidate either meets or does
not; `max_primer` is a cut through a ranking, chosen for the size of the working
set rather than for any property of the candidates below the line. Keeping the
first and dropping the second retains what was excluded arbitrarily and not what
was excluded on a stated rule.

It does not reintroduce the silent zero `all_qc` was built to close. Under
`post_gini` the Gini gate is a hard gate, so the candidate provider never
expands past it and never asks the background index about a candidate it does
not hold.
"""

import json
import pathlib

import pytest

from neoswga.core.candidate_inventory import background_scan_pool

A, C, G, T = "AAAAAAAAAAAA", "CCCCCCCCCCCC", "GGGGGGGGGGGG", "TTTTTTTTTTTT"

CLEARED = [A, C, G, T]
AFTER_GINI = [A, C, G]
SHORTLIST = [A]


def test_the_mode_is_declared_in_the_schema():
    schema = json.loads(pathlib.Path("neoswga/core/schema/params.schema.json").read_text())
    entry = schema["properties"]["candidate_retention"]

    assert set(entry["enum"]) == {"all_qc", "post_gini"}
    assert entry["default"] == "all_qc", "the default must not move silently"


def test_post_gini_indexes_the_evenness_survivors():
    indexed = background_scan_pool(CLEARED, "post_gini", after_gini=AFTER_GINI)

    assert sorted(indexed) == sorted(AFTER_GINI)


def test_it_sits_strictly_below_all_qc():
    """Stated as a containment, not as two separate counts."""
    middle = set(background_scan_pool(CLEARED, "post_gini", after_gini=AFTER_GINI))
    everything = set(background_scan_pool(CLEARED, "all_qc", after_gini=AFTER_GINI))

    assert middle < everything


def test_all_qc_does_not_need_the_post_gini_argument():
    assert sorted(background_scan_pool(CLEARED, "all_qc")) == sorted(CLEARED)


def test_post_gini_without_its_input_is_refused_rather_than_guessed():
    """Falling back to the other mode would be a silent change of policy."""
    with pytest.raises(ValueError, match="post_gini"):
        background_scan_pool(CLEARED, "post_gini")


def test_the_removed_mode_says_what_replaced_it_and_why():
    """A bare enum failure would not tell a user why their config stopped working."""
    with pytest.raises(ValueError) as excinfo:
        background_scan_pool(CLEARED, "legacy")

    message = str(excinfo.value)
    assert "removed" in message
    assert "post_gini" in message
    assert "empty background" in message, "the reason it went is the useful half"


def test_an_unknown_mode_names_the_modes_that_exist():
    with pytest.raises(ValueError) as excinfo:
        background_scan_pool(CLEARED, "whatever")

    message = str(excinfo.value)
    assert "all_qc" in message and "post_gini" in message
    assert "legacy" not in message, "a removed mode must not read as an option"


def test_the_pipeline_passes_the_post_gini_frame_through():
    """The frame before the cut, not the one after it.

    `gini_df` and `filtered_gini_df` differ by exactly the `max_primer` cut, and
    handing the second one to a post-Gini mode would silently index the
    shortlist instead, which is exactly the policy that was just removed.
    """
    import ast
    import inspect
    import textwrap

    from neoswga.core import pipeline

    source = textwrap.dedent(inspect.getsource(pipeline.step2))
    tree = ast.parse(source)
    calls = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and getattr(node.func, "id", None) == "_index_background_by_retention"
    ]

    assert len(calls) == 1, "the retention call site moved"
    names = [getattr(a, "id", None) for a in calls[0].args]
    assert "gini_df" in names, "the post-Gini frame does not reach the retention policy"
    assert "filtered_rate_df" in names
