"""A user who passes a Bloom filter to `optimize` must be told where it is read.

These three options are unimplemented on the optimize path and say so, both in
their help text and at run time through `report_unimplemented_options`. That
convention is deliberate and is not changed here.

What the messages did not say is where the capability the user asked for DOES
live. They pointed at `--no-bg-prefilter`, which controls optimize's own
candidate prefilter, not the pre-built background index the flag names. A user
who spent hours building a filter and handed it to the wrong command was told
their flag does nothing and not where to put it.

The Bloom filter is consumed by `neoswga filter`, before optimization, because
it screens CANDIDATES against the background. Optimize scores panels from
binding positions and has no use for presence without coordinates.
"""

import pytest

from neoswga.cli._common import UNIMPLEMENTED_OPTIONS

BLOOM_OPTIONS = ["background_bloom_path", "background_sampled_path", "use_background_filter"]


@pytest.mark.parametrize("attr", BLOOM_OPTIONS)
def test_the_message_names_the_command_that_reads_the_filter(attr):
    _flag, detail = UNIMPLEMENTED_OPTIONS[attr]
    assert "neoswga filter" in detail, (
        f"{attr} tells the user their flag is unimplemented without naming the "
        f"command that does read a Bloom filter"
    )


@pytest.mark.parametrize("attr", BLOOM_OPTIONS)
def test_the_message_names_the_flag_that_enables_it(attr):
    _flag, detail = UNIMPLEMENTED_OPTIONS[attr]
    assert "--use-bloom-filter" in detail, (
        "both keys are required on that command, so naming the command alone "
        "leaves the user in the half-configured state that screens nothing"
    )


@pytest.mark.parametrize("attr", BLOOM_OPTIONS)
def test_the_named_flag_exists_on_that_command(attr):
    """Advice naming a flag that does not exist is worse than none."""
    from neoswga.cli_unified import create_parser

    parser = create_parser()
    subparsers = [
        action
        for action in parser._actions
        if hasattr(action, "choices") and isinstance(action.choices, dict)
    ][0]
    filter_parser = subparsers.choices["filter"]
    flags = {opt for action in filter_parser._actions for opt in action.option_strings}
    assert "--use-bloom-filter" in flags
    assert "--bloom-filter-path" in flags
