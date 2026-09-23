"""The limitations page must not drift from what the code does.

`docs/LIMITATIONS.md` tells a user what a design run does NOT establish. A page
like that is worse than useless once it is stale: it carries the authority of a
deliberate disclosure while describing a tool that has moved on.

This repository has the defect already. `docs/SCIENCE_CITATIONS.md` states
Klenow processivity as 10,000 bp citing Bambara (1978) while the shipped
registry says 40. The prose was right when written and the code moved, and
CLAUDE.md now tells readers to treat that document as commentary rather than as
the record -- which is a reasonable response to a defect and not a fix for it.

So the quantities the page asserts are pinned against their sources. Only
numbers a reader would ACT on are pinned. Judgements are not: a test asserting
that a sentence still appears would fail on rewording and pass on becoming
wrong, which is worse than no test.
"""

import pathlib
import re

import pytest

PAGE = pathlib.Path(__file__).resolve().parent.parent / "docs" / "LIMITATIONS.md"


@pytest.fixture(scope="module")
def text():
    assert PAGE.exists(), "the READMEs link to this page"
    return PAGE.read_text(encoding="utf-8")


def test_the_two_reaches_are_the_code_s_reaches(text):
    """Conflating them inflates reported coverage severalfold, which is why
    the page names both."""
    from neoswga.core.coverage import polymerase_extension_reach

    assert polymerase_extension_reach("phi29", coverage_metric="realistic") == 3000
    assert polymerase_extension_reach("phi29", coverage_metric="processivity") == 70000

    assert "3 kb" in text and "70 kb" in text


def test_the_evidence_record_count_is_right(text):
    """The page says one of twenty-three states a checkable range."""
    import json

    records = json.loads(
        (
            pathlib.Path(__file__).resolve().parent.parent
            / "neoswga"
            / "core"
            / "registry"
            / "model_evidence.json"
        ).read_text()
    )["records"]

    with_range = sum(1 for r in records if r.get("temperature_range_c"))

    assert len(records) == 23, f"the registry now holds {len(records)} records"
    assert with_range == 1, f"{with_range} records now state a range"
    assert "twenty-three" in text and "one of" in text


def test_the_page_does_not_overstate_its_evidence(text):
    """It cites the number of measurement records as evidence of practice.

    Asserted as an upper bound rather than equality on purpose. Overstating is
    the harmful direction; understating happens every time someone adds a
    record, and a test that fails on that would train people to edit the number
    without reading the page -- which is how the page goes stale in the way
    that matters.
    """
    records = list((PAGE.parent / "validation").glob("*.md"))

    stated = int(re.search(r"keeps (\d+) measurement records", text).group(1))

    assert stated <= len(
        records
    ), f"the page claims {stated} validation records; there are {len(records)}"
    assert stated >= 30, "the claim has drifted so low it no longer says anything"


def test_every_cited_record_exists(text):
    """A limitations page whose evidence 404s is worse than one with none."""
    targets = re.findall(r"\]\((validation/[^)]+)\)", text)

    missing = [t for t in targets if not (PAGE.parent / t).exists()]

    assert targets, "the page should cite its evidence"
    assert not missing, missing


def test_the_retired_model_is_still_retired(text):
    """The page says the amplification model is off the default path. If it is
    ever restored by default, this sentence becomes a false disclaimer."""
    from neoswga.core.parameter import PipelineParameters

    assert hasattr(PipelineParameters, "__dataclass_fields__")
    assert "retired" in text.lower()

    import inspect

    from neoswga.core import pipeline

    source = inspect.getsource(pipeline)
    assert "amp_model" in source, (
        "the --amp-model opt-in has gone, so the page's description of how to "
        "restore the model is stale"
    )


def test_both_readmes_point_at_it():
    """It is only useful if a reader meets it."""
    root = PAGE.parent.parent

    assert "LIMITATIONS.md" in (root / "README.md").read_text()
    assert "LIMITATIONS.md" in (root / "docs" / "README.md").read_text()
