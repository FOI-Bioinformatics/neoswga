"""The conceptual walkthrough's numbers must be the code's numbers.

`docs/guides/why-four-steps.md` explains what each pipeline step is for. Prose
like that goes stale in a particular way: the explanation stays broadly true
while a constant inside it quietly stops matching, and nobody notices because
nothing executes a document.

This repository has been bitten by exactly that. `docs/SCIENCE_CITATIONS.md`
states Klenow processivity as 10,000 bp citing Bambara (1978) while the shipped
registry says 40 bp -- the prose was right when written and the code moved.
CLAUDE.md now tells readers to treat that document as commentary rather than as
the record, which is a reasonable response to a defect and not a fix for it.

So the walkthrough's load-bearing constants are pinned here against their
sources. Only quantities a reader would ACT on are pinned: the two reaches,
because designing to the wrong one inflates coverage severalfold; the evenness
site threshold, because a small-target design must lower it; the per-length
frequency rescaling, because it is the step's least obvious behaviour; and the
retired stage name, because an alias was deliberately not provided.

Prose claims that are judgements rather than constants are NOT pinned. A test
asserting that a sentence still appears would fail on rewording and pass on
becoming wrong, which is worse than no test.
"""

import pathlib
import re

import pytest

DOC = pathlib.Path(__file__).resolve().parent.parent / "docs" / "guides" / "why-four-steps.md"


@pytest.fixture(scope="module")
def text():
    assert DOC.exists(), f"{DOC} is missing; the docs index links to it"
    return DOC.read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# The two reaches
# ---------------------------------------------------------------------------


def test_the_selection_reach_is_the_one_the_document_quotes(text):
    """~3 kb for phi29. Designing to the other one inflates coverage 5-20x."""
    from neoswga.core.coverage import polymerase_extension_reach

    realistic = polymerase_extension_reach("phi29", coverage_metric="realistic")

    assert realistic == 3000, realistic
    assert "~3 kb" in text


def test_the_connectivity_reach_is_the_one_the_document_quotes(text):
    from neoswga.core.coverage import polymerase_extension_reach

    processivity = polymerase_extension_reach("phi29", coverage_metric="processivity")

    assert processivity == 70000, processivity
    assert "~70 kb" in text


def test_the_document_keeps_the_two_reaches_apart(text):
    """The whole point of the table. Conflating them is the reported defect."""
    assert "selecting and scoring coverage" in text
    assert "connectivity" in text


# ---------------------------------------------------------------------------
# Defaults a reader is told to change
# ---------------------------------------------------------------------------


def test_the_evenness_site_threshold_matches_the_default(text):
    """A small-target design must lower this, so the number has to be right."""
    from neoswga.core.primer_attributes import DEFAULT_MIN_GINI_SITES

    assert DEFAULT_MIN_GINI_SITES == 3, DEFAULT_MIN_GINI_SITES
    assert (
        "default\n3)" in text
        or "(default\n3)" in text
        or "default 3" in text
        or ("`min_gini_sites` (default" in text and "3" in text)
    )


def test_the_frequency_rescaling_exponent_matches_the_code(text):
    """`4**(10-k)`: the least obvious thing step 2 does."""
    from neoswga.core.filter import _scale_freq_threshold

    # Ten is the reference length, so the threshold is unchanged there and
    # divides by four for each base beyond it.
    assert _scale_freq_threshold(1e-5, 10) == pytest.approx(1e-5)
    assert _scale_freq_threshold(1e-5, 11) == pytest.approx(1e-5 / 4)
    assert _scale_freq_threshold(1e-5, 12) == pytest.approx(1e-5 / 16)

    assert "4**(10-k)" in text


# ---------------------------------------------------------------------------
# The renamed stage
# ---------------------------------------------------------------------------


def test_the_document_names_the_stage_the_cli_actually_has(text):
    from neoswga.cli_unified import create_parser

    parser = create_parser()
    subparsers = next(a for a in parser._actions if a.__class__.__name__ == "_SubParsersAction")

    assert "prepare-candidates" in subparsers.choices
    assert "score" not in subparsers.choices, (
        "an alias was deliberately not provided, so the walkthrough's claim "
        "that the old name fails is now wrong"
    )
    assert "prepare-candidates" in text


# ---------------------------------------------------------------------------
# Every link resolves
# ---------------------------------------------------------------------------


def test_every_relative_link_resolves(text):
    """A walkthrough is mostly a set of pointers; a dead one wastes a reader."""
    targets = re.findall(r"\]\((?!https?:)([^)#]+)", text)

    missing = [t for t in targets if not (DOC.parent / t).exists()]

    assert not missing, missing
    assert targets, "the document should point somewhere"
