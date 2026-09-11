"""A capability nobody can find is worth little more than one that does not exist.

--auto-size and --show-frontier have been wired since before this audit and
appear in no user-facing document, which is why this project's own GC-tier sweep
was run by hand with repeated -n values.

CLAUDE.md was deliberately absent from DOCS while the plan that added these
flags was running: it carried several hundred lines of the maintainer's
uncommitted work and no plan task was allowed to touch it, so the section it
should gain was collected for hand-over. That hand-over was applied on
2026-09-11 and CLAUDE.md now carries a "Choosing the set size" subsection, so it
joins DOCS and is held to the same three checks as the others.
"""

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent

DOCS = ["README.md", "docs/params-reference.md", "CLAUDE.md"]

# Kept as a separate name because the absence check below once covered a file
# DOCS did not. The two lists are the same now; the distinction costs nothing
# and the absence check is the one that must never narrow.
ALL_DOCS = list(DOCS)


@pytest.mark.parametrize("doc", DOCS)
def test_auto_size_is_documented(doc):
    assert "--auto-size" in (ROOT / doc).read_text(), f"{doc} does not mention --auto-size"


@pytest.mark.parametrize("doc", DOCS)
def test_show_frontier_is_documented(doc):
    assert "--show-frontier" in (ROOT / doc).read_text(), f"{doc} does not mention --show-frontier"


@pytest.mark.parametrize("doc", ALL_DOCS)
def test_no_document_invents_a_pareto_flag(doc):
    """The audit calls it `--pareto`; argparse does not. Documenting a name
    that does not parse is worse than documenting nothing."""
    assert "--pareto" not in (ROOT / doc).read_text(), f"{doc} names a --pareto flag"


@pytest.mark.parametrize("doc", DOCS)
def test_the_docs_state_the_twenty_primer_ceiling(doc):
    """Both flags stop at 20 primers. A reader who does not know that will
    reach for them on a 96- or 160-oligo panel and get a silently clamped
    answer."""
    text = (ROOT / doc).read_text()
    assert "20" in text, f"{doc} does not state the 20-primer ceiling"


def test_the_docs_say_auto_size_does_not_weigh_specificity():
    """recommend_set_size inverts a coverage curve and never reads the pool.
    Describing it as choosing the right size would be wrong."""
    text = (ROOT / "docs" / "params-reference.md").read_text()
    assert "candidate pool" in text.lower()


def test_the_marginal_table_is_documented():
    assert (
        "pp/primer" in (ROOT / "docs" / "params-reference.md").read_text()
    ), "the marginal coverage table optimize now prints is undocumented"


def test_the_documented_flags_actually_parse():
    """The one check that catches a doc drifting away from the parser."""
    from neoswga import cli_unified

    parser = cli_unified.create_parser()
    for argv in (
        ["optimize", "-j", "params.json", "--auto-size"],
        ["optimize", "-j", "params.json", "--show-frontier"],
        ["design", "-j", "params.json", "--auto-size"],
        ["design", "-j", "params.json", "--show-frontier"],
    ):
        parser.parse_args(argv)
