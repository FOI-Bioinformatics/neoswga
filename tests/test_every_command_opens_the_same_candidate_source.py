"""One rule for where a command's candidates come from.

Phase 4 increment 6 of the plan for `docs/validation/pipeline_audit_2026-09-16/`,
third piece: "`optimize` and `expand-primers` onto the same source."

`plan-pool` was put on the inventory in increment 1. `optimize` and
`expand-primers` still read `step3_df.csv` directly, so the 20,670 candidates
the Wolbachia inventory retains could not affect what either of them selected,
and the shortlist they did read is 2,000. Audit finding F1 named all three.

Three copies of "which pool am I searching" would drift, so
`candidate_source.open_source_or_list` holds the rule: the inventory when the
directory has one, the supplied list otherwise, and a line saying which it
took, because which pool a run searched should not have to be inferred.

The frontier opens at exactly the size of the list a command would have read,
which is the same choice increment 1 made for `plan-pool`. It keeps delivered
panels where they were while making the counts visible, and increment 5's refill
is what reaches past it.
"""

import logging

import pandas as pd
import pytest

from neoswga.core.candidate_inventory import record_stage2_inventory
from neoswga.core.candidate_source import (
    InventoryCandidateSource,
    ListCandidateSource,
    open_source_or_list,
)

PRIMERS = ["GCTAAAGACAAT", "TACATAACATAC", "ACGTCAGCACGA", "CAGTCAGGATCA"]
CONDITION = "tm-test:abcdef"


@pytest.fixture
def stocked(tmp_path):
    """A directory with an inventory recorded under one condition."""
    frame = pd.DataFrame(
        {"primer": PRIMERS, "fg_count": [4] * len(PRIMERS), "bg_count": [1] * len(PRIMERS)}
    )
    record_stage2_inventory(
        tmp_path,
        condition_id=CONDITION,
        cleared_hard_gates=frame,
        after_gini=frame,
        shortlisted=frame.iloc[:2],
        indexed=PRIMERS,
    )
    return tmp_path


# -- which source ----------------------------------------------------------


def test_a_directory_with_an_inventory_gives_the_inventory(stocked):
    source = open_source_or_list(stocked, CONDITION, [12], fallback=PRIMERS[:2])

    assert isinstance(source, InventoryCandidateSource)
    assert source.universe_size() == len(PRIMERS)


def test_a_directory_without_one_gives_the_list(tmp_path):
    """A run directory written before the inventory existed still works."""
    source = open_source_or_list(tmp_path, CONDITION, [12], fallback=PRIMERS[:2])

    assert isinstance(source, ListCandidateSource)
    assert source.initial() == PRIMERS[:2]


def test_an_explicit_list_wins_over_the_inventory(stocked):
    """A user who named a pool has said which pool to search."""
    source = open_source_or_list(
        stocked, CONDITION, [12], fallback=PRIMERS[:2], explicit=True
    )

    assert isinstance(source, ListCandidateSource)
    assert source.initial() == PRIMERS[:2]


def test_an_inventory_holding_nothing_for_this_reaction_gives_the_list(stocked):
    """Recorded under another chemistry is the same as not recorded.

    A fingerprint mismatch is exactly how the inventory read as empty when
    `plan-pool` first opened it, and it is worth falling back rather than
    designing over nothing.
    """
    source = open_source_or_list(
        stocked, "tm-test:a-different-reaction", [12], fallback=PRIMERS[:2]
    )

    assert isinstance(source, ListCandidateSource)


def test_the_frontier_opens_at_the_size_of_the_list_it_replaces(stocked):
    """Behaviour-preserving by construction, as in increment 1."""
    source = open_source_or_list(stocked, CONDITION, [12], fallback=PRIMERS[:2])

    assert source.initial() == list(source.frontier())
    assert len(source.frontier()) == 2
    assert source.universe_size() == 4


def test_which_source_was_taken_is_logged(stocked, caplog):
    """"Which pool did this run search" should not have to be inferred."""
    with caplog.at_level(logging.INFO):
        open_source_or_list(stocked, CONDITION, [12], fallback=PRIMERS[:2])

    assert "inventory" in caplog.text.lower()


def test_the_fallback_says_why_it_fell_back(tmp_path, caplog):
    with caplog.at_level(logging.INFO):
        open_source_or_list(tmp_path, CONDITION, [12], fallback=PRIMERS[:2])

    assert "candidate list" in caplog.text.lower()


# -- and every command uses it --------------------------------------------


def _candidate_source_callers():
    """Modules that decide where candidates come from, by AST."""
    import ast
    import pathlib

    package = pathlib.Path(__file__).resolve().parent.parent / "neoswga"
    callers = {}
    for path in sorted(package.rglob("*.py")):
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            name = getattr(node.func, "attr", None) or getattr(node.func, "id", None)
            if name in {"open_source_or_list", "open_candidate_source"}:
                callers.setdefault(path.relative_to(package).as_posix(), set()).add(name)
    return callers


def test_no_command_reads_the_candidate_csv_for_itself():
    """The ratchet. A fourth copy of the decision fails here.

    `dominating_set_optimizer.optimize` is excluded deliberately: it is the
    legacy standalone entry, documented as no longer called from the unified
    CLI and kept for external scripts, so routing it adds risk for no user.
    """
    import ast
    import pathlib

    package = pathlib.Path(__file__).resolve().parent.parent / "neoswga"
    allowed = {
        "core/dominating_set_optimizer.py": "legacy standalone entry, no CLI caller",
        "core/pipeline_qa_integration.py": "re-orders an existing step3_df in place; "
        "it does not choose which candidates a search may reach",
    }
    offenders = {}
    for path in sorted(package.rglob("*.py")):
        relative = path.relative_to(package).as_posix()
        if relative in allowed:
            continue
        source = path.read_text()
        if 'step3_df["primer"]' in source or "step3_df['primer']" in source:
            offenders[relative] = True

    assert not offenders, (
        "These modules pick their own candidate pool out of step3_df.csv instead "
        f"of asking open_source_or_list: {sorted(offenders)}. The inventory holds "
        "every candidate that cleared hard QC and the CSV holds the shortlist, so "
        "reading the CSV directly makes the rest unreachable (audit finding F1)."
    )


def test_the_allowed_list_has_no_stale_entries():
    """It can only shrink, like the other allowlists in this suite."""
    import pathlib

    package = pathlib.Path(__file__).resolve().parent.parent / "neoswga"
    allowed = ["core/dominating_set_optimizer.py", "core/pipeline_qa_integration.py"]
    for relative in allowed:
        source = (package / relative).read_text()
        assert 'step3_df["primer"]' in source or "step3_df[primer_column" in source, (
            f"{relative} no longer reads the candidate CSV. Remove its entry from "
            "the allowed list in this test; a fix is not finished until it drops "
            "its excuse."
        )
