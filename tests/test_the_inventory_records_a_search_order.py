"""Expanding "in rank order" has to have a rank to order by.

Phase 3 of the plan for `docs/validation/pipeline_audit_2026-09-16/`, finding
F10. Not in the audit; found while verifying it.

`CandidateProvider._eligible_in_search_order` reads `step2_rank` off each
candidate's metrics and sorts by it, and its docstring says the traversal is
"Seeded by `step2_rank`, the order stage 2 established". That column is created
in `pipeline.step3`, on the 2,000-row shortlist, AFTER `record_stage2_inventory`
has already run. The inventory therefore holds it for nothing. Every candidate
fell into the same sort bucket and the order collapsed to alphabetical.

Alphabetical is not a neutral default here. Phase 4 hands the search a ranked
frontier and refills it in order, so the rank decides which of 491,836
candidates a design ever sees. Expanding an inventory in alphabetical order is a
lottery with a cost.

The rank is written for every hard-QC candidate rather than only the shortlist,
which is the whole point: a rank that exists only for the 2,000 already selected
cannot order the 489,836 that were not.

The second defect here is cost. The old traversal issued one SELECT and one JSON
parse per candidate, so reading the order was 491,836 round trips. The rank is a
column now and the database sorts.
"""

import pathlib

import pandas as pd

from neoswga.core.candidate_inventory import (
    STAGE2_INVENTORY_NAME,
    CandidateInventory,
    record_stage2_inventory,
)
from neoswga.core.candidate_provider import CandidateProvider

COND = "tm-2026-09-14:abc123"

# Deliberately in an order where alphabetical and ranked disagree: the best
# candidate by background load sorts last by sequence.
BEST = "TTTTAAAACCCC"
MIDDLE = "GGGGAAAACCCC"
WORST = "AAAAGGGGCCCC"


def _cleared():
    """Hard-QC survivors with the counts the ranking is built from."""
    return pd.DataFrame(
        {
            "primer": [WORST, MIDDLE, BEST],
            "fg_count": [100, 100, 100],
            "bg_count": [900, 300, 10],
        }
    )


def _record(tmp_path, shortlisted=()):
    """No shortlist by default, so the ranking alone decides the order.

    The shortlist leads the traversal in the order stage 2 put it in, so a test
    that handed the whole frame over as "the shortlist" would be asserting on
    the frame's row order rather than on the ranking.
    """
    cleared = _cleared()
    keep = cleared[cleared["primer"].isin(list(shortlisted))]
    return record_stage2_inventory(
        tmp_path,
        condition_id=COND,
        cleared_hard_gates=cleared,
        after_gini=cleared,
        shortlisted=keep,
        indexed=list(cleared["primer"]),
    )


def test_the_order_is_not_alphabetical(tmp_path):
    """The defect, stated as the thing that used to happen."""
    path = _record(tmp_path)

    with CandidateInventory(path) as inventory:
        order = list(inventory.iter_eligible(COND, [12]))

    assert order != sorted(order), "the traversal is still alphabetical"


def test_the_least_host_binding_candidate_comes_first(tmp_path):
    """Rank by background load per target site, the quantity stage 2 ranks on."""
    path = _record(tmp_path)

    with CandidateInventory(path) as inventory:
        order = list(inventory.iter_eligible(COND, [12]))

    assert order == [BEST, MIDDLE, WORST]


def test_every_hard_qc_candidate_gets_a_rank_not_just_the_shortlist(tmp_path):
    """The point of the change.

    A rank written only for the candidates already selected cannot order the
    ones that were not, which is exactly the set expansion draws from.
    """
    path = _record(tmp_path, shortlisted=[BEST])

    with CandidateInventory(path) as inventory:
        ranks = [inventory.metrics(s).get("search_rank") for s in (BEST, MIDDLE, WORST)]

    assert all(rank is not None for rank in ranks), "an unranked candidate is unreachable"
    assert ranks == sorted(ranks), "the rank does not follow the ranking"


def test_the_shortlist_leads_the_order(tmp_path):
    """Stage 2 already chose a working set; the search should start there.

    WORST ranks last on its own merits, so if it still comes first when
    shortlisted, the shortlist is leading rather than the raw ranking.
    """
    path = _record(tmp_path, shortlisted=[WORST])

    with CandidateInventory(path) as inventory:
        order = list(inventory.iter_eligible(COND, [12]))

    assert order[0] == WORST
    assert order[1:] == [BEST, MIDDLE]


def test_the_provider_walks_that_order(tmp_path):
    """The consumer the rank exists for."""
    path = _record(tmp_path)

    with CandidateInventory(path) as inventory:
        provider = CandidateProvider(inventory, COND, [12])
        assert provider.initial(2) == [BEST, MIDDLE]


def test_ties_break_on_the_sequence(tmp_path):
    """Equal rank must still give one order, or an unseeded run is not reproducible."""
    cleared = pd.DataFrame(
        {
            "primer": [WORST, MIDDLE, BEST],
            "fg_count": [100, 100, 100],
            "bg_count": [100, 100, 100],
        }
    )
    path = record_stage2_inventory(
        tmp_path,
        condition_id=COND,
        cleared_hard_gates=cleared,
        after_gini=cleared,
        shortlisted=cleared,
        indexed=list(cleared["primer"]),
    )

    with CandidateInventory(path) as inventory:
        order = list(inventory.iter_eligible(COND, [12]))

    assert order == sorted(order)


def test_a_candidate_with_no_counts_still_gets_an_order(tmp_path):
    """Missing counts must not drop a candidate out of the traversal.

    An unranked candidate sorts last rather than disappearing: it cleared the
    hard gates, so it is eligible, and eligibility is not the ranking's to
    revoke.
    """
    cleared = pd.DataFrame({"primer": [BEST, MIDDLE], "fg_count": [100, 0], "bg_count": [10, 5]})
    path = record_stage2_inventory(
        tmp_path,
        condition_id=COND,
        cleared_hard_gates=cleared,
        after_gini=cleared,
        shortlisted=cleared,
        indexed=list(cleared["primer"]),
    )

    with CandidateInventory(path) as inventory:
        order = list(inventory.iter_eligible(COND, [12]))

    assert set(order) == {BEST, MIDDLE}
    assert order[0] == BEST


def test_an_inventory_written_before_ranks_still_reads(tmp_path):
    """Older databases exist; they must not raise."""
    import sqlite3

    path = pathlib.Path(tmp_path) / STAGE2_INVENTORY_NAME
    connection = sqlite3.connect(str(path))
    connection.executescript(
        "CREATE TABLE candidates (sequence TEXT PRIMARY KEY, length INTEGER NOT NULL,"
        " metrics_json TEXT NOT NULL);"
        "CREATE TABLE assessments (sequence TEXT NOT NULL, condition_id TEXT NOT NULL,"
        " policy_version TEXT NOT NULL, passed INTEGER NOT NULL, reasons_json TEXT NOT NULL,"
        " metrics_json TEXT NOT NULL, PRIMARY KEY (sequence, condition_id, policy_version));"
    )
    connection.execute("INSERT INTO candidates VALUES (?, 12, '{}')", (BEST,))
    connection.execute(
        "INSERT INTO assessments VALUES (?, ?, 'qc-2026-09-15', 1, '[]', '{}')", (BEST, COND)
    )
    connection.commit()
    connection.close()

    with CandidateInventory(path) as inventory:
        assert list(inventory.iter_eligible(COND, [12], policy_version="qc-2026-09-15")) == [BEST]


def test_the_traversal_no_longer_sorts_on_a_column_it_cannot_read():
    """The docstring asserted something untrue for as long as the code existed.

    The word may still appear in prose, explaining what used to happen. What
    must not appear is a read: the inventory never held `step2_rank`, so sorting
    on it produced one bucket and an alphabetical walk.
    """
    import ast
    import inspect
    import textwrap

    from neoswga.core.candidate_provider import CandidateProvider as Provider

    tree = ast.parse(textwrap.dedent(inspect.getsource(Provider._eligible_in_search_order)))
    function = tree.body[0]
    body = function.body[1:] if ast.get_docstring(function) else function.body
    literals = {
        node.value
        for statement in body
        for node in ast.walk(statement)
        if isinstance(node, ast.Constant) and isinstance(node.value, str)
    }

    assert (
        "step2_rank" not in literals
    ), "the traversal still reads step2_rank, which the inventory never holds"
