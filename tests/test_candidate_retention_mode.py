"""`max_primer` orders a search; it must not decide what can ever be selected.

Task 3 of the condition-aware pool design plan, second half.

The first half recorded every candidate clearing the declared hard gates. This
half makes them genuinely addressable, which needs one more thing: background
positions. Foreground positions were already scanned for every hard-QC survivor,
because the Gini gate reads them. Background positions were scanned only for the
shortlist, so a candidate the ranking cut had no background index, and a design
that later reached for it would score it against an empty background -- the
silent-zero shape this project has met repeatedly.

`candidate_retention` selects the policy:

``all_qc`` (default for new designs)
    Index every candidate that cleared the hard gates. `max_primer` still
    chooses the working shortlist written to `step2_df.csv`, so the optimizer's
    runtime is unchanged, but nothing the ranking set aside is unreachable.

``legacy``
    Index the shortlist only, reproducing the historical behaviour for
    comparison against an existing result.
"""

import pytest

from neoswga.core.candidate_inventory import CandidateInventory, record_stage2_inventory

A, C, G = "AAAAAAAAAAAA", "CCCCCCCCCCCC", "GGGGGGGGGGGG"
COND = "tm-2026-09-14:abc123"


def test_the_retention_mode_is_declared_in_the_schema():
    import json
    import pathlib

    schema = json.loads(pathlib.Path("neoswga/core/schema/params.schema.json").read_text())
    entry = schema["properties"]["candidate_retention"]
    assert entry["default"] == "all_qc"
    assert set(entry["enum"]) == {"all_qc", "legacy"}


def test_the_mode_binds_a_parameter_global():
    """Otherwise it joins the inert-key class the ratchet exists to prevent."""
    from neoswga.core import parameter

    assert hasattr(parameter, "candidate_retention")


def test_all_qc_indexes_every_hard_qc_survivor():
    """The pool handed to the background scan, by mode."""
    from neoswga.core.candidate_inventory import background_scan_pool as _background_scan_pool

    cleared = [A, C, G]
    shortlist = [A]

    assert sorted(_background_scan_pool(cleared, shortlist, "all_qc")) == sorted(cleared)
    assert list(_background_scan_pool(cleared, shortlist, "legacy")) == shortlist


def test_an_unknown_mode_is_refused_rather_than_guessed():
    from neoswga.core.candidate_inventory import background_scan_pool as _background_scan_pool

    with pytest.raises(ValueError, match="candidate_retention"):
        _background_scan_pool([A], [A], "whatever")


def test_the_shortlist_is_still_what_the_optimizer_reads(tmp_path):
    """Retention is not the same as search budget.

    Indexing everything must not change the pool the optimizer works over, or
    a design's runtime would move with a retention setting rather than with a
    search setting. `step2_df.csv` stays the shortlist.
    """
    import pandas as pd

    cleared = pd.DataFrame({"primer": [A, C, G], "fg_count": [9, 5, 1]})
    shortlisted = cleared[cleared["primer"] == A]

    path = record_stage2_inventory(
        tmp_path,
        condition_id=COND,
        cleared_hard_gates=cleared,
        after_gini=cleared,
        shortlisted=shortlisted,
    )

    with CandidateInventory(path) as inv:
        assert sorted(inv.iter_eligible(COND, [12])) == sorted([A, C, G])
        assert [s for s in (A, C, G) if inv.metrics(s)["shortlisted"]] == [A]


def test_the_counts_distinguish_enumerated_from_examined(tmp_path):
    """Counted, assessed, hard-QC passed and shortlisted are four numbers.

    The funnel reported one of them and labelled it as if it were the pool.
    """
    import pandas as pd

    cleared = pd.DataFrame({"primer": [A, C, G]})
    path = record_stage2_inventory(
        tmp_path,
        condition_id=COND,
        cleared_hard_gates=cleared,
        after_gini=cleared[cleared["primer"] != G],
        shortlisted=cleared[cleared["primer"] == A],
    )

    with CandidateInventory(path) as inv:
        counts = inv.counts()
        assert counts["candidates"] == 3
        assert counts["assessed"] == 3
        assert counts["hard_qc_passed"] == 3
        assert sum(1 for s in (A, C, G) if inv.metrics(s)["shortlisted"]) == 1
