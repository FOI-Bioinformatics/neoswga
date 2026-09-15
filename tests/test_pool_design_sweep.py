"""Comparing chemistries means redesigning the pool, not rescoring one panel.

Task 7 of the condition-aware pool design plan.

`sweep_conditions` scores ONE fixed panel under several reactions. That answers
"how does this panel behave if I change the buffer", which is a useful question
and not the one a design comparison asks. A panel chosen under one chemistry
carries that chemistry's decisions: which candidates passed the Tm window, which
the ranking favoured, which the dimer screen removed. Rescoring it elsewhere
measures the transplant, not the alternative.

`design_sweep` runs a separate design per condition and per oligo length, from
the shared candidate inventory, and compares the results. The two are kept
distinct because conflating them is how a chemistry can appear better simply by
suiting the panel that was already chosen.

Eligibility is per condition, so a candidate admitted under one reaction and not
another appears in one design and not the other. That is the point rather than
an inconsistency.
"""

import pytest

from neoswga.core.candidate_inventory import CandidateInventory
from neoswga.core.pool_design_sweep import design_sweep, nondominated
from neoswga.core.pool_objective import PoolConstraints

A, C, G = "AAAAAAAAAAAA", "CCCCCCCCCCCC", "GGGGGGGGGGGG"
SHORT = "ACGTACGTAC"


class _Condition:
    def __init__(self, name):
        self.name = name

    def fingerprint(self):
        return f"cond:{self.name}"


@pytest.fixture
def inventory(tmp_path):
    """A candidate eligible under one condition only, and one short oligo."""
    with CandidateInventory(tmp_path / "inv.sqlite") as inv:
        for seq in (A, C, G):
            inv.record_candidate(seq, {"step2_rank": 0})
        inv.record_candidate(SHORT, {"step2_rank": 0})
        for seq in (A, C):
            inv.record_assessment(seq, "cond:warm", passed=True, reasons=[], metrics={})
            inv.record_assessment(seq, "cond:cold", passed=True, reasons=[], metrics={})
        # Admitted only under the warmer reaction.
        inv.record_assessment(G, "cond:warm", passed=True, reasons=[], metrics={})
        inv.record_assessment(G, "cond:cold", passed=False, reasons=["tm_window"], metrics={})
        inv.record_assessment(SHORT, "cond:warm", passed=True, reasons=[], metrics={})
        inv.commit()
    with CandidateInventory(tmp_path / "inv.sqlite") as inv:
        yield inv


def _recorder(seen):
    def run_design(**kwargs):
        provider = kwargs["provider"]
        available = provider.initial(limit=100)
        seen.append((kwargs["condition"].fingerprint(), kwargs["length"], tuple(available)))
        return {
            "rows": [],
            "recommendations": [],
            "eligible_candidates": len(available),
        }

    return run_design


def test_each_condition_designs_from_its_own_eligible_set(inventory):
    seen = []
    design_sweep(
        inventory,
        conditions=[_Condition("warm"), _Condition("cold")],
        lengths=[12],
        sizes=[2],
        coverage_targets=[0.9],
        constraints=PoolConstraints(),
        run_design=_recorder(seen),
    )

    warm = next(s for s in seen if s[0] == "cond:warm")
    cold = next(s for s in seen if s[0] == "cond:cold")
    assert G in warm[2], "a candidate admitted by the warmer reaction was not offered"
    assert G not in cold[2], "a candidate the colder reaction rejects was offered anyway"


def test_every_condition_and_length_appears_in_the_sweep(inventory):
    seen = []
    out = design_sweep(
        inventory,
        conditions=[_Condition("warm"), _Condition("cold")],
        lengths=[10, 12],
        sizes=[2],
        coverage_targets=[0.9],
        constraints=PoolConstraints(),
        run_design=_recorder(seen),
    )

    assert {(d["condition"], d["length"]) for d in out["designs"]} == {
        ("cond:warm", 10),
        ("cond:warm", 12),
        ("cond:cold", 10),
        ("cond:cold", 12),
    }


def test_a_length_this_reaction_admits_nothing_at_is_a_result_not_an_error(inventory):
    """Counted but ineligible differs from never counted.

    The 10-mer was enumerated and clears the warm reaction's gates but not the
    cold one's. That is a finding about the chemistry, so it is recorded as a
    design with no eligible candidates rather than raised, and no design is run.
    """
    seen = []
    out = design_sweep(
        inventory,
        conditions=[_Condition("cold")],
        lengths=[10],
        sizes=[2],
        coverage_targets=[0.9],
        constraints=PoolConstraints(),
        run_design=_recorder(seen),
    )

    assert seen == [], "a design was attempted with nothing to design from"
    entry = out["designs"][0]
    assert entry["eligible_candidates"] == 0
    assert entry["result"] is None
    assert "hard gates" in entry["note"]


def test_a_length_nothing_was_counted_at_says_what_to_run(inventory):
    """An empty frontier that looks like a result is worse than a refusal."""
    with pytest.raises(ValueError, match="count-kmers") as excinfo:
        design_sweep(
            inventory,
            conditions=[_Condition("warm")],
            lengths=[18],
            sizes=[2],
            coverage_targets=[0.9],
            constraints=PoolConstraints(),
            run_design=_recorder([]),
        )

    assert "18" in str(excinfo.value)


def test_the_sweep_records_which_design_each_result_came_from(inventory):
    out = design_sweep(
        inventory,
        conditions=[_Condition("warm")],
        lengths=[12],
        sizes=[2],
        coverage_targets=[0.9],
        constraints=PoolConstraints(),
        run_design=_recorder([]),
    )

    assert out["designs"][0]["condition"] == "cond:warm"
    assert out["designs"][0]["length"] == 12
    assert out["designs"][0]["eligible_candidates"] == 3


def test_the_frontier_keeps_only_nondominated_designs():
    """Dominated means worse or equal on every axis, and worse on one."""
    designs = [
        {"label": "small-good", "size": 4, "coverage": 0.90, "background": 10},
        {"label": "big-worse", "size": 8, "coverage": 0.85, "background": 20},
        {"label": "big-better", "size": 8, "coverage": 0.95, "background": 5},
    ]

    kept = {d["label"] for d in nondominated(designs)}

    assert "small-good" in kept, "fewer oligos at high coverage is not dominated"
    assert "big-better" in kept
    assert "big-worse" not in kept


def test_an_equal_design_does_not_dominate_its_twin():
    designs = [
        {"label": "a", "size": 4, "coverage": 0.9, "background": 10},
        {"label": "b", "size": 4, "coverage": 0.9, "background": 10},
    ]

    assert len(nondominated(designs)) == 2
