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

import json

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


def test_a_design_grid_inherits_the_baseline_before_applying_overrides(tmp_path):
    """An override names what CHANGES, not the whole reaction.

    A grid entry listing only `dmso_percent` must keep the resolved buffer,
    salts and concentration from the run's own configuration. Rebuilding a
    condition from the overrides alone would silently compare designs against
    library defaults rather than against the user's reaction.
    """
    from neoswga.core.pool_design_sweep import load_design_grid

    baseline = {"polymerase": "phi29", "temp": 30.0, "na_conc": 75.0, "mg_conc": 10.0}
    grid = {"lengths": [10, 12], "conditions": [{}, {"dmso_percent": 8.0}]}

    lengths, conditions = load_design_grid(grid, baseline)

    assert lengths == [10, 12]
    assert conditions[0].na_conc == 75.0 and conditions[0].dmso_percent == 0.0
    assert conditions[1].na_conc == 75.0, "an override dropped the resolved buffer"
    assert conditions[1].mg_conc == 10.0
    assert conditions[1].dmso_percent == 8.0


def test_a_grid_override_is_validated_by_the_condition_model(tmp_path):
    from neoswga.core.pool_design_sweep import load_design_grid

    with pytest.raises(Exception):
        load_design_grid(
            {"lengths": [12], "conditions": [{"dmso_percent": 500.0}]},
            {"polymerase": "phi29", "temp": 30.0},
        )


def test_a_grid_must_name_at_least_one_length_and_condition():
    from neoswga.core.pool_design_sweep import load_design_grid

    with pytest.raises(ValueError, match="lengths"):
        load_design_grid({"conditions": [{}]}, {"polymerase": "phi29", "temp": 30.0})
    with pytest.raises(ValueError, match="conditions"):
        load_design_grid({"lengths": [12]}, {"polymerase": "phi29", "temp": 30.0})


def test_the_frontier_names_the_smallest_qualifying_pool_per_target():
    from neoswga.core.pool_design_sweep import frontier

    designs = [
        {
            "condition": "cond:warm",
            "length": 12,
            "result": {
                "rows": [
                    {"size": 4, "eligible": True, "coverage": 0.92, "background_sites": 30},
                    {"size": 8, "eligible": True, "coverage": 0.96, "background_sites": 90},
                ]
            },
        },
        {
            "condition": "cond:cold",
            "length": 12,
            "result": {
                "rows": [
                    {"size": 6, "eligible": True, "coverage": 0.93, "background_sites": 10},
                ]
            },
        },
    ]

    out = frontier(designs, coverage_targets=[0.90, 0.95])

    at_90 = out["by_target"][0.90]
    assert {(e["condition"], e["size"]) for e in at_90} == {("cond:warm", 4), ("cond:cold", 6)}
    at_95 = out["by_target"][0.95]
    assert [(e["condition"], e["size"]) for e in at_95] == [("cond:warm", 8)]


def test_the_frontier_omits_designs_that_never_reach_a_target():
    from neoswga.core.pool_design_sweep import frontier

    designs = [
        {
            "condition": "cond:weak",
            "length": 12,
            "result": {
                "rows": [{"size": 4, "eligible": True, "coverage": 0.40, "background_sites": 1}]
            },
        }
    ]

    out = frontier(designs, coverage_targets=[0.90])

    assert out["by_target"][0.90] == []
    assert "cond:weak" in out["unreached"][0.90]


def test_an_ineligible_panel_cannot_qualify():
    """Failing a constraint is not a coverage result."""
    from neoswga.core.pool_design_sweep import frontier

    designs = [
        {
            "condition": "cond:warm",
            "length": 12,
            "result": {
                "rows": [{"size": 4, "eligible": False, "coverage": 0.99, "background_sites": 1}]
            },
        }
    ]

    assert frontier(designs, coverage_targets=[0.90])["by_target"][0.90] == []


def test_plan_pool_exposes_the_design_grid():
    """The grid has to be reachable, or it is a module nothing calls."""
    from neoswga.cli_unified import create_parser

    parser = create_parser()
    action = next(
        a
        for sub in parser._subparsers._group_actions
        for name, sp in sub.choices.items()
        if name == "plan-pool"
        for a in sp._actions
        if a.dest == "design_grid"
    )
    assert action.default is None, "an absent grid must not be a real default"


def test_a_grid_file_is_read_and_resolved(tmp_path):
    from neoswga.cli.plan_pool import load_grid_file

    path = tmp_path / "grid.json"
    path.write_text(json.dumps({"lengths": [12], "conditions": [{}, {"dmso_percent": 5.0}]}))

    lengths, conditions = load_grid_file(
        str(path), {"polymerase": "phi29", "temp": 30.0, "mg_conc": 10.0}
    )

    assert lengths == [12]
    assert [c.dmso_percent for c in conditions] == [0.0, 5.0]
    assert all(c.mg_conc == 10.0 for c in conditions), "the baseline buffer was dropped"


def test_a_grid_file_that_is_not_a_grid_says_so(tmp_path):
    from neoswga.cli.plan_pool import load_grid_file

    path = tmp_path / "grid.json"
    path.write_text(json.dumps({"something": "else"}))

    with pytest.raises(ValueError, match="lengths"):
        load_grid_file(str(path), {"polymerase": "phi29", "temp": 30.0})
