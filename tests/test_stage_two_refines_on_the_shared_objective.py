"""Stage 2 must refine on the quantity the row is judged on.

Phase 2 of the plan for `docs/validation/pipeline_audit_2026-09-16/`, the last
part of audit finding F2.

Task 5 of the earlier plan gave search and acceptance one definition, and Task 6
used it in the greedy. `refine_hybrid_stage2` was missed. It still calls
`refine_by_swaps` with the raw-bin rule -- covered bases first, then background
load -- while `plan_pool` accepts the result on occupancy-weighted coverage and a
selectivity floor. Two rules, one panel, and the stage that picks it used the one
the report does not quote.

The two disagree in a known direction: a primer with many sites and a melting
temperature well below the reaction temperature touches many bins and
contributes little amplification, so the bin rule prefers it and the accepted
metric does not.

The objective rides on the optimizer rather than through `optimize()`'s
signature. `refine_hybrid_stage2` already takes the optimizer as its first
argument, so nothing public changes, and a plain `optimize` run -- which has no
`PoolObjective` and no constraints -- keeps the rule it has always used. Only the
path that declares an objective refines on one.
"""

from types import SimpleNamespace

import pytest

from neoswga.core.base_optimizer import OptimizerConfig
from neoswga.core.pool_objective import PoolConstraints, PoolObjective
from neoswga.core.swap_refinement import refine_hybrid_stage2

WIDE = "AAAACCCCGGGG"  # many bins, poorly bound
NARROW = "CCCCGGGGTTTT"  # fewer bins, well bound
SPARE = "GGGGTTTTAAAA"


class _Optimizer:
    """Only what `refine_hybrid_stage2` reads off its first argument."""

    def __init__(self, regions, background=None, objective=None):
        self._regions = regions
        self.background_pruning = bool(background)
        self.bg_prefixes = ["bg"] if background else []
        self._background = background or {}
        self.max_dimer_bp = 3
        self.swap_max_evaluations = 5000
        self.swap_max_seconds = 5.0
        self.config = OptimizerConfig(max_dimer_bp=3, max_self_dimer_bp=4)
        if objective is not None:
            self.pool_objective = objective

    def _coverage_bins_by_primer(self, pool):
        return {p: self._regions.get(p, set()) for p in pool}

    def _bin_key(self, region):
        return region

    def _count_background_sites(self, primers):
        return sum(self._background.get(p, 0) for p in primers)


def _regions():
    """WIDE touches four bins, NARROW two. The bin rule prefers WIDE."""
    return {
        WIDE: {_Bin(0, 1000), _Bin(1000, 2000), _Bin(2000, 3000), _Bin(3000, 4000)},
        NARROW: {_Bin(0, 1000), _Bin(1000, 2000)},
        SPARE: {_Bin(4000, 5000)},
    }


class _Bin:
    """A coverage region with the width the bin rule weighs it by."""

    def __init__(self, start, end):
        self.start, self.end = start, end

    def __hash__(self):
        return hash((self.start, self.end))

    def __eq__(self, other):
        return (self.start, self.end) == (other.start, other.end)


def _objective(coverages, density=90.0, floor=None):
    """Occupancy-weighted coverage from a table, against an optional floor."""

    def evaluate(primers):
        key = tuple(sorted(primers))
        return SimpleNamespace(
            fg_coverage=coverages[key],
            effective_fg_coverage=coverages[key],
            selectivity_density=density if isinstance(density, float) else density[key],
            total_bg_sites=0,
        )

    return PoolObjective(evaluate, PoolConstraints(min_selectivity_density=floor))


def test_without_an_objective_the_bin_rule_is_unchanged():
    """Plain `optimize` has no objective and must keep the rule it had."""
    opt = _Optimizer(_regions())

    result = refine_hybrid_stage2(opt, [SPARE], [SPARE, WIDE, NARROW], fixed_primers=[])

    assert set(result) == {WIDE}, "the bin rule should take the widest primer"


def test_with_an_objective_the_accepted_metric_decides():
    """The case the bin rule gets wrong.

    WIDE covers twice the bins and half the occupancy-weighted coverage. The
    stage that chooses the panel should choose on the second, because that is
    what the row is accepted on.
    """
    coverages = {
        (SPARE,): 0.05,
        (WIDE,): 0.10,
        (NARROW,): 0.30,
        (NARROW, WIDE): 0.35,
        (SPARE, WIDE): 0.12,
        (NARROW, SPARE): 0.31,
    }
    opt = _Optimizer(_regions(), objective=_objective(coverages))

    result = refine_hybrid_stage2(opt, [SPARE], [SPARE, WIDE, NARROW], fixed_primers=[])

    assert set(result) == {NARROW}, "stage 2 refined on bins, not on the accepted metric"


def test_the_objective_keeps_a_feasible_panel_feasible():
    """A constraint is not a scoring term stage 2 may outbid."""
    coverages = {
        (SPARE,): 0.20,
        (WIDE,): 0.99,
        (NARROW,): 0.25,
        (NARROW, WIDE): 0.99,
        (SPARE, WIDE): 0.99,
        (NARROW, SPARE): 0.30,
    }
    density = {
        (SPARE,): 90.0,
        (WIDE,): 1.0,
        (NARROW,): 90.0,
        (NARROW, WIDE): 1.0,
        (SPARE, WIDE): 1.0,
        (NARROW, SPARE): 90.0,
    }
    objective = _objective(coverages, density=density, floor=10.0)
    opt = _Optimizer(_regions(), objective=objective)

    result = refine_hybrid_stage2(opt, [SPARE], [SPARE, WIDE, NARROW], fixed_primers=[])

    assert WIDE not in result, "stage 2 took a violating primer for the coverage"
    assert not objective.violations(result)


def test_a_fixed_primer_survives_either_rule():
    """The iterative workflow depends on validated oligos being kept."""
    coverages = {
        (SPARE,): 0.05,
        (NARROW, SPARE): 0.31,
        (SPARE, WIDE): 0.12,
    }
    opt = _Optimizer(_regions(), objective=_objective(coverages))

    result = refine_hybrid_stage2(opt, [SPARE], [SPARE, WIDE, NARROW], fixed_primers=[SPARE])

    assert SPARE in result


def test_plan_pool_attaches_the_objective_through_the_wrapper():
    """A fast guard, and no longer a misleading one.

    This used to assert by AST that `plan_pool` contained an assignment to an
    attribute named `pool_objective`. It did, and Stage 2 still received None,
    because the optimizer `plan_pool` is handed is a wrapper and the stage that
    reads the attribute is the delegate inside it. Asserting the assignment
    exists says nothing about where it lands.

    So this now checks that the attachment goes through
    `attach_search_config`, which sets it on both. The path itself is asserted
    end to end, through a real factory-built optimizer, in
    `tests/test_the_objective_reaches_the_stage_that_refines.py`; this remains
    only as a cheap signal that the mechanism was not swapped back.
    """
    import ast
    import inspect
    import textwrap

    from neoswga.core import pool_planner

    source = textwrap.dedent(inspect.getsource(pool_planner.plan_pool))
    attached = {
        node.args[1].value
        for node in ast.walk(ast.parse(source))
        if isinstance(node, ast.Call)
        and getattr(node.func, "id", None) == "attach_search_config"
        and len(node.args) >= 2
        and isinstance(node.args[1], ast.Constant)
    }

    assert "pool_objective" in attached, (
        "plan_pool no longer attaches its objective through attach_search_config. "
        "A plain attribute assignment reaches the wrapper only, which is how "
        "stage 2 came to refine on raw bins while the row was accepted on "
        "occupancy-weighted coverage."
    )


@pytest.mark.parametrize("attribute", ["pool_objective"])
def test_the_refinement_reads_the_attribute_it_is_given(attribute):
    """Names the contract, so renaming one end without the other fails here."""
    import inspect

    from neoswga.core import swap_refinement

    assert attribute in inspect.getsource(swap_refinement.refine_hybrid_stage2)
