"""Step 4 reads `background_aggregate`: per-host figures, and limits that bind.

Phase 7c. 7a plumbed the key and 7b gave it to the step-2 gate; this is the
step-4 half, and it replaces 7a's non-effect test
(`test_background_aggregate_moves_no_step4_figure_before_7c`) with the positive
assertion.

Four claims, each the one a half-wired dispatch would break silently.

**Nothing moves under `sum`.** The pooled figures are compared against
constants RECORDED ON THE PRE-CHANGE TREE (`63875a1`, the tip of 7b) rather
than against a second run of the new code, which would only prove the new path
agrees with itself. `scripts`-free: the script that produced them is in the
report for this phase, and every value below was printed by it under both trees
and diffed.

**`worst-case` on one host is `sum` on one host, exactly.** Not nearly: the
whole `to_dict()` is compared. `_coverage_over_prefixes` returns
`(fraction * length) / length`, which is not bit-identical to `fraction`, so
this holds because `per_host_metrics` carries the pooled figures over for a
single host rather than recomputing them.

**A limit binds per host.** A panel clean on a large host and dense on a small
one passes `min_selectivity_density` pooled and fails it per host, with no new
limit key: `pool_objective._LIMITS` is unchanged and reads its value through
`per_host_metrics.limit_figure`.

**Unknown is refused, not satisfied.** A host with no position data makes the
worst host unknown, so the reduction is `None` and `_unmet` fails the limit. A
worst case over the hosts that happened to answer is not a worst case.
"""

from __future__ import annotations

import json
import random

import h5py
import numpy as np
import pytest

from neoswga.core import parameter, per_host_metrics
from neoswga.core.background_aggregate import SUM, WORST_CASE
from neoswga.core.base_optimizer import OptimizerConfig
from neoswga.core.dominating_set_adapter import DominatingSetAdapter
from neoswga.core.position_cache import PositionCache
from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.selectivity import MAX_SELECTIVITY
from neoswga.core.utility import reverse_complement

K = 12
LENGTHS = {"fg": 40_000, "small": 40_000, "large": 1_720_000}

#: Every pooled step-4 figure for the two-host world below, as `63875a1`
#: printed it. The large host carries no site and the small one carries every
#: site, so pooling dilutes the small host to 2.3 percent of the denominator --
#: which is why the pooled density is 44.0 and the worst host's is 1.0.
POOLED_TWO_HOST = {
    "bg_coverage": 0.02159715909090909,
    "bg_total_length": 1760000,
    "coverage_uniformity": 0.11105,
    "effective_bg_sites": 8.0,
    "effective_fg_coverage": 0.950275,
    "effective_fg_sites": 8.0,
    "extension_reach": 3000,
    "fg_coverage": 0.950275,
    "fg_total_length": 40000,
    "gap_entropy": 0.5032583347756457,
    "gap_gini": 0.11105,
    "headline_coverage": 1.0,
    "headline_reach": 10000,
    "max_gap": 5000,
    "mean_gap": 4444.444444444444,
    "selectivity_density": 44.0,
    "selectivity_mode": "exact",
    "selectivity_ratio": 1.0,
    "total_bg_sites": 8,
    "total_fg_sites": 8,
}

#: The same world with the small host alone, from the same tree. The pooled
#: density is 1.0 here, which is what the two-host worst case reduces to.
POOLED_ONE_HOST = dict(
    POOLED_TWO_HOST,
    bg_coverage=0.950275,
    bg_total_length=40000,
    selectivity_density=1.0,
)


# ----------------------------------------------------------------------
# The world
# ----------------------------------------------------------------------


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    """One 40 kb target, a 40 kb host and a 1.72 Mb host, 43 times its size.

    Every panel oligo binds the small host eight times and the large host
    nowhere, which is the configuration a pooled background hides. `absent` is
    a fourth reference with NO dataset for the panel's oligos: scanned-and-empty
    and never-scanned must not read the same, and this is the second one.
    """
    directory = tmp_path_factory.mktemp("per_host_world")
    rng = random.Random(7001)
    primers = sorted({"".join(rng.choice("ACGT") for _ in range(K)) for _ in range(8)})

    prefixes = {}
    for name, length in LENGTHS.items():
        prefix = str(directory / name)
        prefixes[name] = prefix
        step = length // 8
        with h5py.File(f"{prefix}_{K}mer_positions.h5", "w") as handle:
            for primer in primers:
                hits = [] if name == "large" else [i * step + 11 for i in range(8)]
                handle.create_dataset(primer, data=np.array(hits, dtype=np.int64))
                handle.create_dataset(reverse_complement(primer), data=np.array([], dtype=np.int64))

    # Indexed, but not for THESE oligos: one unrelated key, so the file exists
    # and the prefix is loadable while every panel member is unanswered.
    absent = str(directory / "absent")
    prefixes["absent"] = absent
    with h5py.File(f"{absent}_{K}mer_positions.h5", "w") as handle:
        handle.create_dataset("A" * K, data=np.array([7], dtype=np.int64))

    return {"primers": primers, **prefixes}


def _metrics(world, monkeypatch, aggregate, hosts=("small", "large"), explicit=False):
    """`compute_metrics` for the panel, under one reading of the background.

    `explicit` passes the resolved value as the constructor argument instead of
    leaving it to the `parameter` globals, which is the route
    `unified_optimizer` takes. Both must give the same answer.
    """
    monkeypatch.setattr(
        parameter, "background_aggregate", None if explicit else aggregate, raising=False
    )
    bg = [world[name] for name in hosts]
    extra = {"background_aggregate": aggregate} if explicit else {}
    optimizer = DominatingSetAdapter(
        position_cache=PositionCache([world["fg"], *bg], world["primers"], on_missing="warn"),
        fg_prefixes=[world["fg"]],
        fg_seq_lengths=[LENGTHS["fg"]],
        bg_prefixes=bg,
        bg_seq_lengths=[LENGTHS.get(name, 500_000) for name in hosts],
        config=OptimizerConfig(target_set_size=4),
        conditions=ReactionConditions(temp=30.0, polymerase="phi29"),
        **extra,
    )
    return optimizer, optimizer.compute_metrics(world["primers"])


def _objective(metrics, **limits):
    from neoswga.core.pool_objective import PoolConstraints, PoolObjective

    return PoolObjective(lambda _primers: metrics, PoolConstraints(**limits))


# ----------------------------------------------------------------------
# Nothing moves under the default
# ----------------------------------------------------------------------


@pytest.mark.parametrize("aggregate", [None, SUM, WORST_CASE])
def test_the_pooled_figures_are_what_the_pre_change_tree_printed(world, monkeypatch, aggregate):
    """The pooled fields keep their pooled meaning under EVERY aggregate.

    That is the design decision this phase rests on: rather than overwrite
    `selectivity_density` with the worst host's under `worst-case` -- which
    would make one key mean two things across the summary, the ensemble table
    and the report -- the pooled fields are left alone and the reduction gets
    its own name. So this assertion holds for `worst-case` too, and no output
    has to be re-read to know which figure it is showing.

    The constants come from `63875a1`, not from a second run of this code.
    """
    _optimizer, metrics = _metrics(world, monkeypatch, aggregate)

    for name, expected in POOLED_TWO_HOST.items():
        assert getattr(metrics, name) == expected, name


def test_the_default_run_is_sum_and_says_so(world, monkeypatch):
    """An absent key is `sum`, and the metrics record which reading was used."""
    _optimizer, metrics = _metrics(world, monkeypatch, None)

    assert metrics.background_aggregate == SUM
    assert metrics.to_dict()["background_aggregate"] == SUM


def test_constructing_without_the_argument_does_not_select_worst_case(world, monkeypatch):
    """`BaseOptimizer.__init__`'s default was the OPPOSITE of the key's.

    It read `background_aggregate: str = "worst-case"` while the params.json key
    defaults to `sum`, and nothing passed the argument -- so the first caller to
    pass it would have been the first to notice. Neither default may stand in
    for the other: with nothing passed the value comes from the one seam, which
    answers `sum` for a design that did not set the key.
    """
    monkeypatch.setattr(parameter, "background_aggregate", None, raising=False)
    optimizer, _metrics_ = _metrics(world, monkeypatch, None)

    assert optimizer.background_aggregate == SUM

    import inspect

    from neoswga.core.base_optimizer import BaseOptimizer

    default = inspect.signature(BaseOptimizer.__init__).parameters["background_aggregate"].default
    assert default is None, "a literal default here can disagree with the key"


def test_the_explicit_argument_and_the_globals_agree(world, monkeypatch):
    """`unified_optimizer` passes the resolved value; a step also resolves it.

    Two routes, one definition -- the property `background_aggregate.py`'s
    docstring claims. If they could differ, a run would select under one
    reading and report under the other.
    """
    _o1, by_globals = _metrics(world, monkeypatch, WORST_CASE)
    _o2, by_argument = _metrics(world, monkeypatch, WORST_CASE, explicit=True)

    assert by_globals.to_dict() == by_argument.to_dict()


def test_unified_optimizer_passes_the_resolved_value():
    """The flag reaches the code that decides, rather than being resolved twice.

    Asserted on the source of the one line that does it, because the
    alternative -- letting every optimizer reach into `parameter` -- is Known
    Issue 8 one layer down and would pass any behavioural test.
    """
    import inspect

    from neoswga.core import unified_optimizer

    source = inspect.getsource(unified_optimizer.run_optimization)

    assert 'kwargs["background_aggregate"] = resolve_background_aggregate()' in source


# ----------------------------------------------------------------------
# One host: the two readings are the same design
# ----------------------------------------------------------------------


@pytest.mark.parametrize("host", ["small", "large"])
def test_worst_case_with_one_host_is_sum_with_one_host_exactly(world, monkeypatch, host):
    """Every figure, not the ones that happen to be compared elsewhere.

    `to_dict()` is the whole serialized record, so a field added later is
    covered without this test being touched. Exact equality rather than
    `approx`: with one host the pooled figures ARE the per-host figures, and
    `per_host_metrics` carries them over rather than recomputing a coverage
    fraction that `(f * L) / L` would perturb in the last bit.

    `background_aggregate` itself is the one key that differs, and it has to:
    the record says which reading the run ASKED for, which is a property of the
    request rather than of the panel. Everything derived from it is equal.
    """
    _o1, pooled = _metrics(world, monkeypatch, SUM, hosts=(host,))
    _o2, worst = _metrics(world, monkeypatch, WORST_CASE, hosts=(host,))

    requested = {"background_aggregate"}
    left = {k: v for k, v in pooled.to_dict().items() if k not in requested}
    right = {k: v for k, v in worst.to_dict().items() if k not in requested}

    assert left == right
    assert pooled.background_aggregate == SUM
    assert worst.background_aggregate == WORST_CASE


def test_one_host_still_matches_the_pre_change_tree(world, monkeypatch):
    """And the one-host figures are the ones `63875a1` printed."""
    _optimizer, metrics = _metrics(world, monkeypatch, WORST_CASE, hosts=("small",))

    for name, expected in POOLED_ONE_HOST.items():
        assert getattr(metrics, name) == expected, name


@pytest.mark.parametrize("aggregate", [SUM, WORST_CASE])
def test_a_limit_on_one_host_binds_identically_under_both(world, monkeypatch, aggregate):
    """The limits agree too, which is what `limit_figure` has to guarantee."""
    _optimizer, metrics = _metrics(world, monkeypatch, aggregate, hosts=("small",))

    assert _objective(metrics, min_selectivity_density=0.5).violations(["X"]) == ()
    assert _objective(metrics, min_selectivity_density=2.0).violations(["X"]) == (
        "selectivity below minimum",
    )


# ----------------------------------------------------------------------
# Two hosts: the limit binds per host
# ----------------------------------------------------------------------


def test_the_worst_host_is_measured_on_its_own_length(world, monkeypatch):
    """The per-host table, and the three reductions over it."""
    _optimizer, metrics = _metrics(world, monkeypatch, WORST_CASE)

    assert set(metrics.per_host) == {world["small"], world["large"]}
    small = metrics.per_host[world["small"]]
    large = metrics.per_host[world["large"]]

    assert small["sites"] == 8 and small["length_bp"] == 40_000
    assert small["sites_per_mb"] == pytest.approx(200.0)
    assert small["selectivity_density"] == pytest.approx(1.0)
    # The large host is a measurement, not an absence: it was scanned and binds
    # nowhere, so its zero is a result and its density is the "nothing
    # detected" sentinel rather than a small number.
    assert large["sites"] == 0 and large["unavailable"] is None
    assert large["selectivity_density"] == MAX_SELECTIVITY

    assert metrics.worst_host == world["small"]
    assert metrics.worst_host_selectivity_density == pytest.approx(1.0)
    assert metrics.worst_host_background_sites == 8
    assert metrics.worst_host_bg_coverage == pytest.approx(0.950275)


def test_a_density_floor_passes_pooled_and_fails_per_host(world, monkeypatch):
    """The verification this phase exists for.

    A panel clean on the 1.72 Mb host and dense on the 40 kb one. Pooled, the
    small host is 2.3 percent of the denominator and the design reports a
    density of 44.0; on the small host's own length it is 1.0. A floor of 10
    therefore passes under `sum` and fails under `worst-case`, with the SAME
    limit key and the same panel.
    """
    _o1, pooled = _metrics(world, monkeypatch, SUM)
    _o2, worst = _metrics(world, monkeypatch, WORST_CASE)

    assert _objective(pooled, min_selectivity_density=10.0).violations(["X"]) == ()
    assert _objective(worst, min_selectivity_density=10.0).violations(["X"]) == (
        "selectivity below minimum",
    )
    # And by how much, because a search ranks on the shortfall rather than on
    # the count of violations.
    assert _objective(pooled, min_selectivity_density=10.0).shortfall(["X"]) == 0.0
    assert _objective(worst, min_selectivity_density=10.0).shortfall(["X"]) > 0.0


def test_a_host_coverage_ceiling_binds_per_host(world, monkeypatch):
    """`max_host_coverage`, the one limit that sees host site POSITION.

    Pooled, the panel covers 2.2 percent of the background; on the small host
    it covers 95 percent. A ceiling of 0.5 passes pooled and fails per host.
    """
    _o1, pooled = _metrics(world, monkeypatch, SUM)
    _o2, worst = _metrics(world, monkeypatch, WORST_CASE)

    assert _objective(pooled, max_host_coverage=0.5).violations(["X"]) == ()
    assert _objective(worst, max_host_coverage=0.5).violations(["X"]) == (
        "host coverage above maximum",
    )


def test_a_background_site_ceiling_reduces_to_the_densest_host(world, monkeypatch):
    """`max_background_sites` takes the HIGHEST per-host count, not the lowest.

    Each limit reduces independently, because `worst-case` means the panel must
    clear every limit whichever host turns up -- and the host with the lowest
    density need not be the one with the most sites. Here one host carries all
    eight and the other none, so the reduction is 8 and the pooled sum is 8 as
    well: the assertion that matters is that the reduction is the MAXIMUM over
    hosts rather than a sum that could exceed any single host's ceiling.
    """
    _optimizer, worst = _metrics(world, monkeypatch, WORST_CASE)

    assert worst.worst_host_background_sites == max(row["sites"] for row in worst.per_host.values())
    assert _objective(worst, max_background_sites=4).violations(["X"]) == (
        "background sites above maximum",
    )
    assert _objective(worst, max_background_sites=12).violations(["X"]) == ()


def test_no_limit_key_was_added(world, monkeypatch):
    """One limit, applied per host. The Phase 7 guard, asserted.

    A per-host threshold list is what this must not become: `_LIMITS` and
    `LIMIT_KEYS` carry exactly the rows they carried before, and the per-host
    behaviour comes from which FIGURE each row is read against.
    """
    from neoswga.core.panel_acceptance import LIMIT_KEYS
    from neoswga.core.pool_objective import _LIMITS

    assert not [key for key in LIMIT_KEYS if "host" in key and key != "max_host_coverage"]
    assert not [field for field, _m, _s, _msg, _bg in _LIMITS if field.startswith("per_host")]


def test_the_reduced_metrics_are_exactly_the_background_measured_limits():
    """`AGGREGATE_REDUCED` must track `_LIMITS`, or a new limit reads the pool.

    A fourth background-measured limit added without a reduction here would be
    evaluated on the pooled figure under `worst-case` while every other limit
    was evaluated per host -- the half-wired dispatch this file exists to catch.
    """
    from neoswga.core.pool_objective import _LIMITS

    background_measured = {metric for _f, metric, _s, _msg, needs_bg in _LIMITS if needs_bg}

    assert set(per_host_metrics.AGGREGATE_REDUCED) == background_measured


def test_limit_figure_under_sum_is_the_plain_lookup(world, monkeypatch):
    """What keeps the default path byte-identical: no dispatch at all."""
    _optimizer, metrics = _metrics(world, monkeypatch, SUM)

    for metric in per_host_metrics.AGGREGATE_REDUCED:
        assert per_host_metrics.limit_figure(metrics, metric) == getattr(metrics, metric)
    assert per_host_metrics.limit_figure(metrics, "max_gap") == metrics.max_gap


def test_a_metrics_object_with_no_reading_is_treated_as_sum():
    """A hand-built or legacy metrics object must not acquire a per-host reading.

    `evaluate-set` and several tests construct objects that carry the pooled
    fields and nothing else. Reading them as `worst-case` would silently refuse
    every limit on them.
    """

    class Bare:
        selectivity_density = 7.0

    assert per_host_metrics.limit_figure(Bare(), "selectivity_density") == 7.0


# ----------------------------------------------------------------------
# Unknown is not the worst of the rest
# ----------------------------------------------------------------------


def test_an_unscanned_host_makes_the_worst_case_unknown(world, monkeypatch):
    """A reduction over a set with an unmeasured member is unmeasured.

    `absent` is indexed and holds no entry for any panel oligo, so
    `get_positions` answers with an empty array -- a zero that reads as perfect
    specificity. `has_entry` is the question asked instead, and the three
    reductions come back `None` rather than the worst of the two hosts that did
    answer.
    """
    _optimizer, metrics = _metrics(world, monkeypatch, WORST_CASE, hosts=("small", "absent"))

    assert metrics.worst_host_selectivity_density is None
    assert metrics.worst_host_background_sites is None
    assert metrics.worst_host_bg_coverage is None
    assert metrics.worst_host is None
    row = metrics.per_host[world["absent"]]
    assert row["sites"] is None
    assert "never scanned" in row["unavailable"]


@pytest.mark.parametrize(
    ("limits", "message"),
    [
        ({"min_selectivity_density": 0.001}, "selectivity below minimum"),
        ({"max_background_sites": 10_000_000}, "background sites above maximum"),
        ({"max_host_coverage": 0.999}, "host coverage above maximum"),
    ],
)
def test_a_background_limit_on_an_unscanned_host_is_refused(world, monkeypatch, limits, message):
    """Refused, not satisfied -- and the limits chosen here are ones every
    measured panel clears easily, so a pass could only come from the missing
    measurement being read as a favourable number."""
    _optimizer, worst = _metrics(world, monkeypatch, WORST_CASE, hosts=("small", "absent"))
    _o2, pooled = _metrics(world, monkeypatch, SUM, hosts=("small", "absent"))

    assert _objective(pooled, **limits).violations(["X"]) == ()
    assert _objective(worst, **limits).violations(["X"]) == (message,)
    assert _objective(worst, **limits).shortfall(["X"]) == float("inf")


def test_an_unscanned_host_leaves_the_pooled_figures_alone(world, monkeypatch):
    """Under `sum` the unscanned host is exactly as invisible as it was.

    Phase 7c must not turn `sum` into a stricter gate by the back door. That an
    absent background count passes under `sum` is the documented difference
    from `worst-case`, recorded in the schema for the step-2 gate and the same
    here.
    """
    _o1, with_absent = _metrics(world, monkeypatch, SUM, hosts=("small", "absent"))

    assert with_absent.total_bg_sites == 8
    assert with_absent.selectivity_density == pytest.approx(8 / 40000 / (8 / 540000))


# ----------------------------------------------------------------------
# Known Issue 6: no per-host output invites a ratio comparison
# ----------------------------------------------------------------------


def test_no_per_host_output_carries_a_selectivity_ratio(world, monkeypatch):
    """The count ratio moves with how big a host happens to be.

    So a per-host table of it would invite exactly the comparison Known Issue 6
    records as invalid, and the figure that survives the substitution --
    density -- is the only ratio offered per host. Each row says so itself.
    """
    _optimizer, metrics = _metrics(world, monkeypatch, WORST_CASE)
    serialized = json.dumps(metrics.to_dict()["per_host"])

    assert "selectivity_ratio" not in serialized
    assert "_ratio" not in serialized
    for row in metrics.per_host.values():
        assert row["comparable"] == "selectivity_density"

    # And no reduction over hosts of the pooled ratio exists either.
    assert not [name for name in metrics.to_dict() if name.startswith("worst_host_selectivity_r")]


def test_the_pooled_ratio_is_not_reduced_over_hosts(world, monkeypatch):
    """`selectivity_ratio` is the pooled one under both readings.

    It is also the term `normalized_score` reads, so reducing it over hosts
    would change which method an ensemble run picks -- a decision no
    measurement supports.
    """
    _o1, pooled = _metrics(world, monkeypatch, SUM)
    _o2, worst = _metrics(world, monkeypatch, WORST_CASE)

    assert pooled.selectivity_ratio == worst.selectivity_ratio
    assert pooled.normalized_score() == worst.normalized_score()


# ----------------------------------------------------------------------
# The search sees what the acceptance path will hold it to
# ----------------------------------------------------------------------


def test_the_focused_evaluator_records_the_same_reading(world, monkeypatch):
    """`compute_pool_metrics` is the search evaluator; it must not disagree.

    Without this the bounded repair would chase the pooled density and then be
    rejected on the worst host's -- two definitions of one quantity, which is
    the defect `pool_metrics` was written with a pinning test to prevent.
    """
    optimizer, full = _metrics(world, monkeypatch, WORST_CASE)
    focused = optimizer.compute_pool_metrics(world["primers"])

    assert focused.background_aggregate == full.background_aggregate
    assert focused.worst_host == full.worst_host
    assert focused.worst_host_selectivity_density == full.worst_host_selectivity_density
    assert focused.worst_host_background_sites == full.worst_host_background_sites
    assert per_host_metrics.limit_figure(
        focused, "selectivity_density"
    ) == per_host_metrics.limit_figure(full, "selectivity_density")


def test_the_focused_evaluator_refuses_a_host_coverage_limit_as_it_always_did(world, monkeypatch):
    """It carries no background coverage, so that limit routes elsewhere.

    Unchanged behaviour, asserted because adding `worst_host_bg_coverage` here
    would have made `plan-pool` start enforcing a limit it cannot measure.
    """
    optimizer, _full = _metrics(world, monkeypatch, WORST_CASE)
    focused = optimizer.compute_pool_metrics(world["primers"])

    assert per_host_metrics.limit_figure(focused, "bg_coverage") is None
    assert not hasattr(focused, "worst_host_bg_coverage")


# ----------------------------------------------------------------------
# The search-time host charge
# ----------------------------------------------------------------------


def test_the_host_charge_on_one_host_is_the_pooled_sum_exactly():
    """Integer-identical, so a single-host run cannot move under `worst-case`."""
    for sites in (0, 1, 7, 1234):
        assert per_host_metrics.host_site_charge({"h": sites}, [3_000_000], WORST_CASE) == sites


def test_the_host_charge_has_one_reading_and_refuses_the_other():
    """The search-time host term is the pooled sum, and says so by raising.

    Phase 7c first shipped a `worst-case` branch here that charged the worst
    host's site density carried back onto the pooled length. The ranking was
    right -- counts are not comparable across hosts of different size -- and the
    magnitude was not: it multiplies the charge by
    `sum(host lengths) / worst host length`, 44x in this world and 1029x on this
    project's own three-host panel, while `background_weight` scales this term
    against a coverage term and was chosen against the pooled magnitude.

    So the reduction is withdrawn rather than left behind a note, and the
    function raises rather than pooling silently: a route that skipped the
    up-front refusal must not get a number instead of an error.
    `tests/test_the_host_term_has_one_calibration.py` holds the refusal end to
    end.
    """
    from neoswga.core.exceptions import InvalidDesignRequest

    sites = {"small": 9, "large": 100}
    lengths = [40_000, 1_720_000]

    assert per_host_metrics.host_site_charge(sites, lengths, SUM) == 109

    with pytest.raises(InvalidDesignRequest) as caught:
        per_host_metrics.host_site_charge(sites, lengths, WORST_CASE)
    message = str(caught.value)
    assert "background_weight" in message
    assert "What would settle it" in message


# ----------------------------------------------------------------------
# What the record and the report say
# ----------------------------------------------------------------------


def test_the_panel_assessment_carries_the_hosts_and_the_aggregate(world, monkeypatch):
    """`PanelAssessment` gains `per_host` beside `per_target`.

    Carried from the metrics rather than remeasured: the evaluator measured each
    host during the run, and a second measurement here is how two records of
    one panel come to disagree.
    """
    from neoswga.core.panel_evaluation import evaluate_panel

    _optimizer, metrics = _metrics(world, monkeypatch, WORST_CASE)

    class Request:
        request_hash = "abc"
        target_size = len(world["primers"])
        max_dimer_bp = 20
        coverage_reach = 3000

    assessment = evaluate_panel(Request(), world["primers"], metrics)

    assert set(assessment.per_host) == {world["small"], world["large"]}
    assert assessment.background_aggregate == WORST_CASE
    saved = assessment.as_dict()
    assert saved["background_aggregate"] == WORST_CASE
    assert saved["per_host"][world["small"]]["selectivity_density"] == pytest.approx(1.0)
    assert saved["metrics"]["worst_host_selectivity_density"]["value"] == pytest.approx(1.0)


def test_the_report_loader_keeps_both_keys(tmp_path):
    """A key the loader drops is written and discarded, which 7b found happening.

    `occupancy_ranking_mode` was added, written to `filter_stats.json` and
    silently dropped by a hand-maintained whitelist. So this drives a saved
    assessment through the real loader rather than trusting that the field
    arrives.
    """
    from neoswga.core.report.metrics import (
        _ASSESSMENT_CACHE,
        SpecificityMetrics,
        _apply_panel_assessment,
    )

    saved = {
        "metrics": {"fg_coverage": {"value": 0.5, "unavailable": None}},
        "per_target": {},
        "per_host": {"small": {"sites": 8, "length_bp": 40_000}},
        "background_aggregate": WORST_CASE,
    }
    report = tmp_path / "step4_improved_df_validation.json"
    report.write_text(json.dumps({"ok": True, "issues": [], "assessment": saved}))
    _ASSESSMENT_CACHE.pop(report, None)

    class Coverage:
        total_bases = 0
        overall_coverage = 0.0
        from_optimizer = False
        covered_bases = 0

    class Metrics:
        coverage = Coverage()
        specificity = SpecificityMetrics()
        thermodynamics = None
        uniformity = None
        panel_assessment = None

    metrics = Metrics()
    _apply_panel_assessment(metrics, tmp_path)

    assert metrics.specificity.background_aggregate == WORST_CASE
    assert metrics.specificity.per_host == {"small": {"sites": 8, "length_bp": 40_000}}


def test_the_report_renders_a_per_host_table_only_with_two_or_more_hosts():
    """Following `_render_per_target_coverage`, including its threshold."""
    from neoswga.core.report.metrics import SpecificityMetrics
    from neoswga.core.report.per_host_table import render_per_host_background

    one = SpecificityMetrics(per_host={"small": {"sites": 8, "length_bp": 40_000}})
    assert render_per_host_background(one) == ""
    assert render_per_host_background(None) == ""

    two = SpecificityMetrics(
        background_aggregate=WORST_CASE,
        per_host={
            "small": {
                "length_bp": 40_000,
                "sites": 8,
                "sites_per_mb": 200.0,
                "coverage": 0.95,
                "selectivity_density": 1.0,
                "unavailable": None,
            },
            "absent": {"length_bp": 500_000, "unavailable": "never scanned"},
        },
    )
    html = render_per_host_background(two)

    assert "Per-Host Background" in html
    assert "worst-case" in html
    assert "selectivity density" in html
    # An unmeasured host says so rather than rendering as the best possible
    # result: zero sites and zero coverage are the FAVOURABLE values here.
    assert "not measured: never scanned" in html
    assert ">0</td>" not in html
    # No ratio column, for the Known Issue 6 reason.
    assert "Selectivity ratio" not in html


def test_the_summary_json_names_the_aggregate_beside_every_figure(world, monkeypatch):
    """`to_dict()` is what `step4_improved_df_summary.json` carries.

    Both namings in one record: the pooled figures under the keys they always
    had, and the reductions under names that say what they are, with the
    aggregate beside them. No key means two different things in two outputs.
    """
    _optimizer, metrics = _metrics(world, monkeypatch, WORST_CASE)
    saved = metrics.to_dict()

    assert saved["background_aggregate"] == WORST_CASE
    assert saved["selectivity_density"] == POOLED_TWO_HOST["selectivity_density"]
    assert saved["worst_host_selectivity_density"] == pytest.approx(1.0)
    assert saved["worst_host_background_sites"] == 8
    assert saved["worst_host_bg_coverage"] == pytest.approx(0.950275)
    assert saved["worst_host"] == world["small"]
    assert set(saved["per_host"]) == {world["small"], world["large"]}
    json.dumps(saved)


# ----------------------------------------------------------------------
# The modelled-profile path
# ----------------------------------------------------------------------


def test_the_modelled_profile_panel_reads_the_same_key(world, monkeypatch):
    """One key, one answer to what "several hosts" means.

    `_modelled_site_load` defaulted to `worst-case` through
    `getattr(self, "background_aggregate", "worst-case")` while the key it
    shadowed defaults to `sum`. That second default is gone: a panel of modelled
    profiles now aggregates by the key, which also makes `worst-case` over such
    a panel reachable for the first time. Only a panel of SEVERAL profiles is
    affected, and nothing in the CLI or the schema supplies one.
    """
    from neoswga.core.host_profile import HostProfile

    class Profile(HostProfile):
        pass

    monkeypatch.setattr(parameter, "background_aggregate", WORST_CASE, raising=False)
    optimizer = DominatingSetAdapter(
        position_cache=PositionCache([world["fg"]], world["primers"]),
        fg_prefixes=[world["fg"]],
        fg_seq_lengths=[LENGTHS["fg"]],
        bg_prefixes=[],
        bg_seq_lengths=[],
        config=OptimizerConfig(target_set_size=4),
    )

    assert optimizer.background_aggregate == WORST_CASE

    import inspect

    from neoswga.core.base_optimizer import BaseOptimizer

    source = inspect.getsource(BaseOptimizer._modelled_site_load)
    assert "mode = self.background_aggregate" in source
    assert '"worst-case"' not in source, "the second default is what this removes"
