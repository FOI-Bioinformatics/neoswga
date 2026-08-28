"""Tests for the Gillespie trajectory itself (`GillespieSimulator.simulate`).

The existing stochastic tests construct the simulator and check its parameters;
none of them ran the reaction loop. That loop is where the algorithm actually
lives -- propensity calculation, exponential waiting times, reaction selection,
resource depletion and the termination conditions -- and it is what produces the
enrichment figures used to validate the network heuristics.

Networks here are built on a real `networkx` graph rather than a Mock, so
`_build_extension_network` and the edge-distance logic run for real. Mocking the
graph would leave the extension reactions unreachable, which is most of the
model.
"""

import pytest

nx = pytest.importorskip("networkx")

from neoswga.core.stochastic_simulator import (
    GillespieSimulator,
    ReactionParameters,
    ReactionState,
)


class _Network:
    """Minimal stand-in exposing what the simulator reads: sites and a graph."""

    def __init__(self, sites, spacing=2_000):
        self.binding_sites = list(sites)
        self.graph = nx.DiGraph()
        self.graph.add_nodes_from(self.binding_sites)
        # Chain each site to the next, with a reachable distance.
        for a, b in zip(self.binding_sites, self.binding_sites[1:]):
            self.graph.add_edge(a, b, distance=abs(b - a))


@pytest.fixture
def simulator():
    fg = _Network([1_000, 3_000, 5_000, 7_000])
    bg = _Network([2_000, 6_000])
    return GillespieSimulator(
        primers=["TTGACCATGA", "GCATTACGGT"],
        fg_network=fg,
        bg_network=bg,
        params=ReactionParameters(),
        seed=1234,
    )


# Simulated durations here are fractions of a second, not the 8-hour default.
# That is not a shortcut: exact Gillespie advances one molecular event at a
# time, and at these propensities even a 3-site network needs on the order of
# 10^8 events to reach 8 hours. The loop is what is under test, and it runs
# identically at any duration.
SIM_TIME = 0.05
SAMPLE = 0.01


def _run(sim, max_time=SIM_TIME, sample_interval=SAMPLE):
    return sim.simulate(max_time=max_time, sample_interval=sample_interval)


# ----------------------------------------------------------------------
# The extension network the loop depends on
# ----------------------------------------------------------------------


def test_extension_network_is_built_from_graph_edges(simulator):
    """Without edges there are no extension reactions, and the model reduces to
    binding and unbinding -- which would still 'run' and report nothing."""
    assert simulator.fg_extensions, "no foreground extension edges were built"
    for _site, targets in simulator.fg_extensions.items():
        for _target, distance in targets:
            assert distance > 0


def test_sites_are_indexed_not_used_as_positions(simulator):
    """Extensions are keyed by index into `fg_sites`, not by genome position.

    Confusing the two silently addresses the wrong site.
    """
    assert set(simulator.fg_extensions) <= set(range(len(simulator.fg_sites)))


# ----------------------------------------------------------------------
# The trajectory
# ----------------------------------------------------------------------


def test_simulate_returns_a_time_ordered_history(simulator):
    history = _run(simulator)

    assert len(history) >= 2, "no snapshots beyond the initial state"
    times = [snap["time"] for snap in history]
    assert times == sorted(times), "snapshots are not in time order"
    assert times[0] == 0.0


def test_simulation_advances_time_towards_the_limit(simulator):
    history = _run(simulator)
    assert history[-1]["time"] > 0


def test_longer_runs_produce_more_snapshots(simulator):
    short = _run(simulator, max_time=SIM_TIME / 4)
    long = _run(simulator, max_time=SIM_TIME * 4)
    assert len(long) > len(short)


def test_sample_interval_controls_snapshot_density(simulator):
    coarse = _run(simulator, max_time=SIM_TIME, sample_interval=SIM_TIME / 2)
    fine = _run(simulator, max_time=SIM_TIME, sample_interval=SIM_TIME / 20)
    assert len(fine) > len(coarse)


def test_snapshots_carry_the_fields_the_validation_reads(simulator):
    snap = _run(simulator)[-1]
    for key in (
        "time",
        "fg_molecules",
        "fg_products",
        "total_fg",
        "bg_molecules",
        "bg_products",
        "total_bg",
        "enrichment",
        "fg_amplification",
    ):
        assert key in snap, f"{key} missing from snapshot"


def test_amplification_products_accumulate(simulator):
    """A reaction that never extends is not amplifying anything."""
    history = _run(simulator, max_time=SIM_TIME * 4)
    assert history[-1]["fg_products"] > 0, "no target products were ever synthesised"


def test_molecule_counts_never_go_negative(simulator):
    """Resource depletion and unbinding both decrement counters."""
    for snap in _run(simulator, max_time=SIM_TIME * 4):
        assert snap["fg_products"] >= 0
        assert snap["bg_products"] >= 0
        assert snap["total_fg"] >= snap["fg_molecules"]


def test_enrichment_is_finite_and_positive(simulator):
    """`enrichment` divides by background; a zero denominator must not produce
    inf or nan, since these values are reported and serialised."""
    import math

    for snap in _run(simulator):
        assert math.isfinite(snap["enrichment"])
        assert snap["enrichment"] >= 0


# ----------------------------------------------------------------------
# Reproducibility
# ----------------------------------------------------------------------


def test_same_seed_reproduces_the_trajectory():
    """A stochastic simulator that cannot be replayed cannot be debugged.

    The simulator keeps its own RNG rather than using the global numpy state,
    so an unrelated call elsewhere must not perturb the trajectory.
    """

    def run(seed):
        sim = GillespieSimulator(
            primers=["TTGACCATGA"],
            fg_network=_Network([1_000, 3_000, 5_000]),
            bg_network=_Network([2_000]),
            params=ReactionParameters(),
            seed=seed,
        )
        return [
            (s["time"], s["fg_products"])
            for s in sim.simulate(max_time=SIM_TIME, sample_interval=SAMPLE)
        ]

    assert run(7) == run(7)


def test_different_seeds_diverge():
    def run(seed):
        sim = GillespieSimulator(
            primers=["TTGACCATGA"],
            fg_network=_Network([1_000, 3_000, 5_000]),
            bg_network=_Network([2_000]),
            params=ReactionParameters(),
            seed=seed,
        )
        return [s["time"] for s in sim.simulate(max_time=SIM_TIME, sample_interval=SAMPLE)]

    assert run(1) != run(2), "the seed has no effect; the run is not stochastic"


def test_global_numpy_rng_is_not_disturbed():
    """Using np.random directly would make every other seeded test order-dependent."""
    import numpy as np

    np.random.seed(99)
    before = np.random.random()

    np.random.seed(99)
    sim = GillespieSimulator(
        primers=["TTGACCATGA"],
        fg_network=_Network([1_000, 3_000]),
        bg_network=_Network([2_000]),
        params=ReactionParameters(),
        seed=5,
    )
    sim.simulate(max_time=SIM_TIME, sample_interval=SAMPLE)
    after = np.random.random()

    assert before == after, "the simulation consumed the global numpy RNG"


# ----------------------------------------------------------------------
# Termination
# ----------------------------------------------------------------------


def test_a_network_with_no_sites_terminates_rather_than_hanging(caplog):
    """No binding sites means no reactions; it must stop and say so."""
    import logging

    sim = GillespieSimulator(
        primers=["TTGACCATGA"],
        fg_network=_Network([]),
        bg_network=_Network([]),
        params=ReactionParameters(),
        seed=1,
    )
    with caplog.at_level(logging.WARNING):
        history = sim.simulate(max_time=SIM_TIME, sample_interval=SAMPLE)

    assert history, "no snapshot at all, not even the initial state"
    assert "No reactions possible" in caplog.text or "Zero total propensity" in caplog.text


def test_zero_max_time_returns_only_the_initial_state(simulator):
    history = simulator.simulate(max_time=0.0)
    assert len(history) == 1
    assert history[0]["time"] == 0.0


# ----------------------------------------------------------------------
# The cost of the documented default
# ----------------------------------------------------------------------


def test_iteration_cap_terminates_and_explains(simulator, caplog):
    """Without a bound, an ambitious max_time is indistinguishable from a hang.

    Exact Gillespie advances one event at a time, so the event count scales with
    total propensity and therefore with the number of binding sites. Measured on
    this implementation, reaching the 8-hour default takes roughly 2 wall hours
    at 3 sites, 5 days at 20 and several months at 100 -- and a real design has
    hundreds. The cap does not make such a run accurate; it makes it stop and
    report how far it got.
    """
    import logging

    with caplog.at_level(logging.WARNING):
        history = simulator.simulate(max_time=28_800.0, sample_interval=60.0, max_iterations=2_000)

    assert history, "no trajectory returned"
    assert "Stopped after" in caplog.text
    assert history[-1]["time"] < 28_800.0, "the full duration was somehow reached"


def test_the_cap_can_be_disabled(simulator):
    """A caller who intends to wait must be able to."""
    history = simulator.simulate(max_time=SIM_TIME, sample_interval=SAMPLE, max_iterations=None)
    assert history[-1]["time"] > 0


def test_a_run_under_the_cap_is_unaffected():
    """The bound must not perturb runs that complete within it.

    Fresh simulators per run: the RNG advances with each call, so reusing one
    instance would diverge for reasons unrelated to the cap.
    """

    def run(**kwargs):
        sim = GillespieSimulator(
            primers=["TTGACCATGA"],
            fg_network=_Network([1_000, 3_000, 5_000]),
            bg_network=_Network([2_000]),
            params=ReactionParameters(),
            seed=11,
        )
        return [s["time"] for s in sim.simulate(SIM_TIME, SAMPLE, **kwargs)]

    assert run() == run(max_iterations=None)
