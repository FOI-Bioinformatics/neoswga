"""Coverage weighted by occupancy, so both halves of the objective agree.

Phase 3 made selectivity occupancy-weighted and left coverage a count of
binding sites. That put the two halves of the objective on different footings:
a set could be rewarded for reaching sites it does not actually occupy at the
reaction temperature, and a condition search could spend coverage it could not
see it was spending.

The bias is large and it is not uniform, so it does not cancel in a comparison.
Measured on Prevotella at equiphi29 42 C:

    set     raw coverage   effective coverage   mean occupancy
    k=10        0.520            0.245               0.36
    k=11        0.324            0.220               0.56
    k=12        0.403            0.387               0.86

The raw numbers rank k=10 first; the effective ones rank k=12 first. A design
chosen on the raw figure is chosen on sites that are mostly not bound.

This is also what connects additives to coverage. DMSO and betaine lower
effective Tm, moving every site down its occupancy curve, and against a count
that cost is invisible. With it visible, at equiphi29 42 C:

    k=12  + 5% DMSO + 1M betaine   effective coverage -9%,  selectivity +51%
    k=10  + 5% DMSO + 1M betaine   effective coverage -64%, selectivity +13%

which is the whole design rule in two lines: stringency is affordable when the
primers melt above the reaction temperature and ruinous when they do not.
"""

import random

import pytest

h5py = pytest.importorskip("h5py")

import numpy as np

from neoswga.core.base_optimizer import OptimizerConfig
from neoswga.core.dominating_set_adapter import DominatingSetAdapter
from neoswga.core.position_cache import PositionCache
from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.thermodynamics import reverse_complement

GENOME = 40_000
K = 10


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    """A genome with two primers of deliberately different stability.

    One is GC-rich and stays bound at the reaction temperature; one is AT-rich
    and mostly does not. They bind the same number of sites, so a count-based
    coverage cannot tell them apart and an occupancy-weighted one must.
    """
    rng = random.Random(4)
    directory = tmp_path_factory.mktemp("effcov")
    seq = [rng.choice("ACGT") for _ in range(GENOME)]

    stable = "GCGCGGCCGC"
    weak = "ATATTAAATA"
    for i, primer in enumerate((stable, weak)):
        for j in range(4):
            at = 2_000 + i * 1_000 + j * 9_000
            seq[at : at + K] = list(primer)

    prefix = str(directory / "fg")
    joined = "".join(seq)
    with h5py.File(f"{prefix}_{K}mer_positions.h5", "w") as handle:
        for key in sorted({x for p in (stable, weak) for x in (p, reverse_complement(p))}):
            hits, at = [], joined.find(key)
            while at != -1:
                hits.append(at)
                at = joined.find(key, at + 1)
            if hits:
                handle.create_dataset(key, data=np.array(hits, dtype=np.int32))

    return {"prefix": prefix, "stable": stable, "weak": weak}


def build(world, primers, temp=42.0, **additives):
    cache = PositionCache([world["prefix"]], primers)
    return DominatingSetAdapter(
        position_cache=cache,
        fg_prefixes=[world["prefix"]],
        fg_seq_lengths=[GENOME],
        bg_prefixes=[],
        bg_seq_lengths=[],
        config=OptimizerConfig(target_set_size=len(primers), extension_reach=1000),
        conditions=ReactionConditions(temp=temp, polymerase="equiphi29", **additives),
    ).compute_metrics(primers)


# ----------------------------------------------------------------------
# It measures occupancy, not just reach
# ----------------------------------------------------------------------


def test_effective_coverage_never_exceeds_raw_coverage(world):
    """Occupancy is a probability, so weighting by it can only remove coverage.
    A value above the raw figure would mean the weighting had a sign error."""
    for primers in ([world["stable"]], [world["weak"]], [world["stable"], world["weak"]]):
        m = build(world, primers)
        assert m.effective_fg_coverage <= m.fg_coverage + 1e-9


def test_a_weakly_binding_primer_loses_more_coverage_than_a_stable_one(world):
    """The point of the metric. Both primers bind four sites and reach the same
    distance, so raw coverage cannot separate them."""
    stable = build(world, [world["stable"]])
    weak = build(world, [world["weak"]])

    assert stable.fg_coverage == pytest.approx(weak.fg_coverage, rel=0.2)
    assert weak.effective_fg_coverage < stable.effective_fg_coverage / 2


def test_raising_the_temperature_costs_effective_coverage_only(world):
    """Raw coverage is a property of where primers bind and cannot respond to
    conditions at all -- which is exactly the blindness being fixed."""
    cool = build(world, [world["stable"], world["weak"]], temp=38.0)
    hot = build(world, [world["stable"], world["weak"]], temp=48.0)

    assert cool.fg_coverage == hot.fg_coverage
    assert hot.effective_fg_coverage < cool.effective_fg_coverage


def test_an_additive_costs_effective_coverage(world):
    """DMSO lowers effective Tm, so it moves sites down their occupancy curve.
    Without this term a condition search sees the specificity gain and none of
    the coverage it is paying with."""
    plain = build(world, [world["stable"], world["weak"]])
    dmso = build(world, [world["stable"], world["weak"]], dmso_percent=5.0)

    assert dmso.effective_fg_coverage < plain.effective_fg_coverage


# ----------------------------------------------------------------------
# The arithmetic is the union, not a sum
# ----------------------------------------------------------------------


def test_one_primer_binding_densely_does_not_stack_with_itself(world):
    """Overlapping windows from the SAME primer are one opportunity to prime,
    not several. Multiplying (1 - theta) in per site would understate coverage
    wherever a primer binds densely, which is precisely the good primers."""
    from neoswga.core.occupancy import site_occupancy
    from neoswga.core.thermodynamics import calculate_enthalpy_entropy

    m = build(world, [world["stable"]])

    conditions = ReactionConditions(temp=42.0, polymerase="equiphi29")
    dh, _ = calculate_enthalpy_entropy(world["stable"])
    theta = site_occupancy(dh, conditions.calculate_effective_tm(world["stable"]), 42.0)

    # The primer's sites are spaced 9000 apart with reach 1000, so no two of its
    # windows overlap and effective coverage must be exactly theta * raw. Any
    # self-stacking would multiply (1 - theta) in twice somewhere and drive the
    # ratio below theta.
    assert m.effective_fg_coverage / m.fg_coverage == pytest.approx(theta, rel=1e-4)


def test_two_primers_covering_one_region_combine_as_independent_chances(world):
    """Two chances to prime a region beat one, so the pair must exceed either
    alone -- and stay at or below the raw union."""
    both = build(world, [world["stable"], world["weak"]])
    stable = build(world, [world["stable"]])

    assert both.effective_fg_coverage > stable.effective_fg_coverage
    assert both.effective_fg_coverage <= both.fg_coverage + 1e-9


# ----------------------------------------------------------------------
# Absence of conditions is reported, not invented
# ----------------------------------------------------------------------


def test_without_conditions_it_reports_zero_rather_than_guessing(world):
    """There is no temperature at which to evaluate occupancy, and a fabricated
    number here would be indistinguishable from a measured one."""
    cache = PositionCache([world["prefix"]], [world["stable"]])
    metrics = DominatingSetAdapter(
        position_cache=cache,
        fg_prefixes=[world["prefix"]],
        fg_seq_lengths=[GENOME],
        bg_prefixes=[],
        bg_seq_lengths=[],
        config=OptimizerConfig(target_set_size=1, extension_reach=1000),
        conditions=None,
    ).compute_metrics([world["stable"]])

    assert metrics.effective_fg_coverage == 0.0
    assert metrics.fg_coverage > 0.0


def test_it_is_serialized_beside_the_raw_figure(world):
    """A reader comparing two designs needs both numbers, and the report reads
    this dict."""
    payload = build(world, [world["stable"]]).to_dict()

    assert "effective_fg_coverage" in payload
    assert "fg_coverage" in payload
