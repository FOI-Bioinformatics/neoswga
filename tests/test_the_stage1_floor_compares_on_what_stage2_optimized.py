"""The Stage-1 floor compares on every axis Stage 2 optimized, not raw coverage.

`_network_refine` reduces the Stage-1 set to the requested size, then floors the
result at the Stage-1 greedy prefix so a destructive refinement cannot return
something worse than just taking the prefix. Until 2026-10-26 "worse" was
strictly-lower raw covered bases -- a single number. That vetoed a refinement
which traded a little coverage for a better melting-temperature profile, which
is exactly what a `stage2_tm_weight` run performs, so the Tm term could not
deliver on a mixed-length pool however it was weighted.

The floor now keeps the prefix only when it PARETO-DOMINATES the refinement on
the axes Stage 2 optimized: the accepted coverage (raw, or the attached
objective's), the objective's shortfall when one is attached, and mean Tm
fitness when the Tm term is active. With one axis this is the `>` it replaced,
bit for bit.

These tests drive `_floor_axes` and `_apply_stage1_floor` directly. The
end-to-end delivered panel is deliberately NOT asserted: hybrid Stage 2 is not
deterministic (unseeded `algebraic_connectivity`, hash-ordered set iteration),
so an end-to-end floor assertion would be flaky. The floor's decision, given a
fixed pair of panels, is deterministic, and that is what is pinned here. The
measured end-to-end effect is in docs/validation/2026-10-26-wiring-the-tm-term.md.
"""

import random

import pytest

h5py = pytest.importorskip("h5py")

import numpy as np

from neoswga.core.hybrid_optimizer import HybridOptimizer
from neoswga.core.position_cache import PositionCache
from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.selection_weights import effective_tm, tm_fitness
from neoswga.core.stage2_ranking import prefix_floors_refinement
from neoswga.core.thermodynamics import reverse_complement

GENOME_LENGTH = 30_000
REACTION_TEMP = 30.0


@pytest.fixture
def genome():
    rng = random.Random(4242)
    seq = list("".join(rng.choice("ACGT") for _ in range(GENOME_LENGTH)))
    primers = []
    for i in range(12):
        primer = "".join(rng.choice("ACGT") for _ in range(10))
        if primer in primers:
            continue
        primers.append(primer)
        for j in range(2 + (i % 3)):
            pos = (i * 2_400 + j * 700) % (GENOME_LENGTH - 20)
            seq[pos : pos + 10] = list(primer)
    return {"seq": "".join(seq), "primers": primers}


@pytest.fixture
def cache(tmp_path, genome):
    prefix = str(tmp_path / "target")
    with h5py.File(f"{prefix}_10mer_positions.h5", "w") as f:
        for primer in genome["primers"]:
            for key in {primer, reverse_complement(primer)}:
                positions, i = [], genome["seq"].find(key)
                while i != -1:
                    positions.append(i)
                    i = genome["seq"].find(key, i + 1)
                if positions:
                    f.create_dataset(key, data=np.array(positions, dtype=np.int32))
    return PositionCache([prefix], genome["primers"]), prefix


def _optimizer(cache, **kwargs):
    position_cache, prefix = cache
    return HybridOptimizer(
        position_cache=position_cache,
        fg_prefixes=[prefix],
        fg_seq_lengths=[GENOME_LENGTH],
        coverage_reach=3_000,
        conditions=ReactionConditions(temp=REACTION_TEMP, polymerase="phi29"),
        reaction_temp=REACTION_TEMP,
        **kwargs,
    )


# ---------------------------------------------------------------------------
# The dominance rule, in the abstract
# ---------------------------------------------------------------------------


def test_one_axis_dominance_is_the_strict_comparison_it_replaced():
    """With a single coverage axis the floor must be bit-identical to `>`.

    This is the property the whole change rests on: a run that declares no
    objective and no Tm term keeps the raw-coverage floor it always had.
    """
    assert prefix_floors_refinement((0.86,), (0.82,)) is True  # prefix higher -> floor
    assert prefix_floors_refinement((0.82,), (0.86,)) is False  # refined higher -> keep
    assert prefix_floors_refinement((0.80,), (0.80,)) is False  # tie -> keep refined


def test_a_coverage_for_tm_trade_is_not_dominated():
    """Prefix better on coverage, refinement better on Tm: incomparable.

    The refinement survives. This is the case the old single-number floor got
    wrong, discarding the refinement for covering less.
    """
    # axes order: (coverage, tm_fitness), both higher-is-better.
    prefix = (0.82, 0.34)
    refined = (0.76, 0.51)
    assert prefix_floors_refinement(prefix, refined) is False


def test_a_refinement_worse_on_every_axis_is_floored():
    prefix = (0.82, 0.51)
    refined = (0.76, 0.34)
    assert prefix_floors_refinement(prefix, refined) is True


# ---------------------------------------------------------------------------
# The axes the optimizer actually builds
# ---------------------------------------------------------------------------


def test_floor_axes_is_raw_coverage_alone_without_objective_or_tm(cache, genome):
    opt = _optimizer(cache)  # no stage2_tm_weight, no pool_objective
    axes = opt._floor_axes(genome["primers"][:4], tm_by_primer={})
    assert len(axes) == 1
    assert axes[0] == pytest.approx(opt._calculate_coverage(genome["primers"][:4]))


def test_floor_axes_gains_a_tm_axis_when_the_term_is_active(cache, genome):
    opt = _optimizer(cache, stage2_tm_weight=0.25)
    panel = genome["primers"][:4]
    tm_by = {p: tm_fitness(effective_tm(p, opt.conditions), REACTION_TEMP) for p in panel}
    axes = opt._floor_axes(panel, tm_by_primer=tm_by)
    assert len(axes) == 2
    assert axes[1] == pytest.approx(sum(tm_by.values()) / len(tm_by))


# ---------------------------------------------------------------------------
# The floor decision, end to end through _apply_stage1_floor, deterministically
# ---------------------------------------------------------------------------


def test_the_floor_keeps_a_tm_refinement_the_old_rule_discarded(cache, genome):
    """A refinement that lowers coverage but improves Tm now survives.

    Two panels are chosen so the prefix covers more and the refinement melts
    closer to the optimum. The old raw-coverage rule returned the prefix; the
    dominance floor returns the refinement.
    """
    opt = _optimizer(cache, stage2_tm_weight=0.5)
    primers = genome["primers"]

    # Rank by Tm fitness; the best-Tm half is the "refinement", a
    # higher-coverage set is the "prefix".
    by_tm = sorted(
        primers, key=lambda p: tm_fitness(effective_tm(p, opt.conditions), REACTION_TEMP)
    )
    refinement = by_tm[-4:]  # best Tm
    prefix = by_tm[:4]  # worst Tm

    tm_by = {p: tm_fitness(effective_tm(p, opt.conditions), REACTION_TEMP) for p in primers}

    # Only exercise the intended case: prefix really does cover at least as much
    # while the refinement melts closer. If the random genome does not produce
    # that arrangement, the test has nothing to assert rather than a false pass.
    cov_prefix = opt._calculate_coverage(prefix)
    cov_refined = opt._calculate_coverage(refinement)

    def dev(ps):
        return sum(abs(opt.conditions.calculate_effective_tm(p) - (REACTION_TEMP + 5)) for p in ps)

    if not (cov_prefix >= cov_refined and dev(refinement) < dev(prefix)):
        pytest.skip("random genome did not produce a coverage-for-Tm trade to test")

    kept = opt._apply_stage1_floor(
        current_primers=list(refinement),
        coverage_floor_set=list(prefix) + list(refinement),
        target_count=4,
        tm_by_primer=tm_by,
        verbose=False,
    )
    assert sorted(kept) == sorted(refinement), "the Tm refinement was floored out"

    # And with the term off, the same lower-coverage set is floored, exactly as
    # before the change.
    off = _optimizer(cache)
    kept_off = off._apply_stage1_floor(
        current_primers=list(refinement),
        coverage_floor_set=list(prefix) + list(refinement),
        target_count=4,
        tm_by_primer={},
        verbose=False,
    )
    if cov_prefix > cov_refined + 1e-12:
        assert sorted(kept_off) == sorted(prefix), "single-axis floor stopped matching the old rule"
