"""The requested set size, and who is allowed to change it.

equiphi29 carries `primer_multiplier = 0.85`, which `hybrid_optimizer` applied
to whatever count it was asked for. Two things were wrong with that.

**It depended on the logging flag.** The assignment `final_count =
adjusted_final_count` sat inside `if adjusted_final_count != final_count and
verbose:`, so the multiplier took effect only when verbose logging was on. The
same parameters produced a different primer set depending on whether the run
was chatty, which is not a property a design should have.

**It overrode an explicit request.** At the default target of 6 the multiplier
is a no-op -- `max(6, int(6 * 0.85))` is 6 -- so the only counts it could ever
change were ones above 7, which is to say counts a caller had deliberately
asked for. A params.json requesting 32 primers got 27, and the completion check
then reported "insufficient candidates", advising the operator to relax the
filters and enlarge the pool. Neither was the problem; the tool had reduced its
own target. That is the same failure acca6cc fixed by another route, and the
same misdirection recorded in docs/validation/additive_specificity.md.

An explicitly configured target now wins, in the sense this file tests: the
tool does not quietly lower it before the search starts.

It is NOT a guarantee of the delivered count. `num_primers` was settled as a
REQUEST on 2026-09-14, so a panel may come back shorter than asked for when the
pool cannot supply enough primers that are mutually compatible at
`max_dimer_bp`. The two multiplier tests below therefore assert on the target
`_rescale_for_polymerase` computes, and only bound the delivered panel from
above. `tests/test_delivered_panel_honours_the_dimer_limit.py` holds the other
half: short is allowed, silently short is not.
"""

import logging
import random

import pytest

h5py = pytest.importorskip("h5py")

import numpy as np

from neoswga.core.hybrid_optimizer import HybridOptimizer
from neoswga.core.position_cache import PositionCache
from neoswga.core.thermodynamics import reverse_complement

GENOME_LENGTH = 30_000


@pytest.fixture(scope="module")
def genome():
    """Candidates that survive equiphi29's own thermodynamic filter.

    The window is 37-62 C. Random 10-mers melt around 37 and are mostly
    rejected, which silently collapses any set to the floor of 6 and makes a
    set-size test prove nothing. 12-mers inside the window are what the real
    equiphi29 designs use.
    """
    from neoswga.core.reaction_conditions import ReactionConditions

    conditions = ReactionConditions(temp=42.0, polymerase="equiphi29", mg_conc=10.0)
    rng = random.Random(4242)
    seq = list("".join(rng.choice("ACGT") for _ in range(GENOME_LENGTH)))

    primers, i = [], 0
    while len(primers) < 20 and i < 5_000:
        i += 1
        primer = "".join(rng.choice("ACGT") for _ in range(12))
        if primer in primers or not 42.0 <= conditions.calculate_effective_tm(primer) <= 58.0:
            continue
        slot = len(primers)
        primers.append(primer)
        for j in range(3):
            pos = (slot * 1_400 + j * 430) % (GENOME_LENGTH - 20)
            seq[pos : pos + 12] = list(primer)
    return {"seq": "".join(seq), "primers": primers}


@pytest.fixture
def optimizer(tmp_path, genome):
    prefix = str(tmp_path / "target")
    with h5py.File(f"{prefix}_12mer_positions.h5", "w") as f:
        for primer in genome["primers"]:
            for key in {primer, reverse_complement(primer)}:
                positions, i = [], genome["seq"].find(key)
                while i != -1:
                    positions.append(i)
                    i = genome["seq"].find(key, i + 1)
                if positions:
                    f.create_dataset(key, data=np.array(positions, dtype=np.int32))
    cache = PositionCache([prefix], genome["primers"])
    return HybridOptimizer(
        position_cache=cache,
        fg_prefixes=[prefix],
        fg_seq_lengths=[GENOME_LENGTH],
        bin_size=1_000,
        coverage_reach=3_000,
        polymerase="equiphi29",
    )


def test_multiplier_does_not_depend_on_the_logging_flag(optimizer, genome):
    """A design must not change because logging was switched on."""
    quiet = optimizer.optimize(genome["primers"], final_count=10, verbose=False)
    loud = optimizer.optimize(genome["primers"], final_count=10, verbose=True)

    # Before the fix: 10 quiet, 8 loud. phi29 (multiplier 1.0) matched at both.
    assert len(quiet.primers) == len(loud.primers)


def test_an_explicit_target_is_not_rescaled(optimizer, genome):
    """The tool must not reduce a target the caller named.

    Asserted on the target rather than on the delivered panel. Since dimer
    relaxation defaults off, the greedy stops rather than admitting a pair
    above `max_dimer_bp`, so this fixture delivers fewer than ten whatever the
    target is -- which would make a delivered-count assertion pass for the
    wrong reason, or fail for one. `--num-primers` is a request (decided
    2026-09-14); what this test forbids is the tool quietly lowering the
    request before it even starts.
    """
    assert optimizer._rescale_for_polymerase(10, verbose=False) == 8  # the knob works
    result = optimizer.optimize(
        genome["primers"],
        final_count=10,
        verbose=False,
        apply_polymerase_multiplier=False,
    )

    # Not rescaled to 8 behind the caller's back, and never above the request.
    assert 0 < len(result.primers) <= 10


def test_the_multiplier_still_applies_when_asked_for(optimizer, genome):
    """The knob is preserved for callers that want the polymerase heuristic.

    Also asserted on the target: equiphi29's 0.85 takes 10 to 8, and the
    delivered panel may then be shorter still for pool reasons.
    """
    assert optimizer._rescale_for_polymerase(10, verbose=False) == 8

    result = optimizer.optimize(
        genome["primers"],
        final_count=10,
        verbose=False,
        apply_polymerase_multiplier=True,
    )

    assert 0 < len(result.primers) <= 8


def test_the_multiplier_never_raises_a_target(optimizer):
    """`max(6, ...)` is a floor, so a small request is left alone."""
    assert optimizer._rescale_for_polymerase(6, verbose=False) == 6
    assert optimizer._rescale_for_polymerase(4, verbose=False) == 6


def test_the_adapter_does_not_rescale_by_polymerase(tmp_path, genome, caplog):
    """End of the path a params.json actually travels.

    `num_primers` reaches the optimizer through `OptimizerConfig`, so this is
    the level at which "asked for 32, got 27" was visible to a user.

    Asserted on the adjustment itself rather than on a resulting count. Two
    other things legitimately change how many primers come back -- Stage 1 can
    saturate coverage before reaching the target, and each polymerase applies
    its own Tm window to the candidate pool -- so a count-based assertion here
    would be testing the fixture. What must not happen is the rescaling.
    """
    from neoswga.core.base_optimizer import OptimizerConfig
    from neoswga.core.hybrid_optimizer import HybridBaseOptimizer

    prefix = str(tmp_path / "adapter")
    with h5py.File(f"{prefix}_12mer_positions.h5", "w") as f:
        for primer in genome["primers"]:
            for key in {primer, reverse_complement(primer)}:
                positions, i = [], genome["seq"].find(key)
                while i != -1:
                    positions.append(i)
                    i = genome["seq"].find(key, i + 1)
                if positions:
                    f.create_dataset(key, data=np.array(positions, dtype=np.int32))

    optimizer = HybridBaseOptimizer(
        position_cache=PositionCache([prefix], genome["primers"]),
        fg_prefixes=[prefix],
        fg_seq_lengths=[GENOME_LENGTH],
        config=OptimizerConfig(target_set_size=10, verbose=True),
        polymerase="equiphi29",
    )

    with caplog.at_level(logging.INFO, logger="neoswga.core.hybrid_optimizer"):
        optimizer.optimize(genome["primers"])

    assert "Adjusted target" not in caplog.text
    assert "Target: 10 final primers" in caplog.text


# ----------------------------------------------------------------------
# What the completion check is allowed to claim
# ----------------------------------------------------------------------


def test_shortfall_advice_does_not_assert_a_cause_it_has_not_established():
    """The check sees two numbers. It cannot see why they differ.

    It used to report "PARTIAL result: insufficient candidates" and offer four
    remedies, all of them about enlarging the candidate pool. A shortfall has
    at least one other common cause that the advice actively misdirects from:
    Stage 1 reaching full coverage before it reaches the target, which is the
    optimizer succeeding. Acting on the advice there loosens filters that were
    not the problem, which is how this codebase has been wrong before.
    """
    from neoswga.cli.pipeline import set_size_shortfall_advice

    lines = set_size_shortfall_advice(num_found=8, target_size=10)
    text = " ".join(lines)

    assert "insufficient candidates" not in text
    assert "8" in text and "10" in text
    # Saturation is the benign explanation and must be offered first.
    assert "coverage" in text.lower()


def test_shortfall_advice_still_offers_the_pool_remedies():
    """Enlarging the pool remains a real fix; it just is not the only one."""
    from neoswga.cli.pipeline import set_size_shortfall_advice

    text = " ".join(set_size_shortfall_advice(num_found=8, target_size=10))

    assert "max_primer" in text
