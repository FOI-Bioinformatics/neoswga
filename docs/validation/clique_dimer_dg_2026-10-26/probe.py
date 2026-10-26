"""Probe: clique applied no free-energy floor, and now does.

Run from the repo root:

    PYTHONPATH=. python docs/validation/clique_dimer_dg_2026-10-26/probe.py

Two parts, neither needing a k-mer counter or a genome:

1. The floor reaches the optimizer, and the temperature it is evaluated at
   comes off the attached ReactionConditions.
2. What it does to the compatibility graph: it only ever removes edges the
   base-pair screen had already cleared, and how many depends on the reaction
   temperature -- the same floor removes edges at 30 C and far fewer at 65 C.
"""

import itertools
import logging
import random

import numpy as np

from neoswga.core.clique_optimizer import CliqueOptimizer, build_compatibility_graph
from neoswga.core.dimer import is_dimer_fast, is_dimer_thermodynamic

# The fixture the regression tests use: across this pool nothing exceeds
# max_dimer_bp=3 and nothing self-dimers at 4, and exactly one pair falls below
# a -1.0 kcal/mol floor at 37 C. So the floor is the only discriminator.
FIXED = "GGATCACAGTCT"
PARTNER = "TTCCAGATAGCG"
BENIGN = ["ACATTTTTTCGT", "CTCACAATAGAA", "TTTGAATACTCA"]
FLOOR = -1.0


class _Cache:
    def get_positions(self, *args, **kwargs):
        return np.array([], dtype=np.int64)


def _clique(max_dimer_dg, conditions=None):
    from neoswga.core.base_optimizer import OptimizerConfig

    return CliqueOptimizer(
        position_cache=_Cache(),
        fg_prefixes=["x"],
        fg_seq_lengths=[1000],
        config=OptimizerConfig(target_set_size=4, max_dimer_bp=3, max_dimer_dg=max_dimer_dg),
        conditions=conditions,
    )


def delivered(max_dimer_dg):
    optimizer = _clique(max_dimer_dg)
    # Make the offending candidate the most attractive one, so only the screen
    # can keep it out.
    optimizer._score_set = lambda s: 1.0 if PARTNER in s else 0.0
    result = optimizer.optimize(candidates=[*BENIGN, PARTNER], num_primers=4, fixed_primers=[FIXED])
    return list(result.primers)


def main():
    logging.disable(logging.CRITICAL)

    print("The fixture: no pair exceeds max_dimer_bp=3, and exactly one falls")
    print(f"below a {FLOOR} kcal/mol floor at 37 C.\n")
    pool = [*BENIGN, PARTNER, FIXED]
    bp_pairs = [(a, b) for a, b in itertools.combinations(pool, 2) if is_dimer_fast(a, b, 3)]
    dg_pairs = [
        (a, b)
        for a, b in itertools.combinations(pool, 2)
        if is_dimer_thermodynamic(a, b, delta_g_threshold=FLOOR, temperature=37.0)
    ]
    print(f"   pairs above max_dimer_bp=3 : {bp_pairs}")
    print(f"   pairs below the floor      : {dg_pairs}\n")
    assert not bp_pairs and dg_pairs == [(PARTNER, FIXED)]

    print("1. Is the offending candidate delivered?\n")
    without = delivered(None)
    with_floor = delivered(FLOOR)
    print(f"   max_dimer_dg unset  -> {PARTNER} delivered: {PARTNER in without}  (the control)")
    print(f"   max_dimer_dg {FLOOR}  -> {PARTNER} delivered: {PARTNER in with_floor}")
    assert PARTNER in without, "the no-floor control must still deliver it"
    assert PARTNER not in with_floor, "the configured floor must exclude it"
    print("\n   Before this fix the second row also read True: clique read")
    print("   max_dimer_dg nowhere, so the floor excluded nothing.\n")

    print("2. Does the temperature the floor is evaluated at reach it?\n")
    from neoswga.core.reaction_conditions import ReactionConditions

    for polymerase, temp in (("phi29", 30.0), ("bst", 65.0)):
        got = _clique(FLOOR, ReactionConditions(polymerase=polymerase, temp=temp))._screen_inputs()
        print(f"   conditions {polymerase:6} {temp:5.1f} C -> _screen_inputs = {got}")
        assert got == (FLOOR, temp)

    print("\n3. The floor only ever removes edges, and temperature decides how many:\n")
    random.seed(3)
    primers = list(
        dict.fromkeys("".join(random.choice("ACGT") for _ in range(14)) for _ in range(60))
    )

    def edges(max_dimer_dg, temp):
        return build_compatibility_graph(
            primers,
            max_dimer_bp=7,
            max_self_dimer_bp=8,
            check_self_dimers=True,
            max_dimer_dg=max_dimer_dg,
            temp=temp,
        ).number_of_edges()

    base = edges(None, 30.0)
    print(f"   no floor                      : {base} edges")
    for floor in (-3.0, -6.0):
        for temp in (30.0, 65.0):
            count = edges(floor, temp)
            print(
                f"   floor {floor:5.1f} at {temp:5.1f} C        : {count} edges (-{base - count})"
            )
            assert count <= base, "a floor may only tighten"

    print("\nOK: the floor is opt-in, reaches the screen, follows the reaction, and only tightens.")


if __name__ == "__main__":
    main()
