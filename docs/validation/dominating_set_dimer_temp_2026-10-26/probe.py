"""Probe: the dominating-set free-energy dimer screen ran at 37 C, not the reaction.

Run from the repo root:

    PYTHONPATH=. python docs/validation/dominating_set_dimer_temp_2026-10-26/probe.py

Two parts, neither needing an external tool or a genome:

1. The consequence. For pairs that PASS the base-pair screen (so the free-energy
   floor is the only thing left to reject them), the verdict genuinely differs
   between the phi29 default 30 C and the hardcoded 37 C -- always in the
   PERMISSIVE direction, so screening at 37 C lets through pairs that dimerise
   in the actual reaction.
2. The fix. The reaction temperature now reaches
   `DominatingSetOptimizer.dimer_temp` through the adapter, and with no
   conditions attached the inner's own 37 C default still stands.
"""

import random

import numpy as np

from neoswga.core.dimer import is_dimer_fast, is_dimer_thermodynamic

MAX_DIMER_BP = 7
FLOORS = (-1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0)
TRIALS = 30_000


def flip_envelope():
    """How many random pairs change verdict between 30 C and 37 C, by floor."""
    random.seed(7)
    bases = "ACGT"
    envelope = {}
    permissive_only = True
    for _ in range(TRIALS):
        a = "".join(random.choice(bases) for _ in range(14))
        b = "".join(random.choice(bases) for _ in range(14))
        if is_dimer_fast(a, b, max_dimer_bp=MAX_DIMER_BP):
            continue  # already rejected on length; the floor never runs
        for floor in FLOORS:
            at_30 = is_dimer_thermodynamic(a, b, delta_g_threshold=floor, temperature=30.0)
            at_37 = is_dimer_thermodynamic(a, b, delta_g_threshold=floor, temperature=37.0)
            if at_30 != at_37:
                envelope.setdefault(floor, []).append((a, b, at_30, at_37))
                if not (at_30 and not at_37):
                    permissive_only = False
    return envelope, permissive_only


def forwarding():
    from neoswga.core import unified_optimizer as uo
    from neoswga.core.base_optimizer import OptimizerConfig
    from neoswga.core.optimizer_factory import OptimizerFactory
    from neoswga.core.reaction_conditions import ReactionConditions

    uo._ensure_optimizers_registered()

    class _Cache:
        def get_positions(self, *args, **kwargs):
            return np.array([], dtype=np.int64)

    def build(conditions):
        return OptimizerFactory.create(
            name="dominating-set",
            position_cache=_Cache(),
            fg_prefixes=["fg"],
            fg_seq_lengths=[10_000],
            bg_prefixes=[],
            bg_seq_lengths=[],
            config=OptimizerConfig(target_set_size=6, max_dimer_bp=3, max_dimer_dg=-6.0),
            conditions=conditions,
        )

    out = {}
    for temp, poly in ((30.0, "phi29"), (65.0, "bst")):
        out[temp] = build(ReactionConditions(polymerase=poly, temp=temp))._optimizer.dimer_temp
    out[None] = build(None)._optimizer.dimer_temp
    return out


def main():
    print("1. Does the screen temperature change the verdict?")
    print(f"   {TRIALS} random 14-mer pairs, keeping those that pass max_dimer_bp={MAX_DIMER_BP}")
    print("   (only those reach the free-energy floor at all).\n")
    envelope, permissive_only = flip_envelope()
    print("   floor (kcal/mol)   pairs whose 30 C and 37 C verdicts differ")
    for floor in sorted(envelope, reverse=True):
        rows = envelope[floor]
        a, b, at_30, at_37 = rows[0]
        print(
            f"     {floor:5.1f}            {len(rows):5d}   e.g. {a}/{b}  30C={at_30} 37C={at_37}"
        )
    assert envelope, "expected at least one flip"
    assert permissive_only, "expected every flip to be dimer-at-30 / not-at-37"
    print("\n   Every flip is dimer at 30 C and NOT a dimer at 37 C: screening at")
    print("   37 C under a 30 C reaction is the permissive direction.")
    print(f"   At the Rychlik -6.0 guideline floor: {len(envelope[-6.0])} pairs.\n")

    print("2. Does the reaction temperature reach the inner screen?")
    out = forwarding()
    print(f"   conditions.temp=30.0 (phi29) -> inner.dimer_temp = {out[30.0]}")
    print(f"   conditions.temp=65.0 (bst)   -> inner.dimer_temp = {out[65.0]}")
    print(
        f"   conditions=None              -> inner.dimer_temp = {out[None]} (inner default stands)"
    )
    assert out[30.0] == 30.0 and out[65.0] == 65.0, out
    assert out[None] == 37.0, out
    print("\nOK: the screen follows the reaction, and nothing moves without conditions.")


if __name__ == "__main__":
    main()
