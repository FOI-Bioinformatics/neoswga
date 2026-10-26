"""Probe: hybrid ran every reaction as phi29 at 30 C.

Run from the repo root:

    PYTHONPATH=. python docs/validation/hybrid_reaction_2026-10-26/probe.py

Three parts, none needing a k-mer counter or a genome:

1. What the inner HybridOptimizer resolved to, per configured polymerase,
   through the real factory dispatch and for BOTH hybrid wrappers.
2. The gate that rode on it: `poly_config.thermo_filter` decides whether
   `optimize` runs its Stage-0 thermodynamic pre-screen at all, and phi29's
   preset sets it False while equiphi29's and bst's set it True.
3. What that pre-screen removes when it is allowed to run, under the
   configuration of the SHIPPED `examples/equiphi29_scenario` params.json.
"""

import random

import numpy as np

from neoswga.core.hybrid_optimizer import HybridOptimizer
from neoswga.core.reaction_conditions import ReactionConditions

REACTIONS = (("phi29", 30.0), ("equiphi29", 45.0), ("bst", 65.0))


class _Cache:
    def get_positions(self, *args, **kwargs):
        return np.array([], dtype=np.int64)


def resolved():
    from neoswga.core import unified_optimizer as uo
    from neoswga.core.base_optimizer import OptimizerConfig
    from neoswga.core.optimizer_factory import OptimizerFactory

    uo._ensure_optimizers_registered()
    rows = []
    for method in ("hybrid", "background-aware"):
        for polymerase, temp in REACTIONS:
            optimizer = OptimizerFactory.create(
                name=method,
                position_cache=_Cache(),
                fg_prefixes=["fg"],
                fg_seq_lengths=[10_000],
                bg_prefixes=[],
                bg_seq_lengths=[],
                config=OptimizerConfig(target_set_size=6),
                conditions=ReactionConditions(polymerase=polymerase, temp=temp),
            )
            inner = optimizer._hybrid
            rows.append(
                (
                    method,
                    polymerase,
                    temp,
                    inner.polymerase,
                    inner.reaction_temp,
                    inner.poly_config.thermo_filter,
                )
            )
    return rows


def prescreen_cost():
    """What the Stage-0 pre-screen removes, under the shipped equiphi29 example.

    `optimize` gates the pre-screen on `poly_config.thermo_filter` (see the
    `PRE-STAGE: Thermodynamic Filtering` branch). Resolved as phi29 that flag is
    False, so the whole pre-screen was skipped and every candidate reached
    selection unscreened. This runs the same filter the gate guards.
    """
    random.seed(11)
    candidates = list(
        dict.fromkeys(
            "".join(random.choice("ACGT") for _ in range(random.choice([10, 12, 14, 16, 18])))
            for _ in range(500)
        )
    )
    # The shipped examples/equiphi29_scenario/params.json reaction.
    conditions = ReactionConditions(
        polymerase="equiphi29",
        temp=43.0,
        na_conc=50.0,
        mg_conc=10.0,
        dmso_percent=5.0,
        betaine_m=1.0,
    )
    optimizer = HybridOptimizer(
        position_cache=_Cache(),
        fg_prefixes=["fg"],
        fg_seq_lengths=[10_000],
        min_tm=40,
        max_tm=65,
        conditions=conditions,
    )
    kept = optimizer._thermo_filter_candidates(list(candidates), verbose=False)
    return len(candidates), len(kept), optimizer.poly_config.thermo_filter


def main():
    print("1. What the inner resolved to (real factory dispatch, both wrappers):\n")
    print(
        f"   {'method':17} {'configured':11} {'temp':>6}   {'-> inner':11} {'temp':>6}  thermo_filter"
    )
    for method, polymerase, temp, got_poly, got_temp, flag in resolved():
        print(f"   {method:17} {polymerase:11} {temp:6.1f}   {got_poly:11} {got_temp!s:>6}  {flag}")
    print()
    print("   Before this fix every row read `phi29  30.0  False`, because")
    print("   unified_optimizer passes `conditions` and passes neither")
    print("   `polymerase` nor `reaction_temp`, while both wrappers defaulted")
    print('   polymerase to the literal "phi29".\n')

    total, kept, flag = prescreen_cost()
    print("2+3. The Stage-0 pre-screen, under the shipped equiphi29 example")
    print("     (polymerase equiphi29, reaction_temp 43 C, Tm window 40-65):\n")
    print(f"     equiphi29 preset thermo_filter = {flag} -> `optimize` runs the pre-screen")
    print(f"     pre-screen keeps {kept}/{total} candidates")
    print(f"     resolved as phi29, thermo_filter=False -> pre-screen SKIPPED,")
    print(f"     so {total}/{total} reached selection unscreened ({total - kept} of them")
    print("     ones the configured reaction rejects).")
    assert flag is True
    assert kept < total


if __name__ == "__main__":
    main()
