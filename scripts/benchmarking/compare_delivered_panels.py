"""Run one design twice under two settings and report what moved.

The harness the production-readiness plan needs. Phase 2's consolidation
deletes acceptance rules that sit on the search's stopping condition, so it
moves delivered panels, and this project does not change a delivered panel on
reasoning. Anything of that shape needs a before-and-after on a real design.

## The part that is not obvious

A harness that opens the candidate inventory must compute the SAME reaction
fingerprint the inventory was written under, or the lookup refuses and the run
silently falls back to the supplied list at the list's own frontier -- with no
refills, which is exactly the setting most of these measurements are about.

Reading `params.json` is not enough, and neither is calling `get_params`.
Measured on a Prevotella design on 2026-09-22:

    params.json, fresh process        tm-2026-09-14:2b3536d4305777e8
    after get_params()                tm-2026-09-14:2b3536d4305777e8
    what `filter` actually wrote      tm-2026-09-14:2fb49b12c595c147

The difference is betaine 1.0 M, which the params file never mentions. The
GC-adaptive strategy added it at run time and, as CLAUDE.md records, never
writes back to the file. `effective_conditions` in the run manifest is what
does record it, and is what `export` and `report` already read in preference
to params.json for the same reason.

So this harness reads `effective_conditions` too. Getting that wrong does not
fail loudly: it produces a plausible number from a smaller pool, which is the
shape of the measurement error Known Issue 17 records, where an apparent
density improvement turned out to be the frontier rather than the change under
test.

## Usage

    compare_delivered_panels.py <design_dir> [size] [coverage_target] [density_floor]

`design_dir` must have been through count-kmers, filter and prepare-candidates,
and must still carry the `run_manifest.json` those steps wrote.

Edit `VARIANTS` to say what is being compared. It ships comparing the two
search stopping policies, because that is the measurement outstanding when this
was written.
"""

import json
import pathlib
import sys
import time
from types import SimpleNamespace

import pandas as pd

from neoswga.core import parameter
from neoswga.core import pool_planner as pp
from neoswga.core import search_control as sc
from neoswga.core.base_optimizer import OptimizerConfig
from neoswga.core.optimizer_factory import OptimizerFactory
from neoswga.core.pool_planner import plan_pool
from neoswga.core.position_cache import PositionCache
from neoswga.core.reaction_conditions import build_reaction_conditions
from neoswga.core.run_manifest import read_effective_conditions
from neoswga.core.unified_optimizer import _ensure_optimizers_registered

_ORIGINAL_SEARCH = sc.search_frontiers


def _patch_stop_policy(stop_on_first_qualified):
    """Both names, because `pool_planner` imported the function directly."""

    def patched(*args, **kwargs):
        kwargs["stop_on_first_qualified"] = stop_on_first_qualified
        return _ORIGINAL_SEARCH(*args, **kwargs)

    sc.search_frontiers = pp.search_frontiers = patched


def _restore():
    sc.search_frontiers = pp.search_frontiers = _ORIGINAL_SEARCH


#: name -> setup callable. Each runs before one design and must leave the
#: process able to run the next one.
VARIANTS = {
    "first_feasible": lambda: _patch_stop_policy(True),
    "improve_until_budget": lambda: _patch_stop_policy(False),
}


def conditions_for(design_dir):
    """The reaction the inventory was written under, not the one in the file."""
    effective = read_effective_conditions(str(design_dir))
    if not effective:
        raise SystemExit(
            f"{design_dir}/run_manifest.json records no effective_conditions, so "
            "the reaction the inventory was written under cannot be reproduced "
            "and any refill measurement would silently run over the supplied "
            "list instead. Re-run filter in this directory."
        )
    return build_reaction_conditions(SimpleNamespace(**effective), from_mapping_only=True)


def inventory_source(design_dir, conditions, lengths, frontier):
    """An inventory-backed source, built the way `cli/plan_pool.py` builds one.

    Two mistakes are easy here and both produce a confident number about a
    search that never widened.

    `plan_pool` takes a SOURCE, not a list. Handing it a list wraps it in a
    `ListCandidateSource`, whose `advance()` returns False unconditionally --
    "there is nothing beyond what was handed over" -- so every refill
    measurement reads zero and the run stops with `inventory_exhausted`. The
    CLI passes `InventoryCandidateSource`, and so must anything claiming to
    measure the same thing.

    And the fingerprint must be the one `filter` wrote under, which is what
    `conditions_for` is about.
    """
    from neoswga.core.candidate_source import open_inventory_source

    source = open_inventory_source(
        str(design_dir), conditions.fingerprint(), lengths, frontier=frontier
    )
    if source is None:
        raise SystemExit(f"{design_dir} holds no candidate inventory; run filter there.")
    described = source.describe()
    print(
        f"inventory reachable: {described.get('universe')} candidates, "
        f"frontier opens at {frontier}"
    )
    return source


def run_one(design_dir, params, source, cache, conditions, size, target, floor):
    optimizer = OptimizerFactory.create(
        "hybrid",
        cache,
        params["fg_prefixes"],
        params["fg_seq_lengths"],
        params["bg_prefixes"],
        params["bg_seq_lengths"],
        config=OptimizerConfig(
            max_dimer_bp=params.get("max_dimer_bp", 3),
            max_self_dimer_bp=params.get("max_self_dimer_bp", 4),
            extension_reach=3000,
            fg_circular=params.get("fg_circular", False),
            refinement_method="swap",
            verbose=False,
        ),
        conditions=conditions,
        polymerase=params.get("polymerase", "phi29"),
    )
    started = time.time()
    plan = plan_pool(
        optimizer,
        source,
        [size],
        [target],
        primer_length=params.get("min_k", 12),
        coverage_metric="effective",
        min_selectivity_density=floor,
    )
    row = plan["rows"][0]
    row["_seconds"] = round(time.time() - started, 1)
    return row


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__.split("## Usage")[1].strip())

    design = pathlib.Path(sys.argv[1]).resolve()
    size = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    target = float(sys.argv[3]) if len(sys.argv) > 3 else 0.40
    floor = float(sys.argv[4]) if len(sys.argv) > 4 else 10.0

    params = json.loads((design / "params.json").read_text())

    # The prefixes in params.json are relative to the design directory, and
    # `PositionCache` resolves them against the working directory. Running from
    # anywhere else makes every candidate look unindexed, which the
    # missing-positions guard correctly refuses -- loudly, which is why this is
    # a chdir and not a silently-prefixed path.
    import os

    os.chdir(design)
    parameter.data_dir = "."
    candidates = list(dict.fromkeys(pd.read_csv(design / "step3_df.csv")["primer"]))
    conditions = conditions_for(design)
    lengths = sorted({len(c) for c in candidates})

    source = inventory_source(design, conditions, lengths, len(candidates))

    cache = PositionCache(params["fg_prefixes"] + params["bg_prefixes"], candidates)
    _ensure_optimizers_registered()

    print(
        f"\n{design.name}: pool {len(candidates)}, size {size}, "
        f"coverage target {target}, density floor {floor}\n"
    )

    rows = {}
    for name, setup in VARIANTS.items():
        setup()
        try:
            # A fresh source per variant: the frontier is stateful, and a
            # second run over an already-advanced one would measure the first
            # run's widening rather than its own.
            rows[name] = run_one(
                design,
                params,
                inventory_source(design, conditions, lengths, len(candidates)),
                cache,
                conditions,
                size,
                target,
                floor,
            )
        finally:
            _restore()
        r = rows[name]
        print(
            f"  {name:22} {r['_seconds']:6.1f}s  refills={r.get('frontier_refills')} "
            f"stop={r.get('stop_reason')}"
        )

    print(f"\n{'variant':24}{'coverage':>10}{'density':>10}{'n':>5}")
    for name, r in rows.items():
        print(
            f"{name:24}{r.get('coverage', 0):>10.4f}"
            f"{r.get('selectivity_density', 0):>10.2f}{len(r.get('primers') or []):>5}"
        )

    names = list(rows)
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            pa = set(rows[a].get("primers") or [])
            pb = set(rows[b].get("primers") or [])
            jac = len(pa & pb) / len(pa | pb) if (pa | pb) else 1.0
            verdict = "SAME PANEL" if jac == 1.0 else "PANELS DIFFER"
            print(f"\nJaccard {a} vs {b}: {jac:.3f}  {verdict}")


if __name__ == "__main__":
    main()
