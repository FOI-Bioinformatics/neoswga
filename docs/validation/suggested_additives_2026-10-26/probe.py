#!/usr/bin/env python3
"""Do the additives the tool suggests reach a design decision?

Behind docs/validation/2026-10-26-suggested-additives-reach-design.md. Continues
the chemistry audit (docs/reports/AUDIT_2026-10-26_chemistry_in_pool_design.md),
narrowing to: for every additive the tool accepts, and especially the ones its
suggester and optimizer RECOMMEND, is the effect actually used in design?

"Used in design" here means one of two reachable paths, both genome-free to
test:
  * effective Tm moves  -> the filter Tm gate (candidate eligibility) and the
    occupancy weighting (selectivity/coverage scoring) both see it;
  * require_model_support accepts the design -> the additive is not refused.

An additive that moves no Tm AND is accepted silently reaches neither; one that
is recommended but refused cannot be used at all. Those are the findings.

Run:  PYTHONPATH=. python docs/validation/suggested_additives_2026-10-26/probe.py
"""

from __future__ import annotations

import inspect
import json
import types

from neoswga.core.model_evidence import load_evidence, require_model_support
from neoswga.core.reaction_conditions import ReactionConditions

# Who recommends each additive as a tunable output.
#   suggester = condition_suggester.suggest_conditions (betaine/dmso/trehalose)
#   optimizer = additive_optimizer.AdditiveOptimizer   (+ formamide, glycerol,
#               mg, and an ssb flag on its recommendation)
SUGGESTED_BY = {
    "dmso_percent": ("suggester", "optimizer"),
    "betaine_m": ("suggester", "optimizer"),
    "trehalose_m": ("suggester", "optimizer"),
    "formamide_percent": ("optimizer",),
    "glycerol_percent": ("optimizer",),
    "ssb": ("optimizer",),
}

# A typical protocol dose per additive.
DOSES = {
    "dmso_percent": 5.0,
    "betaine_m": 1.0,
    "trehalose_m": 0.6,
    "formamide_percent": 2.0,
    "glycerol_percent": 5.0,
    "ethanol_percent": 3.0,
    "urea_m": 0.5,
    "tmac_m": 0.05,
    "propanediol_m": 1.0,
    "bsa_ug_ml": 100.0,
    "peg_percent": 5.0,
    "dtt_mm": 5.0,
    "ssb": True,
}

OLIGOS = ["GCGCGGTTCTAC", "AACCTTGGAACC", "CGCGCGCGATAT", "GGCCAATTCCGG"]


def effective_tm(cond, seq):
    fn = cond.calculate_effective_tm
    params = inspect.signature(fn).parameters
    for name in ("primer", "sequence", "seq"):
        if name in params:
            return fn(**{name: seq})
    return fn(seq)


def tm_delta(field, dose):
    base = ReactionConditions(temp=42.0, polymerase="equiphi29")
    ref = [effective_tm(base, s) for s in OLIGOS]
    cond = ReactionConditions(temp=42.0, polymerase="equiphi29", **{field: dose})
    return max(abs(effective_tm(cond, s) - r) for s, r in zip(OLIGOS, ref, strict=True))


def design_outcome(field, dose):
    """What require_model_support does with a design that sets this additive."""
    cond = ReactionConditions(temp=42.0, polymerase="equiphi29", **{field: dose})
    request = types.SimpleNamespace(conditions=cond, polymerase="equiphi29", primer_lengths=(12,))
    try:
        require_model_support(request)
    except Exception as exc:  # noqa: BLE001 - the refusal is the datum
        return f"refused: {type(exc).__name__}"
    # Accepted. Was it even checked? The additive-support map is the gate.
    from neoswga.core.model_evidence import _ADDITIVE_EVIDENCE

    return "accepted" if field in _ADDITIVE_EVIDENCE else "accepted, never checked"


def main() -> int:
    evidence = load_evidence()
    from neoswga.core.model_evidence import _ADDITIVE_EVIDENCE

    rows = []
    for field, dose in DOSES.items():
        delta = tm_delta(field, dose)
        quantity = _ADDITIVE_EVIDENCE.get(field)
        record = evidence.get(quantity) if quantity else None
        rows.append(
            {
                "additive": field,
                "suggested_by": SUGGESTED_BY.get(field, ()),
                "tm_delta_c": round(delta, 4),
                "moves_tm": delta > 1e-9,
                "design_outcome": design_outcome(field, dose),
                "evidence_status": getattr(record, "status", None) if record else None,
            }
        )

    print(f"{'additive':<18}{'suggested':<22}{'dTm (C)':>9}  {'design outcome':<24}{'evidence'}")
    print("-" * 90)
    for r in rows:
        sug = "+".join(r["suggested_by"]) or "-"
        print(
            f"{r['additive']:<18}{sug:<22}{r['tm_delta_c']:>9.4f}  "
            f"{r['design_outcome']:<24}{r['evidence_status'] or '-'}"
        )

    print("\nSuggested additives that do NOT reach a design decision:")
    for r in rows:
        if not r["suggested_by"]:
            continue
        reaches = r["moves_tm"] and not r["design_outcome"].startswith("refused")
        if not reaches:
            why = "refused" if r["design_outcome"].startswith("refused") else "no Tm effect"
            print(f"  {r['additive']:<18} suggested by {'+'.join(r['suggested_by'])}: {why}")

    out = __file__.replace("probe.py", "probe_results.json")
    with open(out, "w") as handle:
        json.dump(rows, handle, indent=2, default=str)
    print(f"\nWrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
