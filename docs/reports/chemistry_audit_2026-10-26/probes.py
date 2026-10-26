#!/usr/bin/env python3
"""Probes behind docs/reports/AUDIT_2026-10-26_chemistry_in_pool_design.md.

Genome-free. Every probe here exercises the chemistry model directly, so it
runs without a k-mer counter, without an HDF5 index and without a reference
genome. That is deliberate: the questions it answers are about which knobs
reach which calculation, and none of them needs a genome to settle.

What it does NOT do: anything about a delivered panel. Whether a chemistry
change moves the selected oligos needs the candidate pool and the position
indexes, and those are not in a fresh clone. Probe 5 states that boundary
rather than guessing past it.

Run:  PYTHONPATH=. python docs/reports/chemistry_audit_2026-10-26/probes.py
"""

from __future__ import annotations

import inspect
import json
import sys

from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.secondary_structure import check_homodimer

BASE = dict(temp=42.0, polymerase="equiphi29")

OLIGOS = ["GCGCGGTTCTAC", "AACCTTGGAACC", "CGCGCGCGATAT", "GGCCAATTCCGG"]

# One typical protocol dose per knob the JSON schema accepts.
DOSES: dict[str, object] = {
    "dmso_percent": 5.0,
    "betaine_m": 1.0,
    "trehalose_m": 0.6,
    "formamide_percent": 2.0,
    "glycerol_percent": 5.0,
    "bsa_ug_ml": 100.0,
    "peg_percent": 5.0,
    "ethanol_percent": 3.0,
    "urea_m": 0.5,
    "tmac_m": 0.05,
    "propanediol_m": 1.0,
    "dtt_mm": 5.0,
    "ssb": True,
}

SALTS = {"na_conc": 100.0, "k_conc": 50.0, "nh4_conc": 10.0, "dntp_conc": 1.0}


def effective_tm(conditions: ReactionConditions, seq: str) -> float:
    """Call `calculate_effective_tm` whatever its parameter happens to be named."""
    fn = conditions.calculate_effective_tm
    params = inspect.signature(fn).parameters
    for name in ("primer", "sequence", "seq"):
        if name in params:
            return fn(**{name: seq})
    return fn(seq)


def mean_tm(conditions: ReactionConditions) -> float:
    return sum(effective_tm(conditions, s) for s in OLIGOS) / len(OLIGOS)


def probe_1_which_knobs_move_tm() -> dict:
    """Which schema additives move an effective Tm, and which move nothing.

    A knob that moves nothing is not thereby wrong -- BSA and SSB have no
    reason to shift a duplex Tm -- but it is accepted, range-validated and
    printed in the condition summary, so a reader cannot tell the difference
    between "modelled as zero" and "not modelled".
    """
    base = ReactionConditions(**BASE)
    reference = [effective_tm(base, s) for s in OLIGOS]
    out = {}
    for knob, dose in {**DOSES, **SALTS}.items():
        try:
            cond = ReactionConditions(**BASE, **{knob: dose})
        except Exception as exc:  # noqa: BLE001 - the refusal is the result
            out[knob] = {"dose": dose, "delta_tm": None, "verdict": f"refused: {type(exc).__name__}"}
            continue
        delta = max(abs(effective_tm(cond, s) - r) for s, r in zip(OLIGOS, reference))
        out[knob] = {
            "dose": dose,
            "delta_tm": delta,
            "verdict": "moves Tm" if delta > 1e-9 else "no Tm term",
        }
    return out


def probe_2_salt_and_concentration() -> dict:
    """How much Tm leverage the buffer and the per-oligo concentration carry.

    `mg_conc` defaults to None and resolves from the polymerase preset, so a
    sweep has to pass explicit values; passing the preset value and seeing no
    change is a default collision, not an absent term.
    """
    return {
        "resolved_mg_conc_when_none": ReactionConditions(**BASE).mg_conc,
        "mg_sweep_mM": {
            mg: mean_tm(ReactionConditions(**BASE, mg_conc=mg))
            for mg in (0.0, 2.0, 5.0, 10.0, 20.0)
        },
        "primer_conc_sweep_M": {
            f"{pc:.0e}": mean_tm(ReactionConditions(**BASE, primer_conc=pc))
            for pc in (5e-8, 5e-7, 5e-6, 5e-5)
        },
    }


def probe_3_tm_is_temperature_invariant() -> dict:
    """Effective Tm must not move with the reaction temperature.

    Tm is a property of the duplex and the buffer. Temperature acts through
    occupancy, not through Tm, and a Tm that moved with `temp` would be double
    counting. This is a guard on the probe above rather than a finding.
    """
    return {t: mean_tm(ReactionConditions(temp=t, polymerase="equiphi29")) for t in (30.0, 42.0, 45.0)}


def probe_4_dimer_model_sees_salt_not_additives() -> dict:
    """Does the primer-primer duplex model carry the additive correction?

    The additive correction is applied to the primer-TEMPLATE duplex. If it is
    not applied to the primer-PRIMER duplex, a run with DMSO and betaine
    destabilises template binding in the model and leaves dimers at their
    no-additive stability, which biases the screen that rejects pairs.
    """
    cases = {
        "baseline": {},
        "mg_0mM": dict(mg_conc=0.0),
        "mg_20mM": dict(mg_conc=20.0),
        "dmso_10pct": dict(dmso_percent=10.0),
        "betaine_2M": dict(betaine_m=2.0),
        "dmso_10pct_betaine_2M": dict(dmso_percent=10.0, betaine_m=2.0),
    }
    out = {}
    for seq in ("GCGCGCGCGCGC", "GCGCGGTTCTAC"):
        rows = {}
        for name, kw in cases.items():
            result = check_homodimer(seq, ReactionConditions(**BASE, **kw))
            rows[name] = {
                "energy_kcal_mol": float(result["energy"]),
                "tm_c": float(result["tm"]),
                "binding_length_bp": int(result["binding_length"]),
                "forms_dimer": bool(result["forms_dimer"]),
            }
        out[seq] = rows
    return out


def probe_5_what_this_clone_cannot_answer() -> dict:
    """State the boundary instead of guessing past it."""
    import importlib.util
    import shutil

    return {
        "h5py_importable": importlib.util.find_spec("h5py") is not None,
        "kmc_on_path": shutil.which("kmc") is not None,
        "jellyfish_on_path": shutil.which("jellyfish") is not None,
        "note": (
            "Whether a chemistry change moves the DELIVERED panel needs step3_df.csv "
            "and the position indexes. Neither is in a fresh clone, and no reference "
            "genome is either, so that question is open here rather than answered."
        ),
    }


def main() -> int:
    results = {
        "probe_1_which_knobs_move_tm": probe_1_which_knobs_move_tm(),
        "probe_2_salt_and_concentration": probe_2_salt_and_concentration(),
        "probe_3_tm_is_temperature_invariant": probe_3_tm_is_temperature_invariant(),
        "probe_4_dimer_model_sees_salt_not_additives": probe_4_dimer_model_sees_salt_not_additives(),
        "probe_5_what_this_clone_cannot_answer": probe_5_what_this_clone_cannot_answer(),
    }

    p1 = results["probe_1_which_knobs_move_tm"]
    print("Probe 1 -- which schema knobs move an effective Tm (equiphi29, 42 C)\n")
    print(f"{'knob':<22}{'dose':>8}{'max |dTm| (C)':>16}  verdict")
    print("-" * 68)
    for knob, row in p1.items():
        delta = "n/a" if row["delta_tm"] is None else f"{row['delta_tm']:.4f}"
        print(f"{knob:<22}{str(row['dose']):>8}{delta:>16}  {row['verdict']}")
    inert = [k for k, v in p1.items() if v["verdict"] == "no Tm term"]
    print(f"\nNo Tm term ({len(inert)}): {', '.join(inert) or 'none'}")

    p2 = results["probe_2_salt_and_concentration"]
    print(f"\nProbe 2 -- mg_conc resolves to {p2['resolved_mg_conc_when_none']} mM when unset")
    for mg, tm in p2["mg_sweep_mM"].items():
        print(f"   Mg {mg:>5} mM -> mean effective Tm {tm:8.4f} C")
    for pc, tm in p2["primer_conc_sweep_M"].items():
        print(f"   primer {pc:>8} M -> mean effective Tm {tm:8.4f} C")

    print("\nProbe 3 -- effective Tm against reaction temperature (must be flat)")
    for t, tm in results["probe_3_tm_is_temperature_invariant"].items():
        print(f"   temp {t:>5} C -> mean effective Tm {tm:8.4f} C")

    print("\nProbe 4 -- homodimer model against salt and additives")
    for seq, rows in results["probe_4_dimer_model_sees_salt_not_additives"].items():
        print(f"   {seq}")
        for name, row in rows.items():
            print(
                f"      {name:<24} energy {row['energy_kcal_mol']:>10.4f}"
                f"   tm {row['tm_c']:>10.4f}   bp {row['binding_length_bp']}"
            )

    print("\nProbe 5 -- boundary of this clone")
    for key, value in results["probe_5_what_this_clone_cannot_answer"].items():
        print(f"   {key}: {value}")

    out = __file__.replace("probes.py", "probe_results.json")
    with open(out, "w") as handle:
        json.dump(results, handle, indent=2, default=str)
    print(f"\nWrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
