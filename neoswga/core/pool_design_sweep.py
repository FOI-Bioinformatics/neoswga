"""Compare chemistries by redesigning the pool, not by rescoring one panel.

Task 7 of the condition-aware pool design plan.

`condition_sweep.sweep_conditions` scores ONE fixed panel under several
reactions. That answers "how does this panel behave if I change the buffer",
which is a useful question and not the one a design comparison asks.

A panel chosen under one chemistry carries that chemistry's decisions: which
candidates passed the Tm window, which the ranking favoured, which the dimer
screen removed. Rescoring it under a different reaction measures the transplant
rather than the alternative, and a chemistry can look better simply by suiting
the panel that was already chosen.

`design_sweep` runs a separate design per condition and per oligo length, each
drawing from the shared candidate inventory through its own provider, and then
compares the delivered designs. Eligibility is per condition, so a candidate
admitted under one reaction and not another appears in one design and not the
other. That is the point, not an inconsistency.

Panels stay single-length. Mixed-length pools are a separate feature: they
change what a dimer screen and a Tm window mean across the panel, and pretending
otherwise would hide that behind a frontier.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import Any

from .candidate_provider import CandidateProvider
from .pool_objective import PoolConstraints


def nondominated(designs: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    """The designs no other design beats on every axis at once.

    Three axes, all "smaller or larger is better" in a fixed direction: fewer
    oligos, higher coverage, lower background load. A design is dominated when
    another is at least as good on all three and strictly better on one.

    Equal designs do not dominate each other, so a tie keeps both rather than
    silently preferring whichever was evaluated first.
    """

    def dominates(a, b):
        not_worse = (
            a["size"] <= b["size"]
            and a["coverage"] >= b["coverage"]
            and a["background"] <= b["background"]
        )
        strictly_better = (
            a["size"] < b["size"]
            or a["coverage"] > b["coverage"]
            or a["background"] < b["background"]
        )
        return not_worse and strictly_better

    return [d for d in designs if not any(dominates(other, d) for other in designs)]


def design_sweep(
    inventory,
    conditions: Sequence[Any],
    lengths: Sequence[int],
    sizes: Sequence[int],
    coverage_targets: Sequence[float],
    constraints: PoolConstraints,
    run_design: Callable[..., dict[str, Any]],
) -> dict[str, Any]:
    """Design once per condition and length, then compare what came back.

    `run_design` receives `condition`, `length`, `sizes`, `coverage_targets`,
    `constraints` and `provider`, and returns the `plan_pool` result dictionary.
    Injecting it keeps this module free of the optimizer construction, so a test
    can exercise the sweep without building a position cache.

    A length nothing was counted at raises rather than contributing an empty
    design. An empty frontier that looks like a result is worse than a refusal:
    it reads as "no panel is possible at this length" when the truth is "nobody
    has counted k-mers at this length yet".
    """
    designs = []
    for condition in conditions:
        fingerprint = condition.fingerprint()
        for length in lengths:
            provider = CandidateProvider(inventory, fingerprint, [length])
            # A length nothing was COUNTED at is a missing input, and reporting
            # an empty frontier for it would read as "no panel is possible" when
            # the truth is "nobody has counted k-mers at this length yet".
            if not inventory.has_length(length):
                raise ValueError(
                    f"No candidates of length {length} were ever counted, so this "
                    f"length cannot be compared. Run 'neoswga count-kmers -j "
                    f"params.json' with min_k/max_k covering {length}, then "
                    f"'neoswga filter'."
                )

            # Counted, but nothing this reaction admits, IS a result: this
            # chemistry has no eligible candidate at this length. Recorded as a
            # design with none rather than raised.
            available = provider.initial(limit=1)
            if not available:
                designs.append(
                    {
                        "condition": fingerprint,
                        "length": length,
                        "eligible_candidates": 0,
                        "result": None,
                        "note": "no candidate of this length clears this reaction's hard gates",
                    }
                )
                continue

            result = run_design(
                condition=condition,
                length=length,
                sizes=sizes,
                coverage_targets=coverage_targets,
                constraints=constraints,
                provider=provider,
            )
            designs.append(
                {
                    "condition": fingerprint,
                    "length": length,
                    "eligible_candidates": result.get(
                        "eligible_candidates", provider.unexamined(excluded=[])
                    ),
                    "result": result,
                }
            )

    return {
        "designs": designs,
        "lengths": list(lengths),
        "coverage_targets": list(coverage_targets),
        "constraints": constraints,
    }


def load_design_grid(grid: dict[str, Any], baseline: dict[str, Any]):
    """Turn a design grid into lengths and fully resolved reaction conditions.

    A grid entry names what CHANGES, not the whole reaction. Each entry is
    applied on top of the run's own resolved configuration, so a grid varying
    only DMSO keeps the user's buffer, salts and per-oligo concentration.
    Rebuilding a condition from the overrides alone would compare designs
    against library defaults rather than against the reaction being run, and the
    comparison would look like a chemistry result.

    Every resulting condition is constructed through `ReactionConditions`, so an
    out-of-range dose is refused by the existing model rather than by a second
    copy of its rules here.
    """
    from .reaction_conditions import ReactionConditions

    lengths = list(grid.get("lengths") or [])
    if not lengths:
        raise ValueError("A design grid must name at least one oligo length under 'lengths'.")
    overrides = grid.get("conditions")
    if not overrides:
        raise ValueError(
            "A design grid must name at least one entry under 'conditions'. Use "
            "an empty object to include the baseline reaction unchanged."
        )

    conditions = []
    for entry in overrides:
        resolved = dict(baseline)
        resolved.update(entry or {})
        conditions.append(ReactionConditions(**resolved))
    return lengths, conditions


def frontier(designs: Sequence[dict[str, Any]], coverage_targets: Sequence[float]):
    """The smallest qualifying pool per coverage target, across designs.

    "Qualifying" means eligible: a panel failing a constraint is not a coverage
    result however high its coverage reads, which is the distinction the old
    report blurred by ranking on coverage alone.

    Per target, the smallest qualifying panel from each design is collected and
    the dominated ones dropped, so what remains is the set of genuinely
    different trade-offs rather than one winner chosen by an arbitrary weighting
    between size, coverage and background load.

    Designs that never reach a target are named under `unreached` rather than
    omitted silently: "this chemistry did not get there" is a result.
    """
    by_target: dict[float, list[dict[str, Any]]] = {}
    unreached: dict[float, list[str]] = {}

    for target in coverage_targets:
        candidates = []
        missed = []
        for design in designs:
            result = design.get("result") or {}
            qualifying = [
                row
                for row in result.get("rows", [])
                if row.get("eligible")
                and row.get("coverage") is not None
                and row["coverage"] >= target
            ]
            if not qualifying:
                missed.append(design["condition"])
                continue
            best = min(qualifying, key=lambda row: (row["size"], -row["coverage"]))
            candidates.append(
                {
                    "condition": design["condition"],
                    "length": design["length"],
                    "size": best["size"],
                    "coverage": best["coverage"],
                    "background": best.get("background_sites") or 0,
                }
            )
        by_target[target] = nondominated(candidates)
        unreached[target] = missed

    return {"by_target": by_target, "unreached": unreached}
