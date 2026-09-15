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

from typing import Any, Callable, Dict, List, Sequence

from .candidate_provider import CandidateProvider
from .pool_objective import PoolConstraints


def nondominated(designs: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
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
    run_design: Callable[..., Dict[str, Any]],
) -> Dict[str, Any]:
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
