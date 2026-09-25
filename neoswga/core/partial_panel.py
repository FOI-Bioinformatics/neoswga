"""Whether a violating partial panel is worth extending.

Task 6 of the condition-aware pool design plan.

A search over partial panels has to decide whether a panel that already violates
a constraint can be abandoned. Getting that wrong in one direction wastes work;
in the other it silently discards panels that would have qualified, which is the
worse error because nothing in the output records that it happened.

The constraints behave differently under additions, and the difference is not a
detail of implementation:

**Background site count is monotonic.** Every primer added contributes its own
sites, so the total never falls. A partial panel already above the cap cannot be
rescued by anything added later, and pruning it is safe.

**Selectivity density is not.** It is a ratio of target to background site
density, and a primer that binds the target well and the background little
RAISES it. A partial panel below the floor may be above it two primers later.
Pruning on density discards solutions.

**A dimer violation is not a partial solution at all.** No addition removes a
pair that is already in the panel, so such a panel is excluded rather than
pruned as unpromising.

Anything unrecognised is not pruned. A constraint nobody has reasoned about here
might well be non-monotonic, and the conservative direction is to keep looking.
"""

from __future__ import annotations

from collections.abc import Sequence

# A violation that no addition can undo, so the panel is excluded outright.
_UNRECOVERABLE = frozenset({"dimer constraint"})

# A violation whose quantity only grows as primers are added, so a panel that
# has it will still have it however it is extended.
_MONOTONIC = "background sites above maximum"


def can_prune(
    violations: Sequence[str],
    max_background_sites: int | None,
    min_selectivity_density: float | None,
) -> bool:
    """Whether this partial panel can be abandoned without losing a solution.

    `min_selectivity_density` is accepted and deliberately unused in the
    decision: its presence is what makes the density violation non-prunable, and
    naming it in the signature keeps that visible at every call site rather than
    leaving a reader to infer why density is absent from the logic.
    """
    if not violations:
        return False

    if any(v in _UNRECOVERABLE for v in violations):
        return True

    # Only prunable when a cap is actually configured; without one there is no
    # monotonic bound to prune against, whatever the violation string says.
    if max_background_sites is not None and _MONOTONIC in violations:
        return True

    return False
