"""What the surviving candidate pool can and cannot do, at this reaction.

Two measurements that `filter` reports and neither of which removes a
candidate:

- whether the pool can tell a true site from a near-miss, which is Known Issue
  17's subject and is a pool-wide mean;
- when the pool mixes oligo lengths, how differently those lengths are bound.

The second is not derivable from the first. Occupancy rises steeply with
length at a fixed temperature, so on a mixed pool the pool-wide mean describes
no length in it -- on the bundled plasmid pool it is 2.405 while the
per-length values run 3.010 down to 1.015.

Extracted from `core/pipeline.py` on 2026-09-21, when adding the second
measurement pushed both that module and its `step2` past their size budgets.
Reporting what a pool is, is its own concern.
"""

from __future__ import annotations

import logging

from neoswga.core.occupancy import log_discrimination_profile

logger = logging.getLogger(__name__)

__all__ = ["report_pool_diagnostics"]


def report_pool_diagnostics(pool, parameter) -> None:
    """Describe `pool` under the reaction `parameter` resolves.

    Best-effort throughout: a diagnostic must never fail the step that
    produced the pool it describes.
    """
    log_discrimination_profile(pool)
    try:
        from neoswga.core.length_occupancy import log_occupancy_spread
        from neoswga.core.reaction_conditions import build_reaction_conditions

        log_occupancy_spread(pool, build_reaction_conditions(parameter), label="candidate pool")
    except Exception as exc:  # pragma: no cover - diagnostic only
        logger.debug(f"Could not report occupancy by length: {exc}")
