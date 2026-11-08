"""`filter_stats.json`: the real per-stage funnel of step 2, and what shaped it.

Extracted from `core/pipeline.py` on 2026-10-02 when recording
`background_aggregate` beside `occupancy_ranking_mode` put that module against
its size ceiling. `pipeline` re-exports the function as `_write_filter_stats`,
the name it had there, so nothing that reads it had to change.

Reports showed a fabricated per-stage estimate before this file existed, which
is why it is written at all.
"""

import json
import logging
import os

from neoswga.core import parameter

logger = logging.getLogger(__name__)

__all__ = ["write_filter_stats"]


def write_filter_stats(funnel, final_candidates, background_aggregate=None):
    """Write the funnel, plus the two keys that are not stage counts.

    `occupancy_ranking_mode` names the mismatch model that ordered the
    candidates. The ranking's own output is the `occupancy_ratio` column of
    step2_df.csv, so without it the pool's ORDER came out of an extrapolated
    table with nothing saying so. It is "occupancy" under the shipped uniform
    model, so no existing figure in this file moves.

    `background_aggregate` names the reading of "the background" the frequency
    gate applied, as the gate reported it rather than as a second lookup would
    answer. It is "sum" on every design that does not set the key, so again no
    existing figure moves.

    Both are keys the report reads. `report.metrics.FilteringStats` is built
    from this file, so a key it does not declare is a TypeError in `report`
    rather than an ignored field, and a declared field the loader's whitelist
    omits is read by nothing; that whitelist is now derived from the dataclass
    for exactly this reason.

    Best-effort: a diagnostic must never fail the step that produced it.
    """
    try:
        funnel["final_candidates"] = final_candidates
        from .mismatch_model import resolve_mismatch_model, site_load_mode

        funnel["occupancy_ranking_mode"] = site_load_mode(resolve_mismatch_model())
        if background_aggregate is not None:
            funnel["background_aggregate"] = background_aggregate
        with open(os.path.join(parameter.data_dir, "filter_stats.json"), "w") as handle:
            json.dump(funnel, handle, indent=2)
    except Exception as exc:  # pragma: no cover - defensive
        logger.debug(f"Could not write filter_stats.json: {exc}")
