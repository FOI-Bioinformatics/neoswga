"""Writing the step 4 results: the CSV, the summary JSON and the audit record.

Split out of `unified_optimizer`, which is about SELECTING a primer set and
sits against the size ratchet in `tests/test_module_size_ratchet.py`. The
serialisation had grown to about a hundred lines and shares nothing with the
optimisation logic beyond the result object.
"""

import json
import logging
import os
from typing import Optional

import pandas as pd

from . import parameter
from .base_optimizer import OptimizationResult, OptimizationStatus
from .dimer import worst_heterodimer
from .step4_audit_trail import write_audit_trail as _write_audit_trail

logger = logging.getLogger(__name__)


def _result_normalized_score(result, application=None):
    """The comparable [0,1] score, under the profile the run actually used.

    `OptimizationResult.normalized_score` is a property and so cannot take the
    profile; it applies the balanced weights. Everywhere the profile matters --
    the CSV column, the summary -- should agree with the ensemble comparison
    table, which is application-weighted.
    """
    if result.status in (OptimizationStatus.ERROR, OptimizationStatus.NO_CONVERGENCE):
        return 0.0
    return float(result.metrics.normalized_score(application=application))


def save_results(
    result: OptimizationResult,
    output_path: str,
    include_all_sets: bool = False,
    application: Optional[str] = None,
    primer_sets=None,
) -> None:
    """
    Save optimization results to CSV.

    Args:
        result: OptimizationResult to save
        output_path: Path for output CSV
        include_all_sets: Whether to include all found sets (not just best)
    """
    if not result.primers:
        logger.warning("No primers to save")
        return

    # Build results DataFrame with set-level aggregate metrics.
    #
    # `set_index` was hardcoded to 0 because a run only ever produced one set.
    # It now numbers the alternatives `max_sets` asked for, best first, which is
    # what the column was there for. Set-level metrics belong to the primary
    # set; the alternatives are offered as candidates, not scored separately.
    metrics = result.metrics
    sets = [tuple(s) for s in (primer_sets or [result.primers]) if s]
    records = []
    for set_index, primers in enumerate(sets):
        for primer in primers:
            records.append(
                {
                    "primer": primer,
                    "set_index": set_index,
                    "score": result.score,
                    # Weighted by the SAME application profile the run was
                    # scored and (for an ensemble) selected under. The property
                    # applies the balanced defaults regardless, so under
                    # `--application clinical` this column disagreed with the
                    # ensemble comparison table for the very same primer set.
                    "normalized_score": _result_normalized_score(result, application),
                    "coverage": metrics.fg_coverage,
                    "bg_coverage": metrics.bg_coverage,
                    "selectivity": metrics.selectivity_ratio,
                    "mean_gap": metrics.mean_gap,
                    "max_gap": metrics.max_gap,
                    "coverage_uniformity": metrics.coverage_uniformity,
                    "gap_gini": metrics.gap_gini,
                    "gap_entropy": metrics.gap_entropy,
                    "dimer_risk_score": metrics.dimer_risk_score,
                    "strand_alternation": metrics.strand_alternation_score,
                    "strand_coverage_ratio": metrics.strand_coverage_ratio,
                    "mean_tm": metrics.mean_tm,
                    "optimizer": result.optimizer_name,
                }
            )

    df = pd.DataFrame(records)
    df.to_csv(output_path, index=False)
    logger.info(f"Results saved to {output_path}")

    # Write summary JSON alongside CSV
    summary_path = os.path.splitext(output_path)[0] + "_summary.json"
    try:
        summary = result.to_dict()
        # Synthesis cost, itemised. Set size is the main lever on coverage, so
        # "how many primers" is really "how many oligos will you buy"; without
        # this the trade-off has a missing axis. Assumes the standard SWGA
        # modification (two 3' PTO bonds against phi29's exonuclease), which is
        # the largest single component for primers this short.
        try:
            from neoswga.core.oligo_cost import cost_per_coverage_point, set_cost

            cost = set_cost(list(result.primers))
            summary["synthesis_cost"] = cost.to_dict()
            summary["synthesis_cost"]["usd_per_coverage_point"] = round(
                cost_per_coverage_point(cost, result.metrics.fg_coverage), 2
            )
        except Exception as exc:  # never lose the result over a cost estimate
            logger.debug(f"Could not estimate synthesis cost: {exc}")

        # The worst pair in the pool being delivered, recorded unconditionally
        # so a reader can see it without re-deriving it. Whether it is a problem
        # depends on the configured max_dimer_bp, which the validator applies;
        # the measurement itself is always worth having.
        try:
            length, pair = worst_heterodimer(list(result.primers))
            summary["worst_heterodimer"] = {
                "length_bp": length,
                "pair": list(pair) if pair else None,
            }
        except Exception as exc:
            logger.debug(f"Could not measure the delivered pool's dimers: {exc}")

        with open(summary_path, "w") as f:
            json.dump(summary, f, indent=2, default=str)
        logger.info(f"Summary saved to {summary_path}")
    except Exception as e:
        logger.warning(f"Could not write summary JSON: {e}")

    # Write audit trail for reproducibility
    _write_audit_trail(output_path, result)
