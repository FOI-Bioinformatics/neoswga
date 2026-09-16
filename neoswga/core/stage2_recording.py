"""What stage 2 durably records, and what it indexes.

Lifted out of `core/pipeline.py`, which sits against a 1600-line budget. These
two belong together and belong apart from the filter step: one decides which
candidates get a background position index, the other writes the verdict and its
provenance. Both answer "what does this run leave behind", which is a different
question from "what does this run compute".

`pipeline` is imported inside the functions rather than at module scope, because
it imports this module and the cycle would otherwise bite at import time.
"""

from __future__ import annotations

import logging

from neoswga.core import parameter

logger = logging.getLogger(__name__)


def _record_run_inventory(cleared_hard_gates, after_gini, shortlisted, indexed):
    """Durably record what this filter run enumerated and what it merely ranked.

    Keyed on the reaction chemistry and on the admission rules actually in
    force, so a re-run under different thresholds opens its own generation
    instead of leaving the previous run's verdicts selectable.
    """
    from neoswga.core.candidate_inventory import record_stage2_inventory
    from neoswga.core.filter import _get_reaction_conditions
    from neoswga.core.qc_policy import resolved_qc_policy

    record_stage2_inventory(
        parameter.data_dir,
        condition_id=_get_reaction_conditions().fingerprint(),
        cleared_hard_gates=cleared_hard_gates,
        after_gini=after_gini,
        shortlisted=shortlisted,
        indexed=indexed,
        qc_policy=resolved_qc_policy(),
    )


def _index_background_by_retention(
    filtered_rate_df, gini_df, filtered_gini_df, bg_prefixes, bg_genomes
):
    """Index by retention policy, search by shortlist.

    The background scan is one Aho-Corasick pass over the background for every
    primer at once, so widening it costs storage rather than proportional time.
    Measured on the bundled plasmid example: filter 1.4 s to 1.6 s, position
    files 4.1 MB to 8.0 MB, and `step2_df.csv` unchanged at 37 rows -- so the
    optimizer's pool, and therefore its runtime, does not move with this setting.
    """
    from neoswga.core.candidate_inventory import background_scan_pool

    retention = getattr(parameter, "candidate_retention", "all_qc")
    to_index = background_scan_pool(
        filtered_rate_df["primer"], retention, after_gini=gini_df["primer"]
    )
    from neoswga.core.pipeline import _scan_background_positions

    _scan_background_positions(to_index, bg_prefixes, bg_genomes)
    # Name the mode that was NOT used. Naming the active one made the sentence
    # read as a hypothetical about what had just happened: under 'post_gini' it
    # said "'post_gini' would index 20670" beside having indexed exactly that.
    alternative = {
        "all_qc": f"'post_gini' would index the {len(gini_df)} that also cleared "
        "the evenness gate",
        "post_gini": f"'all_qc' would index all {len(filtered_rate_df)}",
    }.get(retention, "no other mode is configured")
    logger.info(
        "Background index: %d of %d hard-QC candidates indexed under "
        "candidate_retention=%r; %s. The %d-candidate max_primer shortlist no "
        "longer bounds the index: a candidate the ranking cut used to have no "
        "background index at all, which scores as perfect specificity rather "
        "than as a missing measurement.",
        len(to_index),
        len(filtered_rate_df),
        retention,
        alternative,
        len(filtered_gini_df),
    )
    return to_index
