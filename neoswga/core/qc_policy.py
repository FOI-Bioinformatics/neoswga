"""Which admission rules a run actually enforced.

The candidate inventory keys a verdict on the reaction chemistry and on the QC
policy that produced it. The chemistry half has a fingerprint already; this is
the other half.

It lives apart from `candidate_inventory` because that module is deliberately
free of the `parameter` globals, and apart from `pipeline` because the question
"what did this run enforce" is asked by anything that records or reads a
verdict, not only by the filter step.
"""

from __future__ import annotations

from typing import Any, Dict

# Gates that decide ADMISSION. `max_primer` is deliberately absent: it is a cut
# through a ranking rather than a requirement, so a run that changes it has not
# changed which candidates qualify and must not retire their verdicts.
ADMISSION_THRESHOLDS = (
    "min_fg_freq",
    "max_bg_freq",
    "max_bl_freq",
    "min_tm",
    "max_tm",
    "gc_min",
    "gc_max",
    "max_gini",
    "min_gini_sites",
    "max_homopolymer_run",
    "gc_clamp_window",
    "max_gc_in_clamp",
    "max_self_dimer_bp",
    "excl_threshold",
    "min_k",
    "max_k",
    # The retention mode decides admission under `post_gini`, where the evenness
    # gate stops being a ranking and becomes a requirement. Two modes sharing a
    # digest would let one hand its verdicts to the other.
    "candidate_retention",
)


def resolved_qc_policy(parameter=None) -> Dict[str, Any]:
    """The hard-QC thresholds in force, read after `get_params` has populated them.

    What was ACTUALLY enforced rather than what the file asked for: the
    GC-adaptive strategy and `retune_for_polymerase` both move these at run
    time and never write back to params.json.

    A name that is not set is recorded as `None` rather than filled in with its
    default, because "not configured" and "configured to the default" are
    different statements about what was enforced, and the fingerprint should
    tell them apart.
    """
    if parameter is None:
        from neoswga.core import parameter as parameter_module

        parameter = parameter_module
    return {name: getattr(parameter, name, None) for name in ADMISSION_THRESHOLDS}
