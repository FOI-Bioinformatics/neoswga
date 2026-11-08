"""One normalised score for a panel, whatever optimizer produced it.

Extracted from `base_optimizer` on 2026-10-02, for the reason that module's
size ratchet exists: adding the per-host background figures to
`PrimerSetMetrics` would have pushed the file past its 1600-line ceiling, and
the ratchet's own instruction is to extract rather than raise. Nothing moved
but the text. `PrimerSetMetrics.normalized_score` and
`PrimerSetMetrics.APPLICATION_WEIGHTS` are the same objects reached by the same
names, and `base_optimizer` re-exports `APPLICATION_WEIGHTS` so an importer
reaching for it there still finds it.

Why a normalised score exists at all: raw `score` is whatever the winning
optimizer's own objective returned, and those are not comparable across
methods -- `hybrid` reports a predicted amplification, `dominating-set` a cover
size. `normalized_score` is the one figure the ensemble ranks on, so it has to
be computed the same way for every method.
"""

from __future__ import annotations

import math
from typing import Any

from .selectivity import SELECTIVITY_REFERENCE

__all__ = ["APPLICATION_WEIGHTS", "normalized_panel_score"]

# Application-profile weights for POST-HOC scoring (Phase 14C).
# Different use cases weigh the component metrics differently:
# clinical applications favour selectivity and dimer safety over raw
# coverage, while metagenomics favours coverage above all.
APPLICATION_WEIGHTS = {
    "balanced": {
        "coverage_w": 0.35,
        "selectivity_w": 0.30,
        "dimer_w": 0.15,
        "evenness_w": 0.10,
        "tm_w": 0.10,
    },
    "discovery": {
        "coverage_w": 0.50,
        "selectivity_w": 0.15,
        "dimer_w": 0.10,
        "evenness_w": 0.15,
        "tm_w": 0.10,
    },
    "clinical": {
        "coverage_w": 0.20,
        "selectivity_w": 0.45,
        "dimer_w": 0.20,
        "evenness_w": 0.10,
        "tm_w": 0.05,
    },
    "enrichment": {
        "coverage_w": 0.35,
        "selectivity_w": 0.30,
        "dimer_w": 0.15,
        "evenness_w": 0.10,
        "tm_w": 0.10,
    },
    "metagenomics": {
        "coverage_w": 0.60,
        "selectivity_w": 0.10,
        "dimer_w": 0.10,
        "evenness_w": 0.15,
        "tm_w": 0.05,
    },
}


def normalized_panel_score(
    metrics: Any,
    coverage_w: float = 0.35,
    selectivity_w: float = 0.30,
    dimer_w: float = 0.15,
    evenness_w: float = 0.10,
    tm_w: float = 0.10,
    application: str | None = None,
) -> float:
    """A [0,1] composite score for one panel, independent of its optimizer.

    All component metrics are clamped to [0,1] before weighting.

    `application` (one of 'balanced', 'discovery', 'clinical', 'enrichment',
    'metagenomics') overrides the individual weights with that preset; see
    `APPLICATION_WEIGHTS`.

    `selectivity_ratio` is the term read here, and it is the POOLED one under
    every `background_aggregate`. That is deliberate: the ratio is not
    comparable across hosts of different size (Known Issue 6), so a per-host
    reduction of it would not mean anything, and reducing the score over hosts
    would also change which method an ensemble run picks -- a decision no
    measurement supports. The per-host figures bind limits, not the score.
    """
    if application is not None:
        weights = APPLICATION_WEIGHTS.get(application.lower())
        if weights is not None:
            coverage_w = weights["coverage_w"]
            selectivity_w = weights["selectivity_w"]
            dimer_w = weights["dimer_w"]
            evenness_w = weights["evenness_w"]
            tm_w = weights["tm_w"]

    # Coverage: already [0,1].
    #
    # Prefer the occupancy-weighted figure when conditions were available to
    # compute it. Selecting on the raw count rewards a set for reaching sites it
    # does not occupy at the reaction temperature, and the bias is not uniform
    # across candidates -- on Prevotella at equiphi29 42 C a 10-mer set reports
    # 0.520 raw against 0.245 effective while a 12-mer set reports 0.403 against
    # 0.387, so the two figures rank them in opposite orders. Optimizing the raw
    # one also hides the coverage that additives cost, which is the lever this
    # tool exists to use.
    #
    # Falls back to the raw figure when no conditions are attached, matching how
    # `selectivity_mode` degrades: None means "not computed"; a computed zero
    # remains zero.
    cov = min(metrics.fg_coverage, 1.0)
    if metrics.effective_fg_coverage is not None:
        cov = min(metrics.effective_fg_coverage, 1.0)

    # Selectivity, on a log scale.
    #
    # This was `min(ratio / 100, 1.0)`, calibrated when selectivity was a count
    # of exact matches -- ratios routinely exceeded 100 and saturated the cap, so
    # the term behaved as a near-binary "good enough". Occupancy weighting counts
    # near-match background load too, so ratios now land around 5-15 on the
    # bundled example, where the same linear form contributes 0.05-0.15 of its
    # weight and the term is effectively ignored instead. Saturating and
    # vanishing are the same failure at opposite ends.
    #
    # A ratio spanning orders of magnitude wants a log scale: a design at
    # SELECTIVITY_REFERENCE scores 1.0, and the term stays responsive across
    # everything below it rather than being decided in the first or last few
    # percent of the range.
    sel = 0.0
    if metrics.selectivity_ratio > 0:
        sel = math.log10(1.0 + metrics.selectivity_ratio) / math.log10(1.0 + SELECTIVITY_REFERENCE)
        sel = max(0.0, min(sel, 1.0))

    # Dimer safety: 1 - risk (lower risk is better)
    dimer_safe = 1.0 - min(metrics.dimer_risk_score, 1.0)

    # Evenness: 1 - Gini (lower Gini is more uniform)
    even = 1.0 - min(metrics.gap_gini, 1.0)

    # Tm tightness: narrow Tm spread is better
    tm_spread = abs(metrics.tm_range[1] - metrics.tm_range[0]) if metrics.tm_range else 0.0
    tm_tight = max(0.0, 1.0 - tm_spread / 20.0)  # 0C spread -> 1.0, 20C+ -> 0.0

    return (
        coverage_w * cov
        + selectivity_w * sel
        + dimer_w * dimer_safe
        + evenness_w * even
        + tm_w * tm_tight
    )
