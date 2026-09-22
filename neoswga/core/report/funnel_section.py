"""How the candidate pool narrowed, and how much of it the search looked at.

Two renderers, one subject. The funnel shows the filtering stages down to the
shortlist; the reach note answers the question the funnel stops one step short
of -- how many of the candidates behind that shortlist the search then
examined.

Extracted from `technical_report.py` on 2026-09-22. That module was 3 lines
under its ceiling, which is where `base_optimizer.py` was when three pull
requests each measuring exactly its budget merged to 3 over it. A ceiling with
no headroom fails on the merge rather than on the change;
`tests/test_a_ceiling_needs_headroom.py` now says so before it happens.
"""

from __future__ import annotations

from html import escape as html_escape
from typing import Dict, List, Optional

__all__ = ["render_candidate_reach", "render_funnel"]


def render_candidate_reach(reach: Optional[Dict]) -> str:
    """How much of the available pool the search examined.

    Rendered beneath the funnel because it answers the question the funnel
    stops one step short of. The funnel ends at "shortlisted by max_primer";
    this says how many of the candidates behind that shortlist the search
    then looked at.

    Nothing is rendered when the figure is absent. A run over a named
    candidate list has no universe behind it, and a directory written before
    this field existed has no record either -- and 0% would say the search
    reached none of the pool, which is the favourable-default failure this
    project keeps meeting from the other direction.
    """
    if not reach or not reach.get("universe"):
        return ""

    examined, universe = int(reach["examined"]), int(reach["universe"])
    if reach.get("complete"):
        return f"<p><em>The search examined all {universe:,} available candidates." f"</em></p>"

    return (
        f"<p><em>The search examined <strong>{examined:,} of {universe:,}</strong> "
        f"available candidates ({examined / universe:.1%}). It stops at the first "
        f"frontier that satisfies the configured limits; widening it further was "
        f"measured to raise coverage slightly and cost more specificity, so this "
        f"is a deliberate default rather than an incomplete run. See "
        f"<code>docs/validation/looking_further_costs_specificity_2026-09-22.md"
        f"</code>.</em></p>"
    )


def render_funnel(stages: List[tuple]) -> str:
    """Render filtering funnel visualization."""
    if not stages:
        return (
            "<p><em>Filtering funnel not recorded for this run "
            "(filter_stats.json absent). Re-run <code>neoswga filter</code> to "
            "capture per-stage counts.)</em></p>"
        )

    max_count = max(s[1] for s in stages) if stages else 1
    html = ""

    for i, (label, count) in enumerate(stages):
        width_pct = max(5, (count / max_count) * 100)
        pct_of_total = (count / stages[0][1] * 100) if stages[0][1] > 0 else 0
        safe_label = html_escape(str(label))

        html += f"""
        <div class="funnel-stage">
            <div class="funnel-bar" style="width: {width_pct}%">{count:,}</div>
            <span class="funnel-label">{safe_label}</span>
            <span class="funnel-pct">{pct_of_total:.1f}%</span>
        </div>
        """
    return html
