"""The per-host background table in the technical report.

Its own module because `report/technical_report.py` is at its size budget, and
because the table has one subject: what a panel does to each listed host on
that host's own length, and which reading of "the background" the limits were
applied under.

It is a sibling of `_render_per_target_coverage` in that module and follows it,
including the two-or-more threshold: with one host the pooled figures in the
specificity section already say this, and a one-row table claiming to be a
comparison would be worse than no table.
"""

from __future__ import annotations

from .utils import html_escape

__all__ = ["render_per_host_background"]


def render_per_host_background(specificity) -> str:
    """Render per-host background figures when a run named two or more hosts.

    `selectivity_density` is the column a reader should compare between rows,
    and the note says so. The count ratio is deliberately not here: it carries
    no genome length, so between hosts of different size it moves with the host
    rather than with the panel (Known Issue 6), and a column of it would invite
    exactly the comparison that does not hold.
    """
    if specificity is None:
        return ""
    per_host = getattr(specificity, "per_host", None) or {}
    if len(per_host) < 2:
        return ""
    aggregate = getattr(specificity, "background_aggregate", "sum")
    rows = ""
    for prefix, row in per_host.items():
        reason = row.get("unavailable")
        if reason:
            rows += (
                f"<tr><td>{html_escape(str(prefix))}</td>"
                f"<td>{_format_bp(row.get('length_bp'))}</td>"
                f'<td colspan="4">not measured: {html_escape(str(reason))}</td></tr>'
            )
            continue
        rows += (
            f"<tr><td>{html_escape(str(prefix))}</td>"
            f"<td>{_format_bp(row.get('length_bp'))}</td>"
            f"<td>{_host_cell(row.get('sites'), '{:,.0f}')}</td>"
            f"<td>{_host_cell(row.get('sites_per_mb'), '{:,.1f}')}</td>"
            f"<td>{_host_cell(row.get('coverage'), '{:.1%}')}</td>"
            f"<td>{_host_cell(row.get('selectivity_density'), '{:,.2f}')}</td></tr>"
        )
    applied = (
        "Every limit was applied to the WORST host for that limit."
        if aggregate == "worst-case"
        else "Limits were applied to the pooled background: every host at once, "
        "one summed load over one summed length."
    )
    return (
        '<div class="section"><h2>Per-Host Background '
        '<span class="badge measured">MEASURED</span></h2>'
        "<table><tr><th>Host</th><th>Length</th><th>Sites</th><th>Sites/Mb</th>"
        "<th>Coverage</th><th>Selectivity density</th></tr>"
        f"{rows}</table>"
        f"<p><strong>background_aggregate: {html_escape(str(aggregate))}.</strong> "
        f"{applied} Each host's figures are measured against its OWN length, so "
        "<em>selectivity density</em> is the column to compare between rows; a "
        "count ratio is not comparable across hosts of different size.</p></div>"
    )


def _host_cell(value, template: str) -> str:
    """One per-host number, or an explicit absence.

    `_safe_num` returns 0.0 for None, and zero is the FAVOURABLE value for
    every column here: no sites, no coverage. A host that was not measured must
    not render as the best possible result.
    """
    if value is None:
        return "not measured"
    try:
        return template.format(float(value))
    except (TypeError, ValueError):
        return "not measured"


def _format_bp(value) -> str:
    if value is None:
        return "not recorded"
    try:
        return f"{int(value):,} bp"
    except (TypeError, ValueError):
        return "not recorded"
