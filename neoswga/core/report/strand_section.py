"""The strand figures, all of them, for every genome that was measured.

The quality report rendered two numbers: the alternation score and the
site-count balance of the FIRST target genome. The optimizer writes ten on a
one-host run -- five per genome, in `metrics.strand_stats` of
`step4_improved_df_summary.json` -- so the host's strand structure and three of
the five figures per genome reached disk and no reader. Known Issue 8's shape,
in the report rather than in the search.

Two readings this section has to keep apart:

**A starved strand on the TARGET is bad and on the HOST is good.** Exponential
amplification needs two sites in convergent orientation within the
polymerase's reach; one site primes linearly at best. A target genome whose
reverse-strand sites are missing amplifies linearly, and a host in the same
state amplifies less, which is the outcome a design wants. The role comes from
params.json rather than from a guess, and a prefix params.json does not place
is labelled `unknown`.

**An unmeasurable alternation score is not a zero.** Alternation is the
fraction of ADJACENT site pairs on opposite strands, so a panel with fewer than
two sites on a genome scores 0/0 and the optimizer writes `null`. The report
read it through a float conversion with a literal 0.0 default and rendered
"0.00", which is the FAVOURABLE-direction confusion of Known Issues 5, 6 and
13: "could not be measured" shown as "measured: no alternation at all".

**Nothing here is a threshold.** Known Issue 18 keeps the strand scalars
unconstrainable, and no cut on them has a reference. The two things flagged are
a fact (one strand carries no sites at all) and a mechanism (the widest stretch
with no convergent pair, against the polymerase's reach, which is the distance
the enzyme can actually span). Neither ranks a panel.

This is the gap `docs/validation/2026-09-22-swga-comparative-review.md` records
on this side: COATswga ranks on the MEAN of its two strand coverages, so
0.99/0.30 scores 0.645 and beats a balanced 0.60/0.60, and nothing in it flags
the imbalance. NeoSWGA had the data and did not report it either.
"""

from __future__ import annotations

from html import escape as html_escape
from typing import Any

#: Column label and figure key, in the order the table renders them.
_COLUMNS: tuple[tuple[str, str], ...] = (
    ("Site balance", "strand_coverage_ratio"),
    ("Alternation", "strand_alternation_score"),
    ("Widest convergent gap (bp)", "strand_alternation_gap_max"),
    ("Mean convergent gap (bp)", "strand_alternation_gap_mean"),
    ("Longest same-strand run", "longest_same_strand_run"),
    # swga 2.0's `off_gap_gini` on a host, `on_gap_gini` on a target.
    ("Gap evenness per strand", "gap_gini_per_strand"),
)

_UNMEASURED = "<em>not measurable</em>"


def _figure(key: str, value: Any) -> str:
    """One cell. `None` says so; it never becomes a number."""
    if value is None:
        return _UNMEASURED
    try:
        number = float(value)
    except (TypeError, ValueError):
        return html_escape(str(value))
    if key == "longest_same_strand_run":
        return f"{int(number):,}"
    if key.endswith("_gap_max") or key.endswith("_gap_mean"):
        return f"{number:,.0f}"
    return f"{number:.2f}"


def _headline_rows(uniformity: Any) -> str:
    """The two scalars the section has always carried, for the first target."""
    rows = ""
    alt = uniformity.strand_alternation_score
    ratio = uniformity.strand_coverage_ratio
    if alt is not None or ratio is not None:
        rows += (
            f"<tr><td>Strand alternation score</td><td>{_figure('x', alt)}</td>"
            "<td>fraction of adjacent sites that alternate strands "
            "(higher = better exponential amplification); 0/0 and so not "
            "measurable below two sites</td></tr>"
            f"<tr><td>Strand coverage ratio</td><td>{_figure('x', ratio)}</td>"
            "<td>min(forward, reverse) / max(forward, reverse) over binding "
            "site COUNTS, not per-strand coverage (1.0 = balanced)</td></tr>"
        )
    return rows


def _genome_rows(stats: dict[str, dict[str, Any]], roles: dict[str, str]) -> str:
    """One row per genome the optimizer could measure.

    A prefix absent from `strand_stats` gets no row rather than a row of
    zeros: nobody measured it, which is not the same as measuring zero.
    """
    rows = ""
    for prefix in sorted(stats):
        figures = stats[prefix] or {}
        cells = "".join(f"<td>{_figure(key, figures.get(key))}</td>" for _label, key in _COLUMNS)
        role = roles.get(prefix, "unknown")
        rows += f"<tr><td>{html_escape(prefix)}</td><td>{html_escape(role)}</td>{cells}</tr>"
    return rows


def _starved_strand(prefix: str, role: str, figures: dict[str, Any]) -> str | None:
    """The one fact worth flagging: a genome with sites on one strand only.

    `strand_coverage_ratio` is min/max over the site counts, so 0.0 with any
    site at all means the minority strand is empty. That is a measurement and
    not a threshold, which is why it is reportable where a balance CUT is not.
    """
    ratio = figures.get("strand_coverage_ratio")
    if ratio is None or float(ratio) > 0.0:
        return None
    if role == "host":
        return (
            f"<li><strong>{html_escape(prefix)}</strong> (host) has binding sites on one "
            "strand only. On a host that is the favourable direction: with no convergent "
            "pair the host primes linearly at best.</li>"
        )
    return (
        f"<li><strong>{html_escape(prefix)}</strong> ({html_escape(role)}) has binding sites "
        "on one strand only, so it holds no convergent pair anywhere and amplifies linearly "
        "at best. This is the imbalance a mean of the two strand coverages hides.</li>"
    )


def _gap_against_reach(prefix: str, role: str, figures: dict[str, Any], reach: int) -> str | None:
    """The widest stretch with no convergent pair, in units of the reach.

    The reach is the distance the polymerase can span, so a gap wider than it
    is a mechanism statement rather than a calibrated limit: no pair of
    opposite-strand sites within that stretch can face each other across it.
    No threshold on the ratio is implied, and none has a reference.
    """
    gap = figures.get("strand_alternation_gap_max")
    if gap is None or not reach:
        return None
    multiple = float(gap) / float(reach)
    if multiple <= 1.0:
        return None
    if role == "host":
        return (
            f"<li><strong>{html_escape(prefix)}</strong> (host): widest stretch with no "
            f"convergent pair is {float(gap):,.0f} bp, {multiple:.1f}x the {reach:,} bp reach. "
            "Larger is better here. No reference.</li>"
        )
    return (
        f"<li><strong>{html_escape(prefix)}</strong> ({html_escape(role)}): widest stretch with "
        f"no convergent pair is {float(gap):,.0f} bp, {multiple:.1f}x the {reach:,} bp reach, so "
        "no two sites face each other across it. No reference.</li>"
    )


def _notes(stats: dict[str, dict[str, Any]], roles: dict[str, str], reach: int) -> str:
    """At most one note per genome.

    The starved-strand note wins when both would fire, because the gap figures
    hold the GENOME LENGTH when no opposite-strand pair exists at all. Printing
    that as "423x the reach" would dress a conservative encoding up as a
    measured gap, which is the confusion the whole section exists to avoid.
    """
    items = ""
    for prefix in sorted(stats):
        figures = stats[prefix] or {}
        role = roles.get(prefix, "unknown")
        note = _starved_strand(prefix, role, figures) or _gap_against_reach(
            prefix, role, figures, reach
        )
        if note:
            items += note
    if not items:
        return ""
    return f"<ul>{items}</ul>"


def render_strand_balance(uniformity: Any, reach: int | None = None) -> str:
    """The strand section, or `""` when the optimizer measured nothing."""
    if uniformity is None or not getattr(uniformity, "from_optimizer", False):
        return ""
    stats = {
        str(prefix): dict(figures)
        for prefix, figures in (getattr(uniformity, "strand_stats", None) or {}).items()
        if isinstance(figures, dict)
    }
    roles = dict(getattr(uniformity, "strand_roles", None) or {})
    headline = _headline_rows(uniformity)
    if not headline and not stats:
        return ""

    html = (
        '<div class="section"><h2>Strand Balance <span class="badge measured">MEASURED</span></h2>'
    )
    if headline:
        html += f"<table><tr><th>Metric</th><th>Value</th><th>Meaning</th></tr>{headline}</table>"
    if stats:
        headers = "".join(f"<th>{html_escape(label)}</th>" for label, _key in _COLUMNS)
        html += (
            "<p>Per genome. A target whose strands are starved on one side amplifies "
            "linearly there; a host in the same state amplifies less, which is the "
            "favourable direction. None of these figures is constrainable "
            "(Known Issue 18) and none has a published threshold.</p>"
            f"<table><tr><th>Genome</th><th>Role</th>{headers}</tr>"
            f"{_genome_rows(stats, roles)}</table>"
        )
        html += _notes(stats, roles, int(reach or 0))
    return html + "</div>"
