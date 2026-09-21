"""Portable pool-size report and oligo exports."""

import csv
import html
import json
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

_PLAN_KEYS = {
    "rows",
    "recommendations",
    "primer_length",
    "coverage_metric",
    "extension_reach",
    "background_assessed",
    "min_selectivity_density",
    "max_background_sites",
    "interpretation",
}

# Only what the page and the CSV actually read. Everything else a row may carry
# is optional and fetched with `.get`.
_ROW_KEYS = {"requested_size", "size", "primers", "eligible"}


def _text(value):
    """Escape a value taken from the saved plan for inclusion in the page.

    `report-pool` renders a JSON file this tool did not necessarily write, so
    every interpolated value is treated as text rather than markup.
    """
    return html.escape(str(value))


def _check_plan(plan):
    """Reject an unusable plan before anything is written.

    A malformed row used to escape as a bare `KeyError` or `IndexError`, which
    the CLI reports as a missing *parameter* and answers with advice about
    `neoswga init` and params.json. Neither is the file at fault.
    """
    if not isinstance(plan, dict) or _PLAN_KEYS - plan.keys():
        raise ValueError("Expected saved plan-pool results (pool_plan.json)")
    rows, recommendations = plan["rows"], plan["recommendations"]
    if not isinstance(rows, list) or not isinstance(recommendations, list):
        raise ValueError("Expected saved plan-pool results (pool_plan.json): rows must be a list")
    for row in rows:
        if not isinstance(row, dict) or _ROW_KEYS - row.keys():
            raise ValueError(
                "Expected saved plan-pool results (pool_plan.json): a row is missing "
                + ", ".join(sorted(_ROW_KEYS - (row.keys() if isinstance(row, dict) else set())))
            )
        if not isinstance(row["primers"], list):
            raise ValueError(
                "Expected saved plan-pool results (pool_plan.json): row primers must be a list"
            )
        # The exported FASTA is what reaches a synthesis vendor. `plan_pool`
        # enforces this on the candidates it is given, but a saved plan is
        # rendered without passing through it, and a primer holding a newline
        # split into a second, invented record.
        for primer in row["primers"]:
            if not isinstance(primer, str) or not primer or set(primer) - set("ACGT"):
                raise ValueError(
                    "Expected saved plan-pool results (pool_plan.json): oligos must be "
                    f"non-empty and contain only A, C, G and T, but a row carries {primer!r}"
                )
    for rec in recommendations:
        if not isinstance(rec, dict) or {"target_coverage", "row_index"} - rec.keys():
            raise ValueError(
                "Expected saved plan-pool results (pool_plan.json): a recommendation is incomplete"
            )
        index = rec["row_index"]
        if index is not None and not (isinstance(index, int) and 0 <= index < len(rows)):
            raise ValueError(
                "Expected saved plan-pool results (pool_plan.json): a recommendation points "
                f"at row {index!r}, outside the {len(rows)} rows recorded"
            )


def _concentration_rows(params, plan):
    """Per-oligo and total nominal concentration, both stated.

    They are different experiments. The model assumes a fixed concentration PER
    OLIGO, so a larger panel is a larger total; a fixed total shared across a
    growing panel would lower every oligo's concentration and shift every
    melting temperature, and is not modelled. Reporting only one of the two
    leaves which experiment was assumed to be inferred.
    """
    per_oligo = params.get("primer_conc")
    if per_oligo is None:
        return {}
    sizes = [r.get("size") for r in plan.get("rows", []) if r.get("size")]
    largest = max(sizes) if sizes else None
    rows = {"Concentration per oligo (M)": f"{per_oligo:g}"}
    if largest:
        rows["Total nominal concentration at the largest panel (M)"] = (
            f"{per_oligo * largest:g} ({largest} oligos x {per_oligo:g})"
        )
    return rows


def _accounting_rows(plan):
    """How many candidates existed, and how many the search actually looked at."""
    counts = plan.get("design_counts") or {}
    rows = {}
    if counts:
        rows["Candidates counted / hard-QC passed / examined"] = (
            f"{counts.get('counted', 0):,} / {counts.get('hard_qc_passed', 0):,} / "
            f"{counts.get('examined', 0):,}"
        )
    if plan.get("stop_reason"):
        rows["Search ended because"] = plan["stop_reason"]
    return rows


def _write_figure(plan, rows, output):
    """Draw the coverage and specificity scatter for a saved plan.

    Split out of `write_pool_plan` so that function stays inside the project's
    function-length budget.

    matplotlib is an optional dependency, in the `viz` extra rather than the
    base install, and this module imported it unconditionally. A user who had
    not installed that extra got a `ModuleNotFoundError` from `plan-pool` and
    `report-pool` instead of a report. The rest of the report -- the JSON, the
    CSV, the oligo FASTAs and every number in the page -- needs no plotting, so
    its absence costs one image and nothing else.

    Returns True when the figure was written, so the page can omit the `<img>`
    rather than pointing at a file that is not there.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        from matplotlib import pyplot as plt
    except ImportError:
        logger.info(
            "matplotlib is not installed, so %s holds no figure. "
            "Install the 'viz' extra for it; every other output is unaffected.",
            output / "pool_sizes.png",
        )
        return False

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    fig.suptitle(plan.get("title", "Oligo pool design"))
    for good, color, label in [
        (
            True,
            "#18764a",
            (
                "Meets specificity and dimer limits"
                if plan["background_assessed"]
                else "Dimer limits pass; background unassessed"
            ),
        ),
        (False, "#a84232", "Fails at least one limit"),
    ]:
        subset = [r for r in rows if r["eligible"] == good]
        axes[0].scatter(
            [r["size"] for r in subset],
            [100 * r["coverage"] for r in subset],
            color=color,
            label=label,
        )
    for rec in plan["recommendations"]:
        axes[0].axhline(100 * rec["target_coverage"], color="#777", linestyle=":", linewidth=1)
    axes[0].set(
        xlabel="Oligos in delivered pool",
        ylabel="Estimated target coverage (%)",
        ylim=(0, 102),
        title=f"{plan['primer_length']}-mers; {plan['coverage_metric']} coverage",
    )
    axes[0].legend(fontsize=7)
    known = [r for r in rows if r.get("selectivity_density") is not None]
    if known:
        axes[1].scatter(
            [r["size"] for r in known], [r["selectivity_density"] for r in known], color="#315b8a"
        )
        axes[1].set_yscale("symlog", linthresh=1)
        if plan["min_selectivity_density"] is not None:
            axes[1].axhline(plan["min_selectivity_density"], color="#777", linestyle=":")
    else:
        axes[1].text(0.5, 0.5, "Background not assessed", ha="center", transform=axes[1].transAxes)
    axes[1].set(
        xlabel="Oligos in delivered pool",
        ylabel="Target/background site-density ratio",
        title="Estimated specificity",
    )
    fig.savefig(output / "pool_sizes.png", dpi=160)
    plt.close(fig)
    return True


def _write_sizes_csv(plan, output):
    """One row per evaluated panel size, for a reader who wants the numbers.

    `extrasaction="ignore"` means a plan carrying keys this list does not name
    still writes, which is what lets an older saved plan render.
    """
    fields = [
        "requested_size",
        "size",
        "coverage",
        "raw_coverage",
        "effective_coverage",
        "selectivity_density",
        "background_sites",
        "eligible",
        "status",
        "repaired_by",
        "primers",
    ]
    with (output / "pool_sizes.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in plan["rows"]:
            # A panel the planner had to repair is not the panel the optimizer
            # first returned, and a reader comparing designs should be able to
            # see which rows needed a second attempt.
            repair = row.get("repair") or {}
            writer.writerow(
                {
                    **row,
                    "repaired_by": repair.get("method") or "",
                    "primers": ";".join(row["primers"]),
                }
            )


def write_pool_plan(plan, output):
    """Render a saved plan to a fresh directory, without accessing genome inputs.

    This function is usable from design scripts as well as the CLI. Existing
    reports are preserved to avoid leaving stale oligo exports from another run.
    """
    output = Path(output)
    if output.exists() and not output.is_dir():
        raise ValueError(f"{output} is not a directory; pass a new or empty directory")
    if output.exists() and any(output.iterdir()):
        raise ValueError(
            "Output directory is not empty; use a new directory to preserve earlier designs"
        )
    _check_plan(plan)
    serialized = json.dumps(plan, indent=2, allow_nan=False)
    output.mkdir(parents=True, exist_ok=True)
    (output / "pool_plan.json").write_text(serialized)
    _write_sizes_csv(plan, output)
    targets = []
    for rec in plan["recommendations"]:
        target = rec["target_coverage"]
        if rec["row_index"] is None:
            targets.append(f"<li>{target:.0%} coverage: no qualifying panel found.</li>")
            continue
        row = plan["rows"][rec["row_index"]]
        name = f"target_{target * 100:g}pct_oligos.fasta"
        text = "".join(f">oligo_{i:03d}\n{p}\n" for i, p in enumerate(row["primers"], 1))
        (output / name).write_text(text)
        targets.append(
            f"<li>{target:.0%} coverage: smallest found = {_text(row['size'])} oligos "
            f'(<a href="{name}">FASTA</a>).</li>'
        )
    rows = [r for r in plan["rows"] if r.get("coverage") is not None]
    figure_written = _write_figure(plan, rows, output)
    table = []
    for r in plan["rows"]:
        coverage = f"{r['coverage']:.1%}" if r.get("coverage") is not None else "unavailable"
        density = (
            f"{r['selectivity_density']:.2f}"
            if r.get("selectivity_density") is not None
            else "not assessed"
        )
        # `r["status"]` was the default argument to `.get`, so it was read even
        # when `failed_constraints` was present and a row carrying one but not
        # the other raised.
        failed = r.get("failed_constraints") or ([r["status"]] if "status" in r else [])
        reason = "passes" if r["eligible"] else "; ".join(str(x) for x in failed) or "unknown"
        method = (r.get("repair") or {}).get("method")
        if method:
            reason = f"{reason} (repaired by {method})"
        sites = r["background_sites"] if r.get("background_sites") is not None else "not assessed"
        # Both coverage figures, side by side. One number with a metric name
        # beside it invited the reader to treat a modelling choice as a result;
        # the gap between the two IS the contribution of the temperature and
        # additive model, and it is only visible when both are shown.
        geometric = (
            f"{r['raw_coverage']:.1%}" if r.get("raw_coverage") is not None else "unavailable"
        )
        weighted = (
            f"{r['effective_coverage']:.1%}"
            if r.get("effective_coverage") is not None
            else "unavailable"
        )
        stage_details = ""
        history = r.get("stage_history") or []
        if history:
            items = []
            for stage in history:
                coverage = stage.get("coverage")
                coverage_text = f"{coverage:.1%}" if coverage is not None else "unavailable"
                text = (
                    f"{stage['stage']}: {stage.get('before_size', 0)} to {stage.get('after_size', 0)} oligos; "
                    f"{stage.get('coverage_metric') or 'coverage'} {coverage_text}"
                )
                if stage.get("stop_reason"):
                    text += f"; {stage['stop_reason']}"
                if stage.get("evaluations") is not None:
                    text += f"; {stage['evaluations']} evaluations"
                items.append(f"<li>{_text(text)}</li>")
            stage_details = (
                "<details><summary>Search stages</summary><ul>" + "".join(items) + "</ul></details>"
            )
        table.append(
            f"<tr><td>{_text(r['requested_size'])}</td><td>{_text(r['size'])}</td>"
            f"<td>{_text(geometric)}</td><td>{_text(weighted)}</td>"
            f"<td>{_text(density)}</td><td>{_text(sites)}</td>"
            f"<td>{_text(reason)}{stage_details}</td></tr>"
        )
    density_limit = _text(
        plan["min_selectivity_density"]
        if plan["min_selectivity_density"] is not None
        else "not set"
    )
    background_limit = _text(
        plan["max_background_sites"] if plan["max_background_sites"] is not None else "not set"
    )
    title = _text(plan.get("title", "Oligo pool design"))
    params = plan.get("design_parameters", {})
    inputs = plan.get("inputs", {})

    def reference_names(genome_key, prefix_key):
        paths = params.get(genome_key) or inputs.get(prefix_key) or []
        return ", ".join(Path(p).name for p in paths) or "not recorded"

    context = {
        "Target reference(s)": reference_names("fg_genomes", "foreground_prefixes"),
        "Background reference(s)": (
            reference_names("bg_genomes", "background_prefixes")
            if plan["background_assessed"]
            else "not assessed"
        ),
        "Candidate oligos": plan.get("candidate_count", "not recorded"),
        "Optimizer": plan.get("optimizer", "not recorded"),
        "Polymerase model": params.get("polymerase", "not recorded"),
        "Reaction temperature (C)": params.get("reaction_temp", "not recorded"),
        "Pairwise dimer limit (bp)": plan.get("max_dimer_bp", "not recorded"),
        "Self-dimer limit (bp)": plan.get("max_self_dimer_bp", "not recorded"),
        # The reaction identity, so two reports can be told apart by more than
        # their filenames, and so a coefficient revision is visible.
        "Reaction fingerprint": plan.get("condition_fingerprint", "not recorded"),
        **_concentration_rows(params, plan),
        **_accounting_rows(plan),
    }
    context_html = "".join(
        f"<tr><th>{_text(key)}</th><td>{_text(value)}</td></tr>" for key, value in context.items()
    )
    eligible = [r for r in rows if r["eligible"]]
    best = max(eligible, key=lambda r: (r["coverage"], -r["size"])) if eligible else None
    best_text = (
        f"Highest qualifying estimate among evaluated panels: {best['coverage']:.1%} "
        f"with {_text(best['size'])} oligos."
        if best
        else "No evaluated panel passed all specified constraints."
    )
    figure_html = (
        '<img src="pool_sizes.png" alt="Estimated coverage and specificity ' 'versus oligo count">'
        if figure_written
        else "<p><em>No figure: matplotlib is not installed. Install the "
        "<code>viz</code> extra to include one; every number below is "
        "unaffected.</em></p>"
    )
    sweep = plan.get("reach_sensitivity") or []
    if sweep:
        sweep_cells = []
        for row in sweep:
            reach_text = _text(format(row["reach"], ","))
            coverage_text = _text(format(row["coverage"], ".1%"))
            sweep_cells.append(f"<tr><td>{reach_text}</td><td>{coverage_text}</td></tr>")
        sweep_rows = "".join(sweep_cells)
        sensitivity_html = (
            "<h2>How much of that is the window radius?</h2>"
            "<p>The same delivered panel, scored at other extension reaches. The "
            "reach is a design-density convention taken from published successful "
            "sets, not a measured extension distribution, and coverage is close to "
            "linear in it. Read the recommendation together with this table.</p>"
            "<table><tr><th>Extension reach (bp)</th><th>Estimated coverage</th></tr>"
            f"{sweep_rows}</table>"
        )
    else:
        sensitivity_html = ""

    # What the reported number does and does not account for (audit F2, F7).
    scope_html = (
        "<details><summary>What this coverage figure includes</summary>"
        "<p>It is the union of symmetric windows around exact binding sites, "
        "weighted by an equilibrium occupancy derived from melting temperature "
        "and duplex enthalpy. It is a <strong>proxy</strong>, not a predicted "
        "recovery.</p>"
        "<p>It does not resolve strand-directed extension, competing template "
        "molecules, repeated priming, enzyme activity, reaction duration, "
        "depletion, sequencing depth, or breadth at any depth threshold. An "
        "equilibrium binding fraction is not a probability of amplification or "
        "of read recovery.</p>"
        "<p>Additives act on it only through melting temperature. Glycerol, BSA, "
        "PEG, SSB and DTT have effects in the optional mechanistic model and no "
        "term here, and inhibition of the polymerase by a Tm-active additive is "
        "not represented. A recipe scored well by this metric has been assessed "
        "for binding discrimination, not for enzyme activity or genome recovery."
        "</p></details>"
    )

    page = f"""<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>{title}</title><style>body{{font:16px system-ui;max-width:1100px;margin:40px auto;padding:0 24px;color:#233}}
table{{border-collapse:collapse;width:100%}}td,th{{padding:9px;border-bottom:1px solid #ddd;text-align:left}}img{{width:100%}}a{{color:#176750}}</style>
<h1>{title}: {_text(plan["primer_length"])}-mers</h1>
<p>{_text(plan["interpretation"])}</p>
<p>Coverage metric: {_text(plan["coverage_metric"])}; extension reach: {plan["extension_reach"]:,} bp.
Minimum site-density ratio: {density_limit}; maximum exact background sites: {background_limit}. A density ratio is not a prediction of fold enrichment.</p>
<details><summary>Design inputs and settings</summary><table>{context_html}</table></details>
<h2>Smallest qualifying pools found</h2><ul>{"".join(targets)}</ul>
<p>{best_text}</p>
{figure_html}
{sensitivity_html}
{scope_html}
<h2>Evaluated panels</h2>
<p>Two coverage figures, because they answer different questions. <strong>Geometric</strong> is the union of extension windows around exact binding sites: how much of the target a panel could reach. <strong>Occupancy-weighted</strong> weights each window by how much of the time its site is bound at the reaction temperature: how much it plausibly reaches under this chemistry. A higher modelled density is <strong>not demonstrated enrichment</strong>; every figure here is computed, and <a href="https://github.com/FOI-Bioinformatics/neoswga/blob/main/docs/validation/evidence_matrix_2026-09-15.md">the evidence matrix</a> records which have no compatible observation at all.</p>
<table><tr><th>Requested</th><th>Delivered</th><th>Geometric coverage</th><th>Occupancy-weighted coverage</th><th>Density ratio</th><th>Background sites</th><th>Constraints</th></tr>{"".join(table)}</table>
<p><a href="pool_sizes.csv">All sizes (CSV)</a> · <a href="pool_plan.json">Full results and oligo sequences (JSON)</a></p></html>"""
    (output / "pool_plan.html").write_text(page)
    return output / "pool_plan.html"
