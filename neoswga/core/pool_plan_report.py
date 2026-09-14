"""Portable pool-size report and oligo exports."""

import csv
import html
import json
from pathlib import Path

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
        "primers",
    ]
    with (output / "pool_sizes.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in plan["rows"]:
            writer.writerow({**row, "primers": ";".join(row["primers"])})
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
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    fig.suptitle(plan.get("title", "Oligo pool design"))
    rows = [r for r in plan["rows"] if r.get("coverage") is not None]
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
        sites = r["background_sites"] if r.get("background_sites") is not None else "not assessed"
        table.append(
            f"<tr><td>{_text(r['requested_size'])}</td><td>{_text(r['size'])}</td>"
            f"<td>{_text(coverage)}</td><td>{_text(density)}</td><td>{_text(sites)}</td>"
            f"<td>{_text(reason)}</td></tr>"
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
<img src="pool_sizes.png" alt="Estimated coverage and specificity versus oligo count">
<h2>Evaluated panels</h2><table><tr><th>Requested</th><th>Delivered</th><th>Coverage</th><th>Density ratio</th><th>Background sites</th><th>Constraints</th></tr>{"".join(table)}</table>
<p><a href="pool_sizes.csv">All sizes (CSV)</a> · <a href="pool_plan.json">Full results and oligo sequences (JSON)</a></p></html>"""
    (output / "pool_plan.html").write_text(page)
    return output / "pool_plan.html"
