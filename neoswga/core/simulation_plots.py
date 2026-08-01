#!/usr/bin/env python3
"""
Visualization for SWGA simulation results.

Generates plots and visualizations for simulation analysis including:
- Genome coverage heatmaps
- Primer distribution plots
- Gap visualization
- Enrichment comparison
- Primer contribution charts

Usage:
    from neoswga.core.simulation_plots import generate_plots

    generate_plots(simulation_result, simulator, analysis, output_file)
"""

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

try:
    import matplotlib

    matplotlib.use("Agg")  # Non-interactive backend
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    _MATPLOTLIB_AVAILABLE = True
except ImportError:
    plt = None
    Rectangle = None
    mpatches = None
    _MATPLOTLIB_AVAILABLE = False


def _require_matplotlib():
    if not _MATPLOTLIB_AVAILABLE:
        raise ImportError(
            "matplotlib is required for simulation plotting. "
            "Install with: pip install neoswga[viz]"
        )


def plot_replication_summary(results: dict, output_file: str = "simulation_plots.png") -> None:
    """Plot a summary of the agent-based replication simulation.

    Consumes the dict returned by
    ``replication_simulator.simulate_primer_set`` (keys ``mean_coverage``,
    ``mean_amplification``, ``all_results`` = per-replicate ``SimulationResult``
    objects). Produces a 2x2 PNG: per-replicate coverage, per-replicate
    amplification, the mean coverage profile across the genome, and a metrics
    text panel. Uses only real simulation data (no fabricated values).
    """
    _require_matplotlib()
    all_results = list(results.get("all_results", []) or [])

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("SWGA Replication Simulation Summary", fontsize=14, fontweight="bold")

    # (0,0) per-replicate final coverage
    covs = [getattr(r, "final_coverage_fraction", 0.0) for r in all_results]
    if covs:
        axes[0, 0].bar(range(1, len(covs) + 1), [c * 100 for c in covs], color="#2c7fb8")
    axes[0, 0].set_title("Final coverage per replicate")
    axes[0, 0].set_xlabel("Replicate")
    axes[0, 0].set_ylabel("Coverage (%)")

    # (0,1) per-replicate amplification fold
    amps = [getattr(r, "amplification_fold", 0.0) for r in all_results]
    if amps:
        axes[0, 1].bar(range(1, len(amps) + 1), amps, color="#7fcdbb")
    axes[0, 1].set_title("Amplification fold per replicate")
    axes[0, 1].set_xlabel("Replicate")
    axes[0, 1].set_ylabel("Fold")

    # (1,0) mean coverage profile across the genome (binned)
    arrays = [
        np.asarray(getattr(r, "coverage", None), dtype=float)
        for r in all_results
        if getattr(r, "coverage", None) is not None
    ]
    if arrays:
        mean_cov = np.mean(arrays, axis=0)
        n_bins = min(200, len(mean_cov)) or 1
        edges = np.linspace(0, len(mean_cov), n_bins + 1, dtype=int)
        binned = [
            mean_cov[edges[i] : edges[i + 1]].mean() if edges[i + 1] > edges[i] else 0.0
            for i in range(n_bins)
        ]
        axes[1, 0].fill_between(range(n_bins), [b * 100 for b in binned], color="#41b6c4")
    axes[1, 0].set_title("Mean coverage profile")
    axes[1, 0].set_xlabel("Genome position (binned)")
    axes[1, 0].set_ylabel("Covered (%)")

    # (1,1) metrics text panel
    axes[1, 1].axis("off")
    lines = [
        f"Replicates: {len(all_results)}",
        f"Mean coverage: {results.get('mean_coverage', 0.0):.1%} "
        f"+/- {results.get('std_coverage', 0.0):.1%}",
        f"Mean amplification: {results.get('mean_amplification', 0.0):.1f}x "
        f"+/- {results.get('std_amplification', 0.0):.1f}x",
        f"Mean forks created: {results.get('mean_forks_created', 0.0):.0f}",
        f"Mean fork travel: {results.get('mean_fork_travel', 0.0):,.0f} bp",
        f"Mean displaced strands: {results.get('mean_displaced_strands', 0.0):.0f}",
    ]
    axes[1, 1].text(
        0.02,
        0.95,
        "\n".join(lines),
        va="top",
        ha="left",
        fontsize=12,
        family="monospace",
        transform=axes[1, 1].transAxes,
    )

    plt.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(output_file, dpi=120)
    plt.close(fig)
    logger.info(f"Simulation plots saved to {output_file}")


def generate_plots(result, simulator, analysis=None, output_file="simulation_plots.png"):
    """
    Generate comprehensive visualization of simulation results.

    Args:
        result: SimulationResult from SwgaSimulator
        simulator: SwgaSimulator instance (for access to positions)
        analysis: Optional ComprehensiveAnalysis from SimulationAnalyzer
        output_file: Output file path
    """
    _require_matplotlib()
    logger.info("Generating simulation plots...")

    # Create figure with subplots
    if analysis:
        fig = plt.figure(figsize=(16, 12))
        gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    else:
        fig = plt.figure(figsize=(16, 8))
        gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)

    # Plot 1: Genome coverage heatmap
    ax1 = fig.add_subplot(gs[0, :])
    plot_coverage_heatmap(ax1, simulator.fg_positions, simulator.fg_length)

    # Plot 2: Primer distribution
    ax2 = fig.add_subplot(gs[1, 0])
    plot_primer_distribution(ax2, simulator.fg_positions, simulator.fg_length)

    # Plot 3: Enrichment comparison
    ax3 = fig.add_subplot(gs[1, 1])
    plot_enrichment(ax3, result, simulator)

    # Plot 4: Coverage metrics
    ax4 = fig.add_subplot(gs[1, 2])
    plot_metrics_summary(ax4, result)

    # Additional plots if analysis available
    if analysis:
        # Plot 5: Primer contributions
        ax5 = fig.add_subplot(gs[2, :2])
        plot_primer_contributions(ax5, analysis.primer_contributions[:10])

        # Plot 6: Gap locations
        ax6 = fig.add_subplot(gs[2, 2])
        plot_gaps(ax6, analysis.coverage.gaps)

    plt.suptitle("SWGA Simulation Results", fontsize=16, fontweight="bold")

    # Save figure
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()

    logger.info(f"Plots saved to: {output_file}")


def plot_coverage_heatmap(ax, positions_dict: Dict, genome_length: int):
    """
    Plot genome coverage as heatmap.

    Shows primer binding across the genome.
    """
    # Create bins
    bin_size = 10000
    n_bins = (genome_length + bin_size - 1) // bin_size
    coverage_array = np.zeros(n_bins)

    # Count binding sites per bin
    for primer, pos_data in positions_dict.items():
        for pos in pos_data["+"]:
            bin_idx = int(pos) // bin_size
            if 0 <= bin_idx < n_bins:
                coverage_array[bin_idx] += 1
        for pos in pos_data["-"]:
            bin_idx = int(pos) // bin_size
            if 0 <= bin_idx < n_bins:
                coverage_array[bin_idx] += 1

    # Reshape for heatmap (multiple rows)
    n_cols = 100
    n_rows = (n_bins + n_cols - 1) // n_cols
    heatmap_data = np.zeros((n_rows, n_cols))

    for i in range(n_bins):
        row = i // n_cols
        col = i % n_cols
        heatmap_data[row, col] = coverage_array[i]

    # Plot
    im = ax.imshow(heatmap_data, aspect="auto", cmap="YlOrRd", interpolation="nearest")
    ax.set_title("Genome Coverage Heatmap", fontweight="bold")
    ax.set_xlabel("Genome Position (100 bins per row)")
    ax.set_ylabel("Row")
    ax.set_xticks([])
    ax.set_yticks([])

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, orientation="horizontal", pad=0.05, aspect=30)
    cbar.set_label("Binding Sites per 10kb Bin")

    # Statistics text
    covered = np.sum(coverage_array > 0)
    coverage_pct = covered / n_bins * 100
    mean_sites = np.mean(coverage_array[coverage_array > 0]) if covered > 0 else 0

    stats_text = f"Coverage: {coverage_pct:.1f}%  |  Mean sites/bin: {mean_sites:.1f}"
    ax.text(0.5, -0.15, stats_text, transform=ax.transAxes, ha="center", va="top", fontsize=10)


def plot_primer_distribution(ax, positions_dict: Dict, genome_length: int):
    """
    Plot distribution of primers across genome.

    Shows how primers are spatially distributed.
    """
    # Collect all positions
    all_positions = []
    for primer, pos_data in positions_dict.items():
        all_positions.extend(pos_data["+"])
        all_positions.extend(pos_data["-"])

    all_positions = sorted(all_positions)

    # Calculate inter-site distances
    if len(all_positions) > 1:
        distances = [all_positions[i + 1] - all_positions[i] for i in range(len(all_positions) - 1)]
    else:
        distances = []

    # Plot histogram
    if distances:
        ax.hist(distances, bins=50, color="steelblue", edgecolor="black", alpha=0.7)
        ax.axvline(
            np.median(distances),
            color="red",
            linestyle="--",
            linewidth=2,
            label=f"Median: {np.median(distances):,.0f} bp",
        )
        ax.axvline(
            np.mean(distances),
            color="orange",
            linestyle="--",
            linewidth=2,
            label=f"Mean: {np.mean(distances):,.0f} bp",
        )

    ax.set_title("Inter-Site Distance Distribution", fontweight="bold")
    ax.set_xlabel("Distance (bp)")
    ax.set_ylabel("Frequency")
    ax.legend()
    ax.grid(alpha=0.3)


def plot_enrichment(ax, result, simulator):
    """
    Plot enrichment comparison.

    Bar chart showing target vs background amplification.
    """
    target_amp = result.target_amplification
    bg_amp = result.background_amplification

    # Use log scale for better visualization
    target_log = np.log10(target_amp + 1)
    bg_log = np.log10(bg_amp + 1)

    bars = ax.bar(
        ["Target", "Background"],
        [target_log, bg_log],
        color=["#2ecc71", "#e74c3c"],
        edgecolor="black",
        linewidth=1.5,
    )

    ax.set_title("Amplification Comparison (log scale)", fontweight="bold")
    ax.set_ylabel("log10(Amplification)")
    ax.grid(alpha=0.3, axis="y")

    # Add value labels on bars
    for bar, val, orig_val in zip(bars, [target_log, bg_log], [target_amp, bg_amp]):
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            height,
            f"{orig_val:.1e}",
            ha="center",
            va="bottom",
            fontweight="bold",
        )

    # Add enrichment text
    enrichment_text = f"Enrichment: {result.enrichment:.0f}×"
    ax.text(
        0.5,
        0.95,
        enrichment_text,
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=12,
        fontweight="bold",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
    )


def plot_metrics_summary(ax, result):
    """
    Plot key metrics as dashboard.

    Shows coverage, uniformity, specificity, and recommendation.
    """
    ax.axis("off")

    # Metrics to display
    metrics = [
        ("Coverage", result.target_coverage, 0.5),
        ("Uniformity", result.target_uniformity, 0.5),
        ("Specificity", result.specificity_score, 0.65),
    ]

    y_pos = 0.9
    for name, value, threshold in metrics:
        # Determine color
        if value >= threshold + 0.2:
            color = "#2ecc71"  # Green
            status = "Excellent"
        elif value >= threshold:
            color = "#f39c12"  # Orange
            status = "Good"
        else:
            color = "#e74c3c"  # Red
            status = "Poor"

        # Plot metric
        ax.text(0.05, y_pos, f"{name}:", fontsize=12, fontweight="bold", va="center")
        ax.text(
            0.95,
            y_pos,
            f"{value:.1%}",
            fontsize=12,
            ha="right",
            va="center",
            color=color,
            fontweight="bold",
        )

        # Progress bar
        bar_y = y_pos - 0.05
        bar_width = 0.9
        bar_height = 0.03

        # Background
        ax.add_patch(
            Rectangle(
                (0.05, bar_y), bar_width, bar_height, facecolor="lightgray", edgecolor="black"
            )
        )

        # Filled portion
        ax.add_patch(
            Rectangle(
                (0.05, bar_y), bar_width * value, bar_height, facecolor=color, edgecolor="black"
            )
        )

        y_pos -= 0.25

    # Recommendation
    rec_colors = {"EXCELLENT": "#2ecc71", "GOOD": "#f39c12", "FAIR": "#e67e22", "POOR": "#e74c3c"}
    rec_color = rec_colors.get(result.recommendation, "gray")

    ax.text(0.5, 0.15, "RECOMMENDATION", fontsize=10, ha="center", va="top", fontweight="bold")
    ax.text(
        0.5,
        0.05,
        result.recommendation,
        fontsize=16,
        ha="center",
        va="top",
        fontweight="bold",
        color=rec_color,
        bbox=dict(boxstyle="round", facecolor=rec_color, alpha=0.2, edgecolor=rec_color),
    )


def plot_primer_contributions(ax, primer_contributions: List):
    """
    Plot primer contributions as horizontal bar chart.

    Shows top contributing primers.
    """
    if not primer_contributions:
        ax.text(
            0.5,
            0.5,
            "No primer contribution data",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        ax.axis("off")
        return

    # Extract data
    primers = [p.primer[:12] for p in primer_contributions]  # Truncate long primers
    scores = [p.contribution_score for p in primer_contributions]
    specificities = [min(p.specificity, 1000) for p in primer_contributions]  # Cap for color

    # Create bars
    y_pos = np.arange(len(primers))
    colors = plt.cm.RdYlGn(np.array(specificities) / max(specificities))

    bars = ax.barh(y_pos, scores, color=colors, edgecolor="black", linewidth=0.5)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(primers, fontsize=8)
    ax.set_xlabel("Contribution Score")
    ax.set_title("Top Primer Contributions", fontweight="bold")
    ax.set_xlim(0, 1.0)
    ax.grid(alpha=0.3, axis="x")

    # Add value labels
    for i, (bar, score) in enumerate(zip(bars, scores)):
        width = bar.get_width()
        ax.text(
            width + 0.02,
            bar.get_y() + bar.get_height() / 2,
            f"{score:.2f}",
            ha="left",
            va="center",
            fontsize=7,
        )


def _severity_from_size(length):
    """Bucket a gap by length, using the same thresholds as simulation_analysis.

    Kept in step with `SimulationAnalyzer._analyze_coverage` so an unlabelled
    gap lands in the bucket it would have been given there.
    """
    if length > 100_000:
        return "critical"
    if length > 50_000:
        return "high"
    if length > 20_000:
        return "medium"
    return "low"


def plot_gaps(ax, gaps: List[Dict]):
    """
    Plot gap severity distribution.

    Shows distribution of gaps by severity.
    """
    if not gaps:
        ax.text(
            0.5,
            0.5,
            "No gaps detected\n(excellent coverage!)",
            ha="center",
            va="center",
            fontsize=12,
            fontweight="bold",
            color="green",
            transform=ax.transAxes,
        )
        ax.axis("off")
        return

    # Count by severity.
    #
    # Two gap shapes exist in this project. `simulation_analysis` labels each
    # gap with a severity; `SwgaSimulator._identify_gaps` emits
    # {start, end, length, mean_coverage} and does not. Only the labelled shape
    # is passed here today, but indexing the key directly meant handing over the
    # other one -- the obvious thing to try -- raised a bare KeyError from
    # inside the plotting code, in the panel whose whole job is to report
    # coverage problems. Fall back to bucketing by size on the same thresholds
    # simulation_analysis uses.
    severity_counts = {"critical": 0, "high": 0, "medium": 0, "low": 0}
    for gap in gaps:
        severity = gap.get("severity")
        if severity not in severity_counts:
            severity = _severity_from_size(gap.get("length", gap.get("size", 0)))
        severity_counts[severity] += 1

    # Plot pie chart
    labels = []
    sizes = []
    colors_map = {"critical": "#e74c3c", "high": "#e67e22", "medium": "#f39c12", "low": "#f1c40f"}
    colors = []

    for severity in ["critical", "high", "medium", "low"]:
        if severity_counts[severity] > 0:
            labels.append(f"{severity.capitalize()} ({severity_counts[severity]})")
            sizes.append(severity_counts[severity])
            colors.append(colors_map[severity])

    if sizes:
        wedges, texts, autotexts = ax.pie(
            sizes, labels=labels, colors=colors, autopct="%1.1f%%", startangle=90
        )
        for autotext in autotexts:
            autotext.set_color("white")
            autotext.set_fontweight("bold")

    ax.set_title("Gap Severity Distribution", fontweight="bold")


def generate_interactive_html(
    result, simulator, analysis, output_file="simulation_interactive.html"
):
    """
    Generate interactive HTML visualization.

    Uses plotly for interactive plots.
    """
    logger.info("Generating interactive HTML visualization...")

    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError:
        logger.warning("plotly not installed - skipping interactive HTML generation")
        logger.info("Install with: pip install plotly")
        return

    # Create figure with subplots
    fig = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=("Coverage Heatmap", "Primer Distribution", "Enrichment", "Metrics"),
        specs=[
            [{"type": "heatmap"}, {"type": "histogram"}],
            [{"type": "bar"}, {"type": "indicator"}],
        ],
    )

    # Add plots (implementation details omitted for brevity)
    # ... plotly-specific code ...

    # Save
    fig.write_html(output_file)
    logger.info(f"Interactive HTML saved to: {output_file}")


if __name__ == "__main__":
    # Test with dummy data
    print("This module is meant to be imported, not run directly")
    print("Usage: from neoswga.core.simulation_plots import generate_plots")
