"""
Pareto frontier visualization for SWGA primer set size optimization.

Provides functions to visualize the coverage vs fg/bg ratio tradeoff
and generate reports for the Pareto frontier analysis.

Usage:
    from neoswga.core.pareto_frontier import plot_frontier, generate_frontier_report

    # Plot the frontier
    fig = plot_frontier(frontier_result, application='clinical')
    fig.savefig('pareto_frontier.png')

    # Generate text report
    report = generate_frontier_report(frontier_result, application='clinical')
    print(report)
"""

import logging
from typing import List, Optional, Dict, Any, TYPE_CHECKING

if TYPE_CHECKING:
    from neoswga.core.set_size_optimizer import SetSizeMetrics, FrontierResult

logger = logging.getLogger(__name__)


# Application zone definitions for visualization
APPLICATION_ZONES = {
    'discovery': {
        'label': 'Discovery Zone',
        'color': '#90EE90',  # Light green
        'min_coverage': 0.85,
        'min_ratio': 1.5,
        'description': 'High coverage, lower specificity',
    },
    'clinical': {
        'label': 'Clinical Zone',
        'color': '#FFB6C1',  # Light pink
        'min_coverage': 0.60,
        'min_ratio': 8.0,
        'description': 'High specificity, moderate coverage',
    },
    'enrichment': {
        'label': 'Enrichment Zone',
        'color': '#ADD8E6',  # Light blue
        'min_coverage': 0.75,
        'min_ratio': 4.0,
        'description': 'Balanced coverage and specificity',
    },
    'metagenomics': {
        'label': 'Metagenomics Zone',
        'color': '#FFFACD',  # Light yellow
        'min_coverage': 0.90,
        'min_ratio': 1.0,
        'description': 'Maximum coverage, lowest specificity threshold',
    },
}


def plot_frontier(
    frontier_result: 'FrontierResult',
    application: str = 'enrichment',
    show_zones: bool = True,
    show_all_points: bool = True,
    figsize: tuple = (10, 8),
    title: Optional[str] = None,
) -> Any:
    """
    Create a matplotlib plot of the Pareto frontier.

    Args:
        frontier_result: FrontierResult from ParetoFrontierGenerator
        application: Current application profile (for zone highlighting)
        show_zones: Whether to show application zone overlays
        show_all_points: Whether to show non-Pareto-optimal points
        figsize: Figure size in inches
        title: Optional custom title

    Returns:
        matplotlib Figure object

    Raises:
        ImportError: If matplotlib is not available
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
    except ImportError:
        raise ImportError("matplotlib is required for plotting. Install with: pip install matplotlib")

    fig, ax = plt.subplots(figsize=figsize)

    pareto = frontier_result.pareto_points
    all_points = frontier_result.all_points
    selected = frontier_result.selected_point

    if not pareto:
        ax.text(0.5, 0.5, 'No Pareto points available',
                ha='center', va='center', transform=ax.transAxes)
        return fig

    # Extract data
    pareto_coverage = [p.fg_coverage for p in pareto]
    pareto_ratio = [p.fg_bg_ratio for p in pareto]
    pareto_sizes = [p.set_size for p in pareto]

    # Plot application zones as background regions
    if show_zones:
        max_ratio = max(pareto_ratio) * 1.1 if pareto_ratio else 20
        zone = APPLICATION_ZONES.get(application, APPLICATION_ZONES['enrichment'])
        ax.axhspan(zone['min_ratio'], max_ratio, alpha=0.2, color=zone['color'],
                   label=f"{zone['label']}")
        ax.axvspan(zone['min_coverage'], 1.0, alpha=0.1, color=zone['color'])

    # Plot non-Pareto points (faded)
    if show_all_points and all_points:
        non_pareto = [p for p in all_points if not p.is_pareto_optimal]
        if non_pareto:
            np_coverage = [p.fg_coverage for p in non_pareto]
            np_ratio = [p.fg_bg_ratio for p in non_pareto]
            ax.scatter(np_coverage, np_ratio, c='gray', alpha=0.3, s=50,
                       label='Non-optimal', marker='o')

    # Plot Pareto frontier line
    sorted_pareto = sorted(zip(pareto_coverage, pareto_ratio, pareto_sizes),
                           key=lambda x: x[0])
    cov_sorted = [x[0] for x in sorted_pareto]
    ratio_sorted = [x[1] for x in sorted_pareto]

    ax.plot(cov_sorted, ratio_sorted, 'b-', alpha=0.5, linewidth=2)

    # Plot Pareto points with size labels
    scatter = ax.scatter(pareto_coverage, pareto_ratio, c='blue', s=100,
                         label='Pareto optimal', marker='o', edgecolors='darkblue',
                         linewidths=1.5, zorder=5)

    # Add size labels
    for cov, ratio, size in sorted_pareto:
        ax.annotate(str(size), (cov, ratio),
                    textcoords="offset points", xytext=(5, 5),
                    fontsize=9, color='darkblue')

    # Highlight selected point
    if selected:
        ax.scatter([selected.fg_coverage], [selected.fg_bg_ratio],
                   c='red', s=200, marker='*', label='Selected',
                   edgecolors='darkred', linewidths=2, zorder=10)

    # Labels and title
    ax.set_xlabel('Target Genome Coverage', fontsize=12)
    ax.set_ylabel('Foreground/Background Ratio (fg/bg)', fontsize=12)

    if title:
        ax.set_title(title, fontsize=14)
    else:
        ax.set_title(f'Coverage vs Specificity Pareto Frontier\n({application.capitalize()} Application)',
                     fontsize=14)

    # Set axis limits with some padding
    ax.set_xlim(0, 1.05)
    if pareto_ratio:
        ax.set_ylim(0, max(pareto_ratio) * 1.15)

    # Legend
    ax.legend(loc='upper left')

    # Grid
    ax.grid(True, alpha=0.3)

    # Add explanation text
    explanation = frontier_result.selection_explanation
    if explanation:
        ax.text(0.02, 0.02, explanation,
                transform=ax.transAxes, fontsize=9,
                verticalalignment='bottom', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    return fig


def generate_frontier_report(
    frontier_result: 'FrontierResult',
    application: str = 'enrichment',
    include_all_points: bool = False,
) -> str:
    """
    Generate a text report of the Pareto frontier analysis.

    Args:
        frontier_result: FrontierResult from ParetoFrontierGenerator
        application: Application profile for context
        include_all_points: Whether to include non-Pareto-optimal points

    Returns:
        Formatted text report
    """
    lines = []
    lines.append("=" * 60)
    lines.append("Pareto Frontier Analysis Report")
    lines.append("=" * 60)
    lines.append("")

    # Summary
    pareto = frontier_result.pareto_points
    all_points = frontier_result.all_points
    selected = frontier_result.selected_point

    lines.append(f"Application profile: {application}")
    lines.append(f"Total evaluated sizes: {len(all_points)}")
    lines.append(f"Pareto-optimal points: {len(pareto)}")
    lines.append("")

    # Pareto frontier table
    lines.append("Pareto Frontier (sorted by set size):")
    lines.append("-" * 60)
    lines.append(f"{'Size':>6}  {'Coverage':>10}  {'fg/bg Ratio':>12}  {'fg Sites':>10}")
    lines.append("-" * 60)

    for point in sorted(pareto, key=lambda p: p.set_size):
        marker = " <-- SELECTED" if selected and point.set_size == selected.set_size else ""
        lines.append(
            f"{point.set_size:>6}  {point.fg_coverage:>10.1%}  {point.fg_bg_ratio:>12.1f}  "
            f"{point.fg_binding_sites:>10,}{marker}"
        )

    lines.append("-" * 60)
    lines.append("")

    # Selection explanation
    if frontier_result.selection_explanation:
        lines.append("Selection Decision:")
        lines.append(frontier_result.selection_explanation)
        lines.append("")

    # Selected point details
    if selected:
        lines.append("Selected Set Details:")
        lines.append("-" * 40)
        lines.append(f"  Set size: {selected.set_size} primers")
        lines.append(f"  Target coverage: {selected.fg_coverage:.1%}")
        lines.append(f"  Background coverage: {selected.bg_coverage:.1%}")
        lines.append(f"  fg/bg ratio: {selected.fg_bg_ratio:.1f}")
        lines.append(f"  Foreground binding sites: {selected.fg_binding_sites:,}")
        lines.append(f"  Background binding sites: {selected.bg_binding_sites:,}")

        if selected.primers:
            lines.append("")
            lines.append(f"  Primers ({len(selected.primers)}):")
            for i, primer in enumerate(selected.primers[:10], 1):  # Limit to first 10
                lines.append(f"    {i}. {primer}")
            if len(selected.primers) > 10:
                lines.append(f"    ... and {len(selected.primers) - 10} more")

        lines.append("")

    # Non-Pareto points (optional)
    if include_all_points:
        non_pareto = [p for p in all_points if not p.is_pareto_optimal]
        if non_pareto:
            lines.append("Non-Pareto-Optimal Points (dominated):")
            lines.append("-" * 60)
            lines.append(f"{'Size':>6}  {'Coverage':>10}  {'fg/bg Ratio':>12}")
            lines.append("-" * 60)

            for point in sorted(non_pareto, key=lambda p: p.set_size):
                lines.append(
                    f"{point.set_size:>6}  {point.fg_coverage:>10.1%}  {point.fg_bg_ratio:>12.1f}"
                )
            lines.append("")

    # Application zone context
    zone = APPLICATION_ZONES.get(application, APPLICATION_ZONES['enrichment'])
    lines.append("Application Zone Thresholds:")
    lines.append(f"  Minimum coverage target: {zone['min_coverage']:.0%}")
    lines.append(f"  Minimum fg/bg ratio: {zone['min_ratio']:.1f}")
    lines.append(f"  {zone['description']}")
    lines.append("")

    # Recommendations
    lines.append("Recommendations:")
    lines.append("-" * 40)

    if selected:
        if selected.fg_coverage >= zone['min_coverage'] and selected.fg_bg_ratio >= zone['min_ratio']:
            lines.append("  The selected set meets all application requirements.")
        elif selected.fg_coverage < zone['min_coverage']:
            lines.append(f"  Warning: Coverage ({selected.fg_coverage:.1%}) is below target ({zone['min_coverage']:.0%}).")
            lines.append("  Consider using more primers or different candidates.")
        elif selected.fg_bg_ratio < zone['min_ratio']:
            lines.append(f"  Warning: fg/bg ratio ({selected.fg_bg_ratio:.1f}) is below target ({zone['min_ratio']:.1f}).")
            lines.append("  Consider using fewer, more selective primers.")
    else:
        lines.append("  No set was selected. Review the Pareto frontier manually.")

    lines.append("")
    lines.append("=" * 60)

    return "\n".join(lines)


def generate_frontier_json(
    frontier_result: 'FrontierResult',
    application: str = 'enrichment',
) -> Dict[str, Any]:
    """
    Generate JSON-serializable representation of frontier analysis.

    Args:
        frontier_result: FrontierResult from ParetoFrontierGenerator
        application: Application profile for context

    Returns:
        Dictionary suitable for JSON serialization
    """
    zone = APPLICATION_ZONES.get(application, APPLICATION_ZONES['enrichment'])

    result = {
        'application': application,
        'zone_thresholds': {
            'min_coverage': zone['min_coverage'],
            'min_fg_bg_ratio': zone['min_ratio'],
        },
        'frontier': frontier_result.to_dict(),
        'summary': {
            'total_evaluated': len(frontier_result.all_points),
            'pareto_optimal_count': len(frontier_result.pareto_points),
        },
    }

    if frontier_result.selected_point:
        selected = frontier_result.selected_point
        result['summary']['selected'] = {
            'set_size': selected.set_size,
            'fg_coverage': selected.fg_coverage,
            'fg_bg_ratio': selected.fg_bg_ratio,
            'meets_coverage_target': selected.fg_coverage >= zone['min_coverage'],
            'meets_ratio_target': selected.fg_bg_ratio >= zone['min_ratio'],
        }

    return result


def summarize_frontier_for_cli(
    frontier_result: 'FrontierResult',
    application: str = 'enrichment',
) -> str:
    """
    Generate a concise CLI-friendly summary of the frontier.

    Args:
        frontier_result: FrontierResult from ParetoFrontierGenerator
        application: Application profile for context

    Returns:
        Short summary string suitable for CLI output
    """
    pareto = frontier_result.pareto_points
    selected = frontier_result.selected_point

    if not pareto:
        return "No Pareto-optimal points found."

    # Size range on frontier
    min_size = min(p.set_size for p in pareto)
    max_size = max(p.set_size for p in pareto)

    # Coverage range
    min_cov = min(p.fg_coverage for p in pareto)
    max_cov = max(p.fg_coverage for p in pareto)

    # Ratio range
    min_ratio = min(p.fg_bg_ratio for p in pareto)
    max_ratio = max(p.fg_bg_ratio for p in pareto)

    lines = [
        f"Pareto frontier: {len(pareto)} optimal points",
        f"  Size range: {min_size}-{max_size} primers",
        f"  Coverage range: {min_cov:.0%}-{max_cov:.0%}",
        f"  fg/bg ratio range: {min_ratio:.1f}-{max_ratio:.1f}",
    ]

    if selected:
        lines.append(f"Selected for '{application}': {selected.set_size} primers "
                     f"({selected.fg_coverage:.0%} coverage, {selected.fg_bg_ratio:.1f}x ratio)")

    return "\n".join(lines)
