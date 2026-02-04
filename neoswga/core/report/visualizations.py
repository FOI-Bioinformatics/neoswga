"""
Interactive Plotly visualizations for SWGA reports.

Provides interactive charts that can be embedded in HTML reports:
- Filtering funnel chart
- Quality component radar chart
- Tm/GC distribution plots
- Coverage vs specificity scatter

All functions return empty strings when Plotly is not installed,
allowing graceful degradation.
"""

from typing import List, Optional, Dict, Any, TYPE_CHECKING
import logging

logger = logging.getLogger(__name__)

# Check Plotly availability at import time
_PLOTLY_AVAILABLE = False
try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    _PLOTLY_AVAILABLE = True
except ImportError:
    go = None  # type: ignore
    make_subplots = None  # type: ignore

if TYPE_CHECKING:
    from neoswga.core.report.quality import GradeComponent
    from neoswga.core.report.metrics import PrimerMetrics


def is_plotly_available() -> bool:
    """
    Check if Plotly is installed and available.

    Returns:
        True if Plotly can be imported, False otherwise
    """
    return _PLOTLY_AVAILABLE


# Default color scheme for charts
CHART_COLORS = {
    "primary": "#2c5282",
    "secondary": "#4299e1",
    "success": "#48bb78",
    "warning": "#ecc94b",
    "danger": "#f56565",
    "info": "#4299e1",
    "muted": "#a0aec0",
    # Gradient for funnel
    "funnel_gradient": [
        "#2c5282", "#3182ce", "#4299e1", "#63b3ed", "#90cdf4", "#bee3f8"
    ],
    # Radar chart fill
    "radar_fill": "rgba(66, 153, 225, 0.3)",
    "radar_line": "#2c5282",
    # Rating colors
    "excellent": "#48bb78",
    "good": "#4299e1",
    "acceptable": "#ecc94b",
    "poor": "#f56565",
}


def render_filtering_funnel(
    funnel_data: List[tuple],
    include_plotlyjs: str = 'cdn',
    height: int = 400,
) -> str:
    """
    Render filtering funnel as interactive Plotly funnel chart.

    Shows the progressive filtering of primer candidates through
    the pipeline stages.

    Args:
        funnel_data: List of (stage_name, count) tuples
        include_plotlyjs: How to include Plotly.js ('cdn', 'inline', False)
        height: Chart height in pixels

    Returns:
        HTML string with embedded Plotly chart, or empty string if
        Plotly is not available or funnel_data is empty
    """
    if not _PLOTLY_AVAILABLE:
        logger.debug("Plotly not available, skipping funnel chart")
        return ""

    if not funnel_data:
        return ""

    stages = [item[0] for item in funnel_data]
    values = [item[1] for item in funnel_data]

    # Calculate retention percentages
    initial = values[0] if values else 1
    percentages = [f"{(v / initial * 100):.1f}%" for v in values]

    fig = go.Figure(go.Funnel(
        y=stages,
        x=values,
        textposition="inside",
        textinfo="value+percent initial",
        texttemplate="%{x:,.0f}<br>(%{percentInitial:.1%})",
        marker=dict(
            color=CHART_COLORS["funnel_gradient"][:len(stages)],
        ),
        connector=dict(
            line=dict(color=CHART_COLORS["muted"], width=1)
        ),
    ))

    fig.update_layout(
        title=dict(
            text="Primer Filtering Pipeline",
            font=dict(size=16, color="#2d3748"),
        ),
        showlegend=False,
        height=height,
        margin=dict(l=20, r=20, t=50, b=20),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(family="-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif"),
    )

    return fig.to_html(
        full_html=False,
        include_plotlyjs=include_plotlyjs,
        config={'displayModeBar': False, 'responsive': True}
    )


def render_component_radar(
    components: List['GradeComponent'],
    include_plotlyjs: str = 'cdn',
    height: int = 400,
) -> str:
    """
    Render quality components as interactive radar chart.

    Displays normalized scores for each quality component
    (Coverage, Specificity, Uniformity, etc.) on a polar plot.

    Args:
        components: List of GradeComponent objects from quality assessment
        include_plotlyjs: How to include Plotly.js ('cdn', 'inline', False)
        height: Chart height in pixels

    Returns:
        HTML string with embedded Plotly chart, or empty string if
        Plotly is not available or components is empty
    """
    if not _PLOTLY_AVAILABLE:
        logger.debug("Plotly not available, skipping radar chart")
        return ""

    if not components:
        return ""

    # Extract data from components
    names = [c.name for c in components]
    scores = [c.normalized_score for c in components]
    ratings = [c.rating for c in components]

    # Close the polygon by repeating the first point
    names_closed = names + [names[0]]
    scores_closed = scores + [scores[0]]

    # Create hover text with ratings
    hover_text = [
        f"{name}<br>Score: {score:.2f}<br>Rating: {rating}"
        for name, score, rating in zip(names, scores, ratings)
    ]
    hover_text_closed = hover_text + [hover_text[0]]

    fig = go.Figure()

    fig.add_trace(go.Scatterpolar(
        r=scores_closed,
        theta=names_closed,
        fill='toself',
        fillcolor=CHART_COLORS["radar_fill"],
        line=dict(color=CHART_COLORS["radar_line"], width=2),
        name='Quality Score',
        hovertext=hover_text_closed,
        hoverinfo="text",
    ))

    fig.update_layout(
        polar=dict(
            radialaxis=dict(
                visible=True,
                range=[0, 1],
                tickvals=[0.25, 0.5, 0.75, 1.0],
                ticktext=["Poor", "Fair", "Good", "Excellent"],
                gridcolor="#e2e8f0",
                linecolor="#e2e8f0",
            ),
            angularaxis=dict(
                gridcolor="#e2e8f0",
                linecolor="#e2e8f0",
            ),
            bgcolor="rgba(0,0,0,0)",
        ),
        title=dict(
            text="Quality Component Analysis",
            font=dict(size=16, color="#2d3748"),
        ),
        showlegend=False,
        height=height,
        margin=dict(l=60, r=60, t=50, b=40),
        paper_bgcolor="rgba(0,0,0,0)",
        font=dict(family="-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif"),
    )

    return fig.to_html(
        full_html=False,
        include_plotlyjs=include_plotlyjs,
        config={'displayModeBar': False, 'responsive': True}
    )


def render_tm_gc_distribution(
    primers: List['PrimerMetrics'],
    reaction_temp: float = 30.0,
    include_plotlyjs: str = 'cdn',
    height: int = 350,
) -> str:
    """
    Render Tm histogram and Tm vs GC scatter as side-by-side plots.

    Shows the distribution of melting temperatures and the
    relationship between Tm and GC content.

    Args:
        primers: List of PrimerMetrics objects
        reaction_temp: Reaction temperature for reference line
        include_plotlyjs: How to include Plotly.js ('cdn', 'inline', False)
        height: Chart height in pixels

    Returns:
        HTML string with embedded Plotly chart, or empty string if
        Plotly is not available or primers is empty
    """
    if not _PLOTLY_AVAILABLE:
        logger.debug("Plotly not available, skipping Tm/GC distribution")
        return ""

    if not primers:
        return ""

    # Extract data
    tms = [p.tm for p in primers if p.tm > 0]
    gcs = [p.gc_content * 100 for p in primers if p.tm > 0]
    sequences = [p.sequence for p in primers if p.tm > 0]

    if not tms:
        return ""

    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=('Tm Distribution', 'Tm vs GC Content'),
        horizontal_spacing=0.12,
    )

    # Tm histogram
    fig.add_trace(
        go.Histogram(
            x=tms,
            nbinsx=min(20, len(tms)),
            marker_color=CHART_COLORS["primary"],
            name="Tm",
            hovertemplate="Tm: %{x:.1f}C<br>Count: %{y}<extra></extra>",
        ),
        row=1, col=1
    )

    # Add reaction temperature reference line
    fig.add_vline(
        x=reaction_temp,
        line_dash="dash",
        line_color=CHART_COLORS["danger"],
        annotation_text=f"Reaction: {reaction_temp}C",
        annotation_position="top",
        row=1, col=1
    )

    # Tm vs GC scatter
    hover_text = [
        f"Seq: {seq}<br>Tm: {tm:.1f}C<br>GC: {gc:.0f}%"
        for seq, tm, gc in zip(sequences, tms, gcs)
    ]

    fig.add_trace(
        go.Scatter(
            x=gcs,
            y=tms,
            mode='markers',
            marker=dict(
                color=CHART_COLORS["secondary"],
                size=10,
                line=dict(width=1, color=CHART_COLORS["primary"]),
            ),
            name="Primers",
            hovertext=hover_text,
            hoverinfo="text",
        ),
        row=1, col=2
    )

    fig.update_layout(
        height=height,
        showlegend=False,
        margin=dict(l=50, r=30, t=60, b=50),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(family="-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif"),
    )

    fig.update_xaxes(title_text="Tm (C)", gridcolor="#e2e8f0", row=1, col=1)
    fig.update_yaxes(title_text="Count", gridcolor="#e2e8f0", row=1, col=1)
    fig.update_xaxes(title_text="GC Content (%)", gridcolor="#e2e8f0", row=1, col=2)
    fig.update_yaxes(title_text="Tm (C)", gridcolor="#e2e8f0", row=1, col=2)

    return fig.to_html(
        full_html=False,
        include_plotlyjs=include_plotlyjs,
        config={'displayModeBar': False, 'responsive': True}
    )


def render_coverage_specificity_scatter(
    primers: List['PrimerMetrics'],
    genome_size: int = 0,
    include_plotlyjs: str = 'cdn',
    height: int = 400,
) -> str:
    """
    Render coverage vs specificity scatter with Pareto frontier.

    Visualizes the trade-off between coverage (binding sites) and
    specificity (target/background ratio) for each primer.

    Args:
        primers: List of PrimerMetrics objects
        genome_size: Target genome size for density calculation
        include_plotlyjs: How to include Plotly.js ('cdn', 'inline', False)
        height: Chart height in pixels

    Returns:
        HTML string with embedded Plotly chart, or empty string if
        Plotly is not available or primers is empty
    """
    if not _PLOTLY_AVAILABLE:
        logger.debug("Plotly not available, skipping coverage/specificity scatter")
        return ""

    if not primers:
        return ""

    # Extract data
    # Coverage proxy: foreground binding sites (or density if genome_size given)
    if genome_size > 0:
        coverage_values = [(p.fg_sites / genome_size) * 1_000_000 for p in primers]
        x_label = "Binding Density (sites/Mbp)"
    else:
        coverage_values = [p.fg_sites for p in primers]
        x_label = "Foreground Binding Sites"

    # Specificity: ratio of fg_freq to bg_freq (capped for display)
    specificity_values = [min(p.specificity, 10000) for p in primers]
    sequences = [p.sequence for p in primers]

    # Color by quality (amp_pred if available)
    colors = [p.amp_pred if p.amp_pred > 0 else 0.5 for p in primers]

    # Create hover text
    hover_text = [
        f"Seq: {seq}<br>"
        f"FG Sites: {p.fg_sites}<br>"
        f"Specificity: {p.specificity:.0f}x<br>"
        f"Quality: {p.amp_pred:.2f}"
        for seq, p in zip(sequences, primers)
    ]

    fig = go.Figure()

    # Main scatter plot
    fig.add_trace(go.Scatter(
        x=coverage_values,
        y=specificity_values,
        mode='markers',
        marker=dict(
            size=12,
            color=colors,
            colorscale='Viridis',
            showscale=True,
            colorbar=dict(
                title=dict(text="Quality", font=dict(size=11)),
                tickfont=dict(size=10),
            ),
            line=dict(width=1, color='white'),
        ),
        text=hover_text,
        hoverinfo="text",
        name="Primers",
    ))

    # Calculate and add Pareto frontier
    pareto_points = _calculate_pareto_frontier(
        list(zip(coverage_values, specificity_values))
    )
    if len(pareto_points) > 1:
        pareto_x = [p[0] for p in pareto_points]
        pareto_y = [p[1] for p in pareto_points]
        fig.add_trace(go.Scatter(
            x=pareto_x,
            y=pareto_y,
            mode='lines',
            line=dict(color=CHART_COLORS["success"], width=2, dash='dash'),
            name="Pareto Frontier",
            hoverinfo="skip",
        ))

    fig.update_layout(
        title=dict(
            text="Coverage vs Specificity Trade-off",
            font=dict(size=16, color="#2d3748"),
        ),
        xaxis_title=x_label,
        yaxis_title="Specificity (fg/bg ratio)",
        height=height,
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01,
            bgcolor="rgba(255,255,255,0.8)",
        ),
        margin=dict(l=60, r=30, t=60, b=50),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(family="-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif"),
    )

    fig.update_xaxes(gridcolor="#e2e8f0", zeroline=False)
    fig.update_yaxes(gridcolor="#e2e8f0", zeroline=False, type="log")

    return fig.to_html(
        full_html=False,
        include_plotlyjs=include_plotlyjs,
        config={'displayModeBar': False, 'responsive': True}
    )


def _calculate_pareto_frontier(points: List[tuple]) -> List[tuple]:
    """
    Calculate Pareto frontier for a set of 2D points.

    A point is Pareto-optimal if no other point is better in both dimensions.

    Args:
        points: List of (x, y) tuples

    Returns:
        List of Pareto-optimal points sorted by x coordinate
    """
    if not points:
        return []

    # Sort by x descending
    sorted_points = sorted(points, key=lambda p: (-p[0], -p[1]))

    pareto = []
    max_y = float('-inf')

    for point in sorted_points:
        if point[1] > max_y:
            pareto.append(point)
            max_y = point[1]

    # Sort by x ascending for plotting
    return sorted(pareto, key=lambda p: p[0])


def render_primer_heatmap(
    primers: List['PrimerMetrics'],
    include_plotlyjs: str = 'cdn',
    height: int = 400,
) -> str:
    """
    Render primer metrics as a heatmap for comparison.

    Shows normalized values for key metrics across all primers.

    Args:
        primers: List of PrimerMetrics objects
        include_plotlyjs: How to include Plotly.js ('cdn', 'inline', False)
        height: Chart height in pixels

    Returns:
        HTML string with embedded Plotly chart, or empty string if
        Plotly is not available or primers is empty
    """
    if not _PLOTLY_AVAILABLE:
        logger.debug("Plotly not available, skipping primer heatmap")
        return ""

    if not primers or len(primers) < 2:
        return ""

    # Metrics to display (normalized 0-1)
    metrics = ['GC%', 'Tm', 'Specificity', 'Uniformity', 'Quality']

    # Build data matrix
    z_data = []
    y_labels = []

    for p in primers[:15]:  # Limit to 15 primers for readability
        y_labels.append(p.sequence[:10] + "..." if len(p.sequence) > 10 else p.sequence)

        # Normalize each metric to 0-1 range
        gc_norm = p.gc_content  # Already 0-1
        tm_norm = min(max((p.tm - 20) / 40, 0), 1)  # Normalize 20-60C range
        spec_norm = min(p.specificity / 1000, 1)  # Normalize to 1000x
        uniform_norm = 1 - p.gini  # Invert gini
        quality_norm = p.amp_pred if p.amp_pred > 0 else 0.5

        z_data.append([gc_norm, tm_norm, spec_norm, uniform_norm, quality_norm])

    fig = go.Figure(data=go.Heatmap(
        z=z_data,
        x=metrics,
        y=y_labels,
        colorscale='RdYlGn',
        zmin=0,
        zmax=1,
        hovertemplate="Primer: %{y}<br>Metric: %{x}<br>Score: %{z:.2f}<extra></extra>",
    ))

    fig.update_layout(
        title=dict(
            text="Primer Metrics Comparison",
            font=dict(size=16, color="#2d3748"),
        ),
        height=height,
        margin=dict(l=100, r=30, t=60, b=50),
        paper_bgcolor="rgba(0,0,0,0)",
        font=dict(family="-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif"),
    )

    return fig.to_html(
        full_html=False,
        include_plotlyjs=include_plotlyjs,
        config={'displayModeBar': False, 'responsive': True}
    )


def _calculate_heterodimer_dg(seq1: str, seq2: str) -> float:
    """
    Calculate estimated delta G for heterodimer formation between two sequences.

    Uses a simplified nearest-neighbor approximation based on complementarity,
    focusing on the 3' region which is most critical for primer-dimer formation
    during extension.

    Args:
        seq1: First primer sequence (5' to 3')
        seq2: Second primer sequence (5' to 3')

    Returns:
        Estimated delta G in kcal/mol (negative = favorable binding)
    """
    if not seq1 or not seq2:
        return 0.0

    # Reverse complement mapping
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}

    # Get reverse complement of seq2
    seq2_rc = ''.join(complement.get(b, 'N') for b in reversed(seq2.upper()))
    seq1_upper = seq1.upper()

    # Nearest-neighbor delta G values (kcal/mol) at 37C, 1M Na+
    # Based on SantaLucia unified parameters
    nn_dg = {
        'AA': -1.00, 'TT': -1.00,
        'AT': -0.88, 'TA': -0.58,
        'CA': -1.45, 'TG': -1.45,
        'GT': -1.44, 'AC': -1.44,
        'CT': -1.28, 'AG': -1.28,
        'GA': -1.30, 'TC': -1.30,
        'CG': -2.17, 'GC': -2.24,
        'GG': -1.84, 'CC': -1.84,
    }

    # Find the most stable complementary alignment
    best_dg = 0.0
    len1, len2 = len(seq1_upper), len(seq2_rc)

    # Sliding window to find best alignment
    for offset in range(-len2 + 3, len1 - 2):  # Need at least 3bp overlap
        start1 = max(0, offset)
        end1 = min(len1, offset + len2)
        start2 = max(0, -offset)
        overlap_len = end1 - start1

        if overlap_len < 3:
            continue

        # Calculate delta G for this alignment
        alignment_dg = 0.0
        consecutive_matches = 0
        current_run_dg = 0.0
        best_run_dg = 0.0

        for i in range(overlap_len):
            base1 = seq1_upper[start1 + i]
            base2 = seq2_rc[start2 + i]

            if base1 == base2:
                # Watson-Crick base pair
                dinuc = base1 + base2
                bp_dg = nn_dg.get(dinuc, -1.5)
                current_run_dg += bp_dg
                consecutive_matches += 1
            else:
                # Mismatch breaks the run
                if consecutive_matches >= 3:
                    best_run_dg = min(best_run_dg, current_run_dg)
                consecutive_matches = 0
                current_run_dg = 0.0

        # Check final run
        if consecutive_matches >= 3:
            best_run_dg = min(best_run_dg, current_run_dg)

        alignment_dg = best_run_dg

        # Weight 3' complementarity more heavily (last 6 bases of seq1)
        if start1 >= len1 - 6 and alignment_dg < 0:
            alignment_dg *= 1.3  # 30% penalty for 3' binding

        best_dg = min(best_dg, alignment_dg)

    # Add initiation penalty (+1.96 kcal/mol for duplex initiation)
    if best_dg < 0:
        best_dg += 1.96

    return round(best_dg, 2)


def _calculate_self_dimer_dg(seq: str) -> float:
    """
    Calculate estimated delta G for self-dimer (hairpin/self-complementary) formation.

    Args:
        seq: Primer sequence (5' to 3')

    Returns:
        Estimated delta G in kcal/mol (negative = favorable binding)
    """
    return _calculate_heterodimer_dg(seq, seq)


def _build_dimer_matrix(
    primers: List['PrimerMetrics'],
    max_primers: int = 20,
) -> tuple:
    """
    Build a symmetric matrix of heterodimer delta G values.

    Args:
        primers: List of PrimerMetrics objects
        max_primers: Maximum number of primers to include

    Returns:
        Tuple of (matrix, labels, sequences) or (None, None, None) if invalid
    """
    # Limit primers for performance and readability
    primers_subset = primers[:max_primers]
    n = len(primers_subset)

    if n < 2:
        return None, None, None

    # Create labels (truncated sequences)
    labels = []
    sequences = []
    for p in primers_subset:
        label = p.sequence[:8] + "..." if len(p.sequence) > 8 else p.sequence
        labels.append(label)
        sequences.append(p.sequence)

    # Build symmetric matrix
    matrix = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i, n):
            if i == j:
                # Self-dimer
                dg = _calculate_self_dimer_dg(sequences[i])
            else:
                # Heterodimer
                dg = _calculate_heterodimer_dg(sequences[i], sequences[j])
            matrix[i][j] = dg
            matrix[j][i] = dg  # Symmetric

    return matrix, labels, sequences


def render_dimer_network_heatmap(
    primers: List['PrimerMetrics'],
    include_plotlyjs: str = 'cdn',
    height: int = 500,
    max_primers: int = 20,
    show_values: bool = True,
) -> str:
    """
    Render primer-primer dimer interaction network as an interactive heatmap.

    Shows the estimated delta G of heterodimer formation between all primer
    pairs. More negative values (shown in red) indicate stronger binding
    and higher risk of primer-dimer formation.

    The diagonal shows self-dimer potential for each primer.

    Risk levels:
    - Low (green): dG > -4 kcal/mol
    - Moderate (yellow/orange): -6 < dG <= -4 kcal/mol
    - High (red): dG <= -6 kcal/mol

    Args:
        primers: List of PrimerMetrics objects
        include_plotlyjs: How to include Plotly.js ('cdn', 'inline', False)
        height: Chart height in pixels
        max_primers: Maximum number of primers to display (default 20)
        show_values: If True, show dG values on cells with significant interactions

    Returns:
        HTML string with embedded Plotly chart, or empty string if
        Plotly is not available or primers is empty
    """
    if not _PLOTLY_AVAILABLE:
        logger.debug("Plotly not available, skipping dimer network heatmap")
        return ""

    if not primers or len(primers) < 2:
        return ""

    # Build the dimer matrix
    matrix, labels, sequences = _build_dimer_matrix(primers, max_primers)

    if matrix is None:
        return ""

    n = len(labels)

    # Create hover text matrix
    hover_text = []
    for i in range(n):
        row = []
        for j in range(n):
            dg = matrix[i][j]
            risk = "High" if dg < -6 else "Moderate" if dg < -4 else "Low"

            if i == j:
                text = (
                    f"<b>Self-dimer</b><br>"
                    f"Primer: {sequences[i]}<br>"
                    f"dG: {dg:.1f} kcal/mol<br>"
                    f"Risk: {risk}"
                )
            else:
                text = (
                    f"<b>Heterodimer</b><br>"
                    f"Primer 1: {sequences[i]}<br>"
                    f"Primer 2: {sequences[j]}<br>"
                    f"dG: {dg:.1f} kcal/mol<br>"
                    f"Risk: {risk}"
                )
            row.append(text)
        hover_text.append(row)

    # Custom colorscale: Green (safe) -> Yellow (moderate) -> Red (dangerous)
    # Values map to absolute dG, so higher abs value = more risk = redder
    colorscale = [
        [0.0, '#38a169'],    # Green - dG near 0 (no interaction)
        [0.25, '#68d391'],   # Light green
        [0.4, '#f6e05e'],    # Yellow - moderate risk starts
        [0.55, '#ed8936'],   # Orange
        [0.7, '#e53e3e'],    # Red - high risk
        [1.0, '#9b2c2c'],    # Dark red - very high risk
    ]

    # Normalize matrix for colorscale (map dG to 0-1 range)
    # 0 dG -> 0, -8 or worse -> 1
    z_normalized = []
    for row in matrix:
        z_row = []
        for val in row:
            # Normalize: 0 maps to 0, -8 maps to 1
            normalized = min(1.0, max(0.0, abs(val) / 8.0))
            z_row.append(normalized)
        z_normalized.append(z_row)

    fig = go.Figure()

    # Main heatmap
    fig.add_trace(go.Heatmap(
        z=z_normalized,
        x=labels,
        y=labels,
        colorscale=colorscale,
        zmin=0,
        zmax=1,
        text=hover_text,
        hoverinfo="text",
        showscale=True,
        colorbar=dict(
            title=dict(text="Interaction<br>Strength", font=dict(size=11)),
            tickvals=[0, 0.25, 0.5, 0.75, 1.0],
            ticktext=["None", "Low", "Moderate", "High", "Very High"],
            tickfont=dict(size=10),
            len=0.8,
        ),
    ))

    # Add text annotations for significant interactions
    annotations_list = []
    if show_values:
        for i in range(n):
            for j in range(n):
                dg = matrix[i][j]
                # Only show text for moderate/high risk interactions
                if dg < -3:
                    # Use white text on dark cells, dark text on light cells
                    color = 'white' if dg < -5 else '#2d3748'
                    annotations_list.append(dict(
                        x=labels[j],
                        y=labels[i],
                        text=f"{dg:.0f}",
                        showarrow=False,
                        font=dict(size=9, color=color),
                    ))

    # Count risk levels for subtitle
    high_risk = sum(1 for i in range(n) for j in range(i+1, n) if matrix[i][j] < -6)
    moderate_risk = sum(1 for i in range(n) for j in range(i+1, n) if -6 <= matrix[i][j] < -4)
    total_pairs = n * (n - 1) // 2

    subtitle_text = f"{total_pairs} primer pairs analyzed: {high_risk} high risk, {moderate_risk} moderate risk"

    fig.update_layout(
        title=dict(
            text="Primer-Dimer Interaction Network",
            font=dict(size=16, color="#2d3748"),
            subtitle=dict(
                text=subtitle_text,
                font=dict(size=11, color="#718096"),
            ),
        ),
        xaxis=dict(
            title="Primer",
            tickangle=45,
            side="bottom",
            tickfont=dict(size=10),
        ),
        yaxis=dict(
            title="Primer",
            autorange="reversed",  # Match matrix orientation
            tickfont=dict(size=10),
        ),
        height=height,
        width=height + 80,  # Approximately square with colorbar
        margin=dict(l=100, r=100, t=80, b=100),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(family="-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif"),
        annotations=annotations_list,
    )

    return fig.to_html(
        full_html=False,
        include_plotlyjs=include_plotlyjs,
        config={'displayModeBar': True, 'responsive': True}
    )


def render_dimer_network_graph(
    primers: List['PrimerMetrics'],
    include_plotlyjs: str = 'cdn',
    height: int = 500,
    max_primers: int = 15,
    dg_threshold: float = -4.0,
) -> str:
    """
    Render primer-dimer interactions as a network graph.

    Nodes represent primers, edges represent significant dimer interactions.
    Edge thickness and color indicate interaction strength. Node color
    indicates self-dimer risk.

    This is a complementary view to the heatmap, better for visualizing
    which primers have the most problematic interactions.

    Args:
        primers: List of PrimerMetrics objects
        include_plotlyjs: How to include Plotly.js ('cdn', 'inline', False)
        height: Chart height in pixels
        max_primers: Maximum number of primers to display
        dg_threshold: Only show edges with dG below this threshold

    Returns:
        HTML string with embedded Plotly chart, or empty string if
        Plotly is not available or primers is empty
    """
    if not _PLOTLY_AVAILABLE:
        logger.debug("Plotly not available, skipping dimer network graph")
        return ""

    if not primers or len(primers) < 2:
        return ""

    import math

    # Limit primers
    primers_subset = primers[:max_primers]
    n = len(primers_subset)

    if n < 2:
        return ""

    # Create labels and get sequences
    labels = []
    sequences = []
    for p in primers_subset:
        label = p.sequence[:8] + "..." if len(p.sequence) > 8 else p.sequence
        labels.append(label)
        sequences.append(p.sequence)

    # Calculate positions in a circle
    angles = [2 * math.pi * i / n for i in range(n)]
    node_x = [math.cos(a) for a in angles]
    node_y = [math.sin(a) for a in angles]

    # Calculate all pairwise dG values and filter edges
    edges = []
    for i in range(n):
        for j in range(i + 1, n):
            dg = _calculate_heterodimer_dg(sequences[i], sequences[j])
            if dg < dg_threshold:
                edges.append((i, j, dg))

    fig = go.Figure()

    # Add edges
    for i, j, dg in edges:
        # Normalize weight for line thickness (1-5 pixels)
        weight = min(abs(dg) / 8.0, 1.0) * 4 + 1

        # Color based on risk level
        if dg < -6:
            color = '#c53030'  # Red - high risk
        elif dg < -5:
            color = '#ed8936'  # Orange
        else:
            color = '#ecc94b'  # Yellow - moderate

        fig.add_trace(go.Scatter(
            x=[node_x[i], node_x[j], None],
            y=[node_y[i], node_y[j], None],
            mode='lines',
            line=dict(width=weight, color=color),
            hoverinfo='text',
            hovertext=f"<b>{labels[i]}</b> - <b>{labels[j]}</b><br>dG: {dg:.1f} kcal/mol",
            showlegend=False,
        ))

    # Count interactions per primer for node sizing
    interaction_counts = [0] * n
    for i, j, _ in edges:
        interaction_counts[i] += 1
        interaction_counts[j] += 1

    # Node sizes based on interaction count
    node_sizes = [max(25, 15 + count * 6) for count in interaction_counts]

    # Node colors based on self-dimer risk
    node_colors = []
    self_dimer_dgs = []
    for seq in sequences:
        self_dg = _calculate_self_dimer_dg(seq)
        self_dimer_dgs.append(self_dg)
        if self_dg < -6:
            node_colors.append('#c53030')  # Red
        elif self_dg < -4:
            node_colors.append('#ed8936')  # Orange
        else:
            node_colors.append('#38a169')  # Green

    # Create node hover text
    node_text = []
    for i in range(n):
        risk = "High" if self_dimer_dgs[i] < -6 else "Moderate" if self_dimer_dgs[i] < -4 else "Low"
        text = (
            f"<b>{labels[i]}</b><br>"
            f"Sequence: {sequences[i]}<br>"
            f"Self-dimer dG: {self_dimer_dgs[i]:.1f} kcal/mol ({risk})<br>"
            f"Interactions: {interaction_counts[i]}"
        )
        node_text.append(text)

    # Add nodes
    fig.add_trace(go.Scatter(
        x=node_x,
        y=node_y,
        mode='markers+text',
        marker=dict(
            size=node_sizes,
            color=node_colors,
            line=dict(width=2, color='white'),
        ),
        text=labels,
        textposition="top center",
        textfont=dict(size=9),
        hovertext=node_text,
        hoverinfo='text',
        showlegend=False,
    ))

    # Summary statistics for subtitle
    high_risk = sum(1 for _, _, dg in edges if dg < -6)
    moderate_risk = sum(1 for _, _, dg in edges if dg >= -6)
    subtitle = f"{len(edges)} interactions (dG < {dg_threshold}): {high_risk} high risk, {moderate_risk} moderate"

    fig.update_layout(
        title=dict(
            text="Primer-Dimer Network Graph",
            font=dict(size=16, color="#2d3748"),
            subtitle=dict(
                text=subtitle,
                font=dict(size=11, color="#718096"),
            ),
        ),
        xaxis=dict(
            showgrid=False,
            zeroline=False,
            showticklabels=False,
            range=[-1.5, 1.5],
        ),
        yaxis=dict(
            showgrid=False,
            zeroline=False,
            showticklabels=False,
            range=[-1.5, 1.5],
            scaleanchor="x",
            scaleratio=1,
        ),
        height=height,
        margin=dict(l=40, r=40, t=80, b=40),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(family="-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif"),
        hovermode='closest',
    )

    return fig.to_html(
        full_html=False,
        include_plotlyjs=include_plotlyjs,
        config={'displayModeBar': False, 'responsive': True}
    )
