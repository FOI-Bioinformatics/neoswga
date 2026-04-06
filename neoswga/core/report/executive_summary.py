"""
Executive summary report generation for SWGA primer sets.

Generates a one-page HTML report with:
- Quality grade (A-F)
- Key metrics dashboard
- Primer set table
- Go/no-go recommendation
"""

import logging
from dataclasses import dataclass, field
from datetime import datetime
from html import escape as html_escape
from pathlib import Path
from typing import List, Optional

from neoswga.core.report.metrics import (
    PipelineMetrics,
    PrimerMetrics,
    collect_pipeline_metrics,
)
from neoswga.core.report.quality import (
    QualityGrade,
    QualityAssessment,
    calculate_quality_grade,
    format_grade_display,
)
from neoswga.core.report.utils import (
    escape_format_braces,
    get_grade_colors,
    get_rating_class,
    get_progress_class,
    get_version as _get_version,
    GRADE_DESCRIPTIONS,
    ENRICHMENT_EXCELLENT_THRESHOLD,
)
from neoswga.core.report.visualizations import (
    is_plotly_available,
    render_component_radar,
)

logger = logging.getLogger(__name__)


@dataclass
class ExecutiveSummary:
    """Data for executive summary report."""
    metrics: PipelineMetrics
    quality: QualityAssessment
    generated_at: str
    version: str = field(default_factory=_get_version)


# HTML template with embedded CSS
EXECUTIVE_SUMMARY_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>NeoSWGA Primer Set Report</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}

        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto,
                         'Helvetica Neue', Arial, sans-serif;
            line-height: 1.5;
            color: #1a1a1a;
            background: #f8f9fa;
            padding: 20px;
        }}

        .container {{
            max-width: 900px;
            margin: 0 auto;
            background: white;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
            border-radius: 8px;
            overflow: hidden;
        }}

        /* Header */
        .header {{
            background: linear-gradient(135deg, #2c3e50 0%, #34495e 100%);
            color: white;
            padding: 24px 30px;
        }}

        .header h1 {{
            font-size: 1.75em;
            font-weight: 600;
            margin-bottom: 4px;
        }}

        .header .subtitle {{
            opacity: 0.85;
            font-size: 0.95em;
        }}

        .header .meta {{
            display: flex;
            gap: 20px;
            margin-top: 12px;
            font-size: 0.85em;
            opacity: 0.75;
        }}

        /* Grade Section */
        .grade-section {{
            display: flex;
            align-items: center;
            justify-content: center;
            padding: 30px;
            background: {grade_bg};
            border-bottom: 1px solid #e9ecef;
        }}

        .grade-box {{
            text-align: center;
        }}

        .grade-letter {{
            font-size: 4em;
            font-weight: 700;
            color: {grade_color};
            line-height: 1;
        }}

        .grade-label {{
            font-size: 1.1em;
            font-weight: 600;
            color: {grade_color};
            margin-top: 4px;
        }}

        .grade-score {{
            font-size: 0.9em;
            color: #666;
            margin-top: 8px;
        }}

        /* Metrics Grid */
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(4, 1fr);
            gap: 1px;
            background: #e9ecef;
            border-bottom: 1px solid #e9ecef;
        }}

        .metric-card {{
            background: white;
            padding: 20px;
            text-align: center;
        }}

        .metric-card .label {{
            font-size: 0.75em;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            color: #6c757d;
            margin-bottom: 8px;
        }}

        .metric-card .value {{
            font-size: 1.8em;
            font-weight: 600;
            color: #212529;
        }}

        .metric-card .rating {{
            font-size: 0.8em;
            margin-top: 6px;
            padding: 2px 8px;
            border-radius: 10px;
            display: inline-block;
        }}

        .rating-excellent {{ background: #d4edda; color: #155724; }}
        .rating-good {{ background: #cce5ff; color: #004085; }}
        .rating-acceptable {{ background: #fff3cd; color: #856404; }}
        .rating-poor {{ background: #f8d7da; color: #721c24; }}
        .rating-critical {{ background: #721c24; color: white; }}

        .coverage-badge {{
            font-size: 0.65em;
            padding: 1px 6px;
            border-radius: 8px;
            margin-left: 4px;
            vertical-align: middle;
            font-weight: 500;
        }}

        .badge-measured {{ background: #48bb78; color: white; }}
        .badge-estimated {{ background: #ecc94b; color: #744210; }}

        .gap-analysis {{
            padding: 16px 30px;
            background: #f8f9fa;
            border-bottom: 1px solid #e9ecef;
        }}

        .gap-analysis h3 {{
            font-size: 0.95em;
            font-weight: 600;
            color: #2c3e50;
            margin-bottom: 10px;
        }}

        .gap-metrics {{
            display: flex;
            gap: 24px;
            font-size: 0.9em;
            color: #495057;
        }}

        .gap-metric-item {{
            display: flex;
            flex-direction: column;
        }}

        .gap-metric-item .gap-label {{
            font-size: 0.8em;
            color: #6c757d;
            text-transform: uppercase;
            letter-spacing: 0.3px;
        }}

        .gap-metric-item .gap-value {{
            font-weight: 600;
            color: #212529;
        }}

        /* Progress bars */
        .progress-container {{
            margin-top: 8px;
            height: 6px;
            background: #e9ecef;
            border-radius: 3px;
            overflow: hidden;
        }}

        .progress-bar {{
            height: 100%;
            border-radius: 3px;
        }}

        .progress-excellent {{ background: #28a745; }}
        .progress-good {{ background: #17a2b8; }}
        .progress-acceptable {{ background: #ffc107; }}
        .progress-poor {{ background: #dc3545; }}

        /* Primer Table */
        .section {{
            padding: 24px 30px;
            border-bottom: 1px solid #e9ecef;
        }}

        .section h2 {{
            font-size: 1.1em;
            font-weight: 600;
            color: #2c3e50;
            margin-bottom: 16px;
        }}

        table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 0.9em;
        }}

        th {{
            background: #f8f9fa;
            padding: 10px 12px;
            text-align: left;
            font-weight: 600;
            color: #495057;
            border-bottom: 2px solid #dee2e6;
        }}

        td {{
            padding: 10px 12px;
            border-bottom: 1px solid #e9ecef;
        }}

        tr:hover {{
            background: #f8f9fa;
        }}

        .sequence {{
            font-family: 'SF Mono', 'Monaco', 'Inconsolata', monospace;
            font-size: 0.95em;
            letter-spacing: 0.5px;
        }}

        .quality-stars {{
            color: #ffc107;
            letter-spacing: 1px;
        }}

        /* Recommendation */
        .recommendation {{
            padding: 24px 30px;
            background: {rec_bg};
            border-left: 4px solid {rec_border};
        }}

        .recommendation h2 {{
            font-size: 1em;
            font-weight: 600;
            color: {rec_color};
            margin-bottom: 8px;
        }}

        .recommendation p {{
            color: #495057;
            font-size: 0.95em;
            line-height: 1.6;
        }}

        .considerations {{
            margin-top: 16px;
            padding-top: 16px;
            border-top: 1px solid rgba(0,0,0,0.1);
        }}

        .considerations h3 {{
            font-size: 0.85em;
            font-weight: 600;
            color: #495057;
            margin-bottom: 8px;
        }}

        .considerations ul {{
            margin: 0;
            padding-left: 20px;
        }}

        .considerations li {{
            font-size: 0.9em;
            color: #6c757d;
            margin-bottom: 4px;
        }}

        /* Footer */
        .footer {{
            padding: 16px 30px;
            background: #f8f9fa;
            text-align: center;
            font-size: 0.8em;
            color: #6c757d;
        }}

        /* Responsive */
        @media (max-width: 768px) {{
            .metrics-grid {{
                grid-template-columns: repeat(2, 1fr);
            }}

            .header .meta {{
                flex-direction: column;
                gap: 4px;
            }}
        }}

        @media print {{
            body {{
                background: white;
                padding: 0;
            }}

            .container {{
                box-shadow: none;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <!-- Header -->
        <div class="header">
            <h1>NeoSWGA Primer Set Report</h1>
            <div class="subtitle">{target_name} vs {background_name}</div>
            <div class="meta">
                <span>Generated: {timestamp}</span>
                <span>NeoSWGA v{version}</span>
            </div>
        </div>

        <!-- Quality Grade -->
        <div class="grade-section">
            <div class="grade-box">
                <div class="grade-letter">{grade_letter}</div>
                <div class="grade-label">{grade_description}</div>
                <div class="grade-score">Composite Score: {composite_score:.2f}/1.00</div>
            </div>
        </div>

        <!-- Key Metrics -->
        <div class="metrics-grid">
            <div class="metric-card">
                <div class="label">Coverage <span class="coverage-badge {coverage_badge_class}">{coverage_badge_label}</span></div>
                <div class="value">{coverage_pct:.1f}%</div>
                <div class="progress-container">
                    <div class="progress-bar {coverage_progress_class}"
                         style="width: {coverage_pct:.0f}%"></div>
                </div>
                <div class="rating {coverage_rating_class}">{coverage_rating}</div>
            </div>
            <div class="metric-card">
                <div class="label">Enrichment</div>
                <div class="value">{enrichment:.0f}x</div>
                <div class="progress-container">
                    <div class="progress-bar {enrichment_progress_class}"
                         style="width: {enrichment_bar_pct:.0f}%"></div>
                </div>
                <div class="rating {enrichment_rating_class}">{enrichment_rating}</div>
            </div>
            <div class="metric-card">
                <div class="label">Uniformity</div>
                <div class="value">{uniformity:.2f}</div>
                <div class="progress-container">
                    <div class="progress-bar {uniformity_progress_class}"
                         style="width: {uniformity_bar_pct:.0f}%"></div>
                </div>
                <div class="rating {uniformity_rating_class}">{uniformity_rating}</div>
            </div>
            <div class="metric-card">
                <div class="label">Dimer Risk</div>
                <div class="value">{dimer_risk_label}</div>
                <div class="progress-container">
                    <div class="progress-bar {dimer_progress_class}"
                         style="width: {dimer_bar_pct:.0f}%"></div>
                </div>
                <div class="rating {dimer_rating_class}">{dimer_rating}</div>
            </div>
        </div>

        {gap_analysis_html}

        <!-- Interactive Charts (if available) -->
        {interactive_charts}

        <!-- Primer Table -->
        <div class="section">
            <h2>Primer Set ({primer_count} primers)</h2>
            <table>
                <thead>
                    <tr>
                        <th>#</th>
                        <th>Sequence</th>
                        <th>Length</th>
                        <th>Tm (C)</th>
                        <th>GC%</th>
                        <th>Specificity</th>
                        <th>Quality</th>
                    </tr>
                </thead>
                <tbody>
                    {primer_rows}
                </tbody>
            </table>
        </div>

        <!-- Recommendation -->
        <div class="recommendation">
            <h2>{recommendation}</h2>
            <p>{recommendation_details}</p>
            {considerations_html}
        </div>

        <!-- Footer -->
        <div class="footer">
            NeoSWGA v{version} - Advanced SWGA Primer Design<br>
            Report generated: {timestamp}
        </div>
    </div>
</body>
</html>
"""


def _format_primer_row(idx: int, primer: PrimerMetrics) -> str:
    """Format a single primer as a table row."""
    # Calculate quality stars (1-5)
    quality_score = primer.amp_pred if primer.amp_pred > 0 else 0.5
    stars = min(5, max(1, int(quality_score * 5 + 0.5)))
    quality_stars = "★" * stars + "☆" * (5 - stars)

    # Format specificity
    if primer.specificity > 1000:
        spec_str = f"{primer.specificity/1000:.1f}k"
    else:
        spec_str = f"{primer.specificity:.0f}x"

    # Escape sequence to prevent XSS
    safe_sequence = html_escape(primer.sequence)

    return f"""<tr>
        <td>{idx}</td>
        <td class="sequence">{safe_sequence}</td>
        <td>{primer.length}</td>
        <td>{primer.tm:.1f}</td>
        <td>{primer.gc_content*100:.0f}</td>
        <td>{spec_str}</td>
        <td class="quality-stars">{quality_stars}</td>
    </tr>"""


def _format_considerations(considerations: List[str]) -> str:
    """Format considerations as HTML."""
    if not considerations:
        return ""

    # Escape each consideration to prevent XSS
    items = "\n".join(f"<li>{html_escape(c)}</li>" for c in considerations)
    return f"""<div class="considerations">
        <h3>Considerations:</h3>
        <ul>
            {items}
        </ul>
    </div>"""


def _get_component_by_name(
    components: list,
    name: str
) -> Optional[dict]:
    """Get a component by name."""
    for c in components:
        if c.name == name:
            return c
    return None


def create_executive_summary(
    metrics: PipelineMetrics,
    quality: QualityAssessment,
) -> ExecutiveSummary:
    """
    Create executive summary data structure.

    Args:
        metrics: Pipeline metrics
        quality: Quality assessment

    Returns:
        ExecutiveSummary ready for rendering
    """
    return ExecutiveSummary(
        metrics=metrics,
        quality=quality,
        generated_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    )


def render_executive_summary(summary: ExecutiveSummary, interactive: bool = False) -> str:
    """
    Render executive summary to HTML.

    Args:
        summary: ExecutiveSummary data
        interactive: If True, include interactive Plotly charts (requires plotly)

    Returns:
        HTML string
    """
    metrics = summary.metrics
    quality = summary.quality

    # Get color scheme
    colors = get_grade_colors(quality.grade)

    # Get component data
    coverage_comp = _get_component_by_name(quality.components, "Coverage")
    specificity_comp = _get_component_by_name(quality.components, "Specificity")
    uniformity_comp = _get_component_by_name(quality.components, "Uniformity")
    dimer_comp = _get_component_by_name(quality.components, "Dimer Risk")

    # Calculate display values
    coverage_pct = coverage_comp.raw_value * 100 if coverage_comp else 0
    enrichment = specificity_comp.raw_value if specificity_comp else 0
    uniformity = uniformity_comp.raw_value if uniformity_comp else 0.5
    dimer_risk = dimer_comp.raw_value if dimer_comp else 0

    # Progress bar percentages (capped at 100)
    enrichment_bar = min(100, (enrichment / ENRICHMENT_EXCELLENT_THRESHOLD) * 100)
    uniformity_bar = uniformity * 100
    dimer_bar = (1 - dimer_risk) * 100  # Invert for display

    # Dimer risk label
    if dimer_risk < 0.1:
        dimer_label = "Low"
    elif dimer_risk < 0.3:
        dimer_label = "Moderate"
    else:
        dimer_label = "High"

    # Format primer rows (limit to 20 for readability)
    MAX_PRIMERS_IN_TABLE = 20
    displayed_primers = metrics.primers[:MAX_PRIMERS_IN_TABLE]
    primer_rows = "\n".join(
        _format_primer_row(i + 1, p)
        for i, p in enumerate(displayed_primers)
    )

    # Add truncation notice if needed
    if len(metrics.primers) > MAX_PRIMERS_IN_TABLE:
        truncated_count = len(metrics.primers) - MAX_PRIMERS_IN_TABLE
        primer_rows += f"""<tr>
            <td colspan="7" style="text-align: center; font-style: italic; color: #6c757d;">
                ... and {truncated_count} more primer(s) not shown
            </td>
        </tr>"""

    # Target/background names (escaped to prevent XSS and format string injection)
    target_name = escape_format_braces(html_escape(
        metrics.target_genome.name if metrics.target_genome
        else "Target Genome"
    ))
    background_name = escape_format_braces(html_escape(
        metrics.background_genome.name if metrics.background_genome
        else "Background Genome"
    ))

    # Coverage source badge
    from_optimizer = (
        metrics.coverage is not None and metrics.coverage.from_optimizer
    )
    if from_optimizer:
        coverage_badge_label = "Measured"
        coverage_badge_class = "badge-measured"
    else:
        coverage_badge_label = "Estimated"
        coverage_badge_class = "badge-estimated"

    # Gap analysis section
    gap_analysis_html = ""
    if (metrics.coverage is not None and metrics.coverage.mean_gap > 0):
        mean_gap_kb = metrics.coverage.mean_gap / 1000
        max_gap_kb = metrics.coverage.max_gap / 1000
        gap_gini = metrics.coverage.gap_gini
        gap_analysis_html = f'''
        <div class="gap-analysis">
            <h3>Gap Analysis</h3>
            <div class="gap-metrics">
                <div class="gap-metric-item">
                    <span class="gap-label">Mean Gap</span>
                    <span class="gap-value">{mean_gap_kb:.1f} kb</span>
                </div>
                <div class="gap-metric-item">
                    <span class="gap-label">Max Gap</span>
                    <span class="gap-value">{max_gap_kb:.1f} kb</span>
                </div>
                <div class="gap-metric-item">
                    <span class="gap-label">Gap Gini</span>
                    <span class="gap-value">{gap_gini:.3f}</span>
                </div>
            </div>
        </div>'''

    # Generate interactive charts if requested and Plotly is available
    interactive_charts = ""
    if interactive and is_plotly_available():
        # Render component radar chart
        radar_html = render_component_radar(
            quality.components,
            include_plotlyjs='cdn',
            height=350,
        )
        if radar_html:
            interactive_charts = f'''
        <div class="section">
            <h2>Quality Analysis</h2>
            {radar_html}
        </div>
'''
    elif interactive:
        logger.debug("Interactive charts requested but Plotly not available")

    # Render template
    html = EXECUTIVE_SUMMARY_TEMPLATE.format(
        # Colors
        **colors,
        # Header
        target_name=target_name,
        background_name=background_name,
        timestamp=summary.generated_at,
        version=summary.version,
        # Grade
        grade_letter=quality.grade.value,
        grade_description=GRADE_DESCRIPTIONS[quality.grade],
        composite_score=quality.composite_score,
        # Coverage
        coverage_pct=coverage_pct,
        coverage_rating=coverage_comp.rating if coverage_comp else "N/A",
        coverage_rating_class=get_rating_class(
            coverage_comp.rating if coverage_comp else ""
        ),
        coverage_progress_class=get_progress_class(
            coverage_comp.rating if coverage_comp else ""
        ),
        coverage_badge_label=coverage_badge_label,
        coverage_badge_class=coverage_badge_class,
        # Gap analysis
        gap_analysis_html=gap_analysis_html,
        # Enrichment
        enrichment=enrichment,
        enrichment_bar_pct=enrichment_bar,
        enrichment_rating=specificity_comp.rating if specificity_comp else "N/A",
        enrichment_rating_class=get_rating_class(
            specificity_comp.rating if specificity_comp else ""
        ),
        enrichment_progress_class=get_progress_class(
            specificity_comp.rating if specificity_comp else ""
        ),
        # Uniformity
        uniformity=uniformity,
        uniformity_bar_pct=uniformity_bar,
        uniformity_rating=uniformity_comp.rating if uniformity_comp else "N/A",
        uniformity_rating_class=get_rating_class(
            uniformity_comp.rating if uniformity_comp else ""
        ),
        uniformity_progress_class=get_progress_class(
            uniformity_comp.rating if uniformity_comp else ""
        ),
        # Dimer
        dimer_risk_label=dimer_label,
        dimer_bar_pct=dimer_bar,
        dimer_rating=dimer_comp.rating if dimer_comp else "N/A",
        dimer_rating_class=get_rating_class(
            dimer_comp.rating if dimer_comp else ""
        ),
        dimer_progress_class=get_progress_class(
            dimer_comp.rating if dimer_comp else ""
        ),
        # Primers
        primer_count=metrics.primer_count,
        primer_rows=primer_rows,
        # Recommendation (escaped for XSS and format string injection)
        recommendation=escape_format_braces(html_escape(quality.recommendation)),
        recommendation_details=escape_format_braces(html_escape(quality.recommendation_details)),
        considerations_html=_format_considerations(quality.considerations),
        # Interactive charts
        interactive_charts=interactive_charts,
    )

    return html


def generate_executive_summary(
    results_dir: str,
    output_file: Optional[str] = None,
    interactive: bool = False,
) -> ExecutiveSummary:
    """
    Generate executive summary report from pipeline results.

    Args:
        results_dir: Path to results directory
        output_file: Output HTML file path (optional)
        interactive: If True, include interactive Plotly charts (requires plotly)

    Returns:
        ExecutiveSummary object

    Example:
        summary = generate_executive_summary('results/', 'report.html')
        print(f"Grade: {summary.quality.grade.value}")

        # With interactive charts
        summary = generate_executive_summary('results/', 'report.html', interactive=True)
    """
    logger.info(f"Generating executive summary for {results_dir}")

    # Collect metrics
    metrics = collect_pipeline_metrics(results_dir)

    # Calculate quality grade
    quality = calculate_quality_grade(metrics)

    # Create summary
    summary = create_executive_summary(metrics, quality)

    # Render and save if output file specified
    if output_file:
        html = render_executive_summary(summary, interactive=interactive)
        output_path = Path(output_file).resolve()

        try:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_text(html, encoding='utf-8')
            logger.info(f"Report saved to: {output_path}")
        except (PermissionError, OSError) as e:
            logger.error(f"Failed to write report to {output_path}: {e}")
            raise

    return summary
