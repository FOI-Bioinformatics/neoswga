"""
Technical report generation for SWGA primer sets.

Generates a comprehensive multi-page HTML report with:
- Pipeline execution summary
- Coverage analysis with visualizations
- Specificity analysis
- Thermodynamic properties
- Per-primer detailed analysis
- Recommendations
"""

import logging
from dataclasses import dataclass, field
from datetime import datetime
from html import escape as html_escape
from pathlib import Path
from typing import List, Optional, Dict, Any
import json

from neoswga.core.report.metrics import (
    PipelineMetrics,
    PrimerMetrics,
    FilteringStats,
    collect_pipeline_metrics,
)
from neoswga.core.report.quality import (
    QualityGrade,
    QualityAssessment,
    GradeComponent,
    calculate_quality_grade,
)
from neoswga.core.report.utils import (
    escape_format_braces,
    get_grade_colors,
    get_rating_class,
    get_progress_class,
    get_version as _get_version,
)
from neoswga.core.report.visualizations import (
    is_plotly_available,
    render_filtering_funnel,
    render_component_radar,
    render_tm_gc_distribution,
    render_coverage_specificity_scatter,
    render_primer_heatmap,
    render_dimer_network_heatmap,
    render_dimer_network_graph,
)

logger = logging.getLogger(__name__)


# Note: escape_format_braces is now imported from utils as escape_format_braces


@dataclass
class GapInfo:
    """Information about a coverage gap."""
    start: int
    end: int
    length: int
    chromosome: str = ""
    severity: str = "low"  # low, medium, high, critical
    description: str = ""


@dataclass
class PrimerProfile:
    """Detailed profile for a single primer."""
    sequence: str
    length: int
    gc_content: float
    tm: float
    delta_g: float
    hairpin_dg: float
    self_dimer_dg: float
    fg_sites: int
    bg_sites: int
    specificity: float
    gini: float
    strand_ratio: float
    three_prime_seq: str
    three_prime_stability: float
    unique_coverage: float  # % of genome only this primer covers
    contribution_score: float
    quality_rank: int


@dataclass
class InteractionPair:
    """Primer-primer interaction data."""
    primer1: str
    primer2: str
    delta_g: float
    risk_level: str  # low, moderate, high


@dataclass
class TechnicalReportData:
    """All data needed for the technical report."""
    # Metadata
    generated_at: str
    version: str = field(default_factory=_get_version)

    # Core data
    metrics: Optional[PipelineMetrics] = None
    quality: Optional[QualityAssessment] = None

    # Pipeline summary
    total_runtime: float = 0.0
    step_runtimes: Dict[str, float] = field(default_factory=dict)

    # Filtering funnel
    filtering_stages: List[tuple] = field(default_factory=list)

    # Coverage analysis
    coverage_by_region: Dict[str, float] = field(default_factory=dict)
    gaps: List[GapInfo] = field(default_factory=list)
    binding_distribution: Dict[str, int] = field(default_factory=dict)

    # Primer profiles
    primer_profiles: List[PrimerProfile] = field(default_factory=list)

    # Interactions
    interactions: List[InteractionPair] = field(default_factory=list)
    max_interaction_dg: float = 0.0


def _classify_gap_severity(length: int) -> str:
    """Classify gap severity by length."""
    if length >= 100000:
        return "critical"
    elif length >= 50000:
        return "high"
    elif length >= 20000:
        return "medium"
    else:
        return "low"


def _estimate_binding_dg(sequence: str, tm: float) -> float:
    """
    Estimate binding delta G from Tm using simplified thermodynamics.

    Uses the relationship: dG = dH - T*dS
    At Tm, dG = 0, so we can estimate dG at 37C.

    Args:
        sequence: Primer sequence
        tm: Melting temperature in Celsius

    Returns:
        Estimated delta G in kcal/mol (negative = favorable binding)
    """
    if tm <= 0 or not sequence:
        return 0.0

    # Simplified estimation based on Tm and sequence length
    # dG at 37C approximately: dG = -0.2 * (Tm - 37) - 0.4 * len(sequence)
    # This is a rough heuristic that gives reasonable values
    length = len(sequence)
    dg_estimate = -0.2 * (tm - 37.0) - 0.4 * length

    return max(-25.0, min(0.0, dg_estimate))  # Clamp to reasonable range


def _calculate_primer_profile(primer: PrimerMetrics, rank: int) -> PrimerProfile:
    """Create detailed primer profile."""
    # Get 3' end sequence (last 3 bases)
    three_prime = primer.sequence[-3:] if len(primer.sequence) >= 3 else primer.sequence

    # Estimate contribution score (weighted combination)
    # Clamp individual components to [0, 1] to prevent negative total
    specificity_score = min(primer.specificity / 100, 1.0)
    uniformity_score = max(0.0, 1.0 - primer.gini)  # Gini is 0-1
    amp_score = max(0.0, min(primer.amp_pred, 1.0))  # amp_pred normalized to 0-1 in metrics.py
    # Strand ratio contribution: 1.0 is ideal, penalize deviation but clamp to [0, 1]
    strand_deviation = min(abs(primer.strand_ratio - 1.0), 1.0)
    strand_score = 1.0 - strand_deviation

    contribution = (
        0.4 * specificity_score +
        0.3 * uniformity_score +
        0.2 * amp_score +
        0.1 * strand_score
    )

    # Estimate binding delta G from Tm (replaces hardcoded placeholder)
    delta_g = _estimate_binding_dg(primer.sequence, primer.tm)

    return PrimerProfile(
        sequence=primer.sequence,
        length=primer.length,
        gc_content=primer.gc_content,
        tm=primer.tm,
        delta_g=delta_g,
        hairpin_dg=primer.hairpin_dg,
        self_dimer_dg=primer.self_dimer_dg,
        fg_sites=primer.fg_sites,
        bg_sites=primer.bg_sites,
        specificity=primer.specificity,
        gini=primer.gini,
        strand_ratio=primer.strand_ratio,
        three_prime_seq=three_prime,
        three_prime_stability=primer.three_prime_stability,
        unique_coverage=0.0,  # Would need position data to calculate
        contribution_score=contribution,
        quality_rank=rank,
    )


def _estimate_interactions(primers: List[PrimerMetrics]) -> List[InteractionPair]:
    """
    Estimate primer-primer interactions.

    Uses 3' complementarity check to identify potential dimer formation.
    Optimized with early termination for large primer sets.

    Args:
        primers: List of primer metrics

    Returns:
        Top 10 worst interactions sorted by delta G
    """
    # Limit pair comparisons for very large primer sets
    # n*(n-1)/2 comparisons: 20 primers = 190, 50 primers = 1225, 100 primers = 4950
    MAX_PRIMERS = 50  # Beyond this, only analyze top primers
    if len(primers) > MAX_PRIMERS:
        # Sort by specificity and take top primers
        sorted_primers = sorted(primers, key=lambda p: p.specificity, reverse=True)
        primers = sorted_primers[:MAX_PRIMERS]
        logger.debug(f"Limited interaction analysis to top {MAX_PRIMERS} primers")

    interactions = []
    early_terminate = False

    # Pre-compute reverse complements for efficiency
    rc_table = str.maketrans('ATGC', 'TACG')

    for i, p1 in enumerate(primers):
        if early_terminate:
            break

        seq1 = p1.sequence
        if not seq1:
            continue

        for p2 in primers[i+1:]:
            seq2 = p2.sequence
            if not seq2:
                continue

            # Compute reverse complement of p2
            seq2_rc = seq2[::-1].translate(rc_table)

            # Check for 3' complementarity (most important for dimer formation)
            overlap = 0
            max_check = min(6, len(seq1), len(seq2_rc))

            for k in range(max_check):
                if seq1[-(k+1):] == seq2_rc[:k+1]:
                    overlap = k + 1

            # Only record significant interactions (4+ bp overlap)
            if overlap >= 4:
                dg = -1.5 * overlap  # Rough estimate: ~1.5 kcal/mol per bp

                if dg < -6:
                    risk = "high"
                elif dg < -4:
                    risk = "moderate"
                else:
                    risk = "low"

                interactions.append(InteractionPair(
                    primer1=seq1,
                    primer2=seq2,
                    delta_g=dg,
                    risk_level=risk,
                ))

                # Early termination if we have many high-risk interactions
                high_risk_count = sum(1 for x in interactions if x.risk_level == "high")
                if high_risk_count >= 5:
                    early_terminate = True
                    break

    # Sort by delta G (most negative = worst) and return top 10
    return sorted(interactions, key=lambda x: x.delta_g)[:10]


def collect_technical_report_data(results_dir: str) -> TechnicalReportData:
    """
    Collect all data needed for technical report.

    Args:
        results_dir: Path to results directory

    Returns:
        TechnicalReportData with all collected information
    """
    logger.info(f"Collecting technical report data from {results_dir}")

    # Collect base metrics
    metrics = collect_pipeline_metrics(results_dir)
    quality = calculate_quality_grade(metrics)

    # Initialize report data
    data = TechnicalReportData(
        generated_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        metrics=metrics,
        quality=quality,
    )

    # Build filtering funnel if available
    if metrics.filtering:
        data.filtering_stages = metrics.filtering.as_funnel()
    else:
        # Estimate from primer counts
        data.filtering_stages = [
            ("Initial candidates", metrics.primer_count * 100),  # Estimate
            ("After filtering", metrics.primer_count * 10),
            ("Final set", metrics.primer_count),
        ]

    # Create primer profiles
    data.primer_profiles = [
        _calculate_primer_profile(p, i + 1)
        for i, p in enumerate(metrics.primers)
    ]

    # Estimate interactions
    if metrics.primers:
        data.interactions = _estimate_interactions(metrics.primers)
        if data.interactions:
            data.max_interaction_dg = min(i.delta_g for i in data.interactions)

    return data


# HTML Template for Technical Report
TECHNICAL_REPORT_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>NeoSWGA Technical Report</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}

        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto,
                         'Helvetica Neue', Arial, sans-serif;
            line-height: 1.6;
            color: #1a1a1a;
            background: #f8f9fa;
            padding: 20px;
        }}

        .container {{
            max-width: 1000px;
            margin: 0 auto;
            background: white;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
            border-radius: 8px;
            overflow: hidden;
        }}

        /* Header */
        .header {{
            background: linear-gradient(135deg, #1a365d 0%, #2c5282 100%);
            color: white;
            padding: 30px 40px;
        }}

        .header h1 {{
            font-size: 2em;
            font-weight: 600;
            margin-bottom: 8px;
        }}

        .header .subtitle {{
            opacity: 0.9;
            font-size: 1.1em;
        }}

        .header .meta {{
            display: flex;
            gap: 24px;
            margin-top: 16px;
            font-size: 0.9em;
            opacity: 0.8;
        }}

        /* Table of Contents */
        .toc {{
            background: #f8f9fa;
            padding: 24px 40px;
            border-bottom: 1px solid #e9ecef;
        }}

        .toc h2 {{
            font-size: 1.1em;
            color: #495057;
            margin-bottom: 12px;
        }}

        .toc ul {{
            list-style: none;
            display: flex;
            flex-wrap: wrap;
            gap: 8px 24px;
        }}

        .toc a {{
            color: #2c5282;
            text-decoration: none;
            font-size: 0.9em;
        }}

        .toc a:hover {{
            text-decoration: underline;
        }}

        /* Sections */
        .section {{
            padding: 32px 40px;
            border-bottom: 1px solid #e9ecef;
        }}

        .section:last-child {{
            border-bottom: none;
        }}

        .section h2 {{
            font-size: 1.4em;
            color: #1a365d;
            margin-bottom: 20px;
            padding-bottom: 8px;
            border-bottom: 2px solid #e9ecef;
        }}

        .section h3 {{
            font-size: 1.1em;
            color: #2d3748;
            margin: 24px 0 12px 0;
        }}

        .section p {{
            color: #4a5568;
            margin-bottom: 16px;
        }}

        /* Grade Summary */
        .grade-summary {{
            display: flex;
            align-items: center;
            gap: 32px;
            padding: 24px;
            background: {grade_bg};
            border-radius: 8px;
            margin-bottom: 24px;
        }}

        .grade-circle {{
            width: 80px;
            height: 80px;
            border-radius: 50%;
            background: {grade_color};
            color: white;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 2.5em;
            font-weight: 700;
        }}

        .grade-details h3 {{
            margin: 0 0 4px 0;
            color: {grade_text};
        }}

        .grade-details p {{
            margin: 0;
            color: #4a5568;
        }}

        /* Info boxes */
        .info-box {{
            background: #f8f9fa;
            border-radius: 6px;
            padding: 16px 20px;
            margin: 16px 0;
        }}

        .info-box.parameters {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 12px;
        }}

        .param-item {{
            display: flex;
            justify-content: space-between;
        }}

        .param-label {{
            color: #718096;
            font-size: 0.9em;
        }}

        .param-value {{
            font-weight: 600;
            color: #2d3748;
        }}

        /* Pipeline Flow */
        .pipeline-flow {{
            display: flex;
            flex-direction: column;
            gap: 0;
            margin: 20px 0;
        }}

        .pipeline-step {{
            display: flex;
            align-items: center;
            padding: 12px 16px;
            background: #f8f9fa;
            border-left: 4px solid #2c5282;
            margin-left: 20px;
            position: relative;
        }}

        .pipeline-step::before {{
            content: '';
            position: absolute;
            left: -22px;
            top: 50%;
            transform: translateY(-50%);
            width: 12px;
            height: 12px;
            background: #2c5282;
            border-radius: 50%;
        }}

        .pipeline-step::after {{
            content: '';
            position: absolute;
            left: -16px;
            top: 50%;
            height: calc(100% + 24px);
            width: 2px;
            background: #cbd5e0;
        }}

        .pipeline-step:last-child::after {{
            display: none;
        }}

        .step-name {{
            flex: 1;
            font-weight: 500;
        }}

        .step-value {{
            color: #2c5282;
            font-weight: 600;
        }}

        .step-change {{
            font-size: 0.85em;
            color: #718096;
            margin-left: 8px;
        }}

        /* Filtering Funnel */
        .funnel {{
            margin: 24px 0;
        }}

        .funnel-stage {{
            display: flex;
            align-items: center;
            margin-bottom: 8px;
        }}

        .funnel-bar {{
            height: 32px;
            background: linear-gradient(90deg, #2c5282, #4299e1);
            border-radius: 4px;
            display: flex;
            align-items: center;
            padding: 0 12px;
            color: white;
            font-size: 0.85em;
            font-weight: 500;
            min-width: 60px;
        }}

        .funnel-label {{
            margin-left: 12px;
            color: #4a5568;
            font-size: 0.9em;
        }}

        .funnel-pct {{
            margin-left: auto;
            color: #718096;
            font-size: 0.85em;
        }}

        /* Tables */
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 16px 0;
            font-size: 0.9em;
        }}

        th {{
            background: #f8f9fa;
            padding: 12px;
            text-align: left;
            font-weight: 600;
            color: #2d3748;
            border-bottom: 2px solid #e2e8f0;
        }}

        td {{
            padding: 12px;
            border-bottom: 1px solid #e9ecef;
            color: #4a5568;
        }}

        tr:hover {{
            background: #f8f9fa;
        }}

        .sequence {{
            font-family: 'SF Mono', 'Monaco', 'Inconsolata', monospace;
            font-size: 0.95em;
            letter-spacing: 0.5px;
        }}

        /* Metrics Grid */
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 16px;
            margin: 20px 0;
        }}

        .metric-box {{
            background: #f8f9fa;
            padding: 16px;
            border-radius: 6px;
            text-align: center;
        }}

        .metric-box .value {{
            font-size: 1.8em;
            font-weight: 600;
            color: #2d3748;
        }}

        .metric-box .label {{
            font-size: 0.85em;
            color: #718096;
            margin-top: 4px;
        }}

        .metric-box .rating {{
            font-size: 0.8em;
            padding: 2px 8px;
            border-radius: 10px;
            display: inline-block;
            margin-top: 8px;
        }}

        .rating-excellent {{ background: #c6f6d5; color: #22543d; }}
        .rating-good {{ background: #bee3f8; color: #2a4365; }}
        .rating-acceptable {{ background: #fefcbf; color: #744210; }}
        .rating-poor {{ background: #fed7d7; color: #822727; }}

        /* Progress bars */
        .progress-container {{
            height: 8px;
            background: #e2e8f0;
            border-radius: 4px;
            overflow: hidden;
            margin-top: 8px;
        }}

        .progress-fill {{
            height: 100%;
            border-radius: 4px;
            transition: width 0.3s;
        }}

        .progress-excellent {{ background: #48bb78; }}
        .progress-good {{ background: #4299e1; }}
        .progress-acceptable {{ background: #ecc94b; }}
        .progress-poor {{ background: #f56565; }}

        /* Component breakdown */
        .component-row {{
            display: flex;
            align-items: center;
            padding: 12px 0;
            border-bottom: 1px solid #e9ecef;
        }}

        .component-name {{
            width: 140px;
            font-weight: 500;
            color: #2d3748;
        }}

        .component-bar-container {{
            flex: 1;
            margin: 0 16px;
        }}

        .component-bar {{
            height: 24px;
            background: #e2e8f0;
            border-radius: 4px;
            overflow: hidden;
        }}

        .component-bar-fill {{
            height: 100%;
            display: flex;
            align-items: center;
            padding: 0 8px;
            color: white;
            font-size: 0.8em;
            font-weight: 500;
        }}

        .component-score {{
            width: 60px;
            text-align: right;
            font-weight: 600;
        }}

        .component-weight {{
            width: 50px;
            text-align: right;
            color: #718096;
            font-size: 0.85em;
        }}

        /* Primer cards */
        .primer-card {{
            background: #f8f9fa;
            border-radius: 8px;
            padding: 20px;
            margin: 16px 0;
        }}

        .primer-card-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 16px;
        }}

        .primer-card-header h4 {{
            font-family: monospace;
            font-size: 1.1em;
            color: #2d3748;
        }}

        .primer-card-header .rank {{
            background: #2c5282;
            color: white;
            padding: 4px 12px;
            border-radius: 12px;
            font-size: 0.85em;
        }}

        .primer-card-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(140px, 1fr));
            gap: 12px;
        }}

        .primer-stat {{
            text-align: center;
        }}

        .primer-stat .value {{
            font-size: 1.2em;
            font-weight: 600;
            color: #2d3748;
        }}

        .primer-stat .label {{
            font-size: 0.8em;
            color: #718096;
        }}

        /* Interaction matrix */
        .interaction-warning {{
            background: #fffaf0;
            border-left: 4px solid #ed8936;
            padding: 12px 16px;
            margin: 16px 0;
            border-radius: 0 6px 6px 0;
        }}

        .interaction-warning.low {{
            background: #f0fff4;
            border-color: #48bb78;
        }}

        /* Recommendations */
        .recommendation-box {{
            padding: 20px;
            border-radius: 8px;
            margin: 16px 0;
        }}

        .recommendation-box.proceed {{
            background: #f0fff4;
            border: 1px solid #9ae6b4;
        }}

        .recommendation-box.caution {{
            background: #fffaf0;
            border: 1px solid #fbd38d;
        }}

        .recommendation-box.warning {{
            background: #fff5f5;
            border: 1px solid #feb2b2;
        }}

        .recommendation-box h4 {{
            margin-bottom: 8px;
            color: #2d3748;
        }}

        .recommendation-list {{
            margin: 12px 0 0 20px;
        }}

        .recommendation-list li {{
            margin-bottom: 8px;
            color: #4a5568;
        }}

        /* Footer */
        .footer {{
            background: #f8f9fa;
            padding: 20px 40px;
            text-align: center;
            font-size: 0.85em;
            color: #718096;
        }}

        /* Print styles */
        @media print {{
            body {{
                background: white;
                padding: 0;
            }}

            .container {{
                box-shadow: none;
            }}

            .section {{
                page-break-inside: avoid;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <!-- Header -->
        <div class="header">
            <h1>NeoSWGA Technical Report</h1>
            <div class="subtitle">{target_name} vs {background_name}</div>
            <div class="meta">
                <span>Generated: {timestamp}</span>
                <span>NeoSWGA v{version}</span>
                <span>Primers: {primer_count}</span>
            </div>
        </div>

        <!-- Table of Contents -->
        <div class="toc">
            <h2>Contents</h2>
            <ul>
                <li><a href="#summary">1. Executive Summary</a></li>
                <li><a href="#pipeline">2. Pipeline Execution</a></li>
                <li><a href="#coverage">3. Coverage Analysis</a></li>
                <li><a href="#specificity">4. Specificity Analysis</a></li>
                <li><a href="#thermodynamics">5. Thermodynamics</a></li>
                <li><a href="#primers">6. Primer Profiles</a></li>
                <li><a href="#interactions">7. Interactions</a></li>
                <li><a href="#recommendations">8. Recommendations</a></li>
            </ul>
        </div>

        <!-- Section 1: Executive Summary -->
        <div class="section" id="summary">
            <h2>1. Executive Summary</h2>

            <div class="grade-summary">
                <div class="grade-circle">{grade_letter}</div>
                <div class="grade-details">
                    <h3>{recommendation}</h3>
                    <p>{recommendation_details}</p>
                </div>
            </div>

            <div class="metrics-grid">
                {summary_metrics_html}
            </div>

            <h3>Quality Score Breakdown</h3>
            <div class="info-box">
                {components_html}
            </div>
        </div>

        <!-- Section 2: Pipeline Execution -->
        <div class="section" id="pipeline">
            <h2>2. Pipeline Execution Summary</h2>

            <h3>Input Files</h3>
            <div class="info-box parameters">
                <div class="param-item">
                    <span class="param-label">Target Genome</span>
                    <span class="param-value">{target_name}</span>
                </div>
                <div class="param-item">
                    <span class="param-label">Target Size</span>
                    <span class="param-value">{target_size}</span>
                </div>
                <div class="param-item">
                    <span class="param-label">Background Genome</span>
                    <span class="param-value">{background_name}</span>
                </div>
                <div class="param-item">
                    <span class="param-label">Background Size</span>
                    <span class="param-value">{background_size}</span>
                </div>
            </div>

            <h3>Parameters</h3>
            <div class="info-box parameters">
                {parameters_html}
            </div>

            <h3>Filtering Funnel</h3>
            <p>The pipeline progressively filters candidates through multiple stages:</p>
            <div class="funnel">
                {funnel_html}
            </div>
            {interactive_funnel}
        </div>

        <!-- Section 3: Coverage Analysis -->
        <div class="section" id="coverage">
            <h2>3. Coverage Analysis</h2>

            <p>{coverage_source_note}</p>

            <div class="metrics-grid">
                <div class="metric-box">
                    <div class="value">{coverage_pct:.1f}%</div>
                    <div class="label">Overall Coverage</div>
                    <div class="progress-container">
                        <div class="progress-fill {coverage_class}" style="width: {coverage_pct:.0f}%"></div>
                    </div>
                </div>
                <div class="metric-box">
                    <div class="value">{total_sites:,}</div>
                    <div class="label">Binding Sites</div>
                </div>
                <div class="metric-box">
                    <div class="value">{sites_per_mb:.0f}</div>
                    <div class="label">Sites per Mbp</div>
                </div>
                <div class="metric-box">
                    <div class="value">{mean_spacing}</div>
                    <div class="label">Mean Spacing</div>
                </div>
            </div>

            <h3>Binding Uniformity</h3>
            <p>The Gini coefficient measures how evenly primers bind across the genome.
            A value of 0 indicates perfectly uniform binding, while 1 indicates all binding
            concentrated in one location.</p>

            <div class="info-box">
                <div class="param-item">
                    <span class="param-label">Gini Coefficient</span>
                    <span class="param-value">{gini:.3f}</span>
                </div>
                <div class="param-item">
                    <span class="param-label">Uniformity Score</span>
                    <span class="param-value">{uniformity_score:.1%}</span>
                </div>
                <div class="param-item">
                    <span class="param-label">Assessment</span>
                    <span class="param-value">{uniformity_assessment}</span>
                </div>
            </div>

            {gaps_html}

            {interactive_coverage}
        </div>

        <!-- Section 4: Specificity Analysis -->
        <div class="section" id="specificity">
            <h2>4. Specificity Analysis</h2>

            <p>Specificity measures how selectively primers bind to target DNA versus
            background DNA.</p>

            <div class="metrics-grid">
                <div class="metric-box">
                    <div class="value">{enrichment:.0f}x</div>
                    <div class="label">Enrichment Ratio</div>
                    <div class="rating {enrichment_class}">{enrichment_rating}</div>
                </div>
                <div class="metric-box">
                    <div class="value">{target_density:.0f}</div>
                    <div class="label">Target Sites/Mbp</div>
                </div>
                <div class="metric-box">
                    <div class="value">{bg_density:.1f}</div>
                    <div class="label">Background Sites/Mbp</div>
                </div>
                <div class="metric-box">
                    <div class="value">{bg_sites:,}</div>
                    <div class="label">Background Sites</div>
                </div>
            </div>

            <h3>Interpretation</h3>
            <p>{specificity_interpretation}</p>

            <h3>Per-Primer Specificity</h3>
            <table>
                <thead>
                    <tr>
                        <th>Primer</th>
                        <th>Target Sites</th>
                        <th>Background Sites</th>
                        <th>Specificity</th>
                    </tr>
                </thead>
                <tbody>
                    {specificity_table_html}
                </tbody>
            </table>
        </div>

        <!-- Section 5: Thermodynamics -->
        <div class="section" id="thermodynamics">
            <h2>5. Thermodynamic Analysis</h2>

            <h3>Melting Temperature Profile</h3>
            <div class="info-box parameters">
                <div class="param-item">
                    <span class="param-label">Reaction Temperature</span>
                    <span class="param-value">{reaction_temp}C</span>
                </div>
                <div class="param-item">
                    <span class="param-label">Polymerase</span>
                    <span class="param-value">{polymerase}</span>
                </div>
                <div class="param-item">
                    <span class="param-label">Mean Tm</span>
                    <span class="param-value">{mean_tm:.1f}C</span>
                </div>
                <div class="param-item">
                    <span class="param-label">Tm Range</span>
                    <span class="param-value">{min_tm:.1f} - {max_tm:.1f}C</span>
                </div>
            </div>

            <h3>Secondary Structure Risk</h3>
            <p>Analysis of hairpin and self-dimer formation potential at reaction temperature.</p>

            <div class="info-box">
                <div class="param-item">
                    <span class="param-label">Max Self-Dimer dG</span>
                    <span class="param-value">{max_self_dimer:.1f} kcal/mol</span>
                </div>
                <div class="param-item">
                    <span class="param-label">Assessment</span>
                    <span class="param-value">{secondary_structure_assessment}</span>
                </div>
            </div>

            <h3>Primer Tm Distribution</h3>
            {interactive_tm_gc}
            <table>
                <thead>
                    <tr>
                        <th>Primer</th>
                        <th>Length</th>
                        <th>GC%</th>
                        <th>Tm (C)</th>
                        <th>dTm from Reaction</th>
                    </tr>
                </thead>
                <tbody>
                    {tm_table_html}
                </tbody>
            </table>
        </div>

        <!-- Section 6: Primer Profiles -->
        <div class="section" id="primers">
            <h2>6. Detailed Primer Profiles</h2>

            <p>Comprehensive analysis of each primer in the set.</p>

            {interactive_heatmap}

            {primer_profiles_html}
        </div>

        <!-- Section 7: Interactions -->
        <div class="section" id="interactions">
            <h2>7. Primer-Primer Interactions</h2>

            <p>Analysis of potential primer-dimer formation between primers in the set.</p>

            {interaction_assessment_html}

            {interactive_dimer_heatmap}

            <h3>Top Interactions</h3>
            <table>
                <thead>
                    <tr>
                        <th>Primer 1</th>
                        <th>Primer 2</th>
                        <th>dG (kcal/mol)</th>
                        <th>Risk Level</th>
                    </tr>
                </thead>
                <tbody>
                    {interactions_table_html}
                </tbody>
            </table>
        </div>

        <!-- Section 8: Recommendations -->
        <div class="section" id="recommendations">
            <h2>8. Recommendations</h2>

            {final_recommendation_html}

            {considerations_html}

            <h3>Suggested Experimental Protocol</h3>
            <div class="info-box">
                {protocol_html}
            </div>
        </div>

        <!-- Footer -->
        <div class="footer">
            <p><strong>NeoSWGA Technical Report</strong></p>
            <p>Generated with NeoSWGA v{version}</p>
            <p>{timestamp}</p>
        </div>
    </div>
</body>
</html>
"""


# Note: get_grade_colors is now imported from utils


def _format_size(size: int) -> str:
    """Format genome size for display."""
    if size >= 1_000_000_000:
        return f"{size / 1_000_000_000:.2f} Gbp"
    elif size >= 1_000_000:
        return f"{size / 1_000_000:.2f} Mbp"
    elif size >= 1_000:
        return f"{size / 1_000:.1f} kbp"
    else:
        return f"{size} bp"


def _format_spacing(sites: int, genome_size: int) -> str:
    """Format mean spacing for display."""
    if sites <= 0 or genome_size <= 0:
        return "N/A"
    spacing = genome_size / sites
    if spacing >= 1000:
        return f"{spacing / 1000:.1f} kb"
    else:
        return f"{spacing:.0f} bp"


# Note: get_rating_class and get_progress_class are now imported from utils


def _render_summary_metrics(quality: QualityAssessment) -> str:
    """Render summary metrics grid."""
    html = ""
    for comp in quality.components:
        rating_class = get_rating_class(comp.rating)
        progress_class = get_progress_class(comp.rating)

        # Format value based on component
        if comp.name == "Coverage":
            value_str = f"{comp.raw_value:.1%}"
        elif comp.name == "Specificity":
            value_str = f"{comp.raw_value:.0f}x"
        elif comp.name == "Uniformity":
            value_str = f"{comp.raw_value:.2f}"
        elif comp.name == "Thermodynamics":
            value_str = f"{comp.raw_value:.1f}C"
        else:
            value_str = f"{comp.raw_value:.2f}"

        html += f"""
        <div class="metric-box">
            <div class="value">{value_str}</div>
            <div class="label">{comp.name}</div>
            <div class="progress-container">
                <div class="progress-fill {progress_class}" style="width: {comp.normalized_score * 100:.0f}%"></div>
            </div>
            <div class="rating {rating_class}">{comp.rating}</div>
        </div>
        """
    return html


def _render_components_breakdown(quality: QualityAssessment) -> str:
    """Render component score breakdown."""
    html = ""
    for comp in quality.components:
        progress_class = get_progress_class(comp.rating)
        bar_color = {
            "progress-excellent": "#48bb78",
            "progress-good": "#4299e1",
            "progress-acceptable": "#ecc94b",
            "progress-poor": "#f56565",
        }.get(progress_class, "#a0aec0")

        html += f"""
        <div class="component-row">
            <div class="component-name">{comp.name}</div>
            <div class="component-bar-container">
                <div class="component-bar">
                    <div class="component-bar-fill" style="width: {comp.normalized_score * 100:.0f}%; background: {bar_color};">
                        {comp.rating}
                    </div>
                </div>
            </div>
            <div class="component-score">{comp.normalized_score:.2f}</div>
            <div class="component-weight">{comp.weight:.0%}</div>
        </div>
        """
    return html


def _render_parameters(params: Dict) -> str:
    """Render parameters grid."""
    key_params = [
        ('min_k', 'Min K-mer'),
        ('max_k', 'Max K-mer'),
        ('min_fg_freq', 'Min FG Freq'),
        ('max_bg_freq', 'Max BG Freq'),
        ('max_gini', 'Max Gini'),
        ('min_tm', 'Min Tm'),
        ('max_tm', 'Max Tm'),
        ('polymerase', 'Polymerase'),
        ('reaction_temp', 'Reaction Temp'),
        ('optimization_method', 'Optimization'),
    ]

    html = ""
    for key, label in key_params:
        value = params.get(key, 'N/A')
        if isinstance(value, float):
            if value < 0.01:
                value_str = f"{value:.2e}"
            else:
                value_str = f"{value:.2f}"
        else:
            # Escape string values to prevent XSS and format string injection
            value_str = escape_format_braces(html_escape(str(value)))
        html += f"""
        <div class="param-item">
            <span class="param-label">{label}</span>
            <span class="param-value">{value_str}</span>
        </div>
        """
    return html


def _render_funnel(stages: List[tuple]) -> str:
    """Render filtering funnel visualization."""
    if not stages:
        return "<p>Filtering statistics not available.</p>"

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


def _render_specificity_table(primers: List[PrimerMetrics]) -> str:
    """Render per-primer specificity table."""
    html = ""
    for p in primers[:10]:  # Limit to top 10
        spec_str = f"{p.specificity:.0f}x" if p.specificity < 10000 else f"{p.specificity/1000:.1f}kx"
        safe_sequence = html_escape(p.sequence)
        html += f"""
        <tr>
            <td class="sequence">{safe_sequence}</td>
            <td>{p.fg_sites:,}</td>
            <td>{p.bg_sites:,}</td>
            <td>{spec_str}</td>
        </tr>
        """
    return html


def _render_tm_table(primers: List[PrimerMetrics], reaction_temp: float) -> str:
    """Render Tm distribution table."""
    html = ""
    for p in primers:
        delta_tm = reaction_temp - p.tm
        delta_str = f"+{delta_tm:.1f}" if delta_tm >= 0 else f"{delta_tm:.1f}"
        safe_sequence = html_escape(p.sequence)
        html += f"""
        <tr>
            <td class="sequence">{safe_sequence}</td>
            <td>{p.length}</td>
            <td>{p.gc_content * 100:.0f}%</td>
            <td>{p.tm:.1f}</td>
            <td>{delta_str}C</td>
        </tr>
        """
    return html


def _render_primer_profiles(profiles: List[PrimerProfile]) -> str:
    """Render detailed primer profile cards."""
    html = ""
    for profile in profiles[:10]:  # Limit to top 10
        safe_sequence = html_escape(profile.sequence)
        html += f"""
        <div class="primer-card">
            <div class="primer-card-header">
                <h4>{safe_sequence}</h4>
                <span class="rank">#{profile.quality_rank}</span>
            </div>
            <div class="primer-card-grid">
                <div class="primer-stat">
                    <div class="value">{profile.length}</div>
                    <div class="label">Length (bp)</div>
                </div>
                <div class="primer-stat">
                    <div class="value">{profile.gc_content * 100:.0f}%</div>
                    <div class="label">GC Content</div>
                </div>
                <div class="primer-stat">
                    <div class="value">{profile.tm:.1f}C</div>
                    <div class="label">Tm</div>
                </div>
                <div class="primer-stat">
                    <div class="value">{profile.specificity:.0f}x</div>
                    <div class="label">Specificity</div>
                </div>
                <div class="primer-stat">
                    <div class="value">{profile.gini:.2f}</div>
                    <div class="label">Gini</div>
                </div>
                <div class="primer-stat">
                    <div class="value">{profile.contribution_score:.2f}</div>
                    <div class="label">Contribution</div>
                </div>
            </div>
        </div>
        """
    return html


def _render_interactions_table(interactions: List[InteractionPair]) -> str:
    """Render interactions table."""
    if not interactions:
        return "<tr><td colspan='4'>No significant interactions detected.</td></tr>"

    html = ""
    for inter in interactions:
        risk_class = {
            "high": "color: #c53030;",
            "moderate": "color: #b7791f;",
            "low": "color: #276749;",
        }.get(inter.risk_level, "")

        safe_primer1 = html_escape(inter.primer1)
        safe_primer2 = html_escape(inter.primer2)

        html += f"""
        <tr>
            <td class="sequence">{safe_primer1}</td>
            <td class="sequence">{safe_primer2}</td>
            <td>{inter.delta_g:.1f}</td>
            <td style="{risk_class} font-weight: 600;">{inter.risk_level.upper()}</td>
        </tr>
        """
    return html


def _render_interaction_assessment(data: TechnicalReportData) -> str:
    """Render interaction risk assessment."""
    if not data.interactions:
        return """
        <div class="interaction-warning low">
            <strong>Low Risk:</strong> No significant primer-primer interactions detected.
            The primer set should amplify efficiently without dimer competition.
        </div>
        """

    high_risk = sum(1 for i in data.interactions if i.risk_level == "high")

    if high_risk > 0:
        return f"""
        <div class="interaction-warning">
            <strong>Warning:</strong> {high_risk} high-risk interaction(s) detected.
            These may compete with target amplification at lower DNA concentrations.
            Consider monitoring amplification efficiency.
        </div>
        """
    else:
        return """
        <div class="interaction-warning low">
            <strong>Acceptable:</strong> Only low to moderate interactions detected.
            These are unlikely to significantly impact amplification.
        </div>
        """


def _render_final_recommendation(quality: QualityAssessment) -> str:
    """Render final recommendation box."""
    box_class = {
        QualityGrade.A: "proceed",
        QualityGrade.B: "proceed",
        QualityGrade.C: "caution",
        QualityGrade.D: "warning",
        QualityGrade.F: "warning",
    }.get(quality.grade, "caution")

    # Escape recommendation text for XSS and format string injection
    safe_rec = escape_format_braces(html_escape(quality.recommendation))
    safe_details = escape_format_braces(html_escape(quality.recommendation_details))

    return f"""
    <div class="recommendation-box {box_class}">
        <h4>{safe_rec}</h4>
        <p>{safe_details}</p>
    </div>
    """


def _render_considerations(considerations: List[str]) -> str:
    """Render considerations list."""
    if not considerations:
        return ""

    # Escape each consideration to prevent XSS
    items = "\n".join(f"<li>{html_escape(c)}</li>" for c in considerations)
    return f"""
    <h3>Considerations</h3>
    <ul class="recommendation-list">
        {items}
    </ul>
    """


def _render_protocol(params: Dict, quality: QualityAssessment) -> str:
    """Render suggested protocol."""
    polymerase = escape_format_braces(html_escape(str(params.get('polymerase', 'phi29'))))
    temp = params.get('reaction_temp', 30)

    return f"""
    <div class="param-item">
        <span class="param-label">Primer Concentration</span>
        <span class="param-value">2.5 uM each</span>
    </div>
    <div class="param-item">
        <span class="param-label">Reaction Temperature</span>
        <span class="param-value">{temp}C</span>
    </div>
    <div class="param-item">
        <span class="param-label">Reaction Time</span>
        <span class="param-value">8-16 hours</span>
    </div>
    <div class="param-item">
        <span class="param-label">Polymerase</span>
        <span class="param-value">{polymerase}</span>
    </div>
    <div class="param-item">
        <span class="param-label">dNTP Concentration</span>
        <span class="param-value">1 mM each</span>
    </div>
    """


def render_technical_report(data: TechnicalReportData, interactive: bool = False) -> str:
    """
    Render technical report to HTML.

    Args:
        data: TechnicalReportData with all collected information
        interactive: If True, include interactive Plotly charts (requires plotly)

    Returns:
        HTML string
    """
    metrics = data.metrics
    quality = data.quality

    # Get color scheme
    colors = get_grade_colors(quality.grade)

    # Calculate derived values
    coverage_pct = metrics.coverage.overall_coverage * 100 if metrics.coverage else 0
    if metrics.coverage and metrics.coverage.from_optimizer:
        coverage_source_note = (
            "Coverage source: Measured from primer binding positions."
        )
    else:
        coverage_source_note = (
            "Coverage source: Estimated from primer count (~30kb per primer)."
        )
    total_sites = sum(p.fg_sites for p in metrics.primers)
    # Use 0 as default to indicate missing data, consistent with metrics.py
    fg_size = metrics.parameters.get('fg_size', metrics.parameters.get('foreground_size', 0))
    sites_per_mb = (total_sites / fg_size) * 1_000_000 if fg_size > 0 else 0

    # Get component for coverage rating
    coverage_comp = next((c for c in quality.components if c.name == "Coverage"), None)
    coverage_class = get_progress_class(coverage_comp.rating if coverage_comp else "")

    # Uniformity
    gini = metrics.uniformity.max_gini if metrics.uniformity else 0
    uniformity_score = 1 - gini
    if uniformity_score >= 0.85:
        uniformity_assessment = "Excellent - highly uniform binding"
    elif uniformity_score >= 0.7:
        uniformity_assessment = "Good - mostly uniform binding"
    elif uniformity_score >= 0.55:
        uniformity_assessment = "Acceptable - some clustering"
    else:
        uniformity_assessment = "Poor - significant clustering"

    # Specificity
    enrichment = metrics.specificity.enrichment_ratio if metrics.specificity else 0
    target_density = metrics.specificity.target_density if metrics.specificity else 0
    bg_density = metrics.specificity.background_density if metrics.specificity else 0
    bg_sites = sum(p.bg_sites for p in metrics.primers)

    spec_comp = next((c for c in quality.components if c.name == "Specificity"), None)
    enrichment_class = get_rating_class(spec_comp.rating if spec_comp else "")
    enrichment_rating = spec_comp.rating if spec_comp else "N/A"

    if enrichment >= 500:
        specificity_interpretation = (
            f"The enrichment ratio of {enrichment:.0f}x indicates excellent specificity. "
            f"For every background amplicon, approximately {enrichment:.0f} target amplicons "
            "are expected. This is suitable for most SWGA applications."
        )
    elif enrichment >= 100:
        specificity_interpretation = (
            f"The enrichment ratio of {enrichment:.0f}x indicates good specificity. "
            "This should provide adequate target enrichment for standard applications."
        )
    else:
        specificity_interpretation = (
            f"The enrichment ratio of {enrichment:.0f}x is moderate. "
            "Consider using background-aware optimization for improved specificity."
        )

    # Thermodynamics
    thermo = metrics.thermodynamics
    reaction_temp = thermo.reaction_temp if thermo else 30
    polymerase = thermo.polymerase if thermo else "phi29"
    mean_tm = thermo.mean_tm if thermo else 0
    min_tm = thermo.min_tm if thermo else 0
    max_tm = thermo.max_tm if thermo else 0

    # Use min() to get the worst (most negative) dimer score
    worst_self_dimer = min((p.self_dimer_dg for p in metrics.primers), default=0)
    if worst_self_dimer > -3:
        secondary_structure_assessment = "Low risk - no significant secondary structures"
    elif worst_self_dimer > -5:
        secondary_structure_assessment = "Moderate risk - some self-complementarity"
    else:
        secondary_structure_assessment = "High risk - significant secondary structures possible"

    # Gaps - show real data if available from optimizer, otherwise placeholder
    if metrics.coverage and metrics.coverage.from_optimizer and metrics.coverage.mean_gap > 0:
        mean_gap_kb = metrics.coverage.mean_gap / 1000.0
        max_gap_kb = metrics.coverage.max_gap / 1000.0
        gaps_html = f"""
    <h3>Gap Analysis</h3>
    <p>Gap statistics measured from primer binding positions across the target genome.</p>
    <div class="info-box">
        <div class="param-item">
            <span class="param-label">Mean Gap</span>
            <span class="param-value">{mean_gap_kb:.1f} kb</span>
        </div>
        <div class="param-item">
            <span class="param-label">Max Gap</span>
            <span class="param-value">{max_gap_kb:.1f} kb</span>
        </div>
        <div class="param-item">
            <span class="param-label">Gap Gini</span>
            <span class="param-value">{metrics.coverage.gap_gini:.2f} (lower is more uniform)</span>
        </div>
        <div class="param-item">
            <span class="param-label">Gap Entropy</span>
            <span class="param-value">{metrics.coverage.gap_entropy:.2f} bits</span>
        </div>
    </div>
    """
    else:
        gaps_html = """
    <h3>Coverage Gaps</h3>
    <p>Gap analysis requires primer position data. Run optimization to generate gap metrics.</p>
    """

    # Escape target/background names to prevent XSS and format string injection
    target_name = escape_format_braces(html_escape(
        metrics.target_genome.name if metrics.target_genome else "Target"
    ))
    background_name = escape_format_braces(html_escape(
        metrics.background_genome.name if metrics.background_genome else "Background"
    ))

    # Generate interactive charts if requested and Plotly is available
    interactive_funnel = ""
    interactive_coverage = ""
    interactive_tm_gc = ""
    interactive_heatmap = ""
    interactive_dimer_heatmap = ""

    if interactive and is_plotly_available():
        # Use 'cdn' for first chart, False for subsequent to avoid loading
        # Plotly.js multiple times
        include_js = 'cdn'

        # Filtering funnel chart
        funnel_chart = render_filtering_funnel(
            data.filtering_stages,
            include_plotlyjs=include_js,
            height=350,
        )
        if funnel_chart:
            interactive_funnel = f'<div class="interactive-chart">{funnel_chart}</div>'
            include_js = False  # Only include Plotly.js once

        # Coverage vs specificity scatter
        coverage_chart = render_coverage_specificity_scatter(
            metrics.primers,
            genome_size=fg_size,
            include_plotlyjs=include_js,
            height=400,
        )
        if coverage_chart:
            interactive_coverage = f'''
            <h3>Coverage vs Specificity</h3>
            <div class="interactive-chart">{coverage_chart}</div>
'''
            include_js = False

        # Tm/GC distribution
        tm_gc_chart = render_tm_gc_distribution(
            metrics.primers,
            reaction_temp=reaction_temp,
            include_plotlyjs=include_js,
            height=350,
        )
        if tm_gc_chart:
            interactive_tm_gc = f'<div class="interactive-chart">{tm_gc_chart}</div>'
            include_js = False

        # Primer heatmap
        heatmap_chart = render_primer_heatmap(
            metrics.primers,
            include_plotlyjs=include_js,
            height=400,
        )
        if heatmap_chart:
            interactive_heatmap = f'''
            <h3>Primer Metrics Overview</h3>
            <div class="interactive-chart">{heatmap_chart}</div>
'''
            include_js = False

        # Dimer network heatmap
        dimer_heatmap_chart = render_dimer_network_heatmap(
            metrics.primers,
            include_plotlyjs=include_js,
            height=500,
            max_primers=20,
            show_values=True,
        )
        if dimer_heatmap_chart:
            interactive_dimer_heatmap = f'''
            <h3>Interaction Heatmap</h3>
            <p>Estimated delta G (kcal/mol) for heterodimer formation. More negative values
            indicate stronger binding and higher dimer formation risk.</p>
            <div class="interactive-chart">{dimer_heatmap_chart}</div>
'''
    elif interactive:
        logger.debug("Interactive charts requested but Plotly not available")

    # Render all sections
    html = TECHNICAL_REPORT_TEMPLATE.format(
        # Colors
        **colors,
        # Header
        target_name=target_name,
        background_name=background_name,
        timestamp=data.generated_at,
        version=data.version,
        primer_count=metrics.primer_count,
        # Summary (recommendation escaped for XSS and format string injection)
        grade_letter=quality.grade.value,
        recommendation=escape_format_braces(html_escape(quality.recommendation)),
        recommendation_details=escape_format_braces(html_escape(quality.recommendation_details)),
        summary_metrics_html=_render_summary_metrics(quality),
        components_html=_render_components_breakdown(quality),
        # Pipeline
        target_size=_format_size(metrics.target_genome.size) if metrics.target_genome else "N/A",
        background_size=_format_size(metrics.background_genome.size) if metrics.background_genome else "N/A",
        parameters_html=_render_parameters(metrics.parameters),
        funnel_html=_render_funnel(data.filtering_stages),
        # Coverage
        coverage_pct=coverage_pct,
        coverage_class=coverage_class,
        coverage_source_note=coverage_source_note,
        total_sites=total_sites,
        sites_per_mb=sites_per_mb,
        mean_spacing=_format_spacing(total_sites, fg_size),
        gini=gini,
        uniformity_score=uniformity_score,
        uniformity_assessment=uniformity_assessment,
        gaps_html=gaps_html,
        # Specificity
        enrichment=enrichment,
        enrichment_class=enrichment_class,
        enrichment_rating=enrichment_rating,
        target_density=target_density,
        bg_density=bg_density,
        bg_sites=bg_sites,
        specificity_interpretation=specificity_interpretation,
        specificity_table_html=_render_specificity_table(metrics.primers),
        # Thermodynamics
        reaction_temp=reaction_temp,
        polymerase=polymerase,
        mean_tm=mean_tm,
        min_tm=min_tm,
        max_tm=max_tm,
        max_self_dimer=worst_self_dimer,
        secondary_structure_assessment=secondary_structure_assessment,
        tm_table_html=_render_tm_table(metrics.primers, reaction_temp),
        # Primers
        primer_profiles_html=_render_primer_profiles(data.primer_profiles),
        # Interactions
        interaction_assessment_html=_render_interaction_assessment(data),
        interactions_table_html=_render_interactions_table(data.interactions),
        # Recommendations
        final_recommendation_html=_render_final_recommendation(quality),
        considerations_html=_render_considerations(quality.considerations),
        protocol_html=_render_protocol(metrics.parameters, quality),
        # Interactive charts
        interactive_funnel=interactive_funnel,
        interactive_coverage=interactive_coverage,
        interactive_tm_gc=interactive_tm_gc,
        interactive_heatmap=interactive_heatmap,
        interactive_dimer_heatmap=interactive_dimer_heatmap,
    )

    return html


def generate_technical_report(
    results_dir: str,
    output_file: Optional[str] = None,
    interactive: bool = False,
) -> TechnicalReportData:
    """
    Generate comprehensive technical report from pipeline results.

    Args:
        results_dir: Path to results directory
        output_file: Output HTML file path (optional)
        interactive: If True, include interactive Plotly charts (requires plotly)

    Returns:
        TechnicalReportData object

    Example:
        data = generate_technical_report('results/', 'technical_report.html')
        print(f"Grade: {data.quality.grade.value}")

        # With interactive charts
        data = generate_technical_report('results/', 'report.html', interactive=True)
    """
    logger.info(f"Generating technical report for {results_dir}")

    # Collect all data
    data = collect_technical_report_data(results_dir)

    # Render and save if output file specified
    if output_file:
        html = render_technical_report(data, interactive=interactive)
        output_path = Path(output_file).resolve()

        try:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_text(html, encoding='utf-8')
            logger.info(f"Technical report saved to: {output_path}")
        except (PermissionError, OSError) as e:
            logger.error(f"Failed to write report to {output_path}: {e}")
            raise

    return data
