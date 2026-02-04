"""
Report generation module for NeoSWGA.

Provides comprehensive reporting capabilities:
- Executive summary (one-page quality assessment)
- Technical report (detailed multi-page analysis)
- Interactive Plotly visualizations (optional)
- Data export (CSV, JSON)

Usage:
    from neoswga.core.report import generate_executive_summary, generate_technical_report

    # Quick one-page summary
    generate_executive_summary(
        results_dir='results/',
        output_file='summary.html'
    )

    # Full technical report
    generate_technical_report(
        results_dir='results/',
        output_file='technical_report.html'
    )

    # Interactive report (requires plotly)
    generate_executive_summary(
        results_dir='results/',
        output_file='summary.html',
        interactive=True
    )
"""

from neoswga.core.report.executive_summary import (
    generate_executive_summary,
    ExecutiveSummary,
)
from neoswga.core.report.technical_report import (
    generate_technical_report,
    TechnicalReportData,
)
from neoswga.core.report.metrics import (
    PipelineMetrics,
    PrimerMetrics,
    collect_pipeline_metrics,
)
from neoswga.core.report.quality import (
    QualityGrade,
    calculate_quality_grade,
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

__all__ = [
    # Executive summary
    'generate_executive_summary',
    'ExecutiveSummary',
    # Technical report
    'generate_technical_report',
    'TechnicalReportData',
    # Metrics
    'PipelineMetrics',
    'PrimerMetrics',
    'collect_pipeline_metrics',
    # Quality
    'QualityGrade',
    'calculate_quality_grade',
    # Visualizations
    'is_plotly_available',
    'render_filtering_funnel',
    'render_component_radar',
    'render_tm_gc_distribution',
    'render_coverage_specificity_scatter',
    'render_primer_heatmap',
    'render_dimer_network_heatmap',
    'render_dimer_network_graph',
]
