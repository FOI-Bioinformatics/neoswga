"""
Report generation module for NeoSWGA.

Provides comprehensive reporting capabilities:
- Executive summary (one-page quality assessment)
- Technical report (detailed multi-page analysis)
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
]
