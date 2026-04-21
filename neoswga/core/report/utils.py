"""
Shared utilities for report generation.

Provides common functions for HTML rendering, CSS classes, and format string safety.
"""

from html import escape as html_escape
from typing import Any, Dict, List, Optional

from neoswga.core.report.quality import QualityGrade


# Human-readable labels for validator issue codes emitted by
# OptimizationResult.validate() and the saturation / conditions-init
# extensions in unified_optimizer.run_optimization. Keep this dict in
# sync with the codes produced there — a missing entry falls back to
# the underscore-split code so the banner still renders.
VALIDATION_CODE_LABELS: Dict[str, str] = {
    "duplicate_primers": "Duplicate primers in set",
    "set_size_mismatch": "Primer set size differs from requested target",
    "coverage_below_threshold": "Aggregate foreground coverage below threshold",
    "per_target_coverage_below_threshold": "One or more targets below coverage threshold",
    "blacklist_primer_in_set": "Blacklist primer reached the final set",
    "reaction_conditions_init_failed": "Reaction conditions failed to initialise",
    "coverage_saturated_on_small_genome": "Coverage saturated on a small genome (metric unreliable)",
}


# CSS shared between the executive summary and the technical report for
# the validator-warnings banner. Templates inline this via `{VALIDATION_BANNER_CSS}`
# so tweaks to colour or spacing only live in one place. Each template
# can still override the `padding` / `font-size` on its own `.validation-banner`
# rule after including this block if its layout needs finer tuning.
VALIDATION_BANNER_CSS = """
        .validation-banner {
            padding: 14px 30px;
            border-bottom: 1px solid #e9ecef;
            font-size: 0.9em;
        }

        .validation-banner.level-error {
            background: #f8d7da;
            color: #721c24;
            border-left: 4px solid #dc3545;
        }

        .validation-banner.level-warning {
            background: #fff3cd;
            color: #856404;
            border-left: 4px solid #ffc107;
        }

        .validation-banner h3 {
            font-size: 0.95em;
            font-weight: 600;
            margin-bottom: 6px;
        }

        .validation-banner ul {
            margin: 0;
            padding-left: 20px;
        }

        .validation-banner li {
            margin-bottom: 3px;
            line-height: 1.5;
        }

        .validation-banner .code {
            font-family: 'SF Mono', 'Monaco', 'Inconsolata', monospace;
            font-size: 0.85em;
            font-weight: 600;
        }
"""


def render_validation_banner(issues: List[dict]) -> str:
    """Render validator warnings/errors as a visible HTML banner.

    Used by both executive_summary and technical_report to surface
    `per_target_coverage_below_threshold`, `blacklist_primer_in_set`,
    `coverage_saturated_on_small_genome`, and similar validator findings.
    Previously these lived only in step4_improved_df_validation.json
    and a Grade A display could hide a buried coverage warning.

    Returns an empty string when `issues` is empty so callers can
    unconditionally splice the result into their template without
    producing an empty `<div>`.
    """
    if not issues:
        return ""
    has_error = any(i.get("level") == "error" for i in issues)
    level_class = "level-error" if has_error else "level-warning"
    heading = (
        "Validation errors — review before ordering primers"
        if has_error
        else "Validation warnings"
    )
    items: List[str] = []
    for it in issues:
        code = str(it.get("code", "unknown"))
        detail = str(it.get("detail", ""))
        label = VALIDATION_CODE_LABELS.get(code, code.replace("_", " "))
        items.append(
            f"<li><span class=\"code\">{html_escape(code)}</span> — "
            f"{html_escape(label)}: {html_escape(detail)}</li>"
        )
    return (
        f'<div class="validation-banner {level_class}">'
        f'<h3>{html_escape(heading)}</h3>'
        f'<ul>{"".join(items)}</ul>'
        f'</div>'
    )


def get_version() -> str:
    """Get NeoSWGA version from package metadata."""
    try:
        from neoswga import __version__
        return __version__
    except ImportError:
        return "unknown"


# Color scheme for interactive charts (Plotly visualizations)
CHART_COLORS = {
    "primary": "#2c5282",
    "secondary": "#4299e1",
    "success": "#48bb78",
    "warning": "#ecc94b",
    "danger": "#f56565",
    "info": "#4299e1",
    "muted": "#a0aec0",
    # Gradient for funnel charts
    "funnel_gradient": [
        "#2c5282", "#3182ce", "#4299e1", "#63b3ed", "#90cdf4", "#bee3f8"
    ],
    # Radar chart styling
    "radar_fill": "rgba(66, 153, 225, 0.3)",
    "radar_line": "#2c5282",
    # Rating-based colors
    "excellent": "#48bb78",
    "good": "#4299e1",
    "acceptable": "#ecc94b",
    "poor": "#f56565",
}


def escape_for_html(value: Any) -> str:
    """
    Escape a value for safe HTML rendering.

    Converts value to string and escapes HTML special characters.

    Args:
        value: Any value to escape

    Returns:
        HTML-safe string
    """
    return html_escape(str(value))


def escape_format_braces(text: str) -> str:
    """
    Escape braces in text to prevent format string injection.

    User-controlled content with { or } could cause KeyError or expose
    template variables when used with str.format(). This escapes them
    to {{ and }}.

    Args:
        text: Text that may contain braces

    Returns:
        Text with braces escaped for safe use in format strings
    """
    return text.replace('{', '{{').replace('}', '}}')


def get_grade_colors(grade: QualityGrade) -> Dict[str, str]:
    """
    Get color scheme for a quality grade.

    Returns colors for:
    - grade_bg: Background color for grade section
    - grade_color: Text color for grade letter
    - grade_text: Alternative text color for grade display
    - rec_bg: Background color for recommendation
    - rec_border: Border color for recommendation
    - rec_color: Text color for recommendation

    Args:
        grade: Quality grade (A-F)

    Returns:
        Dictionary of color values (hex codes)
    """
    schemes = {
        QualityGrade.A: {
            "grade_bg": "#d4edda",
            "grade_color": "#155724",
            "grade_text": "#22543d",
            "rec_bg": "#d4edda",
            "rec_border": "#28a745",
            "rec_color": "#155724",
        },
        QualityGrade.B: {
            "grade_bg": "#cce5ff",
            "grade_color": "#004085",
            "grade_text": "#2a4365",
            "rec_bg": "#cce5ff",
            "rec_border": "#17a2b8",
            "rec_color": "#004085",
        },
        QualityGrade.C: {
            "grade_bg": "#fff3cd",
            "grade_color": "#856404",
            "grade_text": "#744210",
            "rec_bg": "#fff3cd",
            "rec_border": "#ffc107",
            "rec_color": "#856404",
        },
        QualityGrade.D: {
            "grade_bg": "#f8d7da",
            "grade_color": "#721c24",
            "grade_text": "#822727",
            "rec_bg": "#f8d7da",
            "rec_border": "#dc3545",
            "rec_color": "#721c24",
        },
        QualityGrade.F: {
            "grade_bg": "#721c24",
            "grade_color": "#fff",
            "grade_text": "#822727",
            "rec_bg": "#f8d7da",
            "rec_border": "#721c24",
            "rec_color": "#721c24",
        },
    }
    return schemes.get(grade, schemes[QualityGrade.C])


def get_rating_class(rating: str) -> str:
    """
    Get CSS class for a rating level.

    Args:
        rating: Rating string ("Excellent", "Good", "Acceptable", "Poor", "Critical")

    Returns:
        CSS class name (e.g., "rating-excellent")
    """
    mapping = {
        "Excellent": "rating-excellent",
        "Good": "rating-good",
        "Acceptable": "rating-acceptable",
        "Poor": "rating-poor",
        "Critical": "rating-critical",
    }
    return mapping.get(rating, "rating-acceptable")


def get_progress_class(rating: str) -> str:
    """
    Get CSS class for a progress bar based on rating.

    Args:
        rating: Rating string ("Excellent", "Good", "Acceptable", "Poor", "Critical")

    Returns:
        CSS class name (e.g., "progress-excellent")
    """
    mapping = {
        "Excellent": "progress-excellent",
        "Good": "progress-good",
        "Acceptable": "progress-acceptable",
        "Poor": "progress-poor",
        "Critical": "progress-poor",
    }
    return mapping.get(rating, "progress-acceptable")


# Grade descriptions for display
GRADE_DESCRIPTIONS = {
    QualityGrade.A: "Excellent",
    QualityGrade.B: "Good",
    QualityGrade.C: "Acceptable",
    QualityGrade.D: "Poor",
    QualityGrade.F: "Critical",
}


def get_grade_description(grade: QualityGrade) -> str:
    """
    Get human-readable description for a grade.

    Args:
        grade: Quality grade

    Returns:
        Description string (e.g., "Excellent" for grade A)
    """
    return GRADE_DESCRIPTIONS.get(grade, "Unknown")


# Threshold for considering enrichment as 100% on progress bar
ENRICHMENT_EXCELLENT_THRESHOLD = 500.0


def calculate_enrichment_bar_percent(enrichment: float) -> float:
    """
    Calculate progress bar percentage for enrichment value.

    Uses ENRICHMENT_EXCELLENT_THRESHOLD (500x) as 100%.

    Args:
        enrichment: Enrichment ratio

    Returns:
        Percentage (0-100)
    """
    return min(100.0, (enrichment / ENRICHMENT_EXCELLENT_THRESHOLD) * 100)
