"""
Shared utilities for report generation.

Provides common functions for HTML rendering, CSS classes, and format string safety.
"""

from html import escape as html_escape
from typing import Any, Dict, Optional

from neoswga.core.report.quality import QualityGrade


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
