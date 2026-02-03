#!/usr/bin/env python3
"""
HTML report generation for SWGA simulation.

Creates comprehensive HTML reports with embedded plots, metrics, and recommendations.

Usage:
    from neoswga.core.simulation_report import generate_html_report

    generate_html_report(simulation_result, output_file)
"""

import base64
from io import BytesIO
from pathlib import Path
from datetime import datetime
import logging

logger = logging.getLogger(__name__)


HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SWGA Simulation Report</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}

        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;
            line-height: 1.6;
            color: #333;
            background: #f5f5f5;
            padding: 20px;
        }}

        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            border-radius: 8px;
            overflow: hidden;
        }}

        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            text-align: center;
        }}

        .header h1 {{
            font-size: 2.5em;
            margin-bottom: 10px;
        }}

        .header p {{
            opacity: 0.9;
            font-size: 1.1em;
        }}

        .recommendation {{
            padding: 30px;
            text-align: center;
            font-size: 1.5em;
            font-weight: bold;
            border-bottom: 3px solid #eee;
        }}

        .recommendation.excellent {{
            background: #d4edda;
            color: #155724;
        }}

        .recommendation.good {{
            background: #fff3cd;
            color: #856404;
        }}

        .recommendation.fair {{
            background: #f8d7da;
            color: #721c24;
        }}

        .recommendation.poor {{
            background: #f8d7da;
            color: #721c24;
        }}

        .metrics {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            padding: 30px;
        }}

        .metric-card {{
            background: #f8f9fa;
            border-radius: 8px;
            padding: 20px;
            border-left: 4px solid #667eea;
        }}

        .metric-card h3 {{
            color: #667eea;
            font-size: 0.9em;
            text-transform: uppercase;
            margin-bottom: 10px;
        }}

        .metric-card .value {{
            font-size: 2.5em;
            font-weight: bold;
            color: #333;
        }}

        .metric-card .label {{
            color: #666;
            font-size: 0.9em;
            margin-top: 5px;
        }}

        .section {{
            padding: 30px;
            border-bottom: 1px solid #eee;
        }}

        .section h2 {{
            color: #667eea;
            margin-bottom: 20px;
            font-size: 1.8em;
        }}

        .section h3 {{
            color: #333;
            margin: 20px 0 10px 0;
            font-size: 1.3em;
        }}

        .section p {{
            margin-bottom: 15px;
            color: #666;
        }}

        .section ul {{
            margin-left: 20px;
            margin-bottom: 15px;
        }}

        .section li {{
            margin-bottom: 8px;
            color: #666;
        }}

        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}

        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}

        th {{
            background: #667eea;
            color: white;
            font-weight: 600;
        }}

        tr:hover {{
            background: #f5f5f5;
        }}

        .plot {{
            margin: 20px 0;
            text-align: center;
        }}

        .plot img {{
            max-width: 100%;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}

        .alert {{
            padding: 15px;
            margin: 15px 0;
            border-radius: 4px;
            border-left: 4px solid;
        }}

        .alert.critical {{
            background: #f8d7da;
            border-color: #dc3545;
            color: #721c24;
        }}

        .alert.warning {{
            background: #fff3cd;
            border-color: #ffc107;
            color: #856404;
        }}

        .alert.info {{
            background: #d1ecf1;
            border-color: #17a2b8;
            color: #0c5460;
        }}

        .alert.success {{
            background: #d4edda;
            border-color: #28a745;
            color: #155724;
        }}

        .footer {{
            background: #f8f9fa;
            padding: 20px;
            text-align: center;
            color: #666;
            font-size: 0.9em;
        }}

        .badge {{
            display: inline-block;
            padding: 4px 8px;
            border-radius: 3px;
            font-size: 0.8em;
            font-weight: bold;
            margin-right: 5px;
        }}

        .badge.success {{
            background: #28a745;
            color: white;
        }}

        .badge.warning {{
            background: #ffc107;
            color: #333;
        }}

        .badge.danger {{
            background: #dc3545;
            color: white;
        }}

        .progress-bar {{
            background: #e9ecef;
            border-radius: 4px;
            height: 20px;
            overflow: hidden;
            margin: 10px 0;
        }}

        .progress-fill {{
            background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
            height: 100%;
            display: flex;
            align-items: center;
            justify-content: center;
            color: white;
            font-size: 0.8em;
            font-weight: bold;
        }}
    </style>
</head>
<body>
    <div class="container">
        <!-- Header -->
        <div class="header">
            <h1>SWGA Simulation Report</h1>
            <p>Comprehensive Primer Set Validation</p>
            <p style="font-size: 0.9em; margin-top: 10px;">Generated: {timestamp}</p>
        </div>

        <!-- Recommendation -->
        <div class="recommendation {rec_class}">
            <div>Overall Assessment: <span style="font-size: 1.2em;">{recommendation}</span></div>
            <div style="font-size: 0.8em; margin-top: 10px;">Composite Score: {composite:.2f}/1.00 | Confidence: {confidence:.0%}</div>
        </div>

        <!-- Key Metrics -->
        <div class="metrics">
            <div class="metric-card">
                <h3>Target Coverage</h3>
                <div class="value">{coverage:.1%}</div>
                <div class="label">Genome coverage fraction</div>
                <div class="progress-bar">
                    <div class="progress-fill" style="width: {coverage_pct:.0%}">{coverage:.1%}</div>
                </div>
            </div>

            <div class="metric-card">
                <h3>Uniformity</h3>
                <div class="value">{uniformity:.1%}</div>
                <div class="label">Coverage uniformity</div>
                <div class="progress-bar">
                    <div class="progress-fill" style="width: {uniformity_pct:.0%}">{uniformity:.1%}</div>
                </div>
            </div>

            <div class="metric-card">
                <h3>Enrichment</h3>
                <div class="value">{enrichment:.0f}×</div>
                <div class="label">Target/background ratio</div>
            </div>

            <div class="metric-card">
                <h3>Runtime</h3>
                <div class="value">{runtime:.1f}s</div>
                <div class="label">Simulation mode: {mode}</div>
            </div>
        </div>

        <!-- Target Genome Performance -->
        <div class="section">
            <h2>Target Genome Performance</h2>

            <h3>Coverage Analysis</h3>
            <ul>
                <li><strong>Coverage:</strong> {coverage:.1%} of genome covered</li>
                <li><strong>Uniformity:</strong> {uniformity:.1%} (higher is better)</li>
                <li><strong>Amplification:</strong> {target_amp:.2e}× estimated fold-amplification</li>
                <li><strong>Gaps:</strong> {n_gaps} under-covered regions identified</li>
            </ul>

            {gaps_html}
        </div>

        <!-- Background Genome Analysis -->
        <div class="section">
            <h2>Background Genome Analysis</h2>

            <ul>
                <li><strong>Coverage:</strong> {bg_coverage:.1%} of background genome covered</li>
                <li><strong>Amplification:</strong> {bg_amp:.2e}× estimated fold-amplification</li>
                <li><strong>Enrichment:</strong> {enrichment:.0f}× target/background ratio</li>
                <li><strong>Specificity Score:</strong> {specificity:.2f}/1.00</li>
            </ul>

            {specificity_alert}
        </div>

        <!-- Detailed Results -->
        <div class="section">
            <h2>Detailed Results</h2>

            <h3>Simulation Parameters</h3>
            <table>
                <tr><th>Parameter</th><th>Value</th></tr>
                <tr><td>Mode</td><td>{mode}</td></tr>
                <tr><td>Primers</td><td>{n_primers}</td></tr>
                <tr><td>Runtime</td><td>{runtime:.2f} seconds</td></tr>
                <tr><td>Confidence</td><td>{confidence:.0%}</td></tr>
            </table>

            {details_html}
        </div>

        <!-- Recommendations -->
        {recommendations_html}

        <!-- Footer -->
        <div class="footer">
            <p><strong>NeoSWGA Simulation Report</strong></p>
            <p>Generated with NeoSWGA v3.0 - Advanced SWGA Primer Design</p>
            <p>Report generated: {timestamp}</p>
        </div>
    </div>
</body>
</html>
"""


def generate_html_report(result, output_file: str, analysis=None):
    """
    Generate comprehensive HTML report.

    Args:
        result: SimulationResult from SwgaSimulator
        output_file: Output HTML file path
        analysis: Optional ComprehensiveAnalysis
    """
    logger.info("Generating HTML report...")

    # Prepare data
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    # Recommendation styling
    rec_class = result.recommendation.lower()

    # Coverage percentages for progress bars
    coverage_pct = result.target_coverage * 100
    uniformity_pct = result.target_uniformity * 100

    # Gaps HTML
    if result.target_gaps:
        gaps_html = "<h3>Largest Gaps</h3><ul>"
        for i, gap in enumerate(sorted(result.target_gaps, key=lambda g: g['length'], reverse=True)[:5], 1):
            gaps_html += f"<li>Gap {i}: {gap['start']:,} - {gap['end']:,} bp ({gap['length']:,} bp)</li>"
        gaps_html += "</ul>"
    else:
        gaps_html = '<div class="alert success">No significant gaps detected - excellent coverage!</div>'

    # Specificity alert
    if result.enrichment < 10:
        specificity_alert = '<div class="alert critical"><strong>Low Specificity:</strong> Enrichment < 10× - primers bind frequently to background genome.</div>'
    elif result.enrichment < 100:
        specificity_alert = '<div class="alert warning"><strong>Moderate Specificity:</strong> Enrichment < 100× - acceptable but could be improved.</div>'
    elif result.enrichment > 1000:
        specificity_alert = '<div class="alert success"><strong>Excellent Specificity:</strong> Enrichment > 1,000× - primers are highly specific!</div>'
    else:
        specificity_alert = '<div class="alert info"><strong>Good Specificity:</strong> Enrichment between 100-1,000× - suitable for most applications.</div>'

    # Details HTML
    details_html = ""
    if hasattr(result, 'details') and result.details:
        details_html = "<h3>Additional Details</h3><ul>"
        for key, value in result.details.items():
            if isinstance(value, (int, float)):
                details_html += f"<li><strong>{key.replace('_', ' ').title()}:</strong> {value}</li>"
        details_html += "</ul>"

    # Recommendations HTML
    recommendations_html = ""
    if analysis and analysis.recommendations:
        recommendations_html = '<div class="section"><h2>Recommendations</h2>'
        for rec in analysis.recommendations:
            if 'excellent' in rec.lower() or 'good' in rec.lower():
                alert_class = 'success'
            elif 'low' in rec.lower() or 'poor' in rec.lower():
                alert_class = 'critical'
            else:
                alert_class = 'info'

            recommendations_html += f'<div class="alert {alert_class}">{rec}</div>'
        recommendations_html += '</div>'

    # Fill template
    html = HTML_TEMPLATE.format(
        timestamp=timestamp,
        recommendation=result.recommendation,
        rec_class=rec_class,
        composite=result.composite_score,
        confidence=result.confidence,
        coverage=result.target_coverage,
        coverage_pct=coverage_pct,
        uniformity=result.target_uniformity,
        uniformity_pct=uniformity_pct,
        enrichment=result.enrichment,
        runtime=result.runtime,
        mode=result.mode.upper(),
        target_amp=result.target_amplification,
        n_gaps=len(result.target_gaps),
        gaps_html=gaps_html,
        bg_coverage=result.background_coverage,
        bg_amp=result.background_amplification,
        specificity=result.specificity_score,
        specificity_alert=specificity_alert,
        n_primers="N/A",  # Will be filled from context
        details_html=details_html,
        recommendations_html=recommendations_html
    )

    # Write file
    with open(output_file, 'w') as f:
        f.write(html)

    logger.info(f"HTML report saved to: {output_file}")


if __name__ == '__main__':
    print("This module is meant to be imported, not run directly")
    print("Usage: from neoswga.core.simulation_report import generate_html_report")
