# Report & Summary Step Evaluation and Improvement

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix coverage estimation inaccuracy, integrate step4_summary.json into reports, add missing test coverage for results_interpreter.py, improve validation, and enhance the technical report with actual set-level metrics.

**Architecture:** The report module (`neoswga/core/report/`) collects metrics from CSV files, grades quality (A-F), and renders HTML reports. The main gap is that coverage is estimated (~30kb per primer) rather than read from actual optimizer output. We just added `step4_summary.json` and expanded CSV columns in the optimization step -- the report module must now consume these real metrics instead of estimating them.

**Tech Stack:** Python 3.11, pytest, dataclasses, HTML templating, Plotly (optional)

---

## Team Structure (5 members + 1 lead)

| Role | Tasks | Focus |
|------|-------|-------|
| **Lead** | Coordinate, code review, final verification | Integration |
| **M1** | Task A (summary JSON ingestion) then Task F (results_interpreter tests) | Data flow + test coverage |
| **M2** | Task B (coverage from real data) | Accuracy improvement |
| **M3** | Task C (validation for new files) | Robustness |
| **M4** | Task D (exec summary enrichment) | Report quality |
| **M5** | Task E (technical report enrichment) | Report quality |

---

## Task A: Ingest step4_summary.json in Metrics Collection [M1, CRITICAL]

**Problem:** `collect_pipeline_metrics()` in `metrics.py:461-535` never reads `step4_summary.json` (which we just added in the optimizer improvement). This JSON contains real coverage, selectivity, gap metrics, and normalized scores from `OptimizationResult.to_dict()`. The report module estimates coverage as `n_primers * 30000 / genome_size` (line 349) which is inaccurate.

**Files:**
- Modify: `neoswga/core/report/metrics.py:461-535` -- load and merge summary JSON
- Modify: `neoswga/core/report/metrics.py:176-190` -- add fields to CoverageMetrics
- Test: `tests/report/test_metrics.py` -- add tests for JSON ingestion

**Steps:**

- [ ] **Step 1: Add `optimizer_metrics` field to CoverageMetrics**

In `metrics.py`, the `CoverageMetrics` dataclass (line 176) currently has limited fields. Add fields to store real optimizer output:

```python
@dataclass
class CoverageMetrics:
    """Coverage analysis metrics."""
    overall_coverage: float = 0.0
    covered_bases: int = 0
    total_bases: int = 0
    num_gaps: int = 0
    mean_gap: float = 0.0
    max_gap: float = 0.0
    gap_gini: float = 0.0
    gap_entropy: float = 0.0
    # New: from step4_summary.json
    from_optimizer: bool = False  # True if coverage comes from real optimizer data
```

- [ ] **Step 2: Add `_load_optimizer_summary()` helper**

Add a function after `_load_params()` (line 304):

```python
def _load_optimizer_summary(results_path: Path) -> Optional[Dict]:
    """Load step4 summary JSON if available."""
    summary_file = results_path / 'step4_improved_df_summary.json'
    if not summary_file.exists():
        return None
    try:
        with open(summary_file, encoding='utf-8') as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError) as e:
        logger.warning(f"Failed to load optimizer summary: {e}")
        return None
```

- [ ] **Step 3: Integrate summary into `collect_pipeline_metrics()`**

After line 515 (after `_calculate_uniformity_metrics`), add:

```python
    # Load optimizer summary for real metrics
    optimizer_summary = _load_optimizer_summary(results_path)
    if optimizer_summary and 'metrics' in optimizer_summary:
        opt_metrics = optimizer_summary['metrics']
        # Override estimated coverage with real optimizer data
        if 'fg_coverage' in opt_metrics:
            metrics.coverage.overall_coverage = opt_metrics['fg_coverage']
            metrics.coverage.from_optimizer = True
            if metrics.coverage.total_bases > 0:
                metrics.coverage.covered_bases = int(
                    opt_metrics['fg_coverage'] * metrics.coverage.total_bases
                )
        # Add gap metrics
        metrics.coverage.mean_gap = opt_metrics.get('mean_gap', 0.0)
        metrics.coverage.max_gap = opt_metrics.get('max_gap', 0.0)
        metrics.coverage.gap_gini = opt_metrics.get('gap_gini', 0.0)
        metrics.coverage.gap_entropy = opt_metrics.get('gap_entropy', 0.0)
        logger.info("Using real optimizer metrics for coverage data")
```

- [ ] **Step 4: Write tests for JSON ingestion**

Add to `tests/report/test_metrics.py`:

```python
class TestOptimizerSummaryIngestion:
    """Test loading and merging of step4_summary.json."""

    def test_load_optimizer_summary_missing_file(self, tmp_path):
        from neoswga.core.report.metrics import _load_optimizer_summary
        assert _load_optimizer_summary(tmp_path) is None

    def test_load_optimizer_summary_valid(self, tmp_path):
        import json
        from neoswga.core.report.metrics import _load_optimizer_summary
        summary = {'metrics': {'fg_coverage': 0.85, 'mean_gap': 5000.0}}
        (tmp_path / 'step4_improved_df_summary.json').write_text(json.dumps(summary))
        result = _load_optimizer_summary(tmp_path)
        assert result['metrics']['fg_coverage'] == 0.85

    def test_load_optimizer_summary_corrupt_json(self, tmp_path):
        from neoswga.core.report.metrics import _load_optimizer_summary
        (tmp_path / 'step4_improved_df_summary.json').write_text('not json{')
        assert _load_optimizer_summary(tmp_path) is None

    def test_coverage_uses_optimizer_data(self, tmp_path):
        """When summary JSON exists, coverage should come from optimizer, not estimate."""
        import json
        # Create minimal step4 CSV
        csv_content = "primer,score,coverage\nATCGATCG,1.0,0.9\n"
        (tmp_path / 'step4_improved_df.csv').write_text(csv_content)
        # Create params.json with genome size
        params = {'fg_genome': 'target.fna', 'fg_size': 100000}
        (tmp_path / 'params.json').write_text(json.dumps(params))
        # Create summary with real coverage
        summary = {
            'metrics': {
                'fg_coverage': 0.72,
                'mean_gap': 3500.0,
                'max_gap': 12000.0,
                'gap_gini': 0.35,
                'gap_entropy': 2.1,
            }
        }
        (tmp_path / 'step4_improved_df_summary.json').write_text(json.dumps(summary))

        from neoswga.core.report.metrics import collect_pipeline_metrics
        metrics = collect_pipeline_metrics(str(tmp_path))
        assert metrics.coverage.overall_coverage == 0.72
        assert metrics.coverage.from_optimizer is True
        assert metrics.coverage.mean_gap == 3500.0
```

- [ ] **Step 5: Run tests**

```bash
conda run -n neoswga3 python -m pytest tests/report/test_metrics.py -v --tb=short
```

---

## Task B: Use Real Coverage in Quality Grading [M2, MEDIUM]

**Problem:** `_calculate_coverage_metrics()` (metrics.py:326-353) estimates coverage as `n_primers * 30000 / genome_size`. After Task A injects real coverage from the summary JSON, the quality grade becomes more accurate. However, we should also handle the case where the step4 CSV itself now contains per-primer `coverage` and `coverage_uniformity` columns (added in the optimizer improvement). The `PrimerMetrics.from_row()` method (line 110) does not parse these new columns.

**Files:**
- Modify: `neoswga/core/report/metrics.py:110-148` -- parse new CSV columns in `from_row()`
- Modify: `neoswga/core/report/metrics.py:326-353` -- improve fallback estimation, label as estimated
- Test: `tests/report/test_metrics.py` -- test new column parsing

**Steps:**

- [ ] **Step 1: Extend PrimerMetrics to parse new step4 columns**

The step4 CSV now includes columns: `normalized_score`, `coverage`, `bg_coverage`, `coverage_uniformity`, `gap_gini`, `gap_entropy`, `dimer_risk_score`, `strand_alternation`, `strand_coverage_ratio`, `mean_tm`, `max_gap`. Add parsing for the ones that are useful per-primer:

In `PrimerMetrics.from_row()`, after the current return statement, add parsing for new fields that may be present in step4:

```python
# Parse additional step4 columns (set-level metrics repeated per primer)
dimer_score_val = _safe_float(row.get('dimer_score', row.get('dimer_risk_score', 0)))
```

Note: Most new columns are set-level (same value for all primers in a set). The `from_row()` method should parse `dimer_risk_score` as `dimer_score` if `dimer_score` is not already present.

- [ ] **Step 2: Improve `_calculate_coverage_metrics()` fallback**

When no optimizer summary is available but step4 CSV has a `coverage` column, use that. Update the function:

```python
def _calculate_coverage_metrics(
    primers: List[PrimerMetrics],
    params: Dict,
) -> CoverageMetrics:
    """Calculate coverage metrics from primer data."""
    coverage = CoverageMetrics()

    fg_size = params.get('fg_size', params.get('foreground_size', 0))
    if fg_size > 0:
        coverage.total_bases = fg_size

    # Check if primers have optimizer-provided coverage (from step4 CSV)
    # All primers in a set share the same coverage value
    if primers and hasattr(primers[0], '_raw_coverage'):
        csv_coverage = primers[0]._raw_coverage
        if csv_coverage > 0:
            coverage.overall_coverage = csv_coverage
            coverage.from_optimizer = True
            if coverage.total_bases > 0:
                coverage.covered_bases = int(csv_coverage * coverage.total_bases)
            return coverage

    # Fallback: estimate from primer count (less accurate)
    if coverage.total_bases > 0 and primers:
        n_primers = len(primers)
        avg_amp_range = 30000  # 30 kb average amplification range
        estimated_covered = min(n_primers * avg_amp_range, coverage.total_bases)
        coverage.overall_coverage = estimated_covered / coverage.total_bases
        coverage.covered_bases = int(estimated_covered)
        # Mark as estimated
        coverage.from_optimizer = False

    return coverage
```

- [ ] **Step 3: Write test for new column parsing**

```python
def test_primer_metrics_from_row_with_step4_columns():
    """PrimerMetrics should parse dimer_risk_score from step4 output."""
    row = {
        'primer': 'ATCGATCG',
        'dimer_risk_score': '0.15',
        'coverage': '0.85',
        'fg_freq': '0.001',
        'bg_freq': '0.0001',
    }
    pm = PrimerMetrics.from_row(row)
    assert pm.dimer_score == 0.15
```

- [ ] **Step 4: Run tests**

```bash
conda run -n neoswga3 python -m pytest tests/report/test_metrics.py -v --tb=short
```

---

## Task C: Validate step4_summary.json in Report Validation [M3, MEDIUM]

**Problem:** `validate_results_directory()` in `validation.py:80-140` checks for `step4_improved_df.csv` and optional files, but does not check for `step4_improved_df_summary.json` (which provides accurate coverage data). It should note its presence as an INFO and warn if missing.

**Files:**
- Modify: `neoswga/core/report/validation.py:122-139` -- add summary JSON to optional files
- Modify: `neoswga/core/report/validation.py:143-218` -- validate summary JSON content
- Test: `tests/report/test_e2e_report.py` or new `tests/report/test_validation.py`

**Steps:**

- [ ] **Step 1: Add summary JSON to optional file checks**

In `validate_results_directory()`, add to the `optional_files` list:

```python
    optional_files = [
        ("step3_df.csv", "Some scoring metrics may be unavailable"),
        ("step2_df.csv", "Filtering statistics may be incomplete"),
        ("filter_stats.json", "Filter funnel data will not be shown"),
        ("params.json", "Parameter information will not be shown"),
        ("step4_improved_df_summary.json", "Coverage will be estimated rather than measured"),
    ]
```

- [ ] **Step 2: Add content validation for summary JSON**

In `validate_metrics()`, add a check for coverage source:

```python
    # Check coverage data source
    if metrics.coverage and not metrics.coverage.from_optimizer:
        result.add_warning(
            "Coverage is estimated (no optimizer summary found). "
            "Run optimization step to generate accurate coverage data.",
            "coverage"
        )
    elif metrics.coverage and metrics.coverage.from_optimizer:
        result.add_info(
            "Using measured coverage from optimizer output",
            "coverage"
        )
```

- [ ] **Step 3: Write test**

```python
def test_validation_warns_estimated_coverage():
    from neoswga.core.report.validation import validate_metrics
    from neoswga.core.report.metrics import PipelineMetrics, CoverageMetrics, PrimerMetrics

    metrics = PipelineMetrics(results_dir='/tmp', generated_at='now')
    metrics.primers = [PrimerMetrics(sequence='ATCG', length=4, gc_content=0.5,
                                     tm=30.0, fg_freq=0.001, bg_freq=0.0001,
                                     fg_sites=10, bg_sites=1, gini=0.3,
                                     specificity=10.0)]
    metrics.coverage = CoverageMetrics(overall_coverage=0.5, from_optimizer=False)

    result = validate_metrics(metrics)
    warning_msgs = [i.message for i in result.warnings]
    assert any('estimated' in m.lower() for m in warning_msgs)
```

- [ ] **Step 4: Run tests**

```bash
conda run -n neoswga3 python -m pytest tests/report/ -v --tb=short
```

---

## Task D: Enrich Executive Summary with Optimizer Metrics [M4, MEDIUM]

**Problem:** The executive summary shows 4 metric cards (Coverage%, Enrichment, Uniformity, Dimer Risk). It should also show gap metrics and normalized score when available from the summary JSON. The coverage card should indicate whether it is measured or estimated.

**Files:**
- Modify: `neoswga/core/report/executive_summary.py` -- add "measured/estimated" badge to coverage card, add normalized score
- Test: `tests/report/test_e2e_report.py` -- verify new content appears

**Steps:**

- [ ] **Step 1: Add coverage source indicator**

In the coverage metric card section of `render_executive_summary()`, add a small badge showing "Measured" (green) or "Estimated" (yellow) next to the coverage value. Find the coverage card HTML and add:

```python
coverage_source = "Measured" if metrics.coverage.from_optimizer else "Estimated"
coverage_badge_color = "#48bb78" if metrics.coverage.from_optimizer else "#ecc94b"
```

Then add a `<span>` badge after the coverage percentage in the metric card HTML.

- [ ] **Step 2: Add normalized score to header**

When summary JSON is available, show the normalized_score (0-1) in the grade section as a secondary metric. This requires passing the optimizer score through.

- [ ] **Step 3: Add gap statistics section**

When gap_gini, gap_entropy, mean_gap, max_gap are available (from Task A), render a small "Gap Analysis" subsection below the metric cards.

- [ ] **Step 4: Write E2E test**

```python
def test_executive_summary_shows_coverage_source(self, bacillus_results_dir):
    """Report should indicate whether coverage is measured or estimated."""
    # This test checks the HTML output contains a coverage source indicator
    from neoswga.core.report.executive_summary import generate_executive_summary
    html = generate_executive_summary(str(bacillus_results_dir))
    assert 'Estimated' in html or 'Measured' in html
```

- [ ] **Step 5: Run tests**

```bash
conda run -n neoswga3 python -m pytest tests/report/test_e2e_report.py -v --tb=short
```

---

## Task E: Enrich Technical Report with Set-Level Metrics [M5, MEDIUM]

**Problem:** The technical report (`technical_report.py`, 1825 lines) does not display gap analysis or set-level optimizer metrics. It defines `GapInfo` (line 57) but never populates it. Now that Task A makes real gap metrics available, the technical report should display them.

**Files:**
- Modify: `neoswga/core/report/technical_report.py` -- add gap analysis section, populate GapInfo
- Test: `tests/report/test_e2e_report.py` -- verify gap section appears

**Steps:**

- [ ] **Step 1: Add gap analysis section to technical report HTML**

Find the coverage analysis section in the technical report template. After it, add a "Gap Analysis" subsection that displays:
- Mean gap size (formatted as kb)
- Max gap size (formatted as kb)
- Gap Gini coefficient
- Gap entropy
- A small table or metric cards for these values

Only render this section when `metrics.coverage.from_optimizer` is True.

- [ ] **Step 2: Add optimizer summary section**

Add a new section "Optimizer Performance" that shows:
- Optimizer name and method
- Normalized score
- Number of primers selected vs target
- Status (SUCCESS/PARTIAL)
- Extension reach used

Load this from the summary JSON if available.

- [ ] **Step 3: Add coverage source annotation**

In the existing coverage section of the technical report, add an annotation indicating whether the coverage value is "measured from primer binding positions" or "estimated from primer count".

- [ ] **Step 4: Write E2E test**

```python
def test_technical_report_shows_gap_analysis(self, bacillus_results_dir):
    """Technical report should include gap analysis when data is available."""
    from neoswga.core.report.technical_report import generate_technical_report
    html = generate_technical_report(str(bacillus_results_dir), interactive=False)
    # Should have coverage source annotation at minimum
    assert 'Estimated' in html or 'Measured' in html
```

- [ ] **Step 5: Run tests**

```bash
conda run -n neoswga3 python -m pytest tests/report/test_e2e_report.py -v --tb=short
```

---

## Task F: Add Tests for results_interpreter.py [M1, after Task A]

**Problem:** `results_interpreter.py` (556 lines, 3 classes, 8 methods) has zero test coverage. This module is used by `neoswga interpret` and provides console-based quality assessment.

**Files:**
- Create: `tests/test_results_interpreter.py`
- Reference: `neoswga/core/results_interpreter.py`

**Steps:**

- [ ] **Step 1: Test QualityRating enum and MetricAssessment**

```python
import pytest
from neoswga.core.results_interpreter import (
    QualityRating, MetricAssessment, ResultsReport, ResultsInterpreter
)

def test_quality_rating_values():
    assert QualityRating.EXCELLENT.value == "EXCELLENT"
    assert QualityRating.CRITICAL.value == "CRITICAL"

def test_metric_assessment_creation():
    ma = MetricAssessment(
        name="Coverage", value=0.85, rating=QualityRating.GOOD,
        unit="%", context="85% of genome covered",
        threshold_info="Good: >85%"
    )
    assert ma.rating == QualityRating.GOOD
```

- [ ] **Step 2: Test ResultsInterpreter with synthetic CSV data**

```python
def test_interpreter_loads_step4(tmp_path):
    """Interpreter should load primers from step4 CSV."""
    csv_content = (
        "primer,score,coverage,bg_coverage,selectivity,mean_gap,optimizer\n"
        "ATCGATCG,1.0,0.85,0.1,8.5,5000,hybrid\n"
        "GCTAGCTA,0.9,0.85,0.1,8.5,5000,hybrid\n"
    )
    (tmp_path / 'step4_improved_df.csv').write_text(csv_content)
    params = {'fg_genome': 'target.fna', 'polymerase': 'phi29', 'reaction_temp': 30}
    import json
    (tmp_path / 'params.json').write_text(json.dumps(params))

    interpreter = ResultsInterpreter(str(tmp_path))
    report = interpreter.analyze()
    assert report is not None
    assert report.primer_count == 2
    assert len(report.assessments) > 0
```

- [ ] **Step 3: Test rating thresholds**

```python
def test_rate_metric_coverage():
    """Test coverage rating thresholds."""
    interpreter = ResultsInterpreter('/nonexistent')
    # Access the private method for unit testing
    rating = interpreter._rate_metric(0.96, 'coverage')
    assert rating == QualityRating.EXCELLENT

    rating = interpreter._rate_metric(0.40, 'coverage')
    assert rating == QualityRating.POOR
```

Note: Check whether `_rate_metric` exists as written or uses a different pattern. Read the interpreter code to verify the method signature.

- [ ] **Step 4: Test recommendation generation**

```python
def test_recommendation_excellent():
    interpreter = ResultsInterpreter('/nonexistent')
    rec = interpreter._generate_recommendation(QualityRating.EXCELLENT)
    assert 'synthesis' in rec.lower() or 'proceed' in rec.lower()

def test_recommendation_critical():
    interpreter = ResultsInterpreter('/nonexistent')
    rec = interpreter._generate_recommendation(QualityRating.CRITICAL)
    assert 'not proceed' in rec.lower() or 'do not' in rec.lower()
```

- [ ] **Step 5: Test print_report does not crash**

```python
def test_print_report_no_crash(tmp_path, capsys):
    """print_report should produce output without crashing."""
    csv_content = (
        "primer,score,coverage,selectivity,mean_gap,optimizer\n"
        "ATCGATCG,1.0,0.85,8.5,5000,hybrid\n"
    )
    (tmp_path / 'step4_improved_df.csv').write_text(csv_content)
    import json
    (tmp_path / 'params.json').write_text(json.dumps({'fg_genome': 'test.fna'}))

    interpreter = ResultsInterpreter(str(tmp_path))
    report = interpreter.analyze()
    if report:
        interpreter.print_report(report)
        captured = capsys.readouterr()
        assert len(captured.out) > 0
```

- [ ] **Step 6: Run tests**

```bash
conda run -n neoswga3 python -m pytest tests/test_results_interpreter.py -v --tb=short
```

---

## Dependencies

```
A (summary JSON ingestion) -- independent, start first
B (real coverage) -- can start in parallel with A
C (validation) -- depends on A (needs from_optimizer field)
D (exec summary) -- depends on A (needs coverage source)
E (technical report) -- depends on A (needs gap metrics)
F (interpreter tests) -- independent of A-E

Timeline:
  A, B, F  (parallel, no dependencies)
       |
       v
  C, D, E  (after A completes, parallel with each other)
```

## Verification Protocol (Lead)

```bash
# After all tasks:
conda run -n neoswga3 python -m pytest tests/report/ tests/test_results_interpreter.py -v --tb=short

# Generate reports on plasmid example to verify output
cd examples/plasmid_example
neoswga report -d . --check               # Validation passes
neoswga report -d . --level full           # Technical report generates
neoswga interpret -d .                     # Interpreter runs without crash

# Verify summary JSON appears in report output
grep -l "Measured\|Estimated" examples/plasmid_example/executive_summary.html
```

## Critical Files

| File | Role |
|------|------|
| `neoswga/core/report/metrics.py` | Metrics collection, summary JSON ingestion |
| `neoswga/core/report/quality.py` | Quality grading (uses coverage data) |
| `neoswga/core/report/executive_summary.py` | One-page report with coverage badge |
| `neoswga/core/report/technical_report.py` | Multi-page report with gap analysis |
| `neoswga/core/report/validation.py` | Input validation for new files |
| `neoswga/core/results_interpreter.py` | Console interpreter (needs tests) |
| `tests/report/test_metrics.py` | Metrics tests |
| `tests/report/test_e2e_report.py` | End-to-end report tests |
| `tests/test_results_interpreter.py` | New interpreter tests |
