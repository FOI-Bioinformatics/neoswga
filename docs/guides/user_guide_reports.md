# NeoSWGA Report Generation User Guide

## Introduction

After running the NeoSWGA pipeline, you can generate comprehensive quality reports to assess your primer set. Reports provide:

- Letter grade (A-F) for overall quality
- Detailed metrics for each quality component
- Per-primer analysis
- Actionable recommendations

## Quick Start

### Generate a Report

```bash
# Navigate to your project directory
cd my_swga_project

# Generate an executive summary
neoswga report -d results/

# The report will be saved as results/report.html
```

Open the generated HTML file in any web browser to view your results.

## Report Types

### Executive Summary (Default)

A one-page overview suitable for quick assessment and sharing.

**Contents:**
- Overall quality grade with color-coded indicator
- Key metrics at a glance (coverage, specificity, uniformity)
- Primer table with essential properties
- Go/No-go recommendation

**Generate:**
```bash
neoswga report -d results/
```

### Technical Report

A comprehensive multi-page report for detailed analysis.

**Contents:**
- Pipeline execution summary
- Filtering funnel visualization
- Coverage analysis with binding site distribution
- Specificity analysis with enrichment calculations
- Thermodynamic properties
- Per-primer detailed profiles
- Primer-primer interaction analysis
- Extended recommendations

**Generate:**
```bash
neoswga report -d results/ --level full
```

### JSON Export

Machine-readable format for programmatic analysis.

**Generate:**
```bash
neoswga report -d results/ --format json
```

## Understanding Your Grade

### Grade Scale

| Grade | Score | Meaning |
|-------|-------|---------|
| A | 85-100% | Excellent - Ready for synthesis |
| B | 70-84% | Good - Minor concerns, likely suitable |
| C | 55-69% | Acceptable - Consider optimization |
| D | 40-54% | Poor - Optimization recommended |
| F | <40% | Critical - Do not proceed |

### Quality Components

Your grade is calculated from five weighted components:

#### 1. Coverage (35% weight)

How well your primers cover the target genome.

| Rating | Coverage | Interpretation |
|--------|----------|----------------|
| Excellent | >95% | Complete genome coverage |
| Good | 85-95% | Minor gaps acceptable |
| Acceptable | 70-84% | Some regions may be underrepresented |
| Poor | 50-69% | Significant coverage gaps |
| Critical | <50% | Inadequate for SWGA |

#### 2. Specificity (30% weight)

How well your primers discriminate target from background.

| Rating | Enrichment | Interpretation |
|--------|------------|----------------|
| Excellent | >500x | Outstanding target enrichment |
| Good | 100-500x | Good enrichment expected |
| Acceptable | 50-100x | Moderate enrichment |
| Poor | 20-50x | Low enrichment |
| Critical | <20x | Minimal target preference |

#### 3. Uniformity (20% weight)

How evenly primers bind across the genome (based on Gini index).

| Rating | Uniformity Score | Interpretation |
|--------|------------------|----------------|
| Excellent | >85% | Very even binding distribution |
| Good | 70-85% | Mostly uniform coverage |
| Acceptable | 55-70% | Some clustering observed |
| Poor | 40-55% | Significant clustering |
| Critical | <40% | Highly uneven amplification expected |

#### 4. Thermodynamics (10% weight)

How consistent the melting temperatures are across primers.

| Rating | Tm Range | Interpretation |
|--------|----------|----------------|
| Excellent | <3C | Highly consistent |
| Good | 3-5C | Good consistency |
| Acceptable | 5-8C | Some variation |
| Poor | 8-12C | Inconsistent |
| Critical | >12C | May cause amplification bias |

#### 5. Dimer Risk (5% weight)

Potential for primer-dimer formation.

| Rating | Risk Level | Interpretation |
|--------|------------|----------------|
| Excellent | <10% | Minimal dimer formation expected |
| Good | 10-20% | Low dimer risk |
| Acceptable | 20-35% | Moderate dimer risk |
| Poor | 35-50% | Significant dimer formation likely |
| Critical | >50% | Severe dimer problems expected |

## Validation Mode

Before generating a report, you can validate your results directory:

```bash
neoswga report -d results/ --check
```

This checks for:
- Required files (step4_improved_df.csv)
- Optional files (params.json, filter_stats.json)
- Data completeness

**Example output:**
```
Validating results directory...
  WARNING: Optional file missing: step3_df.csv

Validation passed with 1 warning(s).
```

## Common Workflows

### Workflow 1: Quick Quality Check

```bash
# Run pipeline
neoswga count-kmers -j params.json
neoswga filter -j params.json
neoswga score -j params.json
neoswga optimize -j params.json

# Generate quick report
neoswga report -d results/

# View in browser
open results/report.html
```

### Workflow 2: Detailed Analysis

```bash
# Validate results first
neoswga report -d results/ --check

# Generate full technical report
neoswga report -d results/ --level full -o analysis_report.html

# Also export JSON for further analysis
neoswga report -d results/ --format json -o metrics.json
```

### Workflow 3: Comparing Multiple Runs

```bash
# Generate reports for each run
neoswga report -d run1/results/ -o run1_report.html
neoswga report -d run2/results/ -o run2_report.html
neoswga report -d run3/results/ -o run3_report.html

# Export JSON for programmatic comparison
neoswga report -d run1/results/ --format json -o run1.json
neoswga report -d run2/results/ --format json -o run2.json
neoswga report -d run3/results/ --format json -o run3.json
```

### Workflow 4: Using Python API

```python
from neoswga.core.report import (
    generate_executive_summary,
    collect_pipeline_metrics,
    calculate_quality_grade,
)

# Generate report programmatically
summary = generate_executive_summary("results/", "report.html")

# Access metrics directly
metrics = collect_pipeline_metrics("results/")
quality = calculate_quality_grade(metrics)

# Print summary
print(f"Grade: {quality.grade.value}")
print(f"Score: {quality.composite_score:.1%}")
print(f"Primers: {metrics.primer_count}")
print(f"Recommendation: {quality.recommendation}")

# Analyze individual primers
for primer in metrics.primers:
    print(f"  {primer.sequence}: Tm={primer.tm:.1f}C, Specificity={primer.specificity:.0f}x")
```

## Troubleshooting

### "Required file missing: step4_improved_df.csv"

The pipeline hasn't completed. Run all steps:

```bash
neoswga count-kmers -j params.json
neoswga filter -j params.json
neoswga score -j params.json
neoswga optimize -j params.json
```

### "No primers found in results"

The optimization step didn't produce any primers. Check:
- Filter settings may be too strict
- Try increasing `max_primer` in params.json
- Check filter statistics for where primers are being eliminated

### Low Coverage Score

Possible causes:
- Insufficient primer count
- Target genome has repetitive regions
- Filter settings eliminating good candidates

Solutions:
- Increase `num_primers` in params.json
- Try longer primer lengths (min_k, max_k)
- Relax filtering thresholds

### Low Specificity Score

Possible causes:
- Target and background genomes are closely related
- Too few unique k-mers in target

Solutions:
- Use `--optimization-method=background-aware`
- Try different primer lengths
- Consider alternative background genomes

### Poor Uniformity

Possible causes:
- Genome has GC-biased regions
- Primer binding is clustered

Solutions:
- Use adaptive GC filtering
- Try `--optimization-method=dominating-set`
- Increase `num_primers` for better distribution

## Output File Locations

| File | Location | Description |
|------|----------|-------------|
| Executive Summary | `results/report.html` | Default HTML report |
| Technical Report | `results/technical_report.html` | Full analysis |
| JSON Export | `results/report.json` | Machine-readable |
| Custom | Specified with `-o` | Your chosen path |

## Getting Help

```bash
# View report command help
neoswga report --help

# Validate installation
neoswga validate --quick
```
