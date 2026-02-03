# NeoSWGA Report Module API Reference

## Module: `neoswga.core.report`

### Public Exports

```python
from neoswga.core.report import (
    # Executive summary
    generate_executive_summary,
    ExecutiveSummary,
    # Technical report
    generate_technical_report,
    TechnicalReportData,
    # Metrics
    PipelineMetrics,
    PrimerMetrics,
    collect_pipeline_metrics,
    # Quality
    QualityGrade,
    calculate_quality_grade,
)
```

---

## neoswga.core.report.metrics

### Functions

#### `collect_pipeline_metrics`

```python
def collect_pipeline_metrics(results_dir: str) -> PipelineMetrics
```

Collect all metrics from a pipeline results directory.

**Parameters:**
| Name | Type | Description |
|------|------|-------------|
| `results_dir` | `str` | Path to directory containing pipeline output files |

**Returns:**
| Type | Description |
|------|-------------|
| `PipelineMetrics` | Collected metrics from all available files |

**Raises:**
| Exception | Condition |
|-----------|-----------|
| `FileNotFoundError` | If results_dir doesn't exist |

**Files Read:**
- `step4_improved_df.csv` (required) - Final primer set
- `step3_df.csv` (optional) - Scored primers
- `step2_df.csv` (optional) - Filtered primers
- `params.json` (optional) - Pipeline parameters
- `filter_stats.json` (optional) - Filtering statistics

---

#### `_safe_float`

```python
def _safe_float(value: Any, default: float = 0.0) -> float
```

Safely convert a value to float, handling edge cases.

**Parameters:**
| Name | Type | Default | Description |
|------|------|---------|-------------|
| `value` | `Any` | - | Value to convert |
| `default` | `float` | `0.0` | Default if conversion fails |

**Handles:**
- `None` and empty strings
- `NaN`, `Inf`, `-Inf` (returns default)
- Invalid string formats

---

#### `_safe_int`

```python
def _safe_int(value: Any, default: int = 0) -> int
```

Safely convert a value to int.

**Parameters:**
| Name | Type | Default | Description |
|------|------|---------|-------------|
| `value` | `Any` | - | Value to convert |
| `default` | `int` | `0` | Default if conversion fails |

---

### Classes

#### `PrimerMetrics`

```python
@dataclass
class PrimerMetrics:
    sequence: str
    length: int
    gc_content: float
    tm: float
    fg_freq: float
    bg_freq: float
    fg_sites: int
    bg_sites: int
    gini: float
    specificity: float
    amp_pred: float = 0.0
    dimer_score: float = 0.0
    hairpin_dg: float = 0.0
    self_dimer_dg: float = 0.0
    three_prime_stability: float = 0.0
    strand_ratio: float = 1.0
    quality_rank: int = 0
```

**Class Methods:**

##### `from_row`

```python
@classmethod
def from_row(cls, row: Dict[str, str]) -> PrimerMetrics
```

Create PrimerMetrics from a CSV row dictionary.

**Column Mappings:**
| Field | Primary Key | Alternative Key |
|-------|-------------|-----------------|
| sequence | `sequence` | `primer` |
| tm | `tm` | `Tm` |
| gc_content | `gc` | `gc_content` |
| fg_sites | `fg_count` | `fg_sites` |
| bg_sites | `bg_count` | `bg_sites` |
| gini | `gini` | `gini_index` |
| amp_pred | `amp_pred` | `on.target.pred` |

---

#### `PipelineMetrics`

```python
@dataclass
class PipelineMetrics:
    primers: List[PrimerMetrics]
    primer_count: int
    target_genome_size: int
    background_genome_size: int
    parameters: Dict[str, Any]
    filtering: Optional[FilteringStats] = None
    coverage: Optional[CoverageMetrics] = None
    specificity: Optional[SpecificityMetrics] = None
    uniformity: Optional[UniformityMetrics] = None
    thermodynamics: Optional[ThermodynamicMetrics] = None
```

---

#### `FilteringStats`

```python
@dataclass
class FilteringStats:
    total_kmers: int
    after_frequency: int
    after_background: int
    after_gini: int
    after_thermodynamic: int
    after_complexity: int
    final_candidates: int
```

**Methods:**

##### `as_funnel`

```python
def as_funnel(self) -> List[Tuple[str, int]]
```

Convert to funnel visualization format.

**Returns:**
```python
[
    ("Total k-mers", 2097152),
    ("After frequency filter", 125000),
    ("After background filter", 8500),
    ...
]
```

---

#### `CoverageMetrics`

```python
@dataclass
class CoverageMetrics:
    overall_coverage: float  # 0-1 scale
    mean_spacing: float      # bp between sites
    max_gap: int             # largest gap in bp
    sites_per_mb: float      # binding site density
```

---

#### `SpecificityMetrics`

```python
@dataclass
class SpecificityMetrics:
    enrichment_ratio: float     # fg_density / bg_density
    target_density: float       # sites per bp in target
    background_density: float   # sites per bp in background
    min_specificity: float      # worst primer specificity
    max_specificity: float      # best primer specificity
```

---

#### `UniformityMetrics`

```python
@dataclass
class UniformityMetrics:
    mean_gini: float    # average Gini index
    max_gini: float     # worst Gini index
    min_gini: float     # best Gini index
```

---

#### `ThermodynamicMetrics`

```python
@dataclass
class ThermodynamicMetrics:
    mean_tm: float
    min_tm: float
    max_tm: float
    tm_range: float     # max_tm - min_tm
    mean_gc: float
```

---

## neoswga.core.report.quality

### Functions

#### `calculate_quality_grade`

```python
def calculate_quality_grade(metrics: PipelineMetrics) -> QualityAssessment
```

Calculate quality grade from pipeline metrics.

**Parameters:**
| Name | Type | Description |
|------|------|-------------|
| `metrics` | `PipelineMetrics` | Collected pipeline metrics |

**Returns:**
| Type | Description |
|------|-------------|
| `QualityAssessment` | Complete quality assessment |

---

### Classes

#### `QualityGrade`

```python
class QualityGrade(Enum):
    A = "A"  # Excellent (>= 0.85)
    B = "B"  # Good (0.70-0.84)
    C = "C"  # Acceptable (0.55-0.69)
    D = "D"  # Poor (0.40-0.54)
    F = "F"  # Critical (< 0.40)
```

---

#### `GradeComponent`

```python
@dataclass
class GradeComponent:
    name: str              # Component name
    weight: float          # Weight in composite (0-1)
    raw_value: float       # Original metric value
    normalized_score: float # Normalized 0-1 score
    rating: str            # "Excellent", "Good", etc.
    description: str       # Human-readable description
```

---

#### `QualityAssessment`

```python
@dataclass
class QualityAssessment:
    grade: QualityGrade
    composite_score: float          # 0-1 weighted average
    components: List[GradeComponent]
    recommendation: str             # Brief recommendation
    recommendation_details: str     # Extended recommendation
    considerations: List[str]       # Additional notes
```

---

### Constants

#### Threshold Definitions

```python
COVERAGE_THRESHOLDS = {
    "excellent": 0.95,
    "good": 0.85,
    "acceptable": 0.70,
    "poor": 0.50,
}

ENRICHMENT_THRESHOLDS = {
    "excellent": 500.0,
    "good": 100.0,
    "acceptable": 50.0,
    "poor": 20.0,
}

UNIFORMITY_THRESHOLDS = {
    "excellent": 0.85,
    "good": 0.70,
    "acceptable": 0.55,
    "poor": 0.40,
}

TM_RANGE_THRESHOLDS = {
    "excellent": 3.0,
    "good": 5.0,
    "acceptable": 8.0,
    "poor": 12.0,
}

DIMER_THRESHOLDS = {
    "excellent": 0.1,
    "good": 0.2,
    "acceptable": 0.35,
    "poor": 0.5,
}
```

---

## neoswga.core.report.validation

### Functions

#### `validate_results_directory`

```python
def validate_results_directory(path: str) -> ValidationResult
```

Validate a results directory before report generation.

**Parameters:**
| Name | Type | Description |
|------|------|-------------|
| `path` | `str` | Path to results directory |

**Returns:**
| Type | Description |
|------|-------------|
| `ValidationResult` | Validation result with issues |

**Checks Performed:**
| Check | Level | Condition |
|-------|-------|-----------|
| Directory exists | ERROR | Path must exist |
| step4_improved_df.csv | ERROR | Required file |
| step3_df.csv | WARNING | Optional file |
| step2_df.csv | WARNING | Optional file |
| params.json | INFO | Optional file |

---

#### `validate_metrics`

```python
def validate_metrics(metrics: PipelineMetrics) -> ValidationResult
```

Validate collected metrics for completeness.

**Parameters:**
| Name | Type | Description |
|------|------|-------------|
| `metrics` | `PipelineMetrics` | Collected metrics |

**Checks Performed:**
| Check | Level | Condition |
|-------|-------|-----------|
| Primers present | ERROR | len(primers) > 0 |
| Genome size | WARNING | target_genome_size > 0 |
| Valid Tm values | WARNING | All Tm > 0 |

---

### Classes

#### `ValidationLevel`

```python
class ValidationLevel(Enum):
    ERROR = "error"      # Blocks report generation
    WARNING = "warning"  # Noted but doesn't block
    INFO = "info"        # Informational only
```

---

#### `ValidationIssue`

```python
@dataclass
class ValidationIssue:
    level: ValidationLevel
    message: str
    field: str = ""
```

---

#### `ValidationResult`

```python
class ValidationResult:
    issues: List[ValidationIssue]

    @property
    def errors(self) -> List[ValidationIssue]

    @property
    def warnings(self) -> List[ValidationIssue]

    @property
    def infos(self) -> List[ValidationIssue]

    @property
    def is_valid(self) -> bool  # True if no errors

    def add_error(self, message: str, field: str = "") -> None
    def add_warning(self, message: str, field: str = "") -> None
    def add_info(self, message: str, field: str = "") -> None
```

---

## neoswga.core.report.utils

### Functions

#### `escape_format_braces`

```python
def escape_format_braces(text: str) -> str
```

Escape braces to prevent format string injection.

**Parameters:**
| Name | Type | Description |
|------|------|-------------|
| `text` | `str` | Text that may contain braces |

**Returns:**
| Type | Description |
|------|-------------|
| `str` | Text with `{` -> `{{` and `}` -> `}}` |

---

#### `get_grade_colors`

```python
def get_grade_colors(grade: QualityGrade) -> Dict[str, str]
```

Get color scheme for a quality grade.

**Returns dict with keys:**
| Key | Description |
|-----|-------------|
| `grade_bg` | Background color for grade display |
| `grade_color` | Text color for grade letter |
| `grade_text` | Alternative text color |
| `rec_bg` | Recommendation background |
| `rec_border` | Recommendation border |
| `rec_color` | Recommendation text |

**Color Schemes:**
| Grade | Background | Text |
|-------|------------|------|
| A | Green (#d4edda) | Dark green (#155724) |
| B | Blue (#cce5ff) | Dark blue (#004085) |
| C | Yellow (#fff3cd) | Dark yellow (#856404) |
| D | Red (#f8d7da) | Dark red (#721c24) |
| F | Dark red (#721c24) | White (#fff) |

---

#### `get_rating_class`

```python
def get_rating_class(rating: str) -> str
```

Get CSS class for a rating level.

**Mappings:**
| Rating | CSS Class |
|--------|-----------|
| Excellent | `rating-excellent` |
| Good | `rating-good` |
| Acceptable | `rating-acceptable` |
| Poor | `rating-poor` |
| Critical | `rating-critical` |

---

#### `get_progress_class`

```python
def get_progress_class(rating: str) -> str
```

Get CSS class for progress bar styling.

**Mappings:**
| Rating | CSS Class |
|--------|-----------|
| Excellent | `progress-excellent` |
| Good | `progress-good` |
| Acceptable | `progress-acceptable` |
| Poor | `progress-poor` |
| Critical | `progress-poor` |

---

## neoswga.core.report.executive_summary

### Functions

#### `generate_executive_summary`

```python
def generate_executive_summary(
    results_dir: str,
    output_file: Optional[str] = None
) -> ExecutiveSummary
```

Generate an executive summary report.

**Parameters:**
| Name | Type | Default | Description |
|------|------|---------|-------------|
| `results_dir` | `str` | - | Path to results directory |
| `output_file` | `str` | `None` | Path to write HTML output |

**Returns:**
| Type | Description |
|------|-------------|
| `ExecutiveSummary` | Summary data with metrics and quality |

---

#### `create_executive_summary`

```python
def create_executive_summary(
    metrics: PipelineMetrics,
    quality: QualityAssessment
) -> ExecutiveSummary
```

Create executive summary from pre-collected data.

---

#### `render_executive_summary`

```python
def render_executive_summary(summary: ExecutiveSummary) -> str
```

Render executive summary to HTML string.

---

### Classes

#### `ExecutiveSummary`

```python
@dataclass
class ExecutiveSummary:
    metrics: PipelineMetrics
    quality: QualityAssessment
    generated_at: str
    version: str = "3.0.0"
```

---

## neoswga.core.report.technical_report

### Functions

#### `generate_technical_report`

```python
def generate_technical_report(
    results_dir: str,
    output_file: Optional[str] = None
) -> TechnicalReportData
```

Generate a comprehensive technical report.

**Parameters:**
| Name | Type | Default | Description |
|------|------|---------|-------------|
| `results_dir` | `str` | - | Path to results directory |
| `output_file` | `str` | `None` | Path to write HTML output |

**Returns:**
| Type | Description |
|------|-------------|
| `TechnicalReportData` | Complete report data |

---

### Classes

#### `TechnicalReportData`

```python
@dataclass
class TechnicalReportData:
    generated_at: str
    version: str = "3.0.0"
    metrics: Optional[PipelineMetrics] = None
    quality: Optional[QualityAssessment] = None
    total_runtime: float = 0.0
    step_runtimes: Dict[str, float] = field(default_factory=dict)
    filtering_stages: List[tuple] = field(default_factory=list)
    coverage_by_region: Dict[str, float] = field(default_factory=dict)
    gaps: List[GapInfo] = field(default_factory=list)
    binding_distribution: Dict[str, int] = field(default_factory=dict)
    primer_profiles: List[PrimerProfile] = field(default_factory=list)
    interactions: List[InteractionPair] = field(default_factory=list)
    max_interaction_dg: float = 0.0
```

---

#### `PrimerProfile`

```python
@dataclass
class PrimerProfile:
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
    unique_coverage: float
    contribution_score: float
    quality_rank: int
```

---

#### `InteractionPair`

```python
@dataclass
class InteractionPair:
    primer1: str
    primer2: str
    delta_g: float
    risk_level: str  # "low", "moderate", "high"
```

---

#### `GapInfo`

```python
@dataclass
class GapInfo:
    start: int
    end: int
    length: int
    chromosome: str = ""
    severity: str = "low"  # "low", "medium", "high", "critical"
    description: str = ""
```
