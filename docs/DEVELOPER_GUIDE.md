# NeoSWGA Developer Guide

Technical documentation for developers contributing to or extending NeoSWGA.

## Table of Contents

1. [Development Setup](#development-setup)
2. [Code Architecture](#code-architecture)
3. [Adding New Optimizers](#adding-new-optimizers)
4. [Extending Thermodynamics](#extending-thermodynamics)
5. [Testing](#testing)
6. [Code Style](#code-style)
7. [Contributing](#contributing)

---

## Development Setup

### Prerequisites

- Python 3.11+
- Git
- Jellyfish k-mer counter

### Environment Setup

```bash
# Clone repository
git clone https://github.com/yourusername/neoswga.git
cd neoswga

# Create virtual environment
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# Install in development mode with all dependencies
pip install -e ".[dev]"

# Install pre-commit hooks (optional)
pre-commit install
```

### Verify Installation

```bash
# Run tests
pytest tests/

# Run quick validation
neoswga validate --quick

# Check code style
black --check neoswga/
flake8 neoswga/
mypy neoswga/
```

---

## Code Architecture

### Package Structure

```
neoswga/
    __init__.py              # Package init, version
    cli_unified.py           # Main CLI entry point
    core/                    # Core functionality
        __init__.py
        # Pipeline
        pipeline.py
        improved_pipeline.py
        # Filtering
        filter.py
        adaptive_filters.py
        # Scoring
        rf_preprocessing.py
        # Optimization
        base_optimizer.py
        optimizer_factory.py
        greedy_optimizer.py
        network_optimizer.py
        genetic_algorithm.py
        # Thermodynamics
        thermodynamics.py
        reaction_conditions.py
        # Performance
        position_cache.py
        # Utilities
        utility.py
        parameter.py
        models/
            random_forest_filter.p
```

### Key Design Patterns

#### Factory Pattern (Optimizers)

```python
# Registration via decorator
@OptimizerRegistry.register('my-optimizer', aliases=['my', 'mo'])
class MyOptimizer(BaseOptimizer):
    """My custom optimizer."""

    def optimize(self, candidates, target_size, **kwargs):
        # Implementation
        pass

# Usage
optimizer = OptimizerFactory.create('my-optimizer', cache=cache, ...)
result = optimizer.optimize(candidates, target_size=6)
```

#### Dataclass Configuration

```python
from dataclasses import dataclass

@dataclass
class MyConfig:
    """Configuration for my component."""
    param1: int = 10
    param2: float = 0.5
    param3: str = 'default'

# Usage with defaults
config = MyConfig()
config = MyConfig(param1=20)
```

#### LRU Caching for Performance

```python
from functools import lru_cache

@lru_cache(maxsize=1_000_000)
def expensive_calculation(seq: str, temp: float) -> float:
    """Cached thermodynamic calculation."""
    # Expensive computation
    return result
```

---

## Adding New Optimizers

### Step 1: Create Optimizer Class

Create a new file `neoswga/core/my_optimizer.py`:

```python
"""
My custom optimizer for primer set selection.

Implements [describe algorithm].
"""

from typing import List, Optional
from dataclasses import dataclass
import numpy as np

from .base_optimizer import BaseOptimizer, OptimizationResult, PrimerSetMetrics
from .optimizer_factory import OptimizerRegistry
from .position_cache import PositionCache


@dataclass
class MyOptimizerConfig:
    """Configuration for MyOptimizer."""
    param1: int = 10
    param2: float = 0.5


@OptimizerRegistry.register(
    'my-optimizer',
    aliases=['my', 'mo'],
    description='My custom optimization algorithm'
)
class MyOptimizer(BaseOptimizer):
    """
    My custom optimizer.

    Algorithm:
    1. First step
    2. Second step
    3. Third step
    """

    def __init__(
        self,
        cache: PositionCache,
        fg_prefixes: List[str],
        fg_seq_lengths: List[int],
        config: Optional[MyOptimizerConfig] = None,
        **kwargs
    ):
        """
        Initialize optimizer.

        Args:
            cache: Position cache for fast lookups
            fg_prefixes: Foreground file prefixes
            fg_seq_lengths: Foreground genome lengths
            config: Optimizer configuration
        """
        super().__init__(cache, fg_prefixes, fg_seq_lengths, **kwargs)
        self.config = config or MyOptimizerConfig()

    def optimize(
        self,
        candidates: List[str],
        target_size: int,
        **kwargs
    ) -> OptimizationResult:
        """
        Find optimal primer set.

        Args:
            candidates: List of candidate primers
            target_size: Desired number of primers

        Returns:
            OptimizationResult with selected primers and metrics
        """
        # Your optimization logic here
        selected = self._my_algorithm(candidates, target_size)

        # Calculate metrics
        metrics = self._calculate_metrics(selected)

        return OptimizationResult(
            primer_set=selected,
            metrics=metrics,
            status='success',
            message=f'Selected {len(selected)} primers'
        )

    def _my_algorithm(
        self,
        candidates: List[str],
        target_size: int
    ) -> List[str]:
        """
        Core optimization algorithm.

        Args:
            candidates: Available primers
            target_size: Number to select

        Returns:
            Selected primer list
        """
        selected = []

        # Implementation
        for primer in candidates[:target_size]:
            selected.append(primer)

        return selected
```

### Step 2: Register with the Factory

Add the import to `neoswga/core/unified_optimizer.py` in `_ensure_optimizers_registered()`:

```python
from . import my_optimizer  # My custom optimization
```

This ensures the `@OptimizerFactory.register()` decorator fires when optimizers are first used.

### Step 3: Write Tests

Create `tests/test_my_optimizer.py`:

```python
"""Tests for MyOptimizer."""

import pytest
import numpy as np
from unittest.mock import MagicMock

from neoswga.core.my_optimizer import MyOptimizer, MyOptimizerConfig
from neoswga.core.optimizer_factory import OptimizerFactory


class TestMyOptimizer:
    """Test suite for MyOptimizer."""

    @pytest.fixture
    def mock_cache(self):
        """Create mock position cache."""
        cache = MagicMock()
        cache.get_positions.return_value = np.array([100, 200, 300])
        return cache

    @pytest.fixture
    def optimizer(self, mock_cache):
        """Create optimizer instance."""
        return MyOptimizer(
            cache=mock_cache,
            fg_prefixes=['test'],
            fg_seq_lengths=[10000]
        )

    def test_optimize_returns_correct_size(self, optimizer):
        """Test that optimizer returns correct number of primers."""
        candidates = ['ATCGATCG', 'GCTAGCTA', 'TGCATGCA', 'AATTCCGG']
        result = optimizer.optimize(candidates, target_size=2)

        assert len(result.primer_set) == 2

    def test_optimize_returns_valid_primers(self, optimizer):
        """Test that returned primers are from candidates."""
        candidates = ['ATCGATCG', 'GCTAGCTA', 'TGCATGCA']
        result = optimizer.optimize(candidates, target_size=2)

        for primer in result.primer_set:
            assert primer in candidates

    def test_factory_registration(self, mock_cache):
        """Test that optimizer is registered with factory."""
        optimizer = OptimizerFactory.create(
            'my-optimizer',
            cache=mock_cache,
            fg_prefixes=['test'],
            fg_seq_lengths=[10000]
        )

        assert isinstance(optimizer, MyOptimizer)

    def test_alias_registration(self, mock_cache):
        """Test that aliases work."""
        optimizer = OptimizerFactory.create(
            'my',
            cache=mock_cache,
            fg_prefixes=['test'],
            fg_seq_lengths=[10000]
        )

        assert isinstance(optimizer, MyOptimizer)
```

### Step 4: Update Documentation

Add to `docs/development/optimizers.md`:

```markdown
### my-optimizer

**Aliases:** `my`, `mo`

**Description:** My custom optimization algorithm.

**When to use:** [Describe use cases]

**Configuration:**
- `param1`: [Description] (default: 10)
- `param2`: [Description] (default: 0.5)

**Example:**
```bash
neoswga optimize -j params.json --optimization-method=my-optimizer
```
```

---

## Extending Thermodynamics

### Adding New Additive

Edit `neoswga/core/reaction_conditions.py`:

```python
# Add to AdditiveConcentrations class or ReactionConditions

# In ReactionConditions.__init__
self.my_additive_m: float = my_additive_m  # Molar concentration

# Add effect to calculate_effective_tm method
def calculate_effective_tm(self, seq: str) -> float:
    tm = base_tm

    # Existing corrections
    tm -= 0.6 * self.dmso_percent
    tm -= 2.3 * self.betaine_m

    # Add new additive effect
    # Effect based on literature: [citation]
    tm -= MY_ADDITIVE_EFFECT_PER_M * self.my_additive_m

    return tm
```

### Adding New Polymerase

Edit `neoswga/core/reaction_conditions.py`:

```python
POLYMERASE_CHARACTERISTICS = {
    # Existing entries...

    'my_polymerase': {
        'name': 'My Polymerase',
        'temp_range': (25.0, 40.0),
        'optimal_temp': 37.0,
        'processivity': 50000,  # bp
        'strand_displacement': True,
        'exonuclease': '3to5',
        'error_rate': 1e-5,
        'requires_primer': True,
        'description': 'My custom polymerase for specific applications'
    }
}
```

---

## Testing

### Running Tests

```bash
# All tests
pytest tests/

# Specific test file
pytest tests/test_thermodynamics_edge_cases.py

# With coverage
pytest --cov=neoswga --cov-report=html tests/

# Verbose output
pytest -v tests/

# Run specific test
pytest tests/test_filter.py::test_gc_filter
```

### Test Organization

```
tests/
    conftest.py              # Shared fixtures
    test_filter.py           # Filter module tests
    test_thermodynamics_edge_cases.py
    test_genetic_algorithm_integration.py
    test_position_cache.py
    integration/             # Integration tests
        phi29_baseline/
        equiphi29_baseline/
```

### Writing Tests

```python
"""Example test module."""

import pytest
import numpy as np
from neoswga.core import thermodynamics as thermo


class TestTmCalculation:
    """Tests for Tm calculation."""

    @pytest.fixture
    def standard_conditions(self):
        """Standard reaction conditions."""
        return {'na_conc': 50.0, 'mg_conc': 0.0}

    def test_basic_tm(self, standard_conditions):
        """Test basic Tm calculation."""
        tm = thermo.calculate_tm_with_salt('ATCGATCG', **standard_conditions)
        assert 20 < tm < 40

    def test_gc_content_affects_tm(self, standard_conditions):
        """Higher GC should give higher Tm."""
        low_gc = thermo.calculate_tm_with_salt('AAAATTTT', **standard_conditions)
        high_gc = thermo.calculate_tm_with_salt('GGGGCCCC', **standard_conditions)

        assert high_gc > low_gc

    @pytest.mark.parametrize('seq,expected_gc', [
        ('ATCG', 0.5),
        ('AAAA', 0.0),
        ('GGCC', 1.0),
    ])
    def test_gc_content(self, seq, expected_gc):
        """Test GC content calculation."""
        assert thermo.gc_content(seq) == pytest.approx(expected_gc)

    def test_invalid_sequence_raises(self):
        """Invalid bases should raise error."""
        with pytest.raises(ValueError):
            thermo.calculate_tm_with_salt('ATCGX')
```

### Integration Tests

```python
"""Integration test example."""

import pytest
import os
import tempfile
from pathlib import Path


@pytest.mark.integration
class TestPipelineIntegration:
    """End-to-end pipeline tests."""

    @pytest.fixture
    def test_data_dir(self):
        """Path to test data."""
        return Path(__file__).parent / 'integration' / 'phi29_baseline'

    def test_full_pipeline(self, test_data_dir):
        """Run complete pipeline on test data."""
        # Setup
        params_file = test_data_dir / 'params.json'

        # Run pipeline steps
        # (Implementation depends on how you want to test)

        # Verify output
        output_file = test_data_dir / 'results' / 'step4_improved_df.csv'
        assert output_file.exists()
```

---

## Code Style

### Python Style

We follow PEP 8 with these tools:
- **black**: Code formatting
- **isort**: Import sorting
- **flake8**: Linting
- **mypy**: Type checking

```bash
# Format code
black neoswga/ tests/
isort neoswga/ tests/

# Check style
flake8 neoswga/ tests/
mypy neoswga/
```

### Configuration

See `pyproject.toml`:

```toml
[tool.black]
line-length = 100
target-version = ['py311', 'py312', 'py313']

[tool.isort]
profile = "black"
line_length = 100

[tool.mypy]
python_version = "3.11"
warn_return_any = true
ignore_missing_imports = true
```

### Docstring Style

Use Google-style docstrings:

```python
def calculate_tm_with_salt(
    seq: str,
    na_conc: float = 50.0,
    mg_conc: float = 0.0
) -> float:
    """
    Calculate melting temperature with salt correction.

    Uses the SantaLucia nearest-neighbor model with Owczarzy
    unified salt correction formula.

    Args:
        seq: DNA sequence (5' to 3')
        na_conc: Na+ concentration in mM
        mg_conc: Mg2+ concentration in mM

    Returns:
        Melting temperature in Celsius

    Raises:
        ValueError: If sequence contains invalid bases

    Example:
        >>> calculate_tm_with_salt('ATCGATCG', na_conc=50)
        32.5

    References:
        - SantaLucia (1998) PNAS 95:1460-1465
        - Owczarzy et al. (2008) Biochemistry 47:5336-5353
    """
```

---

## Contributing

### Workflow

1. **Fork** the repository
2. **Create branch** for your feature:
   ```bash
   git checkout -b feature/my-feature
   ```
3. **Make changes** following code style
4. **Add tests** for new functionality
5. **Run tests** and style checks:
   ```bash
   pytest tests/
   black --check neoswga/
   flake8 neoswga/
   ```
6. **Commit** with clear message:
   ```bash
   git commit -m "Add feature: description of changes"
   ```
7. **Push** and create **Pull Request**

### Commit Messages

Use conventional commits:

```
feat: add new optimizer algorithm
fix: correct Tm calculation for palindromes
docs: update API reference
test: add integration tests for GA
refactor: simplify position cache loading
perf: optimize dimer matrix calculation
```

### Pull Request Checklist

- [ ] Tests pass locally
- [ ] New code has tests
- [ ] Documentation updated
- [ ] Code style checks pass
- [ ] Commit messages are clear
- [ ] PR description explains changes

### Code Review

PRs are reviewed for:
- Correctness
- Test coverage
- Code style
- Documentation
- Performance implications

---

## Development Tasks

### Retrain Random Forest Model

After sklearn updates:

```bash
python scripts/retrain_rf_model.py \
    --training-data data/training_set.csv \
    --output neoswga/core/models/random_forest_filter.p
```

### Generate K-mer Files for Testing

```bash
cd tests/integration/phi29_baseline
neoswga count-kmers -j params.json --min-k 6 --max-k 12
```

### Run Benchmarks

```bash
python scripts/benchmarking/run_benchmarks.py --output benchmarks/results/
```

### Update Documentation

```bash
# Build Sphinx docs
cd docs/
make html

# Preview
open _build/html/index.html
```

---

## See Also

- [API Reference](API_REFERENCE.md)
- [Architecture Diagrams](ARCHITECTURE_DIAGRAMS.md)
- [Module Reference](MODULE_REFERENCE.md)
- [Changelog](CHANGELOG.md)
