"""
Shared fixtures for report module tests.

Uses importlib to avoid loading the full neoswga package.
"""

import csv
import json
import sys
import pytest
from pathlib import Path
from typing import Dict, List
import importlib.util

# Add the project root to sys.path for direct imports
project_root = Path(__file__).parent.parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

# Block the main package from loading to avoid dependency issues
# We only need the report submodule which has no external dependencies
sys.modules.setdefault('neoswga', type(sys)('neoswga'))
sys.modules.setdefault('neoswga.core', type(sys)('neoswga.core'))
sys.modules['neoswga'].__path__ = [str(project_root / 'neoswga')]
sys.modules['neoswga.core'].__path__ = [str(project_root / 'neoswga' / 'core')]

# Now we can import the report module directly
from neoswga.core.report.metrics import (
    PipelineMetrics,
    PrimerMetrics,
    GenomeInfo,
    CoverageMetrics,
    SpecificityMetrics,
    ThermodynamicMetrics,
    UniformityMetrics,
    FilteringStats,
    collect_pipeline_metrics,
)


@pytest.fixture
def fixtures_dir():
    """Path to test fixtures directory."""
    return Path(__file__).parent / "fixtures"


@pytest.fixture
def bacillus_results_dir(fixtures_dir):
    """Path to Bacillus test results directory."""
    return fixtures_dir / "bacillus_results"


@pytest.fixture
def minimal_step4_csv(tmp_path):
    """Create a minimal step4_improved_df.csv for testing."""
    csv_path = tmp_path / "step4_improved_df.csv"
    csv_path.write_text(
        "sequence,score,fg_freq,bg_freq,tm,gini,gc,fg_count,bg_count\n"
        "ATCGATCG,0.75,0.001,0.00001,35.5,0.25,0.5,100,1\n"
        "GCTAGCTA,0.65,0.0008,0.00002,34.2,0.30,0.5,80,2\n"
    )
    return tmp_path


@pytest.fixture
def complete_results_dir(tmp_path):
    """Create a complete results directory with all files."""
    # Create step4_improved_df.csv
    step4 = tmp_path / "step4_improved_df.csv"
    step4.write_text(
        "sequence,score,fg_freq,bg_freq,tm,gini,gc,fg_count,bg_count,"
        "dimer_score,hairpin_dg,self_dimer_dg,three_prime_stability,strand_ratio\n"
        "ATCGATCGATCG,0.85,0.002,0.00001,42.5,0.20,0.5,200,2,-3.5,-1.2,-2.8,0.8,1.2\n"
        "GCTAGCTAGCTA,0.80,0.0018,0.00001,41.8,0.22,0.5,180,2,-2.8,-0.9,-2.2,0.7,0.9\n"
        "TAGCTAGCTAGC,0.78,0.0015,0.00002,40.2,0.25,0.5,150,3,-4.0,-1.5,-3.1,0.9,1.1\n"
        "CGATCGATCGAT,0.72,0.0012,0.00001,39.5,0.28,0.5,120,1,-3.2,-1.0,-2.5,0.6,1.0\n"
        "ATGCATGCATGC,0.70,0.0010,0.00002,38.8,0.30,0.5,100,2,-2.5,-0.8,-2.0,0.5,0.8\n"
        "CATGCATGCATG,0.68,0.0009,0.00001,38.2,0.32,0.5,90,1,-2.2,-0.7,-1.8,0.4,1.3\n"
    )

    # Create params.json
    params = tmp_path / "params.json"
    params.write_text(json.dumps({
        "fg": "/path/to/target.fna",
        "bg": "/path/to/background.fna",
        "fg_size": 4000000,
        "bg_size": 3000000000,
        "min_k": 10,
        "max_k": 12,
        "polymerase": "phi29",
        "reaction_temp": 30.0,
        "num_primers": 6,
    }))

    # Create filter_stats.json
    filter_stats = tmp_path / "filter_stats.json"
    filter_stats.write_text(json.dumps({
        "total_kmers": 100000,
        "after_frequency": 50000,
        "after_background": 10000,
        "after_gini": 5000,
        "after_thermodynamic": 1000,
        "after_complexity": 500,
        "final_candidates": 6,
    }))

    return tmp_path


@pytest.fixture
def sample_primer_metrics():
    """Create sample PrimerMetrics for testing."""
    return [
        PrimerMetrics(
            sequence="ATCGATCGATCG",
            length=12,
            gc_content=0.5,
            tm=42.5,
            fg_freq=0.002,
            bg_freq=0.00001,
            fg_sites=200,
            bg_sites=2,
            gini=0.20,
            specificity=200000.0,
            amp_pred=0.85,
            dimer_score=-3.5,
            hairpin_dg=-1.2,
            self_dimer_dg=-2.8,
            three_prime_stability=0.8,
            strand_ratio=1.2,
        ),
        PrimerMetrics(
            sequence="GCTAGCTAGCTA",
            length=12,
            gc_content=0.5,
            tm=41.8,
            fg_freq=0.0018,
            bg_freq=0.00001,
            fg_sites=180,
            bg_sites=2,
            gini=0.22,
            specificity=180000.0,
            amp_pred=0.80,
            dimer_score=-2.8,
            hairpin_dg=-0.9,
            self_dimer_dg=-2.2,
            three_prime_stability=0.7,
            strand_ratio=0.9,
        ),
    ]


@pytest.fixture
def sample_pipeline_metrics(sample_primer_metrics):
    """Create sample PipelineMetrics for testing."""
    return PipelineMetrics(
        results_dir="/test/results",
        generated_at="2024-01-01T00:00:00",
        pipeline_version="3.6.0",
        target_genome=GenomeInfo(
            name="target",
            size=4000000,
            gc_content=0.45,
        ),
        background_genome=GenomeInfo(
            name="background",
            size=3000000000,
            gc_content=0.40,
        ),
        parameters={
            "fg_size": 4000000,
            "bg_size": 3000000000,
            "polymerase": "phi29",
            "reaction_temp": 30.0,
        },
        primers=sample_primer_metrics,
        primer_count=len(sample_primer_metrics),
        coverage=CoverageMetrics(
            overall_coverage=0.85,
            covered_bases=3400000,
            total_bases=4000000,
        ),
        specificity=SpecificityMetrics(
            enrichment_ratio=500.0,
            target_sites=380,
            background_sites=4,
            target_density=95.0,
            background_density=0.0013,
        ),
        thermodynamics=ThermodynamicMetrics(
            mean_tm=42.15,
            min_tm=41.8,
            max_tm=42.5,
            tm_range=0.7,
            reaction_temp=30.0,
            polymerase="phi29",
            max_heterodimer_dg=-3.5,
            dimer_risk_level="low",
        ),
        uniformity=UniformityMetrics(
            mean_gini=0.21,
            max_gini=0.22,
            strand_ratio=1.05,
        ),
        filtering=FilteringStats(
            total_kmers=100000,
            after_frequency=50000,
            after_background=10000,
            after_gini=5000,
            after_thermodynamic=1000,
            after_complexity=500,
            final_candidates=2,
        ),
    )


@pytest.fixture
def minimal_pipeline_metrics():
    """Create minimal/empty PipelineMetrics for edge case testing."""
    return PipelineMetrics(
        results_dir="/empty/results",
        generated_at="2024-01-01T00:00:00",
        primers=[],
        primer_count=0,
    )


@pytest.fixture
def edge_case_primer_row():
    """CSV row with edge case values for testing."""
    return {
        "sequence": "atcgatcg",  # lowercase
        "fg_freq": "nan",
        "bg_freq": "",
        "tm": "inf",
        "gini": "-0.5",  # invalid negative
        "gc": "",  # missing, should calculate
        "fg_count": "abc",  # invalid
        "bg_count": "12.5",  # float string
        "dimer_score": "None",
        "strand_ratio": "null",
    }


def create_csv_file(path: Path, rows: List[Dict]) -> Path:
    """Helper to create a CSV file from a list of dictionaries."""
    if not rows:
        path.write_text("")
        return path

    fieldnames = list(rows[0].keys())
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path
