"""Phase 15A: unified_optimizer populates per_target_coverage on every
optimizer's result so multi-target runs can surface the "target B at
40%" case instead of a single aggregate.

Before this change the field existed (`PrimerSetMetrics.per_target_coverage
= field(default_factory=dict)`) but was never populated by any optimizer,
rendering the post-optimization validator's `per_target_coverage_below_
threshold` warning unreachable.
"""

import json
import os
import shutil
import tempfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
EXAMPLE_DIR = ROOT / "examples" / "plasmid_example"


def _reset_pipeline_state(params_file):
    import neoswga.core.pipeline as pipeline_mod
    from neoswga.core import parameter
    pipeline_mod._initialized = False
    pipeline_mod.fg_prefixes = None
    pipeline_mod.bg_prefixes = None
    pipeline_mod.fg_genomes = None
    pipeline_mod.bg_genomes = None
    pipeline_mod.fg_seq_lengths = None
    pipeline_mod.bg_seq_lengths = None
    pipeline_mod.fg_circular = None
    pipeline_mod.bg_circular = None
    parameter.json_file = params_file


@pytest.fixture
def primed_workdir(tmp_path):
    if not EXAMPLE_DIR.is_dir():
        pytest.skip("plasmid example not available")
    for fname in os.listdir(EXAMPLE_DIR):
        src = EXAMPLE_DIR / fname
        if src.is_file():
            shutil.copy2(src, tmp_path / fname)
    params_path = tmp_path / "params.json"
    with open(params_path) as fh:
        params = json.load(fh)
    params.setdefault("schema_version", 1)
    params.setdefault("min_amp_pred", 0.0)
    params["min_k"] = 8
    params["max_k"] = 10
    with open(params_path, "w") as fh:
        json.dump(params, fh, indent=2)
    os.chdir(tmp_path)
    _reset_pipeline_state(str(params_path))
    from neoswga.core.pipeline import step2, step3
    step2()
    step3()
    return tmp_path


@pytest.mark.integration
@pytest.mark.slow
def test_hybrid_populates_per_target_coverage(primed_workdir):
    from neoswga.core import unified_optimizer as _uo
    _reset_pipeline_state(str(primed_workdir / "params.json"))
    _uo.optimize_step4(optimization_method="hybrid", verbose=False)

    # The wrapper stashes the most recent raw result for CLI consumers;
    # use it to inspect per_target_coverage.
    last = getattr(_uo, "_LAST_RESULT", None)
    assert last is not None, "unified_optimizer did not stash _LAST_RESULT"
    per_target = last.metrics.per_target_coverage
    assert per_target, (
        f"per_target_coverage should be populated after optimize; got {per_target!r}"
    )
    # Plasmid example has one fg prefix (pcDNA)
    assert "pcDNA" in {str(k).split("/")[-1] for k in per_target}, (
        f"expected pcDNA key in per_target_coverage; got {list(per_target)}"
    )
    for prefix, cov in per_target.items():
        assert 0.0 <= cov <= 1.0, (
            f"per_target coverage for {prefix} out of range: {cov}"
        )


@pytest.mark.integration
@pytest.mark.slow
def test_per_target_coverage_serialised_to_json(primed_workdir):
    """PrimerSetMetrics.to_dict must emit per_target_coverage so the CSV
    summary and downstream consumers can read it."""
    from neoswga.core import unified_optimizer as _uo
    _reset_pipeline_state(str(primed_workdir / "params.json"))
    _uo.optimize_step4(optimization_method="hybrid", verbose=False)
    last = _uo._LAST_RESULT
    metrics_dict = last.metrics.to_dict()
    assert "per_target_coverage" in metrics_dict
    # Must be JSON-serialisable
    json.dumps(metrics_dict["per_target_coverage"])


def test_compute_per_prefix_coverage_empty_input_is_safe():
    """Direct unit test: empty prefix list returns (0.0, {})."""
    from neoswga.core.coverage import compute_per_prefix_coverage

    class DummyCache:
        def get_positions(self, *a, **kw):
            return []

    agg, per = compute_per_prefix_coverage(
        cache=DummyCache(), primers=[], prefixes=[], seq_lengths=[],
    )
    assert agg == 0.0 and per == {}


def test_compute_per_prefix_coverage_marks_occupied_region(tmp_path):
    """Direct unit test: a single primer with one position marks 2*extension
    bases as covered, up to genome boundary."""
    from neoswga.core.coverage import compute_per_prefix_coverage

    class FakeCache:
        def get_positions(self, prefix, primer, strand):
            return [500]  # single position

    # 1000 bp genome, extension 200 -> covers positions [300, 700) = 400 bp
    agg, per = compute_per_prefix_coverage(
        cache=FakeCache(), primers=["ACGT"],
        prefixes=["x"], seq_lengths=[1000], extension=200,
    )
    assert abs(agg - 0.40) < 1e-9
    assert abs(per["x"] - 0.40) < 1e-9
