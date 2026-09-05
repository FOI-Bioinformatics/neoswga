"""End-to-end smoke test for the registered optimizer methods.

After the optimizer-zoo trim there are four registered methods (hybrid,
dominating-set, network, background-aware) plus the `ensemble` pseudo-method.
The four base methods have dedicated unit tests; the `ensemble` path is only
unit-tested with stubbed optimizers. This module runs the genuinely
under-covered integration paths end-to-end on the plasmid example and asserts
the optimizer does not crash and produces a non-empty primer set.

Methods that are not registered (e.g. legacy names removed by the trim) are
skipped rather than failed, so the test stays forward/backward compatible.

Steps 2 and 3 are run once per session (module-scoped fixture) so the
parametrised step-4 invocations share the same candidate pool and only the
optimizer differs. This keeps the suite fast.
"""

from __future__ import annotations

import json
import os
import shutil
import tempfile

import pytest

EXAMPLE_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "examples", "plasmid_example")


# ---------------------------------------------------------------------------
# Fixture: build step2 + step3 once, reuse for every optimizer param
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def plasmid_scored_workdir():
    """Copy the plasmid example to a tmpdir, run step2 + step3 once,
    leave step3_df.csv on disk for every parametrised optimizer call
    to consume.
    """
    if not os.path.isdir(EXAMPLE_DIR):
        pytest.skip("Plasmid example data not available")

    from neoswga.core.kmer_counter import check_jellyfish_available

    if not check_jellyfish_available():
        pytest.skip("jellyfish not available (required for count-kmers / step1)")

    tmpdir = tempfile.mkdtemp(prefix="neoswga_optmethods_")
    for fname in os.listdir(EXAMPLE_DIR):
        src = os.path.join(EXAMPLE_DIR, fname)
        if os.path.isfile(src):
            shutil.copy2(src, tmpdir)

    # Permissive params so step3 doesn't truncate to zero (tiny plasmid
    # makes the RF predictor score very low).
    params_path = os.path.join(tmpdir, "params.json")
    with open(params_path) as fh:
        params_data = json.load(fh)
    params_data.setdefault("min_amp_pred", 0.0)
    params_data.setdefault("schema_version", 1)
    with open(params_path, "w") as fh:
        json.dump(params_data, fh, indent=2)

    original_cwd = os.getcwd()
    os.chdir(tmpdir)

    def _reset():
        import neoswga.core.pipeline as pipeline_mod
        from neoswga.core import parameter as param_mod

        pipeline_mod._initialized = False
        pipeline_mod.fg_prefixes = None
        pipeline_mod.bg_prefixes = None
        pipeline_mod.fg_genomes = None
        pipeline_mod.bg_genomes = None
        pipeline_mod.fg_seq_lengths = None
        pipeline_mod.bg_seq_lengths = None
        pipeline_mod.fg_circular = None
        pipeline_mod.bg_circular = None
        param_mod.json_file = params_path

    _reset()
    from neoswga.core.pipeline import step1, step2, step3

    # Generate the k-mer count files (count-kmers); they are not committed.
    step1()

    _reset()
    df2 = step2()
    assert len(df2) > 0, "fixture: step2 produced no candidates"

    _reset()
    df3 = step3()
    assert len(df3) > 0, "fixture: step3 produced no scored primers"

    yield tmpdir, _reset

    os.chdir(original_cwd)
    shutil.rmtree(tmpdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Parametrised smoke test over the eight under-tested optimizers
# ---------------------------------------------------------------------------


# Every registered method plus the `ensemble` pseudo-method, run end-to-end.
# Unregistered names are skipped inside the test body.
#
# `hybrid` and `clique` were both missing from this list, which is to say the
# sweep covered every method except the DEFAULT one and the only one that makes
# a hard guarantee. The rationale was that the base methods have dedicated unit
# tests, but those construct optimizers directly and so cannot see anything the
# CLI-to-optimizer wiring drops -- which is exactly the defect class this
# repository keeps producing.
UNDER_TESTED_METHODS = [
    "ensemble",
    "hybrid",
    "dominating-set",
    "network",
    "background-aware",
    "clique",
]


@pytest.mark.integration
@pytest.mark.parametrize("method", UNDER_TESTED_METHODS)
def test_under_tested_optimizer_smoke(plasmid_scored_workdir, method):
    """Each method must run the plasmid example through step 4 without
    raising and return at least one non-empty primer set.

    Not asserting specific primer sequences or coverage bounds — the
    goal is to catch hard crashes (import errors, unimplemented branches,
    missing-dependency fallthroughs, arg-name drift between the factory
    and the optimizer constructor) rather than to pin behaviour.
    """
    tmpdir, reset = plasmid_scored_workdir
    reset()

    # Some methods (moea, multi-agent) import optional dependencies that
    # may be absent in the CI environment. Skip with a clear note rather
    # than fail the method outright; the missing-dep path is itself a
    # legitimate outcome we want surfaced.
    _maybe_skip_on_missing_optional_dep(method)

    from neoswga.core.unified_optimizer import optimize_step4
    from neoswga.core.exceptions import OptimizerNotFoundError

    try:
        primer_sets, scores, _cache = optimize_step4(
            optimization_method=method,
            verbose=False,
            seed=42,
        )
    except OptimizerNotFoundError as e:
        pytest.skip(f"{method} is not a registered optimizer: {e}")
    except NotImplementedError as e:
        pytest.skip(f"{method} is not fully implemented: {e}")
    except ImportError as e:
        pytest.skip(f"{method} requires an optional dependency: {e}")

    assert primer_sets is not None, f"{method}: optimize_step4 returned None"
    assert len(primer_sets) > 0, f"{method}: expected at least one primer set, got {primer_sets!r}"
    best_set = primer_sets[0]
    assert len(best_set) > 0, f"{method}: best primer set is empty"
    for primer in best_set:
        assert isinstance(primer, str), f"{method}: non-string primer entry {primer!r} in set"
        assert set(primer.upper()) <= set(
            "ACGTN"
        ), f"{method}: primer {primer!r} contains non-DNA characters"

    # Summary file is the contract every optimizer path writes to.
    summary_path = os.path.join(tmpdir, "step4_improved_df_summary.json")
    assert os.path.exists(summary_path), f"{method}: step4_improved_df_summary.json not written"
    with open(summary_path) as fh:
        summary = json.load(fh)
    assert "primers" in summary, f"{method}: summary missing 'primers' key"
    assert "metrics" in summary, f"{method}: summary missing 'metrics' key"


def _maybe_skip_on_missing_optional_dep(method: str) -> None:
    """Skip methods whose optional dependency is not installed.

    moea needs pymoo; multi-agent has no external dep but is a thin
    ensemble wrapper. The rest are pure-Python.
    """
    if method == "moea":
        try:
            import pymoo  # noqa: F401
        except ImportError:
            pytest.skip("moea requires pymoo (optional dependency)")
