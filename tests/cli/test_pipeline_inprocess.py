"""In-process runs of the four core pipeline commands.

The other pipeline tests drive the CLI through a subprocess, which is realistic
but invisible to coverage and cannot inspect intermediate state. These run the
handlers directly so the argument merging, preset resolution and adaptive-GC
paths inside `neoswga/cli/pipeline.py` are actually exercised -- that module was
the largest untested surface in the package at 18%.

Each workspace is built fresh per test. The pipeline mutates module-level
parameter globals, so sharing one across tests would let an override from one
leak into the next and quietly invalidate the comparison.
"""

import json
import os
import shutil

import pytest


def _needs_jellyfish():
    if not shutil.which("jellyfish"):
        pytest.skip("jellyfish not available (required for count-kmers)")


@pytest.fixture
def workspace(tmp_path, genome_seq):
    """A minimal params.json plus genome, ready for count-kmers."""
    _needs_jellyfish()

    fasta = tmp_path / "target.fasta"
    fasta.write_text(
        ">target\n"
        + "\n".join(genome_seq[i : i + 70] for i in range(0, len(genome_seq), 70))
        + "\n"
    )

    def _params(**overrides):
        params = {
            "fg_genomes": [str(fasta)],
            "bg_genomes": [],
            "fg_prefixes": [str(tmp_path / "target")],
            "bg_prefixes": [],
            "data_dir": str(tmp_path / "results"),
            "min_k": 10,
            "max_k": 10,
            "polymerase": "phi29",
            "reaction_temp": 30.0,
            "min_fg_freq": 1e-6,
            "max_bg_freq": 1.0,
            "max_gini": 1.0,
            "max_primer": 40,
            "min_amp_pred": 0,
            "min_tm": 0,
            "max_tm": 100,
            "gc_min": 0.0,
            "gc_max": 1.0,
            "num_primers": 4,
            "target_set_size": 4,
            "max_sets": 2,
            "iterations": 2,
            "cpus": 1,
            "fg_circular": True,
            "schema_version": 2,
        }
        params.update(overrides)
        path = tmp_path / "params.json"
        path.write_text(json.dumps(params, indent=2))
        return path

    return {"root": tmp_path, "params": _params, "fasta": str(fasta)}


def _run(argv, cwd):
    """Dispatch a pipeline command in-process, from `cwd`."""
    from neoswga import cli_unified

    handlers = {
        "count-kmers": cli_unified.run_step1,
        "filter": cli_unified.run_step2,
        "score": cli_unified.run_step3,
        "optimize": cli_unified.run_step4,
    }
    args = cli_unified.create_parser().parse_args(argv)
    previous = os.getcwd()
    os.chdir(cwd)
    try:
        return handlers[args.command](args)
    finally:
        os.chdir(previous)


@pytest.fixture
def counted(workspace):
    """A workspace with k-mer counts and position files already built."""
    params_file = workspace["params"]()
    _run(["count-kmers", "-j", str(params_file)], workspace["root"])
    return {**workspace, "params_file": params_file}


# ----------------------------------------------------------------------
# count-kmers
# ----------------------------------------------------------------------


def test_count_kmers_writes_the_count_table(counted):
    files = os.listdir(counted["root"])
    assert any(f.endswith("mer_all.txt") for f in files), files


def test_position_index_appears_at_filter_not_count_kmers(counted):
    """Which step writes the HDF5 is a real contract, not an implementation detail.

    Anything reading positions -- PositionCache, and therefore every iterate and
    evaluate command -- needs `filter` to have run, not just `count-kmers`.
    """
    root = counted["root"]
    assert not any(f.endswith("_positions.h5") for f in os.listdir(root))

    _run(["filter", "-j", str(counted["params_file"])], root)
    assert any(f.endswith("_positions.h5") for f in os.listdir(root)), os.listdir(root)


def test_k_range_can_be_overridden_on_the_command_line(workspace):
    """`--min-k/--max-k` is the documented way to build non-default lengths."""
    params_file = workspace["params"](min_k=10, max_k=10)
    _run(
        ["count-kmers", "-j", str(params_file), "--min-k", "8", "--max-k", "8"],
        workspace["root"],
    )
    files = os.listdir(workspace["root"])
    assert any("8mer" in f for f in files), files


# ----------------------------------------------------------------------
# filter
# ----------------------------------------------------------------------


def test_filter_produces_candidates_and_a_funnel(counted):
    _run(["filter", "-j", str(counted["params_file"])], counted["root"])

    data_dir = counted["root"] / "results"
    assert (data_dir / "step2_df.csv").is_file()

    stats = json.loads((data_dir / "filter_stats.json").read_text())
    assert stats, "filter_stats.json is empty; the report renders these counts"


def test_filter_preset_flag_works(counted):
    """`filter --preset` raised NameError: load_preset_conditions was called but
    never imported, and no test passed the flag.

    It is the kind of break that survives a CLI split precisely because the
    parser still accepts the argument.
    """
    _run(
        ["filter", "-j", str(counted["params_file"]), "--preset", "standard_phi29"],
        counted["root"],
    )
    assert (counted["root"] / "results" / "step2_df.csv").is_file()


@pytest.mark.parametrize("preset", ["standard_phi29", "enhanced_equiphi29", "high_gc_genome"])
def test_every_displayed_preset_survives_a_real_filter_run(counted, preset):
    _run(["filter", "-j", str(counted["params_file"]), "--preset", preset], counted["root"])
    assert (counted["root"] / "results" / "step2_df.csv").is_file()


def test_gc_bounds_constrain_the_candidates(counted):
    """A narrow GC window must not keep more primers than a wide one."""
    import pandas as pd

    def survivors(gc_min, gc_max):
        _run(
            [
                "filter",
                "-j",
                str(counted["params_file"]),
                "--gc-min",
                str(gc_min),
                "--gc-max",
                str(gc_max),
            ],
            counted["root"],
        )
        return len(pd.read_csv(counted["root"] / "results" / "step2_df.csv"))

    narrow = survivors(0.45, 0.55)
    wide = survivors(0.0, 1.0)
    assert narrow <= wide, f"narrow GC window kept more primers ({narrow}) than wide ({wide})"


def test_tm_bounds_constrain_the_candidates(counted):
    import pandas as pd

    def survivors(min_tm, max_tm):
        _run(
            [
                "filter",
                "-j",
                str(counted["params_file"]),
                "--min-tm",
                str(min_tm),
                "--max-tm",
                str(max_tm),
            ],
            counted["root"],
        )
        return len(pd.read_csv(counted["root"] / "results" / "step2_df.csv"))

    assert survivors(28, 32) <= survivors(0, 100)


# ----------------------------------------------------------------------
# score
# ----------------------------------------------------------------------


@pytest.fixture
def filtered(counted):
    _run(["filter", "-j", str(counted["params_file"])], counted["root"])
    return counted


def test_score_carries_every_candidate_forward(filtered):
    """`score` prepares the candidate pool; it no longer scores it (finding F0).

    It used to assert that the stage ADDED columns. It now adds none: the
    amplification model is retired from the default path, and what the stage
    contributes is the deterministic order that makes an unseeded run
    reproducible. It must still drop nothing.
    """
    import pandas as pd

    _run(["score", "-j", str(filtered["params_file"])], filtered["root"])

    data_dir = filtered["root"] / "results"
    step2 = pd.read_csv(data_dir / "step2_df.csv")
    step3 = pd.read_csv(data_dir / "step3_df.csv")

    assert len(step3) == len(step2), f"score dropped candidates: {len(step2)} -> {len(step3)}"
    assert set(step3["primer"]) == set(step2["primer"])
    assert step3["gini"].is_monotonic_increasing, "the deterministic order was lost"


# ----------------------------------------------------------------------
# optimize
# ----------------------------------------------------------------------


@pytest.fixture
def scored(filtered):
    _run(["score", "-j", str(filtered["params_file"])], filtered["root"])
    return filtered


@pytest.mark.parametrize("method", ["hybrid", "dominating-set", "network"])
def test_each_optimization_method_produces_a_valid_summary(scored, method):
    """All methods must write the same well-formed summary the report reads."""
    _run(
        ["optimize", "-j", str(scored["params_file"]), "-m", method, "--seed", "1"],
        scored["root"],
    )

    def reject(token):
        raise ValueError(f"non-standard JSON token {token!r}")

    summary_path = scored["root"] / "results" / "step4_improved_df_summary.json"
    summary = json.loads(summary_path.read_text(), parse_constant=reject)

    assert summary.get("primers"), f"{method} selected no primers"
    assert summary["metrics"]["fg_coverage"] is not None


def test_seeded_optimization_is_reproducible(scored):
    def run_once():
        _run(
            ["optimize", "-j", str(scored["params_file"]), "-m", "hybrid", "--seed", "42"],
            scored["root"],
        )
        path = scored["root"] / "results" / "step4_improved_df_summary.json"
        return json.loads(path.read_text())["primers"]

    assert run_once() == run_once()


def test_application_profile_is_accepted_and_recorded(scored):
    """`--application clinical` weights selection differently."""
    _run(
        [
            "optimize",
            "-j",
            str(scored["params_file"]),
            "--application",
            "clinical",
            "--seed",
            "1",
        ],
        scored["root"],
    )
    summary_path = scored["root"] / "results" / "step4_improved_df_summary.json"
    assert json.loads(summary_path.read_text()).get("primers")


# ----------------------------------------------------------------------
# --enable-qa
# ----------------------------------------------------------------------
#
# The flag dispatched to `pipeline_qa_integration.run_step2_with_qa()` and
# `run_step3_with_qa()`, neither of which exists. The AttributeError was
# swallowed by each step's `except Exception` and turned into `sys.exit(1)`, so
# every step the flag is offered on failed the run instead of doing QA.


@pytest.fixture
def filtered_qa(counted):
    _run(["filter", "-j", str(counted["params_file"]), "--enable-qa"], counted["root"])
    return counted


def test_filter_with_enable_qa_does_not_crash(counted):
    import pandas as pd

    _run(["filter", "-j", str(counted["params_file"]), "--enable-qa"], counted["root"])

    step2 = pd.read_csv(counted["root"] / "results" / "step2_df.csv")
    assert "qa_score" in step2.columns, step2.columns.tolist()
    assert (counted["root"] / "results" / "qa_report.txt").is_file()


def test_qa_filtering_only_removes_candidates(counted):
    """QA is a filter: it may reject primers, never invent them."""
    import pandas as pd

    def survivors(*extra):
        _run(["filter", "-j", str(counted["params_file"]), *extra], counted["root"])
        return set(pd.read_csv(counted["root"] / "results" / "step2_df.csv")["primer"])

    plain = survivors()
    with_qa = survivors("--enable-qa")
    assert with_qa <= plain


def test_enable_qa_does_not_stick_to_the_next_run(counted):
    """`parameter.enable_qa` was only ever set True, never cleared.

    In one process -- `design`, the SDK, a test session -- that turned QA on
    for every later step whether or not the flag was given.
    """
    import pandas as pd

    _run(["filter", "-j", str(counted["params_file"]), "--enable-qa"], counted["root"])
    _run(["filter", "-j", str(counted["params_file"])], counted["root"])

    step2 = pd.read_csv(counted["root"] / "results" / "step2_df.csv")
    assert "qa_score" not in step2.columns, "QA ran without --enable-qa"


def test_score_with_enable_qa_combines_rf_and_qa_scores(filtered_qa):
    import pandas as pd

    _run(["score", "-j", str(filtered_qa["params_file"]), "--enable-qa"], filtered_qa["root"])

    step3 = pd.read_csv(filtered_qa["root"] / "results" / "step3_df.csv")
    assert "qa_score" in step3.columns, step3.columns.tolist()
    assert "composite_score" in step3.columns, step3.columns.tolist()
    assert step3["composite_score"].notna().all()
    assert step3["composite_score"].is_monotonic_decreasing, "pool is not ranked by composite score"


def test_score_with_enable_qa_works_without_a_qa_filtered_step2(filtered):
    """`score --enable-qa` after a plain `filter` must score, not crash."""
    import pandas as pd

    _run(["score", "-j", str(filtered["params_file"]), "--enable-qa"], filtered["root"])

    step3 = pd.read_csv(filtered["root"] / "results" / "step3_df.csv")
    assert "composite_score" in step3.columns, step3.columns.tolist()


def test_optimize_with_enable_qa_passes_a_qa_filtered_pool(scored, monkeypatch):
    """Step 4's QA hook is the candidate pre-filter, so the pool must reach it."""
    import pandas as pd

    from neoswga.core import unified_optimizer

    seen = {}
    real = unified_optimizer.optimize_step4

    def capture(**kwargs):
        seen["candidates"] = kwargs.get("candidates")
        return real(**kwargs)

    monkeypatch.setattr(unified_optimizer, "optimize_step4", capture)
    _run(
        ["optimize", "-j", str(scored["params_file"]), "--enable-qa", "--seed", "1"],
        scored["root"],
    )

    pool = set(pd.read_csv(scored["root"] / "results" / "step3_df.csv")["primer"])
    assert seen["candidates"] is not None, "QA pool was never built"
    assert set(seen["candidates"]) <= pool


def test_optimize_without_enable_qa_leaves_the_pool_alone(scored, monkeypatch):
    from neoswga.core import unified_optimizer

    seen = {}
    real = unified_optimizer.optimize_step4

    def capture(**kwargs):
        seen["candidates"] = kwargs.get("candidates")
        return real(**kwargs)

    monkeypatch.setattr(unified_optimizer, "optimize_step4", capture)
    _run(["optimize", "-j", str(scored["params_file"]), "--seed", "1"], scored["root"])

    assert seen["candidates"] is None
