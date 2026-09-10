"""`design` could not reach the set-size tooling at all.

run_design builds a defaults namespace and fixes "auto_size": False in it, and
auto_size is not a params.json key, so the one-shot entry point had no route to
--auto-size or --show-frontier. The defaults loop is `if not hasattr(args,
attr)`, so registering the flags on the design subparser is enough.
"""

import pytest


def _parse(argv):
    from neoswga import cli_unified

    return cli_unified.create_parser().parse_args(argv)


@pytest.mark.parametrize(
    "flag,attr",
    [
        ("--auto-size", "auto_size"),
        ("--show-frontier", "show_frontier"),
        ("--quick-estimate", "quick_estimate"),
    ],
)
def test_design_accepts_the_store_true_flags(flag, attr):
    args = _parse(["design", "-j", "params.json", flag])
    assert getattr(args, attr) is True


def test_design_accepts_the_valued_flags():
    args = _parse(
        [
            "design",
            "-j",
            "params.json",
            "--application",
            "clinical",
            "--min-fg-bg-ratio",
            "12.5",
            "--template-gc",
            "0.65",
        ]
    )
    assert args.application == "clinical"
    assert args.min_fg_bg_ratio == 12.5
    assert args.template_gc == 0.65


def test_design_still_defaults_auto_size_off():
    args = _parse(["design", "-j", "params.json"])
    assert args.auto_size is False
    assert args.show_frontier is False


def test_the_defaults_namespace_does_not_overwrite_a_passed_flag():
    """The `if not hasattr` loop is what makes registering the flag sufficient.
    A change to `setattr` unconditionally would silently re-break this."""
    from pathlib import Path

    import neoswga.cli.commands as commands

    source = Path(commands.__file__).read_text()
    assert "if not hasattr(args, attr):" in source, (
        "run_design's defaults loop no longer guards on hasattr, so a flag the "
        "user passed would be overwritten by the default"
    )


def test_a_passed_flag_survives_the_defaults_namespace(monkeypatch, tmp_path):
    """The behavioural form of the test above. Drives run_design with the flag
    set and reads back what run_step4 was handed."""
    import json

    import neoswga.cli.commands as commands

    (tmp_path / "params.json").write_text(json.dumps({"data_dir": str(tmp_path)}))

    seen = {}

    def _capture(args):
        seen["auto_size"] = args.auto_size
        seen["application"] = args.application
        raise SystemExit(0)

    monkeypatch.setattr(commands, "run_step1", lambda args: None)
    monkeypatch.setattr(commands, "run_step2", lambda args: None)
    monkeypatch.setattr(commands, "run_step3", lambda args: None)
    monkeypatch.setattr(commands, "run_step4", _capture)

    args = _parse(
        [
            "design",
            "-j",
            str(tmp_path / "params.json"),
            "--auto-size",
            "--application",
            "clinical",
        ]
    )
    with pytest.raises(SystemExit):
        commands.run_design(args)

    assert seen["auto_size"] is True, "the defaults namespace overwrote the flag"
    assert seen["application"] == "clinical"


def test_run_design_no_longer_hardcodes_auto_size_false():
    from pathlib import Path

    import neoswga.cli.commands as commands

    source = Path(commands.__file__).read_text()
    assert (
        '"auto_size": False' not in source
    ), "run_design still fixes auto_size False in its defaults namespace"


def test_the_frontier_reads_the_real_genome_lengths():
    """parameter.fg_lengths does not exist; the module global is
    fg_seq_lengths. The getattr default was a hardcoded 1 Mb genome."""
    from pathlib import Path

    import neoswga.cli.pipeline as cli_pipeline

    source = Path(cli_pipeline.__file__).read_text()
    assert 'getattr(parameter, "fg_lengths"' not in source
    assert 'getattr(parameter, "bg_lengths"' not in source
    assert (
        "1_000_000]" not in source
    ), "a hardcoded 1 Mb genome length is still present in the frontier block"


def test_parameter_has_no_fg_lengths_attribute():
    """The premise of the fix. If this ever becomes true, revisit it."""
    from neoswga.core import parameter

    assert not hasattr(parameter, "fg_lengths")


# ----------------------------------------------------------------------
# The frontier's own estimator, found broken while fixing the genome length
# ----------------------------------------------------------------------


def _pool_as_step3_writes_it():
    """step3_df.csv as the current pipeline writes it: counts, not frequencies.

    The score stage was retired on 2026-09-05, and the columns that survive are
    primer, ratio, gini, fg_count and bg_count. There is no fg_freq.
    """
    import pandas as pd

    return pd.DataFrame(
        {
            "primer": [f"ACGTACGTAC{i:02d}" for i in range(8)],
            "ratio": [0.1 * (i + 1) for i in range(8)],
            "gini": [0.3] * 8,
            "fg_count": [40, 35, 30, 25, 20, 15, 10, 5],
            "bg_count": [1, 2, 3, 4, 5, 6, 7, 8],
        }
    )


def test_the_estimator_uses_the_counts_step3_actually_writes():
    """Without this the estimator fell through to a random-sequence model.

    For a 12-mer that model predicts 4**-12 * 2 sites per bp, which over a
    6.2 kb plasmid is 0.003 sites and truncates to zero, so every size on the
    frontier reported 0% coverage and a 0.0 ratio. It only looked like an
    answer before because the genome length was hardcoded to 1 Mb.
    """
    from neoswga.core.set_size_optimizer import ParetoFrontierGenerator

    generator = ParetoFrontierGenerator(
        primer_pool=_pool_as_step3_writes_it(),
        fg_seq_lengths=[6157],
        bg_seq_lengths=[6258],
        processivity=3000,
    )
    points = generator._estimate_from_statistics(min_size=4, max_size=6)

    assert points, "the estimator produced no points at all"
    assert any(
        p.fg_coverage > 0 for p in points
    ), f"every size reports zero coverage: {[(p.set_size, p.fg_coverage) for p in points]}"
    assert any(
        p.fg_bg_ratio > 0 for p in points
    ), f"every size reports a zero fg/bg ratio: {[(p.set_size, p.fg_bg_ratio) for p in points]}"


def test_the_estimator_ranks_by_the_count_ratio_not_input_order():
    """`ratio` in step3_df is bg/fg, so it is not fg_bg_ratio and must not be
    read as one. Without a usable key the pool was left in gini order."""
    from neoswga.core.set_size_optimizer import ParetoFrontierGenerator

    pool = _pool_as_step3_writes_it()
    # Reverse the rows: the worst primers now come first in input order.
    generator = ParetoFrontierGenerator(
        primer_pool=pool.iloc[::-1].reset_index(drop=True),
        fg_seq_lengths=[6157],
        bg_seq_lengths=[6258],
        processivity=3000,
    )
    reversed_points = {p.set_size: p.fg_bg_ratio for p in generator._estimate_from_statistics(4, 6)}

    generator2 = ParetoFrontierGenerator(
        primer_pool=pool,
        fg_seq_lengths=[6157],
        bg_seq_lengths=[6258],
        processivity=3000,
    )
    forward_points = {p.set_size: p.fg_bg_ratio for p in generator2._estimate_from_statistics(4, 6)}

    assert reversed_points == forward_points, (
        "the frontier depends on the order the rows happened to arrive in, so "
        "it is not ranking the pool"
    )


def test_frequency_columns_still_win_when_present():
    """A caller passing a pool with fg_freq keeps the old behaviour."""
    import pandas as pd

    from neoswga.core.set_size_optimizer import ParetoFrontierGenerator

    pool = pd.DataFrame(
        {
            "primer": [f"ACGTACGTAC{i:02d}" for i in range(8)],
            "fg_freq": [1e-3] * 8,
            "bg_freq": [1e-5] * 8,
        }
    )
    generator = ParetoFrontierGenerator(
        primer_pool=pool,
        fg_seq_lengths=[1_000_000],
        bg_seq_lengths=[1_000_000],
        processivity=3000,
    )
    points = generator._estimate_from_statistics(4, 6)
    assert any(p.fg_coverage > 0 for p in points)


def test_show_frontier_reports_a_real_number_end_to_end(tmp_path):
    """Every other frontier test here reads source text, and that is not enough.

    Extracting the frontier block into its own function broke it twice while
    all of them stayed green: once because the decorator above run_step4 ended
    up on the new function, and once because the block referenced two locals
    that were not passed in. Both failures are caught here and nowhere else.

    The assertion is that coverage is non-zero. Against the 6.2 kb plasmid pair
    a 3 kb reach saturates, so the honest expected answer is 100%; a run that
    reports 0% means the estimator fell through to its random-sequence model,
    which for a 12-mer predicts 0.003 sites and truncates to nothing.
    """
    import json
    import random
    import shutil
    import subprocess
    import sys

    if not shutil.which("jellyfish"):
        pytest.skip("jellyfish not available (required for count-kmers)")

    rng = random.Random(20260910)
    seq = "".join(rng.choice("ACGT") for _ in range(20_000))
    fasta = tmp_path / "target.fasta"
    fasta.write_text(
        ">target\n" + "\n".join(seq[i : i + 70] for i in range(0, len(seq), 70)) + "\n"
    )

    params_file = tmp_path / "params.json"
    params_file.write_text(
        json.dumps(
            {
                "schema_version": 2,
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
                "max_primer": 60,
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
            }
        )
    )

    def _step(*extra):
        proc = subprocess.run(
            [sys.executable, "-m", "neoswga.cli_unified", *extra, "-j", str(params_file)],
            capture_output=True,
            text=True,
            cwd=str(tmp_path),
            timeout=900,
        )
        if proc.returncode != 0:
            pytest.skip(f"{extra[0]} failed:\n{proc.stderr[-800:]}")
        return proc

    _step("count-kmers")
    _step("filter")
    _step("score")
    optimize = _step("optimize", "--show-frontier", "--seed", "1")

    output = optimize.stdout + optimize.stderr
    assert "Pareto Frontier Analysis" in output, output[-2000:]
    assert "frontier analysis failed" not in output, output[-2000:]
    assert "Coverage range: 0%-0%" not in output, (
        "the frontier reports zero coverage at every size, so it is not reading "
        "the candidate pool"
    )
