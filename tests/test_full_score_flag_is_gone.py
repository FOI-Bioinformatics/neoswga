"""--full-score cost 126 times as much for no change.

Measured on the 449-candidate E. coli pool: 767.6 s against 6.1 s for
--amp-model alone, for a mean absolute change of 0.0016 on an 11.1 to 19.1
scale, Pearson 1.0000, Spearman 0.9999, and an identical delivered order. The
flag's help text claimed the delta-G features were under 2% of model accuracy;
the measurement puts it under 0.02%. Nothing in tests/ ever exercised
create_augmented_df(skip_delta_g=False).
"""

import subprocess
import sys


def test_the_flag_is_rejected():
    result = subprocess.run(
        [sys.executable, "-m", "neoswga.cli_unified", "score", "--full-score"],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "unrecognized arguments" in result.stderr or "--full-score" in result.stderr


def test_amp_model_survives():
    """The flag that works correctly, and the only thing that makes the
    quality column vary, is kept."""
    result = subprocess.run(
        [sys.executable, "-m", "neoswga.cli_unified", "score", "--help"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert "--amp-model" in result.stdout
    assert "--full-score" not in result.stdout


def test_fast_score_is_the_only_mode(monkeypatch, tmp_path):
    """parameter.fast_score is set True unconditionally now."""
    import json

    import neoswga.cli.pipeline as cli_pipeline
    from neoswga.core import parameter
    from neoswga.core import pipeline as pipeline_mod

    fasta = tmp_path / "t.fasta"
    fasta.write_text(">t\n" + "ACGT" * 500 + "\n")
    (tmp_path / "params.json").write_text(
        json.dumps(
            {
                "schema_version": 2,
                "data_dir": str(tmp_path),
                "fg_genomes": [str(fasta)],
                "fg_prefixes": [str(tmp_path / "t")],
                "bg_genomes": [],
                "bg_prefixes": [],
                "fg_seq_lengths": [2000],
                "cpus": 1,
            }
        )
    )

    monkeypatch.setattr(pipeline_mod, "step3", lambda *a, **k: None)
    monkeypatch.setattr(pipeline_mod, "_initialize", lambda: None)
    monkeypatch.setattr(cli_pipeline, "setup_gpu_acceleration", lambda *a, **k: None)
    monkeypatch.setattr(cli_pipeline, "report_unimplemented_options", lambda *a, **k: None)
    parameter.fast_score = False

    class _Args:
        json_file = str(tmp_path / "params.json")
        quiet = True
        seed = None
        enable_qa = False
        amp_model = False
        min_amp_pred = None
        fast_score = False

    cli_pipeline.run_step3(_Args())

    assert parameter.fast_score is True


def test_no_source_file_still_mentions_the_flag():
    """A removed flag that survives in a docstring or an example is worse than
    one that never existed."""
    import subprocess

    result = subprocess.run(
        ["grep", "-rn", "--include=*.py", "--include=*.md", "full.score", "."],
        capture_output=True,
        text=True,
    )
    offenders = [
        line
        for line in result.stdout.splitlines()
        if "AUDIT_pipeline_four_steps" not in line
        and "2026-09-06-observability-and-startup" not in line
        and "test_full_score_flag_is_gone" not in line
        and "fast_score" not in line
        and "fast-score" not in line
    ]
    assert offenders == [], offenders
