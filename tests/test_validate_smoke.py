"""`validate --quick` opened no genome, no params.json and no pipeline step.

It ran three synthetic tests -- adaptive GC filter, Bloom filter, network
optimization. `validate params -j` ran the schema validator. So there was no
seconds-long check of a configuration before a thirty-second filter, or a
sixteen-minute one against hg38.
"""

import json
import shutil

import pytest


def _needs_jellyfish():
    if not shutil.which("jellyfish"):
        pytest.skip("jellyfish not available (required for count-kmers)")


def test_the_smoke_genomes_are_packaged():
    """examples/ is not package data, so a smoke mode reading from there works
    in a checkout and fails for every pip install."""
    from pathlib import Path

    import neoswga.core.smoke as smoke

    root = Path(smoke.__file__).parent
    assert (root / "pcDNA.fasta").is_file()
    assert (root / "pLTR.fasta").is_file()


def test_pyproject_ships_the_smoke_package():
    from pathlib import Path

    # Anchored on this file, not the working directory: several tests chdir,
    # and under xdist this ran in a worker whose cwd was elsewhere.
    root = Path(__file__).resolve().parent.parent
    text = (root / "pyproject.toml").read_text()
    assert '"neoswga.core.smoke"' in text, "the smoke package is not in `packages`"
    assert "*.fasta" in text, "the smoke FASTAs are not in package-data"


def test_smoke_runs_all_four_steps_and_passes():
    _needs_jellyfish()
    from neoswga.core.validation import smoke_validation

    assert smoke_validation(params_path=None, verbose=False) is True


def test_smoke_applies_the_users_chemistry(tmp_path, caplog):
    """The point is to check the user's configuration, not a canned one."""
    import logging

    from neoswga.core.validation import smoke_validation

    _needs_jellyfish()

    fasta = tmp_path / "t.fasta"
    fasta.write_text(">t\n" + "ACGT" * 500 + "\n")
    params = tmp_path / "params.json"
    params.write_text(
        json.dumps(
            {
                "data_dir": str(tmp_path),
                "fg_genomes": [str(fasta)],
                "fg_prefixes": [str(tmp_path / "t")],
                "polymerase": "bst",
                "reaction_temp": 63.0,
                "min_tm": 40,
                "max_tm": 80,
                "min_k": 10,
                "max_k": 12,
            }
        )
    )

    with caplog.at_level(logging.INFO):
        smoke_validation(params_path=str(params), verbose=True)

    assert "bst" in caplog.text


def test_smoke_reports_an_unknown_key_in_the_users_params(tmp_path, caplog):
    import logging

    from neoswga.core.validation import smoke_validation

    params = tmp_path / "params.json"
    params.write_text(
        json.dumps(
            {
                "data_dir": str(tmp_path),
                "fg_genomes": [],
                "fg_prefixes": [],
                "max_bg_freqency": 5e-6,
            }
        )
    )

    with caplog.at_level(logging.WARNING):
        smoke_validation(params_path=str(params), verbose=True)

    assert "max_bg_freqency" in caplog.text


def test_smoke_reports_a_missing_foreground_genome(tmp_path, caplog):
    """The single most common real failure, and `--quick` never looked."""
    import logging

    from neoswga.core.validation import smoke_validation

    params = tmp_path / "params.json"
    params.write_text(
        json.dumps(
            {
                "data_dir": str(tmp_path),
                "fg_genomes": [str(tmp_path / "absent.fasta")],
                "fg_prefixes": [str(tmp_path / "absent")],
            }
        )
    )

    with caplog.at_level(logging.WARNING):
        ok = smoke_validation(params_path=str(params), verbose=True)

    assert ok is False
    assert "absent.fasta" in caplog.text


def test_the_cli_exposes_the_flag():
    from neoswga import cli_unified

    args = cli_unified.create_parser().parse_args(["validate", "--smoke"])
    assert args.smoke is True


def test_quick_is_unchanged():
    """--quick keeps its historical meaning; --smoke is a separate mode."""
    from neoswga import cli_unified

    args = cli_unified.create_parser().parse_args(["validate", "--quick"])
    assert args.quick is True
    assert getattr(args, "smoke", False) is False


def test_the_exit_code_distinguishes_pass_from_fail(tmp_path):
    """A check whose failure exits 0 is useless in CI, and every assertion
    above reads the return value rather than the process status."""
    import subprocess
    import sys

    _needs_jellyfish()

    bad = tmp_path / "bad.json"
    bad.write_text(
        json.dumps(
            {
                "data_dir": str(tmp_path),
                "fg_genomes": [str(tmp_path / "absent.fasta")],
                "fg_prefixes": [str(tmp_path / "absent")],
            }
        )
    )

    failing = subprocess.run(
        [sys.executable, "-m", "neoswga.cli_unified", "validate", "--smoke", "-j", str(bad)],
        capture_output=True,
        text=True,
    )
    assert failing.returncode != 0, failing.stdout + failing.stderr

    passing = subprocess.run(
        [sys.executable, "-m", "neoswga.cli_unified", "validate", "--smoke"],
        capture_output=True,
        text=True,
    )
    assert passing.returncode == 0, passing.stdout + passing.stderr
