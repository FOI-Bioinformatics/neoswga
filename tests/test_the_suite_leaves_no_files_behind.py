"""Running the suite must leave the working tree as it found it.

Three compounding causes put files in the repository root:

- `parameter.data_dir` has no module-level default. `get_params` creates it and
  it then persists for the whole pytest process, so a test that skips
  `get_params` inherits whatever the last one left. The plasmid example's
  params.json carries `"data_dir": "./"`, a RELATIVE path, so the value that
  leaks resolves against whatever the current directory happens to be.
- Two tests call a `run_stepN` CLI handler with `core.pipeline._initialize`
  stubbed out, so their own `data_dir` is never read.
- `unified_optimizer` wrote `step4_improved_df_validation.json` with an
  `or os.getcwd()` fallback, so it wrote to the current directory even when no
  `data_dir` was configured at all.

`.gitignore` carried a section headed "Root-level pipeline run artifacts",
which is the previous response to this: ignore the files rather than stop
writing them. Those entries are removed, and this test is what keeps them
removed.
"""

import os
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def _git(*args):
    return subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, timeout=60
    ).stdout


def test_the_repository_root_holds_no_pipeline_artifacts():
    """The files this actually produced, by name.

    A run manifest and a step-4 validation report have both appeared here.
    Neither is a source file and neither is ignored any more, so either would
    show up as untracked.
    """
    strays = [
        name
        for name in (
            "run_manifest.json",
            "step4_improved_df_validation.json",
            "step2_df.csv",
            "step3_df.csv",
            "step4_improved_df.csv",
            "filter_stats.json",
            "genome_gc.json",
        )
        if (ROOT / name).exists()
    ]
    assert not strays, (
        f"pipeline artifacts in the repository root: {strays}. A test wrote to "
        f"the process working directory, which under pytest is the repo root."
    )


def test_no_untracked_pipeline_artifact_appeared_anywhere():
    """The general form, so a file this test does not name by hand is caught.

    Scoped to artifact SHAPES rather than to any untracked file: a developer
    mid-change legitimately has untracked source. What must never appear is a
    pipeline output that no run was asked to produce.
    """
    import fnmatch

    patterns = (
        "*_df.csv",
        "*_df_summary.json",
        "*_df_validation.json",
        "run_manifest.json",
        "filter_stats.json",
        "genome_gc.json",
        "*mer_all.txt",
        "*_positions.h5",
    )
    untracked = [
        line[3:] for line in _git("status", "--porcelain").splitlines() if line.startswith("?? ")
    ]
    strays = [
        path
        for path in untracked
        if any(fnmatch.fnmatch(os.path.basename(path.rstrip("/")), pat) for pat in patterns)
    ]
    assert not strays, (
        f"untracked pipeline artifacts: {strays}. Something wrote outside the "
        f"data_dir it was given."
    )


def test_the_gitignore_does_not_paper_over_root_artifacts():
    """The entries that were there instead of a fix.

    Root-anchored ignores for pipeline outputs mean someone saw these appear
    and chose to hide them. Removing the writes is the fix; this stops the
    ignores coming back in place of one.
    """
    text = (ROOT / ".gitignore").read_text()
    for entry in ("/step4_improved_df_validation.json", "/step*_df.csv"):
        assert entry not in text, (
            f"{entry!r} is ignored at the repository root. If a run writes there, "
            f"stop the write rather than hiding the file."
        )


def test_the_validation_report_needs_a_configured_data_dir(tmp_path, monkeypatch):
    """No data_dir means no file, rather than a file in the current directory.

    `or os.getcwd()` made an unset data_dir mean "here", which is the fallback
    that put step4_improved_df_validation.json in the repository root.
    """
    import neoswga.core.unified_optimizer as uo
    from neoswga.core import parameter

    monkeypatch.chdir(tmp_path)
    monkeypatch.delattr(parameter, "data_dir", raising=False)

    uo._write_validation_report({"optimizer": "clique", "ok": True})

    assert list(tmp_path.iterdir()) == [], "a report was written with no data_dir set"


def test_the_validation_report_is_written_where_data_dir_points(tmp_path, monkeypatch):
    import neoswga.core.unified_optimizer as uo
    from neoswga.core import parameter

    target = tmp_path / "results"
    target.mkdir()
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(parameter, "data_dir", str(target), raising=False)

    uo._write_validation_report({"optimizer": "clique", "ok": True})

    assert (target / "step4_improved_df_validation.json").is_file()
