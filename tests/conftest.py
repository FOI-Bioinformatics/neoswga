"""Shared fixtures for the neoswga test suite.

Centralizes commonly duplicated test data: primer sequences, genome strings,
reaction conditions, mock caches, and temporary FASTA files.
"""

import glob
import json
import logging
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pytest

from neoswga.core.reaction_conditions import ReactionConditions

_EXAMPLE_DIR = os.path.join(os.path.dirname(__file__), "..", "examples", "plasmid_example")


def plasmid_example_ready() -> bool:
    """Whether the generated artifacts the dependent tests need are present.

    NOT `os.path.isdir(_EXAMPLE_DIR)`, which nineteen test files used. That
    directory is COMMITTED -- git tracks the README, both FASTAs, params.json
    and step2_df.csv.original -- so it is always there and the question always
    answered yes. What those tests actually need is what
    `_prime_plasmid_example` builds, and that runs only when a counter is on
    PATH. Guarding on the directory meant they failed with a missing-k-mer-file
    error rather than skipping with a reason.

    A table is asked for through `table_exists`, in whichever form the counter
    that ran wrote it. This globbed `*mer_all.txt`, which is jellyfish's
    artifact: once KMC was installed alongside it, priming wrote binary
    databases, no text appeared, and 48 tests skipped as unprepared while the
    directory was in fact ready. A guard that reads "not prepared" when
    everything is prepared costs exactly the coverage it was written to
    protect.
    """
    from neoswga.core import kmer_tables

    # `discover_prefixes` rather than stripping a suffix here: a KMC database
    # is `{prefix}_{k}mer.kmc_pre`, and cutting at the last underscore cuts
    # inside `kmc_pre`, which silently answered "not prepared" for a directory
    # that was.
    counted = any(kmer_tables.discover_prefixes(_EXAMPLE_DIR, k) for k in range(4, 31))
    return counted and os.path.exists(os.path.join(_EXAMPLE_DIR, "step3_df.csv"))


def _run_priming():
    """Run count-kmers, filter and score in the example directory.

    Separated from the fixture so a caller can observe it failing. The fixture
    used to swallow every exception here with a bare `except Exception: pass`,
    which hid a priming failure even when jellyfish WAS present and left
    nineteen files failing for a reason that pointed at the wrong thing.
    """
    import neoswga.core.pipeline as pipeline_mod
    from neoswga.core import parameter

    def _reset():
        for attr in (
            "fg_prefixes",
            "bg_prefixes",
            "fg_genomes",
            "bg_genomes",
            "fg_seq_lengths",
            "bg_seq_lengths",
            "fg_circular",
            "bg_circular",
        ):
            setattr(pipeline_mod, attr, None)
        # `_initialize` is guarded by `_initialized` and by whether
        # `parameter.json_file` changed, and priming changes neither between
        # steps. So nulling the globals without clearing this flag left step 2
        # reading `fg_prefixes + bg_prefixes` as None + None, and priming has
        # failed with a TypeError every time it actually ran. Every other reset
        # helper in the suite clears it; this one did not.
        pipeline_mod._initialized = False
        parameter.json_file = "params.json"

    cwd = os.getcwd()
    try:
        os.chdir(_EXAMPLE_DIR)
        _reset()
        pipeline_mod._initialize()
        pipeline_mod.step1()
        _reset()
        pipeline_mod.step2()
        _reset()
        pipeline_mod.step3()
    finally:
        _reset()
        os.chdir(cwd)


# Priming writes into one shared directory, and this fixture is SESSION scoped,
# so under `-n 8` eight sessions run it at once. Observed on a fresh checkout:
# the workers interleave and leave the directory half-built -- one k-mer table
# of the expected set, no `step3_df.csv` -- after which
# `plasmid_example_ready()` is false and 51 tests skip that would otherwise
# run. An earlier interleaving produced 23 failures and 2,548 errors, all of
# them `validate_step2_prerequisites` reading prefixes that priming had reset.
#
# A directory is the lock because `os.mkdir` is atomic on every platform this
# runs on, and it lives outside the repository so it cannot appear in
# `git status` and trip `test_the_suite_leaves_no_files_behind`.
_PRIMING_LOCK = pathlib.Path(tempfile.gettempdir()) / "neoswga-plasmid-priming.lock"
_PRIMING_TIMEOUT_S = 600.0


def _claim_priming_lock() -> bool:
    """True when this worker is the one that should prime.

    A lock older than the timeout is treated as abandoned, because a worker
    killed mid-prime would otherwise block every future run on this machine.
    """
    try:
        _PRIMING_LOCK.mkdir()
        return True
    except FileExistsError:
        pass
    try:
        age = time.time() - _PRIMING_LOCK.stat().st_mtime
    except OSError:
        return False
    if age > _PRIMING_TIMEOUT_S:
        logging.getLogger(__name__).warning(
            "Priming lock %s is %.0f s old; treating it as abandoned.",
            _PRIMING_LOCK,
            age,
        )
        return True
    return False


def _release_priming_lock() -> None:
    try:
        _PRIMING_LOCK.rmdir()
    except OSError:
        pass


def _wait_for_priming() -> None:
    """Block until the worker holding the lock finishes, or give up.

    Giving up is not a failure: the dependent tests guard on
    `plasmid_example_ready()` and skip with a reason, which is what they do
    without jellyfish too.
    """
    deadline = time.monotonic() + _PRIMING_TIMEOUT_S
    while time.monotonic() < deadline:
        if plasmid_example_ready():
            return
        if not _PRIMING_LOCK.exists():
            return
        time.sleep(0.5)
    logging.getLogger(__name__).warning(
        "Waited %.0f s for another worker to prime examples/plasmid_example; "
        "the tests that need it will skip.",
        _PRIMING_TIMEOUT_S,
    )


@pytest.fixture(scope="session", autouse=True)
def _prime_plasmid_example():
    """Prime the plasmid example with pipeline outputs once per session.

    Several tests depend on generated files in ``examples/plasmid_example/`` that
    are NOT committed (they are gitignored build artifacts):

    - Integration tests copy the directory (incl. the ``*_Nmer_all.txt`` k-mer
      files) into a tmpdir.
    - The swap-primer / contract-set / rescore-set CLI smoke tests run the CLI
      with ``cwd=examples/plasmid_example`` and need ``step2_df.csv`` /
      ``step3_df.csv`` and the ``*_positions.h5`` files.

    On a clean checkout / CI these do not exist, so those tests fail. This
    session-scoped autouse fixture runs ``count-kmers`` -> ``filter`` ->
    ``score`` once in the example directory to generate them (the default 6-12
    k range covers every range the tests use). It is a no-op when the outputs
    already exist (a dev's primed directory) or when jellyfish is unavailable
    (the pipeline's hard dependency); dependent tests then skip or surface the
    missing dependency themselves.
    """
    if not os.path.isdir(_EXAMPLE_DIR):
        return
    if plasmid_example_ready():
        return  # already primed

    try:
        from neoswga.core.kmer_counter import check_jellyfish_available
    except Exception as exc:  # pragma: no cover - import guard
        logging.getLogger(__name__).warning(
            "Cannot check for jellyfish (%s); the plasmid example is unprimed and "
            "the tests that need it will skip.",
            exc,
        )
        return
    if not check_jellyfish_available():
        logging.getLogger(__name__).warning(
            "jellyfish is not on PATH; the plasmid example is unprimed and the "
            "tests that need its generated artifacts will skip."
        )
        return

    if not _claim_priming_lock():
        # Another worker is priming the same directory. Wait for it rather than
        # racing: this fixture is SESSION scoped, and under `-n 8` that means
        # eight sessions, each running it against one shared directory.
        _wait_for_priming()
        return

    try:
        _run_priming()
    except Exception:
        # Reported, not swallowed. This used to be `except Exception: pass`, so
        # a priming failure WITH jellyfish present was invisible and nineteen
        # test files then failed with a missing-k-mer-file error that pointed
        # at the wrong thing.
        logging.getLogger(__name__).warning(
            "Priming examples/plasmid_example failed; the tests that need its "
            "generated artifacts will skip.",
            exc_info=True,
        )
    finally:
        _release_priming_lock()


@pytest.fixture
def plasmid_run(tmp_path):
    """A completed `optimize` over a private copy of the plasmid example.

    Steps 1 to 3 are already done: `_prime_plasmid_example` builds them once per
    session in the shared example directory. This copies that primed directory
    and runs step 4 in the copy, so the artifacts a test inspects belong to it
    alone and nothing is written back into `examples/`.

    A subprocess rather than an in-process call. The steps assign module globals
    in `neoswga.core.pipeline` and `neoswga.core.parameter`, and `_run_priming`
    above exists largely to reset them; a fixture that left them set would
    change what a later test in the same worker resolves.

    `params.json` names its genomes and `data_dir` relatively, so the command
    runs with `cwd` set to the copy and needs no path rewriting.
    """
    workdir = tmp_path / "plasmid"
    shutil.copytree(_EXAMPLE_DIR, workdir)

    completed = subprocess.run(
        [sys.executable, "-m", "neoswga.cli_unified", "optimize", "-j", "params.json"],
        cwd=workdir,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        pytest.fail(f"optimize failed in the copied example:\n{completed.stderr[-2000:]}")
    return workdir


# ---------------------------------------------------------------------------
# Primer sequences
# ---------------------------------------------------------------------------


@pytest.fixture
def sample_primers():
    """List of 8-mer primer sequences spanning different GC contents."""
    return [
        "ATCGATCG",  # 50% GC, balanced
        "GCTAGCTA",  # 50% GC, balanced
        "AATTCCGG",  # 50% GC, clustered
        "TTAACCGG",  # 50% GC, clustered
        "GGCCAATT",  # 50% GC, clustered
    ]


@pytest.fixture
def diverse_primers():
    """Primers spanning a wider range of GC content and lengths."""
    return [
        "ATCGATCG",  # 50% GC, 8-mer
        "GCGCGCGC",  # 100% GC, 8-mer
        "AAAATTTT",  # 0% GC, 8-mer
        "ATCGATCGATCG",  # 50% GC, 12-mer
        "GCTAGCTAGC",  # 50% GC, 10-mer
    ]


# ---------------------------------------------------------------------------
# Genome sequences
# ---------------------------------------------------------------------------


@pytest.fixture
def simple_genome():
    """A short repeating genome string (12 bp)."""
    return "ATCGATCGATCG"


@pytest.fixture
def repeat_genome():
    """A homopolymer genome string (10 bp)."""
    return "AAAAAAAAAA"


# ---------------------------------------------------------------------------
# Reaction conditions
# ---------------------------------------------------------------------------


@pytest.fixture
def phi29_conditions():
    """Standard phi29 reaction conditions (30 C)."""
    return ReactionConditions(temp=30.0, polymerase="phi29", mg_conc=2.5)


@pytest.fixture
def equiphi29_conditions():
    """EquiPhi29 reaction conditions (42 C) with DMSO and betaine."""
    return ReactionConditions(
        temp=42.0,
        polymerase="equiphi29",
        mg_conc=2.5,
        dmso_percent=5.0,
        betaine_m=1.0,
    )


# ---------------------------------------------------------------------------
# Mock position cache
# ---------------------------------------------------------------------------


class MockPositionCache:
    """Lightweight mock for PositionCache used across optimizer tests.

    Returns deterministic positions based on an optional positions dict,
    falling back to evenly spaced sites derived from primer length.
    """

    def __init__(self, positions_dict=None, genome_length=100_000):
        self._positions = positions_dict or {}
        self._genome_length = genome_length

    def get_positions(self, prefix, primer, strand="both"):
        key = (prefix, primer, strand)
        if key in self._positions:
            return np.asarray(self._positions[key])
        # Default: shorter primers bind more frequently
        n_sites = max(5, 20 - len(primer))
        spacing = self._genome_length // (n_sites + 1)
        return np.arange(spacing, spacing * (n_sites + 1), spacing)


@pytest.fixture
def mock_position_cache():
    """A MockPositionCache with default settings."""
    return MockPositionCache()


# ---------------------------------------------------------------------------
# Temporary FASTA files
# ---------------------------------------------------------------------------


@pytest.fixture
def tmp_fasta_file(tmp_path):
    """Create a minimal single-sequence FASTA file. Returns the Path."""
    fasta = tmp_path / "test_genome.fasta"
    fasta.write_text(">seq1\nATCGATCGATCGATCGATCG\nGCTAGCTAGCTAGCTAGCTA\n")
    return fasta


@pytest.fixture
def tmp_fasta_pair(tmp_path):
    """Create a foreground and background FASTA pair. Returns (fg, bg) Paths."""
    fg = tmp_path / "target.fasta"
    fg.write_text(">target_seq1\nATCGATCGATCGATCGATCG\nGCTAGCTAGCTAGCTAGCTA\n")
    bg = tmp_path / "background.fasta"
    bg.write_text(">bg_seq1\nAAAACCCCGGGGTTTTAAAA\nCCCCGGGGTTTTAAAACCCC\n")
    return fg, bg


# ---------------------------------------------------------------------------
# Failure context, for the one flake with no known cause
# ---------------------------------------------------------------------------
#
# On 2026-09-10 a full run gave 6 failures and 4124 passes where two immediate
# re-runs on a byte-identical tree gave none, with 25 concurrent agent
# processes the only unusual condition. The roadmap entry recording it says a
# single test name was not enough to diagnose it, which is why the entry cannot
# say more.
#
# Three known sources of load-dependent failure have since been removed: the
# wall-clock assertions, an over-broad caplog assertion and subprocess tests
# that skipped when a step failed. If anything is left, the next occurrence
# should arrive with its own evidence rather than needing a re-run to notice.
#
# Writes nothing on a green run.

_FAILURE_REPORT = Path(__file__).resolve().parent.parent / ".pytest_failure_context.json"
_FAILURES: list = []


def _load_average():
    try:
        return list(os.getloadavg())
    except (OSError, AttributeError):  # pragma: no cover - platform dependent
        return None


def _write_failure_context():
    """Record what a diagnosis of a load-dependent failure would need."""
    if not _FAILURES:
        return
    payload = {
        "failures": _FAILURES,
        "load_average": _load_average(),
        "workers": os.environ.get("PYTEST_XDIST_WORKER_COUNT"),
        "cpu_count": os.cpu_count(),
    }
    try:
        _FAILURE_REPORT.write_text(json.dumps(payload, indent=2))
    except OSError:  # pragma: no cover - diagnostics must never fail a run
        pass


@pytest.fixture(autouse=True)
def _no_test_leaks_the_working_directory():
    """Restore the working directory after every test, and name the leaker.

    Seven tests failed in a full run and passed alone, all of them reading a
    source path relative to the repository root. The cause was two bare
    `os.chdir(tmpdir)` calls in `test_core_optimize_blacklist_guard.py` with no
    restoration, which left a tmpdir as the working directory for every test
    sorting after that file. Those two only began running once
    `examples/plasmid_example` was primed, so the whole class was invisible on
    an unprimed checkout and on CI, which has no jellyfish.

    Restoring silently would hide the next one. Failing names it, at the test
    that leaked rather than at the seven that inherit it. Use
    `monkeypatch.chdir`, which restores on teardown, rather than `os.chdir`.
    """
    before = os.getcwd()
    yield
    after = os.getcwd()
    if after != before:
        os.chdir(before)
        pytest.fail(
            f"this test left the working directory at {after!r} instead of "
            f"{before!r}. Every later test that reads a repository-relative "
            f"path then fails, and passes when run alone. Use "
            f"`monkeypatch.chdir`, which restores on teardown."
        )


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    """Collect the failures, with the worker that ran each one."""
    outcome = yield
    report = outcome.get_result()
    if report.when == "call" and report.failed:
        _FAILURES.append(
            {
                "test": report.nodeid,
                "worker": os.environ.get("PYTEST_XDIST_WORKER", "master"),
                "duration": round(getattr(report, "duration", 0.0), 3),
            }
        )


def pytest_sessionfinish(session, exitstatus):
    _write_failure_context()
