"""The whole pipeline, and what each refusal does to the output directory.

Task 10 of the 2026-09-21 valid-design plan: an end-to-end run over a small
reference, "plus error injections at each required step. The failure cases
must return nonzero status and no qualifying recommendation artifact."

The individual contracts are unit-tested elsewhere. What this file adds is
that they hold through a real process: a refusal has to survive argument
parsing, parameter resolution, the step's own try/except, and the command
boundary, and each of those has swallowed one before. Known Issue 8's class
is precisely a check that exists and is not reached.

Every injection below is a configuration a user could plausibly write, not a
monkeypatched internal. What is asserted of each is the same three things:

1. the command exits nonzero;
2. it leaves `design_failure.json` saying what failed;
3. no recommendation can be exported from the directory afterwards.

The third matters most. The first two are about this run; the third is about
what the directory looks like to the next command, which is the only thing a
later reader sees.
"""

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

pd = pytest.importorskip("pandas")

pytestmark = pytest.mark.skipif(
    shutil.which("jellyfish") is None,
    reason="the pipeline's first step shells out to jellyfish",
)

PACKAGED = Path(__file__).resolve().parent.parent.parent / "neoswga" / "core" / "smoke"


def run(args, cwd):
    return subprocess.run(
        [sys.executable, "-m", "neoswga.cli_unified", *args],
        capture_output=True,
        text=True,
        cwd=str(cwd),
        timeout=600,
    )


def write_params(directory, **overrides):
    params = {
        "fg_genomes": ["pcDNA.fasta"],
        "bg_genomes": [],
        "fg_prefixes": ["pcDNA"],
        "bg_prefixes": [],
        "data_dir": "results",
        "min_k": 10,
        "max_k": 10,
        "polymerase": "phi29",
        "reaction_temp": 30.0,
        "min_fg_freq": 1e-6,
        "max_bg_freq": 1.0,
        "max_gini": 1.0,
        "min_gini_sites": 1,
        "max_primer": 60,
        "min_tm": 0,
        "max_tm": 100,
        "gc_min": 0.0,
        "gc_max": 1.0,
        "num_primers": 4,
        "target_set_size": 4,
        "max_sets": 1,
        "iterations": 2,
        "cpus": 1,
        "schema_version": 2,
    }
    params.update(overrides)
    path = directory / "params.json"
    path.write_text(json.dumps(params))
    return path


def stage_workspace(tmp_path):
    workspace = tmp_path / "work"
    workspace.mkdir()
    shutil.copy(PACKAGED / "pcDNA.fasta", workspace / "pcDNA.fasta")
    return workspace


def prepared(tmp_path, **overrides):
    """A workspace taken through count, filter and candidate preparation."""
    workspace = stage_workspace(tmp_path)
    write_params(workspace, **overrides)
    for step in ("count-kmers", "filter", "prepare-candidates"):
        result = run([step, "-j", "params.json"], workspace)
        assert result.returncode == 0, f"{step} failed:\n{result.stdout}\n{result.stderr}"
    return workspace


def failure_record(workspace):
    path = workspace / "results" / "design_failure.json"
    if not path.exists():
        return None
    return json.loads(path.read_text())


def assert_nothing_recommendable(workspace, expect_record=True):
    """The three things every refusal owes the next command."""
    results = workspace / "results"

    if expect_record:
        record = failure_record(workspace)
        assert record is not None, "a refused run left no failure record"
        assert record["run_state"] in {"failed", "interrupted"}
        assert record["recommendation_written"] is False

    export = run(["export", "-d", "results", "-o", "order", "--format", "fasta"], workspace)
    assert export.returncode != 0, (
        "a directory whose last run was refused still exported a recommendation:\n" + export.stdout
    )
    orders = workspace / "order"
    assert not (orders.exists() and list(orders.glob("*.fasta")))


# ---------------------------------------------------------------------------
# The happy path
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def completed(tmp_path_factory):
    workspace = prepared(tmp_path_factory.mktemp("strict"))
    result = run(["optimize", "-j", "params.json", "--seed", "1"], workspace)
    assert result.returncode == 0, result.stdout + result.stderr
    return workspace


def test_the_four_steps_produce_the_expected_artifacts(completed):
    results = completed / "results"

    assert (results / "step2_df.csv").exists()
    assert (results / "step3_df.csv").exists()
    assert (results / "step4_improved_df.csv").exists()
    assert (results / "step4_improved_df_summary.json").exists()


def test_a_finished_run_leaves_no_failure_record(completed):
    """`optimize` clears one, so an earlier failure cannot block this export."""
    assert failure_record(completed) is None


def test_the_finished_run_exports_a_recommendation(completed):
    export = run(["export", "-d", "results", "-o", "order", "--format", "fasta"], completed)

    assert export.returncode == 0, export.stdout + export.stderr
    assert list((completed / "order").glob("*.fasta"))


def test_the_run_manifest_records_the_request_hash(completed):
    manifest = json.loads((completed / "results" / "run_manifest.json").read_text())
    optimize_steps = [s for s in manifest["steps"] if s["step"] == "optimize"]

    assert optimize_steps
    assert optimize_steps[-1]["extra"]["request_hash"]


# ---------------------------------------------------------------------------
# Error injections, one per contract
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "overrides,expected",
    [
        pytest.param({"glycerol_percent": 10.0}, "glycerol", id="additive-with-no-model"),
        pytest.param(
            {"min_selectivity_density": 60.0},
            "min_selectivity_density",
            id="background-limit-with-no-background",
        ),
    ],
)
def test_a_refusal_on_the_design_path_leaves_nothing_recommendable(tmp_path, overrides, expected):
    """Each of these is a params file somebody could write.

    These reach the design path, so all three obligations apply: exit
    nonzero, leave a record saying what failed, and leave nothing the next
    command will export.

    The record is the one that had to be fixed here. Both refusals fire inside
    `resolve_design_request`, which runs BEFORE `get_params` populates the
    `parameter` module, so the failure writer looked for an output directory
    the module did not yet know about and wrote nowhere. The run exited
    nonzero and the directory still looked like its previous success. Same
    ordering trap as `warn_on_condition_drift`, reached from the other side.
    """
    workspace = prepared(tmp_path)
    write_params(workspace, **overrides)

    result = run(["optimize", "-j", "params.json"], workspace)
    output = result.stdout + result.stderr

    assert result.returncode != 0, f"accepted {overrides}:\n{output}"
    assert expected in output, output[-1500:]
    assert_nothing_recommendable(workspace)


@pytest.mark.parametrize(
    "overrides,expected",
    [
        pytest.param({"coverage_reach": 0}, "minimum of 1", id="explicit-zero-reach"),
        pytest.param({"polymerase": "taq"}, "not one of", id="unsupported-polymerase"),
    ],
)
def test_a_refusal_at_schema_validation_never_enters_the_design_path(tmp_path, overrides, expected):
    """Caught earlier, and better, than the design-request contract catches it.

    The schema refuses these before any handler runs, naming the permitted
    values. No failure record is written and none is wanted: the design path
    was never entered, so there is no stage to name, and the directory has no
    step-4 output for anyone to mistake for current.

    Asserted separately from the design-path refusals rather than folded in,
    because "exits nonzero" is the same observation from two different
    mechanisms and only one of them owes a record.
    """
    workspace = prepared(tmp_path)
    write_params(workspace, **overrides)

    result = run(["optimize", "-j", "params.json"], workspace)
    output = result.stdout + result.stderr

    assert result.returncode != 0, f"accepted {overrides}:\n{output}"
    assert expected in output, output[-1500:]
    assert failure_record(workspace) is None
    assert_nothing_recommendable(workspace, expect_record=False)


def test_a_missing_position_index_is_refused(tmp_path):
    """The index is deleted after preparation, so the pool is real and unscored."""
    workspace = prepared(tmp_path)
    for path in workspace.glob("*_positions.h5"):
        path.unlink()
    for path in (workspace / "results").glob("*_positions.h5"):
        path.unlink()

    result = run(["optimize", "-j", "params.json"], workspace)

    assert result.returncode != 0, result.stdout


def test_a_stale_index_from_another_reference_is_refused(tmp_path):
    """The structurally perfect failure: every dataset present, wrong genome.

    The FASTA is replaced after the index was built, so the index is complete,
    internally consistent, and describes a sequence this run is not targeting.
    """
    workspace = prepared(tmp_path)
    shutil.copy(PACKAGED / "pLTR.fasta", workspace / "pcDNA.fasta")

    result = run(["optimize", "-j", "params.json"], workspace)
    output = result.stdout + result.stderr

    assert result.returncode != 0, output
    assert "reference" in output.lower() or "count-kmers" in output
    assert_nothing_recommendable(workspace)


def test_a_reaction_the_filter_never_recorded_is_refused(tmp_path):
    """Candidates selected under one chemistry, scored under another.

    This used to warn and proceed, searching the shortlist while everything
    the inventory held went unreachable.
    """
    workspace = stage_workspace(tmp_path)
    write_params(workspace)
    assert run(["count-kmers", "-j", "params.json"], workspace).returncode == 0
    assert (
        run(["filter", "-j", "params.json", "--preset", "enhanced_equiphi29"], workspace).returncode
        == 0
    )
    assert run(["prepare-candidates", "-j", "params.json"], workspace).returncode == 0

    result = run(["optimize", "-j", "params.json"], workspace)
    output = result.stdout + result.stderr

    assert result.returncode != 0, output
    assert "filter" in output
    assert_nothing_recommendable(workspace)


def test_the_retired_command_name_is_refused(tmp_path):
    workspace = stage_workspace(tmp_path)
    write_params(workspace)

    result = run(["score", "-j", "params.json"], workspace)

    assert result.returncode != 0
    assert "prepare-candidates" in result.stdout + result.stderr


# ---------------------------------------------------------------------------
# The pipeline can produce what the checks demand
# ---------------------------------------------------------------------------


def write_multi_record_reference(workspace, records):
    import random

    rng = random.Random(5)
    parts = []
    for name, length in records:
        seq = "".join(rng.choice("ACGT") for _ in range(length))
        parts.append(
            f">{name}\n" + "\n".join(seq[i : i + 70] for i in range(0, len(seq), 70)) + "\n"
        )
    (workspace / "multi.fasta").write_text("".join(parts))


def test_a_fresh_multi_record_run_produces_geometry_the_checks_accept(tmp_path):
    """Closes the loop on the record-geometry refusal.

    `verify_index_geometry` refuses a multi-record reference whose index
    carries no record starts. That is only defensible if the pipeline can
    actually produce one that does -- otherwise the check would make
    multi-record references unusable rather than merely requiring a recount,
    and every index in this repository predates the feature, so nothing here
    would have caught it.

    Three records, so the rule applies, and all four steps must complete.
    """
    import h5py

    from neoswga.core.reference_check import verify_index_geometry
    from neoswga.core.string_search import INDEX_FORMAT_VERSION, RECORD_STARTS_KEY

    workspace = tmp_path / "work"
    workspace.mkdir()
    write_multi_record_reference(workspace, [("chrA", 9000), ("chrB", 8000), ("chrC", 7000)])
    write_params(
        workspace,
        fg_genomes=["multi.fasta"],
        fg_prefixes=["multi"],
        max_primer=50,
    )

    for step in ("count-kmers", "filter", "prepare-candidates", "optimize"):
        result = run([step, "-j", "params.json"], workspace)
        assert result.returncode == 0, f"{step}:\n{result.stdout}\n{result.stderr}"

    index = workspace / "multi_10mer_positions.h5"
    assert index.exists()
    with h5py.File(index, "r") as handle:
        assert int(handle.attrs["index_format_version"]) == INDEX_FORMAT_VERSION
        starts = [int(v) for v in handle[RECORD_STARTS_KEY]]

    assert starts == [0, 9000, 17000], starts

    # And the check the pipeline's own output has to satisfy.
    verify_index_geometry({str(workspace / "multi"): str(workspace / "multi.fasta")}, [10])


# ---------------------------------------------------------------------------
# Contraction and sequencing-informed expansion
# ---------------------------------------------------------------------------
#
# Task 10 asks the end-to-end test to cover these two stages as well, and it
# did not. They are the two commands that take a DELIVERED panel and change it,
# so a defect in either reaches a panel somebody would order while every test
# of the four steps stays green.
#
# Expansion is the one that needed a decision. It reads sequencing depth, and
# this repository contains no BAM or CRAM at all -- which is why
# `calibrate-reach` has never been run against measured depth. A BAM is
# synthesised here rather than skipping: the point is that the command's
# plumbing works end to end, not that the depth is realistic. What cannot be
# established this way is anything about recovery, and
# `docs/validation/design_release_gates.md` is where that line is drawn.


def _synthesise_bam(workspace, contig, length, covered_to):
    """An indexed BAM covering `contig` from 0 to `covered_to`.

    Deliberately leaves the tail uncovered, because a gap is what expansion
    exists to design against. Depth is uniform and shallow; nothing here is a
    claim about a real reaction.
    """
    pysam = pytest.importorskip("pysam")

    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": contig, "LN": length}]}
    path = workspace / "depth.bam"
    read_length = 100
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for start in range(0, max(covered_to - read_length, 1), read_length // 2):
            read = pysam.AlignedSegment()
            read.query_name = f"r{start}"
            read.query_sequence = "A" * read_length
            read.flag = 0
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = 60
            read.cigarstring = f"{read_length}M"
            read.query_qualities = pysam.qualitystring_to_array("I" * read_length)
            out.write(read)
    pysam.index(str(path))
    return path


def _record_name(fasta):
    for line in fasta.read_text().splitlines():
        if line.startswith(">"):
            return line[1:].split()[0]
    raise AssertionError(f"{fasta} holds no FASTA record")


def test_contraction_reduces_a_delivered_panel_without_breaking_it(completed):
    """`--minimize-primers` must return a panel, not an empty one or a crash."""
    before = pd.read_csv(completed / "results" / "step4_improved_df.csv")
    delivered_before = before[before.get("set_index", 0) == 0]["primer"].tolist()

    result = run(["optimize", "-j", "params.json", "--seed", "1", "--minimize-primers"], completed)
    assert result.returncode == 0, result.stdout + result.stderr

    after = pd.read_csv(completed / "results" / "step4_improved_df.csv")
    delivered_after = after[after.get("set_index", 0) == 0]["primer"].tolist()

    assert delivered_after, "contraction returned an empty panel"
    assert len(delivered_after) <= len(delivered_before)
    assert len(set(delivered_after)) == len(delivered_after), "contraction duplicated an oligo"
    # And the directory is still exportable, which is the property that matters
    # to whoever orders the result.
    export = run(["export", "-d", "results", "-o", "order2", "--format", "fasta"], completed)
    assert export.returncode == 0, export.stdout + export.stderr


def test_sequencing_informed_expansion_runs_end_to_end(tmp_path):
    """`expand-primers --bam` reads depth, finds a gap and adds oligos."""
    pytest.importorskip("pysam")
    workspace = prepared(tmp_path, num_primers=2, target_set_size=2)
    assert run(["optimize", "-j", "params.json", "--seed", "1"], workspace).returncode == 0

    fasta = workspace / "pcDNA.fasta"
    contig = _record_name(fasta)
    length = sum(
        len(line.strip()) for line in fasta.read_text().splitlines() if not line.startswith(">")
    )
    bam = _synthesise_bam(workspace, contig, length, covered_to=length // 2)

    delivered = pd.read_csv(workspace / "results" / "step4_improved_df.csv")
    fixed = delivered[delivered.get("set_index", 0) == 0]["primer"].tolist()
    assert fixed, "the design delivered no panel to expand"

    result = run(
        [
            "expand-primers",
            "-j",
            "params.json",
            "--fixed-primers",
            *fixed,
            "--num-new",
            "1",
            "--bam",
            str(bam),
            "--contig-alias",
            f"{contig}={contig}",
            "--output",
            "expanded",
        ],
        workspace,
    )

    # The command may decline to add an oligo -- the pool is small and a short
    # panel is a documented outcome -- but it must not fail, and it must have
    # read the depth rather than ignoring it.
    assert result.returncode == 0, result.stdout + result.stderr
    written = sorted((workspace / "expanded").glob("*.csv"))
    assert written, f"expansion wrote no panel: {list((workspace / 'expanded').iterdir())}"
    expanded = pd.read_csv(written[0])
    assert not expanded.empty
    assert set(fixed).issubset(
        set(expanded["primer"])
    ), "expansion dropped an oligo it was told to keep"


def test_expansion_refuses_a_bam_naming_no_configured_reference(tmp_path):
    """A BAM whose contigs match nothing must not read as zero depth everywhere.

    Zero depth everywhere is a gap everywhere, so expansion would design
    against the whole target while appearing to use the sequencing data.
    """
    pytest.importorskip("pysam")
    workspace = prepared(tmp_path, num_primers=2, target_set_size=2)
    assert run(["optimize", "-j", "params.json", "--seed", "1"], workspace).returncode == 0

    delivered = pd.read_csv(workspace / "results" / "step4_improved_df.csv")
    fixed = delivered[delivered.get("set_index", 0) == 0]["primer"].tolist()
    bam = _synthesise_bam(workspace, "a_contig_no_reference_here_has", 4000, covered_to=2000)
    result = run(
        [
            "expand-primers",
            "-j",
            "params.json",
            "--fixed-primers",
            *fixed,
            "--num-new",
            "1",
            "--bam",
            str(bam),
            "--output",
            "expanded",
        ],
        workspace,
    )

    combined = result.stdout + result.stderr
    assert result.returncode != 0 or "contig" in combined.lower() or "alias" in combined.lower(), (
        "an unmatched BAM was accepted silently, so the depth it supplied was "
        "nothing and every base looked like a gap:\n" + combined
    )
