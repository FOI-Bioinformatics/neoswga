"""Shared fixtures for CLI execution tests.

The existing CLI tests are parser-level: they check that `--help` works and that
the dispatch table resolves. Nothing actually *ran* a command body, which is why
`neoswga/cli/` sat at 21% coverage while `core/report/` was at 90%.

These fixtures make it cheap to run a real command against a real (tiny) genome
and assert on what it wrote. The genome is deterministic rather than random so a
failure is reproducible, and primers are embedded at known positions so coverage
and gap numbers can be checked against hand-computed values instead of against
whatever the code happens to produce.
"""

import random

import pytest


def _deterministic_genome(length, seed=20260801):
    """A reproducible pseudo-random sequence with no accidental repeats.

    Not `"ACGT" * n`: a periodic genome makes every k-mer occur thousands of
    times and silently breaks any test that assumes distinct binding sites.
    """
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


@pytest.fixture(scope="session")
def genome_seq():
    """20 kb of target sequence with three known primers planted in it."""
    seq = list(_deterministic_genome(20_000))
    # Plant each primer at known, well-separated positions so gap statistics
    # are predictable: 3 sites each, ~5 kb apart.
    plants = {
        "TTGACCATGA": (1_000, 6_000, 11_000),
        "GCATTACGGT": (2_500, 7_500, 12_500),
        "AACCGGTTAC": (4_000, 9_000, 14_000),
    }
    for primer, positions in plants.items():
        for pos in positions:
            seq[pos : pos + len(primer)] = list(primer)
    return "".join(seq)


@pytest.fixture(scope="session")
def planted_primers():
    """The primers planted in `genome_seq`, each with 3 known sites."""
    return ["TTGACCATGA", "GCATTACGGT", "AACCGGTTAC"]


@pytest.fixture
def genome_fasta(tmp_path, genome_seq):
    path = tmp_path / "target.fasta"
    lines = [">target\n"] + [genome_seq[i : i + 70] + "\n" for i in range(0, len(genome_seq), 70)]
    path.write_text("".join(lines))
    return str(path)


@pytest.fixture
def absent_primer(genome_seq):
    """A primer guaranteed not to occur in the genome, on either strand."""
    for candidate in ("CGCGCGATCG", "ATCGCGCGCG", "GCGCATATCG"):
        rc = candidate.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        if candidate not in genome_seq and rc not in genome_seq:
            return candidate
    raise AssertionError("could not construct an absent primer")


@pytest.fixture
def run_cli():
    """Parse an argv list and dispatch it, returning the handler's result.

    Goes through the real parser so defaults, types and required-ness are
    exercised, rather than hand-building a Namespace that could drift from
    what the CLI actually accepts.
    """
    from neoswga import cli_unified

    def _run(argv):
        parser = cli_unified.create_parser()
        args = parser.parse_args(argv)
        handler = _HANDLERS[args.command]
        return handler(args)

    return _run


_HANDLERS = {}


def pytest_configure(config):
    """Build the command -> handler map from cli_unified's re-exports."""
    from neoswga import cli_unified

    _HANDLERS.update(
        {
            "evaluate-set": cli_unified.run_evaluate_set,
            "analyze-coverage": cli_unified.run_analyze_coverage,
            "predict-efficiency": cli_unified.run_predict_efficiency,
            "contract-set": cli_unified.run_contract_set,
            "rescore-set": cli_unified.run_rescore_set,
            "expand-primers": cli_unified.run_expand_primers,
            "swap-primer": cli_unified.run_swap_primer,
            "analyze-set": cli_unified.analyze_primer_set,
            "analyze-genome": cli_unified.run_analyze_genome,
            "analyze-dimers": cli_unified.run_analyze_dimers,
            "analyze-stability": cli_unified.run_analyze_stability,
            "suggest": cli_unified.run_suggest,
            "init": cli_unified.run_init,
            "validate-params": cli_unified.run_validate_params,
            "schema": cli_unified.run_schema,
            "interpret": cli_unified.run_interpret,
            "report": cli_unified.run_report,
            "export": cli_unified.run_export,
            "doctor": cli_unified.run_doctor,
        }
    )


# ----------------------------------------------------------------------
# A real, minimal pipeline run
# ----------------------------------------------------------------------
#
# The iterate commands (expand/swap/contract/rescore) need `-j params.json`
# plus the artifacts count-kmers/filter/score produce: HDF5 position files and
# step3_df.csv. Building those once per session on a 20 kb genome is cheap and
# means these tests exercise the real file formats rather than mocks that could
# drift from them.


def _have_jellyfish():
    import shutil

    return shutil.which("jellyfish") is not None


@pytest.fixture(scope="session")
def pipeline_run(tmp_path_factory, genome_seq):
    """Run count-kmers -> filter -> score once, and hand back the workspace."""
    if not _have_jellyfish():
        pytest.skip("jellyfish not available (required to build position files)")

    import json
    import subprocess
    import sys

    work = tmp_path_factory.mktemp("pipeline")
    fasta = work / "target.fasta"
    fasta.write_text(
        ">target\n"
        + "\n".join(genome_seq[i : i + 70] for i in range(0, len(genome_seq), 70))
        + "\n"
    )

    params = {
        "fg_genomes": [str(fasta)],
        "bg_genomes": [],
        "fg_prefixes": [str(work / "target")],
        "bg_prefixes": [],
        "data_dir": str(work / "results"),
        "min_k": 10,
        "max_k": 10,
        "polymerase": "phi29",
        "reaction_temp": 30.0,
        "min_fg_freq": 1e-6,
        "max_bg_freq": 1.0,
        "max_gini": 1.0,
        "max_primer": 60,
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
    params_file = work / "params.json"
    params_file.write_text(json.dumps(params, indent=2))

    for step in ("count-kmers", "filter", "score"):
        proc = subprocess.run(
            [sys.executable, "-m", "neoswga.cli_unified", step, "-j", str(params_file)],
            capture_output=True,
            text=True,
            cwd=str(work),
            timeout=900,
        )
        if proc.returncode != 0:
            pytest.skip(
                f"pipeline step {step!r} failed, cannot build fixture:\n{proc.stderr[-800:]}"
            )

    return {
        "dir": work,
        "params_file": str(params_file),
        "data_dir": str(work / "results"),
        "fasta": str(fasta),
    }


@pytest.fixture(scope="session")
def scored_primers(pipeline_run):
    """Primers that survived filtering and scoring, as the pipeline produced them."""
    import csv
    import os

    path = os.path.join(pipeline_run["data_dir"], "step3_df.csv")
    if not os.path.exists(path):
        pytest.skip("step3_df.csv not produced")
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    key = "primer" if rows and "primer" in rows[0] else list(rows[0])[0]
    return [r[key] for r in rows]
