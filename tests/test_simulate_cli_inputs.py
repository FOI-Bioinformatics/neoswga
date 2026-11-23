"""`simulate` runs the reaction it is asked about.

Three defects (2026-11-23 sweep, tier 2): the `SimulationConfig` carrying
`--duration` was built and never passed, so every run simulated the default
8 h; `--polymerase` defaulted to phi29, so a design optimized for equiphi29 at
42 C was simulated at 30 C; and primers were not upper-cased, so a lower-case
sequence found no sites and reported 0% coverage without an error.
"""

import json

import pytest

from neoswga.cli_unified import create_parser
from neoswga.core.run_manifest import write_manifest

PRIMER = "ACGTACGTAC"


@pytest.fixture
def genome(tmp_path):
    path = tmp_path / "g.fasta"
    path.write_text(">g\n" + ("TTTT" * 50 + PRIMER) * 4 + "\n")
    return path


@pytest.fixture
def captured(monkeypatch):
    seen = {}

    def _fake(**kwargs):
        seen.update(kwargs)
        return {
            "mean_coverage": 0.5,
            "std_coverage": 0.0,
            "mean_forks_created": 1.0,
            "mean_fork_travel": 1.0,
        }

    monkeypatch.setattr("neoswga.core.replication_simulator.simulate_primer_set", _fake)
    return seen


def _run(*argv):
    from neoswga.cli.simulate import run_simulate

    run_simulate(create_parser().parse_args(["simulate", *argv]))


def test_the_duration_reaches_the_simulation(tmp_path, genome, captured):
    _run("--primers", PRIMER, "--genome", str(genome), "-o", str(tmp_path / "o"), "--duration", "10")
    assert captured["config"].duration == 600.0


def test_the_recorded_reaction_is_simulated(tmp_path, genome, captured):
    results = tmp_path / "results"
    results.mkdir()
    (results / "step4_improved_df.csv").write_text(f"primer,set_index\n{PRIMER},0\n")
    write_manifest(
        "optimize",
        str(results),
        effective_conditions={"polymerase": "equiphi29", "reaction_temp": 42.0},
    )
    _run("--from-results", str(results), "--genome", str(genome), "-o", str(tmp_path / "o"))
    assert captured["conditions"].polymerase == "equiphi29"
    assert captured["conditions"].temp == 42.0


def test_an_explicit_polymerase_wins(tmp_path, genome, captured):
    _run(
        "--primers", PRIMER, "--genome", str(genome), "-o", str(tmp_path / "o"),
        "--polymerase", "bst",
    )
    assert captured["conditions"].polymerase == "bst"


def test_a_lower_case_primer_finds_its_sites(tmp_path, genome, captured):
    _run("--primers", PRIMER.lower(), "--genome", str(genome), "-o", str(tmp_path / "o"))
    positions = captured["primer_positions"][PRIMER]
    assert len(positions["forward"]) + len(positions["reverse"]) > 0


def test_a_non_acgt_primer_is_refused(tmp_path, genome, captured):
    with pytest.raises(SystemExit):
        _run("--primers", "ACGTXX", "--genome", str(genome), "-o", str(tmp_path / "o"))
    assert not captured
