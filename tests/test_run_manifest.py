"""Tests for the run_manifest reproducibility writer."""

import json
import os

import pytest

from neoswga.core import run_manifest as rm


def test_write_manifest_creates_file(tmp_path):
    out = rm.write_manifest(
        step="count-kmers",
        data_dir=str(tmp_path),
        params_path=None,
        input_files=None,
        seed=42,
    )

    assert out is not None
    assert os.path.basename(out) == rm.MANIFEST_FILENAME

    with open(out) as f:
        data = json.load(f)

    assert data["steps"][0]["step"] == "count-kmers"
    assert data["steps"][0]["seed"] == 42
    assert data["steps"][0]["neoswga_version"]
    assert "timestamp_utc" in data["steps"][0]
    assert "cli_invocation" in data["steps"][0]


def test_write_manifest_appends_across_steps(tmp_path):
    rm.write_manifest(step="count-kmers", data_dir=str(tmp_path))
    rm.write_manifest(step="filter", data_dir=str(tmp_path))
    rm.write_manifest(step="score", data_dir=str(tmp_path))

    with open(tmp_path / rm.MANIFEST_FILENAME) as f:
        data = json.load(f)

    assert [s["step"] for s in data["steps"]] == ["count-kmers", "filter", "score"]


def test_write_manifest_checksums_inputs(tmp_path):
    fasta = tmp_path / "genome.fna"
    fasta.write_text(">seq\nACGTACGT\n")

    out = rm.write_manifest(
        step="count-kmers",
        data_dir=str(tmp_path),
        input_files=[str(fasta)],
    )

    with open(out) as f:
        data = json.load(f)

    checksums = data["steps"][0]["input_checksums"]
    assert str(fasta) in checksums
    assert len(checksums[str(fasta)]) == 64  # sha256 hex


def test_write_manifest_missing_input_silently_skipped(tmp_path):
    out = rm.write_manifest(
        step="filter",
        data_dir=str(tmp_path),
        input_files=[str(tmp_path / "does_not_exist.csv")],
    )

    with open(out) as f:
        data = json.load(f)

    assert data["steps"][0]["input_checksums"] == {}


def test_write_manifest_handles_corrupt_existing_file(tmp_path):
    (tmp_path / rm.MANIFEST_FILENAME).write_text("not valid json {{{")

    out = rm.write_manifest(step="count-kmers", data_dir=str(tmp_path))
    assert out is not None

    with open(out) as f:
        data = json.load(f)

    assert len(data["steps"]) == 1
    assert data["steps"][0]["step"] == "count-kmers"


def test_write_manifest_loads_params_from_path(tmp_path):
    params_path = tmp_path / "params.json"
    params_path.write_text(json.dumps({"min_k": 12, "max_k": 18}))

    out = rm.write_manifest(
        step="count-kmers",
        data_dir=str(tmp_path),
        params_path=str(params_path),
    )

    with open(out) as f:
        data = json.load(f)

    assert data["steps"][0]["resolved_params"] == {"min_k": 12, "max_k": 18}
    assert data["steps"][0]["params_path"] == str(params_path)


def test_write_manifest_returns_none_without_data_dir():
    assert rm.write_manifest(step="filter", data_dir=None) is None
    assert rm.write_manifest(step="filter", data_dir="") is None


def test_reproducibility_fingerprint_stable(tmp_path):
    """Same seed + same params + same inputs -> identical fingerprint fields."""
    fasta = tmp_path / "g.fna"
    fasta.write_text(">x\nAAACCC\n")
    params = tmp_path / "p.json"
    params.write_text(json.dumps({"min_k": 6}))

    runs = []
    for _ in range(2):
        out_dir = tmp_path / f"run_{_}"
        rm.write_manifest(
            step="count-kmers",
            data_dir=str(out_dir),
            params_path=str(params),
            input_files=[str(fasta)],
            seed=7,
        )
        with open(out_dir / rm.MANIFEST_FILENAME) as f:
            runs.append(json.load(f)["steps"][0])

    # Fingerprint = the fields that must be identical for reproducibility
    def fingerprint(entry):
        return {
            "seed": entry["seed"],
            "resolved_params": entry["resolved_params"],
            "input_checksums": entry["input_checksums"],
            "neoswga_version": entry["neoswga_version"],
        }

    assert fingerprint(runs[0]) == fingerprint(runs[1])


# ---------------------------------------------------------------------------
# `resolved_params` is a verbatim copy of params.json, so the buffer the run
# actually used was never recorded. The GC-adaptive strategy sets betaine and
# DMSO at run time from the target's GC, and neither appears in params.json --
# so a reader of the manifest, or `export`/`report` reconstructing conditions
# from it, gets a different reaction than the one that produced the design.
# ---------------------------------------------------------------------------


def test_manifest_records_the_conditions_the_run_used(tmp_path):
    out = rm.write_manifest(
        step="filter",
        data_dir=str(tmp_path),
        effective_conditions={"polymerase": "equiphi29", "reaction_temp": 42.0, "betaine_m": 2.0},
    )

    with open(out) as f:
        entry = json.load(f)["steps"][0]

    assert entry["effective_conditions"]["betaine_m"] == 2.0
    assert entry["effective_conditions"]["polymerase"] == "equiphi29"


def test_effective_conditions_absent_when_not_supplied(tmp_path):
    out = rm.write_manifest(step="filter", data_dir=str(tmp_path))

    with open(out) as f:
        entry = json.load(f)["steps"][0]

    assert entry["effective_conditions"] is None


def test_read_effective_conditions_returns_the_latest(tmp_path):
    rm.write_manifest(
        step="filter",
        data_dir=str(tmp_path),
        effective_conditions={"polymerase": "phi29", "betaine_m": 0.0},
    )
    rm.write_manifest(
        step="optimize",
        data_dir=str(tmp_path),
        effective_conditions={"polymerase": "equiphi29", "betaine_m": 2.0},
    )

    assert rm.read_effective_conditions(str(tmp_path)) == {
        "polymerase": "equiphi29",
        "betaine_m": 2.0,
    }


def test_read_effective_conditions_without_manifest_is_none(tmp_path):
    assert rm.read_effective_conditions(str(tmp_path)) is None


def test_record_run_manifest_snapshots_the_effective_parameter_state(tmp_path):
    """The betaine the GC-adaptive strategy chose is in the parameter module,
    not in params.json, and this is the only place holding both."""
    from types import SimpleNamespace

    from neoswga.cli._common import _record_run_manifest

    parameter = SimpleNamespace(
        data_dir=str(tmp_path),
        polymerase="equiphi29",
        reaction_temp=42.0,
        na_conc=50.0,
        mg_conc=10.0,
        betaine_m=2.0,
        dmso_percent=0.0,
    )
    _record_run_manifest("filter", SimpleNamespace(json_file=None, seed=7), parameter)

    conditions = rm.read_effective_conditions(str(tmp_path))
    assert conditions is not None
    assert conditions["betaine_m"] == 2.0
    assert conditions["polymerase"] == "equiphi29"
    assert conditions["reaction_temp"] == 42.0
    assert conditions["mg_conc"] == 10.0
