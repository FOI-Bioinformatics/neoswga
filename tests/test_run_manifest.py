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
