"""Integration tests for blacklist handling in params.json loading.

Verifies that bl_seq_lengths is auto-computed when missing or mismatched,
so the frequency-based blacklist filter in pipeline._filter_blacklist_penalty
gets correct denominators whether bl_genomes came from params.json or
--blacklist CLI.
"""

import json
from types import SimpleNamespace

import pytest


def _write_fasta(path, seq: str, header: str = "seq"):
    path.write_text(f">{header}\n{seq}\n")


def _mk_args(json_file):
    """Build an argparse.Namespace with the fields get_params() reads,
    all set to None so params.json values / defaults take effect."""
    attrs = {
        "json_file": str(json_file),
        "data_dir": None, "src_dir": None,
        "min_fg_freq": None, "max_bg_freq": None,
        "min_tm": None, "max_tm": None,
        "max_gini": None, "max_primer": None, "min_amp_pred": None,
        "cpus": None, "max_dimer_bp": None, "max_self_dimer_bp": None,
        "verbose": None, "drop_iterations": None, "iterations": None,
        "top_set_count": None, "retries": None, "max_sets": None,
        "selection_metric": None, "fg_circular": None, "bg_circular": None,
        "min_k": None, "max_k": None,
        "fasta_fore": None, "fasta_back": None,
        "kmer_fore": None, "kmer_back": None,
    }
    return SimpleNamespace(**attrs)


def _minimal_params(tmp_path, bl_fasta=None, bl_lengths=None):
    """Create minimal params.json; omit bl_seq_lengths intentionally."""
    fg_fasta = tmp_path / "fg.fna"
    _write_fasta(fg_fasta, "ACGT" * 250)  # 1000 bp

    params = {
        "data_dir": str(tmp_path),
        "src_dir": str(tmp_path),
        "fg_genomes": [str(fg_fasta)],
        "fg_prefixes": [str(tmp_path / "fg")],
        "fg_seq_lengths": [1000],
        "cpus": 1,
    }
    if bl_fasta is not None:
        params["bl_genomes"] = [str(bl_fasta)]
        params["bl_prefixes"] = [str(tmp_path / "bl")]
        if bl_lengths is not None:
            params["bl_seq_lengths"] = bl_lengths

    params_file = tmp_path / "params.json"
    params_file.write_text(json.dumps(params))
    return params_file


def test_bl_seq_lengths_auto_computed_when_missing(tmp_path):
    """When bl_genomes is present but bl_seq_lengths is not, it should be auto-computed."""
    from neoswga.core import parameter

    bl_fasta = tmp_path / "bl.fna"
    _write_fasta(bl_fasta, "ACGT" * 500)  # 2000 bp

    params_file = _minimal_params(tmp_path, bl_fasta=bl_fasta, bl_lengths=None)

    parameter.reset_to_defaults()
    parameter.get_params(_mk_args(params_file))

    assert parameter.bl_seq_lengths == [2000], (
        f"bl_seq_lengths should be auto-computed; got {parameter.bl_seq_lengths}"
    )


def test_bl_seq_lengths_recomputed_when_mismatched(tmp_path):
    """If bl_seq_lengths count mismatches bl_genomes, it should be recomputed."""
    from neoswga.core import parameter

    bl_fasta = tmp_path / "bl.fna"
    _write_fasta(bl_fasta, "ACGT" * 250)  # 1000 bp

    # Intentionally wrong: two lengths for one genome
    params_file = _minimal_params(tmp_path, bl_fasta=bl_fasta, bl_lengths=[123, 456])

    parameter.reset_to_defaults()
    parameter.get_params(_mk_args(params_file))

    assert parameter.bl_seq_lengths == [1000], (
        f"Mismatched bl_seq_lengths should be recomputed; got {parameter.bl_seq_lengths}"
    )


def test_bl_seq_lengths_preserved_when_correct(tmp_path):
    """If bl_seq_lengths matches bl_genomes count, it should be preserved."""
    from neoswga.core import parameter

    bl_fasta = tmp_path / "bl.fna"
    _write_fasta(bl_fasta, "ACGT" * 250)  # actual 1000

    # User explicitly provides a plausible value (e.g., from a prior run)
    params_file = _minimal_params(tmp_path, bl_fasta=bl_fasta, bl_lengths=[1000])

    parameter.reset_to_defaults()
    parameter.get_params(_mk_args(params_file))

    assert parameter.bl_seq_lengths == [1000]
