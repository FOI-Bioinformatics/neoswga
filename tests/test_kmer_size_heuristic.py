"""Tests for the genome-size + GC + polymerase -> k-mer range heuristic.

This heuristic surfaces as both:
- `condition_suggester.suggest_kmer_range()` -- pure function
- A warning in `get_params()` when the user's explicit range looks wrong
  for their genome size
"""

import json
import logging
from types import SimpleNamespace

import pytest

from neoswga.core.condition_suggester import suggest_kmer_range


# ----- Pure function ---------------------------------------------------------

class TestSuggestKmerRange:
    def test_phi29_small_plasmid(self):
        # 5 kb plasmid, balanced GC
        min_k, max_k = suggest_kmer_range(5_000, gc=0.5, polymerase='phi29')
        assert min_k == 6
        assert max_k <= 12

    def test_phi29_typical_bacterium(self):
        # 2 Mb bacterium, balanced GC
        min_k, max_k = suggest_kmer_range(2_000_000, gc=0.5, polymerase='phi29')
        assert min_k >= 8, "1 Mb+ bacteria need longer min_k for specificity"
        assert max_k == 12

    def test_phi29_large_bacterium(self):
        # 10 Mb
        min_k, max_k = suggest_kmer_range(10_000_000, gc=0.5, polymerase='phi29')
        assert min_k >= 10

    def test_phi29_human_scale_background(self):
        min_k, max_k = suggest_kmer_range(3_000_000_000, gc=0.41, polymerase='phi29')
        assert min_k >= 12
        assert max_k >= min_k + 4

    def test_equiphi29_bacterium(self):
        min_k, max_k = suggest_kmer_range(2_000_000, gc=0.5, polymerase='equiphi29')
        assert min_k >= 10
        assert max_k == 18

    def test_bst_lamp_range(self):
        min_k, max_k = suggest_kmer_range(2_000_000, gc=0.5, polymerase='bst')
        assert min_k >= 15
        assert max_k >= 25

    def test_extreme_at_nudges_min_k_up(self):
        # 20% GC, 2 Mb
        normal = suggest_kmer_range(2_000_000, gc=0.5, polymerase='phi29')
        at_rich = suggest_kmer_range(2_000_000, gc=0.2, polymerase='phi29')
        assert at_rich[0] >= normal[0]

    def test_extreme_gc_nudges_min_k_up(self):
        normal = suggest_kmer_range(2_000_000, gc=0.5, polymerase='phi29')
        gc_rich = suggest_kmer_range(2_000_000, gc=0.72, polymerase='phi29')
        assert gc_rich[0] >= normal[0]

    def test_unknown_polymerase_fallback(self):
        # Defensive fallback
        min_k, max_k = suggest_kmer_range(2_000_000, polymerase='unknown_enzyme')
        assert 6 <= min_k <= max_k <= 12


# ----- get_params warning ----------------------------------------------------

def _mk_args(json_file):
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


def _write_params(tmp_path, seq: str, extra=None):
    fg_fasta = tmp_path / "fg.fna"
    fg_fasta.write_text(">seq\n" + seq + "\n")
    params = {
        "data_dir": str(tmp_path),
        "src_dir": str(tmp_path),
        "fg_genomes": [str(fg_fasta)],
        "fg_prefixes": [str(tmp_path / "fg")],
        "fg_seq_lengths": [len(seq)],
        "cpus": 1,
    }
    if extra:
        params.update(extra)
    path = tmp_path / "params.json"
    path.write_text(json.dumps(params))
    return path


def test_large_bacterium_with_short_max_k_warns(tmp_path, caplog):
    """5 Mb genome with max_k=7 (below suggested 8) should warn."""
    from neoswga.core import parameter
    parameter.reset_to_defaults()

    # Claim a 5 Mb genome in params.json (actual FASTA is smaller but
    # fg_seq_lengths takes precedence for this test)
    seq = "ACGT" * 100  # 400 bp actual
    params_file = _write_params(tmp_path, seq, extra={
        "polymerase": "phi29",
        "min_k": 6,
        "max_k": 7,
        "fg_seq_lengths": [5_000_000],
    })
    with caplog.at_level(logging.WARNING, logger="neoswga.core.parameter"):
        parameter.get_params(_mk_args(params_file))

    # Ensure a size-related k-mer warning was emitted
    assert any("K-mer specificity" in rec.message for rec in caplog.records), (
        f"Expected size-mismatch warning; got: {[r.message for r in caplog.records]}"
    )


def test_tiny_plasmid_with_long_min_k_warns(tmp_path, caplog):
    """5 kb plasmid with min_k=15 produces a warning."""
    from neoswga.core import parameter
    parameter.reset_to_defaults()

    seq = "ACGT" * 100
    params_file = _write_params(tmp_path, seq, extra={
        "polymerase": "bst",  # whose defaults would be 15-25 but is obviously wrong for plasmid
        "min_k": 15,
        "max_k": 18,
        "fg_seq_lengths": [5_000],
    })
    with caplog.at_level(logging.WARNING, logger="neoswga.core.parameter"):
        parameter.get_params(_mk_args(params_file))

    assert any("Too-long primers" in rec.message for rec in caplog.records)
