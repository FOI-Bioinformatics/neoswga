"""Tests for auto-enabled adaptive GC filtering in get_params().

The adaptive filter engages automatically when the target genome has non-standard
GC content. Explicit user values always win: an explicit gc_min/gc_max or
"adaptive_gc": false disables the auto behavior.
"""

import json
from types import SimpleNamespace

import pytest


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


def test_at_rich_genome_triggers_adaptive_window(tmp_path, caplog):
    """25% GC genome should get gc_min ~0.10 (clamped to 0.15) and gc_max ~0.40."""
    from neoswga.core import parameter
    parameter.reset_to_defaults()

    # 25% GC sequence
    seq = "A" * 750 + "G" * 250  # 25% GC
    params_file = _write_params(tmp_path, seq)

    import logging
    with caplog.at_level(logging.INFO, logger="neoswga.core.parameter"):
        parameter.get_params(_mk_args(params_file))

    assert 0.14 <= parameter.gc_min <= 0.16
    assert 0.39 <= parameter.gc_max <= 0.41
    # Extreme-GC message surfaces to the user
    assert any("Extreme GC genome" in rec.message for rec in caplog.records)


def test_gc_rich_genome_triggers_adaptive_window(tmp_path, caplog):
    """70% GC genome should get gc_min ~0.55 and gc_max ~0.85."""
    from neoswga.core import parameter
    parameter.reset_to_defaults()

    # 70% GC sequence
    seq = "G" * 700 + "A" * 300
    params_file = _write_params(tmp_path, seq)

    import logging
    with caplog.at_level(logging.INFO, logger="neoswga.core.parameter"):
        parameter.get_params(_mk_args(params_file))

    assert 0.54 <= parameter.gc_min <= 0.56
    assert 0.84 <= parameter.gc_max <= 0.86
    assert any("Extreme GC genome" in rec.message for rec in caplog.records)


def test_explicit_gc_min_overrides_adaptive(tmp_path):
    """User-specified gc_min in params.json wins over the adaptive value."""
    from neoswga.core import parameter
    parameter.reset_to_defaults()

    # 25% GC, but user pins gc_min=0.30
    seq = "A" * 750 + "G" * 250
    params_file = _write_params(tmp_path, seq, extra={"gc_min": 0.30, "gc_max": 0.55})
    parameter.get_params(_mk_args(params_file))

    assert parameter.gc_min == 0.30
    assert parameter.gc_max == 0.55


def test_adaptive_gc_false_disables_auto(tmp_path):
    """adaptive_gc: false in params.json disables auto-adaptive GC."""
    from neoswga.core import parameter
    parameter.reset_to_defaults()

    # 70% GC, but user opts out
    seq = "G" * 700 + "A" * 300
    params_file = _write_params(tmp_path, seq, extra={"adaptive_gc": False})
    parameter.get_params(_mk_args(params_file))

    # Should fall back to default gc_min/gc_max, NOT the adaptive 0.55-0.85
    assert parameter.gc_min == 0.375
    assert parameter.gc_max == 0.625


def test_balanced_gc_skips_extreme_message(tmp_path, caplog):
    """A 50% GC genome still gets adaptive bounds but no extreme-GC warning."""
    from neoswga.core import parameter
    parameter.reset_to_defaults()

    seq = "AG" * 500  # 50% GC
    params_file = _write_params(tmp_path, seq)

    import logging
    with caplog.at_level(logging.INFO, logger="neoswga.core.parameter"):
        parameter.get_params(_mk_args(params_file))

    assert not any("Extreme GC genome" in rec.message for rec in caplog.records)
