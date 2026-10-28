"""The documented `fg_circular` default (True) must survive resolution.

`params.schema.json` and `docs/params-reference.md` both give `fg_circular` a
default of True and `bg_circular` a default of False. Nothing injects the
schema's defaults at runtime, and `get_params` assigns `data["fg_circular"]`
from `get_value_or_default`, which returns None for an absent key -- so the
"fill absent fields" blocks cannot reach it afterwards. Left unfixed, an
omitted `fg_circular` resolved to None, `bool(None)` is False, and a circular
target (plasmid, bacterial chromosome -- the common SWGA case) silently went
LINEAR: no origin-spanning coverage wrap in the network/coverage geometry, and
no origin-spanning gap merge in `expand-primers --bam`. All shipped examples
and the wizard write the key explicitly, so the defect was latent there; it bit
a hand-written params.json that trusted the documented default.

These tests pin that an OMITTED key takes the documented default while an
EXPLICIT true/false is preserved untouched.
"""

import json
from types import SimpleNamespace

from neoswga.core import parameter


def _mk_args(json_file):
    # Every CLI-overridable key as None: the argparse state for a run that
    # passed only -j, so resolution falls entirely to the file (or its
    # default). fg_circular/bg_circular have no CLI flag at all, so None here
    # is exactly what the real command produces.
    return SimpleNamespace(
        json_file=str(json_file),
        data_dir=None,
        src_dir=None,
        min_fg_freq=None,
        max_bg_freq=None,
        min_tm=None,
        max_tm=None,
        max_gini=None,
        max_primer=None,
        min_amp_pred=None,
        cpus=None,
        max_dimer_bp=None,
        max_self_dimer_bp=None,
        verbose=None,
        drop_iterations=None,
        iterations=None,
        top_set_count=None,
        retries=None,
        max_sets=None,
        selection_metric=None,
        fg_circular=None,
        bg_circular=None,
        min_k=None,
        max_k=None,
        fasta_fore=None,
        fasta_back=None,
        kmer_fore=None,
        kmer_back=None,
    )


def _write_params(tmp_path, extra=None):
    fg_fasta = tmp_path / "fg.fna"
    fg_fasta.write_text(">seq\n" + ("ACGT" * 50) + "\n")
    params = {
        "data_dir": str(tmp_path),
        "src_dir": str(tmp_path),
        "fg_genomes": [str(fg_fasta)],
        "fg_prefixes": [str(tmp_path / "fg")],
        "fg_seq_lengths": [200],
        "cpus": 1,
    }
    if extra:
        params.update(extra)
    path = tmp_path / "params.json"
    path.write_text(json.dumps(params))
    return path


def test_omitted_fg_circular_resolves_to_documented_true(tmp_path):
    """An absent fg_circular must resolve to True, not None/False.

    This is the regression: `bool(None)` is False, so the omitted case used to
    make a circular target linear everywhere circular is consumed by truthiness
    (coverage wrap, `expand-primers --bam` origin-spanning gap merge).
    """
    params_file = _write_params(tmp_path)  # fg_circular absent
    resolved = parameter.get_params(_mk_args(params_file))
    assert resolved["fg_circular"] is True
    # iterate.py / pipeline.py read the module global; it must agree and be truthy.
    assert parameter.fg_circular is True
    assert bool(parameter.fg_circular) is True


def test_explicit_fg_circular_false_is_preserved(tmp_path):
    """A user who writes fg_circular: false means linear; the default must not
    clobber that back to True."""
    params_file = _write_params(tmp_path, extra={"fg_circular": False})
    resolved = parameter.get_params(_mk_args(params_file))
    assert resolved["fg_circular"] is False
    assert parameter.fg_circular is False


def test_explicit_fg_circular_true_is_preserved(tmp_path):
    params_file = _write_params(tmp_path, extra={"fg_circular": True})
    resolved = parameter.get_params(_mk_args(params_file))
    assert resolved["fg_circular"] is True
    assert parameter.fg_circular is True


def test_omitted_bg_circular_resolves_to_documented_false(tmp_path):
    """bg_circular's documented default is False. Its omitted value was already
    falsy (None), so this only pins that the stored value is the real default
    rather than None -- and guards against a future copy-paste that gives it
    True by symmetry with fg_circular."""
    params_file = _write_params(tmp_path)  # bg_circular absent
    resolved = parameter.get_params(_mk_args(params_file))
    assert resolved["bg_circular"] is False
    assert parameter.bg_circular is False
