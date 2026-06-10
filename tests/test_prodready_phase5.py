"""Phase 5 (production-readiness v2): UX, --seed surface, manifest provenance."""

import pytest


def _subparsers():
    from neoswga.cli_unified import create_parser

    parser = create_parser()
    return next(a for a in parser._actions if a.__class__.__name__ == "_SubParsersAction")


@pytest.mark.parametrize(
    "cmd", ["expand-primers", "swap-primer", "contract-set", "multi-genome", "simulate"]
)
def test_seed_flag_present_on_set_producing_commands(cmd):
    sub = _subparsers()
    opts = [o for a in sub.choices[cmd]._actions for o in a.option_strings]
    assert "--seed" in opts


def test_expand_primers_accepts_ensemble():
    sub = _subparsers()
    action = next(
        a
        for a in sub.choices["expand-primers"]._actions
        if "--optimization-method" in a.option_strings
    )
    assert "ensemble" in action.choices


def test_apply_seed_is_reproducible():
    import random

    from neoswga import cli_unified

    class _Args:
        seed = 123

    cli_unified._apply_seed(_Args())
    a = [random.random() for _ in range(5)]
    cli_unified._apply_seed(_Args())
    b = [random.random() for _ in range(5)]
    assert a == b


def test_record_run_manifest_uses_args_seed(monkeypatch):
    from neoswga import cli_unified

    captured = {}

    def fake_write_manifest(**kwargs):
        captured.update(kwargs)

    import neoswga.core.run_manifest as rm

    monkeypatch.setattr(rm, "write_manifest", fake_write_manifest)

    class _Args:
        seed = 777
        json_file = "params.json"

    class _Param:
        data_dir = "."
        seed = None  # the bug: previously read from here -> None

    cli_unified._record_run_manifest("optimize", _Args(), _Param(), input_files=["x.csv"])
    assert captured["seed"] == 777  # comes from args, not parameter


def test_seed_help_text_not_stale():
    """--seed help must not reference the removed genetic/moea optimizers."""
    sub = _subparsers()
    action = next(a for a in sub.choices["optimize"]._actions if "--seed" in a.option_strings)
    assert "genetic" not in action.help and "moea" not in action.help
