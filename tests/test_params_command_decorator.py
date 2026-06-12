"""Tests for the params_command CLI-handler preamble decorator."""

from types import SimpleNamespace

import pytest

from neoswga.cli._common import params_command


def test_validate_merge_seed_run_in_order(monkeypatch):
    calls = []

    monkeypatch.setattr(
        "neoswga.cli._common.validate_params_json_file",
        lambda p: calls.append(("validate", p)),
    )
    monkeypatch.setattr(
        "neoswga.cli._common.merge_args_to_parameter",
        lambda args, parameter, names: calls.append(("merge", tuple(names))),
    )
    monkeypatch.setattr("neoswga.cli._common._apply_seed", lambda args: calls.append(("seed",)))

    @params_command(seed=True)
    def handler(args):
        calls.append(("body",))
        return "result"

    out = handler(SimpleNamespace(json_file="params.json"))

    assert out == "result"
    # validate -> merge(json_file) -> seed -> body, in that order.
    assert calls == [
        ("validate", "params.json"),
        ("merge", ("json_file",)),
        ("seed",),
        ("body",),
    ]


def test_merge_none_skips_merge_and_seed_off_by_default(monkeypatch):
    calls = []
    monkeypatch.setattr(
        "neoswga.cli._common.validate_params_json_file",
        lambda p: calls.append("validate"),
    )
    monkeypatch.setattr(
        "neoswga.cli._common.merge_args_to_parameter",
        lambda *a, **k: calls.append("merge"),
    )
    monkeypatch.setattr("neoswga.cli._common._apply_seed", lambda args: calls.append("seed"))

    @params_command(merge=None)
    def handler(args):
        calls.append("body")

    handler(SimpleNamespace(json_file="p.json"))
    assert calls == ["validate", "body"]  # no merge, no seed


def test_validate_skipped_when_no_json_file(monkeypatch):
    calls = []
    monkeypatch.setattr(
        "neoswga.cli._common.validate_params_json_file",
        lambda p: calls.append("validate"),
    )
    monkeypatch.setattr(
        "neoswga.cli._common.merge_args_to_parameter", lambda *a, **k: calls.append("merge")
    )

    @params_command(validate=True, merge=None)
    def handler(args):
        calls.append("body")

    handler(SimpleNamespace())  # no json_file attribute
    assert calls == ["body"]  # validate skipped (no json_file)


def test_preserves_function_metadata():
    @params_command()
    def run_something(args):
        """Docstring stays."""

    assert run_something.__name__ == "run_something"
    assert run_something.__doc__ == "Docstring stays."


def test_bare_decorator_usage(monkeypatch):
    """@params_command without parentheses should also work."""
    calls = []
    monkeypatch.setattr(
        "neoswga.cli._common.validate_params_json_file", lambda p: calls.append("validate")
    )
    monkeypatch.setattr(
        "neoswga.cli._common.merge_args_to_parameter", lambda *a, **k: calls.append("merge")
    )

    @params_command
    def handler(args):
        calls.append("body")

    handler(SimpleNamespace(json_file="p.json"))
    assert calls == ["validate", "merge", "body"]
