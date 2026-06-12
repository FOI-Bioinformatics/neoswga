"""Coverage for neoswga.core.workflow_selector menu/dispatch helpers."""

import builtins

import pytest

from neoswga.core import workflow_selector as ws


def test_print_menu_renders(capsys):
    ws.print_menu("TITLE", [("Opt A", "does A"), ("Opt B", "does B")])
    out = capsys.readouterr().out
    assert "TITLE" in out and "Opt A" in out and "does B" in out
    assert "q." in out  # quit option shown by default


def test_get_choice_valid(monkeypatch):
    monkeypatch.setattr(builtins, "input", lambda _="": "2")
    assert ws.get_choice(3) == 2


def test_get_choice_quit(monkeypatch):
    monkeypatch.setattr(builtins, "input", lambda _="": "q")
    assert ws.get_choice(3) is None


def test_get_choice_retries_then_valid(monkeypatch):
    answers = iter(["0", "9", "abc", "1"])
    monkeypatch.setattr(builtins, "input", lambda _="": next(answers))
    assert ws.get_choice(3) == 1


def test_execute_command_success(monkeypatch):
    class _Result:
        returncode = 0

    monkeypatch.setattr(ws.subprocess, "run", lambda cmd: _Result())
    assert ws.execute_command(["neoswga", "--help"]) is True


def test_execute_command_failure(monkeypatch):
    class _Result:
        returncode = 3

    monkeypatch.setattr(ws.subprocess, "run", lambda cmd: _Result())
    assert ws.execute_command(["neoswga", "boom"]) is False


def test_execute_command_not_found(monkeypatch):
    def _raise(cmd):
        raise FileNotFoundError()

    monkeypatch.setattr(ws.subprocess, "run", _raise)
    assert ws.execute_command(["neoswga"]) is False


def test_prompt_execute_yes_runs(monkeypatch):
    monkeypatch.setattr(builtins, "input", lambda _="": "y")
    monkeypatch.setattr(ws, "execute_command", lambda cmd: True)
    assert ws.prompt_execute(["neoswga", "x"], "desc") is True


def test_prompt_execute_no_skips(monkeypatch):
    monkeypatch.setattr(builtins, "input", lambda _="": "n")
    called = {"ran": False}
    monkeypatch.setattr(ws, "execute_command", lambda cmd: called.__setitem__("ran", True))
    assert ws.prompt_execute(["neoswga", "x"]) is False
    assert called["ran"] is False


def test_run_workflow_selector_quits_immediately(monkeypatch):
    # 'q' at the main menu exits via SystemExit(0).
    monkeypatch.setattr(builtins, "input", lambda _="": "q")
    with pytest.raises(SystemExit) as exc:
        ws.run_workflow_selector()
    assert exc.value.code == 0
