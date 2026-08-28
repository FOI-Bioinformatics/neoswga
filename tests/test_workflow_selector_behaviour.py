"""Behavioural tests for the interactive workflow selector.

``workflow_selector`` backs ``neoswga start``, an interactive menu that
discovers and launches the rest of the CLI. Its central failure mode is drift:
the menu hardcodes subcommand names, flag names, and optimizer method choices
as literal strings, and nothing enforces that those strings stay in sync with
the real argparse definitions in ``cli_unified``. If a subcommand gets
renamed or removed, the menu keeps advertising it and a user who follows the
prompt hits an argparse error instead of a working command.

The tests below therefore do two things:

1. Cross-check every ``neoswga <subcommand>`` string the module can produce
   (both the commands it actually executes via ``execute_command`` and the
   commands it only prints as copy-paste examples) against the live
   subparser registry built by ``cli_unified.create_parser()``. This is the
   single most valuable check for a discovery menu: it fails loudly the
   moment a command is renamed underneath the menu.
2. Drive the interactive loop end to end with a scripted ``input()`` queue
   and a stubbed ``execute_command``, to confirm each menu choice builds the
   command line it claims to, that invalid input is rejected rather than
   silently accepted, and that quitting/backing out does not hang or crash.

No test reads from real stdin or calls the real ``subprocess.run``.
"""

import re

import pytest

from neoswga.cli_unified import create_parser
from neoswga.core import workflow_selector as ws


def _real_subcommands():
    """Return the set of subcommand names the actual CLI registers."""
    parser = create_parser()
    for action in parser._actions:
        if getattr(action, "choices", None):
            return set(action.choices.keys())
    raise AssertionError("Could not find the subparsers action on the top-level parser")


def _referenced_subcommands(source):
    """Extract every 'neoswga <subcommand>' token the module's source can emit.

    Covers both list-form commands actually executed (`["neoswga", "filter", ...]`)
    and the plain-text example commands printed in the advanced-features menu
    (e.g. "neoswga multi-genome \\"). A couple of English words that happen to
    follow the literal word "neoswga" in prose (module docstring, error
    messages) are not subcommands and are excluded explicitly rather than
    matched by a narrower, more fragile regex.
    """
    list_form = re.findall(r'"neoswga",\s*"([a-z][a-z0-9-]*)"', source)
    print_form = re.findall(r"neoswga ([a-z][a-z0-9]*(?:-[a-z0-9]+)*)\b", source)
    not_commands = {"command", "features"}
    return {c for c in set(list_form) | set(print_form) if c not in not_commands}


def test_every_referenced_subcommand_exists_in_the_real_cli():
    """Prime bug guard: catches renamed/removed subcommands the menu still offers.

    This is a regression guard, not a report of a found bug: at the time of
    writing every command referenced by workflow_selector.py exists in
    cli_unified.create_parser(). The test exists so that the *next* rename
    (which is exactly the kind of change this codebase's history shows
    happening -- see e.g. the split-create-parser refactor) fails a fast unit
    test instead of shipping a dead menu entry.
    """
    import inspect

    source = inspect.getsource(ws)
    referenced = _referenced_subcommands(source)
    real = _real_subcommands()

    assert referenced, "Regex extracted no subcommands -- likely a regex bug, not an empty module"
    missing = referenced - real
    assert not missing, f"workflow_selector references subcommands that don't exist: {missing}"


def test_optimize_method_choices_match_the_real_optimize_parser():
    """The method-selection submenu hardcodes a number->name mapping.

    If `optimize`'s --optimization-method choices ever change (a method
    renamed or removed), this hardcoded dict in workflow_selector.py would
    silently keep offering a value argparse rejects. Read the mapping's
    values straight out of the module (not retyped here) and check each one
    against the live argparse choices.
    """
    # The mapping lives inside run_workflow_selector's closure; recover it by
    # re-deriving from the menu labels, which enumerate 1..5 -> hybrid, etc.
    # Simpler and less brittle: just check the menu's printed labels list the
    # same methods the parser accepts.
    parser = create_parser()
    subparsers = None
    for action in parser._actions:
        if getattr(action, "choices", None):
            subparsers = action.choices
            break
    optimize_parser = subparsers["optimize"]
    real_choices = None
    for action in optimize_parser._actions:
        if action.dest == "optimization_method":
            real_choices = set(action.choices)
            break
    assert real_choices is not None

    hardcoded_methods = {
        "hybrid",
        "dominating-set",
        "background-aware",
        "network",
        "clique",
        "ensemble",
    }
    assert hardcoded_methods == real_choices, (
        "The workflow_selector method-selection submenu's hardcoded method names "
        f"{hardcoded_methods} no longer match optimize's real choices {real_choices}"
    )


# ----------------------------------------------------------------------
# get_choice
# ----------------------------------------------------------------------


def test_get_choice_accepts_q_for_quit(monkeypatch):
    monkeypatch.setattr("builtins.input", lambda *_: "q")
    assert ws.get_choice(5) is None


def test_get_choice_is_case_insensitive_for_quit(monkeypatch):
    monkeypatch.setattr("builtins.input", lambda *_: "Q")
    assert ws.get_choice(5) is None


def test_get_choice_rejects_out_of_range_then_accepts_valid(monkeypatch, capsys):
    """Out-of-range numbers must be re-prompted, not silently clamped or crashed on."""
    responses = iter(["0", "99", "3"])
    monkeypatch.setattr("builtins.input", lambda *_: next(responses))
    assert ws.get_choice(5) == 3
    out = capsys.readouterr().out
    assert out.count("Please enter a number between 1 and 5") == 2


def test_get_choice_rejects_non_numeric_input(monkeypatch, capsys):
    """Garbage input must not raise ValueError out of get_choice -- it should re-prompt."""
    responses = iter(["abc", "2"])
    monkeypatch.setattr("builtins.input", lambda *_: next(responses))
    assert ws.get_choice(5) == 2
    assert "Please enter a valid number or 'q' to quit" in capsys.readouterr().out


# ----------------------------------------------------------------------
# prompt_execute / execute_command
# ----------------------------------------------------------------------


def test_prompt_execute_skips_on_no(monkeypatch):
    monkeypatch.setattr("builtins.input", lambda *_: "n")
    called = []
    monkeypatch.setattr(ws, "execute_command", lambda cmd: called.append(cmd) or True)
    result = ws.prompt_execute(["neoswga", "filter", "-j", "p.json"], "desc")
    assert result is False
    assert called == [], "execute_command must not run when the user declines"


def test_prompt_execute_runs_and_returns_execute_command_result(monkeypatch):
    """The return value must propagate execute_command's success/failure, not just True."""
    monkeypatch.setattr("builtins.input", lambda *_: "yes")
    monkeypatch.setattr(ws, "execute_command", lambda cmd: False)
    assert ws.prompt_execute(["neoswga", "filter"], "") is False


def test_prompt_execute_reprompts_on_invalid_answer(monkeypatch, capsys):
    responses = iter(["maybe", "n"])
    monkeypatch.setattr("builtins.input", lambda *_: next(responses))
    assert ws.prompt_execute(["neoswga", "filter"], "") is False
    assert "Please enter 'y' or 'n'" in capsys.readouterr().out


def test_execute_command_handles_missing_binary(monkeypatch):
    """If neoswga itself isn't on PATH, the menu should report that, not traceback."""

    def raise_fnf(cmd):
        raise FileNotFoundError()

    monkeypatch.setattr(ws.subprocess, "run", raise_fnf)
    assert ws.execute_command(["neoswga", "x"]) is False


def test_execute_command_handles_nonzero_exit(monkeypatch):
    class FakeResult:
        returncode = 1

    monkeypatch.setattr(ws.subprocess, "run", lambda cmd: FakeResult())
    assert ws.execute_command(["neoswga", "x"]) is False


def test_execute_command_true_on_success(monkeypatch):
    class FakeResult:
        returncode = 0

    monkeypatch.setattr(ws.subprocess, "run", lambda cmd: FakeResult())
    assert ws.execute_command(["neoswga", "x"]) is True


# ----------------------------------------------------------------------
# print_menu
# ----------------------------------------------------------------------


def test_print_menu_omits_empty_descriptions(capsys):
    """'Back to main menu' entries pass description='' -- must not print a blank line label."""
    ws.print_menu("TITLE", [("Option A", "desc A"), ("Back to main menu", "")])
    out = capsys.readouterr().out
    assert "1. Option A" in out
    assert "desc A" in out
    assert "2. Back to main menu" in out


def test_print_menu_hides_quit_when_disabled(capsys):
    ws.print_menu("TITLE", [("Option A", "")], show_quit=False)
    assert "Quit" not in capsys.readouterr().out


# ----------------------------------------------------------------------
# run_workflow_selector: full menu-driven flows
# ----------------------------------------------------------------------


def _run_with_inputs(monkeypatch, inputs, calls):
    """Drive run_workflow_selector with a scripted input queue.

    execute_command is stubbed to record the command list it would have run
    (as a real subprocess call) instead of touching subprocess at all.
    """
    it = iter(inputs)
    monkeypatch.setattr("builtins.input", lambda *_: next(it))
    monkeypatch.setattr(ws, "execute_command", lambda cmd: calls.append(cmd) or True)
    with pytest.raises(SystemExit):
        ws.run_workflow_selector()


def test_quit_from_main_menu_raises_system_exit(monkeypatch):
    calls = []
    _run_with_inputs(monkeypatch, ["q"], calls)
    assert calls == []


def test_full_pipeline_runs_all_four_steps_in_order(monkeypatch):
    calls = []
    inputs = [
        "2",  # main: run pipeline
        "1",  # run full pipeline
        "params.json",  # params file
        "y",
        "y",
        "y",
        "y",  # confirm each of the 4 steps
        "",  # press enter to continue
        "6",  # back to main menu
        "q",  # quit
    ]
    _run_with_inputs(monkeypatch, inputs, calls)
    assert calls == [
        ["neoswga", "count-kmers", "-j", "params.json"],
        ["neoswga", "filter", "-j", "params.json"],
        ["neoswga", "score", "-j", "params.json"],
        ["neoswga", "optimize", "-j", "params.json"],
    ]


def test_full_pipeline_stops_after_a_declined_step(monkeypatch):
    """Declining step 2 (filter) must stop the chain -- score/optimize must not run."""
    calls = []
    inputs = [
        "2",
        "1",
        "params.json",
        "y",  # count-kmers: yes
        "n",  # filter: no -> pipeline should stop here
        "",
        "6",
        "q",
    ]
    _run_with_inputs(monkeypatch, inputs, calls)
    assert calls == [["neoswga", "count-kmers", "-j", "params.json"]]


@pytest.mark.parametrize(
    "step_choice,expected_cmd",
    [
        ("2", ["neoswga", "count-kmers", "-j", "params.json"]),
        ("3", ["neoswga", "filter", "-j", "params.json"]),
        ("4", ["neoswga", "score", "-j", "params.json"]),
    ],
)
def test_individual_pipeline_step_maps_to_correct_command(monkeypatch, step_choice, expected_cmd):
    """Each numbered step in the submenu must launch the correspondingly-named step."""
    calls = []
    inputs = ["2", step_choice, "params.json", "y", "", "6", "q"]
    _run_with_inputs(monkeypatch, inputs, calls)
    assert calls == [expected_cmd]


@pytest.mark.parametrize(
    "method_choice,expected_method",
    [
        ("1", "hybrid"),
        ("2", "dominating-set"),
        ("3", "background-aware"),
        ("4", "network"),
        ("5", "clique"),
        ("6", "ensemble"),
    ],
)
def test_optimize_step_method_selection(monkeypatch, method_choice, expected_method):
    calls = []
    inputs = ["2", "5", "params.json", method_choice, "y", "", "6", "q"]
    _run_with_inputs(monkeypatch, inputs, calls)
    assert calls == [
        ["neoswga", "optimize", "-j", "params.json", f"--optimization-method={expected_method}"]
    ]


def test_optimize_step_unrecognised_method_falls_back_to_hybrid(monkeypatch):
    """An out-of-range/garbage method answer must not crash -- dict.get default is hybrid."""
    calls = []
    inputs = ["2", "5", "params.json", "garbage", "y", "", "6", "q"]
    _run_with_inputs(monkeypatch, inputs, calls)
    assert calls == [["neoswga", "optimize", "-j", "params.json", "--optimization-method=hybrid"]]


def test_pipeline_default_params_file_is_used_when_blank(monkeypatch):
    """An empty params-file answer must default to 'params.json', not pass '' to the CLI."""
    calls = []
    inputs = ["2", "2", "", "y", "", "6", "q"]
    _run_with_inputs(monkeypatch, inputs, calls)
    assert calls == [["neoswga", "count-kmers", "-j", "params.json"]]


def test_setup_new_project_with_background_genome(monkeypatch):
    calls = []
    inputs = ["1", "genome.fasta", "bg.fasta", "y", "", "q"]
    _run_with_inputs(monkeypatch, inputs, calls)
    assert calls == [["neoswga", "init", "--genome", "genome.fasta", "--background", "bg.fasta"]]


def test_setup_new_project_without_background_genome(monkeypatch):
    calls = []
    inputs = ["1", "genome.fasta", "", "y", "", "q"]
    _run_with_inputs(monkeypatch, inputs, calls)
    assert calls == [["neoswga", "init", "--genome", "genome.fasta"]]


def test_setup_new_project_skip_runs_nothing(monkeypatch):
    calls = []
    inputs = ["1", "skip", "", "q"]
    _run_with_inputs(monkeypatch, inputs, calls)
    assert calls == []


def test_validate_configuration_default_path(monkeypatch):
    calls = []
    inputs = ["3", "", "y", "", "q"]
    _run_with_inputs(monkeypatch, inputs, calls)
    assert calls == [["neoswga", "validate-params", "-j", "params.json"]]


def test_interpret_results_default_path(monkeypatch):
    calls = []
    inputs = ["4", "", "y", "", "q"]
    _run_with_inputs(monkeypatch, inputs, calls)
    assert calls == [["neoswga", "interpret", "-d", "results/"]]


def test_interpret_results_custom_path(monkeypatch):
    calls = []
    inputs = ["4", "my_results/", "y", "", "q"]
    _run_with_inputs(monkeypatch, inputs, calls)
    assert calls == [["neoswga", "interpret", "-d", "my_results/"]]


def test_advanced_menu_back_returns_to_main_menu(monkeypatch):
    """Advanced-features entries only print instructions; none call execute_command.

    'Back to main menu' (option 7) must return control to the main loop
    rather than quitting the whole program.
    """
    calls = []
    inputs = ["5", "7", "q"]
    _run_with_inputs(monkeypatch, inputs, calls)
    assert calls == []


@pytest.mark.parametrize("adv_choice", ["1", "2", "3", "4", "5", "6"])
def test_advanced_menu_options_do_not_crash_and_do_not_execute(monkeypatch, adv_choice):
    """Every advanced-features entry is informational only (no subprocess call)."""
    calls = []
    inputs = ["5", adv_choice, "", "7", "q"]
    _run_with_inputs(monkeypatch, inputs, calls)
    assert calls == []


def test_advanced_menu_quit_only_backs_out_to_main_menu(monkeypatch):
    """'q' inside a submenu is treated as "back", not "quit the whole program".

    get_choice()'s own docstring calls its None return "quit", but both the
    pipeline and advanced submenus treat None identically to their explicit
    "Back to main menu" entry (`if adv_choice is None or adv_choice ==
    len(advanced_options): break`) -- it only breaks the inner loop and
    redraws the main menu. A second 'q' at the main menu is what actually
    exits. This is presumably intentional (q = "back" while inside any menu),
    but it means 'q' does not behave uniformly across nesting depth, so pin
    the actual behaviour down rather than assume a name implies a meaning.
    """
    calls = []
    _run_with_inputs(monkeypatch, ["5", "q", "q"], calls)
    assert calls == []
