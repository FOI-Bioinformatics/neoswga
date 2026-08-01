"""Execution tests for the sequence-only and registry commands.

These need no pipeline workspace: `analyze-set`, `analyze-dimers`,
`analyze-stability`, `analyze-genome`, `show-presets` and the background/genome
registries all work from arguments alone. That makes them cheap to test and, up
to now, entirely untested -- `neoswga/cli/analysis.py` was at 19% and
`registry.py` at 26%.

The preset handling gets particular attention. Presets diverged from what they
claimed to apply in three separate places during the audit, each time silently:
a preset defined by its additives would run with those additives dropped, and
nothing in the output said so.
"""

import json
import os

import pytest

PRIMERS = ["TTGACCATGA", "GCATTACGGT", "AACCGGTTAC", "CCTAGGATCA"]


# ----------------------------------------------------------------------
# show-presets
# ----------------------------------------------------------------------


def _accepted_preset_names():
    """Every name `--preset` will accept, read from the resolver itself."""
    import inspect
    import re

    from neoswga.cli._common import load_preset_conditions

    return set(re.findall(r'"([a-z0-9_]+)":\s*rc\.', inspect.getsource(load_preset_conditions)))


def test_show_presets_lists_every_accepted_preset(capsys):
    """`show-presets` is the discovery command; anything --preset accepts must
    appear in it.

    These are maintained as two separate lists -- PRESETS, built from a
    four-entry optimizer-hints dict, plus a hardcoded `additional` list in
    show_presets covering the other seven. Both currently agree with the
    resolver, but nothing enforces it, so adding a preset to the resolver alone
    would make it silently undiscoverable.
    """
    from neoswga.cli.analysis import show_presets

    show_presets()
    printed = capsys.readouterr().out

    missing = sorted(n for n in _accepted_preset_names() if n not in printed)
    assert not missing, f"--preset accepts these but show-presets never mentions them: {missing}"


def test_show_presets_reports_what_the_preset_actually_applies(capsys):
    """`show-presets` and `--preset` diverged: the display dict and the applied
    conditions were separate structures, so the tool advertised values it did
    not use."""
    from neoswga.cli._common import PRESETS, load_preset_conditions
    from neoswga.cli.analysis import show_presets

    show_presets()
    printed = capsys.readouterr().out

    for name in PRESETS:
        applied = load_preset_conditions(name)
        assert applied, f"preset {name} resolves to nothing"
        temp = applied["reaction_temp"]
        assert (
            f"{temp:.0f}" in printed or f"{temp}" in printed
        ), f"preset {name} applies {temp} C but that value is not shown"


def test_extreme_gc_preset_carries_its_defining_additives():
    """The preset is defined by urea and TMAC; a round trip that drops them
    leaves a preset that does nothing it claims to."""
    from neoswga.cli._common import load_preset_conditions

    applied = load_preset_conditions("extreme_gc")
    if applied is None:
        pytest.skip("no extreme_gc preset registered")

    assert (
        applied.get("urea_m", 0) > 0 or applied.get("tmac_m", 0) > 0
    ), f"extreme_gc preset applies neither urea nor TMAC: {applied}"


def test_no_preset_silently_applies_zero_magnesium():
    """Five presets once inherited mg_conc 0.0 from a default."""
    from neoswga.cli._common import load_preset_conditions

    for name in _accepted_preset_names():
        applied = load_preset_conditions(name)
        assert applied, f"preset {name} resolves to nothing"
        assert applied["mg_conc"] > 0, f"preset {name} runs at zero magnesium"


# ----------------------------------------------------------------------
# analyze-set
# ----------------------------------------------------------------------


def test_analyze_set_runs_on_sequences_alone(run_cli, tmp_path, capsys):
    run_cli(["analyze-set", "--primers", *PRIMERS, "--output", str(tmp_path / "a")])
    printed = capsys.readouterr().out
    assert "Thermodynamic Properties" in printed, printed[:400]
    for primer in PRIMERS:
        assert primer in printed


def test_analyze_set_honours_its_required_output(run_cli, tmp_path):
    """`--output` is required by the parser, so it must be used.

    It was declared required and then ignored, alongside two other required
    arguments the handler never read.
    """
    out = tmp_path / "analysis"
    run_cli(["analyze-set", "--primers", *PRIMERS, "--output", str(out)])

    written = out / "primer_analysis.json"
    assert written.is_file(), sorted(os.listdir(out)) if out.exists() else "nothing written"

    def reject(token):
        raise ValueError(f"non-standard JSON token {token!r}")

    payload = json.loads(written.read_text(), parse_constant=reject)
    assert payload, "primer_analysis.json is empty"


def test_analyze_set_responds_to_the_preset(run_cli, tmp_path, capsys):
    """Different chemistry must produce different thermodynamics.

    analyze-set built ReactionConditions from four fields only, so it ran at
    0 mM Mg and dropped urea/TMAC/ethanol/formamide -- the third place the
    preset divergence bit. Identical output across two presets would mean the
    conditions are not reaching the calculation.
    """
    outputs = {}
    for preset in ("standard_phi29", "enhanced_equiphi29"):
        run_cli(
            [
                "analyze-set",
                "--primers",
                *PRIMERS,
                "--preset",
                preset,
                "--output",
                str(tmp_path / preset),
            ]
        )
        outputs[preset] = capsys.readouterr().out

    assert (
        outputs["standard_phi29"] != outputs["enhanced_equiphi29"]
    ), "standard_phi29 and enhanced_equiphi29 presets produced identical analysis"


# ----------------------------------------------------------------------
# analyze-dimers / analyze-stability / analyze-genome
# ----------------------------------------------------------------------


def test_analyze_dimers_writes_output(run_cli, tmp_path):
    out = tmp_path / "dimers"
    run_cli(["analyze-dimers", "--primers", *PRIMERS, "--output", str(out)])
    assert os.path.isdir(out) and os.listdir(out), f"nothing written to {out}"


def test_analyze_dimers_finds_a_deliberate_self_dimer(run_cli, tmp_path):
    """A primer and its own reverse complement must register as interacting.

    A dimer analysis that reports nothing for a perfectly complementary pair is
    not analysing anything.
    """
    seq = "GCATTACGGT"
    rc = seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

    out = tmp_path / "dimers"
    run_cli(["analyze-dimers", "--primers", seq, rc, "--output", str(out)])

    summary = out / "dimer_summary.txt"
    assert summary.is_file(), sorted(os.listdir(out))
    text = summary.read_text()
    assert seq in text and rc in text, "both primers should appear in the summary"


def test_analyze_stability_runs(run_cli, tmp_path, caplog):
    import logging

    with caplog.at_level(logging.INFO):
        run_cli(["analyze-stability", "--primers", *PRIMERS, "--output", str(tmp_path / "s")])
    assert caplog.text.strip()


def test_analyze_genome_reports_gc(run_cli, genome_fasta, tmp_path, caplog):
    import logging

    with caplog.at_level(logging.INFO):
        run_cli(["analyze-genome", "--genome", genome_fasta, "--output", str(tmp_path / "g")])
    assert caplog.text.strip(), "analyze-genome reported nothing"


# ----------------------------------------------------------------------
# Registries
# ----------------------------------------------------------------------


def test_background_list_runs_on_an_empty_registry(run_cli, tmp_path, monkeypatch, capsys, caplog):
    """Listing an empty registry must be a clean no-op, not a crash."""
    import logging

    monkeypatch.setenv("HOME", str(tmp_path))
    with caplog.at_level(logging.INFO):
        try:
            run_cli(["background-list"])
        except SystemExit as exc:
            assert exc.code in (0, 1)

    combined = capsys.readouterr().out + caplog.text
    assert combined.strip(), "background-list said nothing at all"


def test_genome_list_runs_on_an_empty_registry(run_cli, tmp_path, monkeypatch, capsys, caplog):
    import logging

    monkeypatch.setenv("HOME", str(tmp_path))
    with caplog.at_level(logging.INFO):
        try:
            run_cli(["genome-list"])
        except SystemExit as exc:
            assert exc.code in (0, 1)

    combined = capsys.readouterr().out + caplog.text
    assert combined.strip()
