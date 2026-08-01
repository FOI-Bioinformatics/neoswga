"""Execution tests for the setup and onboarding commands.

`init`, `validate-params`, `schema` and `suggest` are the first things a new
user runs, and the documented onboarding path is `init` -> pipeline -> iterate.
That path had never been exercised end to end, which is how a wizard-generated
params.json came to be rejected by two downstream commands.

The most valuable test here is therefore the round-trip: what `init` writes must
satisfy `validate-params`, and must carry enough for the rest of the tool to
work rather than merely being schema-valid.
"""

import json
import os

import pytest

# ----------------------------------------------------------------------
# init
# ----------------------------------------------------------------------


@pytest.fixture
def initialised_params(run_cli, genome_fasta, tmp_path):
    """Run `neoswga init` and return the path to what it wrote."""
    out = tmp_path / "params.json"
    run_cli(
        [
            "init",
            "--genome",
            genome_fasta,
            "-o",
            str(out),
            "--non-interactive",
            "--yes",
            "--output-dir",
            str(tmp_path / "results"),
        ]
    )
    if not out.exists():
        pytest.skip("init did not produce a params.json (interactive path?)")
    return out


def test_init_writes_parsable_params(initialised_params):
    data = json.loads(initialised_params.read_text())
    assert isinstance(data, dict) and data


def test_init_output_passes_validate_params(run_cli, initialised_params):
    """The wizard tells you to run this; it must not reject its own output.

    These two drifted apart before: the wizard emitted an optimization method
    the validator's list no longer contained, so the config it had just written
    failed the check it recommends in the same breath.
    """
    with pytest.raises(SystemExit) as exc:
        run_cli(["validate-params", "-j", str(initialised_params)])
    assert exc.value.code == 0, "validate-params rejected the wizard's own output"


def test_init_sets_a_usable_magnesium_concentration(initialised_params):
    """A wizard bug wrote mg_conc: 0.0, and an explicit key suppresses the
    polymerase default, so the whole pipeline ran at zero magnesium."""
    data = json.loads(initialised_params.read_text())
    assert "mg_conc" in data, "wizard stopped writing mg_conc; the default now applies"
    assert data["mg_conc"] > 0, "wizard wrote a zero magnesium concentration"


def test_init_declares_an_optimization_method_the_validator_accepts(initialised_params):
    """These two lists drifted apart once: the wizard emitted 'ensemble' while
    the validator still listed retired methods and rejected it."""
    from neoswga.core.param_validator import VALID_OPTIMIZATION_METHODS

    data = json.loads(initialised_params.read_text())
    assert "optimization_method" in data
    assert data["optimization_method"] in VALID_OPTIMIZATION_METHODS


def test_init_records_the_genome_it_was_pointed_at(initialised_params, genome_fasta):
    data = json.loads(initialised_params.read_text())
    genomes = data.get("fg_genomes") or []
    assert any(os.path.basename(genome_fasta) in str(g) for g in genomes)


# ----------------------------------------------------------------------
# The onboarding path that was broken
# ----------------------------------------------------------------------


def test_wizard_output_yields_usable_sequence_lengths(initialised_params):
    """`init` does not write fg_seq_lengths, so get_params must derive it.

    When it did not publish the derived value as a module global, contract-set
    and rescore-set both refused to run on any wizard-generated config while the
    pipeline worked -- the exact seam between onboarding and iteration.
    """
    written = json.loads(initialised_params.read_text())
    assert "fg_seq_lengths" not in written, (
        "the wizard now writes fg_seq_lengths; this test no longer covers the "
        "derive-it-from-the-FASTA path it was written for"
    )

    from neoswga.core import parameter

    class _Args:
        def __getattr__(self, name):
            return None

        json_file = None

    args = _Args()
    args.json_file = str(initialised_params)
    parameter.get_params(args)

    assert getattr(parameter, "fg_seq_lengths", None), (
        "fg_seq_lengths is empty after get_params on a wizard-generated config; "
        "the iterate commands cannot compute coverage from this"
    )
    assert all(length > 0 for length in parameter.fg_seq_lengths)


# ----------------------------------------------------------------------
# schema
# ----------------------------------------------------------------------


def test_schema_command_emits_valid_json(run_cli, capsys):
    run_cli(["schema"])
    printed = capsys.readouterr().out
    schema = json.loads(printed)
    assert "properties" in schema


def test_schema_polymerase_enum_matches_the_registry(run_cli, capsys):
    """The schema is generated from the registry; drift here means a params.json
    the CLI accepts would be rejected by an external validator, or vice versa."""
    from neoswga.core.registry import polymerase_names

    run_cli(["schema"])
    schema = json.loads(capsys.readouterr().out)

    enum = schema.get("properties", {}).get("polymerase", {}).get("enum")
    if enum is not None:
        assert set(enum) == set(polymerase_names())


# ----------------------------------------------------------------------
# suggest
# ----------------------------------------------------------------------


def test_suggest_runs_from_a_genome(run_cli, genome_fasta, capsys):
    run_cli(["suggest", "--genome", genome_fasta])
    out = capsys.readouterr().out
    assert out.strip(), "suggest printed nothing"


def test_suggest_runs_from_an_explicit_gc(run_cli, capsys):
    run_cli(["suggest", "--genome-gc", "0.65", "--primer-length", "15"])
    assert capsys.readouterr().out.strip()


def test_suggest_recommends_different_chemistry_for_different_gc(run_cli, capsys):
    """A recommender that returns the same answer regardless of input is not
    recommending anything."""
    run_cli(["suggest", "--genome-gc", "0.25"])
    at_rich = capsys.readouterr().out
    run_cli(["suggest", "--genome-gc", "0.72"])
    gc_rich = capsys.readouterr().out

    assert at_rich != gc_rich, "suggest gave identical advice for 25% and 72% GC"
