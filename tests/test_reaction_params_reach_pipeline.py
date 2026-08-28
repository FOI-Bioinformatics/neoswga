"""Guard: a reaction-chemistry parameter that is accepted must actually be used.

Four parameters - `glycerol_percent`, `peg_percent`, `bsa_ug_ml` and `ssb` -
were valid in the JSON schema, exposed as CLI flags, and merged onto the
parameter module, but had no field in `PipelineParameters` and were never read
by `get_params`. Setting them changed nothing, `validate-params` passed clean,
and nothing said so. `ssb` had no schema entry at all.

`parameter.REACTION_PARAM_DEFAULTS` is now the single list these flow from.
These tests keep it in step with every surface that has to agree with it, so the
next added additive cannot be half-wired.
"""

import inspect
import json
from pathlib import Path

import pytest

from neoswga.core import parameter as P
from neoswga.core.reaction_conditions import ReactionConditions

SCHEMA_PATH = (
    Path(__file__).resolve().parents[1] / "neoswga" / "core" / "schema" / "params.schema.json"
)

# Parameters the pipeline records but which are deliberately absent from the Tm
# model - they act on the enzyme or are protocol metadata, not on the duplex.
NO_TM_MODEL = {"glycerol_percent", "peg_percent", "bsa_ug_ml", "dtt_mm", "ssb"}


def test_every_reaction_param_has_a_pipeline_field():
    from dataclasses import fields

    declared = {f.name for f in fields(P.PipelineParameters)}
    missing = set(P.REACTION_PARAM_DEFAULTS) - declared
    assert not missing, f"PipelineParameters is missing {sorted(missing)}"


def test_every_reaction_param_has_a_module_global():
    missing = [n for n in P.REACTION_PARAM_DEFAULTS if not hasattr(P, n)]
    assert not missing, f"parameter module globals missing {sorted(missing)}"


def test_every_reaction_param_is_accepted_by_reaction_conditions():
    """Anything the pipeline carries must be settable on the conditions object."""
    accepted = set(inspect.signature(ReactionConditions.__init__).parameters)
    missing = set(P.REACTION_PARAM_DEFAULTS) - accepted
    assert not missing, (
        f"ReactionConditions cannot accept {sorted(missing)} - the pipeline would "
        f"carry values it can never apply"
    )


def test_every_reaction_param_is_in_the_schema():
    schema = json.loads(SCHEMA_PATH.read_text())
    missing = set(P.REACTION_PARAM_DEFAULTS) - set(schema["properties"])
    assert not missing, f"params.schema.json is missing {sorted(missing)}"


def test_every_numeric_reaction_param_has_a_validator_range():
    from neoswga.core.param_validator import PARAM_RANGES

    numeric = {
        n
        for n, d in P.REACTION_PARAM_DEFAULTS.items()
        if isinstance(d, (int, float)) and not isinstance(d, bool)
    }
    missing = numeric - set(PARAM_RANGES)
    assert not missing, f"param_validator has no range for {sorted(missing)}"


@pytest.mark.parametrize("name", sorted(P.REACTION_PARAM_DEFAULTS))
def test_reaction_param_survives_get_params(name, tmp_path, monkeypatch):
    """A value in params.json must reach the module global, not vanish."""
    genome = tmp_path / "t.fna"
    genome.write_text(">c\n" + "ATGC" * 300 + "\n")

    default = P.REACTION_PARAM_DEFAULTS[name]
    probe = True if isinstance(default, bool) else 1.25

    params = {
        "data_dir": str(tmp_path),
        "src_dir": str(tmp_path),
        "fg_genomes": [str(genome)],
        "fg_prefixes": [str(tmp_path / "t")],
        "cpus": 1,
        name: probe,
    }
    cfg = tmp_path / "params.json"
    cfg.write_text(json.dumps(params))

    class Args:
        def __getattr__(self, _):
            return None

    args = Args()
    args.json_file = str(cfg)

    P.get_params(args)
    assert getattr(P, name) == probe, (
        f"{name} was set to {probe} in params.json but the pipeline global is "
        f"{getattr(P, name)!r} - the value was silently dropped"
    )


def test_params_with_no_tm_model_are_documented_as_such():
    """The no-Tm-model set must be exactly what the Tm layer omits.

    These are real reagents with real bench effects; the point is that the tool
    does not model their contribution to melting temperature, and a user should
    not read 'no Tm shift' as 'no effect'.
    """
    from neoswga.core.mechanistic_params import ADDITIVE_TM_PARAMS

    modelled = set(ADDITIVE_TM_PARAMS)
    for name in NO_TM_MODEL:
        stem = name.replace("_percent", "").replace("_ug_ml", "").replace("_mm", "")
        assert stem not in modelled, (
            f"{name} now has a Tm model - remove it from NO_TM_MODEL and give it "
            f"a real coefficient rather than leaving it silently at 0.0"
        )


def test_real_phi29_buffer_is_expressible():
    """The motivating case: 10 mM MgCl2 / 10 mM (NH4)2SO4 / 4 mM DTT, no sodium.

    Before k_conc/nh4_conc existed the only way to represent this was to
    misreport the ammonium and potassium as na_conc.
    """
    conditions = ReactionConditions(
        temp=30.0,
        polymerase="phi29",
        na_conc=0.0,
        nh4_conc=20.0,
        k_conc=66.0,
        mg_conc=10.0,
        dntp_conc=1.6,
        dtt_mm=4.0,
    )
    assert conditions.na_conc == 0.0
    assert conditions.nh4_conc == 20.0
    assert conditions.dtt_mm == 4.0

    tm = conditions.calculate_effective_tm("ATCGATCGATCG")
    assert tm > 0, "a buffer with no sodium must still give a usable Tm"


def test_monovalent_cations_are_interchangeable_by_ionic_strength():
    from neoswga.core import thermodynamics as td

    seq = "ATCGATCGATCG"
    via_na = td.calculate_tm_with_salt(seq, na_conc=50, mg_conc=10)
    via_k = td.calculate_tm_with_salt(seq, na_conc=0, mg_conc=10, k_conc=50)
    assert via_na == pytest.approx(via_k)


def test_dntp_chelation_lowers_free_magnesium():
    from neoswga.core import thermodynamics as td

    assert td.free_magnesium_mm(10.0, 1.6) == pytest.approx(8.4)
    assert td.free_magnesium_mm(10.0, 0.0) == pytest.approx(10.0)
    # Never negative, even if someone over-specifies dNTPs.
    assert td.free_magnesium_mm(1.0, 5.0) == 0.0


def test_new_buffer_params_default_to_previous_behaviour():
    """Zero defaults must reproduce the pre-change numbers exactly."""
    from neoswga.core import thermodynamics as td

    seq = "GCTAGCTAGCTA"
    assert td.calculate_tm_with_salt(seq, 50, 10) == pytest.approx(
        td.calculate_tm_with_salt(seq, 50, 10, k_conc=0.0, nh4_conc=0.0, dntp_conc=0.0)
    )
