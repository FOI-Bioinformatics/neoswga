"""Every reaction-condition field must reach the chemistry that uses it.

`ReactionConditions` takes 20 parameters -- temperature, polymerase, four
buffer species, and thirteen additives. Nothing constructed it with all of
them. Each of the ~19 call sites hand-listed the subset its author happened to
need, so the list of additives a command respects was decided independently in
nineteen places and had drifted accordingly:

    cli/evaluate.py                 2 fields (temp, polymerase)
    core/thermodynamic_filter.py    4
    core/pipeline.py (score step)   5
    cli/iterate.py (x3)             7
    core/filter.py                  11
    core/unified_optimizer.py       11

Not one passed `k_conc`, `nh4_conc`, `dntp_conc` or `dtt_mm`, so the phi29
buffer work landed earlier reached nothing. `glycerol_percent`, `peg_percent`,
`bsa_ug_ml` and `ssb` were parsed by the CLI and dropped at almost every site.

The failure is silent in the worst way: a dropped additive does not raise, it
just leaves the Tm at its no-additive value. A user who sets 5% DMSO and sees
the melting temperature not move has no way to tell whether the model thinks
DMSO does not matter or never heard about it.

The fix is to stop hand-listing. `build_reaction_conditions` derives the field
list from the constructor signature, so a field added to `ReactionConditions`
is forwarded by every caller without anyone remembering to update them. These
tests pin that property -- the parametrised test below enumerates the signature
rather than a copy of it, so it grows automatically too.
"""

import inspect
from pathlib import Path

import pytest

from neoswga.core import parameter
from neoswga.core.reaction_conditions import ReactionConditions, build_reaction_conditions

# Derived from the constructor, never hand-listed: that is the point.
CONDITION_FIELDS = [p for p in inspect.signature(ReactionConditions.__init__).parameters][1:]

# `temp` is the one field with no same-named global (it reads `reaction_temp`),
# and `polymerase` is a string rather than a concentration.
ADDITIVE_FIELDS = [f for f in CONDITION_FIELDS if f not in ("temp", "polymerase")]

# Distinctive, in-range values. Each has to differ from the field's default or
# the assertion could pass without anything being forwarded.
PROBE_VALUES = {
    "dmso_percent": 4.0,
    "betaine_m": 1.2,
    "trehalose_m": 0.4,
    "formamide_percent": 3.0,
    "glycerol_percent": 6.0,
    "bsa_ug_ml": 120.0,
    "peg_percent": 5.0,
    "ethanol_percent": 2.0,
    "urea_m": 0.7,
    "tmac_m": 0.05,
    "propanediol_m": 0.8,
    "na_conc": 75.0,
    "mg_conc": 7.5,
    "k_conc": 45.0,
    "nh4_conc": 12.0,
    "dntp_conc": 0.9,
    "dtt_mm": 3.0,
    "ssb": True,
}


@pytest.fixture
def clean_parameter(monkeypatch):
    """Set every additive to its zero value so one probe at a time is visible."""
    for field in ADDITIVE_FIELDS:
        monkeypatch.setattr(parameter, field, False if field == "ssb" else 0.0, raising=False)
    monkeypatch.setattr(parameter, "polymerase", "phi29", raising=False)
    monkeypatch.setattr(parameter, "reaction_temp", 30.0, raising=False)
    return parameter


# ----------------------------------------------------------------------
# The property that makes this stay fixed
# ----------------------------------------------------------------------


def test_every_constructor_field_has_a_probe_value():
    """Guards the test itself.

    If `ReactionConditions` gains a field and nobody adds a probe for it, the
    parametrised test below would silently skip it -- the same
    forgot-to-update-the-list failure this whole file exists to prevent.
    """
    missing = [f for f in ADDITIVE_FIELDS if f not in PROBE_VALUES]
    assert not missing, f"no probe value for new ReactionConditions field(s): {missing}"


@pytest.mark.parametrize("field", ADDITIVE_FIELDS)
def test_builder_forwards_every_field(clean_parameter, monkeypatch, field):
    """Set one field on `parameter`; it must arrive on the conditions object."""
    probe = PROBE_VALUES[field]
    monkeypatch.setattr(parameter, field, probe, raising=False)

    conditions = build_reaction_conditions()

    assert getattr(conditions, field) == probe, (
        f"{field} was set to {probe!r} on parameter and did not reach "
        f"ReactionConditions (got {getattr(conditions, field)!r})"
    )


def test_temperature_comes_from_reaction_temp(clean_parameter, monkeypatch):
    """The one field whose global has a different name."""
    monkeypatch.setattr(parameter, "reaction_temp", 37.5, raising=False)
    assert build_reaction_conditions().temp == 37.5


def test_missing_temperature_falls_back_to_the_polymerase_optimum(monkeypatch):
    """`ReactionConditions(temp=None)` raises when it range-checks, so a params
    file without `reaction_temp` must not reach the constructor with None."""
    monkeypatch.setattr(parameter, "reaction_temp", None, raising=False)
    monkeypatch.setattr(parameter, "polymerase", "equiphi29", raising=False)

    conditions = build_reaction_conditions()

    assert conditions.temp is not None
    assert 35.0 <= conditions.temp <= 50.0, "not the equiphi29 operating range"


def test_unset_magnesium_keeps_the_polymerase_default(monkeypatch):
    """`mg_conc=None` means "use the polymerase's buffer value" -- 10 mM for
    phi29. Coercing an absent value to 0.0, as several call sites did, runs the
    whole reaction at zero magnesium and the model then penalises it.
    """
    monkeypatch.setattr(parameter, "mg_conc", None, raising=False)
    monkeypatch.setattr(parameter, "polymerase", "phi29", raising=False)
    monkeypatch.setattr(parameter, "reaction_temp", 30.0, raising=False)

    assert build_reaction_conditions().mg_conc == 10.0


def test_explicit_arguments_win_over_parameter(clean_parameter, monkeypatch):
    """CLI flags have to be able to override the params.json values."""
    monkeypatch.setattr(parameter, "polymerase", "phi29", raising=False)

    conditions = build_reaction_conditions(polymerase="equiphi29", temp=42.0)

    assert conditions.polymerase == "equiphi29"
    assert conditions.temp == 42.0


# ----------------------------------------------------------------------
# The call sites
# ----------------------------------------------------------------------

ROOT = Path(__file__).resolve().parent.parent / "neoswga"

# Modules that build conditions from pipeline configuration. These are the ones
# that were dropping fields; they must go through the builder rather than
# hand-listing a subset again.
PRODUCTION_SITES = [
    "cli/evaluate.py",
    "cli/iterate.py",
    "core/filter.py",
    "core/pipeline.py",
    "core/thermodynamic_filter.py",
    "core/unified_optimizer.py",
]


@pytest.mark.parametrize("relpath", PRODUCTION_SITES)
def test_production_sites_do_not_hand_list_additive_fields(relpath):
    """A hand-written list of additive keywords is how the drift started.

    What is flagged is specifically naming additive or buffer fields at the
    call site, since that is the list that goes stale. Two constructions are
    fine and are not flagged: `ReactionConditions(**kwargs)` rebuilt from a
    complete field mapping (the multiprocessing worker in
    `thermodynamic_filter`), and mentions in comments or type annotations.

    Excluded from the list entirely: `reaction_conditions.py` itself, which
    defines the presets and must name fields.
    """
    import re

    source = (ROOT / relpath).read_text()

    offenders = []
    for match in re.finditer(r"ReactionConditions\((.*?)\)", source, re.DOTALL):
        args = match.group(1)
        named = [f for f in ADDITIVE_FIELDS if re.search(rf"\b{f}\s*=", args)]
        if named:
            offenders.append(named)

    assert not offenders, (
        f"{relpath} hand-lists additive fields {offenders}; use "
        f"build_reaction_conditions() so new fields cannot be dropped by omission"
    )


def test_builder_is_importable_without_the_optimizer_stack():
    """It is called from `filter` and from CLI paths, so it must not drag in
    numpy/scipy-heavy modules just to read a few floats."""
    import subprocess
    import sys

    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "from neoswga.core.reaction_conditions import build_reaction_conditions; "
            "print(build_reaction_conditions(temp=30.0).dmso_percent)",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
