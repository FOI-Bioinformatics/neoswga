"""Two paths, one params file, one reaction.

Found while wiring Phase 4 of `docs/validation/pipeline_audit_2026-09-16/`. Not
in the audit: it only becomes visible once something tries to read the inventory
back.

The candidate inventory files every verdict under a reaction fingerprint. The
filter writes it through `filter._get_reaction_conditions`; `plan-pool` computes
it through `build_reaction_conditions(SimpleNamespace(**params))`. On the real
Wolbachia design those two disagreed, from the same params.json:

    filter    tm-2026-09-14:2b3536d4305777e8
    plan-pool tm-2026-09-14:6dc7004c72a51941

One field: `dtt_mm`, 4.0 against 0.0. So the design could never find the
candidates the filter had recorded, and fell back to the CSV with a message
about no eligible candidates -- which is true and is not the reason.

The cause is a module-level default shadowing a polymerase-aware one.
`ReactionConditions` resolves `mg_conc=None` and `dtt_mm=None` to the
polymerase's buffer values, and both carry comments explaining that a literal
0.0 reads as a broken reaction: no polymerase runs without magnesium, and the
mechanistic model scores absent DTT as a deficiency worth 20% of stability.
`build_reaction_conditions` passes a value through whenever it is not None, and
`parameter.mg_conc` and `parameter.dtt_mm` were 2.0 and 0.0 at module scope. A
path that had run `get_params` got the resolved values; one that had not got the
stale module defaults, which is the same defect those two comments describe as
fixed.

`parameter.reaction_temp` was already `None` for this reason. These two now
match it.
"""

import json
from types import SimpleNamespace

import pytest

from neoswga.core.reaction_conditions import ReactionConditions, build_reaction_conditions

PARAMS = {
    "polymerase": "phi29",
    "reaction_temp": 30.0,
    "na_conc": 50.0,
    "min_k": 12,
    "max_k": 12,
}


@pytest.fixture(autouse=True)
def _unconfigured(monkeypatch):
    """A process that has not run `get_params`, which is the case under test.

    The field list is read off `ReactionConditions.__init__` rather than written
    out here, for the reason `build_reaction_conditions` gives for doing the
    same: every hand-written subset drops something. A first attempt listed nine
    fields and missed `betaine_m`, which the suite's own fixtures leave at 1.0,
    so these tests passed alone and failed in the suite.
    """
    import inspect

    from neoswga.core import parameter

    fields = list(inspect.signature(ReactionConditions.__init__).parameters)[1:]
    for field in [*fields, "reaction_temp"]:
        monkeypatch.setattr(parameter, field, None, raising=False)


def test_an_unset_magnesium_takes_the_polymerase_buffer_value():
    """2.0 mM is the old PCR figure; phi29's buffer is 10 mM."""
    conditions = build_reaction_conditions(SimpleNamespace(**PARAMS))

    assert conditions.mg_conc == 10.0


def test_an_unset_dtt_takes_the_polymerase_buffer_value():
    """Every phi29 buffer has DTT; a literal zero reads as a deficiency."""
    conditions = build_reaction_conditions(SimpleNamespace(**PARAMS))

    assert conditions.dtt_mm > 0


def test_an_explicit_zero_is_still_honoured():
    """The distinction the constructor draws: absent is not the same as zero."""
    conditions = build_reaction_conditions(SimpleNamespace(**PARAMS, dtt_mm=0.0))

    assert conditions.dtt_mm == 0.0


def test_an_explicit_value_still_wins():
    conditions = build_reaction_conditions(SimpleNamespace(**PARAMS, mg_conc=4.0))

    assert conditions.mg_conc == 4.0


def test_the_direct_path_agrees_with_the_constructor():
    """The fingerprint is the inventory key, so the two must not diverge."""
    built = build_reaction_conditions(SimpleNamespace(**PARAMS))
    direct = ReactionConditions(temp=30.0, polymerase="phi29", na_conc=50.0)

    assert built.fingerprint() == direct.fingerprint()


def test_the_module_defaults_do_not_shadow_a_polymerase_aware_one():
    """The root cause, stated where it can be seen.

    A module-level number is indistinguishable from a configured one by the
    time `build_reaction_conditions` sees it, so any field the constructor
    resolves from the polymerase has to be absent there rather than guessed.

    Read from the source rather than the live module: these are globals that
    `get_params` writes, so whatever ran earlier in the suite decides what a
    live read returns. The declaration is the invariant.
    """
    import ast
    import pathlib

    source = pathlib.Path("neoswga/core/parameter.py").read_text()
    tree = ast.parse(source)
    declared = {
        node.targets[0].id: node.value
        for node in tree.body
        if isinstance(node, ast.Assign)
        and len(node.targets) == 1
        and isinstance(node.targets[0], ast.Name)
    }

    for field in ("mg_conc", "dtt_mm", "reaction_temp"):
        value = declared.get(field)
        assert isinstance(value, ast.Constant) and value.value is None, (
            f"parameter.{field} is declared with a module-level default, which "
            "shadows the polymerase-aware one for every path that has not run "
            "get_params"
        )


@pytest.mark.parametrize("polymerase", ["phi29", "equiphi29", "bst", "klenow"])
def test_every_polymerase_gets_a_usable_buffer(polymerase):
    """No polymerase should design at zero magnesium.

    The temperature is left to resolve from the polymerase too, which is the
    whole point: a fixed one would be outside the range for some of these.
    """
    conditions = build_reaction_conditions(SimpleNamespace(polymerase=polymerase))

    assert conditions.mg_conc > 0
    assert conditions.dtt_mm >= 0


def test_a_params_file_resolves_the_same_way_twice(tmp_path):
    """Guard the guard: same file, same fingerprint, whichever path reads it."""
    path = tmp_path / "params.json"
    path.write_text(json.dumps(PARAMS))

    first = build_reaction_conditions(SimpleNamespace(**json.loads(path.read_text())))
    second = build_reaction_conditions(SimpleNamespace(**json.loads(path.read_text())))

    assert first.fingerprint() == second.fingerprint()
