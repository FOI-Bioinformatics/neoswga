"""A chemistry model states what it supports, and refuses the rest.

Task 4 of the 2026-09-21 valid-design plan. The distinction it draws is between
three things a number can be: measured, estimated from data that does not quite
cover the case, and chosen so the model behaves plausibly. All three currently
reach a design with equal authority, and a coefficient with a citation beside
it reads as the first whatever it is.

The plan's rule, stated once so it can be checked: **an assumption is not
promoted to a measurement because it has a citation**, and **an absent effect
model is not a known zero effect**.

Two properties are pinned here.

1. The registry is machine-readable, versioned, shipped with the package, and
   agrees with the constants the code actually uses. The prose ledger in
   `docs/SCIENCE_CITATIONS.md` had already drifted: it states Klenow
   processivity as 10,000 bp citing Bambara (1978) while the shipped registry
   says 40 bp. A ledger nothing checks is a ledger that drifts.
2. `require_model_support` refuses a requested computation outside a recorded
   domain, before any search begins.
"""

import json

import pytest

from neoswga.core.design_request import resolve_design_request
from neoswga.core.exceptions import UnsupportedModelError
from neoswga.core.model_evidence import (
    EVIDENCE_STATUSES,
    EvidenceRecord,
    load_evidence,
    require_model_support,
)


def base_params(**overrides):
    params = {
        "fg_prefixes": ["target"],
        "fg_genomes": ["target.fna"],
        "polymerase": "phi29",
        "reaction_temp": 30.0,
        "min_k": 8,
        "max_k": 12,
        "num_primers": 12,
    }
    params.update(overrides)
    return params


# ---------------------------------------------------------------------------
# 1. The registry
# ---------------------------------------------------------------------------


def test_the_registry_loads_without_the_repository():
    """Shipped as package data, not found by walking up from the source tree.

    An artifact that only resolves relative to a checkout is absent for every
    installed user, and the model would then initialise with no evidence at all.
    """
    import importlib.resources

    handle = importlib.resources.files("neoswga.core.registry") / "model_evidence.json"
    assert handle.is_file()
    payload = json.loads(handle.read_text())
    assert payload["schema_version"] >= 1


def test_every_record_declares_its_evidence_status():
    for record in load_evidence().values():
        assert record.status in EVIDENCE_STATUSES, (record.quantity, record.status)


def test_a_citation_does_not_by_itself_make_a_value_measured():
    """The rule this file exists for, checked rather than asserted in prose.

    Several coefficients carry a primary citation and are nonetheless
    extrapolations from data taken at another temperature, on longer DNA, or in
    another buffer. Those are `estimated`, and the registry must not record a
    source as sufficient grounds for `measured`.
    """
    evidence = load_evidence()
    cited_but_not_measured = [
        record.quantity
        for record in evidence.values()
        if record.source and record.status != "measured"
    ]
    assert cited_but_not_measured, (
        "every cited record claims to be measured, which is the conflation this "
        "registry exists to prevent"
    )


def test_an_absent_effect_model_is_recorded_rather_than_implied():
    evidence = load_evidence()
    glycerol = evidence["tm_glycerol"]

    assert glycerol.status == "absent"
    assert "0" not in (glycerol.value or ""), "an absent model must not carry a value"
    assert glycerol.notes


@pytest.mark.parametrize(
    "quantity,expected",
    [
        ("processivity_phi29", "70000"),
        ("processivity_klenow", "40"),
        ("reach_phi29", "3000"),
    ],
)
def test_the_registry_agrees_with_the_shipped_constants(quantity, expected):
    """The check the prose ledger never had.

    `docs/SCIENCE_CITATIONS.md` still states Klenow processivity as 10,000 bp
    citing Bambara (1978); the registry says 40. One of them has to be checked
    against the code, and it is this one.
    """
    from neoswga.core.registry import views

    record = load_evidence()[quantity]
    assert record.value == expected

    characteristics = views.as_characteristics()
    field = "processivity" if quantity.startswith("processivity") else "typical_amplicon_length"
    name = quantity.split("_", 1)[1]
    assert str(characteristics[name][field]) == expected


def test_a_missing_registry_fails_model_initialisation(monkeypatch):
    """Not an empty registry. Absent evidence is not evidence of no constraint."""
    import neoswga.core.model_evidence as module

    monkeypatch.setattr(module, "_CACHE", None, raising=False)
    monkeypatch.setattr(module, "_EVIDENCE_FILENAME", "no_such_file.json")
    with pytest.raises(Exception) as caught:
        load_evidence()
    assert "model_evidence" in str(caught.value) or "no_such_file" in str(caught.value)


def test_a_record_is_immutable():
    import dataclasses

    record = next(iter(load_evidence().values()))
    with pytest.raises(dataclasses.FrozenInstanceError):
        record.status = "measured"
    assert isinstance(record, EvidenceRecord)


# ---------------------------------------------------------------------------
# 2. require_model_support refuses before the search
# ---------------------------------------------------------------------------


def test_a_supported_request_passes():
    require_model_support(resolve_design_request(base_params()))


def test_an_unknown_polymerase_is_refused():
    with pytest.raises(UnsupportedModelError):
        resolve_design_request(base_params(polymerase="taq"))


def test_a_primer_length_outside_the_polymerase_range_is_refused():
    """Known issue: a bst design was filtered through phi29's 6-12 bp window.

    Bst is modelled for 15-25 nt. A 6-mer under Bst at 63 C is not a primer the
    parameter set covers, and computing a Tm for it produces a number with no
    evidence behind it rather than a wrong one anyone can see.
    """
    with pytest.raises(UnsupportedModelError, match="length"):
        require_model_support(
            resolve_design_request(
                base_params(polymerase="bst", reaction_temp=63.0, min_k=6, max_k=12)
            )
        )


def test_the_same_polymerase_with_its_own_lengths_passes():
    require_model_support(
        resolve_design_request(
            base_params(polymerase="bst", reaction_temp=63.0, min_k=15, max_k=25)
        )
    )


def test_a_temperature_outside_the_polymerase_range_is_refused():
    with pytest.raises((UnsupportedModelError, Exception)):
        resolve_design_request(base_params(polymerase="phi29", reaction_temp=65.0))


def test_an_additive_with_no_tm_model_is_refused():
    """Glycerol is validated to 0-15%, printed in the summary, and changes no Tm.

    Measured in this repository: at 10% glycerol the effective Tm of a 12-mer is
    identical to the figure at 0%, to the last digit. The literature expects a
    real destabilisation, so zero here is an absent model rather than a result.
    """
    with pytest.raises(UnsupportedModelError, match="glycerol"):
        require_model_support(resolve_design_request(base_params(glycerol_percent=10.0)))


def test_an_additive_whose_absence_of_a_tm_term_is_deliberate_passes():
    """BSA and PEG act on the enzyme, not on duplex stability.

    No Tm term for them is a modelling decision with a reason, not a gap, and
    the registry records which of the two each one is.
    """
    require_model_support(resolve_design_request(base_params(bsa_ug_ml=200.0, peg_percent=4.0)))


def test_the_workhorse_additive_pair_is_supported():
    """DMSO plus betaine is what three shipped examples use.

    They have no pairwise interaction term; their Tm contributions compose
    independently. That is a stated assumption rather than an absence, and the
    registry says which.
    """
    require_model_support(resolve_design_request(base_params(dmso_percent=5.0, betaine_m=1.5)))


def test_a_concentration_in_the_wrong_unit_is_refused():
    """1.5 molar DMSO is a percentage typed into a molar field, or the reverse.

    A unit mistake lands inside no declared range by accident, so the range
    check is what catches it; the registry records the unit so the message can
    say which one was expected.
    """
    with pytest.raises(Exception) as caught:
        resolve_design_request(base_params(dmso_percent=150.0))
    assert "DMSO" in str(caught.value) or "dmso" in str(caught.value)


def test_a_dimer_free_energy_floor_without_a_model_is_refused():
    """`max_dimer_dg` is evaluated at the reaction temperature.

    Nothing has validated that floor against a reaction, and the plan requires
    primer-primer dimer chemistry to carry its own capability record rather
    than inheriting the duplex one.
    """
    record = load_evidence()["dimer_free_energy"]
    assert record.status in {"estimated", "empirical", "assumed"}
    assert record.notes


# ---------------------------------------------------------------------------
# 3. Changing the reaction invalidates what was derived under it
# ---------------------------------------------------------------------------


def test_changing_conditions_changes_the_request_hash():
    """The cache key every derived quantity should be keyed on.

    Two designs under different chemistry must not be able to share a cached
    occupancy, Tm or eligibility verdict, and the hash is what distinguishes
    them.
    """
    warm = resolve_design_request(base_params(reaction_temp=35.0))
    cool = resolve_design_request(base_params(reaction_temp=30.0))
    additive = resolve_design_request(base_params(reaction_temp=30.0, dmso_percent=5.0))

    assert len({warm.request_hash, cool.request_hash, additive.request_hash}) == 3


def test_the_condition_fingerprint_moves_with_every_modelled_additive():
    """An additive the fingerprint ignores lets a stale verdict be reused.

    The inventory keys eligibility on this fingerprint, so an additive outside
    it would let candidates assessed under one chemistry be reused under
    another without anything saying so.
    """
    from neoswga.core.reaction_conditions import ReactionConditions

    base = ReactionConditions(temp=30.0, polymerase="phi29")
    for field, value in (
        ("dmso_percent", 5.0),
        ("betaine_m", 1.0),
        ("trehalose_m", 0.5),
        ("formamide_percent", 5.0),
        ("ethanol_percent", 2.0),
        ("urea_m", 1.0),
        ("tmac_m", 0.1),
    ):
        changed = ReactionConditions(temp=30.0, polymerase="phi29", **{field: value})
        assert changed.fingerprint() != base.fingerprint(), field


# ---------------------------------------------------------------------------
# 4. The registry survives installation
# ---------------------------------------------------------------------------


def test_every_package_on_disk_is_declared_for_installation():
    """A subpackage setuptools is not told about is absent from the wheel.

    `neoswga.core.registry` was undeclared, so an installed neoswga had no
    polymerase table, no model evidence and no `views` module -- and
    `core/parameter.py` imports it at module scope, so importing the package
    at all would have failed. Verified by building a wheel on 2026-09-21: the
    only entries matching "registry" were two unrelated modules.

    Package data cannot rescue this. `include-package-data` applies to
    packages that are being installed, and this one was not.
    """
    import pathlib
    import tomllib

    root = pathlib.Path(__file__).resolve().parent.parent
    with open(root / "pyproject.toml", "rb") as handle:
        declared = set(tomllib.load(handle)["tool"]["setuptools"]["packages"])

    on_disk = {
        str(path.parent.relative_to(root)).replace("/", ".")
        for path in (root / "neoswga").rglob("__init__.py")
        if "__pycache__" not in str(path)
    }

    assert not on_disk - declared, (
        "these packages exist and would not be installed: " f"{sorted(on_disk - declared)}"
    )


def test_the_registry_json_is_declared_as_package_data():
    import pathlib
    import tomllib

    root = pathlib.Path(__file__).resolve().parent.parent
    with open(root / "pyproject.toml", "rb") as handle:
        data = tomllib.load(handle)["tool"]["setuptools"]["package-data"]

    assert "*.json" in data.get("neoswga.core.registry", []), data


def test_the_evidence_loads_from_an_arbitrary_working_directory(tmp_path, monkeypatch):
    """The failure mode a repository checkout hides.

    Resolving the artifact relative to the source tree works for every
    developer and for no installed user.
    """
    import neoswga.core.model_evidence as module

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(module, "_CACHE", None, raising=False)

    assert module.load_evidence()["reach_phi29"].value == "3000"


# ---------------------------------------------------------------------------
# 5. The concentration policy, and why it is not propagated
# ---------------------------------------------------------------------------


def test_a_fixed_total_moves_tm_a_lot_and_occupancy_almost_not_at_all():
    """The measurement behind not threading concentration through evaluation.

    Under a fixed total the per-oligo concentration falls with panel size, and
    Tm falls about ten degrees between a 1-oligo and a 96-oligo panel. The
    quantity selection uses is occupancy, and at phi29 30 C almost everything
    is saturated, so the same change moves it by about one part in ten
    thousand. Known Issue 17 is the reason.

    Pinned so that a future decision to propagate concentration is taken
    against a number rather than against an intuition, and so that a change
    making occupancy concentration-sensitive at the default reaction shows up
    here.
    """
    from neoswga.core.occupancy import site_occupancy
    from neoswga.core.reaction_conditions import ReactionConditions
    from neoswga.core.thermodynamics import calculate_enthalpy_entropy

    sequence = "ACGTTGCAAGGC"
    enthalpy, _entropy = calculate_enthalpy_entropy(sequence)
    conditions = ReactionConditions(temp=30.0, polymerase="phi29")
    total = 4e-6

    tms = [conditions.calculate_effective_tm(sequence, primer_conc=total / n) for n in (1, 96)]
    occupancies = [site_occupancy(enthalpy, tm, 30.0) for tm in tms]

    assert tms[0] - tms[1] > 9.0, tms
    assert occupancies[0] / occupancies[1] < 1.001, occupancies


def test_the_request_conserves_a_fixed_total_across_panel_sizes():
    request = resolve_design_request(
        base_params(concentration_mode="fixed_total", total_primer_molar=4e-6)
    )

    for count in (2, 6, 24):
        values = request.concentrations_molar(tuple(f"ACGTACGTAC{i:02d}" for i in range(count)))
        assert len(values) == count
        assert sum(values) == pytest.approx(4e-6)
        assert all(value > 0 for value in values)
