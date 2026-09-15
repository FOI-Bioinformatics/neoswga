"""One oligo, one reaction, one melting temperature.

Task 1 of the condition-aware pool design plan. `ReactionConditions.
calculate_effective_tm` is the canonical calculation, and two optimizer paths
reconstructed their own instead.

`NetworkOptimizer._get_primer_tm` called the salt-corrected Tm with sodium
alone, so magnesium took the low-level default of zero. On a 12-mer at 50 mM Na
that is 38.07 C against 46.96 C at the 10 mM magnesium a phi29 reaction runs at,
nearly nine degrees, on the quantity that weights every edge in its scoring. It
also dropped potassium, ammonium, dNTP and the per-oligo concentration.

`ThermodynamicFilter.analyze_primer` applied salt but no additive correction. A
longer primer brought into the window by DMSO or betaine passes the
additive-aware stage-2 gate and is then rejected again by this secondary screen
on its uncorrected Tm, which undoes the chemistry the user configured.

Cache identity is part of the same problem: a Tm cached under a sequence alone
survives a change of reaction.
"""

from types import SimpleNamespace

import pytest

from neoswga.core.network_optimizer import NetworkOptimizer
from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.thermodynamic_filter import ThermodynamicCriteria, ThermodynamicFilter

PRIMER = "ACGTACGTACGT"

# Each varies one field away from the same base reaction.
VARIATIONS = [
    ("baseline", {}),
    ("magnesium", {"mg_conc": 10.0}),
    ("potassium", {"k_conc": 50.0}),
    ("ammonium", {"nh4_conc": 10.0}),
    ("dntp", {"dntp_conc": 0.8}),
    ("primer_conc", {"primer_conc": 5e-6}),
    ("dmso", {"dmso_percent": 8.0}),
    ("betaine", {"betaine_m": 1.5}),
]


def _conditions(**overrides):
    base = {"temp": 30.0, "polymerase": "phi29", "na_conc": 50.0, "mg_conc": 0.0}
    base.update(overrides)
    return ReactionConditions(**base)


@pytest.mark.parametrize("name,overrides", VARIATIONS, ids=[v[0] for v in VARIATIONS])
def test_the_network_path_agrees_with_the_condition_object(name, overrides):
    conditions = _conditions(**overrides)
    holder = SimpleNamespace(conditions=conditions, _tm_cache={})

    assert NetworkOptimizer._get_primer_tm(holder, PRIMER) == pytest.approx(
        conditions.calculate_effective_tm(PRIMER)
    ), f"network scoring computes a different Tm when {name} is set"


@pytest.mark.parametrize("name,overrides", VARIATIONS, ids=[v[0] for v in VARIATIONS])
def test_the_secondary_screen_agrees_with_the_condition_object(name, overrides):
    conditions = _conditions(**overrides)
    screen = ThermodynamicFilter(
        ThermodynamicCriteria(na_conc=conditions.na_conc, mg_conc=conditions.mg_conc),
        conditions=conditions,
    )

    assert screen.analyze_primer(PRIMER).tm == pytest.approx(
        conditions.calculate_effective_tm(PRIMER)
    ), f"the secondary screen computes a different Tm when {name} is set"


def test_an_additive_admitted_primer_is_not_rejected_by_the_secondary_screen():
    """The concrete case: chemistry admits a primer, then the screen undoes it.

    A strong additive load lowers Tm. A primer sitting above the window on its
    uncorrected Tm is brought inside it by the correction, which is the whole
    point of configuring additives. The screen must judge it on the same Tm the
    gate that admitted it used.
    """
    conditions = _conditions(mg_conc=10.0, dmso_percent=10.0, betaine_m=1.5)
    corrected = conditions.calculate_effective_tm(PRIMER)
    uncorrected = _conditions(mg_conc=10.0).calculate_effective_tm(PRIMER)
    assert corrected < uncorrected - 1.0, "fixture does not exercise the correction"

    # A window that admits the corrected Tm and excludes the uncorrected one.
    criteria = ThermodynamicCriteria(
        na_conc=conditions.na_conc,
        mg_conc=conditions.mg_conc,
        min_tm=corrected - 2.0,
        max_tm=(corrected + uncorrected) / 2.0,
    )
    screen = ThermodynamicFilter(criteria, conditions=conditions)

    result = screen.analyze_primer(PRIMER)
    assert not [r for r in result.failure_reasons if "Tm" in r], (
        f"rejected on Tm despite the configured additives admitting it: "
        f"{result.failure_reasons}"
    )


def test_the_fingerprint_changes_with_chemistry_and_not_otherwise():
    assert _conditions().fingerprint() == _conditions().fingerprint()
    assert _conditions().fingerprint() != _conditions(mg_conc=10.0).fingerprint()
    assert _conditions().fingerprint() != _conditions(dmso_percent=5.0).fingerprint()
    assert _conditions().fingerprint() != _conditions(primer_conc=5e-6).fingerprint()


def test_the_fingerprint_covers_the_model_version():
    """So a coefficient revision cannot reuse a cache computed under the old one."""
    from neoswga.core.reaction_conditions import CONDITION_MODEL_VERSION

    assert CONDITION_MODEL_VERSION in _conditions().fingerprint()


def test_a_cached_tm_does_not_survive_a_change_of_reaction():
    """The cache was keyed by sequence alone."""
    weak, strong = _conditions(), _conditions(dmso_percent=10.0)
    holder = SimpleNamespace(conditions=weak, _tm_cache={})

    first = NetworkOptimizer._get_primer_tm(holder, PRIMER)
    holder.conditions = strong
    second = NetworkOptimizer._get_primer_tm(holder, PRIMER)

    assert second != pytest.approx(first)
    assert second == pytest.approx(strong.calculate_effective_tm(PRIMER))


def test_the_screen_passes_the_resolved_conditions_to_the_filter():
    """Otherwise the canonical Tm added above never reaches the screen."""
    import ast
    import pathlib

    source = pathlib.Path("neoswga/core/hybrid_thermo_screen.py").read_text()
    tree = ast.parse(source)
    construction = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "ThermodynamicFilter"
    ]
    assert construction, "ThermodynamicFilter is no longer constructed here"
    for call in construction:
        assert any(kw.arg == "conditions" for kw in call.keywords), (
            "the screen builds a ThermodynamicFilter without the resolved "
            "conditions, so it judges primers on an uncorrected Tm"
        )


def test_an_empty_screen_result_is_not_replaced_by_its_input():
    """Returning the input on an empty result turns a QC failure into a pass.

    The screen logged a warning and returned the unfiltered candidates, so a
    configuration under which nothing passes thermodynamic QC proceeded with
    every candidate it had just rejected. That is not a screen.
    """
    from neoswga.core.hybrid_thermo_screen import ThermoScreenMixin

    class _NothingPasses(ThermoScreenMixin):
        polymerase = "phi29"
        poly_config = SimpleNamespace(reaction_temp=30.0)
        conditions = _conditions(mg_conc=10.0)
        max_dimer_bp = 3
        _thermo_filter_cache = None

        def _thermo_criteria(self):
            # A window no primer can satisfy.
            return ThermodynamicCriteria(min_tm=999.0, max_tm=1000.0)

    assert _NothingPasses()._thermo_filter_candidates([PRIMER, "CCCCCCCCCCCC"], verbose=False) == []


def test_the_screen_cache_does_not_survive_a_change_of_chemistry():
    from neoswga.core.hybrid_thermo_screen import ThermoScreenMixin

    class _Screen(ThermoScreenMixin):
        polymerase = "phi29"
        poly_config = SimpleNamespace(reaction_temp=30.0)
        conditions = _conditions()
        max_dimer_bp = 3
        _thermo_filter_cache = None

        def _thermo_criteria(self):
            return ThermodynamicCriteria(
                na_conc=self.conditions.na_conc, mg_conc=self.conditions.mg_conc
            )

    screen = _Screen()
    screen._thermo_filter_with_cache([PRIMER], verbose=False)
    assert screen._thermo_filter_cache is not None
    cached_under = screen._thermo_filter_cache[0]

    screen.conditions = _conditions(dmso_percent=10.0)
    screen._thermo_filter_with_cache([PRIMER], verbose=False)

    assert (
        screen._thermo_filter_cache[0] != cached_under
    ), "the screen reused a verdict computed under a different reaction"
