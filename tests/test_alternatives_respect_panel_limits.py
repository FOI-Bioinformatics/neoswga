"""An alternative set that breaks a configured panel limit is not delivered.

`collect_alternative_sets` appended every non-empty, distinct alternative. The
contract it runs through marks a limit violation only as a PARTIAL status and
keeps the panel, so with `min_selectivity_density` set, set 0 could meet it
after repair while set 2 missed it and was still written to
step4_improved_df.csv (2026-11-23 sweep, tier 2). `export` refuses sets other
than 0, which limited the reach; anyone reading the CSV still saw it.
"""

from types import SimpleNamespace

from neoswga.core import search_control

CANDIDATES = [f"P{i:02d}" for i in range(12)]


def _run(monkeypatch, proposals, bad):
    calls = iter(proposals)

    def _alternative(optimizer, pool, target_size, budget, constraints, through_contract):
        return SimpleNamespace(primers=next(calls))

    class _Objective:
        def violations(self, primers):
            return ("min_selectivity_density",) if tuple(primers) in bad else ()

    monkeypatch.setattr(search_control, "_search_alternative", _alternative)
    monkeypatch.setattr(
        "neoswga.core.panel_refinement.objective_for_optimizer", lambda optimizer, *a: _Objective()
    )
    primary = SimpleNamespace(primers=("P00", "P01"))
    return search_control.collect_alternative_sets(
        primary, object(), CANDIDATES, target_size=2, max_sets=4, max_iterations=4
    )


def test_a_violating_alternative_is_not_delivered(monkeypatch):
    sets = _run(monkeypatch, [("P02", "P03"), ("P04", "P05"), ("P06", "P07")], bad={("P04", "P05")})
    assert ("P04", "P05") not in sets
    assert sets == [("P00", "P01"), ("P02", "P03"), ("P06", "P07")]


def test_without_limits_every_alternative_stands(monkeypatch):
    sets = _run(monkeypatch, [("P02", "P03"), ("P04", "P05"), ("P06", "P07")], bad=set())
    assert len(sets) == 4
