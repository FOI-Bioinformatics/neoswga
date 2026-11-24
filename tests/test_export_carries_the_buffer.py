"""Exported and reported Tm are corrected for the recorded monovalent ions.

`export` built its `ReactionConditions` from the polymerase, temperature,
Mg2+ and the protocol additives only, so `na_conc`, `k_conc`, `nh4_conc`,
`dntp_conc` and `dtt_mm`, all in the run manifest's `effective_conditions`,
were dropped: na_conc 300 against 50 moves TGTACTGCCAA from 42.8 C to 45.2 C
and the export printed 42.8. The report and `interpret` did the same, and
`interpret` also substituted 50 mM Na+ and 2.5 mM Mg2+ when params.json said
nothing (2026-11-23 sweep, tier 2).
"""

import json

import pytest

from neoswga.core.run_manifest import write_manifest

PRIMER = "TGTACTGCCAA"


@pytest.fixture
def results(tmp_path):
    directory = tmp_path / "results"
    directory.mkdir()
    (directory / "step4_improved_df.csv").write_text(f"primer,set_index\n{PRIMER},0\n")
    (directory / "params.json").write_text(json.dumps({"polymerase": "phi29"}))
    write_manifest(
        "optimize",
        str(directory),
        effective_conditions={
            "polymerase": "phi29",
            "reaction_temp": 30.0,
            "na_conc": 300.0,
            "k_conc": 20.0,
            "dntp_conc": 1.0,
        },
    )
    return directory


def test_the_exporter_corrects_for_the_recorded_buffer(results):
    from neoswga.core.export import PrimerExporter

    conditions = PrimerExporter.from_results_dir(str(results)).conditions()
    assert conditions.na_conc == 300.0
    assert conditions.k_conc == 20.0
    assert conditions.dntp_conc == 1.0


def test_the_report_reconstructs_the_recorded_buffer(results):
    from neoswga.core.report.metrics import _conditions_from_params

    conditions = _conditions_from_params({}, results)
    assert conditions.na_conc == 300.0


def test_interpret_uses_the_recorded_reaction(results, monkeypatch):
    from neoswga.core import mechanistic_model
    from neoswga.core.results_interpreter import ResultsInterpreter

    seen = {}

    class _Model:
        def __init__(self, conditions):
            seen["conditions"] = conditions

        def calculate_effects(self, *_a, **_k):
            raise RuntimeError("stop after construction")

    monkeypatch.setattr(mechanistic_model, "MechanisticModel", _Model)
    ResultsInterpreter(str(results))._estimate_enrichment_fold([{"primer": PRIMER}])
    assert seen["conditions"].na_conc == 300.0
