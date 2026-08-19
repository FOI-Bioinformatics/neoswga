"""The bench protocol has to list what the design was made for.

`export_protocol` builds the markdown a user takes to the bench. It read three
additives from params.json -- betaine, DMSO, trehalose -- and printed those
three. A design optimized with urea and TMAC produced a protocol whose additives
section read "No additives", and a run set up from that protocol is not the
reaction the primers were selected under.

Of the eleven additives the schema accepts, eight carry a melting-temperature
term and therefore changed which primers were chosen. Silently dropping five of
them between the design and the bench is the one defect in this area that leaves
the software entirely.

The magnesium default was wrong in the same place and for the same reason as
everywhere else: `mg_conc: float = 2.5` is a PCR concentration, retired from the
model when the mis-attributed citation behind it was removed
(registry/INCONSISTENCIES.md #6). A protocol is the last place that should still
print it -- 2.5 mM is a quarter of the standard phi29 buffer, and a reaction set
up from it will underperform for a reason nobody will look for.
"""

import json

import pytest

# Every additive the schema accepts that has a Tm model, with a value in range.
ADDITIVES = {
    "betaine_m": 1.0,
    "dmso_percent": 5.0,
    "trehalose_m": 0.3,
    "formamide_percent": 3.0,
    "ethanol_percent": 2.0,
    "urea_m": 1.0,
    "tmac_m": 0.05,
    "propanediol_m": 1.0,
}

# Accepted, modelled through the enzyme rather than Tm, and equally load-bearing
# at the bench.
NON_TM_ADDITIVES = {
    "glycerol_percent": 5.0,
    "peg_percent": 5.0,
    "bsa_ug_ml": 200.0,
}


@pytest.mark.parametrize("name,value", sorted(ADDITIVES.items()))
def test_every_tm_additive_reaches_the_protocol(name, value):
    from neoswga.core.export import generate_protocol

    text = generate_protocol(["ATGCATGCATGC"], **{name: value})

    assert str(value) in text, f"{name}={value} is missing from the protocol"
    assert "No additives" not in text


@pytest.mark.parametrize("name,value", sorted(NON_TM_ADDITIVES.items()))
def test_every_enzyme_additive_reaches_the_protocol(name, value):
    from neoswga.core.export import generate_protocol

    text = generate_protocol(["ATGCATGCATGC"], **{name: value})

    assert str(value) in text, f"{name}={value} is missing from the protocol"


def test_no_additives_is_still_reported_when_there_are_none():
    from neoswga.core.export import generate_protocol

    assert "No additives" in generate_protocol(["ATGCATGCATGC"])


def test_a_params_file_carries_its_additives_through(tmp_path):
    """End to end, which is the path a user actually takes."""
    from neoswga.core.export import PrimerExporter

    results = tmp_path / "results"
    results.mkdir()
    (results / "step4_improved_df.csv").write_text("primer\nATGCATGCATGC\nGCTAGCTAGCTA\n")
    params = dict(ADDITIVES)
    params.update(NON_TM_ADDITIVES)
    params.update({"polymerase": "equiphi29", "reaction_temp": 42.0, "mg_conc": 10.0})
    (results / "params.json").write_text(json.dumps(params))

    text = PrimerExporter.from_results_dir(str(results)).generate_protocol()

    for name, value in {**ADDITIVES, **NON_TM_ADDITIVES}.items():
        assert str(value) in text, f"{name} did not survive the params.json round trip"


def test_the_protocol_magnesium_default_is_the_isothermal_buffer():
    """2.5 mM is PCR. A protocol printing it would have a user set up a reaction
    at a quarter of the phi29 buffer's magnesium."""
    from neoswga.core.export import generate_protocol
    from neoswga.core.parameter import default_mg_conc

    for polymerase in ("phi29", "equiphi29", "bst", "klenow"):
        text = generate_protocol(["ATGCATGCATGC"], polymerase=polymerase)
        expected = default_mg_conc(polymerase)
        assert f"{expected} mM" in text, polymerase


def test_a_params_file_without_magnesium_gets_the_buffer_value(tmp_path):
    from neoswga.core.export import PrimerExporter
    from neoswga.core.parameter import default_mg_conc

    results = tmp_path / "results"
    results.mkdir()
    (results / "step4_improved_df.csv").write_text("primer\nATGCATGCATGC\n")
    (results / "params.json").write_text(json.dumps({"polymerase": "equiphi29"}))

    exporter = PrimerExporter.from_results_dir(str(results))

    assert exporter.mg_conc == default_mg_conc("equiphi29")
