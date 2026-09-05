"""The audit record must actually be written, and must tell two runs apart.

Two defects, one of which I introduced while fixing the other.

`params_hash` covered thirteen keys and omitted everything that most changes the
delivered set: the optimizer, the dimer thresholds, the coverage reach, the Tm
window, the mismatch allowance. Two runs differing only in
`--optimization-method` produced an identical hash, so the record could not
distinguish the runs it exists to distinguish.

And `write_audit_trail` is wrapped in a bare `try/except`, so when it was split
into its own module without the `parameter` import it uses, every call failed
silently. The CSV and the summary kept being written and the audit file simply
stopped updating -- visible only by comparing file timestamps. This file asserts
the record is written at all, because a whole feature disappearing behind a
swallowed exception is the failure mode that has to be caught by something.
"""

import json
from dataclasses import replace

import pytest

from neoswga.core.base_optimizer import OptimizationResult, OptimizationStatus, PrimerSetMetrics


def _result(optimizer_name):
    return OptimizationResult(
        primers=("AAACCCGGGT", "ACCCGGGTTT"),
        score=1.0,
        metrics=replace(PrimerSetMetrics.empty(), fg_coverage=0.9),
        status=OptimizationStatus.SUCCESS,
        optimizer_name=optimizer_name,
        iterations=1,
    )


def _write(tmp_path, optimizer_name):
    from neoswga.core.step4_audit_trail import write_audit_trail

    out = tmp_path / f"{optimizer_name}.csv"
    write_audit_trail(str(out), _result(optimizer_name))
    return tmp_path / f"{optimizer_name}_audit.json"


def test_the_audit_record_is_written(tmp_path):
    """It stopped being written for a while and nothing noticed."""
    path = _write(tmp_path, "network")
    assert path.exists(), (
        "no audit record was produced. `write_audit_trail` swallows every "
        "exception, so a missing import or a renamed attribute silently "
        "disables it while the rest of the step keeps working."
    )


def test_the_record_names_the_optimizer_that_produced_the_set(tmp_path):
    payload = json.loads(_write(tmp_path, "clique").read_text())
    assert payload["optimizer"] == "clique"
    assert payload["parameters"]["optimizer"] == "clique"


def test_two_methods_do_not_share_a_params_hash(tmp_path):
    """The optimizer is the largest single determinant of the result."""
    a = json.loads(_write(tmp_path, "network").read_text())["params_hash"]
    b = json.loads(_write(tmp_path, "clique").read_text())["params_hash"]

    assert a != b, (
        f"both runs hashed to {a}; the record cannot distinguish a network run "
        "from a clique run on the same params.json"
    )


@pytest.mark.parametrize(
    "attr,value",
    [
        ("max_dimer_bp", 2),
        ("coverage_reach", 8_000),
        ("max_mismatches", 3),
        ("min_tm", 12.0),
    ],
)
def test_a_setting_that_changes_the_set_changes_the_hash(tmp_path, monkeypatch, attr, value):
    """Each of these alters which primers come back, and none was hashed."""
    from neoswga.core import parameter

    before = json.loads(_write(tmp_path, "network").read_text())["params_hash"]
    monkeypatch.setattr(parameter, attr, value, raising=False)
    after = json.loads(_write(tmp_path, "network").read_text())["params_hash"]

    assert before != after, (
        f"changing {attr} to {value} left the params_hash at {before}; the "
        "record would report two materially different runs as the same"
    )
