"""Coverage for neoswga.core.experimental_tracker (prediction/outcome log)."""

from types import SimpleNamespace

import pytest

from neoswga.core.experimental_tracker import (
    ExperimentalOutcome,
    ExperimentalTracker,
)


def _prediction(enrichment=100.0, coverage=0.8, confidence=0.7):
    return SimpleNamespace(
        predicted_enrichment=enrichment,
        predicted_coverage=coverage,
        confidence_score=confidence,
    )


@pytest.fixture
def tracker(tmp_path):
    return ExperimentalTracker(log_path=str(tmp_path / "exp.json"))


def test_record_prediction_returns_id(tracker):
    eid = tracker.record_prediction(["ATCG", "GCTA"], _prediction())
    assert isinstance(eid, str) and eid
    exp = tracker.get_experiment(eid)
    assert exp is not None
    assert exp.predicted_enrichment == 100.0
    assert exp.has_outcome is False


def test_record_outcome_and_errors(tracker):
    eid = tracker.record_prediction(["ATCG"], _prediction(enrichment=100.0, coverage=0.8))
    tracker.record_outcome(experiment_id=eid, actual_enrichment=80.0, actual_coverage=0.6)
    exp = tracker.get_experiment(eid)
    assert exp.has_outcome is True
    # errors are |predicted - actual|-style magnitudes (non-negative)
    assert exp.enrichment_error is not None and exp.enrichment_error >= 0
    assert exp.coverage_error is not None and exp.coverage_error >= 0
    assert len(tracker.get_experiments_with_outcomes()) == 1


def test_calibration_report(tracker):
    for i in range(3):
        eid = tracker.record_prediction([f"PRIMER{i}"], _prediction(enrichment=100.0))
        tracker.record_outcome(experiment_id=eid, actual_enrichment=90.0, actual_coverage=0.7)
    report = tracker.get_calibration_report()
    assert report is not None
    d = report.to_dict()
    assert isinstance(d, dict)


def test_persistence_and_clear(tmp_path):
    path = str(tmp_path / "exp.json")
    t1 = ExperimentalTracker(log_path=path)
    t1.record_prediction(["ATCG"], _prediction())
    # reload sees the experiment
    t2 = ExperimentalTracker(log_path=path)
    assert len(t2.list_experiments()) >= 1
    t2.clear()
    t3 = ExperimentalTracker(log_path=path)
    assert len(t3.list_experiments()) == 0


def test_outcome_dataclass_roundtrip():
    o = ExperimentalOutcome(
        experiment_id="e1",
        primer_set=["ATCG"],
        predicted_enrichment=100.0,
        predicted_coverage=0.8,
        confidence_score=0.7,
    )
    restored = ExperimentalOutcome.from_dict(o.to_dict())
    assert restored.experiment_id == "e1"
    assert restored.primer_set == ["ATCG"]
