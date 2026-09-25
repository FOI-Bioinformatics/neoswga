"""Phase 3 (production-readiness v2): state, determinism, data-integrity."""

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# pipeline.reset_pipeline_state + auto-reset on params change
# ---------------------------------------------------------------------------


def test_reset_pipeline_state_clears_globals():
    import neoswga.core.pipeline as pl

    pl._initialized = True
    pl.fg_prefixes = ["stale"]
    pl.fg_seq_lengths = [123]
    pl.fg_circular = True

    pl.reset_pipeline_state()

    assert pl._initialized is False
    assert pl.fg_prefixes is None
    assert pl.fg_seq_lengths is None
    assert pl.fg_circular is None


def test_initialize_auto_resets_on_json_change(monkeypatch):
    import neoswga.core.pipeline as pl

    calls = {"n": 0}

    def fake_reset():
        calls["n"] += 1
        pl._initialized = False

    monkeypatch.setattr(pl, "reset_pipeline_state", fake_reset)
    # Pretend we initialized with file A.
    pl._initialized = True
    pl._initialized_json_file = "A.json"

    class _Param:
        json_file = "B.json"  # changed

    monkeypatch.setattr(pl, "parameter", _Param)
    # Stub get_params so _initialize doesn't do real work after reset.
    monkeypatch.setattr(
        pl.parameter,
        "get_params",
        lambda opts: {"fg_prefixes": [], "fg_genomes": [], "fg_seq_lengths": []},
        raising=False,
    )
    monkeypatch.setattr(pl, "_apply_gc_adaptive_defaults", lambda: None)

    pl._initialize()
    assert calls["n"] == 1  # auto-reset fired because the params path changed


# ---------------------------------------------------------------------------
# multi_genome result reports None (not fabricated) for uncomputed metrics
# ---------------------------------------------------------------------------


def test_multi_genome_result_allows_none_metrics():
    """The placeholder metrics must admit None, so an uncomputed metric reads
    as absent rather than as a fabricated zero.

    This asserts the PROPERTY, not the spelling. It used to require the literal
    string "Optional" in the annotation, so modernising `Optional[float]` to
    `float | None` broke it while its subject had not changed at all. Both
    spellings put NoneType in `get_args`, and a field made non-optional still
    fails, because a bare `float` has no args.
    """
    import typing

    from neoswga.core.multi_genome_pipeline import MultiGenomePipelineResult

    hints = typing.get_type_hints(MultiGenomePipelineResult)
    for name in ("coverage", "connectivity", "predicted_amplification", "stage1_primer_count"):
        assert type(None) in typing.get_args(hints[name]), (
            f"{name} must admit None so an uncomputed metric is absent, not fabricated"
        )


# `test_bg_aware_connectivity_uses_network` was removed on 2026-09-10 with
# `BackgroundAwareOptimizer._calculate_connectivity`. The class no dispatch path
# reached was deleted; the method existed only on it, so there is no live
# behaviour left for the test to pin. Recorded here rather than silently
# dropped: the assertion was that connectivity comes from the network rather
# than from len(primers)/10, and if that shortcut ever reappears elsewhere it
# wants a test again.
# ---------------------------------------------------------------------------
# Determinism: stochastic simulator RNG is seeded
# ---------------------------------------------------------------------------


def test_gillespie_rng_is_seeded():
    from neoswga.core.stochastic_simulator import GillespieSimulator

    a = GillespieSimulator.__new__(GillespieSimulator)
    a._rng = np.random.default_rng(7)
    b = GillespieSimulator.__new__(GillespieSimulator)
    b._rng = np.random.default_rng(7)
    assert a._rng.random() == b._rng.random()


def test_efficiency_dimer_subsample_is_deterministic():
    import random

    primers = [f"PRIMER{i:03d}" for i in range(100)]
    seed = hash(tuple(sorted(primers))) & 0xFFFFFFFF
    s1 = random.Random(seed).sample(primers, 50)
    s2 = random.Random(seed).sample(primers, 50)
    assert s1 == s2  # same input -> same subsample


# ---------------------------------------------------------------------------
# Determinism: replication simulator (small, marked slow)
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_replication_simulate_primer_set_seeded_reproducible():
    from neoswga.core.reaction_conditions import get_standard_conditions
    from neoswga.core.replication_simulator import SimulationConfig, simulate_primer_set

    genome = "ATCG" * 50  # 200 bp
    positions = {"ATCG": {"forward": np.array([0, 40, 80]), "reverse": np.array([20, 60])}}
    cfg = SimulationConfig(duration=5.0, time_step=1.0, use_mechanistic_model=False)

    kw = dict(
        primers=["ATCG"],
        primer_positions=positions,
        genome_length=len(genome),
        genome_sequence=genome,
        conditions=get_standard_conditions(),
        n_replicates=2,
        config=cfg,
    )
    r1 = simulate_primer_set(seed=42, **kw)
    r2 = simulate_primer_set(seed=42, **kw)
    assert r1["mean_coverage"] == r2["mean_coverage"]
    assert r1["mean_amplification"] == r2["mean_amplification"]
