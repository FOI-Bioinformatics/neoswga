"""The sequence heuristics, made configurable and measured rather than rescaled.

`filter.py` applies three rules inherited from PCR primer design -- a GC clamp
over the last 5 bases, a ban on GGG/CCC at the 3' end, and a homopolymer-run
limit. They were hardcoded constants, and they are length-blind: "the last 5
bases" is 83% of a 6-mer but 42% of a 12-mer, and they come from a regime of
18-25mers at 55-65 C. Neither swga 1.0 nor swga 2.0 applies any of them.

The plan was to scale them to primer length. Measurement says not to.

Rejection rates turn out to be flat across length -- the GC clamp rejects
18-19% of random k-mers at every k from 6 to 12 -- so the length-blindness does
not show up as differential harshness. What does vary is redundancy: 79% of the
clamp's rejections at k=7 would already have failed the GC-content filter,
against 28% at k=12. So at short lengths it mostly duplicates an existing
filter and at long lengths it does independent work, which is the argument for
scaling the window.

Scaling it makes the design worse. Against the 141 pcDNA primers with measured
amplification, a window scaled to `min(5, max(2, k//2))` newly admits 20
primers whose median measured amplification is 1.55x, against a pool median of
10.15x -- only 4 of the 20 beat the pool, and 12 of them are 7-mers. The fixed
5-base window earns its place at short lengths precisely by acting as the extra
constraint that redundancy analysis makes it look like it is not.

So the constants become parameters, with the defaults reproducing today's
behaviour exactly, and the rules keep their shape. A user who wants to test the
other choice can; nothing moves silently for one who does not.
"""

import pytest

from neoswga.core import filter as filter_mod
from neoswga.core import parameter


@pytest.fixture(autouse=True)
def restore():
    """These are module globals; leave them as found."""
    snapshot = {
        name: getattr(parameter, name, None)
        for name in ("max_homopolymer_run", "gc_clamp_window", "max_gc_in_clamp")
    }
    yield
    for name, value in snapshot.items():
        if value is None:
            if hasattr(parameter, name):
                delattr(parameter, name)
        else:
            setattr(parameter, name, value)


# ----------------------------------------------------------------------
# Defaults reproduce today's behaviour
# ----------------------------------------------------------------------


def test_defaults_are_the_values_that_were_hardcoded():
    """Nothing moves for a user who does not ask for it."""
    assert filter_mod.resolve_homopolymer_run() == 5
    assert filter_mod.resolve_gc_clamp() == (5, 3)


def test_a_known_clamp_failure_still_fails():
    """GGCCG in the last five bases is 5 G/C, over the limit of 3."""
    assert filter_mod._fails_gc_clamp("AATTGGCCG")


def test_a_known_clamp_pass_still_passes():
    """AGCTA is 2 G/C in the last five -- inside the 1-3 band.

    The band has a floor as well as a ceiling: a 3' end with no G/C at all is
    too unstable to prime reliably, so 0 fails too.
    """
    assert not filter_mod._fails_gc_clamp("TTTTTAGCTA")
    assert filter_mod._fails_gc_clamp("GGCCCAATTA"), "no G/C in the last five should fail"


# ----------------------------------------------------------------------
# They are configurable
# ----------------------------------------------------------------------


def test_the_clamp_window_can_be_widened(monkeypatch):
    monkeypatch.setattr(parameter, "gc_clamp_window", 3, raising=False)
    monkeypatch.setattr(parameter, "max_gc_in_clamp", 3, raising=False)

    # Only the last three bases are considered now, and TTA carries no G/C.
    assert filter_mod.resolve_gc_clamp() == (3, 3)
    # Last three bases are GTA: one G/C, inside the band.
    assert not filter_mod._fails_gc_clamp("GGCCGGCGTA")


def test_the_clamp_threshold_can_be_tightened(monkeypatch):
    monkeypatch.setattr(parameter, "max_gc_in_clamp", 1, raising=False)

    # AAGGA has two G/C in the last five, which now exceeds the limit of 1.
    assert filter_mod._fails_gc_clamp("TTTTTAAGGA")


def test_the_homopolymer_limit_is_configurable(monkeypatch):
    assert not filter_mod._has_homopolymer_run("AAAATTTT")

    monkeypatch.setattr(parameter, "max_homopolymer_run", 4, raising=False)
    assert filter_mod._has_homopolymer_run("AAAATTTT")


def test_homopolymer_patterns_follow_the_configured_length(monkeypatch):
    """The patterns were a module constant derived once at import, so a
    configured value would have had no effect on them."""
    monkeypatch.setattr(parameter, "max_homopolymer_run", 3, raising=False)

    assert filter_mod._has_homopolymer_run("CAAAC")
    assert not filter_mod._has_homopolymer_run("CAAC")


# ----------------------------------------------------------------------
# Why the shape was not changed
# ----------------------------------------------------------------------


def test_the_clamp_rejects_uniformly_across_primer_length():
    """The length-blindness does not manifest as differential harshness.

    A rule over "the last five bases" rejects the same fraction of random
    sequence whatever surrounds those bases, so scaling it cannot be justified
    by short primers being over-filtered -- they are not.
    """
    import random

    rng = random.Random(7)
    rates = []
    for k in range(6, 13):
        seqs = ["".join(rng.choice("ACGT") for _ in range(k)) for _ in range(4000)]
        rates.append(sum(1 for s in seqs if filter_mod._fails_gc_clamp(s)) / len(seqs))

    assert max(rates) - min(rates) < 0.05, f"rejection rates vary by length: {rates}"


def test_scaling_the_window_would_admit_worse_amplifiers():
    """The measurement that settled it, pinned so the conclusion is checkable.

    Against the only benchmark with per-primer measured amplification, a
    length-scaled clamp newly admits primers amplifying at a median of 1.55x
    where the pool median is 10.15x. The refinement that looks principled makes
    the design worse.
    """
    import json
    import statistics
    from pathlib import Path

    data = (
        Path(__file__).parent / "validation" / "data" / "dwivedi_yu_2023_plasmid_amplification.json"
    )
    payload = json.loads(data.read_text())

    amp = {}
    for row in payload["measurements"]:
        if row["designed_for"] == "pcDNA":
            amp.setdefault(row["primer"], []).append(row["on_target_amp"])
    amp = {p: statistics.mean(v) for p, v in amp.items()}
    pool_median = statistics.median(amp.values())

    def fails_scaled(primer):
        window = min(5, max(2, len(primer) // 2))
        return sum(1 for b in primer[-window:] if b in "GC") > 3

    newly_admitted = [p for p in amp if filter_mod._fails_gc_clamp(p) and not fails_scaled(p)]

    assert newly_admitted, "the scaled window should differ from the fixed one"
    assert statistics.median(amp[p] for p in newly_admitted) < pool_median, (
        "scaling the clamp admits primers that amplify better than the pool, "
        "which would reverse the decision to keep the fixed window"
    )


def test_a_mocked_parameter_module_does_not_silently_reconfigure_the_rules(monkeypatch):
    """The trap these resolvers fell into on first writing.

    Several test modules replace `neoswga.core.filter.parameter` with a
    MagicMock, which answers every attribute with a truthy stub whose `int()`
    is 1. A `getattr(...) or default` therefore resolved the GC clamp to a
    1-base window with a 1-base limit and quietly failed eight unrelated
    filter tests. Values are type-checked, not truthiness-checked.
    """
    from unittest.mock import MagicMock

    monkeypatch.setattr(filter_mod, "parameter", MagicMock())

    assert filter_mod.resolve_gc_clamp() == (5, 3)
    assert filter_mod.resolve_homopolymer_run() == 5


def test_a_nonsense_configured_value_falls_back_rather_than_applying(monkeypatch):
    """A string or a negative from a hand-edited params.json would otherwise
    produce a quietly different filter instead of an error."""
    for bad in ("five", -1, None, 0):
        monkeypatch.setattr(parameter, "gc_clamp_window", bad, raising=False)
        assert filter_mod.resolve_gc_clamp()[0] == 5, f"{bad!r} was applied"
