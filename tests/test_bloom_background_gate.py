"""The Bloom background path must actually filter.

`neoswga build-filter` exists so a host-sized background (human, 3 Gbp) can be
screened without holding an exact k-mer index in memory. On that path
`get_bg_rates_via_bloom` had no counts to report -- a Bloom filter answers
presence, not frequency -- so it emitted a sentinel:

    primer_to_count[primer] = max_bg_matches + 1      # 11 by default

The gate downstream is a FREQUENCY:

    bg_bool = bg_count / bg_total_length < scaled_max_bg

A fixed count divided by a growing genome length shrinks without limit, so the
filter got more permissive as the background got larger -- the exact opposite of
what it is for:

    background      sentinel / length      vs max_bg_freq 5e-6
    6 kb plasmid    11 / 6.2e3 = 1.8e-3    rejected
    3 Gbp human     11 / 3.0e9 = 3.7e-9    PASSES

And since an absent primer scores 0, which also passes, *both* branches passed
on a large background: the background filter did nothing at all, silently, at
precisely the scale it was written for. Any specificity work built on top of it
would have been built on a filter that rejected nothing.

The repair is not a better sentinel. Presence/absence cannot answer a frequency
question, and inventing a number to stand in for one is what caused this. The
`SampledGenomeIndex` alongside it can answer it -- `build-filter` already writes
`bg_sampled.pkl` next to `bg_bloom.pkl` -- so the path now finds that index, and
refuses to run without it rather than reporting a number it does not have.
"""

import pytest

pytest.importorskip("pybloom_live")

import random

from neoswga.core import filter as filter_mod
from neoswga.core import parameter
from neoswga.core.background_filter import BackgroundBloomFilter, SampledGenomeIndex

# A primer that will be planted heavily in the background, and one that will not
# appear in it at all.
COMMON = "ACGTACGTAC"
ABSENT = "TTTTTTTTTT"


@pytest.fixture
def background(tmp_path):
    """A background FASTA with COMMON planted many times, plus the artefacts
    `build-filter` writes: a Bloom filter and a sampled index beside it."""
    rng = random.Random(4242)
    seq = [rng.choice("ACGT") for _ in range(60_000)]
    for i in range(300):
        pos = i * 200
        seq[pos : pos + len(COMMON)] = list(COMMON)
    fasta = tmp_path / "bg.fasta"
    fasta.write_text(">bg\n" + "".join(seq) + "\n")

    bloom = BackgroundBloomFilter(capacity=2_000_000, error_rate=0.01)
    bloom.add_genome(str(fasta), min_k=10, max_k=10)
    bloom_path = tmp_path / "bg_bloom.pkl"
    bloom.save(str(bloom_path))

    # sample_rate=1 so the test's expectations are exact rather than
    # statistical; the extrapolation itself is covered elsewhere.
    sampled = SampledGenomeIndex(sample_rate=1)
    sampled.add_genome(str(fasta), min_k=10, max_k=10)
    sampled.save(str(tmp_path / "bg_sampled.pkl"))

    return {"dir": tmp_path, "bloom": str(bloom_path), "length": 60_000}


@pytest.fixture(autouse=True)
def clean_parameter(background, monkeypatch):
    monkeypatch.setattr(parameter, "bloom_filter_path", background["bloom"], raising=False)
    monkeypatch.setattr(parameter, "sampled_index_path", None, raising=False)
    monkeypatch.setattr(parameter, "bloom_max_bg_matches", 10, raising=False)


# ----------------------------------------------------------------------
# The scale inversion
# ----------------------------------------------------------------------


def test_a_heavily_present_primer_is_not_reported_as_rare(background):
    """The core regression.

    COMMON occurs 300 times in 60 kb. Whatever the path reports has to reflect
    that, not a fixed 11 that shrinks to nothing once divided by a real genome
    length.
    """
    counts = filter_mod.get_bg_rates_via_bloom([COMMON], background["bloom"])

    assert counts[COMMON] > 100, (
        f"reported {counts[COMMON]} background sites for a primer planted 300 "
        f"times; a sentinel cannot represent frequency"
    )


def test_the_scale_inversion_is_real_arithmetic(background):
    """States the failure in the terms the gate uses.

    The same sentinel is restrictive on a small background and permissive on a
    large one, because it is an absolute count fed to a frequency test. That
    inversion -- stricter where filtering is easy, useless where it is hard --
    is what made the bug invisible: every small-scale test of this path passed.
    """
    sentinel = 10 + 1  # bloom_max_bg_matches + 1, as it used to be emitted
    max_bg_freq = 5e-6

    assert sentinel / 6_258 > max_bg_freq, "sentinel rejected on a 6 kb plasmid"
    assert sentinel / 3.0e9 < max_bg_freq, "sentinel passed on a 3 Gbp human background"

    # What the fixed path reports instead scales with the genome it measured.
    counts = filter_mod.get_bg_rates_via_bloom([COMMON], background["bloom"])
    assert counts[COMMON] / background["length"] > max_bg_freq


def test_an_absent_primer_still_reports_zero(background):
    """Bloom filters have no false negatives, so absence is trustworthy and
    zero is the honest answer -- unlike presence, which carries no count."""
    counts = filter_mod.get_bg_rates_via_bloom([ABSENT], background["bloom"])

    assert counts[ABSENT] == 0


# ----------------------------------------------------------------------
# Finding the sampled index, and refusing to guess without it
# ----------------------------------------------------------------------


def test_the_sampled_index_is_found_next_to_the_bloom_filter(background):
    """`build-filter` writes bg_bloom.pkl and bg_sampled.pkl into the same
    directory, but using the second required setting `sampled_index_path` by
    hand -- undocumented, and easy to miss. Missing it was what put the path on
    the sentinel branch, so the common case has to work without it.
    """
    counts = filter_mod.get_bg_rates_via_bloom([COMMON], background["bloom"])

    assert counts[COMMON] > 100


def test_an_explicit_sampled_index_still_wins(background, monkeypatch):
    explicit = str(background["dir"] / "bg_sampled.pkl")
    monkeypatch.setattr(parameter, "sampled_index_path", explicit, raising=False)

    counts = filter_mod.get_bg_rates_via_bloom([COMMON], background["bloom"])
    assert counts[COMMON] > 100


def test_without_any_sampled_index_it_refuses_rather_than_guessing(background, monkeypatch):
    """The heart of the fix.

    A Bloom filter cannot answer "how often", and manufacturing a number to
    stand in for one is what silently disabled the filter. Failing here costs a
    user one clear error; the alternative cost them a primer set screened
    against nothing, with no indication anything was wrong.
    """
    (background["dir"] / "bg_sampled.pkl").unlink()
    monkeypatch.setattr(parameter, "sampled_index_path", None, raising=False)

    with pytest.raises(ValueError, match="sampled index"):
        filter_mod.get_bg_rates_via_bloom([COMMON], background["bloom"])


def test_the_refusal_says_how_to_fix_it(background, monkeypatch):
    (background["dir"] / "bg_sampled.pkl").unlink()
    monkeypatch.setattr(parameter, "sampled_index_path", None, raising=False)

    with pytest.raises(ValueError) as excinfo:
        filter_mod.get_bg_rates_via_bloom([COMMON], background["bloom"])

    message = str(excinfo.value)
    assert "build-filter" in message
    assert "sampled_index_path" in message


# ----------------------------------------------------------------------
# The gate, end to end
# ----------------------------------------------------------------------


def test_the_gate_now_discriminates(background):
    """What the filter is for: a primer saturating the background fails the
    frequency gate while one absent from it passes.

    Both passed before, on any background large enough to matter.
    """
    counts = filter_mod.get_bg_rates_via_bloom([COMMON, ABSENT], background["bloom"])
    length = background["length"]
    max_bg_freq = 5e-6

    assert counts[COMMON] / length >= max_bg_freq, "saturating primer should fail the gate"
    assert counts[ABSENT] / length < max_bg_freq, "absent primer should pass the gate"


# ----------------------------------------------------------------------
# Saying so when the sample is too sparse to resolve the threshold
# ----------------------------------------------------------------------


def test_it_warns_when_the_sample_cannot_resolve_the_threshold(background, monkeypatch, caplog):
    """The sampled index extrapolates by `sample_rate`, so it can only resolve
    counts that survive being divided by it.

    The count at the gate's threshold is `max_bg_freq * genome_size`. On a
    3 Gbp human background at 5e-6 that is 15000 sites, sampled ~150 times at
    rate 100 -- resolvable. On a 6 kb plasmid it is 0.03 sites, sampled
    0.0003 times -- so every primer estimates zero and passes, which is the
    same silent under-filtering the sentinel caused, arriving by a different
    route.

    Reported rather than raised: unlike a missing index, the numbers here are
    real, just too coarse, and the caller may know their background better than
    this heuristic does.
    """
    import logging

    monkeypatch.setattr(parameter, "bg_seq_lengths", [background["length"]], raising=False)
    monkeypatch.setattr(parameter, "max_bg_freq", 5e-6, raising=False)

    with caplog.at_level(logging.WARNING):
        filter_mod.get_bg_rates_via_bloom([COMMON], background["bloom"])

    assert "sample_rate" in caplog.text.lower() or "sparse" in caplog.text.lower()


def test_no_warning_when_the_sample_can_resolve_the_threshold(background, monkeypatch, caplog):
    """A human-scale background at the same rate resolves the threshold
    comfortably, and should pass without noise."""
    import logging

    monkeypatch.setattr(parameter, "bg_seq_lengths", [3_000_000_000], raising=False)
    monkeypatch.setattr(parameter, "max_bg_freq", 5e-6, raising=False)

    with caplog.at_level(logging.WARNING):
        filter_mod.get_bg_rates_via_bloom([COMMON], background["bloom"])

    assert "too sparse" not in caplog.text.lower()
