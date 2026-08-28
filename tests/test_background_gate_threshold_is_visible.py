"""`max_bg_freq` does not mean what a user writes, and nothing said so.

`_scale_freq_threshold` divides both frequency gates by `4^(k - 10)`. For the
foreground FLOOR that is necessary: a 12-mer binds a sixteenth as often as a
10-mer by construction, so applying an unscaled floor would reject every long
primer. Its docstring argues exactly that case, and it is right.

The same function is then applied to the background CEILING, where the argument
does not carry. The harm from a background site does not depend on primer
length -- a 12-mer with ten host sites primes the host ten times, exactly as a
10-mer with ten does -- so scaling the ceiling silently reinterprets the
parameter at every length:

    k     max_bg_freq 5e-7 becomes    sites allowed on chr21 (46.7 Mb)
    10          5.00e-07                       23.35
    12          3.13e-08                        1.46      -> bg_count <= 1
    14          1.95e-09                        0.09      -> bg_count == 0

Measured on the *Prevotella* 12-mer pool with the shipped equiphi29 config:

    pass the foreground floor                      626,984
    bg ceiling as configured (<23.35 sites)        573,182
    bg ceiling as applied    (< 1.46 sites)        143,224
    discarded by the scaling                       429,958   (75.0%)

The discarded primers are not worse candidates: median 2 foreground sites
either way, maximum 140 against 179. The scaling is not selecting for quality,
it is cutting the pool by a factor of four.

Whether the ceiling *should* scale is a design question with a real trade-off --
admitting primers with up to 23 host sites would raise coverage and lower
selectivity — and it is not settled here. What is settled is that the effect
must be visible. At k >= 13 the gate collapses to "zero exact background
matches" and no value of `max_bg_freq` in its documented range loosens it, so a
user tuning that parameter is tuning something that has stopped responding.
"""

import logging

import pytest

from neoswga.core.filter import _scale_freq_threshold, describe_freq_gate

BG_LEN = 46_709_983


def test_the_foreground_floor_scaling_is_unchanged():
    """The case the function was written for and gets right: an unscaled floor
    of 1e-5 demands 31.7 sites, which almost no 12-mer has."""
    assert _scale_freq_threshold(1e-5, 10) == pytest.approx(1e-5)
    assert _scale_freq_threshold(1e-5, 12) == pytest.approx(1e-5 / 16)


@pytest.mark.parametrize(
    "k,expected_sites",
    [(10, 23.35), (12, 1.46), (14, 0.09)],
)
def test_the_background_ceiling_in_sites_is_reported(k, expected_sites):
    """A frequency is not a quantity anyone can reason about; sites are."""
    gate = describe_freq_gate(5e-7, k, BG_LEN)

    assert gate["sites"] == pytest.approx(expected_sites, abs=0.02)
    assert gate["scaled"] == pytest.approx(_scale_freq_threshold(5e-7, k))


def test_a_ceiling_below_one_site_is_flagged():
    """Below one site the gate means "zero exact matches" and stops responding
    to the parameter."""
    assert describe_freq_gate(5e-7, 14, BG_LEN)["degenerate"] is True
    assert describe_freq_gate(5e-7, 12, BG_LEN)["degenerate"] is False


def test_the_gate_is_logged_in_sites_not_only_in_frequency(caplog):
    """The scaling was invisible: nothing printed the applied threshold, so a
    user comparing their params.json against the result had no way to see that
    5e-7 had become 3.13e-8."""
    from neoswga.core.filter import log_freq_gates

    with caplog.at_level(logging.INFO):
        log_freq_gates(
            primer_length=12,
            min_fg_freq=1e-5,
            max_bg_freq=5e-7,
            fg_length=3_168_282,
            bg_length=BG_LEN,
        )

    text = caplog.text
    assert "12" in text
    assert "site" in text.lower()
    assert "5e-07" in text or "5.00e-07" in text


def test_a_degenerate_ceiling_warns(caplog):
    from neoswga.core.filter import log_freq_gates

    with caplog.at_level(logging.WARNING):
        log_freq_gates(
            primer_length=14,
            min_fg_freq=1e-5,
            max_bg_freq=5e-7,
            fg_length=3_168_282,
            bg_length=BG_LEN,
        )

    warnings = [r for r in caplog.records if r.levelno >= logging.WARNING]
    assert warnings, "a gate that admits only zero-match primers must say so"
    assert "max_bg_freq" in caplog.text


def test_the_scaling_docstring_covers_the_ceiling_case():
    """It argued only the floor, which is why the ceiling went unexamined."""
    doc = _scale_freq_threshold.__doc__ or ""

    assert "ceiling" in doc.lower()
