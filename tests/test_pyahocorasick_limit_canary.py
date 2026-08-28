"""Keeps the 2**31 workaround honest as pyahocorasick moves.

`string_search.MAX_SCAN_CHUNK` exists because `Automaton.iter()` indexes with
a 32-bit int: handed a string longer than `2**31 - 1` its scan loop never runs
and it yields nothing, with no exception and no warning. Human is 3.1 Gb and
mouse 2.7 Gb, so without the chunking every mammalian background scores zero
binding sites -- indistinguishable, downstream, from a perfect primer set.

Verified against 2.3.1, the current release, and absent from the upstream
tracker at the time. That leaves a standing question no comment can answer: is
it still true? A note saying "re-check on upgrade" only works if somebody reads
it on the day they upgrade.

So the check is automated in two parts.

`test_the_limit_check_is_still_current` is cheap and always runs. It fails the
moment the installed version leaves the set we have actually tested, which is
the moment the question needs asking again. It follows the convention this repo
already uses for black and isort in CI, and for the model checksum allowlist:
pin, and bump deliberately.

`test_the_upstream_limit_is_still_present` is the answer to that question, and
it is the expensive one -- it allocates over two gigabytes, so it runs only
when explicitly asked. That is the command the failure message names.

If upstream ever fixes this, the second test fails. At that point the chunking
is still correct, merely no longer load-bearing, and the decision to keep or
drop it belongs to a human.
"""

import os

import pytest

pytest.importorskip("ahocorasick")

import ahocorasick  # noqa: E402

from neoswga.core import string_search  # noqa: E402

# Versions whose behaviour past 2**31 we have measured, rather than assumed.
VERIFIED_BROKEN = {"2.3.1"}

REVERIFY = "NEOSWGA_VERIFY_AHOCORASICK_LIMIT"

_NEEDLE = "ACGTACGTACGT"


def _installed_version():
    from importlib.metadata import PackageNotFoundError, version

    try:
        return version("pyahocorasick")
    except PackageNotFoundError:  # pragma: no cover - installed by definition
        pytest.skip("pyahocorasick metadata unavailable")


def test_the_limit_check_is_still_current():
    """Fail when the installed version is one we have not measured.

    Not a claim that newer versions are broken -- a claim that we do not know,
    and that MAX_SCAN_CHUNK is currently justified by a measurement that no
    longer applies to what is installed.
    """
    installed = _installed_version()

    assert installed in VERIFIED_BROKEN, (
        f"pyahocorasick {installed} has not been checked against the 2**31 "
        f"scan limit that string_search.MAX_SCAN_CHUNK works around "
        f"(measured on {sorted(VERIFIED_BROKEN)}).\n"
        f"Re-measure with:\n"
        f"    {REVERIFY}=1 pytest "
        f"tests/test_pyahocorasick_limit_canary.py -k still_present\n"
        f"If the limit is still there, add {installed!r} to VERIFIED_BROKEN. "
        f"If it is gone, the chunking is correct but no longer load-bearing, "
        f"and keeping or dropping it is a deliberate choice -- see Known "
        f"Issues #5 in CLAUDE.md."
    )


def test_the_chunk_stays_under_the_limit():
    """Cheap invariant, independent of what upstream does."""
    assert string_search.MAX_SCAN_CHUNK <= 2**31 - 1


@pytest.mark.slow
@pytest.mark.scale
@pytest.mark.skipif(
    os.environ.get(REVERIFY) != "1",
    reason=f"allocates >2 GB; set {REVERIFY}=1 to run",
)
def test_the_upstream_limit_is_still_present():
    """The measurement the workaround rests on.

    The needle sits at offset 1000, far below any boundary, so a match found
    just under the limit and missed just over it isolates the scan being
    skipped entirely rather than truncated.
    """
    automaton = ahocorasick.Automaton()
    automaton.add_word(_NEEDLE, _NEEDLE)
    automaton.make_automaton()

    def matches(total_len):
        haystack = "A" * 1000 + _NEEDLE + "A" * (total_len - 1000 - len(_NEEDLE))
        return len(list(automaton.iter(haystack)))

    assert matches(2**31 - 10) == 1, "sanity: the needle is findable below the limit"
    assert matches(2**31 + 10) == 0, (
        "pyahocorasick now finds matches past 2**31. The limit "
        "string_search.MAX_SCAN_CHUNK works around appears to be fixed "
        "upstream; revisit Known Issues #5."
    )
