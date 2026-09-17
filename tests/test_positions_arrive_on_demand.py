"""Positions for a frontier that moves, and a provider that cannot read a zero.

Phase 4 increment 3 of the plan for `docs/validation/pipeline_audit_2026-09-16/`.

Two things are pinned here.

`PositionCache` was built once over a fixed candidate list. That was sufficient
while the search never looked past the `max_primer` shortlist, and it stops
being sufficient the moment the frontier advances: a candidate admitted later
has no entry, `get_positions` answers with an empty array, and the coverage
computed for it is a zero rather than a measurement. `load` brings a batch in
and `release` lets the one behind it go, so a moving window does not accumulate
the whole 443 MB index.

`CandidateProvider.ensure_positions` was supposed to catch exactly that, and
could not. Its predicate asked whether the primer had a hit on ANY prefix, so a
candidate with fifty foreground sites and no background entry at all passed,
which is the silent-zero shape of Known Issues 5, 6 and 13 reached once more.
The question is whether there is an ENTRY on EVERY required prefix, and a
genuine zero on one of them is an answer rather than a failure.
"""

import h5py
import numpy as np
import pytest

from neoswga.core.candidate_provider import CandidateProvider
from neoswga.core.position_cache import MissingPositionsError, PositionCache

FOREGROUND = {
    "AAAACCCCGGGG": [100, 400, 900],
    "TTTTCCCCGGGG": [250, 650],
    # Indexed on the foreground, and deliberately absent from the background
    # below: the candidate whose specificity is unknown, not zero.
    "GGGGCCCCAAAA": [500],
}

BACKGROUND = {
    "AAAACCCCGGGG": [7000],
    # Present in the index with no sites. A scanned zero is a measurement and
    # must not be mistaken for an absent entry.
    "TTTTCCCCGGGG": [],
}


def _write_index(path, entries, record_starts=(0,)):
    with h5py.File(path, "w") as handle:
        for primer, positions in entries.items():
            handle.create_dataset(primer, data=np.array(positions, dtype=np.int64))
        handle.create_dataset("#record_starts", data=np.array(list(record_starts), dtype=np.int64))


@pytest.fixture
def prefixes(tmp_path):
    """One foreground and one background prefix, both record-aware."""
    _write_index(tmp_path / "fg_12mer_positions.h5", FOREGROUND)
    _write_index(tmp_path / "bg_12mer_positions.h5", BACKGROUND)
    return str(tmp_path / "fg"), str(tmp_path / "bg")


# -- the cache ------------------------------------------------------------


def test_load_reaches_a_candidate_the_cache_was_not_built_over(prefixes):
    """The frontier advanced; the positions have to follow it."""
    fg, bg = prefixes
    cache = PositionCache([fg, bg], ["AAAACCCCGGGG"], on_missing="warn")

    assert cache.get_positions(fg, "TTTTCCCCGGGG", "both").size == 0

    cache.load(["TTTTCCCCGGGG"])

    assert sorted(cache.get_positions(fg, "TTTTCCCCGGGG", "both")) == [250, 650]


def test_load_honours_the_missing_policy_rather_than_bypassing_it(prefixes):
    """A candidate absent from the index is refused, not admitted as a zero."""
    fg, bg = prefixes
    cache = PositionCache([fg, bg], ["AAAACCCCGGGG"], on_missing="error")

    with pytest.raises(MissingPositionsError):
        cache.load(["CACACACACACA"])


def test_load_reports_which_candidates_it_brought_in(prefixes):
    """Loading an already-held candidate is not work, and says so."""
    fg, bg = prefixes
    cache = PositionCache([fg, bg], ["AAAACCCCGGGG"], on_missing="warn")

    assert cache.load(["AAAACCCCGGGG"]) == []
    assert cache.load(["TTTTCCCCGGGG"]) == ["TTTTCCCCGGGG"]


def test_release_frees_the_candidate_the_search_moved_past(prefixes):
    """A window that only grows is not a window."""
    fg, bg = prefixes
    cache = PositionCache([fg, bg], ["AAAACCCCGGGG", "TTTTCCCCGGGG"], on_missing="warn")
    # Materialise the memoized 'both' key as a real search would.
    cache.get_positions(fg, "TTTTCCCCGGGG", "both")
    held = len(cache.cache)

    cache.release(["TTTTCCCCGGGG"])

    assert len(cache.cache) < held
    assert not any(key[1] == "TTTTCCCCGGGG" for key in cache.cache)


def test_reading_a_released_candidate_raises_rather_than_reading_zero(prefixes):
    """Releasing must not open the hole the release was meant to avoid."""
    fg, bg = prefixes
    cache = PositionCache([fg, bg], ["AAAACCCCGGGG", "TTTTCCCCGGGG"], on_missing="warn")
    cache.release(["TTTTCCCCGGGG"])

    with pytest.raises(MissingPositionsError):
        cache.get_positions(fg, "TTTTCCCCGGGG", "both")


def test_release_leaves_the_rest_of_the_cache_alone(prefixes):
    """Including the record geometry, which is per prefix and not per primer."""
    fg, bg = prefixes
    cache = PositionCache([fg, bg], ["AAAACCCCGGGG", "TTTTCCCCGGGG"], on_missing="warn")

    cache.release(["TTTTCCCCGGGG"])

    assert sorted(cache.get_positions(fg, "AAAACCCCGGGG", "both")) == [100, 400, 900]
    assert cache.get_record_starts(fg) == [0]


def test_a_released_candidate_can_be_loaded_again(prefixes):
    """A frontier can revisit; the release is not a verdict."""
    fg, bg = prefixes
    cache = PositionCache([fg, bg], ["AAAACCCCGGGG", "TTTTCCCCGGGG"], on_missing="warn")
    cache.release(["TTTTCCCCGGGG"])

    cache.load(["TTTTCCCCGGGG"])

    assert sorted(cache.get_positions(fg, "TTTTCCCCGGGG", "both")) == [250, 650]


# -- the provider ---------------------------------------------------------


class _StubInventory:
    """Just enough inventory to build a provider; ordering is tested elsewhere."""

    def __init__(self, sequences):
        self._sequences = list(sequences)

    def current_policy(self, condition_id):
        return None

    def iter_eligible(self, condition_id, lengths, *rest):
        return list(self._sequences)


def _provider(sequences=("AAAACCCCGGGG",)):
    return CandidateProvider(_StubInventory(sequences), "cond", [12])


def test_a_provider_without_a_cache_refuses_to_vouch_for_a_batch():
    """The configuration that produced the silent zero now raises."""
    provider = _provider()

    with pytest.raises(MissingPositionsError):
        provider.ensure_positions(["AAAACCCCGGGG"])


def test_a_candidate_indexed_on_only_one_prefix_is_refused(prefixes):
    """Foreground sites and no background entry is unknown specificity."""
    fg, bg = prefixes
    cache = PositionCache([fg, bg], ["GGGGCCCCAAAA"], on_missing="warn")
    provider = _provider()
    provider.attach_positions(cache)

    with pytest.raises(MissingPositionsError) as raised:
        provider.ensure_positions(["GGGGCCCCAAAA"])

    assert bg in str(raised.value)


def test_a_measured_zero_on_one_prefix_is_accepted(prefixes):
    """An indexed candidate that binds the host nowhere is a good candidate."""
    fg, bg = prefixes
    cache = PositionCache([fg, bg], ["TTTTCCCCGGGG"], on_missing="warn")
    provider = _provider()
    provider.attach_positions(cache)

    provider.ensure_positions(["TTTTCCCCGGGG"])


def test_a_candidate_outside_the_cache_is_loaded_rather_than_refused(prefixes):
    """The batch the frontier just admitted is the normal case, not an error."""
    fg, bg = prefixes
    cache = PositionCache([fg, bg], ["GGGGCCCCAAAA"], on_missing="warn")
    provider = _provider()
    provider.attach_positions(cache)

    provider.ensure_positions(["AAAACCCCGGGG"])

    assert sorted(cache.get_positions(bg, "AAAACCCCGGGG", "both")) == [7000]
