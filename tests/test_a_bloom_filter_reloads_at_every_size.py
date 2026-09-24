"""A saved filter must reload whatever geometry it was built with.

`pybloom_live.make_hashfuncs` selects the hash constructor from the filter's
GEOMETRY -- `num_slices` and `num_bits` -- so which constructor ends up in the
pickle moves with capacity and error rate. The safe-pickle allowlist named only
`_hashlib.openssl_sha256`, which is what one observed filter happened to carry.

At the shipped error rate of 0.01, seven slices:

    capacity below about 3,400        num_bits < 2**15     xxhash.xxh3_128
    up to about 224 million           num_bits < 2**31     sha256
    above that                        num_bits >= 2**31    sha512

Other error rates change the slice count and reach sha384 and sha1. So a
filter over a small background could not be reloaded, and neither could one
for a long-oligo design against a host genome, where the distinct k-mer count
exceeds 224 million. `save()` succeeded and `load()` raised, which is the
asymmetry the bitarray entries in that allowlist were already added for.

The ratchet below asserts the RULE rather than one instance: every constructor
`make_hashfuncs` can select must be covered. A behavioural test could only
reach the sha512 arm by allocating 268 MB.
"""

import hashlib

import pytest

pytest.importorskip("pybloom_live")

from neoswga.core.background_filter import BackgroundBloomFilter
from neoswga.core.safe_pickle import _ALLOWED_CLASSES


def _selectable_hash_constructors():
    """Every hashfn `pybloom_live.make_hashfuncs` can return."""
    import xxhash

    return [
        hashlib.sha512,
        hashlib.sha384,
        hashlib.sha256,
        hashlib.sha1,
        xxhash.xxh128,
    ]


@pytest.mark.parametrize("hashfn", _selectable_hash_constructors())
def test_every_hash_pybloom_can_pick_is_allowlisted(hashfn):
    allowed = _ALLOWED_CLASSES["bloom_filter"]
    key = (hashfn.__module__, hashfn.__qualname__)
    assert key in allowed, (
        f"pybloom can build a filter hashed with {key[0]}.{key[1]}, and such a "
        f"filter would save successfully and refuse to load"
    )


@pytest.mark.parametrize("capacity", [1_000, 50_000])
def test_a_filter_round_trips_on_both_sides_of_the_hash_boundary(tmp_path, capacity):
    """1,000 selects xxh3_128 and 50,000 selects sha256, at error rate 0.01."""
    bloom = BackgroundBloomFilter(capacity=capacity, error_rate=0.01)
    bloom.add("ACGTACGTAC")

    path = tmp_path / f"bloom_{capacity}.pkl"
    bloom.save(str(path))
    reloaded = BackgroundBloomFilter.load(str(path))

    assert reloaded.contains("ACGTACGTAC"), (
        f"a filter at capacity {capacity} did not survive a save/load round trip"
    )
