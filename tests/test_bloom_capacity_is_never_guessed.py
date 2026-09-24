"""No caller may take a default Bloom capacity.

`BackgroundBloomFilter.__init__` declared `capacity: int = 3e9` -- a float on
an int annotation, and a default nobody should take, since pybloom allocates
from it upfront: 3e9 items at 9.59 bits is about 3.6 GB.

`genome_library.add_genome` took it. That path builds a filter automatically
for any reference above 50 Mb, and it was wrong in both directions at once:
3.6 GB is far more than the 26.8 MB a k 6-12 filter needs, and far LESS than
the 14.6 billion distinct k-mers a k 6-18 filter holds, which is the range
that path actually computes. pybloom raises IndexError at capacity, and the
surrounding `except Exception` turned that into "Bloom filter build failed"
with `bloom_path = None`, so the reference was registered with no filter and
the run continued.

Making capacity required is what stops a third call site from appearing with
the same defect.
"""

import inspect

import pytest

pytest.importorskip("pybloom_live")

from neoswga.core.background_filter import BackgroundBloomFilter


def test_capacity_has_no_default():
    signature = inspect.signature(BackgroundBloomFilter.__init__)
    assert signature.parameters["capacity"].default is inspect.Parameter.empty, (
        "a default capacity is an allocation nobody chose"
    )


def test_the_genome_library_sizes_its_filter_from_the_range_it_computed():
    """It passed neither a capacity nor the k range it had just counted."""
    source = inspect.getsource(
        __import__("neoswga.core.genome_library", fromlist=["x"])
    )
    assert "distinct_kmer_capacity" in source, (
        "the library builds a Bloom filter without sizing it to the k-mers it holds"
    )
    assert "BackgroundBloomFilter()" not in source


def test_a_capacity_that_is_too_small_is_reported_not_swallowed():
    """pybloom raises IndexError at capacity; that must stay visible."""
    bloom = BackgroundBloomFilter(capacity=3, error_rate=0.01)
    with pytest.raises(IndexError):
        for kmer in ("AAAAAA", "CCCCCC", "GGGGGG", "TTTTTT", "ACGTAC", "TGCATG"):
            bloom.add(kmer)
