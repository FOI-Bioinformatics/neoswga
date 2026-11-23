"""`filter` and `optimize` read one default k range from one params.json.

`get_params` defaults `min_k`/`max_k` from the polymerase (equiphi29 10-18,
bst 15-25, klenow 8-15); `resolve_design_request` hard-coded 6-12. A file that
set `polymerase: equiphi29` and no k range was filtered as 10-18-mers and then
refused by `optimize` for oligo lengths 6-9 (2026-11-23 sweep, tier 1 item 6).
"""

import pytest

from neoswga.core.design_request import resolve_design_request
from neoswga.core.registry import views

TEMPS = {"phi29": 30.0, "equiphi29": 42.0, "bst": 63.0, "klenow": 37.0}


def _params(polymerase):
    return {
        "fg_prefixes": ["target"],
        "fg_genomes": ["target.fna"],
        "polymerase": polymerase,
        "reaction_temp": TEMPS[polymerase],
        "num_primers": 6,
    }


@pytest.mark.parametrize("polymerase", sorted(TEMPS))
def test_request_k_range_matches_get_params(polymerase):
    low, high = views.primer_length_ranges()[polymerase]
    request = resolve_design_request(_params(polymerase))
    assert request.primer_lengths == tuple(range(low, high + 1))


def test_an_explicit_k_range_still_wins():
    params = _params("phi29")
    params.update(min_k=8, max_k=10)
    assert resolve_design_request(params).primer_lengths == (8, 9, 10)
