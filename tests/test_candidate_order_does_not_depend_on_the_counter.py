"""Which k-mer counter ran must not change the candidates or the panel.

Jellyfish dumps k-mers in hash order and KMC in sorted order. Step 2 sorted by
`ratio` then `fg_count`, and rows tying on both kept their INPUT order -- the
order the counter emitted them in. So the counter decided the order of tied
candidates, which leads step 3 into an order-sensitive optimizer, and through
`[:max_primer]` it decided which tied candidates survived the shortlist.

Measured on the plasmid example before the fix: identical values for every
primer, different row order in step 2 and step 3. After it, step 2, step 3 and
the delivered panel are byte-identical between the two counters, and the KMC
run writes no text table at all.

This was found by running the whole pipeline under each counter and comparing
the files, after the claim that the fallback "changes speed, never a result"
had already been written down. The claim was wrong until this change.
"""

import pandas as pd

from neoswga.core import pipeline


def _tied_frame(order):
    rows = {
        "GATGTACTGC": (0.0, 3),
        "TGGCAGTACA": (0.0, 3),
        "TGTACTGCCAA": (0.0, 3),
        "ATGTACTGCCA": (0.0, 3),
        "CCCATTGAC": (0.0, 5),
        "GGCAACTA": (0.4, 3),
    }
    return pd.DataFrame(
        [{"primer": p, "ratio": rows[p][0], "fg_count": rows[p][1]} for p in order]
    )


HASH_ORDER = ["GATGTACTGC", "TGGCAGTACA", "TGTACTGCCAA", "ATGTACTGCCA", "CCCATTGAC", "GGCAACTA"]
SORTED_ORDER = sorted(HASH_ORDER)


def test_tied_candidates_sort_the_same_whatever_order_they_arrive_in():
    a = _tied_frame(HASH_ORDER).sort_values(**pipeline._RANK_BY_RATIO)
    b = _tied_frame(SORTED_ORDER).sort_values(**pipeline._RANK_BY_RATIO)
    assert list(a.primer) == list(b.primer)


def test_the_shortlist_cut_keeps_the_same_candidates_whatever_the_input_order():
    """The costlier half: a cut inside a run of ties chose different
    candidates depending on emission order, not just a different order."""
    cut = 3
    a = _tied_frame(HASH_ORDER).sort_values(**pipeline._RANK_BY_RATIO)[:cut]
    b = _tied_frame(SORTED_ORDER).sort_values(**pipeline._RANK_BY_RATIO)[:cut]
    assert set(a.primer) == set(b.primer)


def test_the_primary_keys_still_lead():
    """The tie-break decides only ties: a better ratio and a higher count still
    come first."""
    ranked = list(_tied_frame(HASH_ORDER).sort_values(**pipeline._RANK_BY_RATIO).primer)
    assert ranked[0] == "CCCATTGAC", "the higher fg_count at ratio 0 must lead"
    assert ranked[-1] == "GGCAACTA", "the worse ratio must come last"


def test_the_occupancy_ranking_breaks_ties_the_same_way():
    frame = pd.DataFrame(
        [{"primer": p, "occupancy_ratio": 0.1, "fg_count": 3} for p in HASH_ORDER]
    )
    reordered = frame.iloc[::-1]
    a = frame.sort_values(**pipeline._RANK_BY_OCCUPANCY)
    b = reordered.sort_values(**pipeline._RANK_BY_OCCUPANCY)
    assert list(a.primer) == list(b.primer)
