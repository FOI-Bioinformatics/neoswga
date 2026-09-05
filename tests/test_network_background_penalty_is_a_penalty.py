"""Host binding must lower a candidate's score, never raise it.

The network greedy scored a candidate as

    base_score = fg_improvement / (bg_penalty + 1.0)

where `bg_penalty` is the CHANGE IN THE AVERAGE component size of the background
network. Adding a primer whose host sites are scattered and isolated pulls that
average DOWN, so `bg_penalty` goes negative, and once it passes -1 the
denominator is negative too. Measured with a stub cache, one primer contributing
scattered background sites:

    isolated bg sites added   bg_penalty   base_score
                          0        +0.00       +2.000
                          1        -2.50       -1.333
                         11        -4.58       -0.558

Three things follow, all of them about which primers get selected rather than
about what gets reported.

Past the sign change the score RISES with host load, so the background penalty
has become a background reward. A `bg_penalty` of exactly -1.0 is reachable and
yields `+inf`, which `optimize_greedy` accepts immediately against its `-inf`
starting best and which no later candidate can displace. And because the Tm and
dimer modifiers are applied multiplicatively, a negative base score inverts
them: at `dimer_penalty=0.5` a fully dimerising candidate scored -0.258 against
-0.516 for a clean one, so the dimer penalty became a dimer bonus.

This is a sign error, and distinct from the documented limitation that a
multiplicative soft penalty can never cut more than its own weight.
"""

import pytest

pytest.importorskip("networkx")

import numpy as np

from neoswga.core.network_optimizer import AmplificationNetwork, NetworkOptimizer

GENOME_LENGTH = 1_000_000


class _Cache:
    """Positions by (prefix, primer, strand); everything else binds nowhere.

    Both strands matter here. An amplification edge joins a forward site to a
    reverse site within reach, so a stub that reports every site as forward
    produces a graph with no edges at all, every component of size one, and an
    average component size pinned at 1.0 that no addition can move -- which
    would make this whole file pass against the unfixed code.
    """

    def __init__(self, positions):
        self._positions = positions

    def get_positions(self, prefix, primer, strand="both"):
        return np.array(self._positions.get((prefix, primer, strand), []), dtype=np.int64)


def _optimizer(positions):
    return NetworkOptimizer(
        position_cache=_Cache(positions),
        fg_prefixes=["fg"],
        bg_prefixes=["bg"],
        fg_seq_lengths=[GENOME_LENGTH],
        bg_seq_lengths=[GENOME_LENGTH],
        max_extension=10_000,
    )


def _score_with_background_sites(n_isolated):
    """Score one candidate that always helps the target and binds the host n times.

    The host sites are spread far past `max_extension` so each is its own
    component, which is what drags the average component size down.
    """
    positions = {
        # The candidate helps the target the same way every time: a forward and
        # a reverse site close enough to form one amplicon.
        ("fg", "CAND", "forward"): [10_000],
        ("fg", "CAND", "reverse"): [15_000],
        # Its host sites are all forward and 100 kb apart, well past
        # max_extension, so each is an isolated component of one and every
        # additional one drags the average component size further down.
        ("bg", "CAND", "forward"): [500_000 + i * 100_000 for i in range(n_isolated)],
        ("bg", "CAND", "reverse"): [],
        # The background network starts holding one real two-site component, so
        # the average it starts from is above 1 and can be pulled down.
        ("bg", "SEED", "forward"): [1_000],
        ("bg", "SEED", "reverse"): [4_000],
    }

    optimizer = _optimizer(positions)

    fg_network = AmplificationNetwork(max_extension=10_000)
    bg_network = AmplificationNetwork(max_extension=10_000)
    bg_network.add_primer_sites("SEED", np.array([1_000], dtype=np.int64), "+")
    bg_network.add_primer_sites("SEED", np.array([4_000], dtype=np.int64), "-")
    bg_network.build_edges()

    return optimizer._evaluate_primer_addition("CAND", [], fg_network, bg_network)


# ---------------------------------------------------------------------------
# The property
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_isolated", [0, 1, 2, 5, 11, 20])
def test_the_score_is_finite_however_much_host_binding_a_candidate_has(n_isolated):
    """`+inf` is selectable and undisplaceable, so it must not be reachable."""
    score = _score_with_background_sites(n_isolated)
    assert np.isfinite(score), (
        f"a candidate with {n_isolated} isolated host sites scored {score}; "
        "greedy compares with a strict > from -inf, so an infinite score is "
        "taken immediately and can never be replaced"
    )


def test_more_host_binding_never_scores_better():
    """The ordering the term exists to impose.

    Not a strict decrease: the background network's structure can leave two
    counts genuinely equivalent. What must never happen is host load buying a
    candidate a HIGHER score.
    """
    counts = [0, 1, 2, 5, 11, 20]
    scores = [_score_with_background_sites(n) for n in counts]

    for (lo, lo_score), (hi, hi_score) in zip(zip(counts, scores), list(zip(counts, scores))[1:]):
        assert hi_score <= lo_score + 1e-9, (
            f"{hi} host sites scored {hi_score:.4f}, better than {lo} host sites "
            f"at {lo_score:.4f}. The background penalty is acting as a reward."
        )


def test_a_clean_candidate_beats_an_identical_one_that_binds_the_host():
    """The decision the greedy actually makes, stated directly."""
    clean = _score_with_background_sites(0)
    dirty = _score_with_background_sites(11)
    assert clean > dirty, (
        f"a candidate binding the host 11 times scored {dirty:.4f} against "
        f"{clean:.4f} for one binding it not at all"
    )
