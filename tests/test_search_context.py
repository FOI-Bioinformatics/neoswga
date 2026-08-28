"""Tests for the typed search-context containers (pure data, no I/O)."""

import numpy as np
import pytest

from neoswga.core.search_context import (
    BFSConfig,
    BFSSearchContext,
    DimerConstraints,
    EvaluationContext,
    GenomeInfo,
    PositionData,
    SearchResult,
    SearchState,
)

# --- GenomeInfo -------------------------------------------------------------


def test_genome_info_from_lists_and_properties():
    g = GenomeInfo.from_lists(["a", "b"], [100, 200], circular=False)
    assert g.prefixes == ("a", "b")
    assert g.seq_lengths == (100, 200)
    assert g.circular is False
    assert g.total_length == 300
    assert g.num_sequences == 2


def test_genome_info_length_mismatch_raises():
    with pytest.raises(ValueError):
        GenomeInfo(prefixes=("a",), seq_lengths=(1, 2))


# --- PositionData -----------------------------------------------------------


def test_position_data_empty_and_add_and_total():
    pd_ = PositionData.empty(2)
    assert pd_.total_sites() == 0
    pd_.add_positions(0, np.array([1, 2, 3]), np.array([4, 5]))
    pd_.add_positions(1, np.array([7]), np.array([]))
    assert pd_.total_sites() == 3 + 2 + 1 + 0


def test_position_data_from_numpy_and_copy_is_independent():
    pd_ = PositionData.from_numpy([(np.array([1, 2]), np.array([3]))])
    clone = pd_.copy()
    clone.add_positions(0, np.array([99]), np.array([]))
    assert pd_.total_sites() == 3  # original unchanged
    assert clone.total_sites() == 4


def test_position_data_add_out_of_range_raises():
    pd_ = PositionData.empty(1)
    with pytest.raises(IndexError):
        pd_.add_positions(5, np.array([1]), np.array([2]))


# --- DimerConstraints -------------------------------------------------------


def test_dimer_constraints_compatibility_and_readonly_matrix():
    matrix = np.array([[False, True], [True, False]])
    dc = DimerConstraints(matrix=matrix, primer_to_index={"A": 0, "B": 1})
    # matrix made read-only
    assert dc.matrix.flags.writeable is False
    # A and B form a dimer (matrix[0,1] True)
    assert dc.are_compatible(["A"], "B") is False
    # Unknown primer is assumed compatible
    assert dc.are_compatible(["A"], "Z") is True
    # No conflict when partner not present
    assert dc.are_compatible([], "B") is True


# --- BFSConfig --------------------------------------------------------------


def test_bfs_config_validation():
    BFSConfig(max_sets=5, iterations=3, selection_method="softmax")
    with pytest.raises(ValueError):
        BFSConfig(max_sets=0)
    with pytest.raises(ValueError):
        BFSConfig(iterations=0)
    with pytest.raises(ValueError):
        BFSConfig(selection_method="bogus")


# --- SearchState ------------------------------------------------------------


def test_search_state_copy_add_primer_and_key():
    state = SearchState(
        primers=["A"],
        score=1.0,
        fg_positions=PositionData.empty(1),
        bg_positions=PositionData.empty(1),
    )
    new = state.add_primer(
        "B",
        fg_new_positions=[(np.array([10]), np.array([20]))],
        bg_new_positions=[(np.array([30]), np.array([]))],
        new_score=2.5,
    )
    assert state.primers == ["A"]  # original untouched
    assert new.primers == ["A", "B"] and new.score == 2.5
    assert new.fg_positions.total_sites() == 2
    assert new.compressed_key == "A,B"
    assert (
        SearchState(["B", "A"], 0.0, PositionData.empty(0), PositionData.empty(0)).compressed_key
        == "A,B"
    )


# --- BFSSearchContext -------------------------------------------------------


def _ctx(pool=("A", "A", "B", "C")):
    return BFSSearchContext(
        primer_pool=list(pool),
        fg_genome=GenomeInfo.from_lists(["fg"], [1000]),
        dimer_constraints=DimerConstraints(np.zeros((3, 3), dtype=bool), {"A": 0, "B": 1, "C": 2}),
        config=BFSConfig(max_sets=2),
    )


def test_context_dedupes_pool_and_empty_raises():
    ctx = _ctx()
    assert ctx.primer_pool == ["A", "B", "C"]  # deduped, order preserved
    assert ctx.has_background is False
    with pytest.raises(ValueError):
        BFSSearchContext(
            primer_pool=[],
            fg_genome=GenomeInfo.from_lists(["fg"], [1]),
            dimer_constraints=DimerConstraints(np.zeros((0, 0), dtype=bool), {}),
            config=BFSConfig(),
        )


def test_context_available_primers_excludes_banned_and_current():
    ctx = _ctx()
    ctx.banned_primers = {"C"}
    avail = ctx.get_available_primers(["A"])
    assert avail == ["B"]  # A in set, C banned


def test_context_score_cache_roundtrip():
    ctx = _ctx()
    assert ctx.get_cached_score(["B", "A"]) is None
    ctx.cache_score(["B", "A"], 3.14)
    assert ctx.get_cached_score(["A", "B"]) == 3.14  # order-independent key


def test_context_initialize_states_empty_and_from_sets():
    ctx = _ctx()
    ctx.initialize_states()
    assert len(ctx.top_states) == 2
    assert all(s.primers == [] for s in ctx.top_states)

    ctx.initialize_states([["A", "B"], ["C"], ["A"]])  # capped at max_sets=2
    assert len(ctx.top_states) == 2
    assert ctx.top_states[0].primers == ["A", "B"]


# --- SearchResult / EvaluationContext --------------------------------------


def test_search_result_best_and_to_dict():
    s1 = SearchState(["A"], 1.0, PositionData.empty(0), PositionData.empty(0))
    s2 = SearchState(["B", "C"], 5.0, PositionData.empty(0), PositionData.empty(0))
    res = SearchResult(states=[s1, s2], best_score=5.0, iterations_completed=3, converged=True)
    assert res.best_state is s2
    assert res.best_primers == ["B", "C"]
    d = res.to_dict()
    assert d["best_primers"] == ["B", "C"] and d["num_states"] == 2 and d["converged"]

    empty = SearchResult(states=[], best_score=0.0, iterations_completed=0)
    assert empty.best_state is None and empty.best_primers == []


def test_evaluation_context_has_background():
    fg = GenomeInfo.from_lists(["fg"], [10])
    assert EvaluationContext(fg_genome=fg).has_background() is False
    assert (
        EvaluationContext(
            fg_genome=fg, bg_genome=GenomeInfo.from_lists(["bg"], [20])
        ).has_background()
        is True
    )
