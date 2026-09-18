"""Three of the five strand quantities were computed and thrown away.

`PositionCache.compute_strand_alternation_stats` returns five figures. The call
site in `base_optimizer._compute_metrics` read two of them, dropped
`strand_alternation_gap_mean`, `strand_alternation_gap_max` and
`longest_same_strand_run`, and `break`ed after the first foreground genome, so a
multi-target design reported one target's strand structure and the host's was
never computed at all.

`strand_alternation_gap_max` is the discarded one that matters most: exponential
amplification needs two sites in convergent orientation within the polymerase's
reach, so the widest gap between opposite-strand sites is the closest quantity
here to the mechanism. swga 2.0 approximates it with `within_mean_gap_ratio` and
computes it for both genomes. Item 3 of
`docs/validation/getting_ahead_on_spacing_2026-09-18.md`.

The other half of this is telling an uncomputed value from a measured one. The
old call site initialised both scalars to 0.0 and left them there when the cache
could not answer, so a zero meant either "measured zero" or "never asked", which
is the shape of Known Issues 5, 6 and 13.
"""

import pytest

from neoswga.core.strand_metrics import (
    STRAND_KEYS,
    collect_strand_stats,
    headline_strand_scalars,
)


class _Cache:
    """A cache that answers for the prefixes it was given and no others."""

    def __init__(self, answers, explode_on=()):
        self._answers = answers
        self._explode_on = set(explode_on)

    def compute_strand_alternation_stats(self, prefix, primers, genome_length):
        if prefix in self._explode_on:
            raise RuntimeError("no positions for this prefix")
        if prefix not in self._answers:
            raise KeyError(prefix)
        return dict(self._answers[prefix])


def _stats(score=0.7, ratio=0.9, gap_max=41_000.0, gap_mean=6_000.0, run=3):
    return {
        "strand_alternation_gap_mean": gap_mean,
        "strand_alternation_gap_max": gap_max,
        "strand_alternation_score": score,
        "strand_coverage_ratio": ratio,
        "longest_same_strand_run": run,
    }


PRIMERS = ["ACGTACGTAC", "TTGCATGCAT"]


class TestEveryQuantityIsKept:
    def test_all_five_keys_are_declared(self):
        assert set(STRAND_KEYS) == {
            "strand_alternation_gap_mean",
            "strand_alternation_gap_max",
            "strand_alternation_score",
            "strand_coverage_ratio",
            "longest_same_strand_run",
        }

    @pytest.mark.parametrize("key", sorted(STRAND_KEYS))
    def test_the_quantity_reaches_the_collected_stats(self, key):
        cache = _Cache({"target": _stats()})

        stats = collect_strand_stats(cache, ["target"], [3_200_000], PRIMERS)

        assert key in stats["target"]

    def test_the_widest_convergent_gap_is_kept(self):
        """The one the old call site dropped, and the closest of the five to the
        amplification mechanism."""
        cache = _Cache({"target": _stats(gap_max=41_000.0)})

        stats = collect_strand_stats(cache, ["target"], [3_200_000], PRIMERS)

        assert stats["target"]["strand_alternation_gap_max"] == 41_000.0


class TestEveryGenomeIsCovered:
    def test_a_second_foreground_genome_is_not_skipped(self):
        """The old loop `break`ed after the first, so a pan-target design
        reported one target's strand structure as if it were the panel's."""
        cache = _Cache({"a": _stats(score=0.7), "b": _stats(score=0.2)})

        stats = collect_strand_stats(cache, ["a", "b"], [1_000, 2_000], PRIMERS)

        assert set(stats) == {"a", "b"}
        assert stats["b"]["strand_alternation_score"] == 0.2

    def test_the_host_is_computed_too(self):
        """It never was, though the same method would have done it."""
        cache = _Cache({"target": _stats(), "host": _stats(gap_max=900.0)})

        stats = collect_strand_stats(cache, ["target", "host"], [3_000, 3_000_000], PRIMERS)

        assert stats["host"]["strand_alternation_gap_max"] == 900.0

    def test_each_genome_is_measured_against_its_own_length(self):
        """Passing one length for every prefix would scale the host's gaps by
        the target's size."""
        seen = []

        class _Recording(_Cache):
            def compute_strand_alternation_stats(self, prefix, primers, genome_length):
                seen.append((prefix, genome_length))
                return _stats()

        collect_strand_stats(_Recording({}), ["a", "b"], [111, 222], PRIMERS)

        assert seen == [("a", 111), ("b", 222)]


class TestAnUncomputedValueIsDistinguishable:
    def test_a_prefix_the_cache_cannot_answer_is_absent_rather_than_zero(self):
        cache = _Cache({"target": _stats()}, explode_on=["host"])

        stats = collect_strand_stats(cache, ["target", "host"], [3_000, 3_000], PRIMERS)

        assert "target" in stats
        assert "host" not in stats

    def test_a_cache_without_the_method_yields_nothing_at_all(self):
        """`StreamingPositionCache` has no such method, and metrics must still
        be computable rather than raising."""

        class _Bare:
            pass

        assert collect_strand_stats(_Bare(), ["target"], [3_000], PRIMERS) == {}

    def test_no_primers_yields_nothing(self):
        assert collect_strand_stats(_Cache({"t": _stats()}), ["t"], [1], []) == {}

    def test_mismatched_prefixes_and_lengths_yield_nothing(self):
        """Zipping a short length list would silently measure the wrong genomes."""
        cache = _Cache({"a": _stats(), "b": _stats()})

        assert collect_strand_stats(cache, ["a", "b"], [100], PRIMERS) == {}


class TestTheHeadlineScalars:
    """The two fields `PrimerSetMetrics` has carried all along."""

    def test_they_come_from_the_first_foreground_genome(self):
        stats = {"a": _stats(score=0.7, ratio=0.9), "b": _stats(score=0.1, ratio=0.2)}

        score, ratio = headline_strand_scalars(stats, ["a", "b"])

        assert (score, ratio) == (0.7, 0.9)

    def test_they_are_none_when_nothing_was_computed(self):
        """Not 0.0. The old initialiser made "never asked" indistinguishable
        from "measured zero", so a limit or a report could not tell them apart.
        """
        score, ratio = headline_strand_scalars({}, ["a"])

        assert score is None
        assert ratio is None

    def test_a_measured_zero_stays_zero(self):
        """A one-site panel has nothing to alternate, and that IS zero."""
        score, ratio = headline_strand_scalars({"a": _stats(score=0.0, ratio=0.0)}, ["a"])

        assert score == 0.0
        assert ratio == 0.0

    def test_a_host_only_result_does_not_masquerade_as_the_target(self):
        """The scalars describe the target. If only the host could be measured,
        they must stay None rather than report the host's balance."""
        score, ratio = headline_strand_scalars({"host": _stats(score=0.5)}, ["target"])

        assert score is None
        assert ratio is None


class TestItReachesTheMetrics:
    def test_the_metrics_carry_the_per_genome_stats(self):
        from neoswga.core.base_optimizer import PrimerSetMetrics

        assert "strand_stats" in PrimerSetMetrics.__dataclass_fields__

    def test_an_empty_metrics_has_no_strand_stats_rather_than_zeros(self):
        from neoswga.core.base_optimizer import PrimerSetMetrics

        empty = PrimerSetMetrics.empty()

        assert empty.strand_stats == {}
        assert empty.strand_alternation_score is None
        assert empty.strand_coverage_ratio is None

    def test_the_stats_reach_the_serialised_form(self):
        """`to_dict` feeds `step4_improved_df_summary.json` and the report."""
        from neoswga.core.base_optimizer import PrimerSetMetrics

        payload = PrimerSetMetrics.empty().to_dict()

        assert "strand_stats" in payload

    def test_the_collector_is_what_base_optimizer_calls(self):
        """Asserted on the call: both ends existing and the path not is exactly
        the class Known Issue 16 records."""
        import ast
        import inspect

        from neoswga.core import base_optimizer

        tree = ast.parse(inspect.getsource(base_optimizer))
        names = {
            node.func.id
            for node in ast.walk(tree)
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
        }

        assert "collect_strand_stats" in names
        assert "headline_strand_scalars" in names
