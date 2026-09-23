"""Item 6, delivered as the ingredient rather than the recipe.

The ranked list's item 6 was occupancy-weighted spacing: a gap bounded by a
site occupied a tenth of the time is not really bounded, and no other tool can
express that. It was ranked last and explicitly least supported, and two
measurements say it should not be built as proposed.

**No weighting formula can be validated.** The 18 published sets with wet-lab
outcomes carry PUBLISHED gap figures, not binding positions, so a weighted gap
cannot be recomputed for them. The unweighted statistic already fails to
separate their winners at every reach-derived threshold
(`scripts/benchmarking/published_gap_thresholds.py`), and weighting a quantity
that does not separate winners gives a better-founded quantity that still does
not. Choosing a weighting rule here would be exactly the unvalidated scoring
change this project has refused all along.

**It would also be a hot-path rewrite.** `_compute_metrics` pools site
positions per prefix and discards which primer each came from, so a weighted
gap needs the gap computation restructured to carry primer identity.

What IS worth having is the ingredient, which is cheap and actionable: how
weakly the weakest primer in the delivered panel is bound. Measured, occupancy
spans 7 to 9 fold within a panel at equiphi29 42 C and only 1.7 to 3.0 fold at
phi29 30 C, so this is a real quantity on the platform where additives work and
nearly constant on the one where they do not.
"""

import pytest

from neoswga.core.panel_regime import assess_panel


class _Metrics:
    def __init__(self, primer_occupancy=None, **kw):
        self.effective_fg_coverage = kw.get("effective_fg_coverage", 0.85)
        self.fg_coverage = 0.9
        self.selectivity_ratio = 8.0
        self.selectivity_density = 45.0
        self.total_bg_sites = 120
        self.max_gap = 31_000.0
        self.mean_gap = 4_800.0
        self.gap_gini = 0.48
        self.bg_coverage = 0.02
        self.strand_coverage_ratio = 0.9
        self.strand_alternation_score = 0.7
        self.strand_stats = {}
        self.primer_occupancy = {} if primer_occupancy is None else dict(primer_occupancy)


def _assess(metrics=None, **overrides):
    kwargs = dict(
        coverage_target=0.80,
        min_fg_bg_ratio=5.0,
        requested_size=12,
        delivered_size=12,
        coverage_reach=3_000,
        genome_length=3_200_000,
    )
    kwargs.update(overrides)
    return assess_panel(metrics or _Metrics(), **kwargs)


def _named(regime, name):
    for criterion in regime.criteria:
        if criterion.name == name:
            return criterion
    raise AssertionError(f"no criterion named {name!r} in {[c.name for c in regime.criteria]}")


class TestTheMetricsCarryIt:
    def test_the_field_exists(self):
        from neoswga.core.base_optimizer import PrimerSetMetrics

        assert "primer_occupancy" in PrimerSetMetrics.__dataclass_fields__

    def test_an_empty_metrics_has_none_rather_than_zeros(self):
        from neoswga.core.base_optimizer import PrimerSetMetrics

        assert PrimerSetMetrics.empty().primer_occupancy == {}

    def test_it_reaches_the_serialised_form(self):
        from neoswga.core.base_optimizer import PrimerSetMetrics

        assert "primer_occupancy" in PrimerSetMetrics.empty().to_dict()

    def test_base_optimizer_computes_it(self):
        """Asserted on the call, because a field nothing populates is the
        defect class this project keeps closing."""
        import ast
        import inspect

        from neoswga.core import base_optimizer

        tree = ast.parse(inspect.getsource(base_optimizer))
        names = {
            node.func.id
            for node in ast.walk(tree)
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
        }

        assert "panel_occupancy" in names


class TestTheWeakestPrimerIsReported:
    def test_it_appears_when_occupancy_was_computed(self):
        regime = _assess(_Metrics({"AAAA": 0.12, "GGGG": 0.95, "ACGT": 0.60}))

        criterion = _named(regime, "weakest_occupancy")
        assert criterion.value == pytest.approx(0.12)

    def test_it_carries_no_reference(self):
        """There is no validated floor on occupancy, so it is reported and
        never ranked, like the other unreferenced criteria."""
        regime = _assess(_Metrics({"AAAA": 0.12, "GGGG": 0.95}))

        criterion = _named(regime, "weakest_occupancy")
        assert criterion.reference is None
        assert criterion.slack is None
        assert "weakest_occupancy" in regime.unreferenced

    def test_the_note_carries_the_median_and_the_spread(self):
        """The minimum alone cannot say whether the panel is uniformly weak or
        carries one outlier, and that difference is what a user acts on."""
        regime = _assess(_Metrics({"A": 0.10, "B": 0.50, "C": 0.90}))
        note = _named(regime, "weakest_occupancy").note

        assert "median" in note.lower()
        assert "0.5" in note

    def test_it_is_absent_when_occupancy_could_not_be_computed(self):
        """No conditions means no temperature at which to evaluate occupancy,
        and a fabricated zero would be indistinguishable from a measurement."""
        regime = _assess(_Metrics({}))

        assert "weakest_occupancy" not in [c.name for c in regime.criteria]

    def test_limiting_never_names_it(self):
        regime = _assess(_Metrics({"A": 0.001}))

        assert regime.limiting != "weakest_occupancy"

    def test_it_reaches_the_printed_output(self):
        from neoswga.core.panel_regime import format_regime

        regime = _assess(_Metrics({"A": 0.12, "B": 0.95}))
        lines = "\n".join(format_regime(regime))

        assert "weakest_occupancy" in lines


class TestTheHelper:
    def test_it_returns_a_reading_per_primer(self):
        from neoswga.core.reaction_conditions import ReactionConditions
        from neoswga.core.strand_metrics import panel_occupancy

        primers = ["ACGTACGTACGT", "GGGGCCCCGGGG"]
        result = panel_occupancy(primers, ReactionConditions(temp=30.0))

        assert set(result) == set(primers)
        assert all(0.0 <= v <= 1.0 for v in result.values())

    def test_no_conditions_yields_nothing(self):
        from neoswga.core.strand_metrics import panel_occupancy

        assert panel_occupancy(["ACGTACGTACGT"], None) == {}

    def test_no_primers_yields_nothing(self):
        from neoswga.core.reaction_conditions import ReactionConditions
        from neoswga.core.strand_metrics import panel_occupancy

        assert panel_occupancy([], ReactionConditions(temp=30.0)) == {}

    def test_a_gc_rich_primer_is_more_occupied_at_a_cold_reaction(self):
        """Sanity, so the helper is not returning a constant."""
        from neoswga.core.reaction_conditions import ReactionConditions
        from neoswga.core.strand_metrics import panel_occupancy

        result = panel_occupancy(["ATATATATATAT", "GCGCGCGCGCGC"], ReactionConditions(temp=30.0))

        assert result["GCGCGCGCGCGC"] > result["ATATATATATAT"]
