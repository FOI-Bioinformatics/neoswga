"""Validation against published SWGA primer sets with wet-lab outcomes.

Data: Dwivedi-Yu et al. (2023) PLOS Comput Biol 19(4):e1010137, S2 Table (six
primer sets for *Prevotella melaninogenica* with their computed statistics) and
S3 Table (percent of reads mapping to target after amplification). This is a
rare thing to have: designed sets, their in-silico metrics, and a measured
outcome for each.

Two things this suite establishes.

**1. Which metric actually predicts success.** The audit plan - and the guidance
originally written from the paper's own summary - assumed binding-site density
(lowest mean binding distance) was the dominant predictor. The published data
does not support that. Across these six sets mean binding distance has almost no
rank correlation with enrichment (rho = -0.09), and ranking by it alone puts the
two *worst-performing* sets first. What separates the winners is the largest
coverage hole (max gap, rho = -0.83) and binding-site evenness (Gini,
rho = -0.77).

That is mechanistically sensible: a region no primer reaches stays unamplified
no matter how dense the average spacing is elsewhere. Prev04 has the best mean
spacing of all six (1980 bp) and finished third; Prev06 has the second-worst
(4980 bp) and won outright.

**2. That the shipped scoring does rank the winners first - but narrowly.**
``PrimerSetMetrics.normalized_score`` surfaces Prev03 and Prev06 as the top two,
which is correct. The margin over the next set is only a few points of score
against a 60-fold spread in measured performance, so the ordering should not be
read as a confident separation on six sets.

Note ``normalized_score`` does not use ``max_gap`` despite computing it, and
``max_gap`` is the strongest single predictor here. Adding a hole term was tried
and made no difference to the ranking on this data, so it was not shipped: a
change that re-ranks all optimizer output needs to demonstrate a benefit first.

**Caveat, stated plainly: n = 6.** These correlations are suggestive, not
established, and any re-weighting tuned on six points risks fitting noise. The
value here is that the *current* ranking is demonstrably wrong on the only real
outcome data available, not that the specific alternative weights are right.
"""

import json
from pathlib import Path

import pytest

DATA = Path(__file__).parent / "data" / "dwivedi_yu_2023_prevotella.json"
WET_LAB_WINNERS = {"Prev03", "Prev06"}
GENOME_BP = 3_168_282  # P. melaninogenica ATCC 25845
BACKGROUND_BP = 3_100_000_000  # human, the background used in the paper


@pytest.fixture(scope="module")
def published():
    return json.loads(DATA.read_text())


@pytest.fixture(scope="module")
def sets(published):
    return published["sets"]


def _spearman(x, y):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0] * len(v)
        for pos, i in enumerate(order):
            r[i] = pos + 1
        return r

    rx, ry = rank(x), rank(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den else 0.0


# ----------------------------------------------------------------------
# The dataset itself
# ----------------------------------------------------------------------


def test_fixture_matches_the_published_tables(sets):
    """Sanity-check the transcription from S2/S3 before drawing conclusions."""
    assert len(sets) == 6
    for name, s in sets.items():
        assert len(s["primers"]) == s["published"]["size"], name
        assert all(set(p) <= set("ACGT") for p in s["primers"]), name
    # S3 Table figures
    assert sets["Prev06"]["outcome"]["pct_reads_on_target"] == 72.22
    assert sets["Prev03"]["outcome"]["pct_reads_on_target"] == 57.34
    assert sets["Prev05"]["outcome"]["pct_reads_on_target"] == 1.20


def test_success_labels_match_the_papers_grouping(sets):
    successful = {n for n, s in sets.items() if s["outcome"]["successful"]}
    assert successful == WET_LAB_WINNERS | {
        "Prev04"
    }, "the paper reports Prev03/Prev04/Prev06 as the effective sets"


# ----------------------------------------------------------------------
# Which metric predicts the outcome
# ----------------------------------------------------------------------


def test_max_gap_predicts_enrichment(sets):
    """The largest coverage hole is the strongest single predictor here."""
    names = list(sets)
    y = [sets[n]["outcome"]["fold_enrichment_vs_control"] for n in names]
    x = [sets[n]["published"]["fg_max_distance"] for n in names]
    rho = _spearman(x, y)
    assert rho < -0.7, f"expected strong negative correlation, got {rho:+.3f}"


def test_evenness_predicts_enrichment(sets):
    names = list(sets)
    y = [sets[n]["outcome"]["fold_enrichment_vs_control"] for n in names]
    x = [sets[n]["published"]["fg_gini"] for n in names]
    assert _spearman(x, y) < -0.7


def test_mean_binding_distance_does_not_predict_enrichment(sets):
    """Guards against reinstating the assumption this data refutes.

    If a future dataset shows mean distance does predict outcome, this test
    should fail loudly so the claim gets re-examined rather than quietly
    reintroduced into the docs.
    """
    names = list(sets)
    y = [sets[n]["outcome"]["fold_enrichment_vs_control"] for n in names]
    x = [sets[n]["published"]["fg_mean_distance"] for n in names]
    rho = _spearman(x, y)
    assert abs(rho) < 0.4, (
        f"mean binding distance now correlates (rho={rho:+.3f}); revisit the "
        f"guidance in the swga-chemistry skill and the scoring weights"
    )


def test_ranking_by_mean_distance_alone_picks_losers(sets):
    """Concretely: the metric we used to recommend selects the wrong sets."""
    order = sorted(sets, key=lambda n: sets[n]["published"]["fg_mean_distance"])
    assert set(order[:2]) & WET_LAB_WINNERS == set(), (
        "ranking by mean binding distance should NOT surface the winners - it "
        "picks Prev04 and Prev05, which achieved 12.8x and 2.0x"
    )


def test_ranking_by_max_gap_recovers_the_top_three(sets):
    order = sorted(sets, key=lambda n: sets[n]["published"]["fg_max_distance"])
    assert order[:3] == [
        "Prev06",
        "Prev03",
        "Prev04",
    ], f"max-gap ranking should reproduce the measured top three, got {order[:3]}"


# ----------------------------------------------------------------------
# Does neoswga's scoring reproduce the outcome?
# ----------------------------------------------------------------------


def _metrics_from_published(p):
    """Build PrimerSetMetrics from the published statistics.

    Coverage and site counts are derived from the reported spacing rather than
    recomputed from the genome, so this exercises the SCORING function without
    requiring a 3 Mb download at test time.
    """
    from neoswga.core.base_optimizer import PrimerSetMetrics

    fg_sites = max(1, round(GENOME_BP / p["fg_mean_distance"]))
    # The background was the HUMAN genome, so background sites must be counted
    # over 3.1 Gb, not over the 3.2 Mb target. Dividing both by the target size
    # understated background sites ~1000-fold and inverted the selectivity
    # ordering. Note the paper's own "fg/bg ratio" column is a ratio of mean
    # SPACINGS (fg_mean/bg_mean), which is a different quantity from
    # PrimerSetMetrics.selectivity_ratio (fg_sites/bg_sites).
    bg_sites = max(1, round(BACKGROUND_BP / p["bg_mean_distance"]))
    return PrimerSetMetrics(
        fg_coverage=min(1.0, fg_sites * 3000 * 2 / GENOME_BP),
        bg_coverage=0.0,
        coverage_uniformity=p["fg_gini"],
        total_fg_sites=fg_sites,
        total_bg_sites=bg_sites,
        selectivity_ratio=fg_sites / bg_sites,
        mean_tm=30.0,
        tm_range=(28.0, 32.0),
        dimer_risk_score=0.0,
        mean_gap=p["fg_mean_distance"],
        max_gap=p["fg_max_distance"],
        gap_gini=p["fg_gini"],
        gap_entropy=1.0,
        strand_alternation_score=0.5,
        strand_coverage_ratio=0.5,
    )


def test_normalized_score_is_computable_for_every_published_set(sets):
    for name, s in sets.items():
        score = _metrics_from_published(s["published"]).normalized_score()
        assert 0.0 <= score <= 1.0, name


@pytest.mark.parametrize("application", ["balanced", "clinical", "enrichment"])
def test_normalized_score_ranks_the_wet_lab_winners_first(sets, application):
    """The shipped scoring does surface the two wet-lab winners.

    This was briefly recorded as a known failure. That was wrong, and the cause
    was a bug in the reconstruction above, not in the scoring: background sites
    were counted over the 3.2 Mb target instead of the 3.1 Gb human background,
    which understated them ~1000-fold and inverted the selectivity ordering.
    With the background sized correctly, normalized_score ranks Prev03 and
    Prev06 first without any change to its weights.
    """
    scored = sorted(
        sets,
        key=lambda n: -_metrics_from_published(sets[n]["published"]).normalized_score(
            application=application
        ),
    )
    assert (
        set(scored[:2]) == WET_LAB_WINNERS
    ), f"{application}: top-2 was {scored[:2]}, expected {sorted(WET_LAB_WINNERS)}"


def test_normalized_score_barely_discriminates_between_these_sets(sets):
    """Documents the magnitude of the problem, not just its direction.

    Actual performance spans 60-fold (2.0x to 120x). If the score spread is
    tiny, the metric cannot be separating good sets from bad regardless of
    whether the ordering happens to come out right.
    """
    scores = [_metrics_from_published(s["published"]).normalized_score() for s in sets.values()]
    spread = max(scores) - min(scores)
    assert spread < 0.10, (
        f"score spread is {spread:.3f}; if this has widened, the scoring was "
        f"changed and this test should be revisited alongside the xfail above"
    )


def test_max_gap_is_available_to_the_scorer(sets):
    """The fix is reachable: the metric exists, it is simply unused."""
    from dataclasses import fields

    from neoswga.core.base_optimizer import PrimerSetMetrics

    assert "max_gap" in {f.name for f in fields(PrimerSetMetrics)}
    m = _metrics_from_published(sets["Prev06"]["published"])
    assert m.max_gap == 31300
