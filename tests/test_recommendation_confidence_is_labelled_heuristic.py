"""A score threshold is not a confidence level.

Audit finding F8. `AdditiveOptimizer` mapped its own search score to the words
"high", "medium" and "low" confidence with no uncertainty model behind them: a
score above 0.5 was reported as high confidence in a recipe that has never been
compared to a measured reaction.

The number is a useful ranking and is kept. What changes is that it is named for
what it is, so a reader is not told the tool is confident when what it means is
that its own heuristic scored well.
"""

from neoswga.core.additive_optimizer import AdditiveRecommendation


def test_the_field_is_named_as_a_heuristic_band():
    assert hasattr(AdditiveRecommendation, "__dataclass_fields__")
    fields = AdditiveRecommendation.__dataclass_fields__
    assert "heuristic_score_band" in fields, (
        "the score band should not be called 'confidence'; it carries no " "uncertainty model"
    )


def test_the_rendered_recommendation_does_not_claim_confidence():
    rec = AdditiveRecommendation(heuristic_score_band="high")
    summary = rec.summary()
    text = summary if isinstance(summary, str) else "\n".join(summary)

    assert "Confidence:" not in text
    assert "Heuristic score band" in text or "heuristic" in text.lower()
