"""Guards for calibrating a model against sequencing runs.

Task 9 of the 2026-09-21 valid-design plan. A calibration fitted on a set of
experiments and then evaluated on one of them reports an in-sample score and
calls it held out. The failure is silent in the worst way: the number is good,
which is exactly what makes nobody check it.

This module holds the checks that have to run before a fitted model may claim
anything. It deliberately does not hold the fitting: a guard that lives inside
the thing it guards is one refactor away from being skipped.
"""

from __future__ import annotations

from typing import Sequence

__all__ = ["require_disjoint_experiments"]


def require_disjoint_experiments(training: Sequence[str], validation: Sequence[str]) -> None:
    """Refuse a validation set that shares an experiment with the training set.

    Both sets must be non-empty. An empty validation set would pass a
    disjointness check trivially and then produce an "out-of-sample" result
    computed on no samples, which is a stronger claim than an in-sample one
    and rests on less.

    A repeated identifier WITHIN one set is untidy rather than leakage, so it
    is allowed: listing a run twice does not make it inform its own
    evaluation.

    The overlapping identifiers are named. A guard that says only that an
    overlap exists leaves the user to find which of twenty runs it was.
    """
    training_set = {str(item) for item in training}
    validation_set = {str(item) for item in validation}

    if not training_set:
        raise ValueError(
            "The training experiment set is empty; there is nothing to fit a " "calibration on."
        )
    if not validation_set:
        raise ValueError(
            "The validation experiment set is empty. Nothing held out is not the "
            "same as nothing overlapping: an out-of-sample result computed on no "
            "samples is not a weaker claim than an in-sample one, it is an "
            "unsupported one."
        )

    shared = sorted(training_set & validation_set)
    if shared:
        raise ValueError(
            f"Training and validation experiments overlap: {', '.join(shared)}. "
            f"A run that informed the fit cannot evaluate it; the score would be "
            f"in-sample while the report called it held out."
        )
