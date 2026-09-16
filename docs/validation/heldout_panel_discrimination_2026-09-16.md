# Do the optimised quantities separate effective panels from the rest?

Measured 2026-09-16. Plan step 252 of
`docs/superpowers/plans/2026-09-15-condition-aware-pool-design.md`, which asks
for a held-out assessment against simple baselines, held out by panel or study
rather than by splitting sites within one experiment, with negative results
reported.

Reproduce with `scripts/benchmarking/heldout_panel_discrimination.py`.

**The result is negative: no clear transferable advantage was established.** The
best quantity tested is five points above a majority-class baseline on held-out
accuracy, which is one panel out of twenty-one, and two quantities rank panels in
opposite directions in two studies from the same paper.

## What could be asked

The evidence matrix
(`docs/validation/evidence_matrix_2026-09-15.md`) rules out calibrating coverage
or genome recovery: no dataset here reports breadth at a stated depth. Three
datasets do carry a published per-panel site geometry alongside a binary
outcome, so a narrower question is available. Does a panel that scores better on
these quantities tend to be the one that worked?

| study | panels | effective | outcome as recorded |
|---|---|---|---|
| Clarke 2017, *M. tuberculosis* | 12 | 4 | the paper's footnote marking the most effective sets |
| Clarke 2017, *Wolbachia* | 5 | 1 | the authors' effective flag |
| Dwivedi-Yu 2023, *Prevotella* | 6 | 3 | reads on target above the control |

Twenty-three panels, eight of them effective, across three studies. Two of the
twelve *M. tuberculosis* entries are deliberate negative controls, so every
figure below is reported with and without them.

The predictors are the papers' own published numbers, not quantities recomputed
here. That makes this a check on the metric definitions rather than on this
tool's implementation of them, which is the stronger reading for a negative
result and the weaker one for a positive.

## Held-out accuracy

Held out by study: the threshold is chosen on two studies and applied to the
third, for each study in turn. The majority-class baseline is what you get by
predicting the commoner label every time.

| predictor | held-out accuracy | majority baseline |
|---|---|---|
| fg/bg ratio | 0.67 | 0.62 |
| background over foreground spacing | 0.57 | 0.62 |
| site density | 0.62 | 0.62 |
| panel size | 0.52 | 0.62 |
| foreground Gini | 0.52 | 0.62 |
| worst foreground gap | 0.52 | 0.62 |
| reach coverage | 0.48 | 0.62 |

Negative controls excluded. Including them moves every figure by at most 0.05,
and by that much only because the controls are designed to fail. The one
ordering it changes is the reach coverage proxy, which rises from last to a tie
for fifth. No conclusion here turns on it.

Nothing clears the baseline by a margin that 21 panels could support. The best,
the published fg/bg ratio, is five points above it, which is one panel. Saying
flatly that nothing beat the baseline would overstate this in the other
direction; the honest reading is that no advantage was established either way.

## Within-study ranking

Ranking is an easier question than thresholding, because it never has to
transfer a cut between studies. Area under the ROC curve, where 0.5 is no
information and 0.0 is a perfect inversion:

| predictor | Mtb | Wolbachia | Prevotella |
|---|---|---|---|
| fg/bg ratio | 0.00 | 1.00 | 0.89 |
| background over foreground spacing | 1.00 | 0.00 | 0.11 |
| worst foreground gap | 0.83 | 0.75 | 1.00 |
| site density | 1.00 | 0.50 | 0.78 |
| foreground Gini | 0.00 | 0.88 | 0.78 |
| panel size | 0.58 | 0.75 | 0.72 |
| reach coverage | 0.50 | 0.50 | 0.67 |

Negative controls excluded.

The first two rows are the finding. Each ranks one study perfectly and the other
perfectly backwards, and both studies come from the same paper with the same
polymerase. A quantity that behaves that way carries no signal that transfers;
it is reading something specific to each experiment.

The worst foreground gap is the only predictor above 0.7 in all three, and it is
also the one whose held-out accuracy is 0.52. Ranking within a study and
carrying a threshold between studies are different things, and only the second
is what a design tool needs.

## Reach coverage is at its ceiling

The coverage proxy is the weakest predictor in the table, and the reason is
mechanical rather than statistical. At a 3 kb reach, seventeen of the
twenty-three panels have a mean site spacing below 6 kb, so the proxy saturates
at 1.0 and cannot separate them.

| study | panels at the ceiling |
|---|---|
| Clarke *M. tuberculosis* | 10 of 12 |
| Clarke *Wolbachia* | 2 of 5 |
| Dwivedi-Yu *Prevotella* | 5 of 6 |

Published panels are already dense enough that geometric coverage is close to
complete, which agrees with what the marginal coverage table shows on this
tool's own designs: coverage rises monotonically and flattens, while whatever
separates a good panel from a bad one is elsewhere.

## Limitations

- Twenty-three panels and eight positives. Every figure here has a confidence
  interval wide enough to contain the baseline.
- Three studies means three held-out folds, and one of them has a single
  positive panel, so its area under the curve is the rank of one observation.
- The outcome labels are not one endpoint. Two are the authors' own judgement of
  which sets were most effective and one is reads on target against a control.
  Pooling them assumes they mean the same thing, and they do not.
- Predictors are published values from three papers, computed by their authors
  under their own definitions. They are not this tool's measurements of the same
  panels.
- No occupancy-weighted quantity appears here, because none of these datasets
  reports one and recomputing it would need the reaction conditions the papers
  do not fully state.

## Position

This does not show the tool's objective is wrong. It shows the available
observations cannot tell, and that the simple baselines cannot tell either. The
design quantities remain what they were described as: modelled site geometry,
not measured amplification.

The thing that would change this is in the evidence matrix under what would be
needed for calibration, and it is not obtainable by re-analysing what is held
here.
