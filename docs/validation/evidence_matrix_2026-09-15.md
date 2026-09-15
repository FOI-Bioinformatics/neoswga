# What the available observations can and cannot validate

Task 8 of the [condition-aware pool design plan](../superpowers/plans/2026-09-15-condition-aware-pool-design.md).
Compiled by inspecting the datasets already in `tests/validation/data/` before
requesting anything new.

The question this answers is narrow and worth stating plainly: for each
quantity the tool predicts, is there an observation here that could confirm or
refute it? For most of them the answer is no, and the reason is usually that the
observation measures a different endpoint rather than that it is unavailable.

## The datasets

| source | endpoint | polymerase | oligo lengths | conditions available | reference available | depth available | independent unit | supported use | limitations |
|---|---|---|---|---|---|---|---|---|---|
| Dwivedi-Yu 2023, plasmid | qPCR fold amplification, 310 measurements | phi29 | 8-12 | partial, one buffer | yes, two plasmids | no | one oligo | sequence-level amplification | Site geometry is constant by construction, so it cannot separate a coverage or selectivity model from a sequence model |
| Dwivedi-Yu 2023, Prevotella | % reads on target, fold enrichment vs control, 6 sets | phi29 | 10-12 | partial | yes | not per base | one panel | panel-level enrichment ranking | 6 panels is few for held-out evaluation; depth per base is not reported |
| Clarke 2017, M. tuberculosis | effective / not, 12 sets | phi29 | 8-12 | partial | yes | no | one panel | binary discrimination | An ordinal "most effective" flag, not a continuous endpoint |
| Clarke 2017, Wolbachia | effective / not, with tier, 5 sets | phi29 | 8-12 | partial | yes | no | one panel | binary discrimination | Tiers are the authors' judgement, not a measured scale |
| Leichty & Brisson 2014 | fold target and fold background, 3 sets | phi29 | 8-12 | partial | yes | no | one panel | selectivity direction | 3 panels; fold enrichment is not breadth |
| Dwivedi-Yu 2023, Leishmania | effective / not, 4 sets | phi29 | 10-12 | partial | yes | no | one panel | binary discrimination | 4 panels |
| Oyola 2016, P. falciparum | fold enrichment and % reads on target, 1 set | phi29 | 10-12 | partial | yes | no | one panel | a single anchor point | One panel cannot be split into training and held-out halves |

## What each predicted quantity can be checked against

| predicted quantity | observation available | verdict |
|---|---|---|
| Sequence-level amplification | 310 qPCR measurements | Checkable. Primer length and melting temperature already correlate at rho +0.52 and +0.43 |
| Panel-level enrichment ranking | 6 Prevotella panels, 1 P. falciparum | Weakly checkable, as a ranking. Too few for held-out evaluation |
| Binary effective / not | 24 panels across four studies | Checkable as discrimination, not as calibration |
| Occupancy-weighted coverage | none | **Not checkable.** No dataset reports breadth at a stated depth |
| Selectivity density | none | **Not checkable.** Fold enrichment is a ratio of yields, not of binding-site densities |
| Genome recovery at a depth threshold | none | **Not checkable.** No dataset reports depth per base |
| Additive effect at 30 C or 42 C | none | **Not checkable.** Every source is single-condition |
| Extension reach | one breadth proxy, one organism | Fitted once, not validated. The 3 kb figure is a design-density convention |

## The distinction that matters most

**qPCR fold amplification is not breadth of genome.** It is the most abundant
observation here, 310 measurements against single-digit panel counts everywhere
else, and it is tempting to treat it as the calibration set. It measures how
much product one oligo generates from a plasmid whose site geometry is constant
by construction. Coverage, selectivity density and genome recovery are all
properties of where sites fall across a genome, which this dataset holds fixed.

Using it to calibrate a coverage model would produce a number that looked
validated and was not.

## What would be needed to calibrate genome recovery

Not held here, and not obtainable by re-analysing what is:

- observed breadth at declared depth thresholds, per base or per window;
- sequencing effort, so breadth can be separated from depth;
- the reference and mapping definition used to call a base covered;
- input composition and DNA quality, since host fraction bounds what is
  achievable regardless of the panel;
- panel identity, polymerase, and the full buffer and additive composition;
- the concentration convention, per oligo or total pool, which changes every
  melting temperature in the model.

Held-out units must be complete reactions or studies. Splitting a single depth
profile into training and test positions measures interpolation within one
experiment and would report a much better error than the model deserves.

## Position

The software changes in tasks 1 to 7 are correctness work, verified against
their own units and against small fixtures. They stand on that.

Claims of measured genome recovery are withheld, and no threshold in the tool
should be read as a recovery target. That is an evidence limitation, not a
defect in the code, and it is not a reason to delay the correctness work.
