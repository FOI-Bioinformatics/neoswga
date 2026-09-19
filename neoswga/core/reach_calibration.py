"""Estimate the per-primer coverage reach from sequencing depth.

The reach is the single most consequential number in a coverage figure and the
least evidenced. The same primer set measures 0.418 coverage at 3 kb, 0.836 at
10 kb and ~1.0 at 70 kb, so which value is used decides both what gets reported
and, because set-cover selects at it, which primers get chosen.

Neither convention in the literature is measured. swga 2.0 uses phi29's ~70 kb
single-molecule processivity, which is how far one molecule *can* be extended,
not how far a primer's product reaches before a neighbour's strand displacement
truncates it. NeoSWGA's ~3 kb default is derived from the inter-primer
*spacings* that Clarke et al. (2017) and Dwivedi-Yu et al. (2023) designed
for -- but a spacing is a choice a designer made, not a reach anyone observed.

A BAM from a real SWGA run does contain the answer. If a primer's product
reaches R bases, sequencing depth should fall away from binding sites on a
length scale of R. Fitting R to observed depth replaces both conventions with a
measurement for the reaction that produced the data.

The estimate is only as good as its input: it needs a BAM from an SWGA reaction
using the primer set given, mapped to the target. Depth from a uniform WGA
library carries no primer-position signal and will not produce a meaningful
number -- `fit_reach` reports the correlation so a caller can tell the
difference rather than trusting the returned value blindly.

The number it returns is a MODEL PARAMETER fitted to one dataset, not a
measured constant of the polymerase. Treat it as the centre of a sensitivity
analysis: `plausible_reaches` gives the range the data cannot separate, and
`cv_correlation` says how much of the in-sample fit survives being held out.
The reported `correlation` is the maximum over the grid, computed on the same
data that chose the reach, so it is optimistic by construction.

Measured on null depth over twelve seeds, that optimism averages +0.058 --
about a third of the 0.15 informativeness threshold, purely from choosing among
eleven candidates. It is not, however, larger than the fold-to-fold noise in
the held-out estimate (sd 0.052), so on any single run the held-out figure may
land above the in-sample one. Read the two together as a stability check rather
than treating either as a corrected value.
"""

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

logger = logging.getLogger(__name__)

# Reaches tried when fitting. Geometric rather than linear: the quantity spans
# more than an order of magnitude and the interesting differences are
# multiplicative, so a linear grid would waste most of its points above 20 kb.
DEFAULT_REACH_GRID = (500, 1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 35000, 70000)

# Below this Spearman correlation between predicted and observed depth, the fit
# is not evidence about reach. Deliberately not a p-value: with millions of
# bases even a trivial correlation is significant, and significance is not the
# question -- whether the primer positions explain the depth profile is.
MIN_INFORMATIVE_CORRELATION = 0.15

# Fewer bins than this and a Spearman correlation is not evidence about reach.
# Two bins correlate perfectly with anything monotone; eight is the point at
# which a rank correlation has enough distinct orderings to be worth reading.
# Below it the fit reports what it has rather than returning a number, which is
# also what keeps a target shorter than one bin from reaching the reshape that
# used to raise a bare numpy ValueError.
MIN_BINS_TO_FIT = 8

# Folds for the blocked cross-validation, and the smallest block worth scoring.
CV_FOLDS = 5
MIN_CV_BLOCK_BINS = 3

# Two reaches whose correlations differ by less than this are not distinguished
# by the data. Reported as a range rather than resolved into a point, because a
# flat optimum and a sharp one must not read as the same answer.
PLAUSIBLE_MARGIN = 0.02


@dataclass(frozen=True)
class ReachFit:
    """Result of fitting reach to observed depth.

    `correlation` is the IN-SAMPLE figure: the maximum over the grid, on the
    data that chose the reach, and so optimistic by construction.
    `cv_correlation` is the same quantity measured on blocks held out of that
    choice. A large gap between them means the fit is selection noise; agreement
    means the profile carries real primer-position signal.
    """

    best_reach: int
    correlation: float
    by_reach: Dict[int, float]
    informative: bool
    contig: str
    note: str
    cv_correlation: Optional[float] = None
    plausible_reaches: Tuple[int, ...] = ()
    at_grid_edge: bool = False
    n_bins: int = 0

    @property
    def in_sample_correlation(self) -> float:
        """What `correlation` has always been, named for what it is."""
        return self.correlation

    def to_dict(self) -> Dict:
        return {
            "best_reach": self.best_reach,
            "correlation": self.correlation,
            "in_sample_correlation": self.correlation,
            "cv_correlation": self.cv_correlation,
            "plausible_reaches": list(self.plausible_reaches),
            "at_grid_edge": self.at_grid_edge,
            "n_bins": self.n_bins,
            "by_reach": {str(k): v for k, v in self.by_reach.items()},
            "informative": self.informative,
            "contig": self.contig,
            "note": self.note,
        }


def predicted_depth(
    positions: Sequence[int], length: int, reach: int, circular: bool = False
) -> np.ndarray:
    """Expected relative depth if each binding site contributes over +/- reach.

    A triangular kernel rather than a rectangular one: product length is a
    distribution, not a hard cutoff, so depth should taper with distance from
    the site rather than stopping. The shape matters less than the scale -- the
    fit is choosing between reaches an order of magnitude apart -- but a
    rectangular kernel makes every reach predict a flat profile once sites
    overlap, which removes the signal the fit runs on.

    `circular` wraps the kernel instead of clipping it. On a plasmid or a
    bacterial chromosome the clipped version loses the part of each near-origin
    triangle that falls off the end, and it loses MORE of it the larger the
    reach -- a bias against exactly the large reaches the fit exists to weigh.
    It is off by default because a multi-record prefix is a concatenation, and
    there the end of the last record is not adjacent to the start of the first.
    """
    profile = np.zeros(length, dtype=np.float64)
    if reach <= 0 or length <= 0:
        return profile
    # One triangle per site, accumulated. Vectorised per site rather than per
    # base: a 3 Mb genome with 500 sites is 500 slice adds, not 1.5e9 updates.
    ramp = 1.0 - np.abs(np.arange(-reach, reach + 1, dtype=np.float64)) / reach
    for pos in positions:
        lo, hi = int(pos) - reach, int(pos) + reach + 1
        if circular and (lo < 0 or hi > length):
            # `np.add.at` rather than fancy-index assignment because a reach
            # wider than the target wraps more than once and the repeated
            # indices must accumulate rather than overwrite.
            np.add.at(profile, np.arange(lo, hi) % length, ramp)
            continue
        klo, khi = 0, ramp.size
        if lo < 0:
            klo = -lo
            lo = 0
        if hi > length:
            khi -= hi - length
            hi = length
        if hi > lo and khi > klo:
            profile[lo:hi] += ramp[klo:khi]
    return profile


def _bin_starts(length: int, bin_size: int, record_starts: Optional[Sequence[int]]) -> np.ndarray:
    """Offsets at which bins begin, restarting at every record boundary.

    Binning the concatenated array as one run let a bin span the end of one
    molecule and the start of the next, averaging depth across a join no
    polymerase crosses. Restarting at each boundary costs a partial bin per
    record and buys an observation that describes one molecule.

    Every base lands in a bin, including the trailing remainder. `usable =
    n_bins * bin_size` discarded up to a bin of the genome, which on the
    bundled plasmid is most of it.
    """
    bounds = sorted({0} | {int(s) for s in (record_starts or ()) if 0 < int(s) < length})
    starts: List[int] = []
    for start, end in zip(bounds, bounds[1:] + [length]):
        starts.extend(range(start, end, bin_size))
    return np.asarray(starts, dtype=np.int64)


def _binned(values: np.ndarray, starts: np.ndarray, widths: np.ndarray) -> np.ndarray:
    """Mean of `values` over each bin, partial bins included.

    `np.add.reduceat` rather than a reshape, because the bins are not all the
    same width once record boundaries and the trailing remainder are honoured.
    """
    return np.add.reduceat(np.asarray(values, dtype=np.float64), starts) / widths


def _spearman(predicted: np.ndarray, observed: np.ndarray) -> float:
    """Rank correlation, with the degenerate cases reported as no signal.

    Ranks rather than values because absolute depth depends on sequencing
    effort, which says nothing about reach.
    """
    from scipy.stats import spearmanr

    if predicted.size < 3 or observed.size < 3:
        return 0.0
    # A constant array has no ranks to correlate. It arises when the reach is
    # large enough relative to the target that every bin is predicted equally,
    # which is a real answer about that reach: it carries no information.
    if np.allclose(predicted, predicted[0]) or np.allclose(observed, observed[0]):
        return 0.0
    rho = spearmanr(predicted, observed).statistic
    return 0.0 if np.isnan(rho) else float(rho)


def _blocked_cv_correlation(
    binned_by_reach: Dict[int, np.ndarray], observed: np.ndarray
) -> Optional[float]:
    """How much of the in-sample fit survives being held out.

    The reported correlation was the maximum over eleven candidates, computed
    on the data that chose among them, so it is optimistic by construction and
    nothing said so. Here each fold re-runs the whole selection on the
    remaining bins and scores the reach it picks on the block it never saw.
    That measures the PROCEDURE, which is what a caller applies.

    The blocks are contiguous rather than interleaved. Neighbouring bins share
    the same binding sites and the same local mappability, so a random split
    would put near-duplicates on both sides and report a held-out figure barely
    below the in-sample one -- which is the reassurance this exists to avoid.
    """
    n_bins = observed.size
    folds = min(CV_FOLDS, n_bins // MIN_CV_BLOCK_BINS)
    if folds < 2:
        return None

    edges = np.linspace(0, n_bins, folds + 1).astype(int)
    held_out: List[float] = []
    for lo, hi in zip(edges, edges[1:]):
        if hi - lo < MIN_CV_BLOCK_BINS:
            continue
        keep = np.ones(n_bins, dtype=bool)
        keep[lo:hi] = False
        if keep.sum() < MIN_CV_BLOCK_BINS:
            continue
        chosen = max(
            binned_by_reach,
            key=lambda r: _spearman(binned_by_reach[r][keep], observed[keep]),
        )
        held_out.append(_spearman(binned_by_reach[chosen][lo:hi], observed[lo:hi]))

    return float(np.mean(held_out)) if held_out else None


def fit_reach(
    depth: np.ndarray,
    positions: Sequence[int],
    contig: str = "",
    reaches: Optional[Sequence[int]] = None,
    bin_size: int = 1000,
    record_starts: Optional[Sequence[int]] = None,
    circular: bool = False,
) -> ReachFit:
    """Choose the reach whose predicted depth profile best matches the observed one.

    Ranked by Spearman correlation, which cares about the shape of the profile
    rather than its scale -- absolute depth depends on sequencing effort, which
    says nothing about reach.

    Binning to `bin_size` before correlating: per-base depth is dominated by
    mapping noise at a scale far below any candidate reach, and correlating raw
    bases mostly measures that noise. `record_starts` keeps a bin from spanning
    two records of a concatenated prefix.

    The result is a model parameter, not a measured constant. Read
    `cv_correlation` rather than `correlation`, and `plausible_reaches` rather
    than `best_reach` alone.
    """
    reaches = list(reaches if reaches is not None else DEFAULT_REACH_GRID)
    length = len(depth)
    if length == 0 or not len(positions):
        return ReachFit(
            best_reach=0,
            correlation=0.0,
            by_reach={},
            informative=False,
            contig=contig,
            note="No depth or no binding sites; nothing to fit.",
        )

    starts = _bin_starts(length, bin_size, record_starts)
    widths = np.diff(np.append(starts, length)).astype(np.float64)
    n_bins = int(starts.size)

    caveats: List[str] = []
    if circular and n_bins and len(set(int(s) for s in (record_starts or ()))) > 1:
        caveats.append(
            "Circular wrapping was refused: this prefix holds more than one "
            "record, so "
            "the end of the last is not adjacent to the start of the first. The "
            "kernel clips at the ends of the concatenation instead."
        )
        circular = False

    if n_bins < MIN_BINS_TO_FIT:
        return ReachFit(
            best_reach=0,
            correlation=0.0,
            by_reach={},
            informative=False,
            contig=contig,
            n_bins=n_bins,
            note=" ".join(
                [
                    f"Only {n_bins} bin(s) of {bin_size:,} bp: fewer than "
                    f"{MIN_BINS_TO_FIT} is not enough to correlate against, so no "
                    f"reach was fitted. Lower bin_size for a short target.",
                    *caveats,
                ]
            ),
        )

    observed = _binned(depth, starts, widths)

    binned_by_reach: Dict[int, np.ndarray] = {
        reach: _binned(predicted_depth(positions, length, reach, circular=circular), starts, widths)
        for reach in reaches
    }
    scores = {reach: _spearman(binned, observed) for reach, binned in binned_by_reach.items()}

    best_reach = max(scores, key=lambda r: scores[r])
    best = scores[best_reach]
    informative = best >= MIN_INFORMATIVE_CORRELATION

    # Everything the data cannot separate from the winner, not the winner
    # alone: a flat optimum and a sharp one must not read as the same answer.
    plausible = tuple(sorted(r for r, rho in scores.items() if rho >= best - PLAUSIBLE_MARGIN))
    at_grid_edge = best_reach == max(reaches)
    cv = _blocked_cv_correlation(binned_by_reach, observed)

    if not informative:
        note = (
            f"Best correlation {best:.3f} is below {MIN_INFORMATIVE_CORRELATION}: the "
            f"primer positions do not explain this depth profile, so the fitted reach "
            f"is not evidence. Check that the BAM comes from an SWGA reaction using "
            f"these primers, mapped to this target."
        )
    else:
        note = (
            f"Fitted reach {best_reach:,} bp (in-sample Spearman {best:.3f} over "
            f"{n_bins:,} bins"
            + (f", held out {cv:.3f}" if cv is not None else "")
            + "). This is a model parameter fitted to one dataset, not a measured "
            "constant: use it as the centre of a sensitivity analysis over "
            + (
                f"{min(plausible):,}-{max(plausible):,} bp, the range this data " "cannot separate."
                if len(plausible) > 1
                else "the grid."
            )
        )
        if at_grid_edge:
            note += (
                f" Reach {best_reach:,} bp is the largest on the grid, so the true "
                f"value may be higher; widen the grid to bound it."
            )

    return ReachFit(
        best_reach=best_reach,
        correlation=best,
        by_reach=scores,
        informative=informative,
        contig=contig,
        note=" ".join([note, *caveats]),
        cv_correlation=cv,
        plausible_reaches=plausible,
        at_grid_edge=at_grid_edge,
        n_bins=n_bins,
    )


def format_reach_table(fit: ReachFit) -> str:
    """Render the fit so the shape of the optimum is visible, not just its location.

    A flat optimum means the data cannot distinguish 3 kb from 20 kb, which is a
    different situation from a sharp one and should not be reported as the same
    number with the same confidence.
    """
    lines = [f"{'reach (bp)':>12}{'Spearman rho':>16}", "-" * 28]
    for reach in sorted(fit.by_reach):
        if reach == fit.best_reach:
            marker = "  <-- best"
        elif reach in fit.plausible_reaches:
            marker = "  (within margin)"
        else:
            marker = ""
        lines.append(f"{reach:>12,}{fit.by_reach[reach]:>16.3f}{marker}")
    lines.append("")
    # Labelled, because the unqualified word invited the figure to be read as
    # validation when it is the maximum over the grid on the data that chose it.
    lines.append(f"In-sample rho (optimistic, maximum over the grid): {fit.correlation:.3f}")
    if fit.cv_correlation is not None:
        lines.append(f"Held-out rho (blocked cross-validation):          {fit.cv_correlation:.3f}")
    if fit.n_bins:
        lines.append(f"Bins correlated:                                 {fit.n_bins:,}")
    lines.append("")
    lines.append(fit.note)
    if fit.informative:
        lines.append("")
        lines.append(f"Use with: neoswga optimize -j params.json --coverage-reach {fit.best_reach}")
        if len(fit.plausible_reaches) > 1:
            lines.append(
                f"and repeat at {min(fit.plausible_reaches):,} and "
                f"{max(fit.plausible_reaches):,} bp: this data does not separate them."
            )
    return "\n".join(lines)
