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
"""

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

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


@dataclass(frozen=True)
class ReachFit:
    """Result of fitting reach to observed depth."""

    best_reach: int
    correlation: float
    by_reach: Dict[int, float]
    informative: bool
    contig: str
    note: str

    def to_dict(self) -> Dict:
        return {
            "best_reach": self.best_reach,
            "correlation": self.correlation,
            "by_reach": {str(k): v for k, v in self.by_reach.items()},
            "informative": self.informative,
            "contig": self.contig,
            "note": self.note,
        }


def predicted_depth(positions: Sequence[int], length: int, reach: int) -> np.ndarray:
    """Expected relative depth if each binding site contributes over +/- reach.

    A triangular kernel rather than a rectangular one: product length is a
    distribution, not a hard cutoff, so depth should taper with distance from
    the site rather than stopping. The shape matters less than the scale -- the
    fit is choosing between reaches an order of magnitude apart -- but a
    rectangular kernel makes every reach predict a flat profile once sites
    overlap, which removes the signal the fit runs on.
    """
    profile = np.zeros(length, dtype=np.float64)
    if reach <= 0:
        return profile
    # One triangle per site, accumulated. Vectorised per site rather than per
    # base: a 3 Mb genome with 500 sites is 500 slice adds, not 1.5e9 updates.
    ramp = 1.0 - np.abs(np.arange(-reach, reach + 1, dtype=np.float64)) / reach
    for pos in positions:
        lo, hi = int(pos) - reach, int(pos) + reach + 1
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


def fit_reach(
    depth: np.ndarray,
    positions: Sequence[int],
    contig: str = "",
    reaches: Optional[Sequence[int]] = None,
    bin_size: int = 1000,
) -> ReachFit:
    """Choose the reach whose predicted depth profile best matches the observed one.

    Ranked by Spearman correlation, which cares about the shape of the profile
    rather than its scale -- absolute depth depends on sequencing effort, which
    says nothing about reach.

    Binning to `bin_size` before correlating: per-base depth is dominated by
    mapping noise at a scale far below any candidate reach, and correlating raw
    bases mostly measures that noise.
    """
    from scipy.stats import spearmanr

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

    n_bins = max(1, length // bin_size)
    usable = n_bins * bin_size
    observed = depth[:usable].reshape(n_bins, bin_size).mean(axis=1)

    scores: Dict[int, float] = {}
    for reach in reaches:
        predicted = predicted_depth(positions, length, reach)
        binned = predicted[:usable].reshape(n_bins, bin_size).mean(axis=1)
        if np.allclose(binned, binned[0]):
            # Every base predicted equally -- the reach is so large relative to
            # the genome that the profile is flat and carries no information.
            scores[reach] = 0.0
            continue
        rho = spearmanr(binned, observed).statistic
        scores[reach] = 0.0 if np.isnan(rho) else float(rho)

    best_reach = max(scores, key=lambda r: scores[r])
    best = scores[best_reach]
    informative = best >= MIN_INFORMATIVE_CORRELATION

    if not informative:
        note = (
            f"Best correlation {best:.3f} is below {MIN_INFORMATIVE_CORRELATION}: the "
            f"primer positions do not explain this depth profile, so the fitted reach "
            f"is not evidence. Check that the BAM comes from an SWGA reaction using "
            f"these primers, mapped to this target."
        )
    elif best_reach == max(scores):
        note = (
            f"Best reach {best_reach:,} bp is the largest on the grid, so the true "
            f"value may be higher; widen the grid to bound it."
        )
    else:
        note = f"Fitted reach {best_reach:,} bp (Spearman {best:.3f} against observed depth)."

    return ReachFit(
        best_reach=best_reach,
        correlation=best,
        by_reach=scores,
        informative=informative,
        contig=contig,
        note=note,
    )


def format_reach_table(fit: ReachFit) -> str:
    """Render the fit so the shape of the optimum is visible, not just its location.

    A flat optimum means the data cannot distinguish 3 kb from 20 kb, which is a
    different situation from a sharp one and should not be reported as the same
    number with the same confidence.
    """
    lines = [f"{'reach (bp)':>12}{'Spearman rho':>16}", "-" * 28]
    for reach in sorted(fit.by_reach):
        marker = "  <-- best" if reach == fit.best_reach else ""
        lines.append(f"{reach:>12,}{fit.by_reach[reach]:>16.3f}{marker}")
    lines.append("")
    lines.append(fit.note)
    if fit.informative:
        lines.append("")
        lines.append(f"Use with: neoswga optimize -j params.json --coverage-reach {fit.best_reach}")
    return "\n".join(lines)
