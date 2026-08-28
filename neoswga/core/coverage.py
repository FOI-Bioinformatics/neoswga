"""Centralized coverage computation for optimizer post-processing and
ad-hoc rescoring CLIs.

The single `compute_per_prefix_coverage` helper consolidates the numpy-
vectorised union-of-extension-windows coverage calculation that previously
lived inline in `cli_unified.run_contract_set` and `run_rescore_set`
(Phase 12D). Extracted here so `unified_optimizer` (Phase 15A) and the
CLIs share the same code path.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


def compute_per_prefix_coverage(
    cache,
    primers: Sequence[str],
    prefixes: Sequence[str],
    seq_lengths: Sequence[int],
    extension: int = 3000,
    strand: str = "both",
    circular: bool = False,
) -> Tuple[float, Dict[str, float]]:
    """Compute union-of-extension-windows coverage across prefixes.

    For each (primer, prefix) pair, fetch positions from the PositionCache
    and mark occupied[start:end] = True on a bool numpy array; coverage
    fraction is the sum of the occupied array divided by the genome length.

    Args:
        cache: PositionCache instance that responds to
            ``get_positions(prefix, primer, strand)``.
        primers: primer sequences to query.
        prefixes: HDF5 file prefixes (usually fg_prefixes or bg_prefixes).
        seq_lengths: genome lengths matching `prefixes` elementwise.
        extension: extension reach in bp. Default 3000 bp corresponds to
            the effective per-primer reach in a dense phi29 SWGA design
            (Clarke et al. 2017 post-hoc <5 kb inter-primer-site filter;
            Dwivedi-Yu et al. 2023 1/2-5 kbp successful-set densities).
            A position
            at ``p`` marks
            ``[max(0, p - extension), min(length, p + extension))``
            as occupied. Use :func:`polymerase_extension_reach` to pick
            a per-polymerase value; pass 70000 explicitly if you want
            the theoretical phi29 processivity upper bound (a different
            question — "could two primers connect in principle?").
        strand: 'both' / 'forward' / 'reverse'.
        circular: If True, windows that run off either end of the
            sequence wrap around to the other end. Default False matches
            linear chromosomes, fragmented assemblies, and concatenated
            multi-genome inputs. Set True for single closed-circular
            bacterial chromosomes or plasmids (passes `fg_circular` /
            `bg_circular` from params.json through to here).

    Returns:
        (aggregate_coverage, per_prefix_coverage_dict). aggregate is the
        ratio of total-bases-covered to total-bases across every prefix;
        per_prefix maps each prefix to its own coverage fraction.

        Empty inputs or an unusable cache yield ``(0.0, {})``.
    """
    if not prefixes or not seq_lengths or cache is None or not primers:
        return 0.0, {}

    per_prefix: Dict[str, float] = {}
    total_cov = 0
    total_len = 0

    for prefix, length in zip(prefixes, seq_lengths):
        if length <= 0:
            per_prefix[prefix] = 0.0
            continue
        occupied = np.zeros(length, dtype=bool)
        # If the window diameter equals or exceeds the genome on a
        # circular target, every position is covered by any single site.
        if circular and 2 * extension >= length:
            occupied[:] = True
        else:
            for primer in primers:
                try:
                    positions = cache.get_positions(prefix, primer, strand)
                except (KeyError, ValueError):
                    # Genuinely-absent primer in this prefix; skip it. Do NOT
                    # swallow I/O / HDF5 errors here — a systemic cache failure
                    # must surface rather than silently undercount coverage.
                    continue
                for pos in positions:
                    _mark_window(occupied, int(pos), extension, length, circular)
        covered = int(occupied.sum())
        per_prefix[prefix] = covered / length if length else 0.0
        total_cov += covered
        total_len += length

    agg = total_cov / total_len if total_len else 0.0
    return agg, per_prefix


def _mark_window(
    occupied: "np.ndarray",
    pos: int,
    extension: int,
    length: int,
    circular: bool,
) -> None:
    """Mark ``occupied[pos-extension:pos+extension]`` as True.

    On circular sequences the window wraps across the origin when it
    runs past either end. On linear sequences the window is clipped.
    Used by :func:`compute_per_prefix_coverage` and
    :meth:`base_optimizer.BaseOptimizer._compute_coverage` to keep the
    two coverage implementations in sync.
    """
    start = pos - extension
    end = pos + extension
    if circular:
        if start < 0 and end > length:
            # Window exceeds genome from both ends; everything is covered.
            occupied[:] = True
        elif start < 0:
            occupied[0:end] = True
            occupied[length + start : length] = True
        elif end > length:
            occupied[start:length] = True
            occupied[0 : end - length] = True
        else:
            occupied[start:end] = True
    else:
        clipped_start = max(0, start)
        clipped_end = min(length, end)
        if clipped_end > clipped_start:
            occupied[clipped_start:clipped_end] = True


def polymerase_extension_reach(
    polymerase: str,
    default: int = 3000,
    coverage_metric: str = "realistic",
) -> int:
    """Resolve the extension reach (in bp) for a polymerase.

    Phase 16 critical gap #2: distinguish two legitimate meanings of
    "extension reach":

    - ``coverage_metric='realistic'`` (default): returns the effective
      per-primer reach in a dense SWGA design (phi29 ~3 kb, equiphi29
      ~4 kb, bst ~1 kb, klenow ~1.5 kb). In dense multi-primer reactions
      extension is truncated by neighbouring primers' strand-displacement
      products after a few kb — Clarke et al. (2017) filter candidate
      sets to <5 kb mean inter-primer-site spacing; Dwivedi-Yu et al.
      (2023) report successful Prevotella sets at 1/2-5 kbp densities.
      Use this for `fg_coverage` / `per_target_coverage`.
    - ``coverage_metric='processivity'``: returns the theoretical
      single-molecule processivity (phi29 70 kb, equiphi29 80 kb,
      bst 2 kb, klenow 10 kb). Use this when the question is graph-level
      reachability — "can primer A and B connect via one extension event
      in principle?" — not "how much genome does the set cover?".

    Prior to Phase 16 this helper returned processivity unconditionally,
    which inflated `fg_coverage` 5-20x over what the selected primer set
    can actually amplify in a dense SWGA reaction. Callers that still
    need the legacy behaviour pass ``coverage_metric='processivity'``
    explicitly.

    Args:
        polymerase: Polymerase name (phi29 / equiphi29 / bst / klenow).
        default: Fallback when the polymerase is unknown or the helper
            is unavailable. Default switched to 3000 bp (phi29 per-primer
            reach) in Phase 16; pass 70000 explicitly if the legacy
            processivity value is intended.
        coverage_metric: 'realistic' (default, Phase 16+) or 'processivity'.

    Returns:
        Extension reach in bp.
    """
    try:
        if coverage_metric == "processivity":
            from .reaction_conditions import get_polymerase_processivity

            return int(get_polymerase_processivity(polymerase))
        # realistic (default)
        from .reaction_conditions import get_typical_amplicon_length

        return int(get_typical_amplicon_length(polymerase))
    except Exception:
        return default


# Reaches at which coverage is reported alongside the one used for selection.
#
# Coverage is meaningless without the reach it was computed at, and the
# published tools do not agree on one: swga 2.0's `coverage_ratio` uses phi29's
# ~70 kb single-molecule processivity, while NeoSWGA selects on the ~3 kb
# realistic per-primer reach. The same primer set scores 0.418 at 3 kb, 0.836 at
# 10 kb and ~1.0 at 70 kb, so a bare "coverage" number cannot be compared across
# tools at all. Reporting the curve rather than one point makes the comparison
# possible without either convention having to be restated.
REPORTING_REACHES = (3000, 10000, 70000)


def product_reach(polymerase, default=10000):
    """Measured MDA product length -- the reach headline coverage is reported at.

    Distinct from the *selection* reach (`resolve_coverage_reach`), and the two
    are deliberately different because they answer different questions.

    Selection uses the design density that published successful sets have (2-5
    kb inter-primer spacing; Clarke 2017, Dwivedi-Yu 2023). Selecting at the
    physical reach instead would stop the greedy once the genome was covered at
    ~10 kb, producing designs sparser than any set with measured wet-lab
    success.

    Reporting uses what is actually amplified: unselective phi29 MDA product
    peaks near 10 kb (Dean et al. 2002 PNAS 99:5261; Picher et al. 2016 Nat
    Commun 7:13296). Reporting at the selection reach would understate coverage,
    since strand displacement means a downstream primer does not truncate
    extension -- it gets displaced -- so spacing does not bound product length.
    """
    try:
        from .registry import views as _views

        entry = _views.as_characteristics().get(str(polymerase).lower())
        if entry:
            return int(entry.get("product_length") or entry["typical_amplicon_length"])
    except Exception:  # pragma: no cover - registry always present in practice
        pass
    return default


def resolve_coverage_reach(polymerase, override=None, default=3000):
    """The reach used for set-cover selection and reported `fg_coverage`.

    `override` is `params.json`'s `coverage_reach` or `--coverage-reach`. It
    exists because the right value is an empirical property of the reaction --
    how far a primer's product actually extends before a neighbour's strand
    displacement truncates it -- and the literature the default rests on
    reports inter-primer *spacing* chosen by designers, not measured reach.
    `neoswga calibrate-reach --bam` estimates it from sequencing depth.

    Without an override this returns the polymerase's realistic per-primer
    reach, unchanged.
    """
    if override is not None:
        if not isinstance(override, (int, float)) or isinstance(override, bool):
            raise ValueError(
                f"coverage_reach must be a positive number of base pairs, got {override!r}"
            )
        if override <= 0:
            raise ValueError(f"coverage_reach must be positive, got {override}")
        return int(override)
    return polymerase_extension_reach(polymerase, default=default, coverage_metric="realistic")


# Benchmark reference points, from the two published set-level datasets with
# wet-lab outcomes. See docs/validation/published_primer_sets.md.
#
# Dwivedi-Yu et al. (2023) PLOS Comput Biol 19:e1010137 - six Prevotella sets.
# Clarke et al. (2017) Bioinformatics 33:2071 - twelve M. tuberculosis sets.
#
# The two disagree about WHICH statistic predicts success, and that disagreement
# is the reason these are reported side by side rather than folded into one
# score. Clarke pre-filtered their primer pool for even binding, which removed
# most of the variance in evenness; density then separated their winners
# cleanly. Dwivedi-Yu did not, so evenness varied and the worst coverage hole
# separated theirs instead. Whichever property is currently limiting is the one
# that predicts, so the user needs to see all three.
BENCHMARK_MEAN_GAP_BP = (2000, 5000)  # successful published sets, both datasets
# Boundary drawn between the two successful sets (0.533, 0.552) and the four
# that failed (0.566, 0.600, 0.623, 0.638) in Dwivedi-Yu 2023. It is a cut read
# off six labelled points, NOT a measured constant - treat it as a rough
# orientation marker rather than a threshold with meaning. It also does not
# transfer: evenness did not separate winners in Clarke 2017 at all.
BENCHMARK_GINI_GOOD = 0.56


def _unknown_gap_rows():
    """Rows used when the gap statistics are absent or not finite."""
    return [
        {"label": label, "value": "unknown", "verdict": "unknown", "detail": detail}
        for label, detail in (
            ("Mean gap", "No gaps measured."),
            ("Max gap", "No gaps measured; nothing to compare against primer reach."),
            ("Gap evenness (Gini)", "No gaps measured."),
        )
    ]


def interpret_gap_metrics(mean_gap, max_gap, gap_gini, extension_reach=3000):
    """Describe gap statistics against the published benchmarks.

    Returns a list of (label, value, verdict) tuples for display. The verdicts
    reference real datasets rather than invented thresholds, and deliberately
    avoid collapsing to a single pass/fail: the two benchmarks disagree about
    which statistic matters, so the honest output shows where the set sits on
    each and lets the reader judge.

    Args:
        mean_gap: Mean distance between adjacent binding sites (bp).
        max_gap: Largest gap between binding sites (bp).
        gap_gini: Gini coefficient of the gap distribution (0 = perfectly even).
        extension_reach: Realistic per-primer reach (bp). A gap wider than twice
            this contains sequence no primer can reach.

    Returns:
        List of dicts with keys: label, value, verdict, detail.
    """

    # A failed optimization carries inf gaps (see base_optimizer.json_safe), and
    # a summary JSON round-tripped through `null` brings them back as None.
    # Neither is a measurement, so report them as such rather than rendering
    # "About inf kb of this gap sits beyond the reach of any flanking primer",
    # which read as a real figure in a real report.
    def _finite(x):
        return isinstance(x, (int, float)) and math.isfinite(x)

    if not (_finite(mean_gap) and _finite(max_gap) and _finite(gap_gini)):
        return _unknown_gap_rows()

    out = []

    lo, hi = BENCHMARK_MEAN_GAP_BP
    if mean_gap <= 0:
        verdict, detail = "unknown", "No gaps measured."
    elif mean_gap <= hi:
        verdict = "within published range"
        detail = (
            f"Successful published sets run {lo/1000:.0f}-{hi/1000:.0f} kb mean "
            f"spacing (Clarke 2017; Dwivedi-Yu 2023)."
        )
    else:
        verdict = "sparser than published sets"
        detail = (
            f"Successful published sets run {lo/1000:.0f}-{hi/1000:.0f} kb. "
            f"Density separated the winners in Clarke 2017, though not in "
            f"Dwivedi-Yu 2023."
        )
    out.append(
        {
            "label": "Mean gap",
            "value": f"{mean_gap/1000:.1f} kb",
            "verdict": verdict,
            "detail": detail,
        }
    )

    reachable = 2 * max(1, extension_reach)
    if max_gap <= 0:
        verdict, detail = "unknown", "No gaps measured."
    elif max_gap <= reachable:
        verdict = "fully reachable"
        detail = (
            f"No stretch exceeds twice the {extension_reach/1000:.0f} kb per-primer "
            f"reach, so no region is out of reach of some primer."
        )
    else:
        unreachable = max_gap - reachable
        verdict = "leaves an unreachable stretch"
        detail = (
            f"About {unreachable/1000:.1f} kb of this gap sits beyond the reach of "
            f"any flanking primer. The worst gap was the strongest predictor of "
            f"measured enrichment in Dwivedi-Yu 2023 (rho -0.83); a set can look "
            f"dense on average and still fail on one hole."
        )
    out.append(
        {
            "label": "Max gap",
            "value": f"{max_gap/1000:.1f} kb",
            "verdict": verdict,
            "detail": detail,
        }
    )

    # Gini of 0 is legitimate when there is a single gap (or perfectly regular
    # spacing) - it means "perfectly even", not "not measured". Only treat it as
    # unknown when there were no gaps to measure at all.
    if mean_gap <= 0 and max_gap <= 0:
        verdict, detail = "unknown", "No gaps measured."
    elif gap_gini <= BENCHMARK_GINI_GOOD:
        verdict = "even"
        detail = (
            f"At or below {BENCHMARK_GINI_GOOD}, the range of the two successful "
            f"sets in Dwivedi-Yu 2023 (0.533 and 0.552). That cut is drawn from "
            f"six labelled sets, so treat it as orientation rather than a "
            f"threshold."
        )
    else:
        verdict = "uneven"
        detail = (
            f"Above {BENCHMARK_GINI_GOOD}. The four sets that failed in "
            f"Dwivedi-Yu 2023 ran 0.566-0.638, the two that worked 0.533-0.552. "
            f"That cut comes from six labelled sets and does not transfer: "
            f"evenness did NOT separate winners in Clarke 2017, whose pool was "
            f"pre-filtered for even binding."
        )
    out.append(
        {
            "label": "Gap evenness (Gini)",
            "value": f"{gap_gini:.3f}",
            "verdict": verdict,
            "detail": detail,
        }
    )

    return out


def gap_regime_note(gap_gini):
    """One-line guidance on which statistic is likely limiting.

    Included because the composite `normalized_score` does not use `max_gap` at
    all, so a reader relying on the score alone cannot see this distinction.
    """
    if gap_gini > BENCHMARK_GINI_GOOD:
        return (
            "Binding is uneven, so coverage holes are the likely limiting factor - "
            "look at max gap before adding primers to raise average density."
        )
    return (
        "Binding is reasonably even, so density is the likely lever - compare "
        "mean gap against the 2-5 kb range of published successful sets."
    )
