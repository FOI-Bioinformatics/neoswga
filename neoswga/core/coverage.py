"""Centralized coverage computation for optimizer post-processing and
ad-hoc rescoring CLIs.

The single `compute_per_prefix_coverage` helper consolidates the numpy-
vectorised union-of-extension-windows coverage calculation that previously
lived inline in `cli_unified.run_contract_set` and `run_rescore_set`
(Phase 12D). Extracted here so `unified_optimizer` (Phase 15A) and the
CLIs share the same code path.
"""

from __future__ import annotations

import bisect
import math
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .exceptions import ModelEvaluationError, ReferenceDataError, UnsupportedModelError


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

    for prefix, length in zip(prefixes, seq_lengths, strict=True):
        if length <= 0:
            per_prefix[prefix] = 0.0
            continue
        occupied = np.zeros(length, dtype=bool)
        found: List[int] = []
        for primer in primers:
            try:
                positions = cache.get_positions(prefix, primer, strand)
            except (KeyError, ValueError):
                # Genuinely-absent primer in this prefix; skip it. Do NOT
                # swallow I/O / HDF5 errors here — a systemic cache failure
                # must surface rather than silently undercount coverage.
                continue
            found.extend(int(pos) for pos in positions)
        starts = _record_starts_for(cache, prefix)

        # If the window diameter equals or exceeds the genome on a circular
        # target, every position is covered by ANY SINGLE SITE -- which is why
        # the premise has to be checked. Without `found`, an empty cache took
        # this branch and reported 1.0 for a panel that binds nothing.
        if circular and 2 * extension >= length and found:
            occupied[:] = True
        else:
            for pos in found:
                _mark_window(occupied, pos, extension, length, circular, record_starts=starts)
        covered = int(occupied.sum())
        per_prefix[prefix] = covered / length if length else 0.0
        total_cov += covered
        total_len += length

    agg = total_cov / total_len if total_len else 0.0
    return agg, per_prefix


def marginal_coverage_curve(
    cache,
    primers: Sequence[str],
    prefixes: Sequence[str],
    seq_lengths: Sequence[int],
    extension: int = 3000,
    strand: str = "both",
    circular: bool = False,
) -> List[Dict[str, float]]:
    """Cumulative coverage as the delivered primers are added, one at a time.

    Set size is the most consequential choice a user makes and nothing showed
    the shape of the return curve. `--auto-size` estimates a size from a
    closed-form coverage model without reading the candidate pool, and
    `--show-frontier` stops at 20 primers, so a design at 96 or 160 had no way
    to see where the gain flattened short of running a sweep by hand.

    This is measured on ONE delivered set, in its delivered order: entry k is
    the coverage of that set's first k primers. It is not a sweep. A prefix of
    a 160-primer set is not the set an optimizer would return if asked for 160,
    so each value is a lower bound on a re-optimization at that size. Read the
    curve for its shape, not its level.

    Cost is one pass over the positions plus one array sum per primer, on the
    foreground genome only.

    Args:
        cache: PositionCache instance responding to
            ``get_positions(prefix, primer, strand)``.
        primers: primer sequences, in the order they should be added.
        prefixes: HDF5 file prefixes, usually fg_prefixes.
        seq_lengths: genome lengths matching `prefixes` elementwise.
        extension: per-primer reach in bp. Pass the same value the run was
            selected and scored on -- :func:`resolve_coverage_reach` -- or the
            curve will not end at the reported `fg_coverage`.
        strand: 'both' / 'forward' / 'reverse'.
        circular: wrap windows that run off either end.

    Returns:
        One dict per primer with keys ``n`` (1-based count), ``primer``,
        ``coverage`` (cumulative fraction over all prefixes) and
        ``marginal_pp`` (the gain that primer added, in percentage points).
        Empty inputs yield an empty list.
    """
    if not prefixes or not seq_lengths or cache is None or not primers:
        return []

    usable = [(p, int(n)) for p, n in zip(prefixes, seq_lengths, strict=True) if int(n) > 0]
    if not usable:
        return []

    occupied = {prefix: np.zeros(length, dtype=bool) for prefix, length in usable}
    total_len = sum(length for _, length in usable)

    curve: List[Dict[str, float]] = []
    previous = 0.0
    for index, primer in enumerate(primers, start=1):
        for prefix, length in usable:
            marks = occupied[prefix]
            try:
                positions = cache.get_positions(prefix, primer, strand)
            except (KeyError, ValueError):
                # Genuinely-absent primer in this prefix; skip it. Do NOT
                # swallow I/O / HDF5 errors here.
                continue
            sites = [int(pos) for pos in positions]
            if not sites:
                continue
            # Same premise as compute_per_prefix_coverage. The check used to sit
            # before the lookup and inside the per-primer loop, so the FIRST
            # primer marked the whole genome whether or not it bound anything,
            # and every later entry inherited that.
            if circular and 2 * extension >= length:
                marks[:] = True
                continue
            starts = _record_starts_for(cache, prefix)
            for pos in sites:
                _mark_window(marks, pos, extension, length, circular, record_starts=starts)

        covered = sum(int(occupied[prefix].sum()) for prefix, _ in usable)
        coverage = covered / total_len if total_len else 0.0
        curve.append(
            {
                "n": index,
                "primer": primer,
                "coverage": coverage,
                "marginal_pp": (coverage - previous) * 100.0,
            }
        )
        previous = coverage

    return curve


def _record_starts_for(cache, prefix):
    """Record offsets for a prefix, when the cache can supply them.

    `getattr` rather than a hard call: these helpers accept any cache-shaped
    object, several tests pass a stub, and an index written before record
    starts were stored has none. In both cases the window is not confined and
    behaviour is what it was.
    """
    getter = getattr(cache, "get_record_starts", None)
    if getter is None:
        # No getter at all: an index written before record starts were stored,
        # or one of the cache-shaped stubs the geometry helpers accept. The
        # window is not confined to a record and behaviour is what it was.
        return None
    try:
        return getter(prefix) or None
    except Exception as exc:
        # A getter that exists and FAILS is a different event. Losing record
        # boundaries lets a coverage window run across a record edge, so the
        # panel is credited with covering bases on a contig its site is not on.
        # That inflates coverage silently and in proportion to how fragmented
        # the reference is.
        raise ReferenceDataError(
            f"record starts for prefix '{prefix}'",
            str(exc),
            "Re-run `neoswga filter` to rebuild the position index.",
        ) from exc


def _mark_window(
    occupied: "np.ndarray",
    pos: int,
    extension: int,
    length: int,
    circular: bool,
    record_starts: Optional[Sequence[int]] = None,
) -> None:
    """Mark ``occupied[pos-extension:pos+extension]`` as True.

    On circular sequences the window wraps across the origin when it
    runs past either end. On linear sequences the window is clipped.
    Used by :func:`compute_per_prefix_coverage` and
    :meth:`base_optimizer.BaseOptimizer._compute_coverage` to keep the
    two coverage implementations in sync.

    ``record_starts`` gives the offsets at which each FASTA record begins in
    the concatenated sequence, including 0. When supplied, the window is
    confined to the record holding ``pos``: a polymerase extending from a site
    near the end of one contig does not continue into the next one, and the
    only reason it appeared to was that concatenation left no record boundary
    for the window to stop at (audit F3).

    Wrapping is suppressed when records are known, because a multi-record file
    is not one circle -- joining its last record to its first would be the same
    error in another direction. A single-record circular genome still wraps,
    since its one record spans the whole sequence.
    """
    start = pos - extension
    end = pos + extension

    if record_starts:
        index = bisect.bisect_right(record_starts, pos) - 1
        lo = record_starts[index] if index >= 0 else 0
        hi = record_starts[index + 1] if index + 1 < len(record_starts) else length
        single_record = len(record_starts) == 1 and lo == 0 and hi == length
        if not single_record:
            clipped_start = max(lo, start)
            clipped_end = min(hi, end)
            if clipped_end > clipped_start:
                occupied[clipped_start:clipped_end] = True
            return

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


def merged_window_intervals(
    positions: Sequence[int],
    extension: int,
    length: int,
    circular: bool,
) -> List[Tuple[int, int]]:
    """The union of one primer's binding windows, as disjoint half-open spans.

    The same geometry :func:`_mark_window` marks into a boolean array, returned
    as intervals instead. A caller that needs the union of several primers'
    windows weighted differently can then accumulate over segment boundaries
    rather than over bases, which is what
    :meth:`base_optimizer.BaseOptimizer._compute_effective_coverage` does: its
    cost stops depending on the length of the genome and starts depending on
    the number of sites, which is smaller by orders of magnitude.

    Kept next to ``_mark_window`` and tested against it, because two
    descriptions of where a window falls is how this codebase has produced
    disagreeing coverage numbers before.

    Spans are sorted, non-overlapping, and clipped to ``[0, length)``. Adjacent
    spans are merged, since ``[0, 5)`` and ``[5, 9)`` cover the same bases as
    ``[0, 9)``. ``record_starts`` is deliberately not a parameter: neither
    ``_union_coverage`` nor ``_compute_effective_coverage`` passes record starts
    to ``_mark_window`` today, so accepting them here would let a caller
    silently change what those two report.
    """
    if length <= 0 or extension < 0:
        return []

    spans: List[Tuple[int, int]] = []
    for raw in positions:
        pos = int(raw)
        start = pos - extension
        end = pos + extension
        if circular:
            if start < 0 and end > length:
                # The window laps the whole molecule; nothing else can add to it.
                return [(0, length)]
            if start < 0:
                if end > 0:
                    spans.append((0, min(end, length)))
                spans.append((max(0, length + start), length))
                continue
            if end > length:
                spans.append((min(start, length), length))
                if end - length > 0:
                    spans.append((0, min(end - length, length)))
                continue
        clipped_start = max(0, start)
        clipped_end = min(length, end)
        if clipped_end > clipped_start:
            spans.append((clipped_start, clipped_end))

    if not spans:
        return []

    spans.sort()
    merged = [list(spans[0])]
    for start, end in spans[1:]:
        if start <= merged[-1][1]:
            if end > merged[-1][1]:
                merged[-1][1] = end
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


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
      ~4 kb, bst ~1 kb, klenow ~1.5 kb). Use this for `fg_coverage` /
      `per_target_coverage`.

      **Where this figure comes from.** It is a DESIGN-DENSITY convention
      taken from sets with measured wet-lab success -- Clarke et al. (2017)
      filter candidates to <5 kb mean inter-primer-site spacing, Dwivedi-Yu
      et al. (2023) report successful Prevotella sets at 1/2-5 kbp densities.
      It is not a measured extension distribution, and no reading of it
      identifies a physical processivity.

      This docstring previously explained it as extension being "truncated by
      neighbouring primers' strand-displacement products". That mechanism is
      wrong, and `product_reach` in this same module says so: phi29 displaces
      a downstream product rather than being stopped by it, so spacing does
      not bound product length. The two statements contradicted each other;
      corrected 2026-09-14 (audit F2) in favour of the mechanism, leaving the
      empirical spacing convention as the actual justification for 3 kb.

      Coverage is highly sensitive to this choice. On one saved 26-oligo wMel
      panel under identical conditions the same design reads 41.3% at 1 kb,
      80.1% at 3 kb, 93.5% at 5 kb and 99.6% at 10 kb. A single figure quoted
      without its reach carries almost no information.
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
    if coverage_metric not in ("processivity", "realistic"):
        raise UnsupportedModelError("coverage reach", coverage_metric, "realistic | processivity")
    try:
        if coverage_metric == "processivity":
            from .reaction_conditions import get_polymerase_processivity

            value = get_polymerase_processivity(polymerase)
        else:
            from .reaction_conditions import get_typical_amplicon_length

            value = get_typical_amplicon_length(polymerase)
    except (KeyError, LookupError, ValueError) as exc:
        # An unknown polymerase used to return `default`. Coverage scales
        # directly with reach -- the same wMel panel reads 41.3% at 1 kb and
        # 93.5% at 5 kb -- so a substituted reach silently rescales every
        # coverage figure in the report and the user has no way to see it.
        raise UnsupportedModelError(
            "polymerase reach", polymerase, "phi29 | equiphi29 | bst | klenow"
        ) from exc
    except Exception as exc:
        raise ModelEvaluationError("polymerase reach", polymerase, str(exc)) from exc

    if value is None:
        raise UnsupportedModelError(
            "polymerase reach", polymerase, "phi29 | equiphi29 | bst | klenow"
        )
    return int(value)


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
    except Exception as exc:
        raise ReferenceDataError(
            "polymerase registry", str(exc), "Reinstall the package; the registry ships with it."
        ) from exc

    if not entry:
        # This reach is the one headline coverage is REPORTED at, so an
        # unrecognised polymerase silently substituting 10 kb changes the
        # headline number rather than any internal one.
        raise UnsupportedModelError(
            "product reach", polymerase, ", ".join(sorted(_views.as_characteristics()))
        )
    return int(entry.get("product_length") or entry["typical_amplicon_length"])


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


def resolve_extension_reach(kwargs, verbose: bool = False) -> int:
    """The reach a run computes coverage at, resolved from kwargs and params.

    Wraps `resolve_coverage_reach` with the precedence an optimizer run needs:
    an explicit kwarg, then the `parameter` global set from params.json or
    `--coverage-reach`, then the polymerase's realistic per-primer reach.
    Extracted from `unified_optimizer.run_optimization` so the resolution sits
    beside the function it defers to.
    """
    import logging

    from . import parameter

    polymerase = kwargs.get("polymerase") or getattr(parameter, "polymerase", "phi29")
    override = kwargs.get("coverage_reach")
    if override is None:
        override = getattr(parameter, "coverage_reach", None)
    reach = resolve_coverage_reach(polymerase, override=override)
    if override is not None and verbose:
        logging.getLogger(__name__).info(
            f"Coverage reach: {reach:,} bp (explicit; polymerase default would be "
            f"used otherwise). Coverage figures are not comparable across "
            f"different reaches."
        )
    return reach


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
