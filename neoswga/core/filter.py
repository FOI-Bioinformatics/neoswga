"""
Primer filtering module for NeoSWGA.

Implements sequence-based filtering rules to select high-quality primer candidates.
"""

import logging
import multiprocessing
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from neoswga.core import dimer, parameter, primer_attributes
from neoswga.core.exceptions import InvalidDesignRequest, ReferenceDataError
from neoswga.core.parameter import (
    EXTREME_AT_GENOME_GC,
    EXTREME_GC_GENOME_GC,
    default_reaction_temp,
    default_tm_range,
)
from neoswga.core.reaction_conditions import ReactionConditions, build_reaction_conditions

logger = logging.getLogger(__name__)

# Expected number of sampled sites at the background threshold below which
# the sampled index cannot distinguish 'rare' from 'absent'. Five is a floor,
# not a precision claim: below it the estimate is dominated by whether a
# single position happened to land in the sample.
_MIN_RESOLVABLE_SAMPLED_SITES = 5.0

# Module-level reaction conditions (lazily initialized)
_reaction_conditions = None


def _get_reaction_conditions() -> ReactionConditions:
    """
    Get ReactionConditions from parameter settings (cached).

    Creates ReactionConditions once from global parameters to avoid
    repeated object creation during filtering.
    """
    global _reaction_conditions
    if _reaction_conditions is None:
        # Hand-listing the fields here dropped nine of them, including every
        # buffer species and the glycerol/PEG/BSA/SSB group the filter CLI
        # accepts. `build_reaction_conditions` reads the field list off the
        # constructor, so nothing can be left out by omission. It also keeps
        # the reaction_temp fallback this function needed: the global EXISTS
        # as None until get_params runs, so a getattr default never fires and
        # ReactionConditions(temp=None) raises when it range-checks.
        _reaction_conditions = build_reaction_conditions()
    return _reaction_conditions


def reset_reaction_conditions():
    """Reset cached reaction conditions (for testing or parameter changes)."""
    global _reaction_conditions
    _reaction_conditions = None


# =============================================================================
# Filtering Constants
# =============================================================================

# Rule 5: Maximum homopolymer run length (consecutive identical bases)
MAX_HOMOPOLYMER_RUN = 5

# Rule 4: Dinucleotide repeat patterns to reject (10bp = 5 repeats of dinucleotide)
DINUCLEOTIDE_REPEAT_PATTERNS = (
    "ATATATATAT",
    "TATATATATA",  # AT repeats
    "AGAGAGAGAG",
    "GAGAGAGAGA",  # AG repeats
    "ACACACACAC",
    "CACACACACA",  # AC repeats
    "TCTCTCTCTC",
    "CTCTCTCTCT",  # TC repeats
    "GTGTGTGTGT",
    "TGTGTGTGTG",  # GT repeats
    "CGCGCGCGCG",
    "GCGCGCGCGC",  # CG repeats
)

# Rule 4: Minimum primer length to check for dinucleotide repeats
MIN_LENGTH_FOR_DINUCLEOTIDE_CHECK = 10

# Rule 4: Minimum count of a nucleotide to consider for repeat patterns
MIN_NUCLEOTIDE_COUNT_FOR_REPEAT = 5

# Rule 3: Maximum G/C bases allowed in last 5 bases (GC clamp)
MAX_GC_IN_LAST_5_BASES = 3
# Bases at the 3' end the clamp looks at. See resolve_gc_clamp for why this
# is not scaled to primer length.
GC_CLAMP_WINDOW = 5

# Rule 1: All 3 bases at 3' end cannot be G/C
MAX_GC_AT_3PRIME_END = 2  # Max 2 of 3 bases can be G/C

# Homopolymer patterns (generated from MAX_HOMOPOLYMER_RUN)
HOMOPOLYMER_PATTERNS = tuple(base * MAX_HOMOPOLYMER_RUN for base in "ACGT")


def _scale_freq_threshold(
    base_threshold: float, primer_length: int, reference_k: int = 10
) -> float:
    """
    Scale a frequency threshold inversely with primer length.

    Longer primers have exponentially fewer exact match sites in a genome.
    A 15bp primer typically has 5-20 sites in a 5 Mbp genome, while a 10bp
    primer may have thousands. This function adjusts the threshold so that
    long primers are not eliminated by thresholds calibrated for short ones.

    That argument is about the foreground FLOOR (`min_fg_freq`) and is right:
    an unscaled 1e-5 floor demands 31.7 sites on a 3.2 Mb target, which almost
    no 12-mer has, so without this every long primer would be rejected.

    It is also applied to the background CEILING (`max_bg_freq`), where the
    argument does not carry: the harm from a background site does not depend on
    primer length. A 12-mer with ten host sites primes the host ten times,
    exactly as a 10-mer with ten does. Scaling the ceiling therefore reinterprets
    the parameter at every length -- on chr21, 5e-7 admits 23.35 sites at k=10,
    1.46 at k=12 and 0.09 at k=14, so from k=13 the gate means "zero exact
    background matches" and stops responding to the parameter.

    Measured cost on the Prevotella 12-mer pool: of 573,182 candidates the
    configured ceiling admits, the scaled one keeps 143,224 -- 75% discarded,
    and not for quality (median 2 foreground sites either way).

    Whether the ceiling should scale is a genuine trade-off, since admitting
    primers with more host sites raises coverage and lowers selectivity, and it
    is deliberately not settled here. What is not optional is that the applied
    threshold be visible: see `describe_freq_gate` and `log_freq_gates`.

    Args:
        base_threshold: Frequency threshold calibrated for reference_k-mers
        primer_length: Actual primer length
        reference_k: Reference primer length (default: 10)

    Returns:
        Scaled threshold
    """
    if primer_length <= reference_k:
        return base_threshold
    length_diff = reference_k - primer_length
    scale_factor = 4.0**length_diff  # e.g. 4^(-5) = 1/1024 for 15bp
    return base_threshold * scale_factor


def describe_freq_gate(base_threshold: float, primer_length: int, genome_length: int) -> dict:
    """What a frequency gate actually admits, in sites.

    A frequency is not a quantity anyone can reason about against a params.json
    value; a site count is. `degenerate` marks a ceiling that has fallen below
    one site, where the gate means "zero exact matches" and no value of the
    parameter in its documented range loosens it.
    """
    scaled = _scale_freq_threshold(base_threshold, primer_length)
    sites = scaled * genome_length
    return {
        "configured": base_threshold,
        "scaled": scaled,
        "sites": sites,
        "degenerate": sites < 1.0,
    }


def log_freq_gates(
    primer_length: int,
    min_fg_freq: float,
    max_bg_freq: float,
    fg_length: int,
    bg_length: int,
) -> None:
    """Report the thresholds a run is really applying, once per primer length.

    The length scaling was invisible: nothing printed the applied value, so a
    user comparing params.json against the result had no way to see that
    `max_bg_freq: 5e-07` had become 3.13e-08 at k=12.
    """
    fg = describe_freq_gate(min_fg_freq, primer_length, fg_length)
    bg = describe_freq_gate(max_bg_freq, primer_length, bg_length) if bg_length else None

    message = (
        f"k={primer_length} frequency gates: min_fg_freq {fg['configured']:.3g} "
        f"-> {fg['scaled']:.3g} (>= {fg['sites']:.2f} foreground sites)"
    )
    if bg is not None:
        message += (
            f"; max_bg_freq {bg['configured']:.3g} -> {bg['scaled']:.3g} "
            f"(< {bg['sites']:.2f} background sites)"
        )
    logger.info(message)

    if bg is not None and bg["degenerate"]:
        logger.warning(
            f"At k={primer_length} the background gate admits only primers with "
            f"zero exact background matches ({bg['sites']:.2f} sites). "
            f"max_bg_freq={bg['configured']:.3g} is scaled by 4^(10-{primer_length}) "
            f"and no value in its documented range loosens this. Raise it by that "
            f"factor if you meant to allow background binding at this length."
        )


def _warn_if_sample_too_sparse(sampled_index) -> None:
    """Say so when the sampled index cannot resolve the gate's threshold.

    `SampledGenomeIndex` stores every `sample_rate`-th position and extrapolates,
    so it can only resolve counts that survive being divided by that rate. The
    count at the background gate's threshold is `max_bg_freq * genome_size`:

        human, 3 Gbp at 5e-6   ->  15000 sites, sampled ~150 times   resolvable
        plasmid, 6 kb at 5e-6  ->   0.03 sites, sampled 0.0003 times  hopeless

    Below resolution every primer estimates zero and passes, which is the same
    silent under-filtering the old sentinel produced, reached by a different
    route. Warned rather than raised: unlike a missing index the numbers here
    are real, only too coarse, and a caller may know their background better
    than this heuristic. Exact k-mer counting is the answer at these sizes, and
    it is affordable precisely because the genome is small.
    """
    if getattr(sampled_index, "source", None) == "kmer_counts":
        # Exact counts read off a jellyfish table. Nothing was sampled, so
        # there is no sparsity to assess: `sample_rate` is 1 and the
        # extrapolation in `estimate_count` is the identity. Without this the
        # check was skipped only because that route leaves `genome_size` at 0,
        # which is an accident rather than a rule.
        return

    genome_size = getattr(sampled_index, "genome_size", 0) or 0
    sample_rate = getattr(sampled_index, "sample_rate", 1) or 1
    if genome_size <= 0:
        return

    lengths = getattr(parameter, "bg_seq_lengths", None) or [genome_size]
    max_bg_freq = getattr(parameter, "max_bg_freq", None)
    if not max_bg_freq:
        return

    threshold_sites = max_bg_freq * sum(lengths)
    expected_sampled = threshold_sites / sample_rate

    if expected_sampled < _MIN_RESOLVABLE_SAMPLED_SITES:
        logger.warning(
            f"Sampled index is too sparse to resolve the background threshold: "
            f"max_bg_freq={max_bg_freq:g} over {sum(lengths):,} bp is "
            f"{threshold_sites:.3g} sites, which sample_rate={sample_rate} sees "
            f"{expected_sampled:.3g} times. Background counts will read near "
            f"zero and most primers will pass. Use exact k-mer counting for a "
            f"background this size, or rebuild the index with a lower sample rate."
        )


def _load_sampled_index(bloom_path: str):
    """The sampled index that turns Bloom presence into a usable count.

    A Bloom filter answers "is this k-mer in the background", not "how often".
    The background gate downstream is a frequency test
    (`bg_count / bg_total_length < max_bg_freq`), so this path used to
    manufacture a stand-in count -- `bloom_max_bg_matches + 1`, 11 by default --
    for every primer the filter reported present.

    An absolute count fed to a frequency test inverts with scale. 11 over a
    6 kb plasmid is 1.8e-3 and gets rejected; 11 over a 3 Gbp human background
    is 3.7e-9 and passes any sane threshold. Since an absent primer scores 0,
    which also passes, BOTH branches passed on a large background: the filter
    did nothing whatsoever, silently, at exactly the scale Bloom exists for.
    Every small-scale test of it passed, which is why it survived.

    So the stand-in is gone. `SampledGenomeIndex` can answer the frequency
    question -- `build-filter` already writes `bg_sampled.pkl` beside
    `bg_bloom.pkl` -- and it is looked for there when `sampled_index_path` is
    unset, because requiring a second undocumented setting is what put callers
    on the sentinel branch to begin with. If no index can be found this raises,
    which costs a user one clear error instead of a primer set screened against
    nothing.
    """
    import os

    from neoswga.core.background_filter import SampledGenomeIndex

    configured = getattr(parameter, "sampled_index_path", None)
    candidates = [configured] if configured else []

    # Where `build-filter` puts it (cli/pipeline.py, background_filter.py:678).
    directory = os.path.dirname(os.path.abspath(bloom_path))
    candidates.append(os.path.join(directory, "bg_sampled.pkl"))
    # And the matching name for a bloom file the user renamed.
    stem = os.path.basename(bloom_path)
    if "bloom" in stem:
        candidates.append(os.path.join(directory, stem.replace("bloom", "sampled")))

    for path in candidates:
        if path and os.path.exists(path):
            try:
                index = SampledGenomeIndex.load(path)
                logger.info(f"Using sampled index for background counts: {path}")
                return index
            except Exception as e:
                logger.warning(f"Could not load sampled index {path}: {e}")

    raise ValueError(
        f"No sampled index found for Bloom background {bloom_path!r}. A Bloom "
        f"filter reports presence, not frequency, and the background gate is a "
        f"frequency test -- without counts this path cannot screen anything. "
        f"Re-run 'neoswga build-filter' (it writes bg_sampled.pkl next to "
        f"bg_bloom.pkl), or set 'sampled_index_path' in params.json, or drop "
        f"the Bloom filter to use exact k-mer counts."
    )


def _reject_lengths_the_filter_cannot_answer(bloom, primer_list, bloom_path) -> None:
    """A length the filter never indexed reads as absent, not as unmeasured.

    `contains` answers False for a k-mer of a length that was never inserted.
    The background count is then zero, and zero clears any frequency gate, so
    a filter built over k 6-12 screens a 13-mer design by passing all of it.
    The filter is not at fault: absence is the honest answer to a question
    outside the domain it was built over. What was missing was any record of
    that domain and any check against it.

    This is the silent-zero family again: not a scan that found nothing, an
    integer that saturated, a cache asked for what it does not hold or a
    dictionary default that flatters, but a query outside a model's domain
    answered as though it were inside.

    A filter carrying no recorded range predates the record. That is UNKNOWN
    rather than wrong, so it warns and proceeds -- the rule `digest_algorithm`
    established for a k-mer table written before its provenance sidecar.
    """
    if bloom.min_k is None or bloom.max_k is None:
        logger.warning(
            "Bloom filter %s records no k-mer range, so the primer lengths it "
            "covers cannot be checked. A length it does not hold reads as "
            "absent from the background and clears the gate unscreened. "
            "Rebuild with 'neoswga build-filter' to record the range.",
            bloom_path,
        )
        return

    lengths = sorted({len(p) for p in primer_list})
    outside = [k for k in lengths if not bloom.covers_length(k)]
    if not outside:
        return

    raise ReferenceDataError(
        artifact=bloom_path,
        reason=(
            f"the filter indexes k-mers of length {bloom.min_k}-{bloom.max_k}, "
            f"but this design uses {', '.join(str(k) for k in outside)}. A "
            f"length the filter does not hold reads as absent from the "
            f"background, so every such primer would clear the background gate "
            f"unscreened"
        ),
        remediation=(
            f"rebuild it over the design's range with 'neoswga build-filter "
            f"--genome <background> -o <dir> --min-k {min(lengths)} "
            f"--max-k {max(lengths)}'"
        ),
    )


def get_bg_rates_via_bloom(primer_list: List[str], bloom_path: str) -> Dict[str, int]:
    """
    Get background rates using a pre-built Bloom filter.

    Memory-efficient alternative to loading full k-mer files for large genomes
    like human genome (3 Gbp). Uses O(1) lookup per primer instead of loading
    entire k-mer dictionary into memory.

    Args:
        primer_list: List of primer sequences to check
        bloom_path: Path to pre-built Bloom filter (.pkl file)

    Returns:
        Dictionary mapping primer -> estimated count (or None if not found)
    """
    if len(primer_list) == 0:
        logger.info("Bloom filter: no primers to check")
        return {}

    from neoswga.core.background_filter import BackgroundBloomFilter

    logger.info(f"Using Bloom filter for background filtering: {bloom_path}")

    # Load bloom filter
    bloom = BackgroundBloomFilter.load(bloom_path)
    _reject_lengths_the_filter_cannot_answer(bloom, primer_list, bloom_path)
    sampled_index = _load_sampled_index(bloom_path)
    _warn_if_sample_too_sparse(sampled_index)

    primer_to_count = {}

    total = len(primer_list)
    bloom_hits = 0

    for primer in primer_list:
        if bloom.contains(primer):
            bloom_hits += 1
            primer_to_count[primer] = sampled_index.estimate_count(primer)
        else:
            # Definitely not in background: a Bloom filter has no false
            # negatives, so absence is trustworthy and zero is an honest
            # answer. Presence is not -- it carries no count, which is why the
            # sampled index above is required rather than optional.
            primer_to_count[primer] = 0

    logger.info(
        f"Bloom filter results: {bloom_hits}/{total} primers found in background "
        f"({100*bloom_hits/total:.1f}%)"
    )

    return primer_to_count


def _resolve_tm_window() -> Tuple[float, float]:
    """The Tm window this filter applies, in Celsius.

    Two faults lived in the two lines this replaces:

        tm_min = getattr(parameter, "min_tm", None) or 15
        tm_max = getattr(parameter, "max_tm", None) or 55

    The fallback was a fixed 15-55 whatever the enzyme, so a params.json that
    did not mention Tm screened a bst design at 63 C through a window built for
    phi29 -- keeping primers that cannot prime at that temperature and
    rejecting ones that can. And `or` treats a configured 0.0 as absent, the
    same sentinel-versus-value confusion found elsewhere in this codebase; 0 C
    is a legitimate way to ask for no lower bound.

    `get_params` now resolves both from the polymerase when params.json is
    silent, so the globals are normally populated. The fallback here is for
    library callers that reach the filter without going through it, and it
    reads the same registry rather than being a fourth independent answer to
    "what Tm is acceptable".
    """
    polymerase = getattr(parameter, "polymerase", None) or "phi29"
    default_low, default_high = default_tm_range(polymerase)

    tm_min = getattr(parameter, "min_tm", None)
    tm_max = getattr(parameter, "max_tm", None)
    return (
        default_low if tm_min is None else tm_min,
        default_high if tm_max is None else tm_max,
    )


def _count_gc(sequence: str) -> int:
    """Count G and C bases in a sequence."""
    return sum(1 for base in sequence if base in "GC")


def _configured_int(name: str, default: int, allow_zero: bool = False) -> int:
    """An integer setting from `parameter`, or the default.

    Type-checked rather than truthiness-checked. `neoswga.core.filter.parameter`
    is replaced with a MagicMock in several test modules, and a MagicMock
    answers every attribute with a truthy stub whose `int()` is 1 -- so a plain
    `getattr(...) or default` silently resolved the GC clamp to a 1-base window
    with a 1-base limit. The same guard rejects a string or None arriving from
    a hand-edited params.json, where the failure would be a quietly different
    filter rather than an error.
    """
    value = getattr(parameter, name, None)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return default
    value = int(value)
    if value < 0 or (value == 0 and not allow_zero):
        return default
    return value


def resolve_homopolymer_run() -> int:
    """Longest run of one base a primer may contain, from configuration.

    Was a module constant, with the rejection patterns derived from it once at
    import -- so a configured value would have had no effect on them. Read per
    call now; the default reproduces the hardcoded 5.
    """
    return _configured_int("max_homopolymer_run", MAX_HOMOPOLYMER_RUN)


def resolve_gc_clamp() -> Tuple[int, int]:
    """(window, max G/C) for the 3'-end clamp, from configuration.

    The window is length-blind: five bases is 83% of a 6-mer and 42% of a
    12-mer, and the rule comes from PCR design on 18-25mers. Scaling it to
    primer length was measured and rejected -- against the only benchmark with
    per-primer amplification, a scaled window newly admits primers amplifying
    at a median of 1.55x where the pool median is 10.15x. The fixed window
    earns its place at short lengths by acting as the extra constraint that a
    redundancy analysis makes it look like it is not. Configurable so that
    conclusion can be re-tested, with the default reproducing it.
    """
    return (
        _configured_int("gc_clamp_window", GC_CLAMP_WINDOW),
        _configured_int("max_gc_in_clamp", MAX_GC_IN_LAST_5_BASES, allow_zero=True),
    )


def _fails_gc_clamp(primer: str) -> bool:
    """Whether the 3'-end GC clamp rejects this primer.

    Carries the genome-GC adaptation with it: an AT-rich target (~25% GC)
    cannot supply a G/C in every primer's 3' end, and a GC-rich one (~65%)
    cannot avoid it, so a fixed band would reject nearly everything for both.
    """
    window, configured_max = resolve_gc_clamp()
    observed = _count_gc(primer[-window:])

    genome_gc = getattr(parameter, "genome_gc", None)
    if genome_gc is not None and genome_gc < EXTREME_AT_GENOME_GC:
        low, high = 0, configured_max
    elif genome_gc is not None and genome_gc > EXTREME_GC_GENOME_GC:
        low, high = 1, max(configured_max + 1, 4)
    else:
        low, high = 1, configured_max

    if observed > high or observed < low:
        logger.debug(f"GC clamp filter: {primer} GC_last{window}={observed} (allowed {low}-{high})")
        return True
    return False


def _has_homopolymer_run(primer: str) -> bool:
    """Check if primer contains homopolymer runs exceeding threshold."""
    run = resolve_homopolymer_run()
    return any(base * run in primer for base in "ACGT")


def _has_dinucleotide_repeats(primer: str, nucleotide_counts: Counter[str]) -> bool:
    """
    Check if primer contains problematic dinucleotide repeats.

    Only checks if primer is long enough and has sufficient counts of
    the nucleotides involved in each repeat pattern.
    """
    if len(primer) < MIN_LENGTH_FOR_DINUCLEOTIDE_CHECK:
        return False

    # Find nucleotides with high counts (potential repeat participants)
    high_count_nucleotides = {
        nucleo
        for nucleo, count in nucleotide_counts.items()
        if count >= MIN_NUCLEOTIDE_COUNT_FOR_REPEAT
    }

    # Only check if at least 2 nucleotides have high counts
    if len(high_count_nucleotides) < 2:
        return False

    # Check each pattern
    for pattern in DINUCLEOTIDE_REPEAT_PATTERNS:
        if pattern in primer:
            return True

    return False


def filter_extra(primer: str) -> bool:
    """
    Filter primer based on sequence quality rules.

    Applies five filtering rules to ensure primer quality:
    1. No 3 consecutive G/C at 3' end (prevents mispriming)
    2. GC content within acceptable range (40-60% default, adaptive for extreme genomes)
    3. GC clamp: 1-3 G/C in last 5 bases (promotes specific 3' binding)
    4. No excessive dinucleotide repeats (prevents mispriming)
    5. No long homopolymer runs (prevents mispriming)

    Args:
        primer: DNA sequence to evaluate (5' to 3' direction)

    Returns:
        True if primer passes all filters, False otherwise
    """
    # Tm filtering with reaction condition corrections
    # Use effective Tm that accounts for additives (DMSO, betaine, etc.)
    conditions = _get_reaction_conditions()
    primer_tm = conditions.calculate_effective_tm(primer)
    tm_min, tm_max = _resolve_tm_window()
    if not (tm_min <= primer_tm <= tm_max):
        logger.debug(
            f"Tm filter: {primer} effective Tm={primer_tm:.1f} outside [{tm_min}, {tm_max}]"
        )
        return False

    # Rule 5: Check for homopolymer runs
    if _has_homopolymer_run(primer):
        logger.debug(f"Homopolymer filter: {primer}")
        return False

    # Calculate nucleotide counts (used for multiple rules)
    nucleotide_counts = Counter(primer)
    gc_count = nucleotide_counts.get("G", 0) + nucleotide_counts.get("C", 0)

    # Rule 2: GC content filtering
    gc_content = gc_count / len(primer)
    if not (parameter.gc_min <= gc_content <= parameter.gc_max):
        logger.debug(f"GC content filter: {primer} GC={gc_content:.1%}")
        return False

    # Rule 3: GC clamp at the 3' end, adapted to genome GC.
    if _fails_gc_clamp(primer):
        return False

    # Rule 1: 3' end cannot have all 3 bases as G/C
    gc_in_last_3 = _count_gc(primer[-3:])
    if gc_in_last_3 > MAX_GC_AT_3PRIME_END:
        logger.debug(f"3' end GC filter: {primer} GC_last3={gc_in_last_3}")
        return False

    # Rule 4: Check for dinucleotide repeats
    if _has_dinucleotide_repeats(primer, nucleotide_counts):
        logger.debug(f"Dinucleotide repeat filter: {primer}")
        return False

    # Self-dimer check
    if dimer.is_dimer_fast(primer, primer, parameter.max_self_dimer_bp):
        logger.debug(f"Self-dimer filter: {primer}")
        return False

    return True


def _resolve_background_source() -> Tuple[bool, Optional[str]]:
    """Which background source this run screens against, or refuse.

    Resolved before any counting: this is a configuration error, and
    reading the foreground k-mer tables first only delays it behind file
    I/O that the run is about to discard.

    Returns (use_bloom, bloom_path).
    """
    # Both keys are required, and that is not obvious: a user who builds a
    # filter and sets only the path gets exact counting and no explanation,
    # having paid the build cost -- hours on a host genome.
    #
    # The comment here used to say "auto-enable if bg_bloom is specified", and
    # neither half was true. `use_bloom_filter` is a module global that always
    # exists and defaults to False, so the `getattr` fallback that would have
    # enabled it could never be reached; and `bg_bloom` is not a schema key and
    # is assigned nowhere, so that arm of the `or` was permanently dead.
    bloom_path = getattr(parameter, "bloom_filter_path", None)
    use_bloom = getattr(parameter, "use_bloom_filter", False)

    if use_bloom and not bloom_path:
        # The damaging pairing, and it was silent. This flag is documented to
        # leave `bg_prefixes` empty -- parameter.py sets `bg_seq_lengths = []`
        # for exactly that case -- so the exact-counting fallback below has no
        # k-mer files to read. Every background count then comes back absent,
        # and an absent count PASSES the gate. That rule is deliberate, since a
        # k-mer missing from a jellyfish table may still have sites the string
        # search finds, but here it is not one k-mer missing: it is the whole
        # background gate off, on a run the user asked to screen a host with.
        raise InvalidDesignRequest(
            field="bloom_filter_path",
            reason=(
                "use_bloom_filter is true but no Bloom filter path is set. The "
                "flag also leaves bg_prefixes empty, so the exact-counting "
                "fallback has no k-mer files to read: every background count "
                "would be absent, and an absent count passes the gate, which "
                "would screen nothing at all. Set bloom_filter_path to a filter "
                "built with 'neoswga build-filter', or set use_bloom_filter to "
                "false to count background k-mers exactly"
            ),
        )

    if bloom_path and not use_bloom:
        logger.warning(
            'bloom_filter_path is set to %s but "use_bloom_filter" is false, '
            "so the Bloom filter is NOT being used and background k-mers are "
            "counted exactly. Set use_bloom_filter to true to use it.",
            bloom_path,
        )

    return use_bloom, bloom_path


def get_all_rates(
    primer_list: List[str],
    fg_prefixes: List[str],
    bg_prefixes: List[str],
    fg_total_length: int,
    bg_total_length: int,
) -> pd.DataFrame:
    """
    Computes the foreground and background binding site frequencies normalized by their respective genome lengths.

    Args:
        primer_list: The list of primers to compute frequencies for.
        fg_prefixes: The list of foreground path prefixes used for creating the kmer files.
        bg_prefixes: The list of background path prefixes used for creating the kmer files.
        fg_total_length: The total number of base pairs in the foregound genome.
        bg_total_length: The total number of base pairs in the background genome.

    Returns:
        df: A pandas dataframe with the sequence, unnormalized counts, and  columns fg_bool and bg_bool which indicate if the sequence passes the respective filters.
    """

    use_bloom, bloom_path = _resolve_background_source()
    primer_to_fg_count = get_rates_for_one_species(primer_list, fg_prefixes)

    if use_bloom:
        primer_to_bg_count = get_bg_rates_via_bloom(primer_list, bloom_path)
    else:
        primer_to_bg_count = get_rates_for_one_species(primer_list, bg_prefixes)

    results = []

    # When a k-mer is absent from the count dictionary (count is None), it passes
    # the frequency filter. This is intentional: k-mers missing from the jellyfish
    # output may still have binding sites found by downstream string search. The
    # subsequent Gini index and scoring steps provide additional filtering, so
    # retaining these candidates at this stage avoids premature exclusion.
    #
    # Frequency thresholds are scaled by primer length: longer primers have
    # exponentially fewer exact match sites, so fixed thresholds calibrated
    # for short primers would eliminate nearly all long (15-18bp) candidates.
    # Cache frequency thresholds by primer length to avoid redundant computation
    _threshold_cache = {}
    for primer in primer_list:
        primer_len = len(primer)
        if primer_len not in _threshold_cache:
            _threshold_cache[primer_len] = (
                _scale_freq_threshold(parameter.min_fg_freq, primer_len),
                _scale_freq_threshold(parameter.max_bg_freq, primer_len),
            )
            # Once per length, say what is actually being applied. The scaling
            # silently divides both gates by 4^(k-10), and at k>=13 the
            # background ceiling falls below one site.
            log_freq_gates(
                primer_length=primer_len,
                min_fg_freq=parameter.min_fg_freq,
                max_bg_freq=parameter.max_bg_freq,
                fg_length=fg_total_length,
                bg_length=bg_total_length,
            )
        scaled_min_fg, scaled_max_bg = _threshold_cache[primer_len]
        fg_count = primer_to_fg_count.get(primer, None)
        # Guard against zero-length genomes (e.g. a foreground FASTA that failed
        # to load to a real length) so a bad input raises clearly upstream
        # rather than a per-primer ZeroDivisionError. A count with no length is
        # uninterpretable; treat it as failing the frequency gate.
        if fg_count is None or fg_total_length <= 0:
            fg_bool = fg_count is None
        else:
            fg_bool = fg_count / fg_total_length > scaled_min_fg
        bg_count = primer_to_bg_count.get(primer, None)
        if bg_count is None or bg_total_length <= 0:
            bg_bool = bg_count is None
        else:
            bg_bool = bg_count / bg_total_length < scaled_max_bg
        results.append([primer, fg_count, bg_count, fg_bool, bg_bool])

    df = pd.DataFrame(results, columns=["primer", "fg_count", "bg_count", "fg_bool", "bg_bool"])

    return df


def get_rates_for_one_species(primer_list: List[str], fname_prefixes: List[str]) -> Dict[str, int]:
    """
    Computes the binding site frequencies for all ppsth prefixes in fname_prefixes.

    Args:
        primer_list: The list of primers to compute frequencies for.
        fg_prefixes: The list of foreground path prefixes used for creating the kmer files.

    Returns:
        all_primer_to_count: A dictonary of primer to frequency.
    """
    stratified_primer_list = {}

    for primer in primer_list:
        k = len(primer)
        if k not in stratified_primer_list:
            stratified_primer_list[k] = []
        stratified_primer_list[k].append(primer)

    tasks = []

    for fname_prefix in fname_prefixes:
        for k, primer_list_k in stratified_primer_list.items():
            tasks.append((primer_list_k, fname_prefix, k))

    # Use ThreadPoolExecutor for I/O-bound file reads (avoids process creation
    # overhead and serialization costs compared to multiprocessing.Pool)
    with ThreadPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
        results = list(executor.map(_get_rate_for_one_file, tasks))

    all_primer_to_count = {}

    for primer_to_count in results:
        for primer, count in primer_to_count.items():
            if primer not in all_primer_to_count:
                all_primer_to_count[primer] = count
            else:
                all_primer_to_count[primer] += count
    return all_primer_to_count


def _get_rate_for_one_file(task: Tuple[List[str], str, int]) -> Dict[str, int]:
    primer_list, fname_prefix, k = task
    primer_set = set(primer_list)
    primer_to_count = {}
    found = 0
    target = len(primer_set)
    with open(fname_prefix + "_" + str(k) + "mer_all.txt", "r") as f_in:
        for line in f_in:
            parts = line.split()
            if parts[0] in primer_set:
                primer_to_count[parts[0]] = int(parts[1])
                found += 1
                if found == target:
                    break

    return {primer: primer_to_count.get(primer, 0) for primer in primer_list}


def check_gini_stage_kept_something(before_df, after_df):
    """Refuse an empty candidate pool, and name the reason it is empty.

    The Gini index is NaN where evenness cannot be measured, and
    `filter.get_gini`'s `.notna()` guard drops those rows. On a small target
    that can be most of the pool: across the plasmid example's frequency
    survivors, 92 of 10,532 position-file keys carry three or more combined
    sites.

    Writing an empty step2_df.csv defers the failure to step 3, which reports
    it as a missing 'primer' column. Raising here names the parameter to change,
    the same way `apply_qa_filter_to_step2_file` does for the QA pass. It
    lives here rather than in `pipeline` because the rows it is about are the
    ones `get_gini` below drops with its `.notna()` guard; `pipeline` imports
    it under the same name. Since
    `min_gini_sites` became configurable the first suggestion is to lower it,
    because that is the gate that emptied the pool and the user can now change
    it without editing the source.
    """
    if len(before_df) == 0 or len(after_df) > 0:
        return

    configured = getattr(parameter, "min_gini_sites", None)
    threshold = (
        int(configured)
        if isinstance(configured, int) and not isinstance(configured, bool) and configured > 0
        else primer_attributes.DEFAULT_MIN_GINI_SITES
    )
    raise ValueError(
        f"The evenness (Gini) stage removed all {len(before_df)} candidates. "
        f"Every one of them binds the target fewer than min_gini_sites "
        f"({threshold}) times across both strands, so their spacing is not "
        f"measurable.\n"
        f"On a small target this is expected. Either lower the threshold, with "
        f'"min_gini_sites": 1 in params.json or --min-gini-sites 1 on '
        f"neoswga filter, which accepts single-site primers and their "
        f"unmeasurable evenness; or admit more abundant k-mers by lowering "
        f"min_k, raising max_bg_freq, or lowering min_fg_freq."
    )


def get_gini(
    fg_prefixes: List[str],
    fg_genomes: List[str],
    fg_seq_lengths: List[int],
    df: pd.DataFrame,
    circular: bool,
    position_cache: Optional[Dict] = None,
) -> pd.DataFrame:
    """Computes the Gini index of the gap distances between binding sites.

    Args:
        fg_prefixes: List of path prefixes to the kmer files of the foreground genome.
        fg_genomes: List of paths to the foreground fasta files.
        fg_seq_lengths: List of sequence length(s) of the foreground genome(s).
        df: Pandas dataframe with column primer containing the primer sequences.
        circular: Whether the genome is circular.
        position_cache: Optional dict mapping (prefix, primer) -> positions
            from string_search.get_positions(). When provided, avoids the
            HDF5 read-back round-trip for a measurable speedup.

    Returns:
        df: Input dataframe with new column 'gini' for the computed Gini indices.

    """
    df["gini"] = primer_attributes.get_gini_from_txt(
        df["primer"].values,
        fg_prefixes,
        fg_genomes,
        fg_seq_lengths,
        circular,
        position_cache=position_cache,
    )

    if len(df["gini"]) == 0:
        df["gini_bool"] = []
        return df

    # Vectorized boolean operation (10-50x faster than apply with lambda)
    df["gini_bool"] = df["gini"].notna() & (df["gini"] < parameter.max_gini)

    return df[df["gini_bool"]]
