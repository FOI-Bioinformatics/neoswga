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


def _has_homopolymer_run(primer: str) -> bool:
    """Check if primer contains homopolymer runs exceeding threshold."""
    for pattern in HOMOPOLYMER_PATTERNS:
        if pattern in primer:
            return True
    return False


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

    # Rule 3: GC clamp (last 5 bases)
    # Adaptive based on genome GC content to avoid eliminating valid primers
    # for AT-rich (e.g. Plasmodium ~25% GC) or GC-rich (e.g. Mycobacterium ~65% GC) targets
    gc_in_last_5 = _count_gc(primer[-5:])
    genome_gc = getattr(parameter, "genome_gc", None)
    if genome_gc is not None and genome_gc < EXTREME_AT_GENOME_GC:
        # AT-rich genome: allow 0 GC in last 5, reject >3
        gc_clamp_min = 0
        gc_clamp_max = MAX_GC_IN_LAST_5_BASES
    elif genome_gc is not None and genome_gc > EXTREME_GC_GENOME_GC:
        # GC-rich genome: allow up to 4 GC in last 5, require >=1 AT
        gc_clamp_min = 1
        gc_clamp_max = 4
    else:
        # Standard: 1-3 GC in last 5
        gc_clamp_min = 1
        gc_clamp_max = MAX_GC_IN_LAST_5_BASES
    if gc_in_last_5 > gc_clamp_max or gc_in_last_5 < gc_clamp_min:
        logger.debug(
            f"GC clamp filter: {primer} GC_last5={gc_in_last_5} "
            f"(allowed {gc_clamp_min}-{gc_clamp_max})"
        )
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

    primer_to_fg_count = get_rates_for_one_species(primer_list, fg_prefixes)

    # Check if Bloom filter should be used for background filtering
    # Auto-enable if bg_bloom is specified in params (common user config)
    bloom_path = getattr(parameter, "bloom_filter_path", None) or getattr(
        parameter, "bg_bloom", None
    )
    use_bloom = getattr(parameter, "use_bloom_filter", bloom_path is not None)

    if use_bloom and bloom_path:
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
