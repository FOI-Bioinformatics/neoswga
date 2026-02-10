"""
Primer filtering module for NeoSWGA.

Implements sequence-based filtering rules to select high-quality primer candidates.
"""

from collections import Counter
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import melting
import multiprocessing

import neoswga.core.parameter as parameter
import neoswga.core.primer_attributes as primer_attributes
import neoswga.core.dimer as dimer
from neoswga.core.reaction_conditions import ReactionConditions
import logging

logger = logging.getLogger(__name__)

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
        _reaction_conditions = ReactionConditions(
            temp=getattr(parameter, 'reaction_temp', 30.0),
            polymerase=getattr(parameter, 'polymerase', 'phi29'),
            dmso_percent=getattr(parameter, 'dmso_percent', 0.0),
            betaine_m=getattr(parameter, 'betaine_m', 0.0),
            trehalose_m=getattr(parameter, 'trehalose_m', 0.0),
            formamide_percent=getattr(parameter, 'formamide_percent', 0.0),
            ethanol_percent=getattr(parameter, 'ethanol_percent', 0.0),
            urea_m=getattr(parameter, 'urea_m', 0.0),
            tmac_m=getattr(parameter, 'tmac_m', 0.0),
            na_conc=getattr(parameter, 'na_conc', 50.0),
            mg_conc=getattr(parameter, 'mg_conc', 0.0),
        )
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
    'ATATATATAT', 'TATATATATA',  # AT repeats
    'AGAGAGAGAG', 'GAGAGAGAGA',  # AG repeats
    'ACACACACAC', 'CACACACACA',  # AC repeats
    'TCTCTCTCTC', 'CTCTCTCTCT',  # TC repeats
    'GTGTGTGTGT', 'TGTGTGTGTG',  # GT repeats
    'CGCGCGCGCG', 'GCGCGCGCGC',  # CG repeats
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
HOMOPOLYMER_PATTERNS = tuple(base * MAX_HOMOPOLYMER_RUN for base in 'ACGT')


def _scale_freq_threshold(base_threshold: float, primer_length: int,
                          reference_k: int = 10) -> float:
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
    scale_factor = 4.0 ** length_diff  # e.g. 4^(-5) = 1/1024 for 15bp
    return base_threshold * scale_factor


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
    from neoswga.core.background_filter import BackgroundBloomFilter, BackgroundFilter

    logger.info(f"Using Bloom filter for background filtering: {bloom_path}")

    # Load bloom filter
    bloom = BackgroundBloomFilter.load(bloom_path)

    # Optional: load sampled index for count estimation
    sampled_path = getattr(parameter, 'sampled_index_path', None)
    sampled_index = None
    if sampled_path:
        from neoswga.core.background_filter import SampledGenomeIndex
        try:
            sampled_index = SampledGenomeIndex.load(sampled_path)
            logger.info(f"Using sampled index for count estimation: {sampled_path}")
        except Exception as e:
            logger.warning(f"Failed to load sampled index: {e}")

    max_bg_matches = getattr(parameter, 'bloom_max_bg_matches', 10)
    primer_to_count = {}

    total = len(primer_list)
    bloom_hits = 0

    for primer in primer_list:
        if bloom.contains(primer):
            bloom_hits += 1
            # If we have sampled index, use it for count estimation
            if sampled_index:
                estimated_count = sampled_index.estimate_count(primer)
                primer_to_count[primer] = estimated_count
            else:
                # Without sampled index, use max_bg_matches as indicator
                # (primer is present, use threshold value)
                primer_to_count[primer] = max_bg_matches + 1
        else:
            # Definitely not in background (no false negatives in Bloom filter)
            primer_to_count[primer] = 0

    logger.info(f"Bloom filter results: {bloom_hits}/{total} primers found in background "
                f"({100*bloom_hits/total:.1f}%)")

    return primer_to_count


def _count_gc(sequence: str) -> int:
    """Count G and C bases in a sequence."""
    return sum(1 for base in sequence if base in 'GC')


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
        nucleo for nucleo, count in nucleotide_counts.items()
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
    tm_min = getattr(parameter, 'min_tm', 15)
    tm_max = getattr(parameter, 'max_tm', 55)
    if not (tm_min <= primer_tm <= tm_max):
        logger.debug(f"Tm filter: {primer} effective Tm={primer_tm:.1f} outside [{tm_min}, {tm_max}]")
        return False

    # Rule 5: Check for homopolymer runs
    if _has_homopolymer_run(primer):
        logger.debug(f"Homopolymer filter: {primer}")
        return False

    # Calculate nucleotide counts (used for multiple rules)
    nucleotide_counts = Counter(primer)
    gc_count = nucleotide_counts.get('G', 0) + nucleotide_counts.get('C', 0)

    # Rule 2: GC content filtering
    gc_content = gc_count / len(primer)
    if not (parameter.gc_min <= gc_content <= parameter.gc_max):
        logger.debug(f"GC content filter: {primer} GC={gc_content:.1%}")
        return False

    # Rule 3: GC clamp (last 5 bases)
    # Adaptive based on genome GC content to avoid eliminating valid primers
    # for AT-rich (e.g. Plasmodium ~25% GC) or GC-rich (e.g. Mycobacterium ~65% GC) targets
    gc_in_last_5 = _count_gc(primer[-5:])
    genome_gc = getattr(parameter, 'genome_gc', None)
    if genome_gc is not None and genome_gc < 0.30:
        # AT-rich genome: allow 0 GC in last 5, reject >3
        gc_clamp_min = 0
        gc_clamp_max = MAX_GC_IN_LAST_5_BASES
    elif genome_gc is not None and genome_gc > 0.70:
        # GC-rich genome: allow up to 4 GC in last 5, require >=1 AT
        gc_clamp_min = 1
        gc_clamp_max = 4
    else:
        # Standard: 1-3 GC in last 5
        gc_clamp_min = 1
        gc_clamp_max = MAX_GC_IN_LAST_5_BASES
    if gc_in_last_5 > gc_clamp_max or gc_in_last_5 < gc_clamp_min:
        logger.debug(f"GC clamp filter: {primer} GC_last5={gc_in_last_5} "
                     f"(allowed {gc_clamp_min}-{gc_clamp_max})")
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
    if dimer.is_dimer(primer, primer, parameter.max_self_dimer_bp):
        logger.debug(f"Self-dimer filter: {primer}")
        return False

    return True


def get_all_rates(
    primer_list: List[str],
    fg_prefixes: List[str],
    bg_prefixes: List[str],
    fg_total_length: int,
    bg_total_length: int
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
    bloom_path = getattr(parameter, 'bloom_filter_path', None) or getattr(parameter, 'bg_bloom', None)
    use_bloom = getattr(parameter, 'use_bloom_filter', bloom_path is not None)

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
    for primer in primer_list:
        primer_len = len(primer)
        scaled_min_fg = _scale_freq_threshold(parameter.min_fg_freq, primer_len)
        scaled_max_bg = _scale_freq_threshold(parameter.max_bg_freq, primer_len)
        fg_count = primer_to_fg_count.get(primer, None)
        fg_bool = (fg_count is None or fg_count / fg_total_length > scaled_min_fg)
        bg_count = primer_to_bg_count.get(primer, None)
        bg_bool = (bg_count is None or bg_count / bg_total_length < scaled_max_bg)
        results.append([primer, fg_count, bg_count, fg_bool, bg_bool])

    df = pd.DataFrame(results, columns=['primer', 'fg_count', 'bg_count', 'fg_bool', 'bg_bool'])

    return df


def get_rates_for_one_species(
    primer_list: List[str], fname_prefixes: List[str]
) -> Dict[str, int]:
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

    # Use context manager to ensure proper pool cleanup (prevents resource leaks)
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        results = pool.map(_get_rate_for_one_file, tasks)

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
    primer_to_count = {}
    with open(fname_prefix + '_' + str(k) + 'mer_all.txt', 'r') as f_in:
        for line in f_in:
            primer = line.split(" ")[0]
            count = line.split(" ")[1]
            primer_to_count[primer] = int(count)

    all_counts = []
    for primer in primer_list:
        if primer in primer_to_count:
            all_counts.append(primer_to_count[primer])
        else:
            all_counts.append(0)
    return dict(zip(primer_list, all_counts))


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
    df['gini'] = primer_attributes.get_gini_from_txt(df['primer'].values, fg_prefixes, fg_genomes, fg_seq_lengths, circular, position_cache=position_cache)

    if len(df['gini']) == 0:
        df['gini_bool'] = []
        return df

    # Vectorized boolean operation (10-50x faster than apply with lambda)
    df['gini_bool'] = df['gini'].notna() & (df['gini'] < parameter.max_gini)

    return df[df['gini_bool']]
