"""The blacklist frequency gate, and the batching it depends on.

Extracted from `core/pipeline.py` on 2026-09-25 when that module reached its
size budget. `pipeline` re-exports the name, so every existing importer --
`unified_optimizer`, `cli/iterate` and five tests -- is unaffected and gets
the identical object.

The reason it is worth its own module is the rule it has to keep. A count
lookup has a fixed cost per CALL: against a KMC database `counts_for` builds a
database of the query set, intersects it and dumps the result, which is three
processes. Asking about one primer at a time paid that cost once per primer,
and on a four-prefix design it turned this gate into minutes. Group by k, ask
once per prefix, and the cost is paid once per group.
"""

import logging

from neoswga.core import kmer_tables

logger = logging.getLogger(__name__)


def _filter_blacklist_penalty(
    primers: list[str], bl_prefixes: list[str], bl_seq_lengths: list[int], max_bl_freq: float = 0.0
) -> tuple[list[bool], list[float]]:
    """Filter primers by blacklist genome frequency.

    Reads k-mer count files for blacklist genomes and calculates per-primer
    frequency. Primers exceeding max_bl_freq are rejected.

    Args:
        primers: List of primer sequences.
        bl_prefixes: Blacklist genome k-mer file prefixes.
        bl_seq_lengths: Blacklist genome lengths for frequency calculation.
        max_bl_freq: Maximum allowed blacklist frequency (0 = any hit rejects).

    Returns:
        Tuple of (boolean mask, list of bl_freq values).
    """
    if len(bl_seq_lengths) != len(bl_prefixes):
        raise ValueError(
            f"bl_seq_lengths has {len(bl_seq_lengths)} entries but bl_prefixes "
            f"has {len(bl_prefixes)}; each blacklist genome must have a known length. "
            f"Check params.json and the blacklist preparation step."
        )

    # One lookup per prefix and k, not one per primer. `counts_for` has a
    # fixed cost per CALL -- against a KMC database it builds a database of
    # the query set, intersects and dumps, which is three processes -- so
    # asking about one primer at a time paid that cost hundreds of times and
    # turned a filter step into minutes. The batched form is the shape
    # `counts_for` was written for, and it is why `get_rates_for_one_species`
    # groups by k before asking.
    by_k: dict[int, list[str]] = {}
    for primer in primers:
        by_k.setdefault(len(primer), []).append(primer)

    counts: dict[str, int] = dict.fromkeys(primers, 0)
    total_length = sum(bl_seq_lengths)
    for prefix in bl_prefixes:
        for k, group in by_k.items():
            if not kmer_tables.table_exists(prefix, k):
                continue
            try:
                found = kmer_tables.counts_for(prefix, k, group)
            except Exception as e:
                logger.debug(f"Ignored error reading kmer file for blacklist penalty: {e}")
                continue
            for primer, count in found.items():
                counts[primer] += count

    mask = []
    bl_freqs = []
    for primer in primers:
        freq = counts[primer] / total_length if total_length > 0 else 0.0
        bl_freqs.append(freq)
        mask.append(freq <= max_bl_freq)
    return mask, bl_freqs
