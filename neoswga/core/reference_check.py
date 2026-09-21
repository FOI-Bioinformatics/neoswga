"""Does this position index describe the genome the design names?

Separate from `position_cache.require_record_metadata`, which asks whether an
index is of a shape a current design can use. This asks a different question:
whether it was built from the right sequence.

The two are split because they are knowable at different layers. Shape is a
property of the file, so the cache can check it wherever it is used. Identity
is a relation between a prefix and a genome, and only the resolved request
names both. Reading the pairing out of the mutable `parameter` module instead
paired the prefixes a call was GIVEN with whatever genomes the module happened
to hold, which under `pytest -n 8` made a design refuse its own index because
another test had left a different FASTA there.

A stale index is the failure a passing test suite is least likely to catch.
Nothing about the file is wrong: every dataset is present, the geometry is
current, the positions are internally consistent. They describe another
sequence, so every coverage and specificity figure computed from them is about
a genome the user is not designing against.
"""

from __future__ import annotations

import logging
import os
from typing import Mapping, Sequence

from .exceptions import ReferenceDataError
from .string_search import RECORD_STARTS_KEY

logger = logging.getLogger(__name__)

__all__ = ["verify_reference_digests", "verify_index_geometry"]


def verify_reference_digests(manifest: Mapping[str, str], lengths: Sequence[int]) -> None:
    """Refuse any index whose recorded reference is not the one being designed against.

    `manifest` maps a k-mer prefix to the genome its index must have been built
    from. A prefix the manifest does not name is not checked, because an unknown
    expectation cannot be compared and inventing one would refuse a valid index.

    An index with no recorded digest is reported as unknown rather than as
    matching, and the remedy is the one recount that gives it one. That is the
    same distinction the k-mer table sidecar draws through `digest_algorithm`:
    a record that cannot be compared is not a record that agrees.
    """
    if not manifest:
        return

    from neoswga.core import kmer_counter

    problems = []
    for prefix, genome in sorted(manifest.items()):
        if not os.path.exists(str(genome)):
            # The genome itself is missing. That is a different failure with a
            # different remedy, and the pipeline reports it where it is read.
            continue
        paths = [f"{prefix}_{k}mer_positions.h5" for k in lengths]
        existing = [path for path in paths if os.path.exists(path)]
        if not existing:
            continue
        try:
            expected = kmer_counter.genome_fingerprint(str(genome))
        except OSError as exc:
            problems.append(f"{prefix}: reference {genome} could not be read ({exc})")
            continue
        for path in existing:
            recorded = _recorded_digest(path)
            if recorded is None:
                problems.append(f"{path}: no recorded reference digest to compare")
            elif recorded != expected:
                problems.append(f"{path}: built from a reference other than {genome}")

    if not problems:
        return
    raise ReferenceDataError(
        "position index",
        "; ".join(problems),
        "Re-run 'neoswga count-kmers -j params.json' followed by "
        "'neoswga filter -j params.json' so the index is rebuilt from the "
        "reference this design names.",
    )


def _recorded_digest(path: str):
    """The reference digest an index carries, or None when it carries none."""
    import h5py

    try:
        with h5py.File(path, "r") as handle:
            value = handle.attrs.get("reference_digest")
    except OSError as exc:
        raise ReferenceDataError(
            f"position index {path}",
            f"cannot be opened ({exc})",
            "Regenerate it with 'neoswga count-kmers' followed by 'neoswga filter'.",
        ) from exc
    return None if value is None else str(value)


def verify_index_geometry(manifest: Mapping[str, str], lengths: Sequence[int]) -> None:
    """Refuse an index with no record geometry for a multi-record reference.

    Record starts let a coverage window stop at a contig edge. An index
    without them lets the window run into the next record, crediting a panel
    with covering bases on a contig its site is not on, which on a fragmented
    assembly inflates coverage throughout.

    A SINGLE-record reference has no joins, so the geometry is not merely
    absent but unnecessary, and refusing it would force a recount for a defect
    that cannot apply. The shipped Wolbachia example is exactly that pair: its
    wMel index carries no record starts and wMel is one record, while its
    Drosophila index carries none and Drosophila has 1,870.

    This lives beside the digest check rather than in `PositionCache` for the
    same reason: which genome a prefix belongs to is a relation only the
    resolved request knows. Deciding it inside the evaluator means reading a
    mutable global and pairing it with the prefixes the call was given, which
    is how a design came to refuse its own index under `pytest -n 8`.
    """
    if not manifest:
        return

    import h5py

    problems = []
    for prefix, genome in sorted(manifest.items()):
        paths = [f"{prefix}_{k}mer_positions.h5" for k in lengths]
        existing = [path for path in paths if os.path.exists(path)]
        if not existing:
            continue
        records = _count_records(genome)
        if records is None or records <= 1:
            continue
        for path in existing:
            try:
                with h5py.File(path, "r") as handle:
                    has_geometry = RECORD_STARTS_KEY in handle
            except OSError as exc:
                raise ReferenceDataError(
                    f"position index {path}",
                    f"cannot be opened ({exc})",
                    "Regenerate it with 'neoswga count-kmers' then 'neoswga filter'.",
                ) from exc
            if not has_geometry:
                problems.append(f"{path}: no record geometry, and {genome} holds {records} records")

    if not problems:
        return
    raise ReferenceDataError(
        "position index",
        "; ".join(problems),
        "Regenerate with 'neoswga count-kmers -j params.json' followed by "
        "'neoswga filter -j params.json'. Without record starts a coverage "
        "window runs past a contig edge into the next record.",
    )


def _count_records(genome) -> "int | None":
    """Header lines in a FASTA, or None when it cannot be read.

    Counting ">" rather than loading the sequence: the loader would hold 8.5 GB
    for hg38 to answer a question this settles by reading bytes.
    """
    if not genome or not os.path.exists(str(genome)):
        return None
    try:
        count = 0
        with open(genome, "rb") as handle:
            for line in handle:
                if line.startswith(b">"):
                    count += 1
        return count or None
    except OSError:
        return None
