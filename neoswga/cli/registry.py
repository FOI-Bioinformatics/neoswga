"""Registry / library CLI handlers: ``background-*`` and ``genome-*``.

Manage the background-genome bloom-filter registry and the local genome
library. Extracted from cli_unified.py; uses only stdlib + lazy core imports.
"""

import logging
import os
import sys

logger = logging.getLogger(__name__)


def run_background_list(args):
    """List available pre-computed background genomes."""
    from neoswga.core.background_registry import BackgroundRegistry

    quiet = getattr(args, "quiet", False)

    registry = BackgroundRegistry(auto_discover=False)

    # Run discovery if requested
    if args.discover:
        if not quiet:
            logger.info("Discovering backgrounds...")
        n_discovered = registry.discover(verbose=not quiet)
        if not quiet:
            logger.info(f"Discovered {n_discovered} new backgrounds")
            logger.info("")

    # Search or list all
    if args.search:
        entries = registry.search(args.search)
        if not quiet:
            logger.info(f"Search results for '{args.search}':")
    else:
        entries = registry.list_all()
        if not quiet:
            logger.info("Available backgrounds:")

    if not entries:
        logger.info("  No backgrounds found")
        logger.info("")
        logger.info("Add backgrounds with:")
        logger.info(
            "  neoswga background-add --name 'Human GRCh38' --species 'Homo sapiens' --bloom-path /path/to/bloom.pkl"
        )
        logger.info("")
        logger.info("Or run discovery to find existing files:")
        logger.info("  neoswga background-list --discover")
        return

    logger.info("")
    for entry in entries:
        logger.info(f"  {entry.name}")
        if not quiet:
            logger.info(f"    Species: {entry.species}")
            if entry.genome_size > 0:
                size_gb = entry.genome_size / 1e9
                logger.info(f"    Size: {size_gb:.2f} Gbp")
            if entry.has_bloom:
                logger.info(f"    Bloom filter: {entry.bloom_path}")
            if entry.has_kmers:
                logger.info(f"    K-mer files: {entry.kmer_prefix}_*mer_all.txt")
                logger.info(f"    K range: {entry.k_range[0]}-{entry.k_range[1]}")
            logger.info("")

    if not quiet:
        logger.info(f"Total: {len(entries)} backgrounds")


def run_background_add(args):
    """Add a pre-computed background to the registry."""
    from neoswga.core.background_registry import BackgroundRegistry

    registry = BackgroundRegistry(auto_discover=False)

    # Validate paths
    if args.bloom_path and not os.path.exists(args.bloom_path):
        logger.error(f"Bloom filter not found: {args.bloom_path}")
        sys.exit(1)

    if args.kmer_prefix:
        # Check if at least one k-mer file exists
        found = False
        for k in range(args.min_k, args.max_k + 1):
            kmer_file = f"{args.kmer_prefix}_{k}mer_all.txt"
            if os.path.exists(kmer_file):
                found = True
                break
        if not found:
            logger.warning(f"No k-mer files found at {args.kmer_prefix}_*mer_all.txt")

    # Add to registry
    entry = registry.add(
        name=args.name,
        species=args.species,
        genome_size=args.genome_size,
        bloom_path=args.bloom_path,
        kmer_prefix=args.kmer_prefix,
        k_range=(args.min_k, args.max_k),
        description=args.description or "",
        overwrite=args.overwrite,
    )

    logger.info(f"Added background: {entry.name}")
    logger.info("")
    logger.info("Background genome registered. It will be used in subsequent filter runs.")


# =========================================================================
# UTILITY: Genome library management
# =========================================================================


def run_genome_add(args):
    """Add a genome to the pre-calculated library."""
    from neoswga.core.genome_library import GenomeLibrary

    # Parse k-ranges
    k_ranges = []
    for part in args.k_ranges.split(","):
        parts = part.strip().split("-")
        if len(parts) == 2:
            k_ranges.append((int(parts[0]), int(parts[1])))
        else:
            logger.error(f"Invalid k-range format: {part}. Use min-max (e.g., 6-12)")
            sys.exit(1)

    fasta_path = args.fasta
    if not os.path.exists(fasta_path):
        logger.error(f"FASTA file not found: {fasta_path}")
        sys.exit(1)

    library = GenomeLibrary()
    bloom = "no" if args.no_bloom else "auto"
    entry = library.add(
        name=args.name,
        fasta_path=fasta_path,
        role=args.role,
        species=args.species,
        k_ranges=k_ranges,
        build_bloom=bloom,
    )
    print(f"Added '{entry.name}' to genome library")
    print(f"  Size: {entry.genome_size:,} bp")
    print(f"  GC: {entry.gc_content:.1%}")
    print(f"  K-mer ranges: {entry.computed_k_ranges}")
    if entry.bloom_path:
        print(f"  Bloom filter: {entry.bloom_path}")


def run_genome_list(args):
    """List genomes in the pre-calculated library."""
    from neoswga.core.genome_library import GenomeLibrary

    library = GenomeLibrary()
    entries = library.list()

    if not entries:
        print("No genomes in library.")
        print("Add one with: neoswga genome-add <name> <fasta>")
        return

    print(f"Genome library ({len(entries)} entries):")
    print("-" * 60)
    for entry in entries:
        size_str = (
            f"{entry.genome_size / 1e6:.1f} Mbp"
            if entry.genome_size >= 1e6
            else f"{entry.genome_size:,} bp"
        )
        print(f"  {entry.name:<25} {size_str:>12}  GC={entry.gc_content:.1%}  role={entry.role}")
        if getattr(args, "verbose", False):
            print(f"    FASTA: {entry.fasta_path}")
            print(f"    K-mer ranges: {entry.computed_k_ranges}")
            if entry.bloom_path:
                print(f"    Bloom: {entry.bloom_path}")
            print(f"    Created: {entry.created_date}")


def run_genome_remove(args):
    """Remove a genome from the pre-calculated library."""
    from neoswga.core.genome_library import GenomeLibrary

    library = GenomeLibrary()
    if library.remove(args.name):
        print(f"Removed '{args.name}' from genome library")
    else:
        logger.error(f"Genome '{args.name}' not found in library")
        sys.exit(1)
