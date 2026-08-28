#!/usr/bin/env python3
"""Fetch the reference genomes the specificity validation needs.

The published-set benchmarks in `tests/validation/data/` record measured
wet-lab outcomes -- percent of reads on target, fold enrichment -- but not the
genomes those primers were designed against. Without the genomes the tests can
only check numbers copied from the papers' tables against each other. With
them, predicted selectivity can be scored against a measured outcome, which is
the only way to find out whether the model is right.

The genomes are not committed: the smallest useful pair is ~50 Mb and human is
3 Gb. This fetches them on request, into a directory the validation tests look
for and skip without.

    python scripts/fetch_reference_genomes.py --list
    python scripts/fetch_reference_genomes.py prevotella human_chr21
    python scripts/fetch_reference_genomes.py --all

A note on the human background. Four of the seven datasets used whole human
genomic DNA. `human_chr21` is 46 Mb against 3.1 Gb, so it is a subset and not
the background those experiments actually had -- repeat content and GC differ
between chromosomes, and a primer's background load scales with the sequence it
is counted against. It is offered because it makes the ranking question
answerable at a cost that fits in a test run: whether the six Prevotella sets
rank in their measured order is a question about relative load, which a
representative subset can address. Absolute selectivity figures from it are not
comparable to the papers'. `human_full` is there for when the absolute number
matters.
"""

import argparse
import gzip
import hashlib
import os
import shutil
import sys
import urllib.request
from dataclasses import dataclass
from typing import Dict, Optional

DEFAULT_DIR = os.path.join("tests", "validation", "genomes")


@dataclass(frozen=True)
class Reference:
    key: str
    description: str
    url: str
    filename: str
    approx_mb: int
    # SHA-256 of the decompressed FASTA. None until first fetch records one --
    # see --print-checksums. A wrong genome silently producing plausible
    # numbers is the failure this guards against.
    sha256: Optional[str] = None


REFERENCES: Dict[str, Reference] = {
    "prevotella": Reference(
        key="prevotella",
        description="Prevotella melaninogenica ATCC 25845 (target, Dwivedi-Yu 2023)",
        url=(
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/144/405/"
            "GCF_000144405.1_ASM14440v1/GCF_000144405.1_ASM14440v1_genomic.fna.gz"
        ),
        filename="prevotella.fna",
        approx_mb=3,
        sha256="87815b0e869cd5dad3b5e4966945304ece4a6874eb1994507e83248f20e4b0ac",
    ),
    "human_chr21": Reference(
        key="human_chr21",
        description="Human chr21, hg38 (background subset; 46 Mb of 3.1 Gb)",
        url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz",
        filename="human_chr21.fna",
        approx_mb=46,
        sha256="35c71b68436d1a278ecb6a1e875af3ba4020738a028a7feac769a6d62790ae1f",
    ),
    "mtb": Reference(
        key="mtb",
        description="M. tuberculosis H37Rv (target, Clarke 2017)",
        url=(
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/"
            "GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz"
        ),
        filename="mtb.fna",
        approx_mb=4,
    ),
    "wolbachia": Reference(
        key="wolbachia",
        description="Wolbachia pipientis wMel (target, Clarke 2017 / Leichty 2014)",
        url=(
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/025/"
            "GCF_000008025.1_ASM802v1/GCF_000008025.1_ASM802v1_genomic.fna.gz"
        ),
        filename="wolbachia.fna",
        approx_mb=1,
    ),
    "human_full": Reference(
        key="human_full",
        description="Human GRCh38 primary assembly (full background; slow)",
        url=(
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/"
            "GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
        ),
        filename="human_full.fna",
        approx_mb=3100,
    ),
}


def sha256_of(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(reference: Reference, directory: str, force: bool = False) -> str:
    target = os.path.join(directory, reference.filename)
    if os.path.exists(target) and not force:
        print(f"  {reference.key}: already present ({target})")
        return target

    os.makedirs(directory, exist_ok=True)
    archive = target + ".gz"
    print(f"  {reference.key}: downloading ~{reference.approx_mb} Mb ...")
    urllib.request.urlretrieve(reference.url, archive)

    print(f"  {reference.key}: decompressing ...")
    with gzip.open(archive, "rb") as src, open(target, "wb") as dst:
        shutil.copyfileobj(src, dst)
    os.remove(archive)

    checksum = sha256_of(target)
    if reference.sha256 and checksum != reference.sha256:
        raise SystemExit(
            f"{reference.key}: checksum mismatch.\n"
            f"  expected {reference.sha256}\n  got      {checksum}\n"
            f"The assembly at that URL is not the one this benchmark was built "
            f"against. Validating against a different genome would produce "
            f"plausible numbers that mean nothing."
        )
    print(
        f"  {reference.key}: done ({os.path.getsize(target) / 1e6:.0f} Mb, sha256 {checksum[:16]}...)"
    )
    return target


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("names", nargs="*", help="References to fetch")
    parser.add_argument("--all", action="store_true", help="Fetch everything except human_full")
    parser.add_argument("--list", action="store_true", help="List available references")
    parser.add_argument("-d", "--directory", default=DEFAULT_DIR)
    parser.add_argument("--force", action="store_true", help="Re-download if present")
    parser.add_argument(
        "--print-checksums",
        action="store_true",
        help="Print sha256 of what is already downloaded, to paste into REFERENCES",
    )
    args = parser.parse_args()

    if args.list:
        print(f"{'key':<14}{'size':>8}  description")
        for ref in REFERENCES.values():
            print(f"{ref.key:<14}{ref.approx_mb:>6} Mb  {ref.description}")
        return 0

    if args.print_checksums:
        for ref in REFERENCES.values():
            path = os.path.join(args.directory, ref.filename)
            if os.path.exists(path):
                print(f'    "{ref.key}": "{sha256_of(path)}",')
        return 0

    names = list(args.names)
    if args.all:
        names = [k for k in REFERENCES if k != "human_full"]
    if not names:
        parser.error("name a reference, or pass --all or --list")

    unknown = [n for n in names if n not in REFERENCES]
    if unknown:
        parser.error(f"unknown reference(s): {unknown}. Try --list.")

    print(f"Fetching into {args.directory}/")
    for name in names:
        fetch(REFERENCES[name], args.directory, force=args.force)

    print(
        "\nNext: build k-mer indexes for the sizes the benchmark uses, e.g.\n"
        "  neoswga count-kmers -j <params.json pointing at these genomes>"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
