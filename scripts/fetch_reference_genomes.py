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
    # Whether the URL serves a gzip. The NCBI and UCSC paths do; the AWS
    # iGenomes mirror serves plain FASTA.
    compressed: bool = True


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
    "wolbachia_wri": Reference(
        key="wolbachia_wri",
        description="Wolbachia sp. wRi, supergroup A (target panel)",
        url=(
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/022/285/"
            "GCF_000022285.1_ASM2228v1/GCF_000022285.1_ASM2228v1_genomic.fna.gz"
        ),
        filename="wolbachia_wri.fna",
        approx_mb=1,
        sha256="e5a9790e907e2229720f0273873bef291888c16d934f1aa57ea76a16deac8b64",
    ),
    "wolbachia_wpip": Reference(
        key="wolbachia_wpip",
        description=(
            "Wolbachia wPip Pel, Culex quinquefasciatus, supergroup B (target panel)"
        ),
        url=(
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/073/005/"
            "GCF_000073005.1_ASM7300v1/GCF_000073005.1_ASM7300v1_genomic.fna.gz"
        ),
        filename="wolbachia_wpip.fna",
        approx_mb=1,
        sha256="368040a1de241d38350c360584f01a47ebad5f366ca7e53b672a223c2a4ef1c3",
    ),
    "wolbachia_walbb": Reference(
        key="wolbachia_walbb",
        description="Wolbachia pipientis wAlbB, supergroup B (target panel)",
        url=(
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/171/285/"
            "GCF_004171285.1_ASM417128v1/GCF_004171285.1_ASM417128v1_genomic.fna.gz"
        ),
        filename="wolbachia_walbb.fna",
        approx_mb=2,
        sha256="7d4bb3efd77d9b7feddc4d8472a3f76b824980a9fa0659e3684f2493798f26a4",
    ),
    "wolbachia_wbm": Reference(
        key="wolbachia_wbm",
        description=(
            "Wolbachia wBm strain TRS of Brugia malayi, supergroup D (target panel)"
        ),
        url=(
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/385/"
            "GCF_000008385.1_ASM838v1/GCF_000008385.1_ASM838v1_genomic.fna.gz"
        ),
        filename="wolbachia_wbm.fna",
        approx_mb=1,
        sha256="ddf06f1ea90c9f66ba4d08c9ef45f5968f5d5f484ce37d87c1d1e3ad57ffb668",
    ),
    # The genus Lactobacillus was split in 2020, so the organism a paper calls
    # Lactobacillus plantarum is Lactiplantibacillus plantarum in this assembly.
    # The key keeps the older genus name; the description records the current one.
    "lactobacillus": Reference(
        key="lactobacillus",
        description=(
            "Lactiplantibacillus plantarum WCFS1, a Drosophila gut commensal "
            "(small host, host panel)"
        ),
        url=(
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/203/855/"
            "GCF_000203855.3_ASM20385v3/GCF_000203855.3_ASM20385v3_genomic.fna.gz"
        ),
        filename="lactobacillus.fna",
        approx_mb=3,
        sha256="ee5924b0aae7801395d0398979178527b9b5f7c78e72f7b633077daeb0d52e76",
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
    # The three below are served by the AWS iGenomes mirror as plain FASTA.
    # They are here because the canonical hosts are not always reachable: a
    # cloud session's egress policy answered 403 to ftp.ncbi.nlm.nih.gov,
    # hgdownload.soe.ucsc.edu and ftp.ebi.ac.uk, which is every source the
    # entries above use. Each carries its own checksum, measured on download.
    "ecoli_k12": Reference(
        key="ecoli_k12",
        description="E. coli K-12 MG1655, 4,639,675 bp (real bacterial target)",
        url=(
            "https://ngi-igenomes.s3.amazonaws.com/igenomes/"
            "Escherichia_coli_K_12_MG1655/NCBI/2001-10-15/Sequence/"
            "WholeGenomeFasta/genome.fa"
        ),
        filename="ecoli_k12.fna",
        approx_mb=5,
        sha256="f6240859803fd3ebe6332fa14cc4783e29cf010d7c3529195b4afa1688f5b910",
        compressed=False,
    ),
    "human_chr1": Reference(
        key="human_chr1",
        description=(
            "Human chr1, hg38, 248,956,422 bp (background large enough that a "
            "12-mer cannot avoid it; 7.4% N)"
        ),
        url=(
            "https://ngi-igenomes.s3.amazonaws.com/igenomes/"
            "Homo_sapiens/UCSC/hg38/Sequence/Chromosomes/chr1.fa"
        ),
        filename="human_chr1.fna",
        approx_mb=254,
        sha256="04ee6db2e94ccc4daddc168453189d8c17a01454ba67ac0443a73a0401408ee0",
        compressed=False,
    ),
    # A SEPARATE key from `human_chr21`, deliberately. Same 46,709,983 bases of
    # hg38 chr21 in the same 50-column soft-masked layout, but NOT the same
    # bytes as UCSC's current chr21.fa.gz, so it would fail `human_chr21`'s
    # recorded checksum. Anything written against the UCSC file (for instance
    # `scripts/benchmarking/mismatch_model_ranking.py`) keeps asking for
    # `human_chr21.fna` and is not silently handed this one.
    "mtb_igenomes": Reference(
        key="mtb_igenomes",
        description=(
            "M. tuberculosis H37Rv, 4,411,532 bp, GC 65.6% (target, Clarke 2017); "
            "the Ensembl H37Rv.EB1 assembly as served by iGenomes, a SEPARATE key "
            "from `mtb` because that entry fetches NCBI's GCF_000195955.2"
        ),
        url=(
            "https://ngi-igenomes.s3.amazonaws.com/igenomes/"
            "Mycobacterium_tuberculosis_H37RV/Ensembl/H37Rv.EB1/Sequence/"
            "WholeGenomeFasta/genome.fa"
        ),
        filename="mtb_igenomes.fna",
        approx_mb=5,
        sha256="4d4e3f479cc1796de041c9efae898ada795dc0e3b209738bf1c76f85c95f1d9f",
        compressed=False,
    ),
    "saureus": Reference(
        key="saureus",
        description=(
            "S. aureus NCTC 8325, 2,821,361 bp, GC 32.9% (small host, host panel)"
        ),
        url=(
            "https://ngi-igenomes.s3.amazonaws.com/igenomes/"
            "Staphylococcus_aureus_NCTC_8325/NCBI/2006-02-13/Sequence/"
            "WholeGenomeFasta/genome.fa"
        ),
        filename="saureus.fna",
        approx_mb=3,
        sha256="f7750aa5049aebdf4bacc677331007c0b3c49687a3c728f3d93ae8d1c5d7e5f2",
        compressed=False,
    ),
    "human_chr21_igenomes": Reference(
        key="human_chr21_igenomes",
        description="Human chr21, hg38, via AWS iGenomes (see the note above)",
        url=(
            "https://ngi-igenomes.s3.amazonaws.com/igenomes/"
            "Homo_sapiens/UCSC/hg38/Sequence/Chromosomes/chr21.fa"
        ),
        filename="human_chr21_igenomes.fna",
        approx_mb=48,
        sha256="a44cf4cd4f6cb2565f7fc4520f3f3d908ec3dfa05de1f50d4ed655b0f640a207",
        compressed=False,
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
    print(f"  {reference.key}: downloading ~{reference.approx_mb} Mb ...")
    if reference.compressed:
        archive = target + ".gz"
        urllib.request.urlretrieve(reference.url, archive)
        print(f"  {reference.key}: decompressing ...")
        with gzip.open(archive, "rb") as src, open(target, "wb") as dst:
            shutil.copyfileobj(src, dst)
        os.remove(archive)
    else:
        urllib.request.urlretrieve(reference.url, target)

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
