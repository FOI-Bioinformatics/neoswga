"""Pairwise shared-k-mer divergence across the target genomes of a panel.

This is a K-MER MEASURE AND NOT ANI. It counts how many exact canonical k-mers
two genomes have in common, so a single substitution removes up to k k-mers from
the shared set and an indel or a rearrangement removes a different number again.
The figures here order a panel from close to distant; they are not a
phylogenetic distance and they are not an alignment identity. Supergroup
assignments quoted elsewhere come from the literature, not from this script.

What is reported per unordered pair:

  intersection, union  sizes of the two sets of canonical k-mers
  jaccard              |A and B| / |A or B|, which is SYMMETRIC
  fraction_of_a_in_b   |A and B| / |A|
  fraction_of_b_in_a   |A and B| / |B|
  expected_*           the same figures if the two sets were drawn
                       independently from the canonical k-mer space

THE CHANCE FLOOR IS NOT ZERO. At k=12 the canonical space holds 8,390,656
k-mers and a 1.3 Mb genome occupies about a ninth of it, so two unrelated
genomes of this size share a measurable fraction by coincidence alone. The
expected figures are reported beside the measured ones so that a low Jaccard
can be read against the floor rather than against zero. A pair at the floor is
uninformative at this k, and the remedy is a longer k, not a different index.

The Jaccard index is the headline because the genomes differ in size (wBm at
1.08 Mb against wAlbB at 1.48 Mb), so a one-sided shared fraction would report a
different number depending on which genome was named first. Both one-sided
fractions are reported as well, because the gap between them is the size
asymmetry and is informative on its own.

Counts are read only through `kmer_tables.iter_table`, which raises when a
prefix was never counted. Nothing here substitutes a zero for a missing table:
an absent table would make two genomes look as though they share nothing.

MEMORY. Each genome's k-mers are held in a Python set. A Wolbachia genome has
at most about 1.5 million distinct canonical 12-mers, so the whole five-genome
panel is a few hundred megabytes and the exact sets are affordable. THIS DOES
NOT TRANSFER TO A HOST-SIZED GENOME: hg38 at k=12 saturates the 12-mer space
and a 144 Mb Drosophila table does not fit this shape comfortably either. A
host-scale version needs a bitset over the encoded k-mer space, or a sketch,
and a sketch answers a different question.

Usage:

    python scripts/benchmarking/wolbachia_panel_divergence.py \
        --manifest tests/validation/genomes/diversity_panel.json \
        --k 12 --output divergence.json
"""

import argparse
import itertools
import json
import os
import sys
import time

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from neoswga.core import kmer_tables  # noqa: E402

DEFAULT_MANIFEST = os.path.join("tests", "validation", "genomes", "diversity_panel.json")

MEASURE = "shared_exact_canonical_kmer_jaccard"
NOT_ANI = (
    "A k-mer measure, not ANI and not a phylogenetic distance. One substitution "
    "removes up to k shared k-mers; an indel or a rearrangement removes a "
    "different number. Use it to order a panel, not to quote an identity."
)


def canonical_space(k):
    """How many distinct canonical k-mers exist at this k.

    A k-mer and its reverse complement collapse to one entry, except for the
    k-mers that are their own reverse complement. Those exist only for even k,
    where there are 4**(k/2) of them.
    """
    total = 4**k
    palindromes = 4 ** (k // 2) if k % 2 == 0 else 0
    return (total + palindromes) // 2


def load_targets(manifest_path):
    """The target entries of the panel manifest, ordered by key.

    Ordering is by key rather than by file order so that the output does not
    change when the manifest is regenerated.
    """
    with open(manifest_path) as handle:
        manifest = json.load(handle)
    targets = [entry for entry in manifest["panel"] if entry.get("role") == "target"]
    if not targets:
        raise ValueError(f"{manifest_path} lists no entry with role 'target'.")
    return sorted(targets, key=lambda entry: entry["key"])


def prefix_for(entry, repo_root):
    """The k-mer table prefix for a manifest entry.

    A prefix is the genome path with its extension removed, which is the
    convention `count-kmers` writes under when `fg_prefixes` is derived from
    `fg_genomes`. The table itself is never named here; `kmer_tables` decides
    whether it is a database or a text dump.
    """
    path = entry["path"]
    if not os.path.isabs(path):
        path = os.path.join(repo_root, path)
    return os.path.splitext(path)[0]


def kmer_set(prefix, k):
    """Every distinct canonical k-mer in this table.

    Counts are discarded: this measure is presence and absence. A prefix with no
    table raises out of `iter_table` rather than reporting an empty set.
    """
    return {kmer for kmer, _count in kmer_tables.iter_table(prefix, k)}


def measure_panel(targets, k, repo_root):
    space = canonical_space(k)
    genomes = []
    sets = {}
    for entry in targets:
        prefix = prefix_for(entry, repo_root)
        started = time.time()
        kmers = kmer_set(prefix, k)
        sets[entry["key"]] = kmers
        genomes.append(
            {
                "key": entry["key"],
                "path": entry["path"],
                "prefix": os.path.relpath(prefix, repo_root),
                "length": entry.get("length"),
                "records": entry.get("records"),
                "reference": bool(entry.get("reference", False)),
                "distinct_kmers": len(kmers),
                "space_occupancy": len(kmers) / space,
                "read_seconds": round(time.time() - started, 2),
            }
        )

    pairs = []
    for first, second in itertools.combinations(sorted(sets), 2):
        a = sets[first]
        b = sets[second]
        shared = len(a & b)
        union = len(a | b)
        # The floor a pair of unrelated genomes of these sizes would reach.
        expected_shared = len(a) * len(b) / space
        expected_union = len(a) + len(b) - expected_shared
        pairs.append(
            {
                "a": first,
                "b": second,
                "distinct_a": len(a),
                "distinct_b": len(b),
                "intersection": shared,
                "union": union,
                "jaccard": shared / union if union else None,
                "fraction_of_a_in_b": shared / len(a) if a else None,
                "fraction_of_b_in_a": shared / len(b) if b else None,
                "expected_intersection_if_independent": expected_shared,
                "expected_jaccard_if_independent": (
                    expected_shared / expected_union if expected_union else None
                ),
            }
        )

    return {
        "measure": MEASURE,
        "k": k,
        "canonical_space": space,
        "not_ani": NOT_ANI,
        "genomes": genomes,
        "pairs": pairs,
    }


def print_report(result):
    k = result["k"]
    print(f"Shared exact canonical {k}-mers. {result['not_ani']}")
    print(f"Canonical {k}-mer space: {result['canonical_space']:,}")
    print()
    print(
        f"{'genome':<18}{'length':>12}{'distinct ' + str(k) + '-mers':>20}{'space occupancy':>18}"
    )
    for genome in result["genomes"]:
        length = genome["length"]
        length_text = f"{length:,}" if isinstance(length, int) else "unknown"
        print(
            f"{genome['key']:<18}{length_text:>12}{genome['distinct_kmers']:>20,}"
            f"{genome['space_occupancy']:>18.4f}"
        )
    print()
    header = (
        f"{'pair':<36}{'shared':>12}{'union':>12}{'jaccard':>10}"
        f"{'A in B':>10}{'B in A':>10}{'chance J':>10}"
    )
    print(header)
    for pair in result["pairs"]:
        label = f"{pair['a']} vs {pair['b']}"
        print(
            f"{label:<36}{pair['intersection']:>12,}{pair['union']:>12,}"
            f"{pair['jaccard']:>10.4f}{pair['fraction_of_a_in_b']:>10.4f}"
            f"{pair['fraction_of_b_in_a']:>10.4f}"
            f"{pair['expected_jaccard_if_independent']:>10.4f}"
        )
    print()
    ordered = sorted(result["pairs"], key=lambda pair: -pair["jaccard"])
    print("Most to least similar by Jaccard:")
    for rank, pair in enumerate(ordered, start=1):
        print(f"  {rank}. {pair['a']} vs {pair['b']}: {pair['jaccard']:.4f}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--manifest",
        default=os.path.join(REPO_ROOT, DEFAULT_MANIFEST),
        help="Panel manifest. Target entries are the ones measured.",
    )
    parser.add_argument("--k", type=int, default=12, help="k-mer length (default 12)")
    parser.add_argument("--output", help="Write the result as JSON to this path.")
    args = parser.parse_args(argv)

    targets = load_targets(args.manifest)
    result = measure_panel(targets, args.k, REPO_ROOT)
    result["manifest"] = os.path.relpath(os.path.abspath(args.manifest), REPO_ROOT)
    print_report(result)

    if args.output:
        with open(args.output, "w") as handle:
            json.dump(result, handle, indent=2, sort_keys=False)
            handle.write("\n")
        print(f"\nWrote {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
