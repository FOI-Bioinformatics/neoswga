"""The intact-site set from a variant table must equal the one from the sequence.

This is the oracle for the whole variant route. Everything else about it is
bookkeeping; the only thing that can quietly be wrong is the arithmetic that
turns a 1-based contig position into a concatenated offset and asks whether it
falls inside a binding site. An off-by-one there, or a record offset applied to
the wrong record, produces plausible figures rather than obviously broken ones.

So the test does not assert a tendency. SNPs are generated with a seed, written
BOTH as a variant table and as a mutated FASTA, and the two routes must agree
exactly:

    sites the variant route calls intact
      ==
    reference sites still found by scanning the mutated FASTA, at the same
    coordinates

SNPs only, so no coordinate shifts, which is what makes "at the same
coordinates" a fair comparison. A site the mutation CREATES in the mutated
FASTA is outside the comparison by construction, which is limit 1 of the route.
"""

import os
import random

import numpy as np
import pytest

from neoswga.core import string_search
from neoswga.core.position_cache import PositionCache
from neoswga.core.thermodynamics import reverse_complement
from neoswga.core.variant_sites import (
    SiteGeometry,
    evaluate_strain_panel,
    intact_mask,
    split_sites_by_strand,
)
from neoswga.core.variant_table import open_variants
from tests.variant_invariants import assert_site_counts_add_up

K = 12


def _panel(*args, **kwargs):
    """`evaluate_strain_panel`, with the site-count invariants asserted on the result."""
    block = evaluate_strain_panel(*args, **kwargs)
    assert_site_counts_add_up(block)
    return block


# ----------------------------------------------------------------------
# Building a reference, its variants, and the mutated sequence
# ----------------------------------------------------------------------


def _write_fasta(path, records):
    with open(path, "w") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 70):
                handle.write(sequence[start : start + 70] + "\n")


def _random_sequence(rng, length):
    return "".join(rng.choice("ACGT") for _ in range(length))


def _snp_positions(rng, records, count, forbidden=()):
    """`(record_name, 1-based position)` pairs, one per SNP, sorted per record.

    `forbidden` holds (record, position) pairs to leave alone, which is how a
    record's first and last bases are kept out of the way when a case needs
    them untouched.
    """
    chosen = []
    for name, sequence in records:
        positions = sorted(rng.sample(range(1, len(sequence) + 1), count))
        chosen.extend(
            (name, position) for position in positions if (name, position) not in set(forbidden)
        )
    return chosen


def _mutate(records, snps):
    """The mutated sequence per record, and the alternate base per SNP."""
    sequences = {name: list(sequence) for name, sequence in records}
    alternates = {}
    for name, position in snps:
        reference = sequences[name][position - 1]
        alternate = {"A": "C", "C": "G", "G": "T", "T": "A"}[reference]
        sequences[name][position - 1] = alternate
        alternates[(name, position)] = (reference, alternate)
    return {name: "".join(bases) for name, bases in sequences.items()}, alternates


def _write_vcf(path, records, snps, alternates, sample=None):
    lines = ["##fileformat=VCFv4.2"]
    for name, sequence in records:
        lines.append(f"##contig=<ID={name},length={len(sequence)}>")
    lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    if sample:
        header += "\tFORMAT\t" + sample
    lines.append(header)
    for name, position in snps:
        reference, alternate = alternates[(name, position)]
        row = f"{name}\t{position}\t.\t{reference}\t{alternate}\t.\tPASS\t."
        if sample:
            row += "\tGT\t1/1"
        lines.append(row)
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")


def _pick_primers(records, rng, per_record=20):
    """Primers that occur in the reference, half of them on the reverse strand.

    The reverse-strand half matters: a site is stored at the forward-strand
    offset of the k-mer for either strand, so the primer's 3' terminus sits at
    the right-hand end of the site on one strand and the left-hand end on the
    other. A panel of forward-only primers cannot tell the two apart.
    """
    primers = []
    for index, (_name, sequence) in enumerate(records):
        starts = rng.sample(range(50, len(sequence) - K - 50), per_record)
        for offset, start in enumerate(starts):
            word = sequence[start : start + K]
            primers.append(word if (offset + index) % 2 == 0 else reverse_complement(word))
    return sorted(set(primers))


def _cache(prefix, fasta, primers, circular=False):
    return PositionCache(
        [prefix], list(primers), genome_paths=[fasta], circular=circular, on_missing="scan"
    )


def _route_intact(cache, prefix, primers, strain, geometry):
    """Sites the variant route calls intact, per primer, as a set."""
    out = {}
    for primer in primers:
        per_strand = split_sites_by_strand(cache, prefix, primer)
        union = np.unique(np.concatenate([per_strand["forward"], per_strand["reverse"]]))
        mask = intact_mask(union, len(primer), strain.starts, strain.ends, geometry)
        out[primer] = {int(position) for position in union[mask]}
    return out


def _reference_sites(cache, prefix, primers):
    return {
        primer: {int(position) for position in cache.get_positions(prefix, primer)}
        for primer in primers
    }


def _sites_in(fasta, primers, circular=False):
    """Every forward-strand offset of each primer or its complement in a FASTA.

    The same scanner the position index is built with, so the record-join rule
    is the one the pipeline applies rather than a second reading of it.
    """
    wanted = sorted({*primers, *(reverse_complement(primer) for primer in primers)})
    found = string_search.get_all_positions_per_k(wanted, fasta, circular=circular)
    return {
        primer: {int(p) for p in (found.get(primer) or [])}
        | {int(p) for p in (found.get(reverse_complement(primer)) or [])}
        for primer in primers
    }


def _oracle(reference_fasta, mutated_fasta, primers, circular=False):
    """Reference sites that are still present, at the same offset, after mutation."""
    reference = _sites_in(reference_fasta, primers, circular)
    mutated = _sites_in(mutated_fasta, primers, circular)
    return {primer: reference[primer] & mutated[primer] for primer in primers}


def _case(tmp_path, seed, records, snp_count, sample=None, forbidden=()):
    """One reference, its SNPs, the mutated FASTA, and the opened table."""
    rng = random.Random(seed)
    reference_fasta = str(tmp_path / f"reference_{seed}.fna")
    mutated_fasta = str(tmp_path / f"mutated_{seed}.fna")
    vcf = str(tmp_path / f"variants_{seed}.vcf")

    _write_fasta(reference_fasta, records)
    snps = _snp_positions(rng, records, snp_count, forbidden=forbidden)
    mutated, alternates = _mutate(records, snps)
    _write_fasta(mutated_fasta, [(name, mutated[name]) for name, _ in records])
    _write_vcf(vcf, records, snps, alternates, sample=sample)

    primers = _pick_primers(records, random.Random(seed + 1000))
    prefix = str(tmp_path / f"reference_{seed}")
    table = open_variants(vcf, reference_fasta, prefix=prefix)
    cache = _cache(prefix, reference_fasta, primers)
    return {
        "prefix": prefix,
        "primers": primers,
        "table": table,
        "cache": cache,
        "reference_fasta": reference_fasta,
        "mutated_fasta": mutated_fasta,
        "snps": snps,
        "records": records,
    }


# ----------------------------------------------------------------------
# The oracle
# ----------------------------------------------------------------------


@pytest.mark.parametrize("seed", [1, 2, 3, 4, 5])
def test_the_intact_sites_are_exactly_the_sites_the_mutated_sequence_still_holds(tmp_path, seed):
    records = [("chr1", _random_sequence(random.Random(seed), 6000))]
    case = _case(tmp_path, seed, records, snp_count=250)

    strain = case["table"].strains[0]
    geometry = SiteGeometry(length=len(records[0][1]))
    route = _route_intact(case["cache"], case["prefix"], case["primers"], strain, geometry)
    oracle = _oracle(case["reference_fasta"], case["mutated_fasta"], case["primers"])

    assert route == oracle
    # The comparison is only worth anything if the mutations actually destroyed
    # some sites; otherwise every empty answer agrees with every other.
    reference = _reference_sites(case["cache"], case["prefix"], case["primers"])
    total = sum(len(reference[p]) for p in case["primers"])
    intact = sum(len(route[p]) for p in case["primers"])
    assert total - intact > 0, "no site was affected, so the equality above proves nothing"
    assert intact > 0, "every site was affected, so the equality above proves little"


@pytest.mark.parametrize("seed", [11, 12, 13])
def test_the_oracle_holds_on_two_records_with_a_snp_beside_the_join(tmp_path, seed):
    """A two-record reference, and a SNP on each side of the concatenation join.

    The join is where a record offset gets applied to the wrong record, and
    where the scanners manufacture k-mers that exist in neither record. A SNP
    at the last base of record 1 and the first base of record 2 puts both
    directly under test.
    """
    rng = random.Random(seed)
    records = [
        ("chrA", _random_sequence(rng, 3000)),
        ("chrB", _random_sequence(rng, 2500)),
    ]
    case = _case(tmp_path, seed, records, snp_count=120)

    # Rebuild with the two join-adjacent SNPs added, keeping everything else.
    snps = sorted(
        set(case["snps"]) | {("chrA", len(records[0][1])), ("chrB", 1)},
        key=lambda pair: ([name for name, _ in records].index(pair[0]), pair[1]),
    )
    mutated, alternates = _mutate(records, snps)
    mutated_fasta = str(tmp_path / f"mutated_join_{seed}.fna")
    vcf = str(tmp_path / f"variants_join_{seed}.vcf")
    _write_fasta(mutated_fasta, [(name, mutated[name]) for name, _ in records])
    _write_vcf(vcf, records, snps, alternates)

    table = open_variants(vcf, case["reference_fasta"], prefix=case["prefix"])
    geometry = SiteGeometry(length=sum(len(sequence) for _, sequence in records))
    route = _route_intact(
        case["cache"], case["prefix"], case["primers"], table.strains[0], geometry
    )
    oracle = _oracle(case["reference_fasta"], mutated_fasta, case["primers"])

    assert route == oracle


def test_a_strain_with_no_variants_reproduces_the_reference_figures_exactly(tmp_path):
    """Not approximately. The per-strain figures come out of the same two
    functions as the reference ones, so an empty variant set must land on the
    same numbers rather than near them."""
    rng = random.Random(77)
    records = [("chr1", _random_sequence(rng, 6000))]
    reference_fasta = str(tmp_path / "reference.fna")
    _write_fasta(reference_fasta, records)

    vcf = str(tmp_path / "none.vcf")
    _write_vcf(vcf, records, [], {}, sample="strain_with_nothing")

    primers = _pick_primers(records, random.Random(78))
    prefix = str(tmp_path / "reference")
    cache = _cache(prefix, reference_fasta, primers)
    table = open_variants(vcf, reference_fasta, prefix=prefix)

    block = _panel(cache, prefix, primers, len(records[0][1]), table, extension=500, circular=False)
    record = block["per_strain"]["strain_with_nothing"]
    assert_site_counts_add_up(block)

    assert record["status"] == "measured"
    assert record["affected_sites"]["value"] == 0.0
    assert record["intact_sites"]["value"] == float(block["reference_sites"])
    assert record["intact_site_fraction"]["value"] == 1.0
    for key in ("coverage_on_intact_sites", "mean_gap_bp", "max_gap_bp", "gap_gini"):
        reference_key = {"coverage_on_intact_sites": "coverage"}.get(key, key)
        assert record[key]["value"] == block["reference"][_figure_name(reference_key)]["value"]


def _figure_name(key):
    return {"mean_gap_bp": "mean_gap", "max_gap_bp": "max_gap", "gap_gini": "gap_gini"}.get(
        key, key
    )


# ----------------------------------------------------------------------
# The 3'-proximal split, on both strands
# ----------------------------------------------------------------------


def _one_site_case(tmp_path, primer_on_reverse, offset_from_three_prime):
    """A reference with exactly one site, and one SNP at a chosen 3' distance.

    `offset_from_three_prime` counts bases back from the primer's 3' terminal
    base, so 0 is the terminal base itself.
    """
    rng = random.Random(999)
    filler_left = _random_sequence(rng, 400)
    word = "ACGTTGCAGGTA"
    filler_right = _random_sequence(rng, 400)
    sequence = filler_left + word + filler_right
    start = len(filler_left)  # 0-based offset of the site

    primer = reverse_complement(word) if primer_on_reverse else word
    # The 3' terminus is at the right-hand end of the site for a forward-strand
    # primer and at the left-hand end for a reverse-strand one.
    anchor = start if primer_on_reverse else start + len(word) - 1
    variant = (
        anchor + offset_from_three_prime
        if primer_on_reverse
        else (anchor - offset_from_three_prime)
    )

    records = [("chr1", sequence)]
    reference_fasta = str(tmp_path / f"ref_{primer_on_reverse}_{offset_from_three_prime}.fna")
    _write_fasta(reference_fasta, records)
    vcf = str(tmp_path / f"var_{primer_on_reverse}_{offset_from_three_prime}.vcf")
    snps = [("chr1", variant + 1)]
    _mutated, alternates = _mutate(records, snps)
    _write_vcf(vcf, records, snps, alternates)

    prefix = str(tmp_path / f"ref_{primer_on_reverse}_{offset_from_three_prime}")
    cache = _cache(prefix, reference_fasta, [primer])
    table = open_variants(vcf, reference_fasta, prefix=prefix)
    block = _panel(
        cache, prefix, [primer], len(sequence), table, extension=200, circular=False, window=5
    )
    return block["per_strain"][table.strains[0].name]


@pytest.mark.parametrize("on_reverse", [False, True])
def test_a_variant_at_the_three_prime_terminus_is_proximal_on_either_strand(tmp_path, on_reverse):
    record = _one_site_case(tmp_path, on_reverse, offset_from_three_prime=0)
    assert record["affected_sites"]["value"] == 1.0
    assert record["affected_three_prime_proximal"]["value"] == 1.0
    assert record["affected_distal"]["value"] == 0.0


@pytest.mark.parametrize("on_reverse", [False, True])
def test_a_variant_at_the_far_end_of_the_site_is_distal_on_either_strand(tmp_path, on_reverse):
    record = _one_site_case(tmp_path, on_reverse, offset_from_three_prime=11)
    assert record["affected_sites"]["value"] == 1.0
    assert record["affected_three_prime_proximal"]["value"] == 0.0
    assert record["affected_distal"]["value"] == 1.0


def test_the_split_is_measured_from_the_other_end_on_the_two_strands(tmp_path):
    """The same distance from the LEFT edge of a site is proximal for a
    reverse-strand primer and distal for a forward-strand one. If the split
    ignored the strand, the two would agree, which is the defect."""
    forward = _one_site_case(tmp_path, False, offset_from_three_prime=11)
    reverse = _one_site_case(tmp_path, True, offset_from_three_prime=0)
    assert forward["affected_distal"]["value"] == 1.0
    assert reverse["affected_three_prime_proximal"]["value"] == 1.0


# ----------------------------------------------------------------------
# Coordinates and limits
# ----------------------------------------------------------------------


def test_the_intervals_are_int64(tmp_path):
    rng = random.Random(5)
    records = [("chr1", _random_sequence(rng, 2000))]
    reference_fasta = str(tmp_path / "reference.fna")
    _write_fasta(reference_fasta, records)
    snps = [("chr1", 100), ("chr1", 900)]
    _mutated, alternates = _mutate(records, snps)
    vcf = str(tmp_path / "v.vcf")
    _write_vcf(vcf, records, snps, alternates)

    table = open_variants(vcf, reference_fasta)
    strain = table.strains[0]
    assert strain.starts.dtype == np.int64
    assert strain.ends.dtype == np.int64
    assert list(strain.starts) == [99, 899]
    assert list(strain.ends) == [100, 900]


def test_a_second_record_is_offset_by_the_first_records_length(tmp_path):
    records = [("chrA", _random_sequence(random.Random(3), 1500)), ("chrB", "ACGT" * 300)]
    reference_fasta = str(tmp_path / "reference.fna")
    _write_fasta(reference_fasta, records)
    snps = [("chrB", 1)]
    _mutated, alternates = _mutate(records, snps)
    vcf = str(tmp_path / "v.vcf")
    _write_vcf(vcf, records, snps, alternates)

    strain = open_variants(vcf, reference_fasta).strains[0]
    assert list(strain.starts) == [1500]


def test_the_three_limits_are_written_into_the_block(tmp_path):
    rng = random.Random(21)
    records = [("chr1", _random_sequence(rng, 3000))]
    reference_fasta = str(tmp_path / "reference.fna")
    _write_fasta(reference_fasta, records)
    snps = [("chr1", 500)]
    _mutated, alternates = _mutate(records, snps)
    vcf = str(tmp_path / "v.vcf")
    _write_vcf(vcf, records, snps, alternates)

    primers = _pick_primers(records, random.Random(22), per_record=2)
    prefix = str(tmp_path / "reference")
    cache = _cache(prefix, reference_fasta, primers)
    block = _panel(
        cache,
        prefix,
        primers,
        len(records[0][1]),
        open_variants(vcf, reference_fasta, prefix=prefix),
        extension=500,
    )

    limits = " ".join(block["limits"]).lower()
    assert len(block["limits"]) == 3
    assert "gained" in limits
    assert "indel" in limits
    assert "absent from the reference" in limits
    assert block["three_prime_window_nt"] == 5
    assert "reporting split" in block["three_prime_window_basis"]


def test_a_primer_with_no_site_has_no_intact_fraction_rather_than_zero(tmp_path):
    """Zero of zero sites intact is not 0.0 intact. A 0.0 here would read as a
    primer whose sites the variants destroyed."""
    rng = random.Random(31)
    records = [("chr1", _random_sequence(rng, 2000))]
    reference_fasta = str(tmp_path / "reference.fna")
    _write_fasta(reference_fasta, records)
    snps = [("chr1", 700)]
    _mutated, alternates = _mutate(records, snps)
    vcf = str(tmp_path / "v.vcf")
    _write_vcf(vcf, records, snps, alternates)

    absent = "CGCGCGCGATCG"
    assert absent not in records[0][1]
    prefix = str(tmp_path / "reference")
    cache = _cache(prefix, reference_fasta, [absent])
    block = _panel(
        cache,
        prefix,
        [absent],
        len(records[0][1]),
        open_variants(vcf, reference_fasta, prefix=prefix),
        extension=500,
    )
    strain = block["per_strain"][next(iter(block["per_strain"]))]
    assert strain["per_primer_intact_fraction"][absent] is None
    assert block["per_primer_intact_fraction_every_strain"][absent] is None


def test_the_reference_fasta_must_exist(tmp_path):
    from neoswga.core.exceptions import ReferenceDataError

    vcf = tmp_path / "v.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    with pytest.raises(ReferenceDataError):
        open_variants(str(vcf), str(tmp_path / "nothing.fna"))
    assert not os.path.exists(tmp_path / "design_failure.json")


# ----------------------------------------------------------------------
# A palindromic oligo is one site, not two
# ----------------------------------------------------------------------

PALINDROME = "ACGCGCGCGT"  # its own reverse complement


def test_a_palindromic_site_is_counted_once_and_a_clean_strain_keeps_all_of_it(tmp_path):
    """A self-reverse-complementary oligo is stored under both strand keys at
    the same offset. Counting the denominator per strand and the numerator over
    the union reported a strain with no variants as having lost a third of the
    panel."""
    assert reverse_complement(PALINDROME) == PALINDROME
    rng = random.Random(404)
    left = _random_sequence(rng, 300)
    right = _random_sequence(rng, 300)
    ordinary = "ACGTTGCAAG"
    sequence = left + PALINDROME + right[:100] + ordinary + right[100:]
    assert sequence.count(PALINDROME) == 1 and sequence.count(ordinary) == 1
    assert reverse_complement(ordinary) not in sequence

    records = [("chr1", sequence)]
    reference_fasta = str(tmp_path / "reference.fna")
    _write_fasta(reference_fasta, records)
    vcf = str(tmp_path / "none.vcf")
    _write_vcf(vcf, records, [], {}, sample="clean")

    primers = [ordinary, PALINDROME]
    prefix = str(tmp_path / "reference")
    cache = _cache(prefix, reference_fasta, primers)
    assert len(cache.get_positions(prefix, PALINDROME, "forward")) == 1
    assert len(cache.get_positions(prefix, PALINDROME, "reverse")) == 1

    block = _panel(
        cache,
        prefix,
        primers,
        len(sequence),
        open_variants(vcf, reference_fasta, prefix=prefix),
        extension=100,
    )
    record = block["per_strain"]["clean"]
    assert block["reference_sites"] == 2
    assert record["intact_sites"]["value"] == 2.0
    assert record["intact_site_fraction"]["value"] == 1.0
    assert block["worst_strain_intact_fraction"]["value"] == 1.0


# ----------------------------------------------------------------------
# Intervals longer than one base
# ----------------------------------------------------------------------


def test_a_long_interval_that_starts_first_still_covers_a_later_site():
    """A 100 bp deletion starting at 100, then a SNP at 120. A site at
    [150, 162) lies inside the deletion and after the SNP: only the running
    maximum of the interval ends sees it, because the SNP's own end (121) does
    not reach it."""
    starts = np.array([100, 120], dtype=np.int64)
    ends = np.array([200, 121], dtype=np.int64)
    geometry = SiteGeometry(length=1000)
    sites = np.array([150, 300, 115], dtype=np.int64)
    mask = intact_mask(sites, 12, starts, ends, geometry)
    assert mask.tolist() == [False, True, False]


# ----------------------------------------------------------------------
# A circular reference, as `evaluate-set` treats one unless --linear
# ----------------------------------------------------------------------


def _wrap_primers(sequence, rng, count=6):
    """K-mers that cross the end of the ring, on both strands, plus ordinary ones."""
    length = len(sequence)
    primers = []
    for back in range(1, K):
        word = sequence[length - back :] + sequence[: K - back]
        primers.append(word if back % 2 else reverse_complement(word))
    for start in rng.sample(range(50, length - K - 50), count):
        primers.append(sequence[start : start + K])
    return sorted(set(primers))


def _circular_case(tmp_path, seed, records):
    rng = random.Random(seed)
    tag = f"circ_{seed}_{len(records)}"
    reference_fasta = str(tmp_path / f"{tag}_reference.fna")
    mutated_fasta = str(tmp_path / f"{tag}_mutated.fna")
    vcf = str(tmp_path / f"{tag}.vcf")
    _write_fasta(reference_fasta, records)

    first_name, _first = records[0]
    last_name, last = records[-1]
    # SNPs in the first k-1 and last k-1 bases of the ring, plus scattered ones.
    edge = {(first_name, p) for p in rng.sample(range(1, K), 3)}
    edge |= {(last_name, len(last) - p) for p in rng.sample(range(0, K - 1), 3)}
    snps = set(_snp_positions(rng, records, 60)) | edge
    order = [name for name, _ in records]
    snps = sorted(snps, key=lambda pair: (order.index(pair[0]), pair[1]))
    mutated, alternates = _mutate(records, snps)
    _write_fasta(mutated_fasta, [(name, mutated[name]) for name, _ in records])
    _write_vcf(vcf, records, snps, alternates)

    joined = "".join(sequence for _, sequence in records)
    primers = _wrap_primers(joined, random.Random(seed + 1))
    prefix = str(tmp_path / f"{tag}_reference")
    cache = _cache(prefix, reference_fasta, primers, circular=True)
    table = open_variants(vcf, reference_fasta, prefix=prefix)
    return {
        "cache": cache,
        "prefix": prefix,
        "primers": primers,
        "table": table,
        "reference_fasta": reference_fasta,
        "mutated_fasta": mutated_fasta,
        "length": len(joined),
    }


def _wrapping(case):
    reference = _reference_sites(case["cache"], case["prefix"], case["primers"])
    return {p: {s for s in reference[p] if s + K > case["length"]} for p in case["primers"]}


@pytest.mark.parametrize("seed", [21, 22, 23, 24])
def test_the_oracle_holds_on_a_circular_reference_across_the_wrap(tmp_path, seed):
    """The scanner wraps the concatenated sequence once, so a site may start in
    the last k-1 bases and finish in the first ones. Its bases there have to be
    the ones the mask looks at, or a SNP in the wrapped part is missed and the
    site is reported intact when it is not."""
    rng = random.Random(seed)
    case = _circular_case(tmp_path, seed, [("ring", _random_sequence(rng, 3000))])
    geometry = SiteGeometry(length=case["length"], circular=True)

    route = _route_intact(
        case["cache"], case["prefix"], case["primers"], case["table"].strains[0], geometry
    )
    oracle = _oracle(case["reference_fasta"], case["mutated_fasta"], case["primers"], True)
    assert route == oracle

    wrapping = _wrapping(case)
    lost_at_wrap = sum(len(wrapping[p] - route[p]) for p in case["primers"])
    assert lost_at_wrap > 0, "no wrap-around site was affected, so the wrap is untested"

    block = _panel(
        case["cache"],
        case["prefix"],
        case["primers"],
        case["length"],
        case["table"],
        extension=200,
        circular=True,
    )
    assert block["per_strain"][case["table"].strains[0].name]["not_assessed_sites"]["value"] == 0


@pytest.mark.parametrize("seed", [31, 32])
def test_the_oracle_holds_on_a_circular_two_record_reference(tmp_path, seed):
    """The ring is the concatenation, not each record: the scanner joins the
    last record's end to the first record's start, and so must the mask."""
    rng = random.Random(seed)
    records = [("partA", _random_sequence(rng, 1800)), ("partB", _random_sequence(rng, 1400))]
    case = _circular_case(tmp_path, seed, records)
    geometry = SiteGeometry(length=case["length"], circular=True)

    route = _route_intact(
        case["cache"], case["prefix"], case["primers"], case["table"].strains[0], geometry
    )
    oracle = _oracle(case["reference_fasta"], case["mutated_fasta"], case["primers"], True)
    assert route == oracle
    assert any(_wrapping(case).values()), "no site crosses the ring's join"


def test_a_wrap_site_in_a_linear_reading_is_not_assessed_and_never_intact(tmp_path):
    """A position the scanner only produces for a circular reference, read as
    linear, has no definition here. It is counted apart, and every figure that
    would quietly leave it out is unavailable."""
    rng = random.Random(41)
    case = _circular_case(tmp_path, 41, [("ring", _random_sequence(rng, 3000))])
    linear = SiteGeometry(length=case["length"], circular=False)
    empty = np.array([], dtype=np.int64)
    wrapped_seen = 0
    for primer in case["primers"]:
        sites = case["cache"].get_positions(case["prefix"], primer)
        mask = intact_mask(sites, K, empty, empty, linear)
        for site, intact in zip(sites.tolist(), mask.tolist(), strict=True):
            if site + K > case["length"]:
                wrapped_seen += 1
                assert not intact, "a site with no definition was reported intact"
    assert wrapped_seen > 0

    block = _panel(
        case["cache"],
        case["prefix"],
        case["primers"],
        case["length"],
        case["table"],
        extension=200,
        circular=False,
    )
    record = block["per_strain"][case["table"].strains[0].name]
    assert record["not_assessed_sites"]["value"] == float(wrapped_seen)
    assert record["intact_site_fraction"]["value"] is None
    assert "geometry" in record["intact_site_fraction"]["unavailable"]
    assert record["coverage_on_intact_sites"]["value"] is None
    assert block["worst_strain_intact_fraction"]["value"] is None


def test_the_three_prime_end_of_a_wrapped_site_is_measured_along_the_primer():
    """A forward-strand site starting 5 bases before the end of a 300 bp ring
    ends at base 6. A variant there is at the primer's 3' terminus; the same
    variant is 11 bases from the 3' end of a reverse-strand primer."""
    from neoswga.core.variant_sites import nearest_variant_distances, three_prime_offset

    geometry = SiteGeometry(length=300, circular=True)
    sites = np.array([295], dtype=np.int64)
    starts = np.array([6], dtype=np.int64)
    ends = np.array([7], dtype=np.int64)
    assert three_prime_offset(295, 12, "forward", geometry) == 6
    assert three_prime_offset(295, 12, "reverse", geometry) == 295
    assert intact_mask(sites, 12, starts, ends, geometry).tolist() == [False]
    assert nearest_variant_distances(sites, 12, "forward", starts, ends, geometry).tolist() == [0]
    assert nearest_variant_distances(sites, 12, "reverse", starts, ends, geometry).tolist() == [11]


def test_a_length_that_disagrees_with_the_fasta_is_refused(tmp_path):
    """The wrap point and every coverage denominator are that length."""
    from neoswga.core.exceptions import ReferenceDataError

    rng = random.Random(51)
    records = [("chr1", _random_sequence(rng, 1000))]
    reference_fasta = str(tmp_path / "reference.fna")
    _write_fasta(reference_fasta, records)
    vcf = str(tmp_path / "none.vcf")
    _write_vcf(vcf, records, [], {})
    prefix = str(tmp_path / "reference")
    primer = records[0][1][100:112]
    cache = _cache(prefix, reference_fasta, [primer])
    with pytest.raises(ReferenceDataError, match="length"):
        evaluate_strain_panel(
            cache, prefix, [primer], 999, open_variants(vcf, reference_fasta, prefix=prefix)
        )
