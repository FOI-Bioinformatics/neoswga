"""Per-reference evaluation: each genome on its own, and honest reductions.

Phase 2 of `docs/superpowers/plans/2026-10-02-genomic-diversity-and-multi-host.md`.

Three things are pinned here, each one a defect this repository has shipped in
another place.

**A pooled figure hides a target.** Everything the pipeline reports over several
references is summed: counts over summed length in the two frequency gates,
sites and selectivity over flat sums of prefixes. A primer present in one strain
and absent from another therefore passes exactly as one present in both. The
first test builds that case and checks that `per_target` separates them while
the pooled number does not.

**A host nobody could measure must not read as a clean host.** The reduction is
the test that matters: a worst case computed over the members that happened to
answer is not a worst case, and `host_profile.aggregate_loads` returning 0.0 for
an empty list is the same shape. So the assertion is on
`worst_host_selectivity_density` being None, not merely on the one host record
saying `unavailable`.

**The two site routes agree.** A host given by FASTA with no table and no
params.json entry must count the same sites as the same host configured with a
table, or the new route is a second implementation rather than a fallback.
"""

from __future__ import annotations

import random

import pytest

from neoswga.core.reference_panel_evaluation import (
    ROLE_HOST,
    ROLE_TARGET,
    ReferenceSpec,
    evaluate_reference_panel,
)

K = 12


# ----------------------------------------------------------------------
# Synthetic references, written as FASTA so the scan route is exercised
# ----------------------------------------------------------------------


def _sequence(length, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _write_fasta(path, sequence, name="ref"):
    path.write_text(
        f">{name}\n" + "\n".join(sequence[i : i + 70] for i in range(0, len(sequence), 70)) + "\n"
    )
    return str(path)


def _plant(sequence, plants):
    """Write each primer into the sequence at the given offsets."""
    seq = list(sequence)
    for primer, offsets in plants.items():
        for offset in offsets:
            seq[offset : offset + len(primer)] = list(primer)
    return "".join(seq)


def _canonical_counts(sequence, kmers, k):
    """Occurrences of each k-mer and its reverse complement, counted directly.

    An oracle, not a call into the package: the point of the comparison is that
    the table route and the scan route agree, and a shared helper would let them
    agree by construction.
    """
    counts = {}
    for kmer in kmers:
        rc = _reverse_complement(kmer)
        total = 0
        for index in range(len(sequence) - k + 1):
            window = sequence[index : index + k]
            if window in (kmer, rc):
                total += 1
        counts[kmer] = total
    return counts


@pytest.fixture
def shared_primer():
    """A 12-mer planted in both targets."""
    return "ACGTTAGGCATA"


@pytest.fixture
def only_in_a():
    """A 12-mer planted in target A alone."""
    return "TTACGCAGTTCA"


@pytest.fixture
def two_targets(tmp_path, shared_primer, only_in_a):
    """Two 20 kb targets. One primer binds both; one binds only the first."""
    a = _plant(
        _sequence(20_000, seed=11),
        {shared_primer: (1_000, 9_000), only_in_a: (3_000, 12_000, 17_000)},
    )
    b = _plant(_sequence(20_000, seed=22), {shared_primer: (2_000, 11_000)})
    specs = [
        ReferenceSpec(
            prefix=str(tmp_path / "target_a"),
            genome=_write_fasta(tmp_path / "target_a.fasta", a, "a"),
            length=len(a),
            role=ROLE_TARGET,
        ),
        ReferenceSpec(
            prefix=str(tmp_path / "target_b"),
            genome=_write_fasta(tmp_path / "target_b.fasta", b, "b"),
            length=len(b),
            role=ROLE_TARGET,
        ),
    ]
    return specs


# ----------------------------------------------------------------------
# A primer that binds one target of two
# ----------------------------------------------------------------------


def test_per_target_separates_a_primer_that_binds_only_one_target(
    two_targets, shared_primer, only_in_a
):
    panel = evaluate_reference_panel([shared_primer, only_in_a], two_targets, extension=500)

    a, b = panel.targets
    assert a.per_primer_sites[only_in_a] >= 3
    assert b.per_primer_sites[only_in_a] == 0
    assert b.per_primer_status[only_in_a] == "no_sites"
    # Both references answered, so the zero is a measurement rather than an
    # absence of one. That distinction is the whole reason for the status field.
    assert a.status == "measured" and b.status == "measured"


def test_the_pooled_site_count_alone_does_not_show_it(two_targets, shared_primer, only_in_a):
    """The pooled figure is consistent with several different per-target shapes.

    This is the defect, stated as a test: summing over references gives a number
    a panel binding both targets evenly would also produce, so the pooled figure
    cannot be read as evidence about either target.
    """
    panel = evaluate_reference_panel([shared_primer, only_in_a], two_targets, extension=500)

    a, b = panel.targets
    pooled = a.sites.value + b.sites.value
    assert pooled == pytest.approx(a.sites.value + b.sites.value)
    # The asymmetry is visible per target and invisible in the sum.
    assert a.sites.value > b.sites.value
    assert a.per_primer_sites[only_in_a] != b.per_primer_sites[only_in_a]


def test_worst_target_coverage_is_the_lower_of_the_two(two_targets, shared_primer, only_in_a):
    panel = evaluate_reference_panel([shared_primer, only_in_a], two_targets, extension=500)

    a, b = panel.targets
    assert panel.worst_target_coverage.value == pytest.approx(
        min(a.coverage.value, b.coverage.value)
    )
    assert b.prefix in panel.worst_target_coverage.basis


# ----------------------------------------------------------------------
# A host by FASTA with no table, against the same host with one
# ----------------------------------------------------------------------


def test_a_host_by_fasta_counts_the_same_sites_as_a_host_with_a_table(
    tmp_path, shared_primer, only_in_a
):
    """The fallback must agree with the route it replaces, or it is a second answer.

    Both sides go through `kmer_tables.counts_for`: one prefix carries a text
    table and no genome, the other carries a genome and no table, so the second
    is answered by `query_scan`. A third, independent count is computed in this
    file as the oracle.
    """
    sequence = _plant(
        _sequence(30_000, seed=33), {shared_primer: (4_000, 20_000), only_in_a: (9_000,)}
    )
    fasta = _write_fasta(tmp_path / "host.fasta", sequence, "host")
    primers = [shared_primer, only_in_a]
    expected = _canonical_counts(sequence, primers, K)

    # The configured host: a prefix with a table beside it and no FASTA.
    configured = tmp_path / "configured" / "host"
    configured.parent.mkdir()
    table = configured.parent / f"{configured.name}_{K}mer_all.txt"
    table.write_text("".join(f"{p}\t{expected[p]}\n" for p in primers))

    # The host given by path: a prefix in a directory with no table at all.
    by_path = tmp_path / "bypath" / "host"
    by_path.parent.mkdir()

    def sites(spec):
        panel = evaluate_reference_panel(primers, [spec])
        record = panel.hosts[0]
        assert record.status == "measured", record.unavailable
        assert record.source == "counts"
        return record.sites.value, record.per_primer_sites

    table_sites, table_per_primer = sites(
        ReferenceSpec(
            prefix=str(configured), genome=None, length=len(sequence), role=ROLE_HOST, scan=False
        )
    )
    path_sites, path_per_primer = sites(
        ReferenceSpec(
            prefix=str(by_path), genome=fasta, length=len(sequence), role=ROLE_HOST, scan=False
        )
    )

    assert path_sites == table_sites
    assert path_per_primer == table_per_primer
    assert path_sites == sum(expected.values())


def test_a_counted_host_reports_coverage_as_unavailable_not_zero(tmp_path, shared_primer):
    """Counting answers the site question and cannot answer a positional one."""
    sequence = _plant(_sequence(10_000, seed=44), {shared_primer: (1_000, 6_000)})
    fasta = _write_fasta(tmp_path / "counted.fasta", sequence, "counted")

    panel = evaluate_reference_panel(
        [shared_primer],
        [
            ReferenceSpec(
                prefix=str(tmp_path / "counted"),
                genome=fasta,
                length=len(sequence),
                role=ROLE_HOST,
                scan=False,
            )
        ],
    )
    host = panel.hosts[0]
    assert host.sites.value == 2
    assert host.coverage.value is None
    assert "counted rather than located" in host.coverage.unavailable
    for measurement in (host.mean_gap, host.max_gap, host.gap_gini):
        assert measurement.value is None and measurement.unavailable


def test_scanning_a_host_gives_the_same_sites_as_counting_it(tmp_path, shared_primer):
    """One quantity, two routes. A disagreement here is a wrong answer somewhere."""
    sequence = _plant(_sequence(10_000, seed=55), {shared_primer: (1_000, 6_000, 8_000)})
    fasta = _write_fasta(tmp_path / "host.fasta", sequence, "host")
    prefix = str(tmp_path / "host")

    counted = evaluate_reference_panel(
        [shared_primer],
        [
            ReferenceSpec(
                prefix=prefix, genome=fasta, length=len(sequence), role=ROLE_HOST, scan=False
            )
        ],
    ).hosts[0]
    scanned = evaluate_reference_panel(
        [shared_primer],
        [
            ReferenceSpec(
                prefix=prefix, genome=fasta, length=len(sequence), role=ROLE_HOST, scan=True
            )
        ],
    ).hosts[0]

    assert scanned.source == "positions"
    assert scanned.sites.value == counted.sites.value == 3
    assert scanned.coverage.value is not None


# ----------------------------------------------------------------------
# The reduction over a set with one unavailable member
# ----------------------------------------------------------------------


def _host_that_cannot_be_measured(tmp_path):
    """A host with no table, no index and no readable genome.

    Not an error at construction: the reference is named, and the answer about it
    is unknown. That is the case the reductions have to survive.
    """
    return ReferenceSpec(
        prefix=str(tmp_path / "missing" / "nowhere"),
        genome=str(tmp_path / "missing" / "nowhere.fasta"),
        length=1_000_000,
        role=ROLE_HOST,
        scan=True,
    )


@pytest.fixture
def one_good_host_and_one_unmeasurable(tmp_path, shared_primer):
    sequence = _plant(_sequence(10_000, seed=66), {shared_primer: (1_000,)})
    target = _plant(_sequence(10_000, seed=77), {shared_primer: (500, 5_000)})
    specs = [
        ReferenceSpec(
            prefix=str(tmp_path / "target"),
            genome=_write_fasta(tmp_path / "target.fasta", target, "t"),
            length=len(target),
            role=ROLE_TARGET,
        ),
        ReferenceSpec(
            prefix=str(tmp_path / "good_host"),
            genome=_write_fasta(tmp_path / "good_host.fasta", sequence, "h"),
            length=len(sequence),
            role=ROLE_HOST,
            scan=True,
        ),
        _host_that_cannot_be_measured(tmp_path),
    ]
    return specs


def test_an_unmeasurable_host_is_unavailable_with_a_reason(
    one_good_host_and_one_unmeasurable, shared_primer
):
    panel = evaluate_reference_panel([shared_primer], one_good_host_and_one_unmeasurable)

    bad = panel.record(one_good_host_and_one_unmeasurable[2].prefix)
    assert bad.status == "unavailable"
    assert bad.sites.value is None
    assert bad.unavailable
    assert bad.site_density.value is None and bad.coverage.value is None


def test_the_worst_host_reduction_is_none_not_the_best_of_the_rest(
    one_good_host_and_one_unmeasurable, shared_primer
):
    """The single most important assertion in this phase.

    The measurable host HAS a density, so a reduction that skipped the
    unmeasurable one would return a number, and that number would read as the
    panel's worst case against its hosts. It is not one.
    """
    panel = evaluate_reference_panel([shared_primer], one_good_host_and_one_unmeasurable)

    measured = [pair for pair in panel.pairs if pair.density.value is not None]
    unmeasured = [pair for pair in panel.pairs if pair.density.value is None]
    assert measured and unmeasured, "the premise needs one of each"

    worst = panel.worst_host_density
    assert worst.value is None
    assert worst.value != min(pair.density.value for pair in measured)
    assert "unknown rather than the lowest of the rest" in worst.unavailable
    assert unmeasured[0].host in worst.unavailable


def test_the_worst_target_reduction_is_none_when_a_target_is_unavailable(tmp_path, shared_primer):
    """Same rule on the foreground, where the missing member is a strain."""
    good = _plant(_sequence(10_000, seed=88), {shared_primer: (1_000, 5_000)})
    specs = [
        ReferenceSpec(
            prefix=str(tmp_path / "good_target"),
            genome=_write_fasta(tmp_path / "good_target.fasta", good, "g"),
            length=len(good),
            role=ROLE_TARGET,
        ),
        ReferenceSpec(
            prefix=str(tmp_path / "absent" / "target"),
            genome=str(tmp_path / "absent" / "target.fasta"),
            length=10_000,
            role=ROLE_TARGET,
        ),
    ]
    panel = evaluate_reference_panel([shared_primer], specs)

    measured = [r for r in panel.targets if r.coverage.value is not None]
    assert measured, "the premise needs a target that did measure"
    assert panel.worst_target_coverage.value is None
    assert "unknown rather than the lowest of the rest" in panel.worst_target_coverage.unavailable


def test_a_primer_outside_the_index_makes_the_reference_unavailable(tmp_path, shared_primer):
    """A half-answered site union is not a measurement, and it propagates.

    The index here holds one of the two primers and scanning is not allowed, so
    the other primer's sites on this reference are unknown. Coverage computed
    from the sites that did answer is a lower bound; reporting it as the
    reference's coverage is the shape that scored outside oligo sets at 0%.
    """
    from neoswga.core import position_index

    outside = "TTACGCAGTTCA"
    prefix = str(tmp_path / "target")
    position_index.write_entries(
        position_index.index_path(prefix, K),
        {shared_primer: [100, 2_000]},
        record_starts=[0],
    )

    panel = evaluate_reference_panel(
        [shared_primer, outside],
        [ReferenceSpec(prefix=prefix, genome=None, length=5_000, role=ROLE_TARGET, scan=False)],
    )

    target = panel.targets[0]
    assert target.status == "unavailable"
    assert target.sites.value is None
    assert "unknown rather than zero" in target.unavailable
    # The primer that WAS indexed still reports its own sites, so the record says
    # what is known as well as what is not.
    assert target.per_primer_sites[shared_primer] == 2
    assert target.per_primer_sites[outside] is None
    assert target.per_primer_status[outside] == "not_indexed"
    assert panel.worst_target_coverage.value is None


def test_counts_answer_a_reference_that_is_not_scanned(tmp_path, shared_primer):
    """scan=False with a genome is the counts route, not an unknown.

    Worth pinning: the cost rule for a host-sized genome is "count, do not hold
    it in memory", and that must still yield a site measurement. Only the
    positional figures are unavailable.
    """
    sequence = _plant(_sequence(5_000, seed=99), {shared_primer: (100, 2_500)})
    panel = evaluate_reference_panel(
        [shared_primer],
        [
            ReferenceSpec(
                prefix=str(tmp_path / "counted"),
                genome=_write_fasta(tmp_path / "counted.fasta", sequence, "c"),
                length=len(sequence),
                role=ROLE_TARGET,
                scan=False,
            )
        ],
    )
    target = panel.targets[0]
    assert target.status == "measured" and target.source == "counts"
    assert target.sites.value == 2
    assert target.coverage.value is None
    # And therefore the worst-target reduction is unknown, not the sites figure.
    assert panel.worst_target_coverage.value is None


# ----------------------------------------------------------------------
# The cross table
# ----------------------------------------------------------------------


def test_density_uses_each_references_own_length(tmp_path, shared_primer):
    """Known Issue 6: the ratio moves with host size and the density does not.

    The same host sequence is presented twice, once with its true length and once
    with a length ten times larger, standing in for the substitution the issue
    describes. The count ratio is identical in both cases -- it carries no length
    -- and the density moves by the length ratio, which is the information the
    pooled figure throws away.
    """
    host_seq = _plant(_sequence(10_000, seed=123), {shared_primer: (1_000, 4_000)})
    target_seq = _plant(_sequence(10_000, seed=321), {shared_primer: (500, 5_000, 9_000)})
    target = ReferenceSpec(
        prefix=str(tmp_path / "t"),
        genome=_write_fasta(tmp_path / "t.fasta", target_seq, "t"),
        length=len(target_seq),
        role=ROLE_TARGET,
    )
    host_path = _write_fasta(tmp_path / "h.fasta", host_seq, "h")

    small = evaluate_reference_panel(
        [shared_primer],
        [
            target,
            ReferenceSpec(
                prefix=str(tmp_path / "h"),
                genome=host_path,
                length=len(host_seq),
                role=ROLE_HOST,
                scan=True,
            ),
        ],
    ).pairs[0]
    large = evaluate_reference_panel(
        [shared_primer],
        [
            target,
            ReferenceSpec(
                prefix=str(tmp_path / "h"),
                genome=host_path,
                length=len(host_seq) * 10,
                role=ROLE_HOST,
                scan=True,
            ),
        ],
    ).pairs[0]

    assert small.ratio.value == pytest.approx(large.ratio.value)
    assert large.density.value == pytest.approx(small.density.value * 10)
    assert "not comparable across hosts" in small.ratio.basis


def test_a_host_with_no_sites_is_recorded_rather_than_divided_by(tmp_path, shared_primer):
    absent = "GGCGCGATATCG"
    target_seq = _plant(_sequence(8_000, seed=7), {shared_primer: (1_000,)})
    host_seq = _sequence(8_000, seed=8)
    panel = evaluate_reference_panel(
        [shared_primer],
        [
            ReferenceSpec(
                prefix=str(tmp_path / "t"),
                genome=_write_fasta(tmp_path / "t.fasta", target_seq, "t"),
                length=len(target_seq),
                role=ROLE_TARGET,
            ),
            ReferenceSpec(
                prefix=str(tmp_path / "h"),
                genome=_write_fasta(tmp_path / "h.fasta", host_seq, "h"),
                length=len(host_seq),
                role=ROLE_HOST,
                scan=True,
            ),
        ],
    )
    assert absent not in panel.primers
    pair = panel.pairs[0]
    assert pair.zero_host_sites is True
    # Unbounded, and said so, rather than reported as a measurement of 1e6.
    assert "unbounded" in pair.density.basis
    assert panel.hosts[0].sites.value == 0


def test_the_json_form_is_serialisable_and_keeps_every_reason(
    one_good_host_and_one_unmeasurable, shared_primer
):
    import json

    panel = evaluate_reference_panel([shared_primer], one_good_host_and_one_unmeasurable)
    text = json.dumps(panel.as_dict(), allow_nan=False)
    restored = json.loads(text)

    assert set(restored) >= {
        "per_target",
        "per_host",
        "target_host_pairs",
        "worst_target_coverage",
        "worst_host_selectivity_density",
    }
    assert restored["worst_host_selectivity_density"]["value"] is None
    assert restored["worst_host_selectivity_density"]["unavailable"]


def test_a_reference_with_no_length_is_unavailable_rather_than_empty(tmp_path, shared_primer):
    """An empty FASTA has no denominator, so it has no measurement.

    Reporting zero sites for it would read as a clean genome, and every density,
    coverage figure and gap is a ratio against the length it does not have.
    """
    empty = tmp_path / "empty.fasta"
    empty.write_text(">empty\n\n")

    panel = evaluate_reference_panel(
        [shared_primer],
        [
            ReferenceSpec(
                prefix=str(tmp_path / "empty"),
                genome=str(empty),
                length=0,
                role=ROLE_HOST,
                scan=True,
            )
        ],
    )
    host = panel.hosts[0]
    assert host.status == "unavailable"
    assert host.sites.value is None
    assert "no length" in host.unavailable


def test_a_duplicate_reference_prefix_is_refused(tmp_path, shared_primer):
    """Two references under one key would overwrite each other's record."""
    spec = ReferenceSpec(
        prefix=str(tmp_path / "same"), genome=None, length=1_000, role=ROLE_TARGET, scan=False
    )
    with pytest.raises(ValueError, match="appears twice"):
        evaluate_reference_panel([shared_primer], [spec, spec])
