"""A host given by path must count what the design's own index counted.

Phase 2 of `docs/superpowers/plans/2026-10-02-genomic-diversity-and-multi-host.md`.

`--background FASTA` exists so a delivered set can be scored against a host that
was never in the design's params.json. It is only worth having if it answers the
same question: the host site count it reports is the quantity the optimizer sums
into `total_bg_sites`, and a route that disagreed would be a second
implementation rather than a way in.

The two routes are genuinely different code. The configured host is read from
`drosophila_12mer_positions.h5` through `PositionCache`; the host by path has no
table and no index at its prefix, so `kmer_tables.counts_for` falls back to
`query_scan` and reads the 144 Mb reference once. Nothing is shared but the
answer.

Measured on the prepared Wolbachia design on 2026-10-02, 10 candidates from its
`step3_df.csv`: 578 target sites and 50 Drosophila sites by both routes,
selectivity 11.56 and density 1,311 by both. The by-path run took 11.9 s.

Skipped where the prepared design is absent, which is every checkout that has
not run `examples/wolbachia_pool_design/prepare.py`: the references and indexes
are gitignored. The plan asked for this check against the design's delivered set
0; that set does not exist on disk -- the prepared directory holds no
`step4_improved_df.csv` -- and `optimize` refuses its indexes because they carry
no recorded reference digest, so the comparison is made on candidates from the
same design's `step3_df.csv` instead. The quantity compared is the same one.
"""

from __future__ import annotations

import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DESIGN = ROOT / "examples" / "wolbachia_pool_design"
WORK = DESIGN / "work"
TARGET = DESIGN / "input" / "wmel.fna"
HOST = DESIGN / "input" / "drosophila.fna"
HOST_INDEX = WORK / "drosophila_12mer_positions.h5"
CANDIDATES = WORK / "step3_df.csv"

#: Enough oligos for the two routes to disagree if either were wrong, and few
#: enough that the single pass over the host stays around ten seconds.
HOW_MANY = 10

pytestmark = pytest.mark.skipif(
    not (TARGET.exists() and HOST.exists() and HOST_INDEX.exists() and CANDIDATES.exists()),
    reason=(
        "the prepared Wolbachia design is absent (gitignored); run "
        "examples/wolbachia_pool_design/prepare.py to build it"
    ),
)


@pytest.fixture(scope="module")
def candidates():
    import csv

    with open(CANDIDATES) as handle:
        rows = list(csv.DictReader(handle))
    return [row["primer"] for row in rows[:HOW_MANY]]


@pytest.fixture(scope="module")
def host_length():
    import json

    manifest = json.loads((WORK / "run_manifest.json").read_text())
    del manifest  # read only to confirm the directory is a finished run
    from neoswga.core.utility import get_all_seq_lengths

    return get_all_seq_lengths(fname_genomes=[str(HOST)], cpus=1)[0]


@pytest.fixture(scope="module")
def target_length():
    from neoswga.core.utility import get_all_seq_lengths

    return get_all_seq_lengths(fname_genomes=[str(TARGET)], cpus=1)[0]


def _sites_from_the_index(primers):
    """What the design itself counts: positions in its own host index."""
    from neoswga.core.position_cache import PositionCache

    prefix = str(WORK / "drosophila")
    cache = PositionCache([prefix], list(primers), on_missing="warn")
    unresolved = {primer for _, primer in cache.missing_primers}
    assert not unresolved, f"these candidates are not in the design's host index: {unresolved}"
    return sum(len(cache.get_positions(prefix, primer)) for primer in primers)


def _sites_by_path(primers, host_length, tmp_path):
    """What `--background` counts: one pass over the FASTA, no table, no index."""
    from neoswga.core.reference_panel_evaluation import (
        ROLE_HOST,
        ReferenceSpec,
        evaluate_reference_panel,
    )

    spec = ReferenceSpec(
        prefix=str(tmp_path / "drosophila"),
        genome=str(HOST),
        length=host_length,
        role=ROLE_HOST,
        scan=False,
    )
    panel = evaluate_reference_panel(primers, [spec])
    record = panel.hosts[0]
    assert record.status == "measured", record.unavailable
    assert record.source == "counts", "the premise is that no table exists at this prefix"
    return int(record.sites.value)


def test_the_two_routes_count_the_same_host_sites(candidates, host_length, tmp_path):
    from_index = _sites_from_the_index(candidates)
    by_path = _sites_by_path(candidates, host_length, tmp_path)

    assert from_index > 0, "a host with no sites would make the comparison vacuous"
    assert by_path == from_index


def test_the_per_host_record_carries_the_hosts_own_length(
    candidates, host_length, target_length, tmp_path
):
    """The density is formed from each reference's own length, not a pooled one.

    Drosophila is 113 times wMel here, so a pooled denominator and a per-host one
    differ by two orders of magnitude -- which is the size of the error Known
    Issue 6 describes.
    """
    from neoswga.core.reference_panel_evaluation import (
        ROLE_HOST,
        ROLE_TARGET,
        ReferenceSpec,
        evaluate_reference_panel,
    )

    panel = evaluate_reference_panel(
        candidates,
        [
            ReferenceSpec(
                prefix=str(WORK / "wmel"),
                genome=str(TARGET),
                length=target_length,
                role=ROLE_TARGET,
                circular=True,
                scan=False,
            ),
            ReferenceSpec(
                prefix=str(WORK / "drosophila"),
                genome=str(HOST),
                length=host_length,
                role=ROLE_HOST,
                scan=False,
            ),
        ],
    )
    target, host = panel.targets[0], panel.hosts[0]
    pair = panel.pairs[0]

    expected_density = (target.sites.value / target_length) / (host.sites.value / host_length)
    assert pair.density.value == pytest.approx(expected_density)
    # The count ratio is the figure that moves with host size; both are reported
    # and only one is comparable.
    assert pair.ratio.value == pytest.approx(target.sites.value / host.sites.value)
    assert pair.density.value > pair.ratio.value
