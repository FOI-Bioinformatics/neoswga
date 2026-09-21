"""A missing answer, a stale one and a recorded zero are three different things.

Task 3 of the 2026-09-21 valid-design plan. The silent-zero class this closes
has already been reached by four separate routes in this project: a scan that
found nothing past 2 Gb (Known Issue 5), an integer that saturated (6), a cache
asked for a prefix it did not hold (13), and a guard whose predicate asked the
wrong question (15). Each time the symptom was the same, because every route
ends at the same number: zero background sites, which reads as a perfectly
specific design.

Five situations are separated here, and exactly one of them is a measurement:

| Situation | Meaning | Design may proceed |
|---|---|---|
| recorded empty site array | this primer occurs nowhere on this reference | yes |
| no entry for a primer or prefix | nobody asked | no |
| index missing record geometry | concatenation-era; windows cross contigs | no |
| index format older than this version | unknown geometry | no |
| reference digest does not match | the index describes a different genome | no |

The last is the one a passing test suite is least likely to catch, because a
stale index is structurally perfect. It simply describes something else.
"""

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

from neoswga.core import string_search
from neoswga.core.exceptions import ReferenceDataError
from neoswga.core.position_cache import MissingPositionsError, PositionCache

PRIMER = "ACGTACGTACGT"
ABSENT = "TTTTTTTTTTTT"


def write_index(prefix, *, record_starts=(0, 60), genome=None, sites=(5, 70)):
    """An index written the modern way, or the old way when starts are None."""
    with h5py.File(f"{prefix}_12mer_positions.h5", "w"):
        pass
    string_search.write_to_h5py(
        {PRIMER: list(sites)}, prefix, record_starts=record_starts, genome_fname=genome
    )


@pytest.fixture
def genome(tmp_path):
    fasta = tmp_path / "two.fna"
    fasta.write_text(">one\n" + "A" * 60 + "\n>two\n" + "C" * 60 + "\n")
    string_search.clear_genome_cache()
    return fasta


@pytest.fixture
def indexed(tmp_path, genome):
    prefix = str(tmp_path / "fg")
    write_index(prefix, genome=str(genome))
    return prefix


# ---------------------------------------------------------------------------
# The one case that is a measurement
# ---------------------------------------------------------------------------


def test_a_recorded_empty_array_is_a_measurement_and_is_allowed(tmp_path, genome):
    """Zero sites, recorded, is the design working. It must not be refused.

    This is the case every other check in this file has to be careful not to
    catch. A candidate indexed against a host it binds nowhere is the most
    specific candidate available, and refusing it would reject the best
    answers the pool has.
    """
    prefix = str(tmp_path / "empty")
    with h5py.File(f"{prefix}_12mer_positions.h5", "w"):
        pass
    string_search.write_to_h5py(
        {PRIMER: [5], ABSENT: []}, prefix, record_starts=[0, 60], genome_fname=str(genome)
    )

    cache = PositionCache([prefix], [PRIMER, ABSENT])
    cache.require_entries([PRIMER, ABSENT])
    assert len(cache.get_positions(prefix, ABSENT, "both")) == 0
    assert cache.has_entry(prefix, ABSENT) is True


def test_an_empty_array_and_a_missing_entry_are_not_the_same_answer(indexed):
    cache = PositionCache([indexed], [PRIMER])

    assert cache.has_entry(indexed, PRIMER) is True
    assert cache.has_entry(indexed, "GGGGGGGGGGGG") is False


# ---------------------------------------------------------------------------
# The four that are not
# ---------------------------------------------------------------------------


def test_a_primer_with_no_entry_anywhere_is_refused(indexed):
    with pytest.raises(MissingPositionsError, match="GGGGGGGGGGGG"):
        PositionCache([indexed], [PRIMER]).require_entries([PRIMER, "GGGGGGGGGGGG"])


def test_a_prefix_the_cache_never_indexed_is_refused(indexed):
    """Known Issue 13: an unindexed prefix answered every lookup with zero."""
    cache = PositionCache([indexed], [PRIMER])

    with pytest.raises(MissingPositionsError, match="host"):
        cache.get_positions("host", PRIMER, "both")


def test_an_index_without_record_geometry_is_refused(tmp_path):
    prefix = str(tmp_path / "old")
    write_index(prefix, record_starts=None)

    with pytest.raises(ReferenceDataError, match="record geometry"):
        PositionCache([prefix], [PRIMER]).require_record_metadata([prefix])


def test_an_index_of_an_older_format_is_refused(tmp_path, genome):
    """A future format bump must refuse the previous one, not read it hopefully.

    The attribute is checked rather than inferred from which datasets happen to
    be present, because the next geometry change may add nothing a reader can
    notice by looking.
    """
    prefix = str(tmp_path / "v0")
    write_index(prefix, genome=str(genome))
    with h5py.File(f"{prefix}_12mer_positions.h5", "a") as handle:
        handle.attrs["index_format_version"] = string_search.INDEX_FORMAT_VERSION - 1

    with pytest.raises(ReferenceDataError, match="format"):
        PositionCache([prefix], [PRIMER]).require_record_metadata([prefix])


def test_an_index_built_from_a_different_reference_is_refused(tmp_path, genome):
    """The structurally perfect failure: a stale index describing another genome.

    Nothing about the file is wrong. Every dataset is present, the geometry is
    current, the positions are internally consistent. They simply belong to a
    different sequence, so every coverage figure computed from them describes
    something the user is not designing against.
    """
    prefix = str(tmp_path / "stale")
    write_index(prefix, genome=str(genome))

    other = tmp_path / "other.fna"
    other.write_text(">one\n" + "G" * 60 + "\n>two\n" + "T" * 60 + "\n")
    string_search.clear_genome_cache()

    with pytest.raises(ReferenceDataError, match="different reference"):
        PositionCache([prefix], [PRIMER]).require_record_metadata(
            [prefix], genomes={prefix: str(other)}
        )


def test_the_matching_reference_passes_the_same_check(tmp_path, genome):
    prefix = str(tmp_path / "current")
    write_index(prefix, genome=str(genome))

    PositionCache([prefix], [PRIMER]).require_record_metadata(
        [prefix], genomes={prefix: str(genome)}
    )


def test_an_index_with_no_recorded_digest_is_refused_when_one_is_expected(tmp_path, genome):
    """Unknown is not "matches". An index predating the digest cannot be vouched for."""
    prefix = str(tmp_path / "nodigest"),
    prefix = prefix[0]
    write_index(prefix, genome=None)
    with h5py.File(f"{prefix}_12mer_positions.h5", "a") as handle:
        handle.attrs["index_format_version"] = string_search.INDEX_FORMAT_VERSION

    with pytest.raises(ReferenceDataError, match="no recorded reference"):
        PositionCache([prefix], [PRIMER]).require_record_metadata(
            [prefix], genomes={prefix: str(genome)}
        )


def test_a_corrupt_index_is_refused_rather_than_read_as_empty(tmp_path, genome):
    prefix = str(tmp_path / "corrupt")
    write_index(prefix, genome=str(genome))
    path = f"{prefix}_12mer_positions.h5"
    with open(path, "r+b") as handle:
        handle.seek(0)
        handle.write(b"not-an-hdf5-file")

    with pytest.raises(ReferenceDataError):
        PositionCache([prefix], [PRIMER]).require_record_metadata([prefix])


# ---------------------------------------------------------------------------
# Every design path enforces it
# ---------------------------------------------------------------------------


def _calls_require_record_metadata(module_name, function_name):
    """Whether a function reaches `require_record_metadata`, transitively.

    Asserted on the PATH rather than on the presence of a call somewhere in the
    module. `tests/test_the_objective_reaches_the_stage_that_refines.py` records
    what the other shape costs: two tests each confirmed one end of a
    connection that did not exist.
    """
    import ast
    import importlib
    import inspect

    seen = set()
    queue = [(module_name, function_name)]
    while queue:
        module_name, function_name = queue.pop()
        if (module_name, function_name) in seen:
            # A name already walked, or a recursive call. Skip it and keep
            # going: returning here would abandon the rest of the queue and
            # report "does not reach" for a path that does.
            continue
        seen.add((module_name, function_name))
        try:
            module = importlib.import_module(module_name)
            tree = ast.parse(inspect.getsource(module))
        except (ImportError, OSError, SyntaxError):
            continue
        for node in ast.walk(tree):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            if node.name != function_name:
                continue
            for inner in ast.walk(node):
                if not isinstance(inner, ast.Call):
                    continue
                name = getattr(inner.func, "attr", None) or getattr(inner.func, "id", None)
                if name == "require_record_metadata":
                    return True
                if name:
                    queue.append((module_name, name))
    return False


@pytest.mark.parametrize(
    "module_name,function_name",
    [
        ("neoswga.core.unified_optimizer", "run_optimization"),
        ("neoswga.core.primer_expansion", "expand_primers"),
        ("neoswga.cli.plan_pool", "run_plan_pool"),
    ],
)
def test_every_design_entry_point_checks_the_index(module_name, function_name):
    """Only `plan-pool` used to. `optimize` and `expand-primers` did not.

    A check one command performs is not a contract. The two commands that
    skipped it are the two that most users run.
    """
    assert _calls_require_record_metadata(module_name, function_name), (
        f"{module_name}.{function_name} designs a panel without checking that the "
        "position index it scores against is current and describes this reference"
    )


# ---------------------------------------------------------------------------
# Reference identity, checked where both halves are named
# ---------------------------------------------------------------------------


def test_the_manifest_pairs_a_prefix_only_with_its_own_genome():
    """Positional pairing, and only where both halves are present.

    A prefix with no genome is absent rather than paired with a guess. Pairing
    it with whatever a mutable global held is what made a design refuse its own
    index under `pytest -n 8`: the prefix a call was given was checked against
    another test's FASTA.
    """
    from neoswga.core.design_request import resolve_design_request

    request = resolve_design_request(
        {
            "fg_prefixes": ["target", "second"],
            "fg_genomes": ["target.fna"],
            "bg_prefixes": ["host"],
            "bg_genomes": ["host.fna"],
            "polymerase": "phi29",
        }
    )

    assert request.reference_manifest() == {"target": "target.fna", "host": "host.fna"}


def test_an_index_from_another_genome_is_refused_at_the_request_layer(tmp_path, genome):
    from neoswga.core.reference_check import verify_reference_digests

    prefix = str(tmp_path / "stale")
    write_index(prefix, genome=str(genome))
    other = tmp_path / "other.fna"
    other.write_text(">one\n" + "G" * 60 + "\n>two\n" + "T" * 60 + "\n")
    string_search.clear_genome_cache()

    with pytest.raises(ReferenceDataError, match="other than"):
        verify_reference_digests({prefix: str(other)}, [12])


def test_the_matching_genome_passes_at_the_request_layer(tmp_path, genome):
    from neoswga.core.reference_check import verify_reference_digests

    prefix = str(tmp_path / "current")
    write_index(prefix, genome=str(genome))

    verify_reference_digests({prefix: str(genome)}, [12])


def test_a_prefix_the_manifest_does_not_name_is_not_checked(tmp_path, genome):
    """An unknown expectation cannot be compared, so it is not invented."""
    from neoswga.core.reference_check import verify_reference_digests

    prefix = str(tmp_path / "unnamed")
    write_index(prefix, genome=str(genome))

    verify_reference_digests({}, [12])


def test_an_index_with_no_recorded_digest_is_unknown_not_matching(tmp_path, genome):
    from neoswga.core.reference_check import verify_reference_digests

    prefix = str(tmp_path / "nodigest")
    write_index(prefix, genome=None)

    with pytest.raises(ReferenceDataError, match="no recorded reference digest"):
        verify_reference_digests({prefix: str(genome)}, [12])
