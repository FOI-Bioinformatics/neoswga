"""A new design must not run on concatenation-era index geometry.

Task 2 of the condition-aware pool design plan. Record starts were added to the
position index on 2026-09-14 so coverage windows stop at record edges and
k-mers are not matched across joins. An index written before that carries no
record geometry, and nothing stopped a new design from using one: windows would
reach across contigs again, silently, on an index that looks complete.

The same applies to a stale index. The provenance sidecar already detects a
reference whose content changed under the same filename, but that check governs
rescanning, not whether a pool design may proceed.

Reading is deliberately unaffected. `report-pool` renders a saved plan without
touching an index at all, and must keep doing so: a historical report describes
what was computed then, and refusing to display it would not make it more
accurate.
"""

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

from neoswga.core import string_search
from neoswga.core.exceptions import ReferenceDataError
from neoswga.core.position_cache import PositionCache

PRIMER = "ACGTACGTACGT"


def _write_index(prefix, *, record_starts, genome=None):
    """An index written the modern way, or the old way when starts are None."""
    with h5py.File(f"{prefix}_12mer_positions.h5", "w"):
        pass
    string_search.write_to_h5py(
        {PRIMER: [5, 70]}, prefix, record_starts=record_starts, genome_fname=genome
    )


@pytest.fixture
def two_records(tmp_path):
    fasta = tmp_path / "two.fna"
    fasta.write_text(">one\n" + "A" * 60 + "\n>two\n" + "C" * 60 + "\n")
    string_search.clear_genome_cache()
    return fasta


def test_a_modern_index_satisfies_the_requirement(tmp_path, two_records):
    prefix = str(tmp_path / "modern")
    _write_index(prefix, record_starts=[0, 60], genome=str(two_records))

    PositionCache([prefix], [PRIMER]).require_record_metadata([prefix])


def test_an_index_without_record_geometry_is_refused(tmp_path):
    """The concatenation-era case this exists for.

    The decision moved layers on 2026-09-21 and the refusal did not go away.
    Whether an index NEEDS record starts depends on how many records its
    reference holds -- a single-record reference has no joins for a window to
    cross -- and only the resolved request pairs a prefix with a genome.
    `PositionCache` cannot answer that without reading a mutable global, which
    is the defect that made a design refuse its own index under `pytest -n 8`.
    """
    from neoswga.core.reference_check import verify_index_geometry

    prefix = str(tmp_path / "old")
    _write_index(prefix, record_starts=None)
    fasta = tmp_path / "joined.fna"
    fasta.write_text(">a\n" + "A" * 60 + "\n>b\n" + "C" * 60 + "\n")

    with pytest.raises(ReferenceDataError, match="record geometry") as excinfo:
        verify_index_geometry({prefix: str(fasta)}, [12])

    message = str(excinfo.value)
    assert "old" in message, "the message must name the prefix at fault"
    assert (
        "count-kmers" in message or "filter" in message
    ), "the message must name the command that regenerates it"


def test_every_offending_prefix_is_named_at_once(tmp_path, two_records):
    """A user regenerating indexes should not discover them one run at a time."""
    from neoswga.core.reference_check import verify_index_geometry

    good, bad_one, bad_two = (str(tmp_path / n) for n in ("good", "bad1", "bad2"))
    _write_index(good, record_starts=[0, 60], genome=str(two_records))
    _write_index(bad_one, record_starts=None)
    _write_index(bad_two, record_starts=None)

    manifest = {p: str(two_records) for p in (good, bad_one, bad_two)}
    with pytest.raises(ReferenceDataError) as excinfo:
        verify_index_geometry(manifest, [12])

    message = str(excinfo.value)
    assert "bad1" in message and "bad2" in message
    assert "good" not in message


def test_a_missing_index_is_refused_rather_than_passed(tmp_path):
    cache = PositionCache([str(tmp_path / "absent")], [PRIMER])
    with pytest.raises(ReferenceDataError):
        cache.require_record_metadata([str(tmp_path / "absent")])


def test_a_single_record_reference_is_accepted(tmp_path):
    """One record is legitimate geometry, not missing geometry."""
    fasta = tmp_path / "one.fna"
    fasta.write_text(">only\n" + "ACGT" * 30 + "\n")
    string_search.clear_genome_cache()
    prefix = str(tmp_path / "single")
    _write_index(prefix, record_starts=[0], genome=str(fasta))

    PositionCache([prefix], [PRIMER]).require_record_metadata([prefix])


def test_the_index_records_its_format_version_and_reference(tmp_path, two_records):
    prefix = str(tmp_path / "meta")
    _write_index(prefix, record_starts=[0, 60], genome=str(two_records))

    with h5py.File(f"{prefix}_12mer_positions.h5", "r") as db:
        assert db.attrs["index_format_version"] >= 1
        assert db.attrs["reference_digest"]


def test_plan_pool_validates_geometry_before_evaluating_a_panel():
    """The requirement has to be invoked, not merely available.

    A capability nothing calls is the failure this project has recorded before,
    when a background genome was read and never acted on.
    """
    import ast
    import pathlib

    source = pathlib.Path("neoswga/cli/plan_pool.py").read_text()
    tree = ast.parse(source)
    calls = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == "require_record_metadata"
    ]
    assert calls, "plan-pool builds a PositionCache without validating its geometry"


def test_report_pool_still_renders_without_any_index(tmp_path):
    """Reading a saved plan must not require an index at all."""
    import ast
    import pathlib

    source = pathlib.Path("neoswga/cli/plan_pool.py").read_text()
    tree = ast.parse(source)
    report = next(
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.FunctionDef) and node.name == "run_report_pool"
    )
    names = {
        node.func.id
        for node in ast.walk(report)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    assert (
        "PositionCache" not in names
    ), "rendering a historical report must not open a position index"
