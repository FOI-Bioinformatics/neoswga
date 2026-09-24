"""Advice that names a real thing must also name a thing that works.

`test_advice_names_something_that_exists.py` is the mechanical floor: it
checks that a named flag, command or script exists. It cannot check whether
following the advice resolves the stated problem, and that is the failure that
has actually cost something in this repository twice.

Each test here pins one case where it did not, found 2026-09-24.
"""

import ast
import json
import pathlib
from dataclasses import fields

import pytest

REPO = pathlib.Path(__file__).resolve().parent.parent


def test_the_qa_failure_does_not_offer_a_knob_that_does_not_exist():
    """It said "Relax the QA stringency or drop --enable-qa". Only the second
    half was possible.

    All ten `QAFilterConfig` fields are neither a schema key nor a `parameter`
    global nor a CLI flag, and the only production caller passes no config, so
    every run uses the "moderate" defaults. The library HAS the lever --
    `create_three_prime_analyzer(stringency="lenient")` -- and nothing can ask
    for it, which is Known Issue 16's shape: a parameter no caller supplies.
    """
    from neoswga.core import parameter
    from neoswga.core.pipeline_qa_integration import QAFilterConfig

    schema = set(
        json.loads((REPO / "neoswga/core/schema/params.schema.json").read_text())["properties"]
    )
    reachable = [
        f.name for f in fields(QAFilterConfig) if f.name in schema or hasattr(parameter, f.name)
    ]

    source = (REPO / "neoswga/core/pipeline_qa_integration.py").read_text()
    if reachable:
        pytest.fail(
            f"{reachable} became configurable; the message may now offer it, "
            "but say which key rather than 'relax the QA stringency'"
        )
    assert "Relax the QA stringency or drop" not in source
    assert "not configurable" in source, "it must say the stringency cannot be changed"


def test_the_bloom_advice_names_what_completes_it():
    """ "Consider pre-building a Bloom filter: neoswga build-filter ..." cost
    hours on hg38 and changed nothing: `filter` never looks in the output
    directory, and setting `bloom_filter_path` alone does not enable it
    either."""
    source = (REPO / "neoswga/core/pipeline.py").read_text()
    index = source.index("Bloom filter")
    window = source[index - 400 : index + 700]

    assert "bloom_filter_path" in window, "name the key that points at the file"
    assert "use_bloom_filter" in window, "name the key that switches it on"


def test_a_path_without_the_switch_says_so_rather_than_counting_exactly():
    """The user who built the filter and set only the path. They used to get
    exact counting and silence, having paid the build cost."""
    source = (REPO / "neoswga/core/filter.py").read_text()

    # Only the warning is asserted here. Whether the dead key is still read is
    # a separate question and `test_the_dead_bloom_key_is_gone` answers it by
    # AST, which does not depend on how the call happens to be formatted.
    assert "is NOT being used" in source, "it must warn when the path is set and the switch is not"


def test_the_dead_bloom_key_is_gone():
    """`bg_bloom` was read by `filter` and is not a schema key, is assigned
    nowhere, and so could never be anything but None. The comment above it
    claimed it auto-enabled the filter."""
    schema = set(
        json.loads((REPO / "neoswga/core/schema/params.schema.json").read_text())["properties"]
    )
    assert "bg_bloom" not in schema, "if it becomes real, this test should change"

    tree = ast.parse((REPO / "neoswga/core/filter.py").read_text())
    reads = [
        n
        for n in ast.walk(tree)
        if isinstance(n, ast.Call)
        and isinstance(n.func, ast.Name)
        and n.func.id == "getattr"
        and len(n.args) >= 2
        and isinstance(n.args[1], ast.Constant)
        and n.args[1].value == "bg_bloom"
    ]

    assert not reads, "filter.py still reads a key nothing can set"


def test_an_optional_dependency_is_not_called_a_corrupted_installation():
    """The handler wraps a whole pipeline step, so it catches an ImportError
    from any optional dependency -- networkx, pysam, pybloom_live, matplotlib
    -- each of which already says what to install. It then added "This may
    indicate a corrupted installation. Try: pip install -e . --force-reinstall"
    underneath: correct advice followed by wrong advice, wrong one last.
    Reinstalling does not add an extra."""
    import logging

    from neoswga.cli._failure import report_import_failure

    def messages(module_name):
        records = []
        handler = logging.Handler()
        handler.emit = records.append
        log = logging.getLogger("neoswga.cli._failure")
        log.addHandler(handler)
        try:
            try:
                __import__(module_name)
            except ImportError as exc:
                report_import_failure(exc)
        finally:
            log.removeHandler(handler)
        return " ".join(r.getMessage() for r in records)

    third_party = messages("networkx_definitely_not_installed")
    assert "optional dependency" in third_party
    assert "force-reinstall" not in third_party, "a missing extra is not a corrupt install"

    our_own = messages("neoswga.core.not_a_real_module")
    assert "force-reinstall" in our_own, "a genuinely broken install must still say so"
