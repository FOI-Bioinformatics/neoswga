"""The resolved request must reach the stage that would use it.

`cli/pipeline.py` resolves a `DesignRequest` before the search begins, logs its
hash and runs the two reference checks against it, then drops it: the local was
named `_request`. So `run_optimization` had no request, and
`panel_evaluation.evaluate_panel`, which takes one, had no production caller at
all -- only its own two test files imported it.

Nothing reads the request yet, so this changes no delivered panel and no
verdict. What it removes is the reason the assessment could not be wired.

The end-to-end test is the one that matters. `attach_search_config` had two
tests of the weaker kind -- one asserting by AST that an attribute of that name
was assigned, one asserting by source text that something read one -- and the
refinement still received None on every run, because both ends existed and the
path between them did not.
"""

import ast
import json
from pathlib import Path

import pytest

from tests.conftest import plasmid_example_ready

PIPELINE = Path(__file__).resolve().parent.parent / "neoswga" / "cli" / "pipeline.py"


def _call_names(node):
    func = node.func
    return getattr(func, "id", None) or getattr(func, "attr", None)


def test_the_command_hands_the_request_to_the_optimizer():
    """A source check on the seam, so a broken wiring names the seam.

    Not sufficient alone -- the end-to-end test below is what says the object
    arrives -- but without it a broken seam surfaces only as an assessment
    silently absent from an artifact, which reads like a feature that was never
    added.
    """
    tree = ast.parse(PIPELINE.read_text())
    forwarded = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and _call_names(node) == "optimize_step4"
        and any(keyword.arg == "design_request" for keyword in node.keywords)
    ]

    assert forwarded, "optimize_step4 is called without design_request"


@pytest.mark.skipif(not plasmid_example_ready(), reason="needs jellyfish and the plasmid example")
def test_the_request_arrives_with_the_hash_the_run_resolved(plasmid_run, monkeypatch):
    """The object itself, through the real command, named by the same hash.

    The hash is the assertion rather than "a request arrived": a request built
    from different parameters would satisfy the weaker claim while describing a
    different design.
    """
    import neoswga.core.unified_optimizer as unified
    from neoswga.cli import pipeline
    from neoswga.cli_unified import create_parser
    from neoswga.core.design_request import resolve_design_request

    params_path = plasmid_run / "params.json"
    expected = resolve_design_request(json.loads(params_path.read_text())).request_hash

    captured = {}
    original = unified.run_optimization

    def recording(*args, **kwargs):
        captured["request"] = kwargs.get("design_request")
        return original(*args, **kwargs)

    monkeypatch.setattr(unified, "run_optimization", recording)
    monkeypatch.chdir(plasmid_run)

    # The real parser, so the namespace carries every flag the handler reads
    # and a flag added later cannot make this test pass by being absent.
    args = create_parser().parse_args(["optimize", "-j", "params.json"])
    pipeline.run_step4(args)

    assert captured["request"] is not None, "run_optimization was given no request"
    assert captured["request"].request_hash == expected
