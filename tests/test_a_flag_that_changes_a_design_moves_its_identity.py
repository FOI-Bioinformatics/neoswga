"""A flag that changes the delivered panel must change `request_hash`.

Review finding F1 of Phase 6c. `min_per_target_coverage` became the first panel
limit with a CLI flag, and `--min-per-target-coverage` steered the search while
`resolve_design_request` built `DesignRequest.constraints` from the params FILE
alone. Measured on the packaged plasmid example, in ONE work directory so the
hashed path fields are identical:

    noflag    hash a00a4b95…  primers CCCATTGACG, AAACCTGTCG, …  ok True
    flag 0.5  hash a00a4b95…  primers CCATTGACG, GGCAACTA, …     ok False

Two runs, one hash, different primers and opposite qualification verdicts, from
a field whose own docstring is "SHA-256 over a canonical form of every setting
that affects a design". `run_manifest.json` had the same hole: `resolved_params`
is a verbatim copy of the file, so the limit that chose the flag run's primers
appeared nowhere in the record of it.

The fix is the mechanism rather than the one key: `design_request.cli_overrides`
collects every flag that names a params.json key and carries `None` as its
argparse default, and both request entry points overlay it onto the file.
`_CLI_OVERLAY_KEYS` and `CLI_OVERLAY_EXCLUDED` between them must account for
every such flag, which is what the ratchet below enforces, so the next option
with this shape cannot be forgotten.

After the fix, the same three runs:

    noflag    hash a00a4b95…  ok True
    flag 0.5  hash b8159951…  ok False   cli_overrides names the limit
    json 0.5  hash b8159951…  ok False   same effective request, same identity
"""

import argparse
import ast
import inspect
import json
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
SCHEMA = json.loads((ROOT / "neoswga" / "core" / "schema" / "params.schema.json").read_text())


def _optimize_parser():
    """The real `optimize` subparser, not a stand-in for it."""
    from neoswga.cli_unified import create_parser

    parser = create_parser()
    actions = [a for a in parser._actions if isinstance(a, argparse._SubParsersAction)]
    return actions[0].choices["optimize"]


class _Args:
    """An argparse namespace carrying only what a run set."""

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


class TestEveryFlagWithThisShapeIsAccountedFor:
    """The ratchet. A flag naming a params.json key, defaulting to `None`, is a
    flag that means "this invocation asked for it" -- so it either overlays or
    is listed with a reason."""

    def test_no_such_flag_is_unaccounted_for(self):
        from neoswga.core.design_request import (
            _CLI_OVERLAY_KEYS,
            _REQUEST_ONLY_KEYS,
            CLI_OVERLAY_EXCLUDED,
        )

        known = set(_CLI_OVERLAY_KEYS) | set(CLI_OVERLAY_EXCLUDED)
        unaccounted = sorted(
            action.dest
            for action in _optimize_parser()._actions
            if action.dest != "help"
            and action.default is None
            and (action.dest in SCHEMA["properties"] or action.dest in _REQUEST_ONLY_KEYS)
            and action.dest not in known
        )

        assert not unaccounted, (
            "these optimize flags name a params.json key and default to None, so "
            "they override the file; add each to _CLI_OVERLAY_KEYS or to "
            f"CLI_OVERLAY_EXCLUDED with its reason: {unaccounted}"
        )

    def test_the_exclusions_are_not_stale(self):
        """An exclusion is a promise about a real flag."""
        from neoswga.core.design_request import CLI_OVERLAY_EXCLUDED

        dests = {action.dest for action in _optimize_parser()._actions}

        assert set(CLI_OVERLAY_EXCLUDED) <= dests

    def test_the_overlay_keys_are_all_declared_somewhere(self):
        from neoswga.core.design_request import _CLI_OVERLAY_KEYS, _REQUEST_ONLY_KEYS

        for key in _CLI_OVERLAY_KEYS:
            assert key in SCHEMA["properties"] or key in _REQUEST_ONLY_KEYS, key

    def test_the_per_target_floor_is_one_of_them(self):
        from neoswga.core.design_request import _CLI_OVERLAY_KEYS

        assert "min_per_target_coverage" in _CLI_OVERLAY_KEYS


class TestCollectingWhatTheInvocationSet:
    def test_an_absent_flag_contributes_nothing(self):
        """The Known Issue 8 rule, applied to identity: no flag, no override."""
        from neoswga.core.design_request import cli_overrides

        assert cli_overrides(_Args(min_per_target_coverage=None, polymerase=None)) == {}

    def test_an_empty_namespace_contributes_nothing(self):
        from neoswga.core.design_request import cli_overrides

        assert cli_overrides(_Args()) == {}

    def test_a_set_flag_is_collected_under_its_params_json_name(self):
        from neoswga.core.design_request import cli_overrides

        assert cli_overrides(_Args(min_per_target_coverage=0.5)) == {"min_per_target_coverage": 0.5}

    def test_the_size_flag_is_not_collected_here(self):
        """`effective_design_request` overlays the effective size afterwards,
        under both its names, and must win."""
        from neoswga.core.design_request import cli_overrides

        assert cli_overrides(_Args(num_primers=40)) == {}


class TestTheRequestCarriesTheFlag:
    """The defect itself: the hash must move, and two spellings of the same
    effective request must agree."""

    @pytest.fixture
    def params_file(self, tmp_path):
        params = {
            "schema_version": 2,
            "fg_prefixes": [str(tmp_path / "t")],
            "fg_seq_lengths": [100000],
            "fg_genomes": [str(tmp_path / "t.fna")],
            "polymerase": "phi29",
            "reaction_temp": 30.0,
            "num_primers": 6,
            "target_set_size": 6,
        }
        path = tmp_path / "params.json"
        path.write_text(json.dumps(params))
        return path

    def test_the_flag_moves_the_request_hash(self, params_file):
        from neoswga.core.design_request import design_request_for_run

        without = design_request_for_run(_Args(json_file=str(params_file)), None)
        with_flag = design_request_for_run(
            _Args(json_file=str(params_file), min_per_target_coverage=0.5), None
        )

        assert without.request_hash != with_flag.request_hash

    def test_the_flag_reaches_the_constraints_on_the_request(self, params_file):
        from neoswga.core.design_request import design_request_for_run

        request = design_request_for_run(
            _Args(json_file=str(params_file), min_per_target_coverage=0.5), None
        )

        assert request.constraints is not None
        assert request.constraints.min_per_target_coverage == 0.5

    def test_the_flag_and_the_file_resolve_to_one_identity(self, params_file):
        """Same effective request, same hash: the invariant the fix establishes
        in both directions."""
        from neoswga.core.design_request import design_request_for_run

        params = json.loads(params_file.read_text())
        params["min_per_target_coverage"] = 0.5
        in_file = params_file.parent / "with_floor.json"
        in_file.write_text(json.dumps(params))

        by_flag = design_request_for_run(
            _Args(json_file=str(params_file), min_per_target_coverage=0.5), None
        )
        by_file = design_request_for_run(_Args(json_file=str(in_file)), None)

        assert by_flag.request_hash == by_file.request_hash

    def test_the_effective_request_carries_it_too(self, params_file):
        """`effective_request_for_manifest` goes through this one, and it is the
        request whose hash the manifest records."""
        from neoswga.core.design_request import effective_design_request

        without = effective_design_request(_Args(json_file=str(params_file)), target_size=6)
        with_flag = effective_design_request(
            _Args(json_file=str(params_file), min_per_target_coverage=0.5), target_size=6
        )

        assert without.request_hash != with_flag.request_hash

    def test_the_size_overlay_still_wins_over_the_flag(self, params_file):
        from neoswga.core.design_request import effective_design_request

        request = effective_design_request(
            _Args(json_file=str(params_file), num_primers=40), target_size=12
        )

        assert request.target_size == 12

    def test_the_manifest_records_the_overrides(self):
        """`resolved_params` is a verbatim copy of the file by design, so the
        overrides ride in `extra` beside `effective_set_size`."""
        text = (ROOT / "neoswga" / "cli" / "pipeline.py").read_text()

        assert '"cli_overrides": _cli_overrides(args)' in text


class TestTheThreadingIntoTheSearchIsPinned:
    """Review finding F3: removing `kwargs` from the two `_dispatch_optimizer`
    calls reverted the delivered panel while 294 tests across ten option and
    ratchet files still passed. The flag stopped steering the search and nothing
    said so, because the "Configured limits" table still printed."""

    def test_every_call_in_the_optimize_path_passes_the_run_flags(self):
        from neoswga.core import unified_optimizer

        tree = ast.parse(inspect.getsource(unified_optimizer))
        calls = [
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "constraints_from_parameter"
        ]

        assert calls, "constraints_from_parameter is no longer called here"
        bare = [node.lineno for node in calls if len(node.args) + len(node.keywords) < 2]
        assert not bare, (
            "constraints_from_parameter is called with one argument at "
            f"unified_optimizer.py:{bare}, so the run's flags do not reach the "
            "limits built there and --min-per-target-coverage stops steering "
            "the search silently"
        )

    def test_all_six_call_sites_are_still_there(self):
        """A site deleted rather than un-threaded would pass the test above."""
        from neoswga.core import unified_optimizer

        tree = ast.parse(inspect.getsource(unified_optimizer))
        calls = [
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "constraints_from_parameter"
        ]

        assert len(calls) == 6, f"expected six call sites, found {len(calls)}"
