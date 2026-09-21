"""A params.json a reader copies out of the docs must not warn.

The user guide shipped a sample configuration containing `top_sets_count`,
which was never a schema key, and `min_amp_pred`, which is real but inert since
the amplification gate left the default path on 2026-09-05. Anyone who copied
that block got an "unrecognised parameter" warning on every run and a gate that
did nothing, and nothing in the suite could see it because a documented
configuration is never executed.

That is the Known Issue 8 class arriving by a route none of its five ratchets
watch. Those check that a schema key binds something, that a CLI flag is read,
and that a capability is reachable. None asks whether the configuration this
project TELLS people to write is one it accepts.

So this walks every JSON block in the user-facing guides that looks like a
params file and puts it through the same validator `neoswga validate params`
uses. Unknown keys fail. Genome paths are not checked, because a documented
example names files that do not exist here and should.
"""

import json
import pathlib
import re

import pytest

DOCS = pathlib.Path(__file__).resolve().parent.parent / "docs"

#: A block is a params file if it carries one of these. Grid files, report
#: fragments and API payloads share the fence and are not configurations.
MARKERS = ("fg_genomes", "fg_prefixes", "polymerase", "min_k")


def params_blocks():
    """(path, index, mapping) for every documented params configuration."""
    roots = [DOCS / "guides", DOCS]
    seen = set()
    for root in roots:
        for path in sorted(root.glob("*.md")):
            if path in seen or "archive" in path.parts or "validation" in path.parts:
                continue
            seen.add(path)
            text = path.read_text(encoding="utf-8")
            for index, block in enumerate(re.findall(r"```json\n(.*?)```", text, re.S)):
                try:
                    loaded = json.loads(block)
                except ValueError:
                    # A deliberately elided example ("..."), not a config.
                    continue
                if isinstance(loaded, dict) and any(m in loaded for m in MARKERS):
                    yield path, index, loaded


def test_the_docs_actually_contain_configurations():
    """Guard the guard: a broken extractor would make every check vacuous."""
    found = list(params_blocks())

    assert found, "no documented params configuration was found to check"


@pytest.mark.parametrize(
    "path,index,config",
    [pytest.param(p, i, c, id=f"{p.name}#{i}") for p, i, c in params_blocks()],
)
def test_a_documented_config_uses_only_current_keys(path, index, config):
    """Retired and misspelled keys are what this catches.

    Checked against the schema directly rather than by running the validator,
    so the failure names the key rather than a process exit status.
    """
    schema = json.loads(
        (
            pathlib.Path(__file__).resolve().parent.parent
            / "neoswga"
            / "core"
            / "schema"
            / "params.schema.json"
        ).read_text(encoding="utf-8")
    )
    known = set(schema["properties"])

    # A leading underscore marks a comment and is accepted by the loader.
    unknown = sorted(k for k in config if not k.startswith("_") and k not in known)

    assert not unknown, (
        f"{path.relative_to(DOCS.parent)} block {index} documents key(s) the "
        f"schema does not define: {unknown}. A reader copying this block gets "
        "an 'unrecognised parameter' warning."
    )


@pytest.mark.parametrize(
    "path,index,config",
    [pytest.param(p, i, c, id=f"{p.name}#{i}") for p, i, c in params_blocks()],
)
def test_a_documented_config_does_not_set_a_retired_gate(path, index, config):
    """`min_amp_pred` is in the schema and does nothing without `--amp-model`.

    Setting it in an example teaches a filter that is not applied. The schema
    cannot refuse it -- the key is real and `--amp-model` restores it -- so the
    rule belongs here, where "what the docs recommend" is the subject.
    """
    assert "min_amp_pred" not in config, (
        f"{path.relative_to(DOCS.parent)} block {index} sets min_amp_pred. The "
        "amplification gate was retired from the default path on 2026-09-05 "
        "and runs only under --amp-model, so this configures nothing. Mention "
        "it in prose with that caveat instead."
    )
