"""A run manifest that cannot be parsed is kept, not erased.

`write_manifest` read the existing file, and when it did not parse it started a
fresh `{"steps": []}` and truncated the file. Every earlier step's entry,
including the `effective_conditions` `export` and `report` read, was lost
without a word (2026-11-23 sweep, tier 2).
"""

import json

import pytest

from neoswga.core.run_manifest import write_manifest


@pytest.mark.parametrize("content", ["{ broken", "[]", '{"steps": "not a list"}'])
def test_an_unusable_manifest_is_preserved(tmp_path, content):
    (tmp_path / "run_manifest.json").write_text(content)
    write_manifest("filter", str(tmp_path))

    kept = list(tmp_path.glob("run_manifest.json.corrupt-*"))
    assert len(kept) == 1
    assert kept[0].read_text() == content
    fresh = json.loads((tmp_path / "run_manifest.json").read_text())
    assert [step["step"] for step in fresh["steps"]] == ["filter"]


def test_a_good_manifest_is_appended_without_a_copy(tmp_path):
    write_manifest("filter", str(tmp_path))
    write_manifest("optimize", str(tmp_path))
    assert not list(tmp_path.glob("run_manifest.json.corrupt-*"))
    steps = json.loads((tmp_path / "run_manifest.json").read_text())["steps"]
    assert [step["step"] for step in steps] == ["filter", "optimize"]
