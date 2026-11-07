"""Every documentation reference to a Markdown file must point at one.

Three stale entries and a dead link spanning five documents were found by hand
in October 2026, all incidentally, across 124 files and ~44,500 lines. One of
the stale entries had been falsified by a commit made two days earlier in the
same session. Correctness of prose there depends on someone happening to read
the right line, which is what this test replaces.

Two shapes are checked, because the dead reference took both and a link-only
check finds one of them:

- a Markdown link, ``[text](path/to/file.md)``
- an inline code span naming a path, ``` `docs/path/to/file.md` ```

Four of the five documents citing
`docs/superpowers/plans/2026-09-15-condition-aware-pool-design.md` -- a plan
that never existed in git -- used the second shape, so the sweep that found the
first missed them.

**Plan, spec and archive documents are exempt from the code-span rule.** A plan
that says "create `docs/validation/model_domain_review.md`" is specifying a
deliverable, and a completed cleanup plan that lists
`tests/integration/test_plan.md` as "Outdated" is recording a file it had
deleted. Neither is rot; treating them as rot would mean a plan could not name
what it proposes. Markdown LINKS are still checked everywhere, including in
plans, because a link promises something a reader can follow now.

**Orphans are deliberately NOT an error.** `docs/README.md` states that
`validation/` holds one record per claim and that its curated table lists only
the ones most often cited, so an unlinked record is the documented norm rather
than a defect.
"""

import pathlib
import re

import pytest

_ROOT = pathlib.Path(__file__).resolve().parent.parent

#: ``[text](target)``, target captured. Angle-bracket and title forms are not
#: used anywhere in this tree, so the simple form is enough.
_LINK = re.compile(r"\[[^\]]*\]\(([^)\s]+)")

#: A single-backtick code span whose whole content looks like a path to a
#: Markdown file. Anchored on the extension so prose in backticks is not swept
#: in, and it requires a separator so a bare ``README.md`` -- which is ambiguous
#: about which directory it means -- is left alone.
_CODE_PATH = re.compile(r"`([A-Za-z0-9_./-]+/[A-Za-z0-9_.-]+\.md)`")

#: Directories whose documents may name files that do not exist: see the module
#: docstring. Scoped to directories rather than to individual files so a new
#: plan is covered without an edit here.
_FORWARD_REFERENCE_DIRS = ("superpowers", "design", "archive")


def _markdown_files():
    return sorted(
        p for p in _ROOT.rglob("*.md") if ".git" not in p.parts and "node_modules" not in p.parts
    )


def _resolve(source: pathlib.Path, target: str) -> pathlib.Path | None:
    """Where `target`, written in `source`, should point.

    A code-span path is resolved relative to the repository root as well as to
    the file's own directory, because both conventions are in use: records cite
    `docs/validation/x.md` from inside `docs/validation/`.
    """
    for candidate in ((source.parent / target), (_ROOT / target)):
        if candidate.exists():
            return candidate
    return None


def _references(path: pathlib.Path):
    text = path.read_text(encoding="utf-8", errors="replace")
    forward_ok = any(part in _FORWARD_REFERENCE_DIRS for part in path.parts)
    for match in _LINK.finditer(text):
        target = match.group(1).split("#")[0].strip()
        if target and not target.startswith(("http://", "https://", "mailto:")):
            if target.endswith(".md"):
                yield target
    if forward_ok:
        return
    for match in _CODE_PATH.finditer(text):
        target = match.group(1)
        # A literal ellipsis is prose, not a path: "docs/...md" appears in a
        # sentence describing the layout.
        if "..." not in target:
            yield target


@pytest.mark.parametrize("path", _markdown_files(), ids=lambda p: str(p.relative_to(_ROOT)))
def test_markdown_references_resolve(path):
    broken = sorted({target for target in _references(path) if _resolve(path, target) is None})
    assert not broken, (
        f"{path.relative_to(_ROOT)} references Markdown files that do not "
        f"exist: {broken}. Either fix the path or, if the document was never "
        f"written, name it in prose instead of citing a path."
    )


def test_the_sweep_sees_both_reference_shapes():
    """Guard the regexes themselves.

    A test whose patterns silently stop matching passes forever, and the dead
    reference that survived in five documents was in the code-span shape.
    """
    index = _ROOT / "docs" / "README.md"
    assert any(
        t.endswith("LIMITATIONS.md") for t in _references(index)
    ), "the link form matched no .md target in docs/README.md"

    sample = "see `docs/validation/tool_comparison.md` for the table"
    assert _CODE_PATH.findall(sample) == [
        "docs/validation/tool_comparison.md"
    ], "the code-span form no longer matches a path in backticks"
    # A bare filename in backticks is intentionally not swept in: it does not
    # say which directory it means.
    assert _CODE_PATH.findall("the `README.md` at the root") == []
