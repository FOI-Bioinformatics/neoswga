"""Each advertised capability names its evidence, and the evidence exists.

`docs/EVIDENCE.md` places every capability this tool advertises on one of five
ordered tiers, from "the code exists" to "a pool was synthesised and
sequenced". It is generated from `docs/capability_evidence.json` so that a
claim and the evidence for it cannot drift apart in prose -- which is the
failure `docs/SCIENCE_CITATIONS.md` already has, stating Klenow processivity as
10,000 bp while the shipped registry says 40.

What these tests check is narrow and worth stating plainly. They cannot verify
that a measurement supports the tier it is cited for; a human read it. They CAN
verify the three things that make the page decay silently:

- every file named exists, so the evidence does not 404;
- the tier is consistent with what is named beside it, so nothing claims to be
  measured while citing no measurement;
- the page matches the metadata, so editing the prose alone fails.

And one claim is asserted outright: **nothing is at the prospective tier.** If
that ever changes it must be a deliberate edit to this test, not to a JSON
file, because it is the strongest claim the project could make and the one
`docs/LIMITATIONS.md` exists to deny.
"""

import json
import pathlib
import subprocess
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
METADATA = ROOT / "docs" / "capability_evidence.json"
PAGE = ROOT / "docs" / "EVIDENCE.md"

#: Weakest first. A tier's requirements are cumulative going up.
TIERS = [
    "implemented",
    "connected",
    "oracle_tested",
    "retrospectively_evaluated",
    "prospectively_validated",
]


@pytest.fixture(scope="module")
def data():
    return json.loads(METADATA.read_text())


@pytest.fixture(scope="module")
def capabilities(data):
    return data["capabilities"]


# ---------------------------------------------------------------------------
# The evidence exists
# ---------------------------------------------------------------------------


def test_every_named_test_file_exists(capabilities):
    """A capability citing a test that has been renamed cites nothing."""
    missing = [
        f"{item['capability']} -> {path}"
        for item in capabilities
        for path in item["tests"]
        if not (ROOT / path).exists()
    ]

    assert not missing, missing


def test_every_named_measurement_exists(capabilities):
    missing = [
        f"{item['capability']} -> {path}"
        for item in capabilities
        for path in item["measurements"]
        if not (ROOT / path).exists()
    ]

    assert not missing, missing


# ---------------------------------------------------------------------------
# The tier is consistent with what is named beside it
# ---------------------------------------------------------------------------


def test_every_tier_is_one_of_the_five(capabilities, data):
    declared = [tier["name"] for tier in data["tiers"]]

    assert declared == TIERS, "the declared tiers no longer match this test's order"
    for item in capabilities:
        assert item["tier"] in TIERS, f"{item['capability']} sits on tier {item['tier']}"


def test_a_capability_claiming_measurement_cites_one(capabilities):
    """The check that stops the tiers becoming adjectives."""
    unsupported = [
        item["capability"]
        for item in capabilities
        if TIERS.index(item["tier"]) >= TIERS.index("retrospectively_evaluated")
        and not item["measurements"]
    ]

    assert not unsupported, (
        "these are placed at or above 'retrospectively evaluated' and name no "
        f"measurement: {unsupported}"
    )


def test_a_capability_claiming_a_test_has_one(capabilities):
    untested = [
        item["capability"]
        for item in capabilities
        if TIERS.index(item["tier"]) >= TIERS.index("oracle_tested") and not item["tests"]
    ]

    assert not untested, untested


def test_every_capability_says_what_reaches_it(capabilities):
    """`connected` is a tier of its own because reachability is not an
    assumption here: six capabilities were built, tested and reachable by no
    command at all, found in one audit."""
    for item in capabilities:
        assert item["reached_by"].strip(), f"{item['capability']} names no command"
        assert item["note"].strip(), f"{item['capability']} carries no caveat"


def test_nothing_claims_to_have_been_validated_prospectively(capabilities):
    """The strongest claim the project could make, and it cannot make it.

    Asserted here rather than left to the metadata on purpose: promoting
    something to this tier must be a deliberate edit to a test, because it
    contradicts `docs/LIMITATIONS.md` and every measurement's caveats.
    """
    claimed = [
        item["capability"] for item in capabilities if item["tier"] == "prospectively_validated"
    ]

    assert not claimed, (
        f"{claimed} claim prospective validation. No pool this tool designed has "
        "been synthesised, run and sequenced, and there is no BAM or CRAM in "
        "this repository to evaluate one against."
    )


# ---------------------------------------------------------------------------
# The page matches the metadata
# ---------------------------------------------------------------------------


def test_the_page_is_not_stale():
    """Editing the prose alone must fail. That is the whole reason to generate it."""
    completed = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "generate_evidence_page.py"), "--check"],
        capture_output=True,
        text=True,
        cwd=ROOT,
    )

    assert completed.returncode == 0, (
        "docs/EVIDENCE.md no longer matches docs/capability_evidence.json. "
        "Edit the JSON and re-run scripts/generate_evidence_page.py.\n" + completed.stderr
    )


def test_the_page_leads_with_what_is_not_established():
    text = PAGE.read_text()

    assert "prospectively validated" in text
    assert "LIMITATIONS.md" in text, "the page must point at the limitations it summarises"


def test_the_readme_points_at_it():
    assert "EVIDENCE.md" in (ROOT / "README.md").read_text()
