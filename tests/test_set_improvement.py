"""`improve-set`: diagnose an existing oligo set and propose small edits to it.

Phase 2b of `docs/superpowers/plans/2026-10-02-genomic-diversity-and-multi-host.md`.

Each test here is one of the plan's synthetic verification items, and each pins
a way this kind of command goes wrong.

**The pooled figure hides the target that needs help.** With a hole on one
target of two, the candidate that raises the pooled coverage most is the wrong
one. The proposal must be the one that fills the hole.

**A proposal that was not checked bounds nothing.** A proposed add is asserted
clean against every oligo that stays, with `dimer.is_dimer_fast` as the oracle
and over several seeds, on the proposals themselves.

**The prediction must be the measurement.** Applying the top proposal and
evaluating the result afresh reproduces the predicted figures exactly. If it
does not, two code paths disagree and that is the defect.

**Unknown is not zero.** An oligo whose sites could not be established on a
target is not "binding nowhere", and must not be proposed for dropping as if it
were.
"""

from __future__ import annotations

import random

import pytest

from neoswga.core.dimer import is_dimer_fast
from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.reference_panel_evaluation import (
    ROLE_HOST,
    ROLE_TARGET,
    ReferenceSpec,
    evaluate_reference_panel,
)
from neoswga.core.set_improvement import (
    KIND_ADD,
    KIND_DROP,
    KIND_SWAP,
    KIND_TRADE_OFF,
    ImprovementSettings,
    improve_set,
    report_lines,
)

#: Per-primer reach for these 20 kb references.
EXTENSION = 500

#: A 12-mer that the fixtures assert is absent from every reference.
NOWHERE = "GATTACAGATTC"


# ----------------------------------------------------------------------
# Synthetic references
# ----------------------------------------------------------------------


def _sequence(length, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _oligo(seed, length=12):
    return _sequence(length, seed=10_000 + seed)


def _reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _plant(sequence, plants):
    seq = list(sequence)
    for primer, offsets in plants.items():
        for offset in offsets:
            seq[offset : offset + len(primer)] = list(primer)
    return "".join(seq)


def _write_fasta(path, sequence, name="ref"):
    path.write_text(
        f">{name}\n" + "\n".join(sequence[i : i + 70] for i in range(0, len(sequence), 70)) + "\n"
    )
    return str(path)


def _target(tmp_path, name, sequence):
    return ReferenceSpec(
        prefix=str(tmp_path / name),
        genome=_write_fasta(tmp_path / f"{name}.fasta", sequence, name),
        length=len(sequence),
        role=ROLE_TARGET,
    )


def _host(tmp_path, name, sequence):
    """A host given by FASTA with no table: sites by count, no coverage."""
    return ReferenceSpec(
        prefix=str(tmp_path / name),
        genome=_write_fasta(tmp_path / f"{name}.fasta", sequence, name),
        length=len(sequence),
        role=ROLE_HOST,
        scan=False,
    )


def _absent(oligo, *sequences):
    rc = _reverse_complement(oligo)
    return all(oligo not in s and rc not in s for s in sequences)


@pytest.fixture
def conditions():
    return ReactionConditions(temp=30.0, polymerase="phi29")


def _settings(**overrides):
    """Wide limits, so each test is about the one rule it names."""
    values = dict(extension=EXTENSION, max_dimer_bp=7, min_tm=-100.0, max_tm=200.0)
    values.update(overrides)
    return ImprovementSettings(**values)


def _pooled_coverage(assessment):
    """Length-weighted coverage over the targets: the figure a pooled view ranks on."""
    total = sum(r.length for r in assessment.targets)
    return sum(r.coverage.value * r.length for r in assessment.targets) / total


# ----------------------------------------------------------------------
# A hole on one target of two
# ----------------------------------------------------------------------


@pytest.fixture
def hole_on_b(tmp_path):
    """Target B is the worse covered. One candidate helps B; another helps A more.

    `in_set` covers the first half of A and the first quarter of B. `fills_b`
    binds only B, over 3 kb the set does not reach. `helps_a` binds only A, over
    the 10 kb the set does not reach, so it raises the POOLED coverage by more.
    """
    in_set, fills_b, helps_a = _oligo(1), _oligo(2), _oligo(3)
    a = _plant(
        _sequence(20_000, seed=101),
        {in_set: range(200, 10_000, 400), helps_a: range(10_200, 19_800, 400)},
    )
    b = _plant(
        _sequence(20_000, seed=202),
        {in_set: range(200, 5_000, 400), fills_b: range(5_200, 8_000, 400)},
    )
    host = _plant(_sequence(30_000, seed=303), {in_set: (1_000, 7_000), fills_b: (12_000,)})
    assert _absent(NOWHERE, a, b, host)
    assert _absent(fills_b, a) and _absent(helps_a, b)
    specs = [
        _target(tmp_path, "target_a", a),
        _target(tmp_path, "target_b", b),
        _host(tmp_path, "host", host),
    ]
    return {"specs": specs, "set": [in_set], "fills_b": fills_b, "helps_a": helps_a}


def test_the_add_fills_the_hole_on_the_worst_target(hole_on_b, conditions):
    """Not the candidate that raises the pooled coverage more."""
    case = hole_on_b
    pool = [case["helps_a"], case["fills_b"]]

    report = improve_set(case["set"], case["specs"], conditions, _settings(), pool)

    # The premise, measured rather than assumed: the other candidate is the
    # better one by the pooled figure, and B is the target that needs help.
    def fresh(extra):
        return evaluate_reference_panel(
            case["set"] + [extra], case["specs"], conditions, extension=EXTENSION
        )

    assert _pooled_coverage(fresh(case["helps_a"])) > _pooled_coverage(fresh(case["fills_b"]))
    assert "target_b" in report.current.worst_target_coverage.basis

    top = report.proposals[0]
    assert top.kind == KIND_ADD
    assert top.add == (case["fills_b"],)
    assert top.worst_target_gain > 0
    b_prefix = case["specs"][1].prefix
    assert top.coverage_change[b_prefix]["change"] > 0
    # The pooled favourite does nothing for the worst target, so it is not an
    # option at all rather than a lower-ranked one.
    assert all(case["helps_a"] not in p.add for p in report.proposals)


def test_applying_the_top_proposal_reproduces_the_prediction_exactly(hole_on_b, conditions):
    """Prediction and measurement are one code path, or this fails."""
    case = hole_on_b
    report = improve_set(
        case["set"], case["specs"], conditions, _settings(), [case["helps_a"], case["fills_b"]]
    )
    top = report.proposals[0]

    edited = [p for p in case["set"] if p not in top.drop] + list(top.add)
    assert tuple(edited) == top.resulting_set
    measured = evaluate_reference_panel(edited, case["specs"], conditions, extension=EXTENSION)

    assert measured.as_dict() == top.predicted.as_dict()
    # And the host is part of what was compared, by count, with no table.
    host = measured.hosts[0]
    assert host.sites.value == 3.0
    assert host.coverage.value is None


def test_the_current_figures_are_the_evaluation_itself(hole_on_b, conditions):
    """Step 1 calls the Phase 2 evaluation; it does not restate it."""
    case = hole_on_b
    report = improve_set(case["set"], case["specs"], conditions, _settings())
    measured = evaluate_reference_panel(case["set"], case["specs"], conditions, extension=EXTENSION)

    assert report.current.as_dict() == measured.as_dict()
    assert report.as_dict()["candidate_pool"]["unavailable"]
    assert not [p for p in report.proposals if p.kind != KIND_DROP]


def test_a_candidate_outside_the_tm_window_is_not_an_option(hole_on_b, conditions):
    case = hole_on_b
    tm = conditions.calculate_effective_tm(case["fills_b"])
    narrow = _settings(min_tm=tm + 5.0, max_tm=tm + 10.0)

    report = improve_set(case["set"], case["specs"], conditions, narrow, [case["fills_b"]])

    assert not report.proposals
    assert report.candidate_pool["rejected_tm"] == 1


def test_a_missed_panel_limit_is_reported_and_removes_nothing(hole_on_b, conditions):
    """The configured limits are advisory here: named, not enforced."""
    case = hole_on_b
    seen = []

    def limit_check(oligos):
        seen.append(tuple(oligos))
        return {"evaluated": True, "violations": ["worst hole above maximum"], "values": []}

    report = improve_set(
        case["set"],
        case["specs"],
        conditions,
        _settings(),
        [case["fills_b"]],
        limit_check=limit_check,
    )

    top = report.proposals[0]
    assert top.add == (case["fills_b"],)
    assert top.advisory["violations"] == ["worst hole above maximum"]
    assert top.resulting_set in seen
    assert any("misses a configured limit" in line for line in report_lines(report.as_dict()))


# ----------------------------------------------------------------------
# Drops
# ----------------------------------------------------------------------


@pytest.fixture
def one_binds_nothing(tmp_path):
    """Two oligos that cover a target, and one that binds no reference."""
    first, second = _oligo(11), _oligo(12)
    a = _plant(_sequence(20_000, seed=404), {first: (1_000, 9_000), second: (5_000, 15_000)})
    assert _absent(NOWHERE, a)
    return {"specs": [_target(tmp_path, "target_a", a)], "set": [first, second, NOWHERE]}


def test_an_oligo_that_binds_nothing_is_proposed_for_dropping(one_binds_nothing, conditions):
    case = one_binds_nothing

    report = improve_set(case["set"], case["specs"], conditions, _settings())

    drops = [p for p in report.proposals if p.kind == KIND_DROP]
    assert [p.drop for p in drops] == [(NOWHERE,)]
    assert "binds no target" in drops[0].reason
    # What is lost is stated per target, and here it is nothing.
    assert all(change["change"] == 0 for change in drops[0].coverage_change.values())

    row = next(r for r in report.attribution if r["oligo"] == NOWHERE)
    target = row["targets"][case["specs"][0].prefix]
    assert target["sites"] == 0 and target["status"] == "no_sites"
    assert target["marginal_coverage"]["value"] == 0
    assert target["sole_cover_of_some_region"] is False
    # The oligos that do cover something are the only cover of their regions.
    others = [r for r in report.attribution if r["oligo"] != NOWHERE]
    assert all(r["targets"][case["specs"][0].prefix]["sole_cover_of_some_region"] for r in others)


def test_a_fixed_oligo_is_never_proposed_for_dropping(one_binds_nothing, conditions):
    """Already ordered or validated: diagnosed like the rest, and left alone."""
    case = one_binds_nothing

    report = improve_set(
        case["set"], case["specs"], conditions, _settings(fixed_oligos=(NOWHERE.lower(),))
    )

    assert all(NOWHERE not in p.drop for p in report.proposals)
    row = next(r for r in report.attribution if r["oligo"] == NOWHERE)
    assert row["fixed"] is True
    assert row["targets"][case["specs"][0].prefix]["sites"] == 0


def test_sites_that_could_not_be_established_are_not_binding_nowhere(tmp_path, conditions):
    """One target cannot be read. Its zero is unknown, and nothing is dropped."""
    first = _oligo(21)
    a = _plant(_sequence(20_000, seed=505), {first: (1_000, 9_000)})
    assert _absent(NOWHERE, a)
    unreadable = ReferenceSpec(
        prefix=str(tmp_path / "target_b"),
        genome=str(tmp_path / "does_not_exist.fasta"),
        length=20_000,
        role=ROLE_TARGET,
    )
    specs = [_target(tmp_path, "target_a", a), unreadable]

    report = improve_set([first, NOWHERE], specs, conditions, _settings())

    assert not report.proposals
    assert "worst target is unknown" in report.proposals_unavailable
    row = next(r for r in report.attribution if r["oligo"] == NOWHERE)
    on_a, on_b = row["targets"][specs[0].prefix], row["targets"][unreadable.prefix]
    assert on_a["sites"] == 0 and on_a["status"] == "no_sites"
    assert on_b["sites"] is None and on_b["status"] == "unavailable"
    assert on_b["unavailable"]
    assert on_b["marginal_coverage"]["value"] is None
    assert on_b["sole_cover_of_some_region"] is None
    text = "\n".join(report_lines(report.as_dict()))
    assert "target_b: unavailable" in text


# ----------------------------------------------------------------------
# The dimer screen, on the proposals themselves
# ----------------------------------------------------------------------


def _dimer_case(tmp_path, seed):
    """A set of three, a hole on target B, and a pool planted across the hole.

    One candidate is built to dimerise with the first oligo of the set and is
    planted more densely than any other, so a ranking that skipped the screen
    would put it first.
    """
    rng = random.Random(seed)
    oligos = ["".join(rng.choice("ACGT") for _ in range(12)) for _ in range(3)]
    partner = _reverse_complement(oligos[0])[:8] + "".join(rng.choice("ACGT") for _ in range(4))
    pool = ["".join(rng.choice("ACGT") for _ in range(12)) for _ in range(30)]

    a = _plant(_sequence(20_000, seed=seed * 7 + 1), {o: range(200, 19_000, 900) for o in oligos})
    plants = {o: range(200, 6_000, 900) for o in oligos}
    for index, candidate in enumerate(pool):
        plants[candidate] = (7_000 + 400 * index,)
    plants[partner] = range(6_300, 19_500, 450)
    b = _plant(_sequence(20_000, seed=seed * 7 + 2), plants)
    specs = [_target(tmp_path, f"a{seed}", a), _target(tmp_path, f"b{seed}", b)]
    return oligos, [partner, *pool], specs, partner


@pytest.mark.parametrize("seed", [1, 2, 3, 4, 5])
def test_a_proposed_add_never_dimerises_with_an_oligo_that_stays(tmp_path, conditions, seed):
    oligos, pool, specs, partner = _dimer_case(tmp_path, seed)
    settings = _settings(max_dimer_bp=3, max_edits=50)

    report = improve_set(oligos, specs, conditions, settings, pool)

    adding = [p for p in report.proposals if p.add]
    assert adding, "the pool was built to hold at least one candidate that helps"
    for proposal in adding:
        retained = [o for o in proposal.resulting_set if o not in proposal.add]
        for added in proposal.add:
            for kept in retained:
                assert not is_dimer_fast(added, kept, max_dimer_bp=3), (added, kept)
                assert not is_dimer_fast(kept, added, max_dimer_bp=3), (added, kept)

    # The screen did something: the candidate built to conflict is in no
    # proposal beside the oligo it conflicts with, and some were rejected.
    assert is_dimer_fast(partner, oligos[0], max_dimer_bp=3)
    assert all(not (partner in p.add and oligos[0] in p.resulting_set) for p in report.proposals)
    pool_record = report.candidate_pool
    assert pool_record["rejected_dimer"] + pool_record["eligible_swap"] > 0

    # A swap drops exactly the oligo its add could not sit beside.
    for proposal in report.proposals:
        if proposal.kind == KIND_SWAP:
            assert len(proposal.drop) == 1 and proposal.drop[0] in proposal.depends_on
            assert is_dimer_fast(proposal.add[0], proposal.drop[0], max_dimer_bp=3)

    # Every reported prediction, of every kind, is what a fresh evaluation gives.
    for proposal in report.proposals:
        measured = evaluate_reference_panel(
            list(proposal.resulting_set), specs, conditions, extension=EXTENSION
        )
        assert measured.as_dict() == proposal.predicted.as_dict()


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_a_fixed_oligo_is_not_swapped_out_either(tmp_path, conditions, seed):
    oligos, pool, specs, _partner = _dimer_case(tmp_path, seed)
    settings = _settings(max_dimer_bp=3, max_edits=50, fixed_oligos=tuple(oligos))

    report = improve_set(oligos, specs, conditions, settings, pool)

    assert all(not p.drop for p in report.proposals)
    assert report.candidate_pool["eligible_swap"] == 0


def test_the_ranking_is_worst_target_gain_then_host_load_then_size(tmp_path, conditions):
    """Two adds with the same gain: the one lighter on the host comes first."""
    in_set, light, heavy = _oligo(41), _oligo(42), _oligo(43)
    a = _plant(_sequence(20_000, seed=808), {in_set: range(200, 19_000, 400)})
    b = _plant(
        _sequence(20_000, seed=909),
        {in_set: range(200, 5_000, 400), light: (9_000,), heavy: (14_000,)},
    )
    host = _plant(_sequence(30_000, seed=111), {heavy: range(1_000, 20_000, 1_000)})
    specs = [
        _target(tmp_path, "target_a", a),
        _target(tmp_path, "target_b", b),
        _host(tmp_path, "host", host),
    ]

    report = improve_set([in_set], specs, conditions, _settings(), [heavy, light])

    adds = [p for p in report.proposals if p.kind == KIND_ADD]
    assert [p.add for p in adds] == [(light,), (heavy,)]
    assert adds[0].worst_target_gain == adds[1].worst_target_gain
    assert adds[0].host_load.value < adds[1].host_load.value

    capped = improve_set(
        [in_set], specs, conditions, _settings(max_edits=1), [heavy, light]
    ).as_dict()
    section = capped["sections"][KIND_ADD]
    assert [p["add"] for p in section["proposals"]] == [[light]]
    assert section["considered"] == 2
    assert section["shown"] == 1


# ----------------------------------------------------------------------
# Sections: what is an improvement, what is a trade-off, and what is shown
# ----------------------------------------------------------------------


def test_an_ordinary_set_is_not_told_to_drop_a_member_for_host_load(tmp_path, conditions):
    """Six oligos alike on the target, with 5, 5, 5, 5, 5 and 6 host sites.

    Removing the sixth raises the target-against-host density, as removing the
    member with the lowest ratio always does. That is a trade-off with a
    coverage cost, and it must be reported as one: in its own section, with the
    two densities and the cost, and not among the drops that cost nothing.
    """
    oligos = [_oligo(60 + index) for index in range(6)]
    plants = {oligo: (500 + 3_000 * index,) for index, oligo in enumerate(oligos)}
    target = _plant(_sequence(20_000, seed=1_001), plants)
    host_plants = {
        oligo: [300 + 700 * index + 4_000 * copy for copy in range(6 if index == 5 else 5)]
        for index, oligo in enumerate(oligos)
    }
    host = _plant(_sequence(30_000, seed=1_002), host_plants)
    specs = [_target(tmp_path, "target", target), _host(tmp_path, "host", host)]

    report = improve_set(oligos, specs, conditions, _settings())
    payload = report.as_dict()

    host_sites = [row["hosts"][specs[1].prefix]["sites"] for row in report.attribution]
    assert host_sites == [5, 5, 5, 5, 5, 6], "the fixture is the near-uniform set it claims"
    target_sites = [row["targets"][specs[0].prefix]["sites"] for row in report.attribution]
    assert target_sites == [1] * 6

    assert report.section(KIND_DROP).considered == 0
    trade = payload["sections"][KIND_TRADE_OFF]
    assert trade["considered"] == 1
    assert [p["drop"] for p in trade["proposals"]] == [[oligos[5]]]
    # The facts a reader needs, carried as data.
    assert trade["is_improvement"] is False
    assert trade["lowest_ratio_member_always_qualifies"] is True
    assert trade["each_entry_evaluated_alone"] is True
    assert trade["jointly_applicable"] is False
    assert payload["sections"][KIND_DROP]["jointly_applicable"] is False

    # What is measured: two densities, and the coverage the removal costs.
    entry = trade["proposals"][0]
    density = entry["worst_host_selectivity_density"]
    assert density["after"]["value"] > density["before"]["value"]
    assert entry["per_target_coverage"][specs[0].prefix]["change"] < 0
    assert entry["worst_target_coverage_gain"] < 0
    assert entry["kind"] == KIND_TRADE_OFF

    text = "\n".join(report_lines(payload))
    assert "host load high relative" not in text
    assert "host load high relative" not in str(payload)
    assert trade["title"] in text and trade["note"] in text


def test_a_drop_that_costs_coverage_is_a_trade_off_with_its_cost_stated(tmp_path, conditions):
    """An oligo with most of the host sites and little of the target."""
    good, costly = _oligo(31), _oligo(32)
    a = _plant(
        _sequence(20_000, seed=606),
        {good: range(200, 16_000, 400), costly: (18_000,)},
    )
    host = _plant(_sequence(30_000, seed=707), {good: (1_000,), costly: range(2_000, 28_000, 500)})
    specs = [_target(tmp_path, "target_a", a), _host(tmp_path, "host", host)]

    report = improve_set([good, costly], specs, conditions, _settings())

    assert not report.section(KIND_DROP).proposals
    (drop,) = report.section(KIND_TRADE_OFF).proposals
    assert drop.drop == (costly,)
    assert drop.worst_target_gain < 0
    assert drop.coverage_change[specs[0].prefix]["change"] < 0
    assert drop.density_after.value > drop.density_before.value
    assert drop.host_load.value < report.current_host_load.value


def test_helpful_candidates_do_not_push_the_useless_oligo_off_the_report(tmp_path, conditions):
    """One covering oligo, one that binds nothing, and a pool of eight that help.

    `max_edits` bounds each section on its own. Applied to one ranked list, the
    five best adds filled the report and the oligo that binds nothing was never
    shown, which is the first thing this command exists to say.
    """
    covering = _oligo(71)
    pool = [_oligo(80 + index) for index in range(8)]
    plants = {covering: range(200, 6_000, 400)}
    for index, candidate in enumerate(pool):
        plants[candidate] = (7_000 + 1_500 * index,)
    target = _plant(_sequence(20_000, seed=1_101), plants)
    assert _absent(NOWHERE, target)
    specs = [_target(tmp_path, "target", target)]

    report = improve_set([covering, NOWHERE], specs, conditions, _settings(), pool)

    assert report.settings.max_edits == 5, "the default, not a value chosen for the test"
    assert [p.drop for p in report.section(KIND_DROP).proposals] == [(NOWHERE,)]
    adds = report.section(KIND_ADD)
    assert (len(adds.proposals), adds.considered) == (5, 8)
    text = "\n".join(report_lines(report.as_dict()))
    assert "shown 5 of 8 considered" in text
    assert "shown 1 of 1 considered" in text

    # Asking for none shows none and still says how many there were.
    silent = improve_set(
        [covering, NOWHERE], specs, conditions, _settings(max_edits=0), pool
    ).as_dict()
    assert all(section["shown"] == 0 for section in silent["sections"].values())
    assert silent["sections"][KIND_ADD]["considered"] == 8
    text = "\n".join(report_lines(silent))
    assert "shown 0 of 8 considered" in text
    assert "no single edit was found" not in text
    assert text.count("none found") == 2, "only the two sections with nothing considered"


# ----------------------------------------------------------------------
# Unknown is not success: candidates, and oligos an index does not hold
# ----------------------------------------------------------------------


def _write_index(prefix, entries, k=12):
    """A position index holding exactly `entries`, written through its one door."""
    from neoswga.core.position_index import index_path, write_entries

    full = {}
    for oligo, positions in entries.items():
        full[oligo] = list(positions)
        full.setdefault(_reverse_complement(oligo), [])
    write_entries(index_path(str(prefix), k), full, record_starts=[0])


def test_a_candidate_with_unknown_host_binding_is_not_an_option(tmp_path, conditions):
    """The host index holds the set's oligo and not the candidate, and is not scanned.

    The candidate fills the hole on the target. Its host binding is unknown,
    which is not a low host load: it is not proposed, it is counted as
    unmeasured, and the report says why.
    """
    in_set, candidate = _oligo(91), _oligo(92)
    target = _plant(
        _sequence(20_000, seed=1_201),
        {in_set: range(200, 6_000, 400), candidate: range(8_000, 12_000, 400)},
    )
    host_prefix = tmp_path / "indexed_host"
    _write_index(host_prefix, {in_set: [1_000, 9_000]})
    host = ReferenceSpec(
        prefix=str(host_prefix), genome=None, length=30_000, role=ROLE_HOST, scan=False
    )
    specs = [_target(tmp_path, "target", target), host]

    report = improve_set([in_set], specs, conditions, _settings(), [candidate])

    assert report.current.hosts[0].sites.value == 2.0, "the host is measured for the set itself"
    assert not report.proposals
    pool = report.candidate_pool
    assert pool["eligible_add"] == 1 and pool["unmeasured"] == 1
    (example,) = pool["unmeasured_examples"]
    assert example["candidate"] == candidate
    assert example["references"] == [host.prefix]
    assert "absent from the position index" in example["reason"]
    assert candidate in "\n".join(report_lines(report.as_dict()))


def test_an_oligo_the_index_does_not_hold_is_unavailable_not_absent(tmp_path, conditions):
    """A target index holding one of two oligos, with nothing to scan.

    The indexed oligo keeps its count. The other is unavailable with the
    reason, and no edit is proposed: in particular it is not dropped as an
    oligo that binds nothing.
    """
    indexed, outside = _oligo(95), _oligo(96)
    prefix = tmp_path / "indexed_target"
    _write_index(prefix, {indexed: [500, 6_000, 12_000]})
    target = ReferenceSpec(
        prefix=str(prefix), genome=None, length=20_000, role=ROLE_TARGET, scan=False
    )

    report = improve_set([indexed, outside], [target], conditions, _settings())

    assert not report.proposals
    assert "worst target is unknown" in report.proposals_unavailable
    rows = {row["oligo"]: row["targets"][target.prefix] for row in report.attribution}
    assert rows[indexed]["sites"] == 3 and rows[indexed]["status"] == "ok"
    assert rows[outside]["sites"] is None and rows[outside]["status"] == "unavailable"
    assert "absent from the position index" in rows[outside]["unavailable"]
    assert rows[outside]["sole_cover_of_some_region"] is None


# ----------------------------------------------------------------------
# The shared reading must not answer for the wrong reference or primer
# ----------------------------------------------------------------------


def test_two_counted_hosts_keep_their_own_counts(tmp_path, conditions):
    """One oligo, two hosts counted from FASTA, different site counts.

    A memo keyed on the k-mer alone would hand the first host's count to the
    second. Each host must report its own, and the whole evaluation must equal
    one made with nothing shared.
    """
    in_set, other = _oligo(101), _oligo(102)
    target = _plant(
        _sequence(20_000, seed=1_301),
        {in_set: range(200, 9_000, 400), other: range(10_000, 15_000, 400)},
    )
    first = _plant(_sequence(30_000, seed=1_302), {in_set: (1_000, 9_000)})
    second = _plant(_sequence(30_000, seed=1_303), {in_set: range(1_000, 26_000, 5_000)})
    specs = [
        _target(tmp_path, "target", target),
        _host(tmp_path, "first_host", first),
        _host(tmp_path, "second_host", second),
    ]

    report = improve_set([in_set, other], specs, conditions, _settings())

    assert [host.sites.value for host in report.current.hosts] == [2.0, 5.0]
    fresh = evaluate_reference_panel([in_set, other], specs, conditions, extension=EXTENSION)
    assert report.current.as_dict() == fresh.as_dict()
    by_host = {row["oligo"]: row["hosts"] for row in report.attribution}
    assert by_host[in_set][specs[1].prefix]["sites"] == 2
    assert by_host[in_set][specs[2].prefix]["sites"] == 5


def test_a_primer_arriving_late_is_looked_up_not_answered_from_the_held_cache(tmp_path, conditions):
    """The shared position cache is rebuilt when asked about a primer it lacks.

    Returning the held cache would answer "no sites" for the newcomer, which
    is the silent zero the evaluation exists to avoid.
    """
    from neoswga.core.set_improvement import MemoisedSources

    early, late = _oligo(111), _oligo(112)
    target = _plant(
        _sequence(20_000, seed=1_401),
        {early: range(200, 6_000, 400), late: range(9_000, 15_000, 400)},
    )
    host = _plant(_sequence(30_000, seed=1_402), {late: (1_000, 2_000, 3_000)})
    specs = [_target(tmp_path, "target", target), _host(tmp_path, "host", host)]
    sources = MemoisedSources()

    evaluate_reference_panel([early], specs, conditions, extension=EXTENSION, sources=sources)
    shared = evaluate_reference_panel(
        [early, late], specs, conditions, extension=EXTENSION, sources=sources
    )
    fresh = evaluate_reference_panel([early, late], specs, conditions, extension=EXTENSION)

    assert shared.targets[0].per_primer_sites[late] == 15
    assert shared.hosts[0].per_primer_sites[late] == 3
    assert shared.as_dict() == fresh.as_dict()


def test_a_stability_floor_with_no_temperature_is_refused(tmp_path, conditions):
    """`max_dimer_dg` is a free energy at a temperature; none is assumed."""
    first, second = _oligo(121), _oligo(122)
    target = _plant(_sequence(20_000, seed=1_501), {first: (1_000,), second: (9_000,)})
    specs = [_target(tmp_path, "target", target)]
    floor = _settings(max_dimer_dg=-6.0)

    with pytest.raises(ValueError, match="no temperature"):
        improve_set([first, second], specs, None, floor)

    # With a reaction the floor is evaluated, and the report is produced.
    report = improve_set([first, second], specs, conditions, floor)
    assert report.settings.max_dimer_dg == -6.0
