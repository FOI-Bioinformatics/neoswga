"""The step-2 background frequency gate, under either reading of "the background".

`core/background_aggregate.py` defines what `background_aggregate` means;
this module is what applies it to the gate `filter.get_all_rates` runs. The
two live apart because the vocabulary is shared with Phase 7c's step-4 figures
while the arithmetic here is step 2's alone, and because `filter.py` and
`pipeline.py` are both near their size budgets.

What the two readings do, stated as the gate sees them:

- `sum`: one count over one length. The count is the sum over every background
  prefix and the length is the sum of every background genome's length, and the
  candidate passes when that pooled frequency is below the scaled
  `max_bg_freq`. This is what the pipeline has always computed and
  `pooled_verdict` is the rule transcribed unchanged.
- `worst-case`: one count per host over that host's own length. The candidate
  passes only when EVERY host's own frequency is below the same scaled ceiling.
  Nothing is pooled, so a host that is a small share of the pooled length is no
  longer diluted by the large ones.

Why the difference is not a matter of taste. Against the *Wolbachia* reference
with three hosts pooled (3.3 Mb, 144 Mb and 3.30 Gb), hg38 is 95.7 percent of
the pooled length, so the pooled ceiling sits within 5 percent of hg38's own
ceiling while the numerator still charges a candidate for its sites on the two
small hosts. A candidate frequent on a small host and clean on the large one
passes the pooled gate on the strength of the large host's length
(`docs/validation/2026-10-02-diversity-baseline.md`, Q4 and Q5). This is
Known Issue 6 -- a figure that moves with background size -- inside the gate.

**An absent per-host count blocks under `worst-case`, and that is a deliberate
difference from today.** Under `sum` an absent POOLED count passes
(`filter.get_all_rates`, the comment above the loop: a k-mer missing from a
counter's table may still have sites the string search finds, and the later
stages filter again). That argument does not carry to a worst case: the worst
of a set one of whose members was never measured is not a measured worst, and
reporting it as a pass would report a specificity the panel has not been shown
to have. So under `worst-case` an absent or unpairable per-host count is a
rejection, and under `sum` the old rule stands untouched. `candidate_retention`
and the Gini stage still see whatever survives, exactly as before.

Three configurations refuse rather than guess, all under `worst-case` only:

- prefixes and lengths that cannot be paired one-to-one, because a count would
  then be divided by another host's length;
- a host whose length is missing or not positive, because there is no ceiling
  to compare against;
- a Bloom filter background, because the filter answers for the pooled
  background and cannot be split per host (`filter.get_bg_rates_via_bloom`).

Per-host columns are written to `step2_df.csv` only when there are at least two
hosts. With one host the column would repeat `bg_count` value for value, and
withholding it keeps every single-host design's table byte-identical, which is
what makes `worst-case` with one host the same design as `sum` with one host
rather than merely the same verdicts.
"""

import os
from collections.abc import Mapping
from dataclasses import dataclass

from .background_aggregate import SUM, WORST_CASE, resolve_background_aggregate
from .exceptions import InvalidDesignRequest

__all__ = [
    "BackgroundHost",
    "COLUMN_PREFIX",
    "MIN_HOSTS_FOR_COLUMNS",
    "columns_are_written",
    "host_column_names",
    "host_panel",
    "per_host_columns",
    "pool_counts",
    "pooled_verdict",
    "refuse_a_pooled_only_background",
    "worst_case_verdict",
]

#: Per-host count columns in `step2_df.csv` are named `bg_count_<basename>`.
COLUMN_PREFIX = "bg_count_"

#: Below this many hosts the per-host column carries nothing `bg_count` does
#: not already carry, so it is not written.
MIN_HOSTS_FOR_COLUMNS = 2


@dataclass(frozen=True)
class BackgroundHost:
    """One host of the background panel, paired with its OWN length.

    The pairing is made once, where the prefix list and the length list are
    both in hand, and travels together afterwards. A prefix carries no record
    of the genome it was counted from, so looking the length up later would be
    a guess, and the wrong guess divides one host's count by another's length.
    """

    prefix: str
    length: int
    column: str


def host_column_names(prefixes: list[str]) -> dict[str, str]:
    """A deterministic `step2_df.csv` column name per background prefix.

    The basename, because a prefix is a path and the directory is the same for
    every host of a run. Two prefixes with the same basename in different
    directories would collide, so the second and later get their position
    appended; the order is the configured order, never a set's.
    """
    names: dict[str, str] = {}
    taken: set[str] = set()
    for index, prefix in enumerate(prefixes):
        base = os.path.basename(str(prefix)) or f"host{index}"
        column = f"{COLUMN_PREFIX}{base}"
        if column in taken:
            column = f"{COLUMN_PREFIX}{base}__{index}"
        taken.add(column)
        names[str(prefix)] = column
    return names


def host_panel(
    bg_prefixes: list[str] | None,
    bg_seq_lengths: list[int] | None,
    aggregate: str | None = None,
) -> tuple[BackgroundHost, ...]:
    """The background hosts, each paired with its own length.

    Empty when there is no background at all, which is a configuration and not
    an unmeasured host: a design that lists no host has no background gate to
    apply under either reading, and `sum` passes every candidate there today.

    Under `sum` an unpairable pair of lists gives an empty panel, because the
    pooled gate needs only the summed length and the panel is then used for
    the report columns alone. Under `worst-case` the same input refuses: the
    gate cannot divide a host's count by its own length if it cannot tell which
    length is its own.
    """
    resolved = resolve_background_aggregate(aggregate)
    prefixes = [str(prefix) for prefix in (bg_prefixes or [])]
    lengths = list(bg_seq_lengths or [])

    if not prefixes:
        return ()

    unusable = len(prefixes) != len(lengths) or any(
        length is None or int(length) <= 0 for length in lengths
    )
    if unusable:
        if resolved == WORST_CASE:
            raise InvalidDesignRequest(
                "background_aggregate",
                (
                    f"'{WORST_CASE}' judges each host on its own length, so every "
                    f"background prefix needs one positive length of its own. Got "
                    f"{len(prefixes)} prefixes and {len(lengths)} lengths "
                    f"({lengths!r}). Set bg_seq_lengths with one entry per entry of "
                    f"bg_genomes, in the same order, or use '{SUM}', which needs "
                    f"only the summed length"
                ),
                resolved,
            )
        return ()

    columns = host_column_names(prefixes)
    return tuple(
        BackgroundHost(prefix=prefix, length=int(length), column=columns[prefix])
        for prefix, length in zip(prefixes, lengths, strict=True)
    )


def refuse_a_pooled_only_background(aggregate: str, source: str, detail: str) -> None:
    """Refuse `worst-case` on a background that exists only as one pooled answer.

    The Bloom filter path is the case: it records membership for the whole
    background and `parameter` deliberately leaves `bg_prefixes` empty for it,
    so there is no per-host count to take a worst case over. Passing every
    candidate there would be the "unknown reported as success" shape, and
    silently pooling would report a figure the user did not ask for.
    """
    if aggregate != WORST_CASE:
        return
    raise InvalidDesignRequest(
        "background_aggregate",
        (
            f"'{WORST_CASE}' needs one count per host, and {source} answers for the "
            f"pooled background only ({detail}). Use '{SUM}', which is what a "
            f"pooled background supports, or count the hosts into per-host k-mer "
            f"tables so each can be judged on its own length"
        ),
        aggregate,
    )


def pool_counts(by_prefix: Mapping[str, Mapping[str, int]]) -> dict[str, int]:
    """One count per primer, summed over the prefixes: today's pooled numerator.

    Separate from the per-prefix dictionaries rather than instead of them. The
    per-(prefix, k) results already arrive apart in
    `filter.get_rates_for_one_species`; keeping them is free and summing them
    here reproduces what that function returned before, key order included,
    because the prefixes are visited in the configured order.
    """
    pooled: dict[str, int] = {}
    for counts in by_prefix.values():
        for primer, count in counts.items():
            if primer in pooled:
                pooled[primer] += count
            else:
                pooled[primer] = count
    return pooled


def pooled_verdict(count: int | None, pooled_length: int, scaled_max_bg: float) -> bool:
    """Today's background gate, transcribed unchanged.

    An absent count passes: a k-mer missing from a counter's table may still
    have sites the later string search finds, and the evenness and ranking
    stages filter again. A count with no length to divide by is uninterpretable
    and fails, which is the guard that turned a zero-length background FASTA
    from a per-primer ZeroDivisionError into a clear rejection.
    """
    if count is None or pooled_length <= 0:
        return count is None
    return count / pooled_length < scaled_max_bg


def worst_case_verdict(
    by_prefix: Mapping[str, Mapping[str, int]],
    primer: str,
    hosts: tuple[BackgroundHost, ...],
    scaled_max_bg: float,
) -> bool:
    """Pass only when every host's OWN frequency is below the ceiling.

    With no hosts this passes, for the reason `host_panel` gives: no host
    listed is a configuration, not an unmeasured host.

    The reduction runs over every host, not over the hosts that answered. A
    host with no entry, or an entry of None, is unknown, and the worst case
    over a set with an unmeasured member is unmeasured, so it rejects. That is
    the one place this gate is stricter than the pooled one.

    "This host was never asked" and "this host answered None" are tested
    apart, because they are two different ways of not knowing and one
    `.get(primer)` would collapse them. Both reject, and neither is read as a
    count.
    """
    for host in hosts:
        if host.prefix not in by_prefix:
            return False
        counts = by_prefix[host.prefix]
        if primer not in counts:
            return False
        count = counts[primer]
        if count is None:
            return False
        if host.length <= 0:
            return False
        if not count / host.length < scaled_max_bg:
            return False
    return True


def per_host_columns(
    by_prefix: Mapping[str, Mapping[str, int]],
    primer: str,
    hosts: tuple[BackgroundHost, ...],
) -> list[int | None]:
    """This primer's count on each host, in panel order, `None` where unknown.

    `None` and not 0: a host that was never asked and a host that answered zero
    are different measurements, and the column is read by a person deciding
    whether the gate had an answer to apply.
    """
    values: list[int | None] = []
    for host in hosts:
        counts = by_prefix.get(host.prefix)
        count = None if counts is None else counts.get(primer)
        values.append(None if count is None else int(count))
    return values


def columns_are_written(hosts: tuple[BackgroundHost, ...]) -> bool:
    """Whether `step2_df.csv` gains the per-host columns for this panel."""
    return len(hosts) >= MIN_HOSTS_FOR_COLUMNS
