"""Does a position-dependent mismatch model rank the measured outcomes better?

Phase 4b, validation item 3. The model in `core/mismatch_model.py` is only
worth having if it predicts better than the uniform one, so this scores the
published primer sets in `tests/validation/data/` both ways and compares how
well each ranks the wet-lab outcome the paper reported.

Methodology is the 2026-09-16 record's
(`docs/validation/heldout_panel_discrimination_2026-09-16.md`) and the controls
are the same ones:

- area under the ROC curve within a study, rank-based, ties counted as half;
- leave-one-study-out thresholded accuracy against the majority-class baseline;
- every figure reported with and without the deliberate negative controls;
- the paper's own published fg/bg ratio carried along as a reference
  predictor, because it was the best one in that record. A new quantity that
  does not beat it has not beaten anything.

**What is different from that record, and it matters.** Those predictors were
the papers' own published numbers, which made it a check on the metric
definitions. These are recomputed here from k-mer tables, so this is a check on
THIS implementation against those outcomes, and it needs the genomes. A study
whose genome or k-mer table is not on this machine is reported as not measured,
with the reason; it is never dropped silently and never scored against a
substitute without the substitution being named in the output.

Usage:

    python scripts/benchmarking/mismatch_model_ranking.py \\
        --table-dir /tmp/neoswga-mismatch-ranking-tables --build-tables
    python scripts/benchmarking/mismatch_model_ranking.py \\
        --table-dir /tmp/neoswga-mismatch-ranking-tables --study prevotella

`--table-dir` is where tables counted by this script are written and looked for
first; the repository's own genome directory is searched after it, so a table
that already exists is never recounted. It defaults to a directory under the
system temp dir and REFUSES a path inside the repository: counting writes `.jf`
tables and a `genome_lengths.json`, and a run has to leave
`git status --porcelain` unchanged. Nothing is ever written beside a reference
genome either.

Saturation check (validation item 4) is a separate mode:

    python scripts/benchmarking/mismatch_model_ranking.py --saturation \\
        --pool step3_df.csv --host-prefix PREFIX --k 12
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
import pathlib
import statistics
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))


def _data_root():
    """Where the reference genomes live.

    They are gitignored, so a git WORKTREE has none of them while the main
    checkout does. Running this from a worktree and silently finding no genome
    would read as "this machine has no data", which is a different finding from
    "this checkout has no data". `--genome-root` overrides.
    """
    marker = ROOT / ".git"
    if marker.is_file():
        text = marker.read_text().strip()
        if text.startswith("gitdir:"):
            gitdir = pathlib.Path(text.split(":", 1)[1].strip())
            parts = gitdir.parts
            if ".git" in parts:
                main = pathlib.Path(*parts[: parts.index(".git")])
                if (main / "tests" / "validation" / "genomes").is_dir():
                    return main
    return ROOT


DATA = ROOT / "tests" / "validation" / "data"
DATA_ROOT = _data_root()
GENOMES = DATA_ROOT / "tests" / "validation" / "genomes"
EXAMPLE_INPUT = DATA_ROOT / "examples" / "wolbachia_pool_design" / "input"

#: Reaction the loads are computed under. The papers do not state conditions
#: fully enough to reconstruct each one, so ONE reaction is used for every
#: panel and named here rather than silently varied. All three studies used a
#: phi29-class isothermal amplification.
REACTION = {"temp": 30.0, "polymerase": "phi29"}

#: Highest mismatch distance summed over. The shipped default. Distance 2 costs
#: about 19x the time (measured; see the validation record) and is reachable
#: through `max_mismatches`, not through a flag here.
MAX_MISMATCHES = 1


class Study:
    """One dataset, its genomes, and what stands in for what.

    `background_note` is not decoration: a background that is not the one the
    paper used changes what the numbers mean, and Known Issue 6 records that a
    partial background moves `selectivity_ratio`. It is printed with every
    figure the study contributes.
    """

    def __init__(
        self,
        key,
        label_key,
        foreground,
        background,
        *,
        background_note="",
        foreground_note="",
    ):
        self.key = key
        self.label_key = label_key
        self.foreground = foreground
        self.background = background
        self.background_note = background_note
        self.foreground_note = foreground_note


STUDIES = (
    Study(
        "clarke_2017_mtb",
        "most_effective",
        ("mtb", GENOMES / "mtb.fna"),
        ("human_chr21", GENOMES / "human_chr21.fna"),
        background_note=(
            "the paper designed against the whole human genome; human chr21 "
            "(46.7 Mb of 3.30 Gb) stands in because no hg38 k-mer table "
            "exists on this machine and counting one was out of scope. "
            "Known Issue 6: a partial background is not a specific design"
        ),
    ),
    Study(
        "clarke_2017_wolbachia",
        "effective",
        ("wolbachia", GENOMES / "wolbachia.fna"),
        ("drosophila", EXAMPLE_INPUT / "drosophila.fna"),
        foreground_note="wMel, the reference strain in this repository",
    ),
    Study(
        "dwivedi_yu_2023_prevotella",
        "successful",
        ("prevotella", GENOMES / "prevotella.fna"),
        ("human_chr21", GENOMES / "human_chr21.fna"),
        background_note=(
            "the paper designed against the whole human genome; human chr21 " "stands in, as above"
        ),
    ),
)

MODELS = ("uniform", "position-dependent", "position-dependent-3prime")


# ----------------------------------------------------------------------
# Tables: found, or counted into a directory of our own
# ----------------------------------------------------------------------


def _search_dirs(table_dir):
    return [pathlib.Path(table_dir), GENOMES, EXAMPLE_INPUT]


def find_prefix(name, k, table_dir):
    """A prefix with a counted table for `k`, or None.

    Goes through `kmer_tables.table_exists`, never a filename glob, so a KMC
    database answers as well as a text dump.
    """
    from neoswga.core import kmer_tables

    for directory in _search_dirs(table_dir):
        prefix = str(directory / name)
        if kmer_tables.table_exists(prefix, k):
            return prefix
    return None


def build_table(name, genome, k, table_dir):
    """Count `genome` at `k` into `table_dir`, through the project's counter.

    Writes only into `table_dir`. The reference genome directory is read-only
    as far as this script is concerned.
    """
    from neoswga.core.kmer_backend import select_backend

    out = pathlib.Path(table_dir)
    out.mkdir(parents=True, exist_ok=True)
    backend = select_backend("jellyfish" if not _have("kmc") else None)
    started = time.perf_counter()
    backend.count(str(genome), k, str(out / name))
    return time.perf_counter() - started


def _have(tool):
    import shutil

    return shutil.which(tool) is not None


def genome_length(path, cache):
    """Total residues, measured by streaming. Cached in `cache` between runs.

    A length is needed for `selectivity_density`, which is the figure Known
    Issue 6 says to compare instead of the ratio.
    """
    key = str(path)
    if key in cache:
        return cache[key]

    from neoswga.core.genome_io import GenomeLoader

    total = sum(len(chunk) for chunk in GenomeLoader().load_genome_streaming(str(path)))
    cache[key] = total
    return total


# ----------------------------------------------------------------------
# Loads
# ----------------------------------------------------------------------


def panel_loads(primers, fg_prefix, bg_prefix, model):
    """`(fg_load, bg_load)` for one panel under one model, or None.

    None when a table for some oligo length is missing: a panel of mixed
    lengths needs a table per length, and a load summed over the lengths that
    happened to have one is a load over an unknown subset of the panel. That is
    the "ran silently over the measured members" defect, so it is refused here.
    """
    from neoswga.core.occupancy import weighted_site_load
    from neoswga.core.reaction_conditions import ReactionConditions

    conditions = ReactionConditions(**REACTION)
    fg = bg = 0.0
    for k in sorted({len(p) for p in primers}):
        same_k = [p for p in primers if len(p) == k]
        try:
            fg += weighted_site_load(
                same_k, [fg_prefix(k)], conditions, MAX_MISMATCHES, None, model
            )
            bg += weighted_site_load(
                same_k, [bg_prefix(k)], conditions, MAX_MISMATCHES, None, model
            )
        except (FileNotFoundError, OSError, TypeError):
            return None
    return fg, bg


def auc(scores, labels, higher_is_better):
    """Rank-based, threshold-free, ties counted as half.

    Written out rather than imported so a reader can see that a tie contributes
    0.5, which matters on samples this small. Identical to the 2026-09-16
    script's function, deliberately.
    """
    pos = [s for s, y in zip(scores, labels, strict=True) if y]
    neg = [s for s, y in zip(scores, labels, strict=True) if not y]
    if not pos or not neg:
        return None
    wins = sum(
        1.0 if (a > b) == higher_is_better and a != b else 0.5 if a == b else 0.0
        for a, b in itertools.product(pos, neg)
    )
    return wins / (len(pos) * len(neg))


def best_threshold(rows, extract, higher_is_better):
    """The cut maximising accuracy on these rows. Ties go to the lower cut."""
    values = sorted({extract(r) for r in rows})
    best, best_acc = None, -1.0
    for cut in values:
        correct = sum(((extract(r) >= cut) == higher_is_better) == r["label"] for r in rows)
        if correct / len(rows) > best_acc:
            best, best_acc = cut, correct / len(rows)
    return best, best_acc


# ----------------------------------------------------------------------
# Scoring
# ----------------------------------------------------------------------


def score_study(study, table_dir, build, length_cache):
    """Rows for one study, or a reason it could not be scored."""
    from neoswga.core.selectivity import (
        selectivity_density_from_loads,
        selectivity_from_loads,
    )

    payload = json.loads((DATA / f"{study.key}.json").read_text())
    panels = {
        set_id: entry
        for set_id, entry in payload["sets"].items()
        if entry.get("primers") and entry.get("published")
    }
    if not panels:
        return None, "the dataset records no panel with both primers and published metrics"

    lengths = sorted({len(p) for entry in panels.values() for p in entry["primers"]})

    prefixes = {}
    for name, genome in (study.foreground, study.background):
        if not pathlib.Path(genome).exists():
            return None, f"{name}: {genome} is not on this machine"
        for k in lengths:
            prefix = find_prefix(name, k, table_dir)
            if prefix is None and build:
                took = build_table(name, genome, k, table_dir)
                print(f"    counted {name} at k={k} in {took:.1f} s", flush=True)
                prefix = find_prefix(name, k, table_dir)
            if prefix is None:
                return None, (
                    f"no k-mer table for {name} at k={k} "
                    f"(the panels use oligo lengths {lengths}); "
                    f"re-run with --build-tables to count one"
                )
            prefixes[(name, k)] = prefix

    fg_name = study.foreground[0]
    bg_name = study.background[0]
    fg_length = genome_length(study.foreground[1], length_cache)
    bg_length = genome_length(study.background[1], length_cache)

    rows = []
    for set_id, entry in sorted(panels.items()):
        primers = list(entry["primers"])
        row = {
            "id": set_id,
            "label": bool(entry["outcome"].get(study.label_key)),
            "control": bool(entry["outcome"].get("negative_control")),
            "published_fg_bg_ratio": entry["published"]["fg_bg_ratio"],
            "size": len(primers),
        }
        for model in MODELS:
            loads = panel_loads(
                primers,
                lambda k, n=fg_name: prefixes[(n, k)],
                lambda k, n=bg_name: prefixes[(n, k)],
                model,
            )
            if loads is None:
                return None, f"panel {set_id}: a load could not be computed under {model}"
            fg, bg = loads
            row[f"{model}:fg_load"] = fg
            row[f"{model}:bg_load"] = bg
            row[f"{model}:selectivity"] = selectivity_from_loads(fg, bg)
            row[f"{model}:density"] = selectivity_density_from_loads(fg, fg_length, bg, bg_length)
        rows.append(row)

    return {
        "study": study.key,
        "rows": rows,
        "fg": fg_name,
        "bg": bg_name,
        "fg_length": fg_length,
        "bg_length": bg_length,
        "notes": [n for n in (study.foreground_note, study.background_note) if n],
    }, None


def predictors():
    """What is being ranked, and which direction is better.

    `density` is the figure Known Issue 6 says to compare: a count ratio moves
    with how much background sequence was supplied, and two of these studies
    are scored against a stand-in host of a different size.
    """
    out = {}
    for model in MODELS:
        out[f"density[{model}]"] = (lambda r, m=model: r[f"{m}:density"], True)
        out[f"selectivity[{model}]"] = (lambda r, m=model: r[f"{m}:selectivity"], True)
    # The controls from the 2026-09-16 record.
    out["published fg/bg ratio (control)"] = (lambda r: r["published_fg_bg_ratio"], True)
    out["panel size (control)"] = (lambda r: r["size"], True)
    return out


def report(scored, failures):
    names = predictors()

    for study, reason in failures:
        print(f"\nNOT MEASURED  {study}: {reason}")

    if not scored:
        print("\nNothing could be scored.")
        return

    for include_controls in (True, False):
        tag = "with" if include_controls else "without"
        print(f"\n===== {tag} the negative controls =====")
        selected = {}
        for result in scored:
            rows = [r for r in result["rows"] if include_controls or not r["control"]]
            selected[result["study"]] = rows
            positives = sum(r["label"] for r in rows)
            print(f"  {result['study']}: {len(rows)} panels, {positives} effective")
            for note in result["notes"]:
                print(f"      note: {note}")

        keys = [r["study"] for r in scored]
        print("\n-- Within-study ranking (AUC; 0.5 is no information, 0.0 a perfect inversion)")
        header = f"{'predictor':<38}" + "".join(f"{k.split('_')[-1][:11]:>12}" for k in keys)
        print(header + f"{'mean':>12}")
        for pname, (extract, higher) in names.items():
            cells, values = "", []
            for key in keys:
                rows = selected[key]
                value = auc([extract(r) for r in rows], [r["label"] for r in rows], higher)
                cells += f"{'n/a':>12}" if value is None else f"{value:>12.2f}"
                if value is not None:
                    values.append(value)
            mean = statistics.fmean(values) if values else float("nan")
            print(f"{pname:<38}{cells}{mean:>12.2f}")

        if len(keys) > 1:
            print("\n-- Leave-one-study-out thresholded accuracy")
            print(f"{'predictor':<38}{'held-out acc':>14}{'majority base':>15}")
            everything = [r for key in keys for r in selected[key]]
            majority = max(
                sum(r["label"] for r in everything),
                len(everything) - sum(r["label"] for r in everything),
            ) / len(everything)
            for pname, (extract, higher) in names.items():
                correct = total = 0
                for held in keys:
                    train = [r for key in keys if key != held for r in selected[key]]
                    cut, _ = best_threshold(train, extract, higher)
                    for r in selected[held]:
                        correct += ((extract(r) >= cut) == higher) == r["label"]
                        total += 1
                print(f"{pname:<38}{correct / total:>14.2f}{majority:>15.2f}")
        else:
            print(
                "\n-- Leave-one-study-out thresholded accuracy: NOT MEASURED, "
                "one study scored. Holding out a fold needs at least two."
            )

    print("\n-- Per-panel loads")
    for result in scored:
        print(f"\n  {result['study']}  fg={result['fg']} bg={result['bg']}")
        print(
            f"    {'panel':<18}{'label':>6}{'ctl':>5}"
            + "".join(f"{m.replace('position-dependent', 'pd'):>22}" for m in MODELS)
        )
        for r in result["rows"]:
            cells = "".join(f"{r[f'{m}:fg_load']:10.1f}/{r[f'{m}:bg_load']:<11.1f}" for m in MODELS)
            print(f"    {r['id']:<18}{int(r['label']):>6}{int(r['control']):>5}{cells}")


# ----------------------------------------------------------------------
# Saturation (validation item 4)
# ----------------------------------------------------------------------


def read_pool(path, k):
    """Oligos of length `k` from a step3_df.csv, a step2_df.csv or a plain list."""
    text = pathlib.Path(path).read_text()
    oligos = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        field = line.split(",")[0].strip().strip('"').upper()
        if len(field) == k and set(field) <= set("ACGT"):
            oligos.append(field)
    seen, out = set(), []
    for oligo in oligos:
        if oligo not in seen:
            seen.add(oligo)
            out.append(oligo)
    return out


def saturation(pool_path, host_prefix, k, limit):
    """Spread of the position-dependent host load across a candidate pool.

    Known Issue 6 records that 99.7% of canonical 12-mers occur in the human
    genome, so allowing one mismatch may make nearly every 12-mer's
    neighbourhood present and the host load may stop separating candidates.
    This reports whether it does. A pool whose loads span little is a finding
    about k=12 against that host, not a defect in the model.
    """
    from neoswga.core.occupancy import weighted_site_load
    from neoswga.core.reaction_conditions import ReactionConditions

    oligos = read_pool(pool_path, k)[:limit]
    if not oligos:
        print(f"NOT MEASURED: no {k}-mers found in {pool_path}")
        return
    conditions = ReactionConditions(**REACTION)

    print(f"\n===== saturation of the host load, {len(oligos)} {k}-mers =====")
    print(f"  pool {pool_path}")
    print(f"  host {host_prefix}")
    summary = {}
    for model in MODELS:
        started = time.perf_counter()
        loads = []
        for oligo in oligos:
            loads.append(
                weighted_site_load([oligo], [host_prefix], conditions, MAX_MISMATCHES, None, model)
            )
        took = time.perf_counter() - started
        loads.sort()
        mean = statistics.fmean(loads)
        spread = (loads[-1] / loads[0]) if loads[0] > 0 else float("inf")
        summary[model] = loads
        print(
            f"\n  {model}"
            f"\n    min {loads[0]:.3f}  p10 {loads[len(loads) // 10]:.3f}"
            f"  median {statistics.median(loads):.3f}"
            f"  p90 {loads[-1 - len(loads) // 10]:.3f}  max {loads[-1]:.3f}"
            f"\n    mean {mean:.3f}  sd {statistics.pstdev(loads):.3f}"
            f"  coefficient of variation {statistics.pstdev(loads) / mean:.3f}"
            f"  max/min {spread:.1f}"
            f"\n    {took / len(oligos) * 1000:.2f} ms per oligo"
        )

    base = summary["uniform"]
    modelled = summary["position-dependent"]
    print(
        "\n  Read this as: a coefficient of variation near zero means the host "
        "load has stopped separating candidates at this k, whatever the model. "
        f"uniform {statistics.pstdev(base) / statistics.fmean(base):.3f} against "
        f"position-dependent {statistics.pstdev(modelled) / statistics.fmean(modelled):.3f}."
    )


# ----------------------------------------------------------------------


def _default_table_dir() -> str:
    """Where to put counted tables when nobody said.

    Not the current directory. Run from the repository root as the usage block
    shows, that wrote `mtb_8mer.jf`, `genome_lengths.json` and the rest into the
    checkout, against this project's habit that a run leaves
    `git status --porcelain` unchanged.
    """
    import tempfile

    return os.path.join(tempfile.gettempdir(), "neoswga-mismatch-ranking-tables")


def _refuse_a_table_dir_inside_the_repository(table_dir: str) -> None:
    """A table directory inside the checkout is refused, not silently used.

    Both roots: the script may run from a git worktree, whose own tree is a
    different directory from the one holding the reference genomes.
    """
    resolved = pathlib.Path(table_dir).expanduser().resolve()
    for root in {ROOT.resolve(), DATA_ROOT.resolve()}:
        if resolved == root or root in resolved.parents:
            raise SystemExit(
                f"--table-dir {resolved} is inside the repository at {root}. "
                f"Counting writes .jf tables and genome_lengths.json, and a run "
                f"must leave `git status --porcelain` unchanged. Point it at a "
                f"scratch directory instead; the default is "
                f"{_default_table_dir()}."
            )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--table-dir",
        default=os.environ.get("MISMATCH_RANKING_TABLES") or _default_table_dir(),
        help=(
            "where tables counted here are written, and looked for first. "
            "Defaults to a directory under the system temp dir; a path inside "
            "the repository is refused"
        ),
    )
    parser.add_argument(
        "--build-tables",
        action="store_true",
        help="count any missing table into --table-dir instead of reporting it absent",
    )
    parser.add_argument("--study", action="append", help="restrict to these studies (substring)")
    parser.add_argument("--json", help="write the per-panel rows here")
    parser.add_argument(
        "--saturation",
        action="store_true",
        help="run the host-load saturation check instead of the ranking",
    )
    parser.add_argument("--pool", help="candidate pool for --saturation")
    parser.add_argument("--host-prefix", help="host k-mer table prefix for --saturation")
    parser.add_argument("--k", type=int, default=12, help="oligo length for --saturation")
    parser.add_argument(
        "--limit", type=int, default=2000, help="how many oligos of the pool to score"
    )
    args = parser.parse_args(argv)
    _refuse_a_table_dir_inside_the_repository(args.table_dir)

    if args.saturation:
        if not args.pool or not args.host_prefix:
            parser.error("--saturation needs --pool and --host-prefix")
        saturation(args.pool, args.host_prefix, args.k, args.limit)
        return 0

    cache_path = pathlib.Path(args.table_dir) / "genome_lengths.json"
    length_cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}

    scored, failures = [], []
    for study in STUDIES:
        if args.study and not any(s in study.key for s in args.study):
            continue
        print(f"\nscoring {study.key} ...", flush=True)
        result, reason = score_study(study, args.table_dir, args.build_tables, length_cache)
        if result is None:
            failures.append((study.key, reason))
        else:
            scored.append(result)

    pathlib.Path(args.table_dir).mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(length_cache, indent=1))

    report(scored, failures)

    if args.json:
        pathlib.Path(args.json).write_text(
            json.dumps({"scored": scored, "not_measured": failures}, indent=1)
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
